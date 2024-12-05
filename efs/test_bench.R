#' test efs methods
#' Run: `Rscript efs/test_bench.R` (from project root)

renv::load() # load
suppressPackageStartupMessages({
library(mlr3)
library(mlr3extralearners)
library(mlr3proba)
library(mlr3fselect)
library(mlr3tuning)
library(rpart)
library(ranger)
library(aorsf)
library(xgboost)
library(glmnet)
library(mboost)
library(CoxBoost)
library(progressr)
library(tictoc)
})
source("efs/helpers.R")

# SETUP CONFIG ----
# Parallelization of XGBoost
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(OMP_THREAD_LIMIT = 1)
Sys.setenv(MKL_NUM_THREADS = 1) # MKL is an Intel-specific thing
# Package-specific settings
try(data.table::setDTthreads(1))
try(RhpcBLASctl::blas_set_num_threads(1))
try(RhpcBLASctl::omp_set_num_threads(1))

# Logging
lgr::get_logger("bbotk")$set_threshold("warn")
lgr::get_logger("mlr3" )$set_threshold("warn")

# progress bars
options(progressr.enable = TRUE)
handlers(on_missing = "ignore", global = TRUE)
handlers("progress")

# TASK ----
task = readRDS(file = "data/wissel2023/task_list.rds")$gex # GEX - TCGA
#task = readRDS(file = "data/osipov2024/surv_task_list.rds")$snv # SNV - Wissel
# task = tsk("gbcs") # German Breast Cancer data for testing

# stratify by status
suppressWarnings({
  task$set_col_roles(cols = "status", add_to = "stratum")
})

# PARAMETERS ----
# Tuning parameters for learners
n_trees = 500 # RSFs
nrounds = 500 # xgboost, glmboost and cv_coxboost
folds = 5 # for inner resampling
inner_rsmp = rsmp("cv", folds = folds)
evals = 25 # random search evaluations (glmboost)
eta = 0.1 # learning rate
nrounds_early = 42 # early stopping rounds (for xgboost, should be < nrounds)
max_depth = 6 # for xgboost

# efs parameters
repeats = 100 # how many subsamples
ratio = 0.8 # % of data for training/tuning/fs
init_rsmp = rsmp("subsampling", repeats = repeats, ratio = ratio)
terminator = trm("none") # RFE => never stop until you've reached `n_features`
n_features = 2 # RFE => stopping criterion (run up to this number of features)
feature_fraction = 0.8 # RFE => keep this % of features in each RFE iteration
# feature subset sizes for RFE, given we go from n => 2, keeping 80% each time
n = task$n_features
subset_sizes = unique(floor(cumprod(c(n, rep(feature_fraction, log(n_features / n) / log(feature_fraction))))))
num_subsets = 15 # RFE iterations => we don't want to do more than this due to RAM issues
if (length(subset_sizes) > num_subsets) {
  # thin the feature subsets to 15 values in total
  indx = unique(round(seq.int(1, length(subset_sizes), length.out = num_subsets)))
  subset_sizes = subset_sizes[indx]
}
rfe = fs("rfe", subset_sizes = subset_sizes)
# for testing
# rfe = fs("rfe", subset_sizes = c(as.integer((task$n_features)/2), as.integer((task$n_features)/4), 2))
# rfe = fs("rfe", feature_number = 1, n_features = 1)
store_bmr = FALSE

# Measure used to optimize and evaluate the learners during the inner resampling process of the training sets
inner_measure = msr("surv.cindex")
# Measure used to evaluate the learners on the test sets generated during the ensemble feature selection process.
measure = msr("surv.cindex")

# print some basic parameters
cat("# Task:", task$id, "\n")
cat("# ", folds, "-fold CV\n", sep = "")
cat("# Random search evals =", evals, "\n")
cat("#", repeats, "subsamples\n")
cat("#", length(subset_sizes), "feature subset sizes (RFE): ", subset_sizes, "\n")

# efs parallelization (number of cores for each `ensemble_fselect()` or `embedded_ensemble_fselect()`)
workers = 10
future::plan("multisession", workers = workers)
cat("#", workers, "workers\n")
cat("---------------------\n")

if (FALSE){
# EFS (WRAPPED-BASED FS) ----
## Tree ----
# minbucket = round(minsplit/3) = 5 obs/leaf
tree_lrn = list(
  lrn("surv.rpart", id = "tree", minsplit = 15)
)

cat("# Wrapper-based efs with SURVIVAL TREE")
tic()
set.seed(42) # reproduce: same subsampling
suppressWarnings({
  efs_tree = ensemble_fselect(
    fselector = rfe,
    task = task,
    learners = tree_lrn,
    init_resampling = init_rsmp,
    inner_resampling = inner_rsmp,
    inner_measure = measure,
    measure = measure,
    terminator = terminator,
    callbacks = list(
      tree = clbks("mlr3fselect.one_se_rule")
    ),
    store_benchmark_result = store_bmr,
    store_models = FALSE
  )
})
toc()
rm_imp(efs_tree)
saveRDS(efs_tree, file = paste0("efs/", task$id, "/efs_tree.rds"))
}
## RSFs ----
rsf_lrns = list(
  lrn("surv.ranger", id = "rsf_logrank", num.trees = n_trees,
      importance = "permutation", splitrule = "logrank"),
  lrn("surv.ranger", id = "rsf_maxstat", num.trees = n_trees,
      importance = "permutation", splitrule = "maxstat"),
  lrn("surv.aorsf", id = "aorsf", n_tree = n_trees, control_type = "fast",
      importance = "permute")
)

rsf_clbks = list(
  rsf_logrank = clbks("mlr3fselect.one_se_rule"),
  rsf_maxstat = clbks("mlr3fselect.one_se_rule"),
  aorsf = clbks("mlr3fselect.one_se_rule")
)

cat("# Wrapper-based efs with RANDOM SURVIVAL FORESTS -", length(rsf_lrns), "learners\n")
tic()
set.seed(42) # reproduce: same subsampling
suppressWarnings({
  efs_rsf = ensemble_fselect(
    fselector = rfe,
    task = task,
    learners = rsf_lrns,
    init_resampling = init_rsmp,
    inner_resampling = rsmp("insample"), # use all training data
    inner_measure = msr("oob_error"), # 1 - C-index
    measure = measure,
    terminator = terminator,
    callbacks = rsf_clbks,
    store_benchmark_result = store_bmr,
    store_models = FALSE
  )
})
toc()
rm_imp(efs_rsf)
oob_to_cindex_convert(efs_rsf)

saveRDS(efs_rsf, file = paste0("efs/", task$id, "/efs_rsf.rds"))

## XGBOOST ----
xgb_params = list(
  eta = eta,
  max_depth = max_depth,
  nrounds = nrounds,
  nrounds_early = nrounds_early
)

xgb_lrns = list(
  create_xgb_lrn(id = "xgb_cox", params = xgb_params),
  #create_xgb_lrn(id = "xgb_aft_norm", params = xgb_params, aft_loss = "normal"),
  #create_xgb_lrn(id = "xgb_aft_extreme", params = xgb_params, aft_loss = "extreme"),
  create_xgb_lrn(id = "xgb_aft_log", params = xgb_params, aft_loss = "logistic")
)

# search space for internal tuning
internal_ss = ps(
  nrounds = p_int(upper = nrounds, aggr = function(x) as.integer(mean(unlist(x))))
)
# per learner
xgb_clbks = list(
  clbk("mlr3fselect.one_se_rule"),
  clbk("mlr3fselect.internal_tuning", internal_search_space = internal_ss)
)

cat("# Wrapper-based efs with XGBOOST -", length(xgb_lrns), "learners\n")
# ALL TOGETHER => NEEDS MEMORY
# tic()
# set.seed(42) # reproduce: same subsampling
# suppressWarnings({
#   efs_xgb = ensemble_fselect(
#     fselector = rfe,
#     task = task,
#     learners = xgb_lrns,
#     init_resampling = init_rsmp,
#     inner_resampling = inner_rsmp,
#     inner_measure = measure,
#     measure = measure,
#     terminator = terminator,
#     callbacks = list(
#       xgb_cox = xgb_clbks,
#       xgb_aft_norm = xgb_clbks,
#       xgb_aft_extreme = xgb_clbks,
#       xgb_aft_log = xgb_clbks
#     ),
#     store_benchmark_result = store_bmr,
#     store_models = FALSE
#   )
# })
# toc()
# rm_imp(efs_xgb)

efs_list = list()
for (learner in xgb_lrns) {
  cat("#", learner$id, "\n")

  xgb_clbk = list(xgb_clbks)
  names(xgb_clbk) = learner$id

  tic()
  set.seed(42) # reproduce: same subsampling
  efs_xgb = ensemble_fselect(
    fselector = rfe,
    task = task,
    learners = list(learner),
    init_resampling = init_rsmp,
    inner_resampling = inner_rsmp,
    inner_measure = measure,
    measure = measure,
    terminator = terminator,
    callbacks = xgb_clbk,
    store_benchmark_result = store_bmr,
    store_models = FALSE
  )
  toc()
  rm_imp(efs_xgb)

  # add to the xgb list
  efs_list[[learner$id]] = efs_xgb
}

# combine xgboost results
efs_xgb = cmb_efs(efs_list, features = task$feature_names, measure = measure, inner_measure = inner_measure)

saveRDS(efs_xgb, file = paste0("efs/", task$id, "/efs_xgb.rds"))

# EFS (EMBEDDED FS) ----
## glmboost ----
# search space for tuning via random search
glmb_ss = ps(
  mstop = p_int(10, nrounds), # boosting rounds
  nu = p_dbl(0, eta) # learning rate
)

glmb_params = list(
  search_space = glmb_ss,
  inner_rsmp = inner_rsmp,
  measure = measure,
  evals = evals
)

# we ran them individually below
glmb_lrns = list(
  create_glmb_at(id = "glmb_cox", family = "coxph", params = glmb_params),
  create_glmb_at(id = "glmb_loglog", family = "loglog", params = glmb_params)#,
  #create_glmb_at(id = "glmb_weibull", family = "weibull", params = glmb_params),
  #create_glmb_at(id = "glmb_lognorm", family = "lognormal", params = glmb_params)
)

cat("# Embedded efs with GLMBOOST -", length(glmb_lrns), "learners\n")
efs_list = list()
for (learner in glmb_lrns) {
  cat("#", learner$id, "\n")

  tic()
  set.seed(42) # reproduce: same subsampling
  efs_glmb = embedded_ensemble_fselect(
    task = task,
    learners = list(learner),
    init_resampling = init_rsmp,
    measure = measure,
    store_benchmark_result = store_bmr
  )
  toc()
  # remove rows which had 0 features selected (due to whatever reason) from embedded efs
  rm_zero_feat(efs_glmb)

  # add to the glmb list
  efs_list[[learner$id]] = efs_glmb
}

# combine glmboost results
efs_glmb = cmb_efs(efs_list, features = task$feature_names, measure = measure, inner_measure = NULL)
saveRDS(efs_glmb, file = paste0("efs/", task$id, "/efs_glmb.rds"))

## CoxBoost ----
# Does internal CV tuning of the nrounds
coxboost = lrn("surv.cv_coxboost", id = "coxboost", standardize = TRUE,
               return.score = FALSE, penalty = "optimCoxBoostPenalty",
               maxstepno = nrounds, K = folds)

cat("# Embedded efs with COXBOOST learner\n")
tic()
set.seed(42) # reproduce: same subsampling
efs_coxb = embedded_ensemble_fselect(
  task = task,
  learners = list(coxboost),
  init_resampling = init_rsmp,
  measure = measure,
  store_benchmark_result = store_bmr
)
toc()
# remove rows which had 0 features selected (due to whatever reason) from embedded efs
rm_zero_feat(efs_coxb)
saveRDS(efs_coxb, file = paste0("efs/", task$id, "/efs_coxb.rds"))

## CoxLasso ----
# Does internal CV tuning of lambda (s => which lambda, lambda.min => more features)
coxlasso = lrn("surv.cv_glmnet", id = "coxlasso", standardize = TRUE, alpha = 1,
               nfolds = folds, type.measure = "C", s = "lambda.min")

cat("# Embedded efs with COXLASSO learner\n")
tic()
set.seed(42) # reproduce: same subsampling
efs_coxlasso = embedded_ensemble_fselect(
  task = task,
  learners = list(coxlasso),
  init_resampling = init_rsmp,
  measure = measure,
  store_benchmark_result = store_bmr
)
toc()
# remove rows which had 0 features selected (due to whatever reason) from embedded efs
rm_zero_feat(efs_coxlasso)
saveRDS(efs_coxlasso, file = paste0("efs/", task$id, "/efs_coxlasso.rds"))

# Combine ALL efs results ----
efs_list = list(
  efs_rsf,
  efs_xgb,
  efs_glmb,
  efs_coxb,
  efs_coxlasso
)

efs_all = cmb_efs(efs_list, features = task$feature_names, measure = measure, inner_measure = inner_measure)
saveRDS(efs_all, file = paste0("efs/", task$id, "/efs_all.rds"))
