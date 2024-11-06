# test efs methods on mRNA data by Wissel
renv::load()
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
# mRNA only
task = readRDS(file = "data/wissel2023/task_list.rds")$gex
# stratify by status
task$set_col_roles(cols = "status", add_to = "stratum")

# PARAMETERS ----
# Tuning parameters for learners
n_trees = 500 # for RSFs
nrounds = 500 # for xgboost, glmboost and cv_coxboost
folds = 5 # for inner resampling
inner_rsmp = rsmp("cv", folds = folds)
evals = 50 # how many tuning evaluations to do with random search (for glmboost)
eta = 0.1 # learning rate
nrounds_early = 42 # early stopping rounds (for xgboost, should be < nrounds)
max_depth = 6 # for xgboost

# efs parameters
repeats = 100 # how many subsamples
ratio = 0.9 # % of data for training/tuning/fs
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
print(paste0(length(subset_sizes), " feature subset sizes (RFE): ",  paste0(subset_sizes, collapse = ",")))
rfe = fs("rfe", subset_sizes = subset_sizes)
# for testing
# rfe = fs("rfe", subset_sizes = c(as.integer((task$n_features)/2), 2))
# measure used in
## 1) RFE optimization
## 2) score the test sets of `init_rsmp`
## 3) during fs via RFE (on the train set's of `init_rsmp`)
measure = msr("surv.cindex")
store_bmr = FALSE

# efs parallelization (number of cores for each `ensemble_fselect()` or `embedded_ensemble_fselect()`)
workers = 10
future::plan("multicore", workers = workers)
sprintf("%s workers", workers)

# EFS runs ----
if (FALSE) {
## Tree ----
# minbucket = round(minsplit/3) = 5 obs/leaf
tree_lrn = list(
  lrn("surv.rpart", id = "tree", minsplit = 15)
)

set.seed(42) # reproduce: same subsampling
sprintf("# Wrapper-based efs with %s - %i learner(s)", "SURVIVAL TREE", length(tree_lrn))
tic()
suppressWarnings({
  efs_tree = ensemble_fselect(
    fselector = rfe,
    task = task,
    learners = tree_lrn,
    init_resampling = init_rsmp,
    inner_resampling = inner_rsmp,
    measure = measure,
    terminator = terminator,
    callbacks = list(
      clbks("mlr3fselect.one_se_rule")
    ),
    store_benchmark_result = store_bmr,
    store_models = FALSE
  )
})
toc()
rm_imp(efs_tree)

saveRDS(efs_tree, file = paste0("efs/", task$id, "/efs_tree.rds"))

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
  clbks("mlr3fselect.one_se_rule"),
  clbks("mlr3fselect.one_se_rule"),
  clbks("mlr3fselect.one_se_rule")
)

sprintf("# Wrapper-based efs with %s - %i learner(s)", "RANDOM SURVIVAL FORESTS", length(rsf_lrns))
tic()
set.seed(42) # reproduce: same subsampling
suppressWarnings({
  efs_rsf = ensemble_fselect(
    fselector = rfe,
    task = task,
    learners = rsf_lrns,
    init_resampling = init_rsmp,
    inner_resampling = rsmp("insample"), # use all training data
    measure = msr("oob_error"), # 1 - C-index
    terminator = terminator,
    callbacks = rsf_clbks,
    store_benchmark_result = store_bmr,
    store_models = FALSE
  )
})
toc()
# hack: C-index = 1 - OOB_ERROR
oob_to_cindex_convert(efs_rsf)
rm_imp(efs_rsf)

saveRDS(efs_rsf, file = paste0("efs/", task$id, "/efs_rsf.rds"))

## XGBOOST ----
xgb_params = list(
  eta = eta,
  max_depth = max_depth,
  nrounds = nrounds,
  nrounds_early = nrounds_early
)

xgb_lrns = list(
  create_xgb_lrn(id = "xgb_cox"),
  create_xgb_lrn(id = "xgb_aft_norm", aft_loss = "normal"),
  create_xgb_lrn(id = "xgb_aft_extreme", aft_loss = "extreme"),
  create_xgb_lrn(id = "xgb_aft_log", aft_loss = "logistic")
)

# search space for internal tuning
internal_ss = ps(
  nrounds = p_int(upper = nrounds, aggr = function(x) as.integer(mean(unlist(x))))
)
# per learner
xgb_clbk = list(
  clbk("mlr3fselect.one_se_rule"),
  clbk("mlr3fselect.internal_tuning", internal_search_space = internal_ss)
)
xgb_clbks = list(xgb_clbk, xgb_clbk, xgb_clbk, xgb_clbk)

sprintf("# Wrapper-based efs with %s - %i learner(s)", "XGBOOST", length(xgb_lrns))
tic()
set.seed(42) # reproduce: same subsampling
suppressWarnings({
  efs_xgb = ensemble_fselect(
    fselector = rfe,
    task = task,
    learners = xgb_lrns,
    init_resampling = init_rsmp,
    inner_resampling = inner_rsmp,
    measure = measure,
    terminator = terminator,
    callbacks = xgb_clbks,
    store_benchmark_result = store_bmr,
    store_models = FALSE
  )
})
toc()
rm_imp(efs_xgb)

saveRDS(efs_xgb, file = paste0("efs/", task$id, "/efs_xgb.rds"))
}
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

glmb_lrns = list(
  create_glmb_at(id = "glmb_cox", family = "coxph", params = glmb_params),
  create_glmb_at(id = "glmb_loglog", family = "loglog", params = glmb_params),
  create_glmb_at(id = "glmb_weibull", family = "weibull", params = glmb_params),
  create_glmb_at(id = "glmb_lognorm", family = "lognormal", params = glmb_params)
)

sprintf("# Embedded efs with %s - %i learner(s)", "GLMBOOST", length(glmb_lrns))
tic()
set.seed(42) # reproduce: same subsampling
efs_glmb = embedded_ensemble_fselect(
  task = task,
  learners = glmb_lrns,
  init_resampling = init_rsmp,
  measure = measure,
  store_benchmark_result = store_bmr
)
toc()
# remove rows which had 0 features selected (due to whatever reason) from embedded efs
rm_zero_feat(efs_glmb)
saveRDS(efs_glmb, file = paste0("efs/", task$id, "/efs_glmb.rds"))

## CoxBoost ----
# Does internal CV tuning of the nrounds
coxboost = lrn("surv.cv_coxboost", id = "coxboost", standardize = TRUE,
               return.score = FALSE, penalty = "optimCoxBoostPenalty",
               maxstepno = nrounds, K = folds)

sprintf("# Embedded efs with %s - %i learner(s)", "COXBOOST", 1)
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
# Does internal CV tuning of lambda (s => which lambda)
coxlasso = lrn("surv.cv_glmnet", id = "coxlasso", standardize = TRUE, alpha = 1,
               nfolds = folds, type.measure = "C", s = "lambda.min")

sprintf("# Embedded efs with %s - %i learner(s)", "COXLASSO", 1)
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

# Combine results ----
res_list = list(
  #efs_tree$result,
  #efs_rsf$result,
  #efs_xgb$result,
  efs_glmb$result,
  efs_coxb$result,
  efs_coxlasso$result
)
result = data.table::rbindlist(l = res_list, fill = TRUE)

efs_all = EnsembleFSResult$new(
  result = result,
  features = task$feature_names,
  measure_id = "surv.cindex",
  minimize = FALSE
)
saveRDS(efs_all, file = paste0("efs/", task$id, "/efs_all.rds"))
