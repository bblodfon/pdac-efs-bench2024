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
n_features = 2 # RFE => stopping criterion (run up to this number of features)
feature_fraction = 0.8 # RFE => keep 80% of features in each RFE iteration
rfe = fs("rfe", n_features = n_features, feature_fraction = feature_fraction)
rfe_test = fs("rfe", subset_sizes = c(as.integer((task$n_features)/2), 2)) # for testing
measure = msr("surv.cindex") # measure used in 1) RFE optimization 2) scoring the
# test sets of `init_rsmp` 3) for tuning, using the train set's of `init_rsmp`
terminator = trm("none") # RFE => never stop until you've reached `n_features`
store_bmr = FALSE

# efs parallelization (number of cores for each `ensemble_fselect()` or `embedded_ensemble_fselect()`)
workers = 10
future::plan("multicore", workers = workers)
sprintf("%s workers", workers)

# LEARNERS ----
## Tree ----
# minbucket = round(minsplit/3) = 5 obs/leaf
tree_lrn = list(
  lrn("surv.rpart", id = "tree", minsplit = 15)
)

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

## XGBOOST ----
xgb_params = list(
  eta = eta,
  max_depth = max_depth,
  nrounds = nrounds,
  nrounds_early = nrounds_early
)

# convenience function to create an xgboost learner
create_xgb_lrn = function(id, params = xgb_params, aft_loss = NULL) {
  checkmate::assert_list(params, types = "numeric", len = 4, any.missing = FALSE, null.ok = FALSE)
  checkmate::check_subset(
    x = names(params),
    choices = c("eta", "max_depth", "rounds", "nrounds_early"),
    empty.ok = FALSE
  )

  eta = params$eta
  max_depth = params$max_depth
  nrounds = params$nrounds
  nrounds_early = params$nrounds_early

  if (is.null(aft_loss)) {
    lrn("surv.xgboost.cox", id = id,
        booster = "gbtree", tree_method = "hist",
        eta = eta, max_depth = max_depth, nrounds = nrounds,
        early_stopping_rounds = nrounds_early, validate = "test")
  } else {
    lrn("surv.xgboost.aft", id = id,
        booster = "gbtree", tree_method = "hist",
        eta = eta, max_depth = max_depth, nrounds = nrounds,
        early_stopping_rounds = nrounds_early, validate = "test",
        aft_loss_distribution_scale = 1, # sigma = 1
        aft_loss_distribution = aft_loss)
  }
}

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

# convenience function to create a glmboost AutoTuner learner
create_glmb_at = function(id, family = "coxph", params = glmb_params) {
  checkmate::assert_list(params, len = 4, any.missing = FALSE, null.ok = FALSE)
  checkmate::check_subset(
    x = names(params),
    choices = c("search_space", "inner_rsmp", "measure", "evals")
  )

  search_space = params$search_space
  inner_rsmp = params$inner_rsmp
  measure = params$measure
  evals = params$evals

  learner = lrn("surv.glmboost", id = id, family = family, center = TRUE)

  suppressWarnings({
    auto_tuner(
      learner = learner,
      search_space = search_space,
      resampling = inner_rsmp,
      measure = measure,
      tuner = tnr("random_search"),
      term_evals = evals,
      store_tuning_instance = FALSE, # archive
      store_benchmark_result = FALSE, # benchmark result of inner resamplings
      store_models = FALSE,
      callbacks = NULL # no 1se rule exists yet
    )
  })
}

glmb_lrns = list(
  create_glmb_at(id = "glmb_cox", family = "coxph"),
  create_glmb_at(id = "glmb_loglog", family = "loglog"),
  create_glmb_at(id = "glmb_weibull", family = "weibull"),
  create_glmb_at(id = "glmb_lognorm", family = "lognormal")
)

## CoxBoost ----
# Does internal CV tuning of the nrounds
coxboost = lrn("surv.cv_coxboost", id = "coxboost", standardize = TRUE,
               return.score = FALSE, penalty = "optimCoxBoostPenalty",
               maxstepno = nrounds, K = folds)
## CoxLasso ----
# Does internal CV tuning of lambda (choose lambda based on the `1se` rule)
coxlasso = lrn("surv.cv_glmnet", id = "coxlasso", standardize = TRUE, alpha = 1,
               nfolds = folds, type.measure = "C", s = "lambda.1se")

# EFS runs ----
if (FALSE) {
## Tree ----
set.seed(42) # reproduce: same subsampling
sprintf("# Wrapper-based efs with %s - %i learner(s)", "SURVIVAL TREE", length(tree_lrn))
tic()
suppressWarnings({
efs_tree = ensemble_fselect(
  #fselector = rfe_test,
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
saveRDS(efs_tree, paste0(task$id, "_efs_tree.rds"))

## RSFs ----
sprintf("# Wrapper-based efs with %s - %i learner(s)", "RANDOM SURVIVAL FORESTS", length(rsf_lrns))
tic()
set.seed(42) # reproduce: same subsampling
suppressWarnings({
efs_rsf = ensemble_fselect(
  #fselector = rfe_test,
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
efs_rsf$.__enclos_env__$private$.result[, surv.cindex := 1 - oob_error]
efs_rsf$.__enclos_env__$private$.result$oob_error = NULL
efs_rsf$.__enclos_env__$private$.measure_id = "surv.cindex"
efs_rsf$.__enclos_env__$private$.minimize = FALSE

saveRDS(efs_rsf, paste0(task$id, "_efs_rsf.rds"))

## XGBOOST ----
sprintf("# Wrapper-based efs with %s - %i learner(s)", "XGBOOST", length(xgb_lrns))
tic()
set.seed(42) # reproduce: same subsampling
suppressWarnings({
efs_xgb = ensemble_fselect(
  #fselector = rfe_test,
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
saveRDS(efs_xgb, paste0(task$id, "_efs_xgb.rds"))
}
## glmboost ----
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
saveRDS(efs_glmb, paste0(task$id, "_efs_glmb.rds"))

## CoxBoost ----
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
saveRDS(efs_coxb, paste0(task$id, "_efs_coxb.rds"))

## CoxLasso ----
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
saveRDS(efs_coxlasso, paste0(task$id, "_efs_coxlasso.rds"))

# Combine results ----
res_list = list(
  efs_tree$result,
  efs_rsf$result,
  efs_xgb$result,
  efs_glmb$result,
  efs_coxb$result,
  efs_coxlasso$result
)
result = data.table::rbindlist(l = res_list, fill = TRUE)
result$importance = NULL # no need to keep this

efs_res = EnsembleFSResult$new(
  result = result,
  features = task$feature_names,
  measure_id = "surv.cindex",
  minimize = FALSE
)
saveRDS(efs_res, paste0(task$id, "_efs_res.rds"))
