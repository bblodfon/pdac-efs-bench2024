#' test efs methods
#' Run: `Rscript efs/test_bench.R` (from project root)

renv::load() # load
suppressPackageStartupMessages({
library(mlr3)
library(mlr3extralearners)
library(mlr3proba)
library(mlr3fselect)
library(mlr3tuning)
library(mlr3pipelines)
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
source("efs/callbacks.R")

# SETUP CONFIG ----
cfg = config::get(file = "efs/config.yml")

# dir to save results
save_dir = file.path(cfg$results_dir, cfg$dataset_id)

# Parallelization of XGBoost
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(OMP_THREAD_LIMIT = 1)
Sys.setenv(MKL_NUM_THREADS = 1) # MKL is an Intel-specific thing
# Package-specific settings
try(data.table::setDTthreads(1))
try(RhpcBLASctl::blas_set_num_threads(1))
try(RhpcBLASctl::omp_set_num_threads(1))

# progress bars
options(progressr.enable = TRUE)
handlers(on_missing = "ignore", global = TRUE)
handlers("progress")

# TASK ----
# tasks = readRDS(file = "data/wissel2023/task_list.rds")
# #tasks = readRDS(file = "data/osipov2024/task_list.rds")
# for (task in tasks) {
#   task$id
# }

# test one task
task = cfg$task

# stratify by status
suppressWarnings({
  task$set_col_roles(cols = "status", add_to = "stratum")
})

# CONFIG ----
# Tuning parameters for learners
n_trees = cfg$learner_params$n_trees # RSFs
nrounds = cfg$learner_params$nrounds # xgboost, glmboost and cv_coxboost
folds = cfg$learner_params$folds # for inner resampling
inner_rsmp = rsmp("cv", folds = folds) # for xgboost and glmboost
evals = cfg$learner_params$evals # random search evaluations (glmboost)
eta = cfg$learner_params$eta # learning rate
nrounds_early = cfg$learner_params$nrounds_early # early stopping rounds (for xgboost, should be < nrounds)
max_depth = cfg$learner_params$max_depth # for xgboost

# efs parameters
rfe = fs("rfe")
# rfe = fs("rfe", subset_sizes = c(as.integer((task$n_features)/2), 2)) # for testing
terminator = trm("none") # RFE => never stop iters, defined by `subset_sizes`
store_bmr = cfg$efs$store_bmr
repeats = cfg$efs$repeats # how many subsamples per learner
ratio = cfg$efs$ratio # % of data for training/tuning/fs
init_rsmp = rsmp("subsampling", repeats = repeats, ratio = ratio)
ss_clbk = clbk("mlr3fselect.rfe_subset_sizes") # random subset sizes
#' `size` => how many RFE iters (subset sizes)
#' `shape2` => beta distr parameter, > 1, the larger => more skewed towards less features
#' We assume that `size` < number of task features
#' Some rules of thump for selecting (`size`, `shape2`) for different #features
if (task$n_features > 10000) {
  # e.g. 20000 => (15, 50)
  ss_clbk$state$size = 15
  ss_clbk$state$shape2 = 50
} else if (task$n_features > 5000) {
  # e.g. 7000 => (15, 20)
  ss_clbk$state$size = 15
  ss_clbk$state$shape2 = 30
} else if (task$n_features > 1500) {
  # e.g. 2000 => (15, 20)
  ss_clbk$state$size = 15
  ss_clbk$state$shape2 = 20
} else if (task$n_features > 500) {
  # e.g. 1300 or 800 => (15, 15)
  ss_clbk$state$size = 15
  ss_clbk$state$shape2 = 15
} else if (task$n_features > 100) {
  # e.g. 300 => (10, 5)
  ss_clbk$state$size = 10
  ss_clbk$state$shape2 = 5
} else if (task$n_features > 10) {
  # e.g. 60 => (10, 3)
  ss_clbk$state$size = 10
  ss_clbk$state$shape2 = 3
} else {
  # e.g. 8 => (n-1, 1) => generates set c(n-1,n-2,...,2)
  ss_clbk$state$size = task$n_features - 1
  ss_clbk$state$shape2 = 1
}

# Generate some example sizes
# gen_sizes(n_features = task$n_features, size = ss_clbk$state$size, shape2 = ss_clbk$state$shape2)

# Measure used to optimize and evaluate the learners during the inner resampling process of the training sets
inner_measure = msr("surv.cindex")
# Measure used to evaluate the learners on the test sets generated during the ensemble feature selection process.
measure = msr("surv.cindex")

# print some basic info
cat("# Dataset:", cfg$dataset_id, "\n")
cat("# Task:", task$id, "\n")
cat("# Inner Resampling: ", folds, "-fold CV\n", sep = "")
cat("# Random search #evaluations: ", evals, "\n")
cat("# #Subsamples: ", repeats, "\n")
cat("# Example feature subset sizes (RFE): ", gen_sizes(n_features = task$n_features, size = ss_clbk$state$size, shape2 = ss_clbk$state$shape2), "\n")

# efs parallelization (number of cores for each `ensemble_fselect()` or `embedded_ensemble_fselect()`)
workers = cfg$workers
future::plan("multisession", workers = workers)
cat("#", workers, "workers\n")
cat("---------------------\n")

# EFS (WRAPPED-BASED FS) ----
## RSFs ----
if (cfg$use$RSF) {
  aorsf_lrn =
    po("removeconstants") %>>%
    lrn("surv.aorsf", n_tree = n_trees, control_type = "fast", importance = "permute") |>
    as_learner()
  aorsf_lrn$id = "aorsf"

  rsf_lrns = list(
    lrn("surv.ranger", id = "rsf_logrank", num.trees = n_trees,
        importance = "permutation", splitrule = "logrank"),
    lrn("surv.ranger", id = "rsf_maxstat", num.trees = n_trees,
        importance = "permutation", splitrule = "maxstat"),
    aorsf_lrn
  )

  rsf_clbks = list(
    rsf_logrank = list(clbk("mlr3fselect.one_se_rule"), ss_clbk),
    rsf_maxstat = list(clbk("mlr3fselect.one_se_rule"), ss_clbk),
    aorsf = list(clbk("mlr3fselect.one_se_rule"), ss_clbk)
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

  saveRDS(efs_rsf, file = file.path(save_dir, task$id, "efs_rsf.rds"))
}

## XGBOOST ----
if (cfg$use$XGBoost) {
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
    clbk("mlr3fselect.one_se_rule"), ss_clbk,
    clbk("mlr3fselect.internal_tuning", internal_search_space = internal_ss)
  )

  cat("# Wrapper-based efs with XGBOOST -", length(xgb_lrns), "learners\n")
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

    # add to the xgb list
    efs_list[[learner$id]] = efs_xgb
  }

  # combine xgboost results
  efs_xgb = do.call(c, efs_list)
  saveRDS(efs_xgb, file = file.path(save_dir, task$id, "efs_xgb.rds"))
}

# EFS (EMBEDDED FS) ----
## GLMBoost ----
# search space for tuning via random search
if (cfg$use$GLMBoost) {
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
  efs_glmb = do.call(c, efs_list)
  saveRDS(efs_glmb, file = file.path(save_dir, task$id, "efs_glmb.rds"))
}

## CoxBoost ----
if (cfg$use$CoxBoost) {
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
  saveRDS(efs_coxb, file = file.path(save_dir, task$id, "efs_coxboost.rds"))
}

## CoxLasso ----
if (cfg$use$CoxLasso) {
  # Does internal CV tuning of lambda
  # `s` => which `lambda` to use, `lambda.min` => more features, higher C-index
  # `lambda.1se` => less features, smaller C-index
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
  saveRDS(efs_coxlasso, file = file.path(save_dir, task$id, "efs_coxlasso.rds"))
}

# Combine ALL efs results ----
efs_list = purrr::compact(list(
  if (cfg$use$RSF) efs_rsf,
  if (cfg$use$XGBoost) efs_xgb,
  if (cfg$use$GLMBoost) efs_glmb,
  if (cfg$use$CoxBoost) efs_coxb,
  if (cfg$use$CoxLasso) efs_coxlasso
))

efs_all = do.call(c, efs_list)
saveRDS(efs_all, file = file.path(save_dir, task$id, "efs_all.rds"))
