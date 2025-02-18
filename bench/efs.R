#' This script runs the ensemble feature selection procedure, for a specific
#' dataset, omic and resampling id (train set)
#'
#' Execute: `Rscript bench/efs.R <dataset_id> <omic_id> <rsmp_id>` (from project root)
#' e.g. `Rscript bench/efs.R wissel2023 gex 42`

# use `renv`
#renv::load()

# CMD args ----
args = commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript bench/efs.R <dataset_id> <omic_id> <rsmp_id>")
}

dataset_id = args[1]
omic_id = args[2]
rsmp_id = as.integer(args[3])

# perform checks
library(checkmate)
dataset_path = file.path("data", dataset_id)
assert_directory(dataset_path)

task_list = readRDS(file = file.path(dataset_path, "task_list.rds"))
omic_ids = names(task_list)
assert_subset(omic_id, omic_ids)

subsampling = readRDS(file = file.path(dataset_path, "subsampling.rds"))
assert_number(rsmp_id, lower = 1, upper = subsampling$iters)

# make directory for results if it doesn't already exist
res_path = file.path("bench", "efs", dataset_id, omic_id)
if (!test_directory_exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
}

# LIBRARIES ----
suppressPackageStartupMessages({
  library(mlr3)
  library(mlr3extralearners)
  library(mlr3proba)
  library(mlr3fselect)
  library(mlr3tuning)
  library(mlr3pipelines)
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

# TASK ----
task = task_list[[omic_id]]$clone()
# run efs only using the train set of this particular subsampling
task$filter(rows = subsampling$train_set(rsmp_id))
task$set_col_roles(cols = "status", add_to = "stratum") # stratify by status

# CONFIG ----
# get config parameters from file
cfg = config::get(file = "bench/config.yml")

# verbose message printing
print_msg = function(...) {
  if (cfg$verbose > 0) cat(...)
}

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

# result efs file with ALL results
efs_path = file.path(res_path, paste0("efs_", rsmp_id, ".rds"))

# if results exist, don't run again
if (file.exists(efs_path) && !cfg$overwrite) {
  cat("Exiting... results already exist and we don't overwrite\n")
  quit()
}

# print some basic info about this run
print_msg("Dataset:", dataset_id,
          "\nOmic:", omic_id,
          "\nResampling Iteration:", rsmp_id,
          "\n")

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

# Measure used to optimize and evaluate the learners during the inner resampling
# process in the training sets
inner_measure = msr("surv.cindex")
# Measure used to evaluate the learners on the test sets generated during the
# ensemble feature selection process
measure = msr("surv.cindex")

# Parallelization
workers = cfg$workers
future::plan("multicore", workers = workers)

# keep execution times
efs_times = tibble::tibble(id = character(), time = numeric())

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

  print_msg("# Wrapper-based efs with RANDOM SURVIVAL FORESTS -", length(rsf_lrns), "learners\n")
  start_time = Sys.time()
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
  stop_time = Sys.time()
  time_diff = round(as.double(stop_time - start_time, units = "secs"), digits = 2)
  print_msg(time_diff, "secs\n")
  efs_times = efs_times |> tibble::add_row(id = "rsf", time = time_diff)
  # saveRDS(efs_rsf, file = file.path(res_path, "efs_rsf.rds"))
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

  print_msg("# Wrapper-based efs with XGBOOST -", length(xgb_lrns), "learners\n")
  efs_list = list()
  for (learner in xgb_lrns) {
    print_msg("#", learner$id, "\n")

    xgb_clbk = list(xgb_clbks)
    names(xgb_clbk) = learner$id

    start_time = Sys.time()
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
    stop_time = Sys.time()
    time_diff = round(as.double(stop_time - start_time, units = "secs"), digits = 2)
    print_msg(time_diff, "secs\n")
    efs_times = efs_times |> tibble::add_row(id = learner$id, time = time_diff)

    # add to the xgb list
    efs_list[[learner$id]] = efs_xgb
  }

  # combine xgboost results
  efs_xgb = do.call(c, efs_list)
  # saveRDS(efs_xgb, file = file.path(res_path, "efs_xgb.rds"))
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
    create_glmb_at(id = "glmb_loglog", family = "loglog", params = glmb_params)
  )

  print_msg("# Embedded efs with GLMBOOST -", length(glmb_lrns), "learners\n")
  efs_list = list()
  for (learner in glmb_lrns) {
    print_msg("#", learner$id, "\n")

    start_time = Sys.time()
    set.seed(42) # reproduce: same subsampling
    suppressWarnings({
      efs_glmb = embedded_ensemble_fselect(
        task = task,
        learners = list(learner),
        init_resampling = init_rsmp,
        measure = measure,
        store_benchmark_result = store_bmr
      )
    })
    stop_time = Sys.time()
    time_diff = round(as.double(stop_time - start_time, units = "secs"), digits = 2)
    print_msg(time_diff, "secs\n")
    efs_times = efs_times |> tibble::add_row(id = learner$id, time = time_diff)

    # remove rows which had 0 features selected (due to whatever reason) from embedded efs
    rm_zero_feat(efs_glmb)

    # add to the glmb list
    efs_list[[learner$id]] = efs_glmb
  }

  # combine glmboost results
  efs_glmb = do.call(c, efs_list)
  # saveRDS(efs_glmb, file = file.path(res_path, "efs_glmb.rds"))
}

## CoxBoost ----
if (cfg$use$CoxBoost) {
  # Does internal CV tuning of the nrounds
  coxboost = lrn("surv.cv_coxboost", id = "coxboost", standardize = TRUE,
                 return.score = FALSE, penalty = "optimCoxBoostPenalty",
                 maxstepno = nrounds, K = folds)

  print_msg("# Embedded efs with COXBOOST learner\n")
  start_time = Sys.time()
  set.seed(42) # reproduce: same subsampling
  suppressWarnings({
    efs_coxb = embedded_ensemble_fselect(
      task = task,
      learners = list(coxboost),
      init_resampling = init_rsmp,
      measure = measure,
      store_benchmark_result = store_bmr
    )
  })
  stop_time = Sys.time()
  time_diff = round(as.double(stop_time - start_time, units = "secs"), digits = 2)
  print_msg(time_diff, "secs\n")
  efs_times = efs_times |> tibble::add_row(id = "coxboost", time = time_diff)

  # remove rows which had 0 features selected (due to whatever reason) from embedded efs
  rm_zero_feat(efs_coxb)
  #saveRDS(efs_coxb, file = file.path(res_path, "efs_coxboost.rds"))
}

## CoxLasso ----
if (cfg$use$CoxLasso) {
  # Does internal CV tuning of lambda
  # `s` => which `lambda` to use, `lambda.min` => more features, higher C-index
  # `lambda.1se` => less features, smaller C-index
  coxlasso = lrn("surv.cv_glmnet", id = "coxlasso", standardize = TRUE, alpha = 1,
                 nfolds = folds, type.measure = "C", s = "lambda.min")

  print_msg("# Embedded efs with COXLASSO learner\n")
  start_time = Sys.time()
  set.seed(42) # reproduce: same subsampling
  efs_coxlasso = embedded_ensemble_fselect(
    task = task,
    learners = list(coxlasso),
    init_resampling = init_rsmp,
    measure = measure,
    store_benchmark_result = store_bmr
  )
  stop_time = Sys.time()
  time_diff = round(as.double(stop_time - start_time, units = "secs"), digits = 2)
  print_msg(time_diff, "secs\n")
  efs_times = efs_times |> tibble::add_row(id = "coxlasso", time = time_diff)

  # remove rows which had 0 features selected (due to whatever reason) from embedded efs
  rm_zero_feat(efs_coxlasso)
  # saveRDS(efs_coxlasso, file = file.path(res_path, "efs_coxlasso.rds"))
}

# SAVE results ----
# combine all efs results in one object
efs_list = purrr::compact(list(
  if (cfg$use$RSF) efs_rsf,
  if (cfg$use$XGBoost) efs_xgb,
  if (cfg$use$GLMBoost) efs_glmb,
  if (cfg$use$CoxBoost) efs_coxb,
  if (cfg$use$CoxLasso) efs_coxlasso
))

efs_all = do.call(c, efs_list)
saveRDS(efs_all, file = efs_path)

# save timings
times_file = file.path(res_path, paste0("times_", rsmp_id, ".csv"))
readr::write_csv(efs_times, file = times_file)

# report efs total time
cat("Total time:", sum(efs_times$time), "secs\n")
