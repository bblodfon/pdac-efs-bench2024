#' FEATURE SELECTION (PER OMIC AND SUBSAMPLING)
#'
#' This script performs feature selection per omic and per subsampling iteration
#' (100 in total) - either using CoxLasso or the ensemble feature selection (efs)
#' method on two PDAC datasets.
#'
#' The efs, being computationally intensive, has been executed in a separate
#' `run_efs.sh` script, so we just load the result object here and decide on the
#' number of features via the estimated Pareto front).
#'
#' Execute: `Rscript bench/run_fs.R` (from project root)

# use `renv`
#renv::load()

# Load required libraries
suppressPackageStartupMessages({
  library(mlr3)
  library(mlr3fselect)
  library(mlr3extralearners)
  library(mlr3proba)
  library(fastVoteR)
  library(glmnet)
  library(checkmate)
  library(readr)
  library(dplyr)
  library(tibble)
  library(future.apply)
  library(progressr)
})
source("bench/helpers.R")

# Set parallel execution
plan("multicore", workers = 30)

# Enable progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers("progress")

# Define datasets
dataset_ids = c("wissel2023", "osipov2024")

# Construct a parameter grid (all combinations of dataset_id, omic_id, rsmp_id)
grid_list = lapply(dataset_ids, function(dataset_id) {
  dataset_path = file.path("data", dataset_id)
  assert_directory(dataset_path)

  # Load omic task IDs
  omic_ids = readr::read_csv(file.path(dataset_path, "omic_ids.csv"),
                             col_names = FALSE, show_col_types = FALSE)[[1L]]

  # Load subsampling iterations
  subsampling = readRDS(file.path(dataset_path, "subsampling.rds"))
  rsmp_ids = seq(subsampling$iters)

  expand.grid(dataset_id = dataset_id, omic_id = omic_ids, rsmp_id = rsmp_ids,
              stringsAsFactors = FALSE)
})

# Flatten the list into a dataframe
grid_df = do.call(rbind, grid_list)

# Parallelized function for feature selection
feature_selection = function(params, p) {
  set.seed(42)
  dataset_id = params$dataset_id
  omic_id = params$omic_id
  rsmp_id = params$rsmp_id
  print_params = list(dataset_id = dataset_id, omic_id = omic_id, rsmp_id = rsmp_id)

  #' Notify progress via `p = progressr::progressor()`
  p(sprintf("Dataset: %s, Omic: %s, Subsampling Iter: %i", dataset_id, omic_id, rsmp_id))

  # Load ensemble feature selection (EFS) results
  file_name = file.path("bench", "efs", dataset_id, omic_id, paste0("efs_", rsmp_id, ".rds"))
  efs = readRDS(file_name)

  # SELECT FEATURES via EFS (all models)
  ## get number of features via the Pareto method
  efs_nfeats = get_nfeats(efs, lrn_ids = NULL, type = "estimated",
                          upper_bound = "max_efs", print_params = print_params)

  efs_feats = efs$feature_ranking(
    method = "sav",
    use_weights = TRUE,
    committee_size = efs_nfeats
  )[["feature"]]

  # SELECT FEATURES via EFS (only CoxLasso model)
  efs_coxlasso_nfeats = get_nfeats(efs, lrn_ids = "coxlasso", type = "estimated",
                                   upper_bound = "max_efs", print_params = print_params)

  efs_coxlasso_feats = efs$feature_ranking(
    method = "sav",
    use_weights = TRUE,
    committee_size = efs_coxlasso_nfeats
  )[["feature"]]

  # SELECT FEATURES via EFS (only RSF models)
  rsf_ids = c("rsf_logrank.fselector", "rsf_maxstat.fselector", "aorsf.fselector")
  efs_rsf_nfeats = get_nfeats(efs, lrn_ids = rsf_ids, type = "estimated",
                              upper_bound = "max_efs", print_params = print_params)

  efs_rsf_feats = efs$feature_ranking(
    method = "sav",
    use_weights = TRUE,
    committee_size = efs_rsf_nfeats
  )[["feature"]]

  # SELECT FEATURES via COXLASSO
  task_list = readRDS(file.path("data", dataset_id, "task_list.rds"))
  task = task_list[[omic_id]]$clone() # get single-omic task

  coxlasso = lrn("surv.cv_glmnet", id = "coxlasso", standardize = TRUE, alpha = 1,
                 nfolds = 5, type.measure = "C", s = "lambda.min")
  subsampling = readRDS(file = file.path("data", dataset_id, "subsampling.rds"))
  coxlasso$train(task, row_ids = subsampling$train_set(rsmp_id))
  coxlasso_feats = coxlasso$selected_features()

  # if zero features, take the next lambda
  if (length(coxlasso_feats) == 0) {
    cat(sprintf("[WARNING]: CoxLasso results in 0 selected features, Dataset: %s, Omic: %s, Subsampling Iter: %i\n",
        dataset_id, omic_id, rsmp_id))
    for (lambda in coxlasso$model$model$lambda) {
      coxlasso_feats = coxlasso$selected_features(lambda = lambda)
      if (length(coxlasso_feats) > 0) break
    }
  }

  # Return result as a tibble
  tibble(
    dataset_id = dataset_id,
    omic_id = omic_id,
    rsmp_id = rsmp_id,
    # efs: all models (9)
    efs_all_feats = list(efs_feats),
    efs_all_nfeats = efs_nfeats,
    # efs: just coxlasso (1)
    efs_coxlasso_feats = list(efs_coxlasso_feats),
    efs_coxlasso_nfeats = efs_coxlasso_nfeats,
    # efs: RSF (3)
    efs_rsf_feats = list(efs_rsf_feats),
    efs_rsf_nfeats = efs_rsf_nfeats,
    # Simple coxlasso for feature selection
    coxlasso_nfeats = length(coxlasso_feats),
    coxlasso_feats = list(coxlasso_feats)
  )
}

# execute function of feature selection
execute_fs = function() {
  row_seq = seq_len(nrow(grid_df))

  # Progress tracking
  p = progressr::progressor(along = row_seq)

  data_list = future_lapply(row_seq, function(i) {
    feature_selection(grid_df[i, ], p)
  }, future.seed = TRUE)

  # Combine results into a single dataframe
  bind_rows(data_list)
}

fs_data = execute_fs()

# Save results
saveRDS(fs_data, file = "bench/fs.rds")
