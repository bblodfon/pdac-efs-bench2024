#' Multi-omics benchmark on PDAC datasets
#'
#' This script runs a multi-omics benchmark on the two PDAC datasets:
#' 1) First we perform feature selection per omic and per subsampling iteration
#' (100 in total) - either using CoxLasso or the ensmelbe feature selection method.
#' The latter, being computationally intensive, has been executed in a separate
#' `run_efs.sh` script, so we just load the result object here and decide on the
#' number of features via the estimated Pareto front method).
#' 2) We perform *late integration*: i.e. for each subsampling iteration, we combine
#' the omic-specific features found by CoxLasso and efs, to make two multiomics
#' datasets (one guided by the coxlasso embedded fs and the other by the efs),
#' which we then fit to a CoxLasso model on the train set and predict using the test set.
#'
#' Execute: `Rscript bench/run_mm_bench.R` (from project root)

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

# Set parallel execution
plan("multicore", workers = 15)

# Enable progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers("progress")

# Define datasets
dataset_ids = c("wissel2023", "osipov2024")

cat("FEATURE SELECTION\n")

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

  #' Notify progress via `p = progressr::progressor()`
  p(sprintf("Dataset: %s, #Omic: %s, #Subsampling Iter: %i", dataset_id, omic_id, rsmp_id))

  # Load ensemble feature selection (EFS) results
  file_name = file.path("bench", "efs", dataset_id, omic_id, paste0("efs_", rsmp_id, ".rds"))
  efs = readRDS(file_name)

  # SELECT FEATURES via EFS
  pf_nfeats = efs$pareto_front()$n_features
  if (length(unique(pf_nfeats)) == 1) {
    cat(sprintf("[WARNING]: All Pareto front points have the same number of features (%i),
                Dataset: %s, #Omic: %s, #Subsampling Iter: %i\n", pf_nfeats[1L],
                dataset_id, omic_id, rsmp_id))
    n_feats = pf_nfeats[1L]
  } else {
    # Choose max features from empirical ePF for upper limit of the estimated PF
    # Use 20, if PF doesn't have points with more features than that
    max_nfeats = max(max(pf_nfeats), 20)
    n_feats = efs$knee_points(type = "estimated", max_nfeatures = max_nfeats)$n_features
  }

  efs_feats = efs$feature_ranking(
    method = "sav",
    use_weights = TRUE,
    committee_size = n_feats
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
    cat("[WARNING]: CoxLasso results in zero selected features\n")
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
    efs_n_feats = length(efs_feats),
    efs_feats = list(efs_feats),
    coxlasso_n_feats = length(coxlasso_feats),
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

# LATE INTEGRATION ----
# measure = msr("surv.cindex")
