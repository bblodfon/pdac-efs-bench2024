#' Multi-omics benchmark on PDAC datasets
#'
#' This script runs a multi-omics benchmark on the two PDAC datasets:
#' 1) First the feature selection per omic and per subsampling iteration
#' has been performed in the `run_fs.R` script. We just load the results here.
#' 2) We perform *late integration*: i.e. for each subsampling iteration, we combine
#' the omic-specific features to make multiomics datasets, for which we then
#' fit to a CoxLasso model on the train set and predict using the test set.
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
plan("multicore", workers = 30)

# Enable progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers("progress")

# load feature selection results
fs = readRDS(file = "bench/fs.rds")

# Define datasets
dataset_ids = c("wissel2023", "osipov2024")

# LATE INTEGRATION ----
# Construct a parameter grid (all combinations of dataset_id, all omic_ids, rsmp_id)
grid_list = lapply(dataset_ids, function(dataset_id) {
  dataset_path = file.path("data", dataset_id)
  assert_directory(dataset_path)

  # Load omic task IDs
  omic_ids = readr::read_csv(file.path(dataset_path, "omic_ids.csv"),
                             col_names = FALSE, show_col_types = FALSE)[[1L]]
  omic_ids = c("clinical", omic_ids)

  # Load subsampling iterations
  subsampling = readRDS(file.path(dataset_path, "subsampling.rds"))
  rsmp_ids = seq(subsampling$iters)

  expand.grid(dataset_id = dataset_id, omic_id = list(omic_ids), rsmp_id = rsmp_ids,
              stringsAsFactors = FALSE)
})

grid_df = do.call(rbind, grid_list)

# Parallelized function for multi-omics benchmark
mm_bench = function(params, p) {
  set.seed(42)
  dataset_id = params$dataset_id
  omic_id = params$omic_id
  rsmp_id = params$rsmp_id
  print_params = list(dataset_id = dataset_id, omic_id = omic_id, rsmp_id = rsmp_id)

  #' Notify progress via `p = progressr::progressor()`
  p(sprintf("Dataset: %s, Omic: %s, Subsampling Iter: %i", dataset_id, omic_id, rsmp_id))

  # make multi-omics task
  task_list = readRDS(file.path("data", dataset_id, "task_list.rds"))
}

# execute function of benchmark
execute_bench = function() {
  row_seq = seq_len(nrow(grid_df))

  # Progress tracking
  p = progressr::progressor(along = row_seq)

  data_list = future_lapply(row_seq, function(i) {
    mm_bench(grid_df[i, ], p)
  }, future.seed = TRUE)

  # Combine results into a single dataframe
  bind_rows(data_list)
}

bench_res = execute_bench()

# Save results
saveRDS(bench_res, file = "bench/bench_res.rds")

measure = msr("surv.cindex")
