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
# how many different feature selection methods?
fs_methods = c("efs_all_feats", "coxlasso_nfeats")
stopifnot(all(fs_methods %in% names(fs)))

# Define datasets
dataset_ids = c("wissel2023", "osipov2024")

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

  expand.grid(
    dataset_id = dataset_id,
    fs_method = fs_methods,
    rsmp_id = rsmp_ids,
    stringsAsFactors = FALSE
  )
})

grid_df = do.call(rbind, grid_list)

# Parallelized function for multi-omics benchmark
mm_bench = function(params, p) {
  set.seed(42)
  dataset_id = params$dataset_id
  fs_method_id = params$fs_method
  rsmp_id = params$rsmp_id
  #print_params = list(dataset_id = dataset_id, fs_method_id = fs_method, rsmp_id = rsmp_id)

  #' Notify progress via `p = progressr::progressor()`
  p(sprintf("Dataset: %s, FS-method: %s, Subsampling Iter: %i",
            dataset_id, fs_method, rsmp_id))

  # get all omics tasks for this dataset
  task_list = readRDS(file.path("data", dataset_id, "task_list.rds"))

  # subset each omic task to the selected features
  for (omic_id in names(task_list)) {
    if (omic_id != "clinical") {
      result = fs |>
        filter(dataset_id == !!dataset_id,
               omic_id == !!omic_id,
               rsmp_id == !!rsmp_id
        ) |> pull(fs_method_id)
      selected_features = result[[1L]]
      task_list[[omic_id]]$select(cols = selected_features)
    }
  }

  # combine omics to a combined multi-omics dataset
  all_data = mlr3misc::map_dtc(names(task_list), function(omic_id) {
    task_list[[omic_id]]$data()
  })

  # remove auxiliary target columns
  all_data = all_data[, !grepl("^time\\.|^status\\.", names(all_data)), with = FALSE]

  # make multi-omics task
  mm_task = as_task_surv(x = all_data, time = "time", event = "status")

  # standardize data (mean = 0, sd = 1)
  pos = po("scale")
  mm_task = pos$train(list(mm_task))[[1L]]

  # train/test split
  subsampling = readRDS(file = file.path("data", dataset_id, "subsampling.rds"))
  train_set = subsampling$train_set(rsmp_id)
  test_set = subsampling$test_set(rsmp_id)

  # train CoxLasso model
  #' `standardize = FALSE` as data (train and test set) is already standardized
  coxlasso = lrn("surv.cv_glmnet", id = "coxlasso", standardize = FALSE, alpha = 1,
                 nfolds = 5, type.measure = "C", s = "lambda.min")
  coxlasso$train(mm_task, row_ids = train_set)
  coxlasso_feats = coxlasso$selected_features()

  # predict with CoxLasso model
  p = coxlasso$predict(mm_task, row_ids = test_set)

  # measure performance via C-index
  harrel_c = p$score(msr("surv.cindex"))
  uno_c = p$score(msr("surv.cindex", weight_meth = "G2"),
                  task = mm_task, train_set = train_set)
  dcalib = p$score(msr("surv.dcalib", truncate = 16)) # see doc for truncate value
  ibrier = p$score(msr("surv.graf", times = c(6, 12, 24), ERV = TRUE), # IBS
                  task = mm_task, train_set = train_set)

  # Return result as a tibble
  tibble(
    dataset_id = dataset_id,
    fs_method_id = fs_method_id,
    rsmp_id = rsmp_id,
    model = "LI-coxlasso", # for now only this
    harrel_c = harrel_c,
    uno_c = uno_c,
    dcalib = dcalib,
    ibrier = ibrier
  )
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
