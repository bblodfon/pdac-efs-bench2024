#' Multi-omics (MM) benchmark on PDAC datasets
#'
#' This script runs a multi-omics benchmark on the two PDAC datasets:
#'
#' 1) **Feature Selection:**
#'    - Feature selection per omic and per subsampling iteration has already
#'    been performed in `run_fs.R`.
#'    - The results are loaded here.
#'
#' 2) **Late Integration:**
#'    - For each subsampling iteration, per-omic-selected features are combined
#'    to form multi-omics datasets.
#'    - These datasets are then used to train a model on the training set and
#'    predict on the test set.
#'    - Several model options are available (Coxlasso, RSF, Cox).
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
  library(mlr3pipelines)
  library(mlr3misc)
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
# define/get the different feature selection methods
#fs_method_ids = c("efs_all_feats", "coxlasso_feats")
fs_method_ids = colnames(fs)[endsWith(colnames(fs), "_feats")]

# Define the models we will use for benchmarking and the data these will be applied on
# This is after feature selection per omic is done and late-integration of features
# is performed
model_data_configs = c(
  # Model: CoxLasso, Data: ALL => Clinical + OMICS
  "coxlasso-all",
  # Model: RSF, Data: ALL => Clinical + OMICS
  "rsf-all",
  # Model: CoxLasso, Data: Clinical + GEX only (Reference for Wissel2023/TCGA data)
  "coxlasso-clinical+gex",
  # Model: Cox, Data: Clinical (Reference, for both datasets)
  "cox-clinical",
  # Model: RSF, Data: Clinical (Reference, for both datasets)
  "rsf-clinical",
  # Model: RSF, Data: Clinical + GEX only (Reference for Wissel2023/TCGA data)
  "rsf-clinical+gex"

  # TO TEST
  # Model: CoxLasso, Data: Clinical (mandatory covariates) + MM dataset
  # Model: CoxLasso|RSF, Data: Clinical + MM (remove Mutation + CNV)
)

# Define datasets
dataset_ids = c("wissel2023", "osipov2024")

# Construct a parameter grid
grid_list = lapply(dataset_ids, function(dataset_id) {
  dataset_path = file.path("data", dataset_id)
  assert_directory(dataset_path)

  # Load subsampling iterations
  subsampling = readRDS(file.path(dataset_path, "subsampling.rds"))
  rsmp_ids = seq(subsampling$iters)

  expand.grid(
    dataset_id = dataset_id,
    fs_method_id = fs_method_ids,
    model_data_config = model_data_configs,
    rsmp_id = rsmp_ids,
    stringsAsFactors = FALSE
  )
})

grid_df = do.call(rbind, grid_list)

# remove un-needed configurations
grid_df = grid_df |>
  # Osipov dataset doesn't have GEX data in our benchmark
  filter(!(endsWith(model_data_config, "+gex") & dataset_id == "osipov2024")) |>
  # {cox|rsf, clinical} and {coxlasso|rsf, clinical+gex} don't depend on the fs method
  group_by(rsmp_id, model_data_config) |>
  mutate(fs_method_id = ifelse(
    model_data_config == "cox-clinical" | model_data_config == "coxlasso-clinical+gex" |
      model_data_config == "rsf-clinical" | model_data_config == "rsf-clinical+gex",
    NA,
    fs_method_id
  )) |>
  distinct(dataset_id, fs_method_id, model_data_config, rsmp_id, .keep_all = TRUE) |>
  ungroup()

# Parallelized function for multi-omics benchmark
mm_bench = function(params, p) {
  set.seed(42)
  dataset_id = params$dataset_id
  fs_method_id = params$fs_method_id
  model_data_config = params$model_data_config
  config = strsplit(model_data_config, split = "-")[[1L]]
  model = config[1]
  data = config[2]
  rsmp_id = params$rsmp_id

  #' Notify progress via `p = progressr::progressor()`
  p(sprintf("Dataset: %s, FS-method: %s, Model: %s, Data: %s, Subsampling Iter: %i",
            dataset_id, fs_method_id, model, data, rsmp_id))

  # get clinical + ALL omics tasks for this dataset
  task_list = readRDS(file.path("data", dataset_id, "task_list.rds"))

  # CREATE TASK
  if (data == "clinical") {
    # just take the clinical data, no standardization
    task = task_list$clinical
  } else if (data == "clinical+gex") {
    use_fs = FALSE
    gex_features = if (use_fs) {
      # CoxLasso selected features for GEX data
      result = fs |>
        filter(dataset_id == !!dataset_id,
               omic_id == "gex",
               rsmp_id == !!rsmp_id
        ) |> pull(coxlasso_feats)
      result[[1L]]
    } else {
      # no feature selection for GEX data
      task_list[["gex"]]$feature_names
    }

    gex_data = task_list[["gex"]]$data(cols = gex_features)
    clinical_data = task_list$clinical$data() # time, status included

    # combine clinical + gex data
    all_data = cbind(clinical_data, gex_data)

    # make clinical + gex task
    task = as_task_surv(x = all_data, time = "time", event = "status")

    # standardize data (mean = 0, sd = 1)
    pos = po("scale")
    task = pos$train(list(task))[[1L]]
  } else {
    # combine all omics to a combined multi-omics dataset
    all_data = map_dtc(names(task_list), function(omic_id) {
      if (omic_id == "clinical") {
        # no fs for clinical features
        return(task_list[[omic_id]]$data())
      }

      # fs for omic features
      result = fs |>
        filter(dataset_id == !!dataset_id,
               omic_id == !!omic_id,
               rsmp_id == !!rsmp_id
        ) |> pull(fs_method_id)
      selected_features = result[[1L]]
      task_list[[omic_id]]$data(cols = selected_features)
    })

    # make multi-omics task
    task = as_task_surv(x = all_data, time = "time", event = "status")

    # standardize data (mean = 0, sd = 1)
    pos = po("scale")
    task = pos$train(list(task))[[1L]]
  }

  # GET RESAMPLING: train/test split
  subsampling = readRDS(file = file.path("data", dataset_id, "subsampling.rds"))
  train_set = subsampling$train_set(rsmp_id)
  test_set = subsampling$test_set(rsmp_id)

  # CREATE MODEL
  if (model == "coxlasso") {
    # CoxLasso
    #' `standardize = FALSE` as data (train and test set) is already standardized
    # same config as in Wissel et al. 2023
    learner = lrn("surv.cv_glmnet", id = model, standardize = FALSE, alpha = 1,
                  nfolds = 5, type.measure = "deviance", grouped = TRUE, s = "lambda.min")
  } else if (model == "rsf") {
    # Random Survival Forest
    learner = lrn("surv.ranger", id = model, importance = "none",
                  num.trees = 2000, splitrule = "logrank")
  } else if (model == "cox") {
    # Simple Cox PH
    learner = lrn("surv.coxph")
  } else {
    stopf("Model %s not implemented in this benchmark", model)
  }

  # train model on train set
  learner$train(task, row_ids = train_set)

  # predict with model on test set
  p = learner$predict(task, row_ids = test_set)

  # measure performance via C-index
  harrell_c = p$score(msr("surv.cindex"))
  uno_c = p$score(msr("surv.cindex", weight_meth = "G2"),
                  task = task, train_set = train_set)
  # investigate IBS at different time points and integrated
  ibrier = p$score(msr("surv.graf", times = c(6, 12, 24), ERV = TRUE), # IBS(3 time points)
                   task = task, train_set = train_set)
  brier_at6 = p$score(msr("surv.graf", times = 6, integrated = FALSE, ERV = TRUE), # BS(t = 6 months)
                   task = task, train_set = train_set)
  brier_tmax_12 = p$score(msr("surv.graf", t_max = 12, ERV = TRUE), # IBS(t_max = 12 months)
                   task = task, train_set = train_set)
  brier_tmax_24 = p$score(msr("surv.graf", t_max = 24, ERV = TRUE), # IBS(t_max = 12 months)
                    task = task, train_set = train_set)

  # Return result as a tibble
  tibble(
    dataset_id = dataset_id,
    fs_method_id = fs_method_id,
    rsmp_id = rsmp_id,
    model_data_config = model_data_config,
    # include the train task here?
    # task = task$filter(rows = train_set),
    task_nfeats = task$n_features, # how many multi-omics features were used
    task_feats = list(task$feature_names),
    harrell_c = harrell_c,
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

result = execute_bench()

# Save results
saveRDS(result, file = "bench/result.rds")
