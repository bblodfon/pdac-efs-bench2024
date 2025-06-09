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
  library(stringr)
  library(future.apply)
  library(progressr)
})
source("bench/blockForest.R")

# Set parallel execution
plan("multicore", workers = 50)

# Enable progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers("progress")

# load feature selection results
fs = readRDS(file = "bench/fs.rds")
# define/get the different feature selection methods
fs_method_ids = colnames(fs)[endsWith(colnames(fs), "_feats")]

#' Pattern: a-b-c, where:
#' `a`: integration model for late fusion (or baseline)
#' `b`: data type(s)
#' `c`: FS method (if absent all FS methods are considered)
model_data_configs = c(
  # Integration Model: CoxLasso, Data: ALL => Clinical + OMICS
  "coxlasso-all",
  # Integration Model: RSF, Data: ALL => Clinical + OMICS
  "rsf-all",
  # (Baseline) Model: CoxPH, Data: Clinical (Reference, for both datasets)
  "cox-clinical",
  # (Baseline) Model: RSF, Data: Clinical (Reference, for both datasets)
  "rsf-clinical",
  # Integration Model: RSF, Data: Clinical + GEX, FS method for GEX: hEFS (9 models)
  "rsf-clinical+gex-efs_all_feats",
  # Integration Model: RSF, Data: GEX, FS method for GEX: hEFS (9 models)
  "rsf-gex-efs_all_feats",
  # Integration Model: CoxLasso, Data: Clinical + GEX, FS method for GEX: CoxLasso
  "coxlasso-clinical+gex-coxlasso_feats",
  # Integration Model: BlockForest, Data: ALL
  "blockforest-all",
  # Integration Model: BlockForest, Data: Clinical + GEX, FS method for GEX: CoxLasso
  "blockforest-clinical+gex-efs_all_feats"
)

# Define datasets
dataset_ids = c("wissel2023", "osipov2024")

# Construct a parameter grid
grid_list = lapply(dataset_ids, function(dataset_id) {
  dataset_path = file.path("data", dataset_id)
  checkmate::assert_directory(dataset_path)

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

# filter configurations
grid_df_filtered = grid_df |>
  # Step 1: osipov2024 dataset doesn't have GEX data
  filter(!(dataset_id == "osipov2024" & str_detect(model_data_config, "gex"))) |>
  # Step 2: for configs with a 3rd part (pattern: a-b-c), ensure 3rd part == fs_method_id
  filter({
    third_part = str_match(model_data_config, "^[^-]+-[^-]+-(.+)$")[,2]
    is.na(third_part) | third_part == fs_method_id
  })

# Clinical-only config DO NOT do any feature selection
clinical_only = grid_df_filtered |>
  filter(model_data_config %in% c("cox-clinical", "rsf-clinical")) |>
  distinct(dataset_id, model_data_config, rsmp_id) |>
  mutate(fs_method_id = NA)

# Extract the rest
rest = grid_df_filtered |>
  filter(!model_data_config %in% c("cox-clinical", "rsf-clinical"))

if (nrow(clinical_only) > 0) {
  grid_df = bind_rows(clinical_only, rest)
} else {
  grid_df = rest
}

# how many configs to run?
nrow(grid_df)

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
  } else if (data == "gex") {
    # Feature selection for GEX data
    result = fs |>
      filter(dataset_id == !!dataset_id,
             omic_id == "gex",
             rsmp_id == !!rsmp_id
      ) |> pull(fs_method_id)
    gex_features = result[[1L]]
    all_data = task_list[["gex"]]$data(cols = c("time", "status", gex_features))
    task = as_task_surv(x = all_data, time = "time", event = "status")

    # standardize data (mean = 0, sd = 1)
    pos = po("scale")
    task = pos$train(list(task))[[1L]]
  } else if (data == "clinical+gex") {
    # Feature selection for GEX data
    result = fs |>
      filter(dataset_id == !!dataset_id,
             omic_id == "gex",
             rsmp_id == !!rsmp_id
      ) |> pull(fs_method_id)
    gex_features = result[[1L]]
    # gex_features = task_list[["gex"]]$feature_names # if no FS for GEX data

    gex_data = task_list[["gex"]]$data(cols = gex_features)
    clinical_data = task_list$clinical$data() # time, status included

    # combine clinical + gex data
    all_data = cbind(clinical_data, gex_data)

    # make clinical + gex task
    task = as_task_surv(x = all_data, time = "time", event = "status")

    # standardize data (mean = 0, sd = 1)
    pos = po("scale")
    task = pos$train(list(task))[[1L]]
  } else { # `all`
    # combine all omics to a combined multi-omics dataset
    all_data = map_dtc(names(task_list), function(omic_id) {
      # investigate: remove mutation omic (default FALSE: DON'T DO THIS)
      if (FALSE && dataset_id == "wissel2023" && omic_id == "mutation") {
        return(data.table())
      }

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
  } else if (model == "blockforest") {
    omic_prefixes = setdiff(names(task_list), "clinical")
    blocks = get_block_indices(feature_names = task$feature_names,
                               omic_prefixes = omic_prefixes)
    learner = lrn("surv.blockforest", blocks = blocks, splitrule = "logrank",
                  num.trees = 2000, nsets = 500, num.trees.pre = 500, num.threads = 8)
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
  # investigate IBS (t_max = 24 months)
  brier_tmax24 = p$score(msr("surv.graf", t_max = 24, ERV = TRUE),
                         task = task, train_set = train_set)

  # Return result as a tibble
  tibble(
    dataset_id = dataset_id,
    fs_method_id = fs_method_id,
    rsmp_id = rsmp_id,
    model_data_config = model_data_config,
    # include the train task here if needed
    # task = task$filter(rows = train_set),
    task_nfeats = task$n_features, # how many multi-omics features were used
    task_feats = list(task$feature_names),
    harrell_c = harrell_c,
    uno_c = uno_c,
    brier_tmax24 = brier_tmax24
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
saveRDS(result, file = "bench/result_bf.rds")

