#' Multi-omics (MM) benchmark on PDAC datasets
#'
#' This script runs a multi-omics benchmark on the available PDAC datasets:
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
#'    - Several model options are available (Coxlasso, RSF, CoxPH, BlockForest)
#'    as well as which data types to integrate.
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
#' `a`: Integration model for late fusion (or baseline) => {`coxlasso`, `rsf`, `blockforest`}
#' `b`: Data type(s) => {`all`, `clinical`, `gex`, `clinical+gex`}.
#'      Note: `all` means clinical + ALL omics available
#' `c`: FS method = {`none`, `coxlasso_all_feats`, etc.}.
#'      Note" If absent, all available FS methods (`fs_method_ids`) are considered.
model_data_fs_configs = c(
  # Integration Model: CoxLasso, Data: ALL => Clinical + OMICS
  "coxlasso-all",
  # Integration Model: CoxLasso, Data: ALL => Clinical + OMICS, no FS method
  "coxlasso-all-none",
  # Integration Model: RSF, Data: ALL => Clinical + OMICS
  "rsf-all",
  # Integration Model: RSF, Data: ALL => Clinical + OMICS, no FS method
  "rsf-all-none",
  # (Baseline) Model: CoxPH, Data: Clinical (Reference, for both datasets)
  "cox-clinical",
  # (Baseline) Model: RSF, Data: Clinical (Reference, for both datasets)
  "rsf-clinical",
  # Integration Model: RSF, Data: Clinical + GEX, FS method for GEX: hEFS (9 models)
  "rsf-clinical+gex-efs_all_feats",
  # Integration Model: RSF, Data: GEX, FS method for GEX: hEFS (9 models)
  "rsf-gex-efs_all_feats",
  # Integration Model: CoxLasso, Data: Clinical + GEX, FS method for GEX: CoxLasso
  #"coxlasso-clinical+gex-coxlasso_feats",
  # Integration Model: BlockForest, Data: ALL
  "blockforest-all",
  # Integration Model: BlockForest, Data: ALL, no FS method
  "blockforest-all-none",
  # Integration Model: BlockForest, Data: Clinical + GEX, FS method for GEX: CoxLasso
  "blockforest-clinical+gex-efs_all_feats"
)

# Define datasets
dataset_ids = c("wissel2023", "osipov2024", "cao2021")

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
    model_data_fs_config = model_data_fs_configs,
    rsmp_id = rsmp_ids,
    stringsAsFactors = FALSE
  )
})

grid_df = do.call(rbind, grid_list)

# filter configurations
grid_df_filtered = grid_df |>
  # Step 1: osipov2024 dataset doesn't have GEX data
  filter(!(dataset_id == "osipov2024" & str_detect(model_data_fs_config, "gex"))) |>
  # Step 2: for configs with a 3rd part (pattern: a-b-c), ensure 3rd part == fs_method_id
  filter({
    third_part = str_match(model_data_fs_config, "^[^-]+-[^-]+-(.+)$")[,2]
    is.na(third_part) | third_part == fs_method_id | third_part == "none"
  })

# Clinical-only and "none" fs config DO NOT do any feature selection
nofs_df = grid_df_filtered |>
  filter(str_ends(model_data_fs_config, "clinical") |
         str_ends(model_data_fs_config, "none")) |>
  distinct(dataset_id, model_data_fs_config, rsmp_id) |>
  mutate(fs_method_id = NA)

# Extract the rest
rest = grid_df_filtered |>
  filter(!str_ends(model_data_fs_config, "clinical") &
         !str_ends(model_data_fs_config, "none"))

if (nrow(nofs_df) > 0) {
  grid_df = bind_rows(nofs_df, rest)
} else {
  grid_df = rest
}

# how many configs to run?
cat(sprintf("--- %i configs to run in total ----", nrow(grid_df)))

# helper function to select features
select_features = function(dataset_id, fs_method_id, rsmp_id, omic_id) {
  if (fs_method_id == "none" || is.na(fs_method_id)) return(NULL)

  result = fs |>
    filter(dataset_id == !!dataset_id,
           omic_id == !!omic_id,
           rsmp_id == !!rsmp_id) |>
    pull(fs_method_id)

  if (length(result) == 0) return(NULL)

  result[[1L]]
}

# Parallelized function for multi-omics benchmark
mm_bench = function(params, p) {
  set.seed(42)

  dataset_id = params$dataset_id
  fs_method_id = params$fs_method_id %||% "none" # none or NA means: no FS
  model_data_fs_config = params$model_data_fs_config
  rsmp_id = params$rsmp_id

  # parse config
  parts = strsplit(model_data_fs_config, split = "-")[[1L]]
  model = parts[1]
  data = parts[2]

  #' Notify progress via `p = progressr::progressor()`
  p(sprintf("Dataset: %s, FS-method: %s, Model: %s, Data: %s, SubSmp iter: %i",
            dataset_id, fs_method_id, model, data, rsmp_id))

  # Load tasks and sub-sampling (train/test split)
  task_list = readRDS(file.path("data", dataset_id, "task_list.rds"))
  subsampling = readRDS(file.path("data", dataset_id, "subsampling.rds"))
  train_set = subsampling$train_set(rsmp_id)
  test_set = subsampling$test_set(rsmp_id)

  # === Create Task ===
  task = switch(data,
    "clinical" = {
      task_list$clinical
    },
    "gex" = {
      gex_features = select_features(dataset_id, fs_method_id, rsmp_id, omic_id = "gex")
      gex_data = task_list$gex$data(
        cols = c("time", "status", gex_features %||% task_list$gex$feature_names)
      )
      task = as_task_surv(x = gex_data, time = "time", event = "status")

      # standardize data (mean = 0, sd = 1)
      pos = po("scale")
      pos$train(list(task))[[1L]]
    },
    "clinical+gex" = {
      gex_features = select_features(dataset_id, fs_method_id, rsmp_id, omic_id = "gex")
      gex_data = task_list$gex$data(
        cols = gex_features %||% task_list$gex$feature_names
      )
      clinical_data = task_list$clinical$data()
      all_data = cbind(clinical_data, gex_data)
      task = as_task_surv(x = all_data, time = "time", event = "status")

      # standardize data (mean = 0, sd = 1)
      pos = po("scale")
      pos$train(list(task))[[1L]]
    },
    "all" = {
      # combine all omics to a combined multi-omics dataset
      omic_data = map_dtc(names(task_list), function(omic_id) {
        # investigate: remove mutation omic (default FALSE: DON'T DO THIS)
        if (FALSE && dataset_id == "wissel2023" && omic_id == "mutation") {
          return(data.table())
        }

        features = select_features(dataset_id, fs_method_id, rsmp_id, omic_id)

        if (omic_id == "clinical") {
          task_list[[omic_id]]$data()
        } else {
          task_list[[omic_id]]$data(cols = features %||% task_list[[omic_id]]$feature_names)
        }
      })
      task = as_task_surv(x = omic_data, time = "time", event = "status")

      # standardize data (mean = 0, sd = 1)
      pos = po("scale")
      pos$train(list(task))[[1L]]
    },
    stopf("Unsupported data config: %s", data)
  )

  # === Learner ===
  learner = switch(model,
     # CoxLasso
     #' `standardize = FALSE` as data (train and test set) is already standardized
     #' same config as in Wissel et al. (2023)
     "coxlasso" = lrn("surv.cv_glmnet", id = model, standardize = FALSE,
                      alpha = 1, nfolds = 5, type.measure = "deviance",
                      grouped = TRUE, s = "lambda.min"),
     # Random Survival Forest
     "rsf" = lrn("surv.ranger", id = model, importance = "none",
                 num.trees = 2000, splitrule = "logrank", num.threads = 3),
     # Cox Proportional Hazards
     "cox" = lrn("surv.coxph"),
     # Block Forest
     "blockforest" = {
       omic_prefixes = setdiff(names(task_list), "clinical")
       blocks = get_block_indices(feature_names = task$feature_names,
                                  omic_prefixes = omic_prefixes)
       lrn("surv.blockforest", blocks = blocks, splitrule = "logrank",
           num.trees = 2000, nsets = 300, num.trees.pre = 1000, num.threads = 8)
     },
     stopf("Model %s not implemented", model)
  )

  # === Fit and Predict ===
  learner$train(task, row_ids = train_set)
  pred = learner$predict(task, row_ids = test_set)

  # === Score ===
  tibble(
    dataset_id = dataset_id,
    fs_method_id = fs_method_id,
    rsmp_id = rsmp_id,
    model_data_config = model_data_fs_config,
    task_nfeats = task$n_features, # how many multi-omics features were used
    task_feats = list(task$feature_names),
    harrell_c = pred$score(msr("surv.cindex")),
    uno_c = pred$score(msr("surv.cindex", weight_meth = "G2"),
                       task = task, train_set = train_set),
    brier_tmax24 = pred$score(msr("surv.graf", t_max = 24, ERV = TRUE),
                              task = task, train_set = train_set)
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
