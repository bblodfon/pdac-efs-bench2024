#' This script runs a LASSO-based embedded feature selection, for a specific
#' dataset, omic and resampling id (train set)
#'
#' Execute: `Rscript bench/run_lasso_fs.R <dataset_id> <omic_id> <rsmp_id>` (from project root)
#' e.g. `Rscript bench/run_lasso_fs.R wissel2023 gex 42`

# use `renv`
#renv::load()
#renv::deactivate()

# LIBRARIES ----
suppressPackageStartupMessages({
  library(mlr3)
  library(mlr3extralearners)
  library(mlr3proba)
  library(glmnet)
  library(checkmate)
})

# CMD args ----
args = commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript bench/run_lasso_fs.R <dataset_id> <omic_id> <rsmp_id>")
}

dataset_id = args[1]
omic_id = args[2]
rsmp_id = as.integer(args[3])

# perform checks
dataset_path = file.path("data", dataset_id)
assert_directory(dataset_path)

task_list = readRDS(file = file.path(dataset_path, "task_list.rds"))
omic_ids = names(task_list)
assert_subset(omic_id, omic_ids)

subsampling = readRDS(file = file.path(dataset_path, "subsampling.rds"))
assert_number(rsmp_id, lower = 1, upper = subsampling$iters)

# make directory for results if it doesn't already exist
res_path = file.path("bench", "fs", dataset_id, omic_id)
if (!test_directory_exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
}

# TASK ----
task = task_list[[omic_id]]$clone()
# run cox lasso only using the train set of this particular subsampling
task$filter(rows = subsampling$train_set(rsmp_id))
task$set_col_roles(cols = "status", add_to = "stratum") # stratify by status

# CONFIG ----
cfg = config::get(file = "bench/config.yml")

# verbose message printing
print_msg = function(...) {
  if (cfg$verbose > 0) cat(...)
}

# result lasso selected features
lasso_fs_path = file.path(res_path, paste0("lasso_", rsmp_id, ".rds"))

# if results exist, don't run again
if (file.exists(lasso_fs_path) && !cfg$overwrite) {
  cat("Exiting... results already exist and we don't overwrite\n")
  quit()
}

# print some basic info about this run
print_msg("Dataset:", dataset_id,
          "\nOmic:", omic_id,
          "\nResampling Iteration:", rsmp_id,
          "\n")

measure = msr("surv.cindex")

# Run CoxLasso model
# Does internal CV tuning of lambda
#' `s` => which `lambda` to use, `lambda.min` => more features, higher C-index
#' `lambda.1se` => less features, smaller C-index
coxlasso = lrn("surv.cv_glmnet", id = "coxlasso", standardize = TRUE, alpha = 1,
               nfolds = cfg$learner_params$folds, type.measure = "C", s = "lambda.min")
set.seed(42)
coxlasso$train(task)

saveRDS(coxlasso, file = lasso_fs_path)
# coxlasso$selected_features()
# unname(coxlasso$timings["train"])
