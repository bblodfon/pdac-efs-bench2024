suppressPackageStartupMessages({
  library(mlr3)
  library(mlr3proba)
  library(mlr3fselect)
  library(mlr3tuning)
  library(rpart)
  library(progressr)
  library(tictoc)
})

source("efs/helpers.R")
source("efs/callbacks.R")

# Logging
lgr::get_logger("bbotk")$set_threshold("warn")
lgr::get_logger("mlr3" )$set_threshold("warn")

# progress bars
options(progressr.enable = TRUE)
handlers(on_missing = "ignore", global = TRUE)
handlers("progress")

# parallelization
future::plan("multisession", workers = 10)

# Get RFE iters/subset_sizes
n_features = 2 # RFE => stopping criterion (run up to this number of features)
feature_fraction = 0.8 # RFE => keep this % of features in each RFE iteration
# feature subset sizes for RFE, given we go from n => 2, keeping 80% each time
n = task$n_features
subset_sizes = unique(floor(cumprod(c(n, rep(feature_fraction, log(n_features / n) / log(feature_fraction))))))
num_subsets = 15 # RFE iterations => we don't want to do more than this due to RAM issues
if (length(subset_sizes) > num_subsets) {
  # thin the feature subsets to 15 values in total
  indx = unique(round(seq.int(1, length(subset_sizes), length.out = num_subsets)))
  subset_sizes = subset_sizes[indx]
}
rfe = fs("rfe", subset_sizes = subset_sizes)

rfe = fs("rfe")
ss_clbk = clbk("mlr3fselect.rfe_subset_sizes")
# set these parameters in the callback
ss_clbk$state$size = 15 # how many RFE iters (subset sizes)
ss_clbk$state$shape2 = 15 # beta distr parameter, > 1 => more skewed towards lower number of features

cat("# Wrapper-based efs with SURVIVAL TREE")
tic()
set.seed(42) # reproduce: same subsampling
suppressWarnings({
  efs_tree_no_calbk = ensemble_fselect(
    fselector = rfe,
    task = readRDS(file = "data/wissel2023/task_list.rds")$gex,
    learners = list(lrn("surv.rpart", id = "tree", minsplit = 15)),
    init_resampling = rsmp("subsampling", repeats = 100, ratio = 0.8),
    inner_resampling = rsmp("cv", folds = 5),
    inner_measure = msr("surv.cindex"),
    measure = msr("surv.cindex"),
    terminator = trm("none"),
    callbacks = list(
      tree = list(clbk("mlr3fselect.one_se_rule"), ss_clbk)
    ),
    store_benchmark_result = store_bmr,
    store_models = FALSE
  )
})
toc()

# saveRDS(efs_tree, file = paste0("efs/", task$id, "/efs_tree.rds"))
