library(mlr3)
library(mlr3extralearners)
library(blockForest) # 0.2.6 version
library(R6)
library(paradox)

#' Full documentation => `BlockForest::blockfor`.
LearnerSurvBlockForest = R6::R6Class("LearnerSurvBlockForest",
  inherit = mlr3proba::LearnerSurv,
  public = list(
    initialize = function() {
      ps = paradox::ps(
        blocks = p_uty(tags = c("train", "required")),
        block.method = p_fct(c("BlockForest", "RandomBlock", "BLockVarSel", "VarProb", "SplitWeights"), default = "BlockForest", tags = "train"),
        num.trees = p_int(1, 2000, default = 2000, tags = "train"),
        mtry = p_uty(default = NULL, tags = "train"),
        nsets = p_int(1, 300, default = 300, tags = "train"),
        num.trees.pre = p_int(1, 1500, default = 1500, tags = "train"),
        splitrule = p_fct(c("logrank", "extratrees", "C", "maxstat"), default = "extratrees", tags = "train"),
        always.select.block = p_int(0, 1, default = 0, tags = "train"),
        importance = p_fct(c("none", "impurity", "impurity_corrected", "permutation"), tags = "train"),
        num.threads = p_int(1L, default = 1L, tags = c("train", "predict", "threads"))
      )

      ps$values = list(
        block.method = "BlockForest",
        num.trees = 2000,
        mtry = NULL,
        nsets = 300,
        num.trees.pre = 1500,
        splitrule = "extratrees",
        always.select.block = 0,
        num.threads = 1
      )

      super$initialize(
        id = "surv.blockforest",
        param_set = ps,
        predict_types = c("crank", "distr"),
        feature_types = c("logical", "integer", "numeric", "character", "factor", "ordered"),
        properties = c("weights", "importance"),
        packages = c("mlr3extralearners", "blockForest"),
        label = "Block Forests: Random Forests for Blocks of Clinical and Omics Covariate Data"
      )
    },

    #' @description
    #' The importance scores are extracted from the model slot `variable.importance`.
    #' @return Named `numeric()`.
    importance = function() {
      if (is.null(self$model)) {
        stopf("No model stored")
      }
      if (self$model$forest$importance.mode == "none") {
        stopf("No importance stored")
      }

      sort(self$model$forest$variable.importance, decreasing = TRUE)
    }
  ),

  private = list(
    .train = function(task) {
      pv = self$param_set$get_values(tags = "train")

      mlr3misc::invoke(blockForest::blockfor,
        X = task$data(cols = task$feature_names),
        y = task$truth(),
        case.weights = task$weights$weight,
        .args = pv
      )
    },

    .predict = function(task) {
      pv = self$param_set$get_values(tags = "predict")
      newdata = mlr3extralearners:::ordered_features(task, self)
      prediction = mlr3misc::invoke(predict, object = self$model$forest, data = newdata, .args = pv)
      mlr3proba::.surv_return(times = prediction$unique.death.times, surv = prediction$survival)
    }
  )
)

x = utils::getFromNamespace("mlr_learners", ns = "mlr3")
x$add("surv.blockforest", LearnerSurvBlockForest)

# helper function to derive the `blocks` parameter for the learner,
# given a vector of features from clinical + multi-omics task
get_block_indices = function(feature_names, omic_prefixes) {
  block_indices = list()

  for (prefix in omic_prefixes) {
    block_indices[[prefix]] = which(startsWith(feature_names, paste0(prefix, "_")))
  }

  # Clinical = everything not matched to an omic prefix
  matched = unlist(block_indices)
  block_indices[["clinical"]] = setdiff(seq_along(feature_names), matched)

  # Remove empty blocks
  block_indices = Filter(length, block_indices)

  # Sanity check
  total_indices = sum(lengths(block_indices))
  if (total_indices != length(feature_names)) {
    stop(sprintf("Sanity check failed: %d indices vs %d feature names",
                 total_indices, length(feature_names)))
  }

  block_indices
}
