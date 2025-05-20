# convenience function to create an xgboost learner
create_xgb_lrn = function(id, params, aft_loss = NULL) {
  checkmate::assert_list(params, types = "numeric", len = 4, any.missing = FALSE, null.ok = FALSE)
  checkmate::check_subset(
    x = names(params),
    choices = c("eta", "max_depth", "rounds", "nrounds_early"),
    empty.ok = FALSE
  )

  eta = params$eta
  max_depth = params$max_depth
  nrounds = params$nrounds
  nrounds_early = params$nrounds_early

  if (is.null(aft_loss)) {
    lrn("surv.xgboost.cox", id = id,
        booster = "gbtree", tree_method = "hist",
        eta = eta, max_depth = max_depth, nrounds = nrounds,
        early_stopping_rounds = nrounds_early, validate = "test")
  } else {
    lrn("surv.xgboost.aft", id = id,
        booster = "gbtree", tree_method = "hist",
        eta = eta, max_depth = max_depth, nrounds = nrounds,
        early_stopping_rounds = nrounds_early, validate = "test",
        aft_loss_distribution_scale = 1, # sigma = 1
        aft_loss_distribution = aft_loss)
  }
}

# convenience function to create a glmboost `AutoTuner` learner
create_glmb_at = function(id = "glmb_cox", family = "coxph", params) {
  checkmate::assert_list(params, len = 4, any.missing = FALSE, null.ok = FALSE)
  checkmate::check_subset(
    x = names(params),
    choices = c("search_space", "inner_rsmp", "measure", "evals")
  )

  search_space = params$search_space
  inner_rsmp = params$inner_rsmp
  measure = params$measure
  evals = params$evals

  learner = lrn("surv.glmboost", family = family, center = TRUE)

  suppressWarnings({
    at = auto_tuner(
      learner = learner,
      search_space = search_space,
      resampling = inner_rsmp,
      measure = measure,
      tuner = tnr("random_search"),
      term_evals = evals,
      store_tuning_instance = FALSE, # archive
      store_benchmark_result = FALSE, # benchmark result of inner resamplings
      store_models = FALSE,
      callbacks = clbk("mlr3tuning.one_se_rule")
    )
    at$id = id
  })

  at
}

# hack: remove rows which had 0 features selected (due to whatever reason)
# from an `EnsembleFSResult` result (mostly applies to the embedded efs methods)
rm_zero_feat = function(efs) {
  n_zeros = sum(efs$.__enclos_env__$private$.result$n_features == 0)
  if (n_zeros > 0) {
    cat(n_zeros, "model(s) on some subsamples selected 0 features! These are effectively removed\n")
    efs$.__enclos_env__$private$.result = efs$.__enclos_env__$private$.result[n_features > 0]
  }
}

# hack: put as measure the C-index = 1 - OOB_ERROR for RSFs
oob_to_cindex_convert = function(efs) {
  efs$.__enclos_env__$private$.result[, surv.cindex_inner := 1 - oob_error_inner]
  efs$.__enclos_env__$private$.result$oob_error_inner = NULL
  efs$.__enclos_env__$private$.inner_measure = msr("surv.cindex")
}

#' @title Get Number of Features from Pareto Front
#'
#' @description
#' Identifies the optimal number of features using the knee point method
#' on the Pareto front (either empirical or estimated).
#'
#' @param efs `EnsembleFSResult` - The ensemble feature selection result object.
#' @param lrn_ids `character()` - Learner IDs to include (default: all learners).
#' @param type `character(1)` - Type of Pareto front to use: `"empirical"` or
#' `"estimated"` (default).
#' @param upper_bound `character(1)` - Defines upper bound for estimated Pareto front: `"max_efs"` (default, max number of features in the efs object) or `"max_pf"` (max
#' number of features in the empirical PF).
#' This number influences the knee point estimation as larger upper bounds result
#' in larger returned number of features.
#' @param print_params `list()` - Named list with dataset identifiers
#' (`dataset_id`, `omic_id`, `rsmp_id`).
#'
#' @return `integer` - The (estimated) number of features at the knee point.
get_nfeats = function(efs, lrn_ids = NULL, type = "estimated", upper_bound = "max_efs",
                      print_params = NULL) {
  assert_choice(type, c("empirical", "estimated"))
  assert_choice(upper_bound, c("max_efs", "max_pf"))

  efs_copy = efs$clone()

  # Filter learners if specified
  if (!is.null(lrn_ids) && length(lrn_ids) > 0) {
    efs_copy$.__enclos_env__$private$.result =
      efs_copy$.__enclos_env__$private$.result[learner_id %in% lrn_ids]
  }

  # Compute Pareto front
  pf = efs_copy$pareto_front()
  pf_nfeats = pf$n_features

  # Handle case where all Pareto front points have the same feature count
  if (length(unique(pf_nfeats)) == 1) {
    if (!is.null(print_params)) {
      cat(sprintf("[WARNING]: All Pareto front points have %i features. Dataset: %s, Omic: %s, Iter: %i\n",
                  pf_nfeats[1L], print_params$dataset_id, print_params$omic_id, print_params$rsmp_id))
    }
    return(pf_nfeats[1L])
  }

  # Optional diagnostic check for upper limit of estimated Pareto front
  if (FALSE) {
    nfeat_candidates = efs_copy$result |>
      select(learner_id, n_features) |>
      group_by(learner_id) |>
      summarize(max_nfeats = max(n_features)) |>
      arrange(desc(max_nfeats)) |>
      pull()
    pf_estimated = efs_copy$pareto_front(type = "estimated", max_nfeatures = 100)
    pf0001 = pf_estimated[c(TRUE, diff(surv.cindex) > 0.0001)]
    pf001  = pf_estimated[c(TRUE, diff(surv.cindex) > 0.001)]
    cat(sprintf("[CHECK]: %i,%i,%i,%i,%i\n",
                max(efs_copy$result$n_features), nfeat_candidates[2], max(pf_nfeats),
                max(pf0001$n_features), max(pf001$n_features)))
  }

  # Determine upper bound for estimated Pareto front
  if (type == "estimated") {
    max_nfeats = ifelse(upper_bound == "max_efs",
                        max(efs_copy$result$n_features),
                        max(pf_nfeats))

    # Return knee point for estimated Pareto front
    n_feats = efs_copy$knee_points(type = "estimated", max_nfeatures = max_nfeats)$n_features
  } else {
    n_feats = efs_copy$knee_points()$n_features
  }

  n_feats
}
