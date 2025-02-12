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
