gen_sizes = function(n_features, size, shape2) {
  bias = stats::dbeta(seq(1, 0.001, length.out = n_features), shape1 = 0.5, shape2 = shape2)
  # if (length(bias) != n_features) {
  #   stop("Try different `shape2` or `size` parameters...")
  # }
  bb = bias/sum(bias) # normalize to probabilities of number selection n:1

  bb[1] = 1 # always include the largest number of features (n)
  bb = head(bb, -1) # do not include just 1 feature
  # randomly sample from n:2 vector with (biased/skewed) probabilities `bb`
  sort(sample(x = n_features:2, size = size, replace = FALSE, prob = bb), decreasing = TRUE)
}

load_callback_subset_sizes = function() {
  gen_sizes = function(n_features, size, shape2) {
    bias = stats::dbeta(seq(1, 0.001, length.out = n_features), shape1 = 0.5, shape2 = shape2)
    # if (length(bias) != n_features) {
    #   stop("Try different `shape2` or `size` parameters...")
    # }
    bb = bias/sum(bias) # normalize to probabilities of number selection n:1

    bb[1] = 1 # always include the largest number of features (n)
    bb = head(bb, -1) # do not include just 1 feature
    # randomly sample from n:2 vector with (biased/skewed) probabilities `bb`
    sort(sample(x = n_features:2, size = size, replace = FALSE, prob = bb), decreasing = TRUE)
  }

  callback_batch_fselect("mlr3fselect.rfe_subset_sizes",
     label = "Random subset sizes",
     man = "mlr3fselect::mlr3fselect.rfe_subset_sizes",

     on_optimization_begin = function(callback, context) {
       if (context$optimizer$id != "rfe") {
         stop("`mlr3fselect.rfe_subset_sizes` callback works only with 'rfe' fselector")
       }

       # get task features number
       n_features = context$instance$objective$task$n_features
       # how many subset sizes to generate
       size = as.integer(checkmate::assert_number(callback$state$size, lower = 1, na.ok = FALSE, null.ok = FALSE))
       checkmate::assert_true(size < n_features)

       # beta argument for Beta distribution (the larger than 1 => more skewed towards lower features)
       shape2 = checkmate::assert_number(callback$state$shape2, na.ok = FALSE, null.ok = FALSE)

       # generate the RFE subset sizes
       sizes = gen_sizes(n_features, size, shape2)

       # change `subset_sizes` from fselector `rfe`
       context$optimizer$param_set$set_values(
         .values = list(subset_sizes = sizes, feature_fraction = NULL, feature_number = NULL)
       )
     }
  )
}
x = utils::getFromNamespace("mlr_callbacks", ns = "mlr3misc")
x$add("mlr3fselect.rfe_subset_sizes", load_callback_subset_sizes())
