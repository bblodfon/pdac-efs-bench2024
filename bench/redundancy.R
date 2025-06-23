#' **Redundancy** analysis of feature selection methods applied to multi-omic datasets.
#' We do this analysis in a separate script since it takes some time even if parallelized.
#'
#' We calculate the pair-wise Pearson, Spearman and XI correlation (Chatterjee 2021) coefficients
#' for every selected feature set and provide the following summary scores:
#' - Mean absolute correlation between all unique feature pairs (Redundancy Rate or RR).
#' - The Adjusted Redundancy Rate (chance-corrected, i.e. the adjusted score is
#' independent of the number of chosen features so FS methods with different number
#' of selected features can be compared)
#' - Proportion of FDR-significant pairs (Significant Redundancy Proportion or SRP).
#' i.e. "What proportion of pairs are significantly redundant?"
#'
#' Execute: `Rscript bench/redundancy.R` (from project root)

suppressPackageStartupMessages({
  library(XICOR)
  library(dplyr)
  library(tidyr)
  library(progressr)
  library(future.apply)
  library(mlr3misc)
})

# progress bars
options(progressr.enable = TRUE)
handlers(on_missing = "ignore", global = TRUE)
handlers("progress")

#' Helper function to compute redundancy scores in a given `mlr3` task
#' `alpha`: FDR threshold for SRP
#' `correct_for_chance`: whether to quantify how much the observed RR differs
#' from what would be expected under random feature selection (same feature set size)
#' `N`: number of random feature subsets to consider to correct for chance, i.e.
#' `N` random feature sets of the same size as the selected features are generated
compute_redundancy = function(task, train_set, selected_features, alpha = 0.05,
                              correct_for_chance = FALSE, N = 1000) {
  stopifnot(N >= 10)
  p = length(selected_features)

  if (p <= 1) {
    redundancy = tibble(
      rr_pearson = NA, srp_pearson = NA,
      rr_spearman = NA, srp_spearman = NA,
      rr_xicor = NA, srp_xicor = NA,
      arr_xicor = base::switch(correct_for_chance, NA)
    )
  } else {
    data = task$data(rows = train_set, cols = selected_features)

    combs = utils::combn(ncol(data), 2) # unique pairs only
    n_pairs = ncol(combs)

    pearson_coeff = numeric(n_pairs)
    pearson_pvalues = numeric(n_pairs)
    spearman_coeff = numeric(n_pairs)
    spearman_pvalues = numeric(n_pairs)
    xicor_coeff = numeric(n_pairs)
    xicor_pvalues = numeric(n_pairs)

    for (i in 1:n_pairs) {
      x = data[[combs[1, i]]]
      y = data[[combs[2, i]]]

      ptest = cor.test(x, y, method = "pearson")
      pearson_coeff[i] = ptest$estimate
      pearson_pvalues[i] = ptest$p.value

      stest = cor.test(x, y, method = "spearman", exact = FALSE)
      spearman_coeff[i] = stest$estimate
      spearman_pvalues[i] = stest$p.value

      # choose the max value of xicor, since it's not a symmetric measure of dependence
      # see (Chatterjee 2021) - Remark No. 1
      xc = list(xicor(x, y, pvalue = TRUE), xicor(y, x, pvalue = TRUE))
      indx = which.max(map_dbl(xc, "xi"))
      xicor_coeff[i] = xc[[indx]]$xi
      xicor_pvalues[i] = xc[[indx]]$pval
    }

    xi_obs = mean(abs(xicor_coeff))

    if (correct_for_chance) {
      # Estimate null distribution
      all_data = task$data(rows = train_set)
      xi_null = numeric(N)
      for (i in seq_len(N)) {
        samp = sample(ncol(all_data), size = p)
        combs_null = utils::combn(samp, 2)
        xi_vals = numeric(ncol(combs_null))
        for (j in seq_along(xi_vals)) {
          x = all_data[[combs_null[1, j]]]
          y = all_data[[combs_null[2, j]]]
          xi_scores = c(xicor(x, y), xicor(y, x))
          xi_scores = xi_scores[is.finite(xi_scores)] # remove Inf values
          xi_vals[j] = max(c(xi_scores, -Inf)) # just in case both scores are Inf and were removed!
        }
        xi_null[i] = mean(abs(xi_vals[is.finite(xi_vals)])) # only non-Inf values
      }

      xi_expected = mean(xi_null)
      xi_max = max(xi_null)
      adj_xi = (xi_obs - xi_expected) / (xi_max - xi_expected)
    }

    redundancy = tibble(
      #' `rr` as in "redundancy rate"
      #' `arr` as in "adjusted redundancy rate"
      #' `srp` as in "significant redundancy proportion"

      #' Pearson
      rr_pearson = mean(abs(pearson_coeff)),
      srp_pearson = sum(p.adjust(pearson_pvalues, "fdr") < alpha) / n_pairs,
      #' Spearman
      rr_spearman = mean(abs(spearman_coeff)),
      srp_spearman = sum(p.adjust(spearman_pvalues, "fdr") < alpha) / n_pairs,
      #' Xicor
      rr_xicor = xi_obs,
      srp_xicor = sum(p.adjust(xicor_pvalues, "fdr") < alpha) / n_pairs,
      #' Chance-corrected Xicor
      arr_xicor = base::switch(correct_for_chance, adj_xi)
    )
  }

  redundancy
}

# load the feature selection results per omic
fs = readRDS(file = "bench/fs.rds")

fs_long = fs |>
  pivot_longer(
    cols = ends_with("_feats"),
    names_to = "fs_method_id",
    values_to = "feats"
  ) |>
  mutate(fs_method_id = sub("_feats$", "", fs_method_id)) |> # cleanup fs method names
  select(dataset_id, omic_id, rsmp_id, fs_method_id, feats)

# get task lists beforehand for memory efficiency
dataset_ids = unique(fs$dataset_id)

task_lists = lapply(dataset_ids, function(id) {
  readRDS(file.path("data", id, "task_list.rds"))
})
names(task_lists) = dataset_ids

# same for subsamplings
subsamplings = lapply(dataset_ids, function(id) {
  readRDS(file.path("data", id, "subsampling.rds"))
})
names(subsamplings) = dataset_ids

# for reproducibility
set.seed(42)

# parallelization
plan(multicore, workers = 15) # change to e.g. 100 if `correct_for_chance = TRUE`

with_progress({
  total_rows = nrow(fs_long)
  p = progressor(along = 1:total_rows)

  fs_red = future.apply::future_lapply(seq_len(total_rows), function(i) {
    data = fs_long[i, ]

    dataset_id = data$dataset_id
    omic_id = data$omic_id
    rsmp_id = data$rsmp_id
    fs_method_id = data$fs_method_id
    selected_features = data$feats[[1]]

    p(sprintf("Dataset: %s, Omic: %s, SubSmp iter: %i, FS method: %s",
              dataset_id, omic_id, rsmp_id, fs_method_id))

    # get task
    task = task_lists[[dataset_id]][[omic_id]]

    # get train set
    train_set = subsamplings[[dataset_id]]$train_set(i = rsmp_id)

    # get redundancy scores
    r = compute_redundancy(task, train_set, selected_features)

    bind_cols(data, r) |> select(!feats)
  }, future.seed = TRUE) |> bind_rows()
})

# save result to disk
saveRDS(fs_red, file = "bench/fs_red.rds")
