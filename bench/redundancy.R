#' **Redundancy** analysis of feature selection methods applied to multi-omic datasets.
#' We do this analysis in a separate script since it takes some time even if parallelized.
#'
#' We calculate the pair-wise Pearson, Spearman and XI correlation (Chatterjee 2021) coefficients
#' for every selected feature set and provide the following summary scores:
#' - Mean absolute correlation between all unique feature pairs (Redundancy Rate or RR).
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

#' `alpha`: FDR threshold for SRP
#' `N`: number of random feature subsets to consider to correct for chance, i.e.
#' `N` random feature sets of the same size as the selected features are generated
compute_redundancy_per_method = function(row, task_list, omic_id, train_set, method,
                                         alpha = 0.05, N = 1000) {
  selected_features = row[[method]][[1]]
  p = length(selected_features)
  task = task_list[[omic_id]]

  if (p <= 1) {
    redundancy = c(
      rrate_pearson = NA, srp_pearson = NA,
      rrate_spearman = NA, srp_spearman = NA,
      rrate_xicor = NA, srp_xicor = NA, adj_xicor = NA
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

      test = cor.test(x, y, method = "pearson")
      pearson_coeff[i] = test$estimate
      pearson_pvalues[i] = test$p.value

      test = cor.test(x, y, method = "spearman", exact = FALSE)
      spearman_coeff[i] = test$estimate
      spearman_pvalues[i] = test$p.value

      # choose the max value of xicor, since it's not a symmetric measure of dependence
      # see (Chatterjee 2021) - Remark No. 1
      xc = list(xicor(x, y, pvalue = TRUE), xicor(y, x, pvalue = TRUE))
      indx = which.max(map_dbl(xc, "xi"))
      xicor_coeff[i] = xc[[indx]]$xi
      xicor_pvalues[i] = xc[[indx]]$pval
    }

    # Estimate null distribution
    all_data = task$data(rows = train_set)
    xi_null = numeric(N)
    for (i in seq_len(N)) {
      samp = sample(ncol(all_data), p)
      combs_null = utils::combn(samp, 2)
      xi_vals = numeric(ncol(combs_null))
      for (j in seq_along(xi_vals)) {
        x = all_data[[combs_null[1, j]]]
        y = all_data[[combs_null[2, j]]]
        xi_vals[j] = max(xicor(x, y), xicor(y, x))
      }
      xi_null[i] = mean(abs(xi_vals))
    }

    xi_expected = mean(xi_null)
    xi_max = max(xi_null)
    xi_obs = mean(abs(xicor_coeff))
    adj_xi = (xi_obs - xi_expected) / (xi_max - xi_expected)

    redundancy = c(
      #' `rrate` as in "redundancy rate"
      #' Pearson
      rrate_pearson = mean(abs(pearson_coeff)),
      srp_pearson = sum(p.adjust(pearson_pvalues, "fdr") < alpha) / n_pairs,
      #' Spearman
      rrate_spearman = mean(abs(spearman_coeff)),
      srp_spearman = sum(p.adjust(spearman_pvalues, "fdr") < alpha) / n_pairs,
      #' Xicor
      rrate_xicor = xi_obs,
      srp_xicor = sum(p.adjust(xicor_pvalues, "fdr") < alpha) / n_pairs,
      adj_xicor = adj_xi
    )
  }

  redundancy
}

compute_redundancy = function(row, task_list, subsampling, methods){
  redundancy = lapply(methods, function(method) {
    result = compute_redundancy_per_method(
      row = row,
      task_list = task_list,
      omic_id = row[["omic_id"]],
      train_set = subsampling$train_set(i = row[["rsmp_id"]]),
      method = method)
    names(result) = paste0(gsub("_feats", "", method), "_", names(result))

    result
  })

  unlist(redundancy)
}

# load the feature selection results per omic
fs = readRDS(file = "bench/fs.rds")
methods = colnames(select(fs, ends_with("_feats")))
dataset_ids = unique(fs$dataset_id)
fs_red = list()

# for reproducibility
set.seed(42)

# parallelization
plan(multicore, workers = 15)

fs_red = list()

with_progress({
  total_rows = nrow(fs)
  p = progressor(along = 1:total_rows)

  ## loop over datasets (outer loop sequential, inner parallel)
  for (dset in dataset_ids) {

    subsampling = readRDS(file.path("data", dset, "subsampling.rds"))
    task_list   = readRDS(file.path("data", dset, "task_list.rds"))
    fs_dset     = filter(fs, dataset_id == dset)

    red_dset = future.apply::future_lapply(seq_len(nrow(fs_dset)), function(i) {
      row = fs_dset[i, ]

      r = compute_redundancy(
        row         = row,
        task_list   = task_list,
        subsampling = subsampling,
        methods     = methods
      )

      p(sprintf("%s row %d/%d", dset, i, nrow(fs_dset)))
      as.data.frame(as.list(r))
    }, future.seed = TRUE) |> bind_rows()

    fs_red[[dset]] = bind_cols(fs_dset, red_dset)
  }
})

fs_red = bind_rows(fs_red) |> select(-ends_with("feats"), -ends_with("coxlasso_train_time"))

# save result to disk
saveRDS(fs_red, file = "bench/fs_red.rds")
