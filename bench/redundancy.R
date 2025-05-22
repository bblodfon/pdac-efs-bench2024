#' **Redundancy** analysis of feature selection methods applied to multi-omic datasets.
#' We do this analysis in a separate script since it takes some time even if parallelized.
#'
#' We calculate the pair-wise Pearson, Spearman and XI correlation (Chatterjee 2021) coefficients
#' for every selected feature set and provide the following summary scores:
#' - Mean absolute correlation between all unique feature pairs (Redundancy Rate).
#' - Proportion of significant pairwise correlations, after multiple testing correction (FDR).
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

compute_redundancy_per_method = function(row, task_list, omic_id, train_set, method) {
  selected_features = row[[method]][[1]]

  if (length(selected_features) == 1) {
    redundancy = c(
      values_pearson = NA,
      sign_pearson = NA,
      prop_sign_pearson = NA,
      values_spearman = NA,
      sign_spearman = NA,
      prop_sign_spearman = NA,
      values_xicor = NA,
      sign_xicor = NA,
      prop_sign_xicor = NA
    )
  } else {
    data = task_list[[omic_id]]$data(rows = train_set, cols = selected_features)

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

      xc = list(
        xicor(x, y, pvalue = TRUE),
        xicor(y, x, pvalue = TRUE)
      )

      # choose the max value of xicor, see (Chatterjee 2021) - Remark No. 1
      indx = which.max(map_dbl(xc, "xi"))
      xicor_coeff[i] = xc[[indx]]$xi
      xicor_pvalues[i] = xc[[indx]]$pval
    }

    redundancy = c(
      #' `rrate` as in "redundancy rate"
      #' Pearson
      rrate_pearson = mean(abs(pearson_coeff)),
      prop5_pearson = sum(p.adjust(pearson_pvalues, "fdr") < 0.05) / n_pairs,
      prop1_pearson = sum(p.adjust(pearson_pvalues, "fdr") < 0.01) / n_pairs,
      #' Spearman
      rrate_spearman = mean(abs(spearman_coeff)),
      prop5_spearman = sum(p.adjust(spearman_pvalues, "fdr") < 0.05) / n_pairs,
      prop1_spearman = sum(p.adjust(spearman_pvalues, "fdr") < 0.01) / n_pairs,
      #' Xicor
      rrate_xicor = mean(abs(xicor_coeff)),
      prop5_xicor = sum(p.adjust(xicor_pvalues, "fdr") < 0.05) / n_pairs,
      prop1_xicor = sum(p.adjust(xicor_pvalues, "fdr") < 0.01) / n_pairs
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
        methods     = methods)

      p(sprintf("%s row %d/%d", dset, i, nrow(fs_dset)))
      as.data.frame(as.list(r))
    }, future.seed = TRUE) |> bind_rows()

    fs_red[[dset]] = bind_cols(fs_dset, red_dset)
  }
})

fs_red = bind_rows(fs_red) |> select(-ends_with("feats"), -ends_with("coxlasso_train_time"))

# save result to disk
saveRDS(fs_red, file = "bench/fs_red.rds")
