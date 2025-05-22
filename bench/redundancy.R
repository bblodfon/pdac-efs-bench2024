#' **Redundancy** analysis of feature selection methods applied to multi-omic datasets.
#'
#' We calculate the pair-wise Pearson, Spearman and XI correlation (Chatterjee 2021) coefficients
#' for every selected feature set and provide the following summary scores:
#' - Average of the absolute values of the coefficients.
#' - Number of significantly correlated variables (FDR-adjusted p-values, alpha = 0.05).
#' - Proportion of significantly correlated features divided by the total number of selected features.
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

compute_redundancy_per_method = function(row, task_list, omic_id, train_set,
                                         method, thres = 0.05, symmetric_xicor = TRUE) {
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

    combs = utils::combn(ncol(data), 2)
    pearson_coeff = numeric(ncol(combs))
    pearson_pvalues = numeric(ncol(combs))
    spearman_coeff = numeric(ncol(combs))
    spearman_pvalues = numeric(ncol(combs))
    xicor_coeff = numeric()
    xicor_pvalues = numeric()

    for (i in 1:ncol(combs)) {
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

      if (symmetric_xicor) {
        # choose the max value of xicor for directions according to
        # (Chatterjee 2021) paper Remark No. 1
        indx = which.max(map_dbl(xc, "xi"))
        xicor_coeff = c(xicor_coeff, xc[[indx]]$xi)
        xicor_pvalues = c(xicor_pvalues, xc[[indx]]$pval)
      } else {
        # add both (x,y) and (y,x) values to the result vector for hypothesis testing
        xicor_coeff = c(xicor_coeff, map_dbl(xc, "xi"))
        xicor_pvalues = c(xicor_pvalues, map_dbl(xc, "pval"))
      }
    }

    redundancy = c(
      values_pearson = mean(abs(pearson_coeff)),
      sign_pearson = sum(p.adjust(pearson_pvalues, "fdr") < thres),
      prop_sign_pearson = sum(p.adjust(pearson_pvalues, "fdr") < thres) / length(pearson_pvalues),
      values_spearman = mean(abs(spearman_coeff)),
      sign_spearman = sum(p.adjust(spearman_pvalues, "fdr") < thres),
      prop_sign_spearman = sum(p.adjust(spearman_pvalues, "fdr") < thres) / length(spearman_pvalues),
      values_xicor = mean(abs(xicor_coeff)),
      sign_xicor = sum(p.adjust(xicor_pvalues, "fdr") < thres),
      prop_sign_xicor = sum(p.adjust(xicor_pvalues, "fdr") < thres) / length(xicor_pvalues)
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
