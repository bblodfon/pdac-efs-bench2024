#' helper function to impute omics data, refactored from Wissel's script:
#' https://github.com/BoevaLab/Multi-omics-noise-resistance/blob/main/noise_resistance/R/prep/prepare_tcga_data.R
#' Imputation strategy: >10% NA → drop feature; else median-impute (after NA removal)
impute = function(df) {
  # Step 1: Drop features (columns) with >10% missing values across samples
  missing_rate = colMeans(is.na(df))
  too_many_missing = sum(missing_rate > 0.1)
  if (too_many_missing > 0) {
    cat(sprintf("%d features with more than 10%% NAs across samples were removed.\n",
                too_many_missing))
    df = df[, missing_rate <= 0.1, drop = FALSE]
  }

  # Step 2: Impute remaining missing values using column medians
  na_counts = colSums(is.na(df)) # feature-wise, how many NAs we have?
  feats_to_impute = sum(na_counts > 0)
  if (feats_to_impute > 0) {
    cat(sprintf("Imputing %d features with <=10%% missing values using median feature values.\n",
                feats_to_impute))
    for (j in which(na_counts > 0)) { # get feature indexes that have NAs
      df[is.na(df[, j]), j] = median(df[[j]], na.rm = TRUE)
    }
  }

  df
}
