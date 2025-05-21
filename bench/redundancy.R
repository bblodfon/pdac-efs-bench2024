library("XICOR")
library("dplyr")
library("tidyr")


compute_redundancy_per_method = function(row, task_list, omic_id, train_set, method){
  selected_features <- row[[method]][[1]]
  
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
  } else{
    data = task_list[[omic_id]]$data(rows = train_set, cols = selected_features)
    
    combs <- combn(ncol(data), 2)
    pearson_coeff <- numeric(ncol(combs))
    pearson_pvalues <- numeric(ncol(combs))
    spearman_coeff <- numeric(ncol(combs))
    spearman_pvalues <- numeric(ncol(combs))
    xicor_coeff <- numeric()
    xicor_pvalues <- numeric()
    
    for (i in 1:ncol(combs)) {
      x <- data[[combs[1, i]]]
      y <- data[[combs[2, i]]]
      test <- cor.test(x, y, method = "pearson")
      pearson_coeff[i] <- test$estimate
      pearson_pvalues[i] <- test$p.value
      test <- cor.test(x, y, method = "spearman", exact=FALSE)
      spearman_coeff[i] <- test$estimate
      spearman_pvalues[i] <- test$p.value
      test <- xicor(x, y, pvalue = TRUE)
      xicor_coeff <- c(xicor_coeff, test$xi)
      xicor_pvalues <- c(xicor_pvalues, test$pval)
      test <- xicor(y, x, pvalue = TRUE)
      xicor_coeff <- c(xicor_coeff, test$xi)
      xicor_pvalues <- c(xicor_pvalues, test$pval)
    }
    
    redundancy = c(
      values_pearson = mean(abs(pearson_coeff)),
      sign_pearson = sum(p.adjust(pearson_pvalues, "fdr") < 0.05),
      prop_sign_pearson = sum(p.adjust(pearson_pvalues, "fdr") < 0.05) / length(pearson_pvalues),
      values_spearman = mean(abs(spearman_coeff)),
      sign_spearman = sum(p.adjust(spearman_pvalues, "fdr") < 0.05),
      prop_sign_spearman = sum(p.adjust(spearman_pvalues, "fdr") < 0.05) / length(spearman_pvalues),
      values_xicor = mean(abs(xicor_coeff)),
      sign_xicor = sum(p.adjust(xicor_pvalues, "fdr") < 0.05),
      prop_sign_xicor = sum(p.adjust(xicor_pvalues, "fdr") < 0.05) / length(xicor_pvalues)
    )
    
  }
  
  return(redundancy)
}


compute_redundancy = function(row, task_list, subsampling, methods){
  redundancy = lapply(methods, function(method){
    result = compute_redundancy_per_method(
      row = row,
      task_list = task_list,
      omic_id = row[["omic_id"]],
      train_set = subsampling$train_set(i = row[["rsmp_id"]]),
      method = method)
    names(result) <- paste0(gsub("_feats", "", method), "_", names(result))
    return(result)
  })
  redundancy = unlist(redundancy)
  return(redundancy)
}


fs = readRDS(file = "bench/fs.rds")
methods = colnames(select(fs, ends_with("_feats")))
dataset_ids <- unique(fs$dataset_id)
fs_final <- list()

for (dataset_id in dataset_ids) {
  subsampling = readRDS(file = file.path("data", dataset_id, "subsampling.rds"))
  task_list = readRDS(file.path("data", dataset_id, "task_list.rds"))
  fs_dataset = filter(fs, dataset_id == !!dataset_id)
  redundancy = bind_rows(lapply(1:nrow(fs_dataset), function(i) {
    row = fs_dataset[i, ]
    result = compute_redundancy(
      row = row,
      task_list = task_list,
      subsampling = subsampling,
      methods = methods
    )
    as.data.frame(as.list(result))
  }))
  
  fs_final[[dataset_id]] <- bind_cols(fs_dataset, redundancy)
}
fs_final = bind_rows(fs_final)

saveRDS(fs_final, file = "bench/fs_red.rds")
