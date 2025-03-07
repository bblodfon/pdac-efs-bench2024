# some test figures
suppressPackageStartupMessages({
  library(mlr3viz)
  library(data.table)
  library(mlr3misc)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)  # For creating custom color scale
})

# per learner class
efs_tree = readRDS(file = "efs/gex/efs_tree.rds")
efs_coxlasso = readRDS(file = "efs/gex/efs_coxlasso.rds")
efs_glmb = readRDS(file = "efs/gex/efs_glmb.rds")
efs_rsf = readRDS(file = "efs/gex/efs_rsf.rds")

# All results combined
efs_all = readRDS(file = "efs/gex/efs_all.rds") # GEX TCGA
efs_all = readRDS(file = "efs/snv/efs_all.rds") # SNV Osipov

## pareto
autoplot(efs_all, type = "pareto") +
  scale_color_brewer(palette = "Set1")
autoplot(efs_all, type = "pareto", pareto_front = "estimated") +
  scale_color_brewer(palette = "Set1")
## how many features?
efs_all$knee_points()
efs_all$knee_points(type = "estimated")

## performance
autoplot(efs_all, type = "performance") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed")
## nfeatures
autoplot(efs_all, type = "n_features") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
## stability
autoplot(efs_all, type = "stability", stability_measure = "nogueira", stability_args = list(p = length(efs_all$.__enclos_env__$private$.features))) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## compare feature ranking methods
nfeats_cutoff = efs_all$knee_points(type = "estimated")$n_features
methods = list("av", "sav", "seq_pav", "seq_phragmen")

features_list = lapply(methods, function(method) {
  efs_all$feature_ranking(method = method, use_weights = TRUE, committee_size = nfeats_cutoff)[["feature"]]
})
names(features_list) = methods
n_methods = length(methods)

# Jaccard similarity between feature ranking methods
jaccard_matrix = matrix(0, nrow = n_methods, ncol = n_methods)
rownames(jaccard_matrix) = methods
colnames(jaccard_matrix) = methods

for (i in 1:n_methods) {
  for (j in i:n_methods) {
    score = stabm::stabilityJaccard(features = list(features_list[[i]], features_list[[j]]))
    jaccard_matrix[i, j] = score
    jaccard_matrix[j, i] = score  # Fill the symmetric part
  }
}

breaks = c(min(jaccard_matrix), (max(jaccard_matrix) + min(jaccard_matrix))/2, max(jaccard_matrix))
# breaks = c(0.5, 0.75, 1)
col_fun = colorRamp2(breaks = breaks, colors = c("#FF6B6B", "white", "#4CAF50"))
ComplexHeatmap::Heatmap(
  matrix = jaccard_matrix,
  name = "Jaccard Index",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.5),
  column_names_rot = 45,
  heatmap_legend_param = list(title = "Jaccard")
)
