library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

# feature selection data
fs = readRDS(file = "bench/fs.rds")

hue_colors = scale_color_hue()$palette(n = 4)

# Density plots of #features ----
fs |>
  select(efs_all_nfeats, efs_coxlasso_nfeats, efs_rsf_nfeats, coxlasso_nfeats) |>
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Features") |>
  filter(Features <= 50) |>
  ggplot(aes(x = Features, fill = Method, color = Method)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, alpha = 0.3,
                   position = "identity") +
    #geom_density(alpha = 0.5, linewidth = 1) +
    theme_minimal() +
    labs(title = "Feature Selection Distributions (<= 50 features)",
         x = "Number of Features",
         y = "Density") +
    scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
    theme(legend.title = element_blank(),
          legend.position = "top")

# Sparsity Boxplots ----
custom_colors = c(
  "efs_all_nfeats" = hue_colors[1],
  "efs_coxlasso_nfeats" = hue_colors[2],
  "efs_rsf_nfeats" = hue_colors[3],
  "coxlasso_nfeats" = hue_colors[4]
)

long_fs = fs |>
  select(dataset_id, omic_id, ends_with("nfeats")) |>
  pivot_longer(
    cols = starts_with("efs") | starts_with("coxlasso"),
    names_to = "method",
    values_to = "n_feats"
  )

# Reorder omic_id within each dataset based on median number of features
long_fs2 = long_fs |>
  group_by(dataset_id, omic_id) |>
  mutate( #  Use median for ordering
    omic_id = fct_reorder(omic_id, n_feats, .fun = median, .desc = TRUE)
  ) |>
  ungroup()

long_fs2 |>
  ggplot(aes(x = omic_id, y = n_feats, fill = method)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7,
                 position = position_dodge2(preserve = "single")) +
    geom_jitter(aes(color = method, size = n_feats), show.legend = FALSE,
                position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
                alpha = 0.7, size = 0.1) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = custom_colors) +
    theme_minimal() +
    facet_wrap(~ dataset_id, scales = "free_x") + # One plot per dataset_id
    labs(
      x = "Omics",
      y = "Number of Features (Sparsity)",
      fill = "Feature Selection Method",
      title = "Feature Selection Sparsity across Omic Types and Methods"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels for readability
    )

# remove coxlasso
long_fs3 = long_fs |>
  filter(method != "coxlasso_nfeats") |>
  group_by(dataset_id, omic_id) |>
  mutate( #  Use median for ordering
    omic_id = fct_reorder(omic_id, n_feats, .fun = median, .desc = TRUE)
  ) |>
  ungroup()

long_fs3 |>
  ggplot(aes(x = omic_id, y = n_feats, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  geom_jitter(aes(color = method, size = n_feats), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 0.01) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  geom_hline(yintercept = 10, color = "red", linetype = "dashed", size = 0.8) +
  theme_minimal() +
  facet_wrap(~ dataset_id, scales = "free_x") + # One plot per dataset_id
  labs(
    x = "Omics",
    y = "Number of Features (Sparsity)",
    fill = "Feature Selection Method",
    title = "Feature Selection Sparsity across Omic Types and Methods"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels for readability
  )

# Performance Boxplots ----
result = readRDS(file = "bench/result.rds")
sort(unlist(mlr3misc::map(result$coxlasso_feats, length))) # no zeros ok

library(ggplot2)
library(dplyr)
library(tidyr)

# Convert data to long format for ggplot
result_long = result |>
  pivot_longer(cols = c(harrel_c, uno_c, dcalib, ibrier),
               names_to = "measure", values_to = "value") |>
  filter(!measure %in% c("dcalib", "ibrier")) # Exclude these metrics

result_long |>
  ggplot(aes(x = fs_method_id, y = value, fill = fs_method_id)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(aes(color = fs_method_id), width = 0.2, alpha = 0.7, size = 0.1) +
    scale_fill_manual(values = c("efs_all_feats" = "#1f77b4", "coxlasso_feats" = "#ff7f0e")) +
    scale_color_manual(values = c("efs_all_feats" = "#1f77b4", "coxlasso_feats" = "#ff7f0e")) +
    geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", size = 0.8) +
    theme_minimal() +
    facet_grid(measure ~ dataset_id, scales = "free_y") +
    labs(
      x = "",
      y = "Performance Score",
      fill = "FS Method",
      color = "FS Method"
    ) +
    ylim(c(0, 1)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
