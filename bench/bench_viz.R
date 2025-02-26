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

# Boxplots of #features
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
    # dodge ensures boxplots for different methods stay side by side
    geom_boxplot(outlier.shape = NA, alpha = 0.7,
                 position = position_dodge2(preserve = "single")) +
    scale_fill_manual(values = custom_colors) +
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
  # dodge ensures boxplots for different methods stay side by side
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = custom_colors) +
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
