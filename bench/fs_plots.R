#' **Sparsity and stability** analysis of feature selection methods applied to
#' multi-omic datasets.
#' We visualize how many features were selected per omic type and how consistent
#' the selections were across resamplings.
#'
#' Execute: `Rscript bench/fs_plots.R` (from project root)

suppressPackageStartupMessages({
  library(ggplot2)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(stringr)
  library(stabm)
  library(mlr3misc)
})

# load the feature selection results per omic ----
fs = readRDS(file = "bench/fs.rds")

## Per-Omic Sparsity ----
# hue_colors = scale_color_hue()$palette(n = 5)
# more manual hue colors:
# hue_colors = c("#F8766D", "#A3A500", "#00BB4E", "#E76BF3", "#35A2FF", "#9590FF")

# RColorBrewer::brewer.pal(n = 9, name = "Set1")
set1_colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")

custom_colors = c(
  "hEFS (9 models)" = set1_colors[1],
  "EFS (CoxLasso)" = set1_colors[2],
  "hEFS (3 RSFs)" = set1_colors[3],
  "CoxLasso" = set1_colors[4]
)

# dataset labels
osipov_abbrev = "Osipov et al. (2024)"
wissel_abbrev = "Wissel et al. (2023)"
dataset_labels = c(osipov2024 = osipov_abbrev, wissel2023 = wissel_abbrev)

fs_long = fs |>
  select(dataset_id, omic_id, ends_with("nfeats")) |>
  pivot_longer(
    cols = starts_with("efs") | starts_with("coxlasso"),
    names_to = "method",
    values_to = "n_feats"
  ) |>
  mutate(method = case_when(
    method == "coxlasso_nfeats" ~ "CoxLasso",
    method == "efs_all_nfeats" ~ "hEFS (9 models)",
    method == "efs_coxlasso_nfeats" ~ "EFS (CoxLasso)",
    method == "efs_rsf_nfeats" ~ "hEFS (3 RSFs)",
    TRUE ~ method  # Keep other values unchanged
  ))

# Reorder omic_id within each dataset based on median number of features
fs_long |>
  group_by(dataset_id, omic_id) |>
  mutate(method = fct_reorder(method, n_feats, .fun = median, .desc = TRUE)) |>
  ungroup() |>
  ggplot(aes(x = omic_id, y = n_feats, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  geom_jitter(aes(color = method), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 0.1) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  facet_wrap(
    ~ dataset_id,
    scales = "free_x",  # One plot per dataset_id
    labeller = as_labeller(dataset_labels)
  ) +
  labs(
    x = "Omics",
    y = "Number of Features",
    fill = "Feature Selection\nMethod",
    title = "Feature Selection Sparsity across Omic Types"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("bench/img/sparsity.png", width = 7, height = 5, dpi = 600, bg = "white")

# Same plot, without coxlasso fs
fs_long |>
  filter(method != "CoxLasso") |>
  group_by(dataset_id, omic_id) |>
  mutate(method = fct_reorder(method, n_feats, .fun = median, .desc = TRUE)) |>
  ungroup() |>
  ggplot(aes(x = omic_id, y = n_feats, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  geom_jitter(aes(color = method), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 0.1) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  facet_wrap(
    ~ dataset_id,
    scales = "free_x",  # One plot per dataset_id
    labeller = as_labeller(dataset_labels)
  ) +
  labs(
    x = "Omics",
    y = "Number of Features",
    fill = "Feature Selection Method",
    title = "Feature Selection Sparsity across Omic Types"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("bench/img/sparsity_no_coxlasso.png", width = 7, height = 5, dpi = 600, bg = "white")

## Per-Omic Stability ----
# Long format per method
fs_long = fs |>
  select(dataset_id, omic_id, rsmp_id, ends_with("_feats")) |>
  pivot_longer(
    cols = ends_with("_feats"),
    names_to = "method",
    values_to = "features"
  ) |>
  mutate(method = case_when(
    method == "coxlasso_feats" ~ "CoxLasso",
    method == "efs_all_feats" ~ "hEFS (9 models)",
    method == "efs_coxlasso_feats" ~ "EFS (CoxLasso)",
    method == "efs_rsf_feats" ~ "hEFS (3 RSFs)",
    TRUE ~ method  # Keep other values unchanged
  ))

### Jaccard ----
# Assess stability across all 100 resamplings - 1 value per (dataset, omic) combo
stab_summary_jacc = fs_long |>
  group_by(dataset_id, omic_id, method) |>
  summarize(
    stability = stabm::stabilityJaccard(features),
    .groups = "drop"
  )

stab_summary_jacc |>
  group_by(dataset_id, omic_id) |>
  mutate(method = fct_reorder(method, stability, .fun = median, .desc = TRUE)) |>
  ungroup() |>
  ggplot(aes(x = omic_id, y = stability, fill = method)) +
    geom_col(position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = custom_colors) +
    theme_minimal() +
    facet_wrap(
      ~ dataset_id,
      scales = "free_x",
      labeller = as_labeller(dataset_labels)
    ) +
    ylim(c(0, 0.65)) +
    labs(
      x = "Omics",
      y = "Jaccard Similarity",
      fill = "Feature Selection\nMethod",
      title = "Feature Selection Stability across Omic Types"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("bench/img/stability_jaccard.png", width = 7, height = 5, dpi = 600, bg = "white")

### Nogueira ----
# make a nice table with (dataset, omic, #features) since the Nogueira measure needs `p`
dataset_ids = unique(fs$dataset_id)
p_lookup = mlr3misc::map_dtr(dataset_ids, function(dataset_id) {
  task_list = readRDS(file.path("data", dataset_id, "task_list.rds"))
  omic_ids = readr::read_csv(
    file.path("data", dataset_id, "omic_ids.csv"),
    col_names = FALSE, show_col_types = FALSE, progress = FALSE
  )[[1L]]

  tibble::tibble(
    dataset_id = dataset_id,
    omic_id = omic_ids,
    p = mlr3misc::map_int(omic_ids, function(omic_id) task_list[[omic_id]]$n_features)
  )
})

stab_summary_nog = fs_long |>
  left_join(p_lookup, by = c("dataset_id", "omic_id")) |>
  group_by(dataset_id, omic_id, method) |>
  summarize(
    stability = stabm::stabilityNogueira(features, p = p[1]),
    .groups = "drop"
  )

stab_summary_nog |>
  group_by(dataset_id, omic_id) |>
  mutate(method = fct_reorder(method, stability, .fun = median, .desc = TRUE)) |>
  ungroup() |>
  ggplot(aes(x = omic_id, y = stability, fill = method)) +
  geom_col(position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  facet_wrap(
    ~ dataset_id,
    scales = "free_x",
    labeller = as_labeller(dataset_labels)
  ) +
  ylim(c(0, 0.65)) +
  labs(
    x = "Omics",
    y = "Nogueira Similarity",
    fill = "Feature Selection\nMethod",
    title = "Feature Selection Stability across Omic Types"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("bench/img/stability_nogueira.png", width = 7, height = 5, dpi = 600, bg = "white")

### Variability - Jaccard ----
# We repeatedly subsample the 100 resamplings (e.g., draw 50 at random) and
# compute the stability each time

# Set number of subsamples and size of each subsample
n_subsamples = 50
subsample_size = 50

set.seed(42)
stab_jacc_rsmp = fs_long |>
  group_by(dataset_id, omic_id, method) |>
  summarize(
    # Do n_subsamples stability computations per group
    stability = list(
      replicate(n_subsamples, {
        sampled_feats = sample(features, subsample_size)
        stabm::stabilityJaccard(sampled_feats)
      }, simplify = TRUE)
    ),
    .groups = "drop"
  ) |>
  unnest_longer(stability)

stab_jacc_rsmp |>
  group_by(dataset_id, omic_id) |>
  mutate(method = fct_reorder(method, stability, .fun = median, .desc = TRUE)) |>
  ungroup() |>
  ggplot(aes(x = omic_id, y = stability, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  geom_jitter(aes(color = method), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 0.1) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  facet_wrap(
    ~ dataset_id,
    scales = "free_x",
    labeller = as_labeller(dataset_labels)
  ) +
  labs(
    x = "Omics",
    y = "Jaccard Similarity",
    fill = "Feature Selection\nMethod",
    title = "Feature Selection Stability across Omic Types"
  ) +
  ylim(c(0, 0.65)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("bench/img/stability_jaccard_var.png", width = 7, height = 5, dpi = 600, bg = "white")

### Variability - Noqueira ----
set.seed(42)
stab_nog_rsmp = fs_long |>
  left_join(p_lookup, by = c("dataset_id", "omic_id")) |>
  group_by(dataset_id, omic_id, method) |>
  summarize(
    stability = list(
      replicate(n_subsamples, {
        sampled_feats = sample(features, subsample_size)
        stabm::stabilityNogueira(sampled_feats, p = p[1])
      }, simplify = TRUE)
    ),
    .groups = "drop"
  ) |>
  unnest_longer(stability)

stab_nog_rsmp |>
  group_by(dataset_id, omic_id) |>
  mutate(method = fct_reorder(method, stability, .fun = median, .desc = TRUE)) |>
  ungroup() |>
  ggplot(aes(x = omic_id, y = stability, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  geom_jitter(aes(color = method), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 0.1) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  facet_wrap(
    ~ dataset_id,
    scales = "free_x",
    labeller = as_labeller(dataset_labels)
  ) +
  labs(
    x = "Omics",
    y = "Noqueira Similarity",
    fill = "Feature Selection\nMethod",
    title = "Feature Selection Stability across Omic Types"
  ) +
  ylim(c(0, 0.67)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave("bench/img/stability_nogueira_var.png", width = 7, height = 5, dpi = 600, bg = "white")
