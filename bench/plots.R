suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(stringr)
})

# FEATURE SELECTION RESULTS (PER OMIC) ----
fs = readRDS(file = "bench/fs.rds")

# hue_colors = scale_color_hue()$palette(n = 5)
# more manual hue colors:
hue_colors = c(
  "#F8766D",
  "#A3A500",
  "#00BB4E",
  "#E76BF3",
  "#35A2FF",
  "#9590FF"
)

# RColorBrewer::brewer.pal(n = 9, name = "Set1")
set1_colors = c(
  "#E41A1C",
  "#377EB8",
  "#4DAF4A",
  "#984EA3",
  "#FF7F00",
  "#FFFF33",
  "#A65628",
  "#F781BF",
  "#999999"
)

# Per-Omic Sparsity ----
custom_colors = c(
  "EFS (9 models)" = set1_colors[1],
  "EFS (CoxLasso)" = set1_colors[2],
  "EFS (3 RSFs)" = set1_colors[3],
  "CoxLasso" = set1_colors[4]
)

fs_long = fs |>
  select(dataset_id, omic_id, ends_with("nfeats")) |>
  pivot_longer(
    cols = starts_with("efs") | starts_with("coxlasso"),
    names_to = "method",
    values_to = "n_feats"
  ) |>
  mutate(method = case_when(
    method == "coxlasso_nfeats" ~ "CoxLasso",
    method == "efs_all_nfeats" ~ "EFS (9 models)",
    method == "efs_coxlasso_nfeats" ~ "EFS (CoxLasso)",
    method == "efs_rsf_nfeats" ~ "EFS (3 RSFs)",
    TRUE ~ method  # Keep other values unchanged
  ))

# Reorder omic_id within each dataset based on median number of features
fs_long |>
  group_by(dataset_id, omic_id) |>
  mutate(
    method = fct_reorder(method, n_feats, .fun = median, .desc = TRUE)
  ) |>
  ungroup() |>
  #ggplot(aes(x = method, y = n_feats, fill = omic_id)) +
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
      labeller = as_labeller(
        c(osipov2024 = "Osipov et. al (2024)",
          wissel2023 = "Wissel et. al (2023)")
      )
    ) +
    labs(
      x = "Omics",
      y = "Number of Features",
      fill = "Feature Selection\nMethod",
      title = "Feature Selection Sparsity across Omic Types"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
ggsave("bench/img/omic_sparsity.png", width = 7, height = 5, dpi = 600, bg = "white")

# Same plot, without coxlasso fs
fs_long |>
  filter(method != "CoxLasso") |>
  group_by(dataset_id, omic_id) |>
  mutate(
    method = fct_reorder(method, n_feats, .fun = median, .desc = TRUE)
  ) |>
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
    labeller = as_labeller(
      c(osipov2024 = "Osipov et. al (2024)",
        wissel2023 = "Wissel et. al (2023)")
    )
  ) +
  labs(
    x = "Omics",
    y = "Number of Features",
    fill = "Feature Selection Method",
    title = "Feature Selection Sparsity across Omic Types"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave("bench/img/omic_sparsity_no_coxlasso.png", width = 7, height = 5, dpi = 600, bg = "white")

# BENCHMARK RESULTS ----
result = readRDS(file = "bench/result.rds") # all omics
#result = readRDS(file = "bench/result_no_mut.rds") # no mutation (Wissel)
#result = readRDS(file = "bench/result_no_mut_or_cnv.rds") # no mutation or CNV (Wissel)
#result = readRDS(file = "bench/result_with_lymph.rds") # + extra clinical variable

# Convert data to long format for ggplot
result_long = result |>
  # rename FS methods
  mutate(fs_method_id = case_when(
    fs_method_id == "coxlasso_feats" ~ "CoxLasso",
    fs_method_id == "efs_all_feats" ~ "EFS (9 models)",
    fs_method_id == "efs_coxlasso_feats" ~ "EFS (CoxLasso)",
    fs_method_id == "efs_rsf_feats" ~ "EFS (3 RSFs)",
    # `NA` fs method means that the 4 above FS methods were not used at all
    is.na(fs_method_id) ~ "Baseline",
    TRUE ~ fs_method_id  # Keep other values unchanged
  )) |>
  # rename model-data-configs
  rename(model = model_data_config) |>
  mutate(model = case_when(
    model == "cox-clinical" ~ "Cox (Clinical)",
    model == "rsf-clinical" ~ "RSF (Clinical)",
    model == "coxlasso-clinical+gex" ~ "CoxLasso (Clinical + GEX)",
    model == "rsf-clinical+gex" ~ "RSF (Clinical + GEX)",
    model == "coxlasso-all" ~ "CoxLasso (ALL)",
    model == "rsf-all" ~ "RSF (ALL)",
    TRUE ~ model  # Keep other values unchanged
  )) |>
  pivot_longer(cols = c(harrell_c, uno_c, brier_t12, brier_t24, brier_tmax24),
               names_to = "measure", values_to = "value")

# Multi-omics (MM) Sparsity ----
custom_colors = c(
  "EFS (9 models)" = set1_colors[1],
  "EFS (CoxLasso)" = set1_colors[2],
  "EFS (3 RSFs)" = set1_colors[3],
  "CoxLasso" = set1_colors[4]
)

result_long |>
  select(dataset_id, fs_method_id, rsmp_id, task_nfeats) |> # `model` doesn't play any role here
  filter(!fs_method_id %in% "Baseline") |> # remove Baseline
  distinct(dataset_id, fs_method_id, rsmp_id, .keep_all = TRUE) |> # remove duplicates
  #group_by(fs_method_id) |>
  mutate(
    fs_method_id = fct_reorder(fs_method_id, task_nfeats, .fun = median, .desc = TRUE)
  ) |>
  ungroup() |>
  ggplot(aes(x = fs_method_id, y = task_nfeats, fill = fs_method_id)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(aes(color = fs_method_id), show.legend = FALSE,
                position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.75),
                alpha = 0.7, size = 0.01) +
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values = custom_colors) +
    theme_minimal(base_size = 14) +
    facet_wrap(
      ~ dataset_id,
      scales = "free_x",  # One plot per dataset_id
      labeller = as_labeller(
        c(osipov2024 = "Osipov et. al (2024)",
          wissel2023 = "Wissel et. al (2023)")
      )
    ) +
    theme(
      axis.text.x = element_blank()
      #axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels for readability
    ) +
    labs(
      x = "",
      y = "Number of Features",
      #title = "Sparsity of Multi-omics Datasets",
      fill = "Feature Selection\nMethod"
    )
ggsave("bench/img/mm_sparsity.png", width = 7, height = 5, dpi = 600, bg = "white")

# Multi-omics FS contribution ----
extract_omics = function(features) {
  case_when(
    str_starts(features, "cnv_") ~ "CNV",
    str_starts(features, "snv_") ~ "SNV",
    str_starts(features, "indel") ~ "Insertions/Deletions",
    str_starts(features, "path_") ~ "Pathology",
    str_starts(features, "gex_") ~ "Gene Expression",
    str_starts(features, "meth_") ~ "Methylation",
    str_starts(features, "mutation_") ~ "Mutation",
    str_starts(features, "rppa_") ~ "RPPA",
    TRUE ~ "Clinical" # everything else
  )
}

ctr_res = result |>
  # remove duplicated or configs that don't have multi-omics fs incorporated
  filter(model_data_config == "coxlasso-all") |>
  select(dataset_id, fs_method_id, rsmp_id, task_nfeats, task_feats) |>
  unnest(task_feats) |> # smart, every row is per feature
  mutate(omics_type = extract_omics(task_feats)) |>
  group_by(dataset_id, fs_method_id, rsmp_id, omics_type, task_nfeats) |>
  summarise(count = n(), .groups = "drop") |>
  mutate(percentage = count / task_nfeats) |>
  group_by(dataset_id, fs_method_id, omics_type) |>
  summarise(avg_percentage = mean(percentage), avg_count = round(mean(count)), .groups = "drop") |>
  mutate(fs_method_id = case_when(
    fs_method_id == "coxlasso_feats" ~ "CoxLasso",
    fs_method_id == "efs_all_feats" ~ "EFS (9 models)",
    fs_method_id == "efs_coxlasso_feats" ~ "EFS (CoxLasso)",
    fs_method_id == "efs_rsf_feats" ~ "EFS (3 RSFs)",
    TRUE ~ fs_method_id  # there shouldn't be any other category here
  ))

# Stacked percentage plot
ctr_res |>
  ggplot(aes(x = reorder(fs_method_id, -avg_percentage, sum),
             y = avg_percentage, fill = omics_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = set1_colors) +
  facet_wrap(~ dataset_id,
    labeller = as_labeller(
      c(osipov2024 = "Osipov et. al (2024)",
        wissel2023 = "Wissel et. al (2023)")
    )
  ) +
  labs(
    x = "Feature Selection Method",
    y = "Average Contribution (%)",
    fill = "Omic Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )
ggsave("bench/img/omic_contribution_stacked.png", width = 8, height = 5, dpi = 600, bg = "white")

# Average counts barplot
## order by average total number of features per omic
ctr_res2 = ctr_res |>
  filter(omics_type != "Clinical") |>
  group_by(dataset_id, omics_type) |>
  summarize(avg_total = sum(avg_count), .groups = "drop") |>
  arrange(dataset_id, desc(avg_total)) |>
  mutate(omics_type = factor(omics_type, levels = unique(omics_type))) |>
  select(-avg_total) |>
  left_join(ctr_res, by = c("dataset_id", "omics_type")) |>
  mutate(omics_type = factor(omics_type, levels = unique(omics_type)))

ctr_res2 |>
  ggplot(aes(x = reorder(fs_method_id, -avg_count, sum),
             y = avg_count, fill = omics_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = set1_colors) +
  facet_wrap(
    ~ dataset_id,
    scales = "free_x",
    labeller = as_labeller(
      c(osipov2024 = "Osipov et. al (2024)",
        wissel2023 = "Wissel et. al (2023)")
    )
  ) +
  labs(
    x = "Feature Selection Method",
    y = "Average Number of Features",
    fill = "Omic Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )
ggsave("bench/img/omic_contribution_counts.png", width = 8, height = 5, dpi = 600, bg = "white")

# C-index ----
# some stats
result_long |>
  filter(measure == "harrell_c") |>
  group_by(dataset_id, fs_method_id, model) |>
  summarize(cindex = median(value)) |>
  arrange(dataset_id, fs_method_id, model) |>
  print(n = 22)

# Harrell's C
result_long |>
  filter(measure == "harrell_c") |>
  mutate(fs_method_id = factor(
    fs_method_id,
    levels = c("Baseline", "CoxLasso", "EFS (9 models)", "EFS (3 RSFs)", "EFS (CoxLasso)"))
  ) |>
  ggplot(aes(x = fs_method_id, y = value, fill = model)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(aes(color = model), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 0.01) +
    #scale_fill_manual(values = custom_colors) +
    #scale_color_manual(values = custom_colors) +
    geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", linewidth = 0.8) +
    theme_minimal(base_size = 14) +
    # `measure ~ dataset_id` if more measures
    facet_wrap(~ dataset_id,
      labeller = as_labeller(
        c(osipov2024 = "Osipov et. al (2024)",
          wissel2023 = "Wissel et. al (2023)")
      )
    ) +
    labs(
      x = "Feature Selection Method",
      y = "C-index",
      fill = "Integration Model (Data)"
    ) +
    ylim(c(0, 1)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
ggsave("bench/img/cindex.png", width = 9, height = 5, dpi = 600, bg = "white")

# Uno's C
result_long |>
  filter(measure == "uno_c") |>
  mutate(fs_method_id = factor(
    fs_method_id,
    levels = c("Baseline", "CoxLasso", "EFS (9 models)", "EFS (3 RSFs)", "EFS (CoxLasso)"))
  ) |>
  ggplot(aes(x = fs_method_id, y = value, fill = model)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = model), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 0.01) +
  #scale_fill_manual(values = custom_colors) +
  #scale_color_manual(values = custom_colors) +
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", linewidth = 0.8) +
  theme_minimal(base_size = 14) +
  # `measure ~ dataset_id` if more measures
  facet_wrap(
    ~ dataset_id,
    labeller = as_labeller(
      c(osipov2024 = "Osipov et. al (2024)",
        wissel2023 = "Wissel et. al (2023)")
      )
  ) +
  labs(
    x = "Feature Selection Method",
    y = "Uno's C-index",
    fill = "Integration Model (Data)"
  ) +
  ylim(c(0, 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )
ggsave("bench/img/uno_cindex.png", width = 9, height = 5, dpi = 600, bg = "white")

# IBS ----
result_long |>
  #filter(grepl("brier", measure)) |>
  filter(measure == "brier_tmax24") |> # `brier_t24`, `brier_tmax24`
  ggplot(aes(x = fs_method_id, y = value, fill = model)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = model), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 0.01) +
  #scale_fill_manual(values = custom_colors) +
  #scale_color_manual(values = custom_colors) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.8) +
  theme_minimal(base_size = 14) +
  facet_grid(
    measure ~ dataset_id,
    labeller = as_labeller(
      c(brier_t12 = "t = 1 year",
        brier_t24 = "t = 2 years",
        brier_tmax24 = "Up to t = 2 years",
        osipov2024 = "Osipov et. al (2024)",
        wissel2023 = "Wissel et. al (2023)"
      )
    ),
  ) +
  labs(
    x = "Feature Selection Method",
    y = "IBS ERV",
    fill = "Integration Model (Data)"
  ) +
  ylim(c(-1, 0.5)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )
ggsave("bench/img/ibs.png", width = 9, height = 5, dpi = 600, bg = "white")
