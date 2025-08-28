#' **Sparsity, stability, runtime, redundancy** analysis of feature selection methods
#' applied to multi-omic datasets.
#' For more details on the redundancy analysis, see `bench/redundancy.R`.
#'
#' This script visualizes and compares feature selection methods in terms of:
#' - **Sparsity**: Number of features selected per omic type
#' - **Stability**: Consistency of selected features across resampling iterations
#' - **Execution Time**: Computational time required by each method
#' - **Redundancy**: Degree of correlation or dependency among the selected features
#'
#' Execute: `Rscript bench/fs_plots.R` (from project root)

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(stringr)
  library(stabm)
  library(mlr3misc)
})

# hue_colors = scale_color_hue()$palette(n = 5)
# more manual hue colors:
# hue_colors = c("#F8766D", "#A3A500", "#00BB4E", "#E76BF3", "#35A2FF", "#9590FF")

# RColorBrewer::brewer.pal(n = 9, name = "Set1") # plus two more I added
set1_colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
                "#A65628", "#F781BF", "#999999", "#17BECF","#6A5ACD")

custom_colors = c(
  "hEFS (9 models)" = set1_colors[1],
  "EFS (CoxLasso)" = set1_colors[2],
  "hEFS (3 RSFs)" = set1_colors[3],
  "CoxLasso" = set1_colors[4]
)

# dataset labels
#osipov_abbrev = "Osipov et al. (2024)"
#wissel_abbrev = "Wissel et al. (2023)"
#cao_abbrev    = "Cao et al. (2021)"
osipov_abbrev = "MolTwin"
wissel_abbrev = "TCGA"
cao_abbrev    = "CPTAC"
dataset_labels = c(
  osipov2024 = osipov_abbrev,
  wissel2023 = wissel_abbrev,
  cao2021    = cao_abbrev
)

rename_omics = function(id) {
  case_when(
    id == "cnv" ~ "CNV",
    id == "snv" ~ "SNV",
    id == "indel" ~ "INDELs",
    id == "path" ~ "Pathology",
    id == "gex" ~ "GEX",
    id == "meth" ~ "Methylation",
    id == "mutation" ~ "Mutation",
    id == "rppa" ~ "RPPA",
    id == "prot" ~ "Proteomics",
    id == "phos" ~ "PhosphoP",
    id == "ngly" ~ "N-GlycoP",
    TRUE ~ "" # omic id that is not mapped yet
  )
}

# load the feature selection results per omic
fs = readRDS(file = "bench/fs.rds")

# Sparsity (Per-Omic) ----
cat("Sparsity")

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
p_sparse = fs_long |>
  group_by(dataset_id, omic_id) |>
  mutate(method = fct_reorder(method, n_feats, .fun = median, .desc = TRUE)) |>
  ungroup() |>
  mutate(omic_id = rename_omics(omic_id)) |>
  ggplot(aes(x = omic_id, y = n_feats, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  geom_jitter(aes(color = method), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 0.1) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 14) +
  facet_wrap(
    ~ dataset_id,
    scales = "free_x",  # One plot per dataset_id
    labeller = as_labeller(dataset_labels)
  ) +
  labs(
    x = "Omics",
    y = "Number of Selected Features",
    fill = "FS Method",
    title = "Feature Selection Sparsity across Omic Types"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    text = element_text(family = "Arial")
  )

p_sparse2 = ggdraw(p_sparse) +
  draw_text("Lower is better", x = 0.17, y = 0.6, angle = 90, size = 10, hjust = 0) +
    draw_line(
      x = c(0.145, 0.145),
      y = c(0.8, 0.55),
      arrow = arrow(length = unit(0.03, "npc")),  # arrowhead at end
      colour = "black",
      size = 0.8
    )
ggsave("bench/img/sparsity.png", plot = p_sparse2, width = 7, height = 5,
       dpi = 600, bg = "white")

# Same plot, without coxlasso fs
p_sparse_nocoxlasso = fs_long |>
  filter(method != "CoxLasso") |>
  group_by(dataset_id, omic_id) |>
  mutate(method = fct_reorder(method, n_feats, .fun = median, .desc = TRUE)) |>
  ungroup() |>
  mutate(omic_id = rename_omics(omic_id)) |>
  ggplot(aes(x = omic_id, y = n_feats, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  geom_jitter(aes(color = method), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 0.1) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 14) +
  facet_wrap(
    ~ dataset_id,
    scales = "free_x",
    labeller = as_labeller(dataset_labels)
  ) +
  labs(
    x = NULL, # Omics
    y = "Number of Selected Features",
    fill = "FS Method",
    #title = "Feature Selection Sparsity across Omic Types"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    text = element_text(family = "Arial"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )
ggsave("bench/img/sparsity_no_coxlasso.png", plot = p_sparse_nocoxlasso,
       width = 8, height = 5, dpi = 600, bg = "white")

# Stability (Per-Omic) ----
cat("Stability")

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

## Jaccard ----
method_lvls = c("hEFS (9 models)", "EFS (CoxLasso)", "hEFS (3 RSFs)", "CoxLasso")

# Assess stability across all 100 resamplings - 1 value per (dataset, omic) combo
stab_summary_jacc = fs_long |>
  group_by(dataset_id, omic_id, method) |>
  summarize(
    stability = stabm::stabilityJaccard(features),
    .groups = "drop"
  )

p_jacc = stab_summary_jacc |>
  #group_by(dataset_id, omic_id) |>
  #mutate(method = fct_reorder(method, stability, .fun = median, .desc = TRUE)) |>
  #ungroup() |>
  mutate(method = factor(method, levels = method_lvls)) |>
  mutate(omic_id = rename_omics(omic_id)) |>
  ggplot(aes(x = omic_id, y = stability, fill = method)) +
    geom_col(position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = custom_colors) +
    theme_minimal(base_size = 14) +
    facet_wrap(
      ~ dataset_id,
      scales = "free_x",
      labeller = as_labeller(dataset_labels)
    ) +
    ylim(c(0, 0.6)) +
    labs(
      x = "Omics",
      y = "Jaccard Similarity",
      fill = "FS Method",
      title = "Feature Selection Stability across Omic Types"
    ) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      text = element_text(family = "Arial")
    )

p_jacc2 =
  ggdraw(p_jacc) +
  draw_text("Higher is better", x = 0.14, y = 0.6, angle = 90, size = 10, hjust = 0) +
  draw_line(
    x = c(0.115, 0.115),
    y = c(0.57, 0.82),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black",
    size = 0.8
  )
ggsave("bench/img/stability_jaccard.png", plot = p_jacc2, width = 7, height = 5, dpi = 600, bg = "white")

## Nogueira ----
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

p_nog = stab_summary_nog |>
  mutate(method = factor(method, levels = method_lvls)) |>
  #group_by(dataset_id, omic_id) |>
  #mutate(method = fct_reorder(method, stability, .fun = median, .desc = TRUE)) |>
  #ungroup() |>
  mutate(omic_id = rename_omics(omic_id)) |>
  ggplot(aes(x = omic_id, y = stability, fill = method)) +
  geom_col(position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal(base_size = 14) +
  facet_wrap(
    ~ dataset_id,
    scales = "free_x",
    labeller = as_labeller(dataset_labels)
  ) +
  ylim(c(0, 0.62)) +
  labs(
    x = "Omics",
    y = "Nogueira Similarity",
    fill = "FS Method",
    title = "Feature Selection Stability across Omic Types"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    text = element_text(family = "Arial")
  )

p_nog2 =
  ggdraw(p_nog) +
  draw_text("Higher is better", x = 0.14, y = 0.65, angle = 90, size = 10, hjust = 0) +
  draw_line(
    x = c(0.115, 0.115),
    y = c(0.6, 0.85),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black",
    size = 0.8
  )
ggsave("bench/img/stability_nogueira.png", plot = p_nog2, width = 8, height = 5,
       dpi = 600, bg = "white")

## Phi ----
# this stability measure is the average pair-wise Pearson correlation between
# binary vectors of length p, where features get 1 if selected, otherwise 0.
stab_summary_phi = fs_long |>
  left_join(p_lookup, by = c("dataset_id", "omic_id")) |>
  group_by(dataset_id, omic_id, method) |>
  summarize(
    stability = stabm::stabilityPhi(features, p = p[1]),
    .groups = "drop"
  )

p_phi = stab_summary_phi |>
  mutate(method = factor(method, levels = method_lvls)) |>
  mutate(omic_id = rename_omics(omic_id)) |>
  ggplot(aes(x = omic_id, y = stability, fill = method)) +
  geom_col(position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal(base_size = 14) +
  facet_wrap(
    ~ dataset_id,
    scales = "free_x",
    labeller = as_labeller(dataset_labels)
  ) +
  ylim(c(0, 0.62)) +
  labs(
    x = "Omics",
    y = "Phi Similarity",
    fill = "FS Method",
    title = "Feature Selection Stability across Omic Types"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    text = element_text(family = "Arial")
  )
ggsave("bench/img/stability_phi.png", plot = p_phi, width = 8, height = 5,
       dpi = 600, bg = "white")

## Jaccard (Variability) ----
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

p_jacc_rsmp = stab_jacc_rsmp |>
  mutate(method = factor(method, levels = method_lvls)) |>
  mutate(omic_id = rename_omics(omic_id)) |>
  ggplot(aes(x = omic_id, y = stability, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  geom_jitter(aes(color = method), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75),
              alpha = 0.5, size = 0.01) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 14) +
  facet_wrap(
    ~ dataset_id,
    scales = "free_x",
    labeller = as_labeller(dataset_labels)
  ) +
  labs(
    x = NULL, # Omics
    y = "Jaccard Similarity",
    fill = "FS Method",
    #title = "Feature Selection Stability across Omic Types"
  ) +
  ylim(c(0, 0.6)) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    text = element_text(family = "Arial"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )
ggsave("bench/img/stability_jaccard_var.png", plot = p_jacc_rsmp,
       width = 8, height = 5, dpi = 600, bg = "white")

## Noqueira (Variability) ----
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

p_nog_rsmp = stab_nog_rsmp |>
  mutate(method = factor(method, levels = method_lvls)) |>
  mutate(omic_id = rename_omics(omic_id)) |>
  ggplot(aes(x = omic_id, y = stability, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  geom_jitter(aes(color = method), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75),
              alpha = 0.5, size = 0.01) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 14) +
  facet_wrap(
    ~ dataset_id,
    scales = "free_x",
    labeller = as_labeller(dataset_labels)
  ) +
  labs(
    x = NULL, # Omics
    y = "Noqueira Similarity",
    fill = "FS Method",
    #title = "Feature Selection Stability across Omic Types"
  ) +
  ylim(c(0, 0.65)) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    text = element_text(family = "Arial"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )
ggsave("bench/img/stability_nogueira_var.png", plot = p_nog_rsmp, width = 8,
       height = 5, dpi = 600, bg = "white")

# Runtime ----
cat("Runtime")

#' Compare Coxlasso vs hEFS variants
timings = fs |>
  select(dataset_id, omic_id, rsmp_id, coxlasso_train_time) |>
  rename(coxlasso = coxlasso_train_time)
summary(timings$coxlasso) # 0.5 - 6.3 secs => CoxLasso is super fast

# function to get the efs timings from the stored .csv files
read_efs_times = function(dataset_id, omic_id, rsmp_id) {
  file_path = file.path("bench", "efs", dataset_id, omic_id, paste0("times_", rsmp_id, ".csv"))
  if (!file.exists(file_path)) {
    return(NULL)
  }

  df = readr::read_csv(file_path, show_col_types = FALSE, progress = FALSE) # (id, time), see `efs.R`
  has_aorsf = "aorsf" %in% df$id

  df |>
    bind_rows(
      tibble::tibble(id = "all", time = sum(df$time, na.rm = TRUE))
    ) |>
    pivot_wider(names_from = "id", values_from = "time") |>
    rename_with(~ paste0("efs_", .)) |>
    mutate(
      # the 3 RSF models aggregated time (`efs_rsf` is RSF logrank + maxstat)
      efs_rsf_total = ifelse(has_aorsf, efs_rsf + efs_aorsf, efs_rsf),
      # the RSF using `logrank` splitrule is ~3x slower that using `maxstat` splitrule
      # based on some previous benchmarks I had done => this is an estimation ofc
      # as I ran these two together in this benchmark (see `efs.R`), but it is good enough!
      efs_rsf_logrank = efs_rsf * (3/4),
      efs_rsf_maxstat = efs_rsf * (1/4),
    ) |>
    select(!efs_rsf)
}

# add the efs execution times
## per learned id (EFS) and aggregated across some learners (hEFS)
timings = timings |>
  rowwise() |>
  mutate(read_efs_times(dataset_id, omic_id, rsmp_id)) |>
  ungroup()
# timings = readRDS("bench/timings.rds") # stored, to avoid running the above

timings_long = timings |>
  pivot_longer(
    cols = starts_with("efs") | starts_with("coxlasso"),
    names_to = "method",
    values_to = "time"
  ) |>
  mutate(method = case_when(
    method == "coxlasso" ~ "CoxLasso",
    method == "efs_rsf_logrank" ~ "RSF-logrank",
    method == "efs_rsf_maxstat" ~ "RSF-maxstat",
    method == "efs_aorsf" ~ "AORSF",
    method == "efs_rsf_total" ~ "hEFS (3 RSFs)",
    method == "efs_xgb_cox" ~ "XGBoost-Cox",
    method == "efs_xgb_aft_log" ~ "XGBoost-AFT",
    method == "efs_glmb_cox" ~ "GLMBoost-Cox",
    method == "efs_glmb_loglog" ~ "GLMBoost-AFT",
    method == "efs_coxboost" ~ "CoxBoost",
    method == "efs_coxlasso" ~ "EFS (CoxLasso)",
    method == "efs_all" ~ "hEFS (9 models)",
    TRUE ~ method  # there shouldn't be other categories here
  ))

fs_method_ids = c("CoxLasso", "EFS (CoxLasso)", "hEFS (3 RSFs)", "hEFS (9 models)")
timing_stats = timings_long |>
  group_by(method) |> #dataset_id, method) |> #omic_id, method) |>
  summarise(
    mean_time = mean(time, na.rm = TRUE),
    sd_time   = sd(time, na.rm = TRUE)
    #.groups = "drop"
  ) |>
  mutate(summary = sprintf("%.2f ± %.2f", mean_time, sd_time))
timing_stats |>
  arrange(desc(mean_time)) #|> filter(method %in% fs_method_ids)

# the 4 fs methods as in the previous per-omic plots
p_times = timings_long |>
  filter(method %in% fs_method_ids) |>
  mutate(time = time/60) |> # convert to min
  mutate(method = fct_reorder(method, time, .fun = median, .desc = TRUE)) |>
  mutate(omic_id = rename_omics(omic_id)) |>
  ggplot(aes(x = omic_id, y = time, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  geom_jitter(aes(color = method), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75),
              alpha = 0.3, size = 0.01) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 14) +
  facet_wrap(
    ~ dataset_id,
    scales = "free_x",
    labeller = as_labeller(dataset_labels)
  ) +
  labs(
    x = NULL, # Omics
    y = "Execution Time (minutes)",
    fill = "FS Method",
    #title = "Execution Time of FS Methods across Omic Types"
  ) +
  ylim(c(0, 15)) + # 15 min upper bound
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    text = element_text(family = "Arial"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )
ggsave("bench/img/timings.png", plot = p_times, width = 8, height = 5,
       dpi = 600, bg = "white")

efs_colors = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
               "#1D91C0", "#BABABA", "#000000")

# remove the hybrid methods (that combine multiple learners) + simple CoxLasso,
# so leaving only the homogeneous ensemble methods, i.e. method + 100 subsamplings
p_times2 = timings_long |>
  filter(!method %in% c("CoxLasso", "hEFS (3 RSFs)", "hEFS (9 models)")) |>
  # remove the EFS from CoxLasso, here we keep only the model names
  mutate(method = case_match(method, "EFS (CoxLasso)" ~ "CoxLasso", .default = method)) |>
  mutate(method = fct_reorder(method, time, .fun = median, .na_rm = TRUE, .desc = TRUE)) |>
  mutate(time = time/60) |> # convert to min
  mutate(omic_id = rename_omics(omic_id)) |>
  ggplot(aes(x = omic_id, y = time, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  geom_jitter(aes(color = method), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75),
              alpha = 0.3, size = 0.01) +
  scale_fill_manual(values = efs_colors) +
  scale_color_manual(values = efs_colors) +
  theme_minimal(base_size = 14) +
  facet_wrap(
    ~ dataset_id,
    scales = "free_x",
    labeller = as_labeller(dataset_labels)
  ) +
  labs(
    x = NULL, # Omics
    y = "Execution Time (minutes)",
    fill = "Model",
    #title = "Execution Time of FS Methods across Omic Types"
  ) +
  ylim(0, 5) + # 5 min upper bound
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    text = element_text(family = "Arial"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )
ggsave("bench/img/timings_per_learner.png", plot = p_times2, width = 8, height = 5,
       dpi = 600, bg = "white")

# Redundancy (Per-Omic) ----
res = readRDS("bench/fs_red.rds")

res_long = res |>
  pivot_longer(
    cols = c(rr_pearson, srp_pearson, rr_spearman, srp_spearman, rr_xicor, srp_xicor), # arr_xicor
    names_to = c("type", "metric"),
    names_pattern = "(rr|srp|arr)_(pearson|spearman|xicor)",
    values_to = "value"
  ) |>
  rename(method = fs_method_id) |>
  mutate(method = case_when(
    method == "coxlasso" ~ "CoxLasso",
    method == "efs_rsf" ~ "hEFS (3 RSFs)",
    method == "efs_coxlasso" ~ "EFS (CoxLasso)",
    method == "efs_all" ~ "hEFS (9 models)",
    TRUE ~ method  # there shouldn't be other categories here
  )) |>
  mutate(omic_id = rename_omics(omic_id))

metric_labels = c(
  pearson = "Pearson",
  spearman = "Spearman",
  xicor = "ξ Correlation"
)

method_lvls = c("CoxLasso", "EFS (CoxLasso)", "hEFS (3 RSFs)", "hEFS (9 models)")

## Redundancy Rate ----
p_rrate_all = res_long |>
  filter(type == "rr", !is.na(value)) |> # all metrics
  mutate(method = factor(method, levels = method_lvls)) |>
  # Reorder method based on overall median redundancy
  #mutate(method = fct_reorder(method, value, .fun = median, .na_rm = TRUE, .desc = FALSE)) |>
  ggplot(aes(x = omic_id, y = value, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  geom_jitter(aes(color = method), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 0.1) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 14) +
  facet_grid(
    metric ~ dataset_id,
    scales = "free_x",
    labeller = labeller(
      metric = as_labeller(metric_labels),
      dataset_id = as_labeller(dataset_labels))
  ) +
  labs(
    x = NULL, #Omics
    y = "Redundancy Rate",
    fill = "FS Method",
    #title = "Redundancy across Omic Types and Metrics"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    text = element_text(family = "Arial"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )
ggsave("bench/img/redundancy_rate_all_metrics.png", plot = p_rrate_all,
       width = 10, height = 8, dpi = 600, bg = "white")

p_rrate_xicor = res_long |>
  filter(type == "rr", metric == "xicor", !is.na(value)) |>
  mutate(method = factor(method, levels = method_lvls)) |>
  ggplot(aes(x = omic_id, y = value, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  geom_jitter(aes(color = method), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 0.1) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 14) +
  facet_grid(
    ~ dataset_id,
    scales = "free_x",
    labeller = labeller(dataset_id = as_labeller(dataset_labels))
  ) +
  labs(
    x = NULL, # Omics
    y = "Redundancy Rate", # (ξ Correlation)
    fill = "FS Method",
    #title = "Redundancy across Omic Types"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    text = element_text(family = "Arial"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )
ggsave("bench/img/redundancy_rate_xicor.png", plot = p_rrate_xicor,
       width = 10, height = 5, dpi = 600, bg = "white")

## Significant Redundancy Proportion ----
p_srp_all = res_long |>
  filter(type == "srp", !is.na(value)) |> # all metrics
  mutate(method = factor(method, levels = method_lvls)) |>
  #mutate(method = fct_reorder(method, value, .fun = median, .na_rm = TRUE, .desc = FALSE)) |>
  ggplot(aes(x = omic_id, y = value, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  geom_jitter(aes(color = method), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 0.1) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 14) +
  facet_grid(
    metric ~ dataset_id,
    scales = "free_x",
    labeller = labeller(
      metric = as_labeller(metric_labels),
      dataset_id = as_labeller(dataset_labels))
  ) +
  labs(
    x = NULL, # Omics
    y = "Significant Redundancy Proportion",
    fill = "FS Method",
    #title = "Redundancy across Omic Types and Metrics"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    text = element_text(family = "Arial"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )
ggsave("bench/img/srp5_all_metrics.png", plot = p_srp_all, width = 10, height = 8,
       dpi = 600, bg = "white") # srp5 => due to `alpha = 0.05`

p_srp5_xicor = res_long |>
  filter(type == "srp", metric == "xicor", !is.na(value)) |>
  mutate(method = factor(method, levels = method_lvls)) |>
  ggplot(aes(x = omic_id, y = value, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  geom_jitter(aes(color = method), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 0.1) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 14) +
  facet_grid(
    ~ dataset_id,
    scales = "free_x",
    labeller = labeller(dataset_id = as_labeller(dataset_labels))
  ) +
  labs(
    x = NULL, # Omics
    y = "Significant Redundancy Proportion", #\n(ξ Correlation)",
    fill = "FS Method",
    #title = "Redundancy across Omic Types"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    text = element_text(family = "Arial"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )
ggsave("bench/img/srp5_xicor.png", plot = p_srp5_xicor, width = 10, height = 5,
       dpi = 600, bg = "white")

## Adjusted Redundancy Rate ----
if (FALSE) {
# there are NA's in all scores due to feature sizes of length 1, but the following
# are extra and are due to null distribution/expected score being `NaN`
res |> filter(!is.na(rr_pearson), is.na(arr_xicor))

res_long1 = res_long |>
  filter(type == "arr") |> # only "xicor" metric here
  filter(!is.na(value), !is.infinite(value)) |>
  mutate(method = factor(method, levels = method_lvls))

p_arr = res_long1 |>
  ggplot(aes(x = omic_id, y = value, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  geom_jitter(aes(color = method), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 0.1) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 14) +
  facet_grid(
    ~ dataset_id,
    scales = "free_x",
    labeller = labeller(dataset_id = as_labeller(dataset_labels))
  ) +
  labs(
    x = "Omics",
    y = "Adjusted Redundancy Rate",
    fill = "FS Method",
    title = "Redundancy across Omic Types and Metrics"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    text = element_text(family = "Arial")
  )
ggsave("bench/img/redundancy_adjusted_xicor.png", plot = p_arr, width = 7, height = 5,
       dpi = 600, bg = "white")

p_arr2 = res_long1 |>
  filter(!omic_id == "mutation") |> # due to very different scale of values
  ggplot(aes(x = omic_id, y = value, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge2(preserve = "single")) +
  geom_jitter(aes(color = method), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 0.1) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 14) +
  facet_grid(
    ~ dataset_id,
    scales = "free_x",
    labeller = labeller(dataset_id = as_labeller(dataset_labels))
  ) +
  labs(
    x = "Omics",
    y = "Adjusted Redundancy Rate",
    fill = "FS Method",
    title = "Redundancy across Omic Types and Metrics"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    text = element_text(family = "Arial")
  )
ggsave("bench/img/redundancy_adjusted_xicor_no_mut.png", plot = p_arr2,
       width = 7, height = 5, dpi = 600, bg = "white")
}
