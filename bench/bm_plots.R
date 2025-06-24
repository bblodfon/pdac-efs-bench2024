suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(stringr)
})

# Benchmark results ----
result = readRDS(file = "bench/result.rds")

# dataset labels
osipov_abbrev = "Osipov et al. (2024)"
wissel_abbrev = "Wissel et al. (2023)"
dataset_labels = c(osipov2024 = osipov_abbrev, wissel2023 = wissel_abbrev)

# colors
# RColorBrewer::brewer.pal(n = 9, name = "Set1") # plus two more I added
set1_colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
                "#A65628", "#F781BF", "#999999", "#17BECF","#6A5ACD")

custom_colors = c(
  "hEFS (9 models)" = set1_colors[1],
  "EFS (CoxLasso)" = set1_colors[2],
  "hEFS (3 RSFs)" = set1_colors[3],
  "CoxLasso" = set1_colors[4]
)

# rename data
result = result |>
  # rename model-data-configs
  rename(model = model_data_config) |>
  mutate(model = case_when(
    model %in% c("coxlasso-all", "coxlasso-all-none") ~ "CoxLasso (ALL)",
    model %in% c("rsf-all", "rsf-all-none") ~ "RSF (ALL)",
    model == "cox-clinical" ~ "Cox (Clinical)",
    model == "rsf-clinical" ~ "RSF (Clinical)",
    model %in% c("blockforest-all", "blockforest-all-none") ~ "BlockForest (ALL)",
    model == "rsf-clinical+gex-efs_all_feats" ~ "RSF (Clinical + GEX)",
    model == "rsf-gex-efs_all_feats" ~ "RSF (GEX)",
    model == "coxlasso-clinical+gex-coxlasso_feats" ~ "CoxLasso (Clinical + GEX)",
    model == "blockforest-clinical+gex-efs_all_feats" ~ "BlockForest (Clinical + GEX)",
    TRUE ~ model
  )) |>
  # rename FS methods
  mutate(fs_method_id = case_when(
    fs_method_id == "coxlasso_feats" ~ "CoxLasso",
    fs_method_id == "efs_all_feats" ~ "hEFS (9 models)",
    fs_method_id == "efs_coxlasso_feats" ~ "EFS (CoxLasso)",
    fs_method_id == "efs_rsf_feats" ~ "hEFS (3 RSFs)",
    # `NA` fs method means that the 4 above FS methods were not used at all
    # distinguish clinical as Baselines from other no-fs
    is.na(fs_method_id) & model %in% c("Cox (Clinical)", "RSF (Clinical)")  ~ "Baseline",
    is.na(fs_method_id) ~ "no FS",
    TRUE ~ fs_method_id  # Keep other values unchanged
  ))

# check data
result |> group_by(dataset_id, model, fs_method_id) |> count() |> print(n = 100)

# Convert data to long format for ggplot
result_long = result |>
  pivot_longer(cols = c(harrell_c, uno_c, brier_tmax24),
               names_to = "measure", values_to = "value")

# useful for filtering later - without "no FS" or "Baselines"
fs_method_ids = c("CoxLasso", "hEFS (9 models)", "hEFS (3 RSFs)", "EFS (CoxLasso)")

# MM Sparsity ----
# MM => Multi-Modal
p_mmsparse = result |>
  # keep one model config that has multi-omics fs (ALL), as we are interested
  # only in number of features after late integration
  filter(model == "CoxLasso (ALL)", fs_method_id != "no FS") |>
  mutate(
    fs_method_id = fct_reorder(fs_method_id, task_nfeats, .fun = median, .desc = TRUE)
  ) |>
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
      labeller = as_labeller(dataset_labels)
    ) +
    theme(
      axis.text.x = element_blank(),
      # legend.position = "bottom",
      # legend.box = "horizontal",
      # legend.margin = margin(t = -20),
      text = element_text(family = "Arial")
    ) +
    #guides(fill = guide_legend(nrow = 2)) +
    labs(
      x = "",
      y = "Number of Total Features",
      title = "Sparsity of Multi-omics Datasets",
      fill = "FS Method"
    )

p_mmsparse2 = ggdraw(p_mmsparse) +
  draw_text("Lower is better", x = 0.17, y = 0.6, angle = 90, size = 10, hjust = 0) +
  draw_line(
    x = c(0.145, 0.145),
    y = c(0.8, 0.55),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black",
    size = 0.8
  )
ggsave("bench/img/mm_sparsity.png", plot = p_mmsparse2, width = 7, height = 5,
       dpi = 600, bg = "white")

# MM contribution ----
# MM => Multi-Modal
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
  # keep one model config that has multi-omics fs (ALL), as we are interested
  # only in number of features after late integration
  filter(model == "CoxLasso (ALL)", fs_method_id != "no FS") |>
  select(dataset_id, fs_method_id, rsmp_id, task_nfeats, task_feats) |>
  tidyr::unnest(task_feats) |> # smart, every row is per feature
  mutate(omics_type = extract_omics(task_feats)) |>
  group_by(dataset_id, fs_method_id, rsmp_id, omics_type, task_nfeats) |>
  summarise(count = n(), .groups = "drop") |>
  mutate(percentage = count / task_nfeats) |>
  group_by(dataset_id, fs_method_id, omics_type) |>
  summarise(avg_percentage = mean(percentage), avg_count = round(mean(count)), .groups = "drop")

# Stacked percentage plot
p_omic_per = ctr_res |>
  # ggplot(aes(x = reorder(fs_method_id, -avg_percentage, sum),
  #            y = avg_percentage, fill = omics_type)) +
  mutate(fs_method_id = factor(fs_method_id, levels = fs_method_ids)) |>
  ggplot(aes(x = fs_method_id, y = 100*avg_percentage, fill = omics_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = set1_colors) +
  facet_wrap(~ dataset_id, labeller = as_labeller(dataset_labels)) +
  labs(
    x = "FS Method",
    y = "Average Contribution",
    fill = "Data Type"
  ) +
  theme_minimal(base_size = 14) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "top",
    text = element_text(family = "Arial")
  )
ggsave("bench/img/omic_contribution_stacked.png", plot = p_omic_per,
       width = 9, height = 5, dpi = 600, bg = "white")

# Average counts barplot
## order by average total number of features per omic
if (FALSE) {
  ctr_res2 = ctr_res |>
    filter(omics_type != "Clinical") |>
    group_by(dataset_id, omics_type) |>
    summarize(avg_total = sum(avg_count), .groups = "drop") |>
    arrange(dataset_id, desc(avg_total)) |>
    mutate(omics_type = factor(omics_type, levels = unique(omics_type))) |>
    select(-avg_total) |>
    left_join(ctr_res, by = c("dataset_id", "omics_type")) |>
    mutate(omics_type = factor(omics_type, levels = unique(omics_type)))

  p_omic_counts = ctr_res2 |>
    ggplot(aes(x = reorder(fs_method_id, -avg_count, sum),
               y = avg_count, fill = omics_type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = set1_colors) +
    facet_wrap(
      ~ dataset_id,
      scales = "free_x",
      labeller = as_labeller(dataset_labels)
    ) +
    labs(
      x = "Feature Selection Method",
      y = "Average Number of Features",
      fill = "Omic Type"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      text = element_text(family = "Arial")
    )
  ggsave("bench/img/omic_contribution_counts.png", plot = p_omic_counts,
         width = 8, height = 5, dpi = 600, bg = "white")
}

# Predictivity (C-index) ----
# some stats for Harrell's C-index per dataset
## Wissel
result_long |>
  filter(dataset_id == "wissel2023", measure == "harrell_c", grepl("^(RSF|BlockF)", model)) |>
  group_by(dataset_id, model, fs_method_id) |>
  summarize(cindex = median(value), .groups = "drop") |>
  arrange(desc(cindex))

result_long |>
  filter(dataset_id == "wissel2023", measure == "harrell_c", grepl("^Cox", model)) |>
  group_by(dataset_id, model, fs_method_id) |>
  summarize(cindex = median(value), .groups = "drop") |>
  arrange(desc(cindex))

## Osipov
result_long |>
  filter(dataset_id == "osipov2024", measure == "harrell_c", grepl("^(RSF|BlockF)", model)) |>
  group_by(dataset_id, model, fs_method_id) |>
  summarize(cindex = median(value), .groups = "drop") |>
  arrange(desc(cindex))

result_long |>
  filter(dataset_id == "osipov2024", measure == "harrell_c", grepl("^Cox", model)) |>
  group_by(dataset_id, model, fs_method_id) |>
  summarize(cindex = median(value), .groups = "drop") |>
  arrange(desc(cindex))

measures = c("harrell_c", "uno_c") # note: Uno's results same as Harrell's C
fs_colors = c(custom_colors, "Baseline" = set1_colors[9], "no FS" = set1_colors[5])

## CoxLasso (integration model) ----
for (meas in measures) {
  for (keep_nofs in c(TRUE, FALSE)) {
    res_cox =
      result_long |>
      filter(measure == !!meas) |>
      filter(model %in% c("CoxLasso (ALL)", "Cox (Clinical)")) |>
      (\(df) if (keep_nofs) df else filter(df, fs_method_id != "no FS"))() |>
      mutate(
        fs_method_id = factor(fs_method_id, levels = c("Baseline", fs_method_ids, "no FS"))
      )

    p_cox =
      res_cox |>
      ggplot(aes(x = fs_method_id, y = value, fill = fs_method_id)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(aes(color = fs_method_id), show.legend = FALSE,
                  position = position_jitterdodge(jitter.width = 1.2, dodge.width = 0.75),
                  alpha = 0.7, size = 0.01) +
      geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", linewidth = 0.8) +
      facet_wrap(
        ~ dataset_id,
        labeller = as_labeller(dataset_labels)
      ) +
      scale_fill_manual(
        values = fs_colors,
        breaks = c(fs_method_ids, "no FS")
      ) +
      scale_color_manual(
        values = fs_colors,
        breaks = c(fs_method_ids, "no FS")
      ) +
      ylim(c(0, 1)) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.margin = margin(t = -20),
        text = element_text(family = "Arial")
      ) +
      labs(
        x = "",
        y = paste0(ifelse(meas == "harrell_c", "Harrell's", "Uno's"), " C-index"),
        #title = "CoxLasso Integration vs Cox Clinical Baseline",
        fill = "FS Method"
      ) +
      geom_vline(xintercept = 1.5, linetype = "dotted", color = "black", linewidth = 0.8) +
      annotate("text", x = 1, y = 0.95, label = "CoxPH\n(Clinical)",
               size = 4, hjust = 0.5, family = "Arial") +
      annotate("text", x = 3.5, y = 0.95, label = "CoxLasso\n(Clinical + Multi-omics)",
               size = 4, hjust = 0.5, family = "Arial")

    p_cox2 = ggdraw(p_cox) +
      draw_text("Higher is better", x = 0.117, y = 0.15, angle = 90, size = 10, hjust = 0) +
      draw_line(
        x = c(0.1, 0.1),
        y = c(0.15, 0.4),
        arrow = arrow(length = unit(0.03, "npc")),
        colour = "black",
        size = 0.8
      )
    ggsave(paste0("bench/img/", meas, "_coxlasso", ifelse(keep_nofs, "_nofs", ""), ".png"),
           plot = p_cox2, width = 9, height = 5, dpi = 600, bg = "white")
  }
}

## RSF (integration model) ----
for (meas in measures) {
  for (keep_nofs in c(TRUE, FALSE)) {
    res_rsf =
      result_long |>
      filter(measure == !!meas) |>
      filter(model %in% c("RSF (ALL)", "RSF (Clinical)")) |>
      (\(df) if (keep_nofs) df else filter(df, fs_method_id != "no FS"))() |>
      mutate(
        fs_method_id = factor(fs_method_id, levels = c("Baseline", fs_method_ids, "no FS"))
      )

    p_rsf =
      res_rsf |>
      ggplot(aes(x = fs_method_id, y = value, fill = fs_method_id)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(aes(color = fs_method_id), show.legend = FALSE,
                  position = position_jitterdodge(jitter.width = 1.2, dodge.width = 0.75),
                  alpha = 0.7, size = 0.01) +
      geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", linewidth = 0.8) +
      facet_wrap(
        ~ dataset_id,
        labeller = as_labeller(dataset_labels)
      ) +
      scale_color_manual(
        values = fs_colors,
        breaks = c(fs_method_ids, "no FS")
      ) +
      scale_fill_manual(
        values = fs_colors,
        breaks = c(fs_method_ids, "no FS")
      ) +
      ylim(c(0, 1)) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.margin = margin(t = -20),
        text = element_text(family = "Arial")
      ) +
      labs(
        x = "",
        y = paste0(ifelse(meas == "harrell_c", "Harrell's", "Uno's"), " C-index"),
        #title = "RSF Integration vs RSF Clinical Baseline",
        fill = "FS Method"
      ) +
      geom_vline(xintercept = 1.5, linetype = "dotted", color = "black", linewidth = 0.8) +
      annotate("text", x = 1, y = 0.95, label = "RSF\n(Clinical)",
               size = 4, hjust = 0.5, family = "Arial") +
      annotate("text", x = 3.5, y = 0.95, label = "RSF\n(Clinical + Multi-omics)",
               size = 4, hjust = 0.5, family = "Arial")

    p_rsf2 = ggdraw(p_rsf) +
      draw_text("Higher is better", x = 0.117, y = 0.15, angle = 90, size = 10, hjust = 0) +
      draw_line(
        x = c(0.1, 0.1),
        y = c(0.15, 0.4),
        arrow = arrow(length = unit(0.03, "npc")),
        colour = "black",
        size = 0.8
      )
    ggsave(paste0("bench/img/", meas, "_rsf", ifelse(keep_nofs, "_nofs", ""), ".png"),
           plot = p_rsf2, width = 9, height = 5, dpi = 600, bg = "white")
  }
}

## BlockForest (integration model) ----
for (meas in measures) {
  for (keep_nofs in c(TRUE, FALSE)) {
    res_bf =
      result_long |>
      filter(measure == !!meas) |>
      filter(model %in% c("BlockForest (ALL)", "RSF (Clinical)")) |>
      (\(df) if (keep_nofs) df else filter(df, fs_method_id != "no FS"))() |>
      mutate(
        fs_method_id = factor(fs_method_id, levels = c("Baseline", fs_method_ids, "no FS"))
      )

    p_bf =
      res_bf |>
      ggplot(aes(x = fs_method_id, y = value, fill = fs_method_id)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(aes(color = fs_method_id), show.legend = FALSE,
                  position = position_jitterdodge(jitter.width = 1.2, dodge.width = 0.75),
                  alpha = 0.7, size = 0.01) +
      geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", linewidth = 0.8) +
      facet_wrap(
        ~ dataset_id,
        labeller = as_labeller(dataset_labels)
      ) +
      scale_color_manual(
        values = fs_colors,
        breaks = c(fs_method_ids, "no FS")
      ) +
      scale_fill_manual(
        values = fs_colors,
        breaks = c(fs_method_ids, "no FS")
      ) +
      ylim(c(0, 1)) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.margin = margin(t = -20),
        text = element_text(family = "Arial")
      ) +
      labs(
        x = "",
        y = paste0(ifelse(meas == "harrell_c", "Harrell's", "Uno's"), " C-index"),
        #title = "BlockForest Integration vs RSF Clinical Baseline",
        fill = "FS Method"
      ) +
      geom_vline(xintercept = 1.5, linetype = "dotted", color = "black", linewidth = 0.8) +
      annotate("text", x = 1, y = 0.95, label = "RSF\n(Clinical)",
               size = 4, hjust = 0.5, family = "Arial") +
      annotate("text", x = 3.5, y = 0.95, label = "BlockForest\n(Clinical + Multi-omics)",
               size = 4, hjust = 0.5, family = "Arial")

    p_bf2 = ggdraw(p_bf) +
      draw_text("Higher is better", x = 0.117, y = 0.15, angle = 90, size = 10, hjust = 0) +
      draw_line(
        x = c(0.1, 0.1),
        y = c(0.15, 0.4),
        arrow = arrow(length = unit(0.03, "npc")),
        colour = "black",
        size = 0.8
      )
    ggsave(paste0("bench/img/", meas, "_bf", ifelse(keep_nofs, "_nofs", ""), ".png"),
           plot = p_bf2, width = 9, height = 5, dpi = 600, bg = "white")
  }
}

## Baselines only ----
for (meas in measures) {
  p_base = result_long |>
    filter(measure == meas, fs_method_id == "Baseline") |>
    mutate(model = case_when(
      model == "Cox (Clinical)" ~ "CoxPH",
      model == "RSF (Clinical)" ~ "RSF",
      TRUE ~ model # there shouldn't be any other category here
    )) |>
    ggplot(aes(x = model, y = value, fill = model)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(aes(color = model), width = 0.2, alpha = 0.5, size = 0.7, show.legend = FALSE) +
      geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", linewidth = 0.8) +
      facet_wrap(
        ~ dataset_id,
        labeller = as_labeller(dataset_labels)
      ) +
      scale_color_brewer(palette = "Set1") +
      scale_fill_brewer(palette = "Set1") +
      ylim(c(0, 1)) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_blank(),
        text = element_text(family = "Arial")
      ) +
      labs(
        x = "",
        y = paste0(ifelse(meas == "harrell_c", "Harrell's", "Uno's"), " C-index"),
        title = "Compare baseline models (only clinical data)",
        fill = "Baseline Model"
      )
  ggsave(paste0("bench/img/", meas, "_baselines.png"), plot = p_base,
         width = 7, height = 5, dpi = 600, bg = "white")
}

## ALL vs GEX+Clinical vs Clinical (Wissel/TCGA) ----
included_models = c(
  "RSF (ALL)", "RSF (Clinical)", "RSF (GEX)", "RSF (Clinical + GEX)",
  "BlockForest (ALL)",  "BlockForest (Clinical + GEX)"
)

p_cmps = result_long |>
  # Harrell's C-index and Wissel dataset only
  filter(dataset_id == "wissel2023", measure == "harrell_c") |>
  filter(model %in% included_models) |>
  filter(fs_method_id %in% c("hEFS (9 models)", "Baseline")) |>
  mutate(model = fct_reorder(model, value, .fun = median, .desc = TRUE)) |>
  ggplot(aes(x = model, y = value, fill = model)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = model), width = 0.2, alpha = 0.5, size = 0.7, show.legend = FALSE) +
    geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", linewidth = 0.8) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    ylim(c(0, 1)) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_blank(),
      text = element_text(family = "Arial")
    ) +
    labs(
      x = "",
      y = "Harrell's C-index",
      title = "Single-omic vs Multi-omics",
      subtitle = "hEFS applied for GEX and ALL omics",
      fill = "Model"
    )
ggsave("bench/img/harrell_c_cmps.png", plot = p_cmps, width = 7, height = 5,
       dpi = 600, bg = "white")

# Predictivity (ISBS) ----
p_isbs = result_long |>
  filter(measure == "brier_tmax24") |>
  ggplot(aes(x = fs_method_id, y = value, fill = model)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = model), show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 0.01) +
  scale_fill_manual(values = set1_colors) +
  scale_color_manual(values = set1_colors) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.8) +
  theme_minimal(base_size = 14) +
  facet_grid(
    measure ~ dataset_id,
    labeller = as_labeller(
      c(brier_tmax24 = "Up to t = 2 years",
        osipov2024 = "Osipov et al. (2024)",
        wissel2023 = "Wissel et al. (2023)"
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
ggsave("bench/img/isbs.png", plot = p_isbs, width = 11, height = 6, dpi = 600, bg = "white")
