#' here we combine figures from `fs_plots.R` and `bm_plots.R` scripts
#' so these must be ran beforehand so that the plots are available here
suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

# Figure 3 ----
p_sparse3 = p_sparse +
  labs(title = NULL, x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.7),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    text = element_text(family = "Arial"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(face = "bold", size = 14),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )

p_mmsparse3 = p_mmsparse +
  labs(title = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    legend.position = "none",
    text = element_text(family = "Arial"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(face = "bold", size = 14),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )

p_omic_per3 = p_omic_per +
  labs(x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "top",
    text = element_text(family = "Arial"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )

top_row = plot_grid(p_sparse3, labels = "a")
bottom_row = plot_grid(
  p_mmsparse3,
  p_omic_per3,
  rel_widths = c(1.5, 2),
  labels = c("b","c"),
  ncol = 2
)

# Combine the two rows with two plots above and one below
fig3 = plot_grid(
  top_row,
  bottom_row,
  nrow = 2,
  rel_heights = c(1, 1),
  label_fontfamily = "Arial"
)
ggsave("bench/img/fig3.png", plot = fig3, width = 13, height = 10, dpi = 600, bg = "white")
ggsave("bench/img/fig3.pdf", plot = fig3, width = 13, height = 10, dpi = 600, bg = "white", device = cairo_pdf)

# Sup. Fig. 4 ----
# Load ensemble feature selection (EFS) results
library(mlr3)
library(mlr3fselect)
library(mlr3viz)

# 'good' example
file_name = file.path("bench", "efs", "wissel2023", "cnv", "efs_7.rds")
efs1 = readRDS(file_name)
kp1 = efs1$knee_points(type = "estimated") # 15

res1 = efs1$.__enclos_env__$private$.result

# Rename learner_id
rename_map = c(
  "rsf_logrank.fselector" = "RSF-logrank",
  "rsf_maxstat.fselector" = "RSF-maxstat",
  "aorsf.fselector" = "AORSF",
  "xgb_cox.fselector" = "XGBoost-Cox",
  "xgb_aft_log.fselector" = "XGBoost-AFT",
  "glmb_cox" = "GLMBoost-Cox",
  "glmb_loglog" = "GLMBoost-AFT",
  "coxboost" = "CoxBoost",
  "coxlasso" = "CoxLasso"
)

# Apply the mapping
res1[, learner_id := rename_map[learner_id]]

# colors for learners
set1_colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                "#A65628", "#F781BF", "#17BECF","#6A5ACD")

# estimated Pareto front (ePF)
pfe1 = efs1$pareto_front(type = "estimated")

p1 = autoplot(efs1, pareto_front = "stepwise", theme = theme_minimal(base_size = 16)) +
  scale_color_manual(values = set1_colors) +
  # add ePF
  geom_line(data = pfe1, mapping = aes(
            x = .data[["n_features"]],
            y = .data[["surv.cindex"]]),
            color = "grey50", linetype = "dashed", linewidth = 0.7) +
  # add knee point based on ePF
  geom_point(data = kp1, aes(x = n_features, y = surv.cindex),
             shape = 4, size = 5, color = "red", stroke = 1.5) +
  labs(y = "Harrell's C-index") +
  theme(
    legend.position = "none",
    text = element_text(family = "Arial")
  )

# 'bad' example
file_name = file.path("bench", "efs", "wissel2023", "cnv", "efs_91.rds")
efs2 = readRDS(file_name)
kp2 = efs2$knee_points(type = "estimated") # 45

res2 = efs2$.__enclos_env__$private$.result
res2[, learner_id := rename_map[learner_id]]
# estimated Pareto front (ePF)
pfe2 = efs2$pareto_front(type = "estimated")

p2 = autoplot(efs2, pareto_front = "stepwise", theme =  theme_minimal(base_size = 16)) +
  scale_color_manual(values = set1_colors) +
  # add ePF
  geom_line(data = pfe2, mapping = aes(
            x = .data[["n_features"]],
            y = .data[["surv.cindex"]]),
            color = "grey50", linetype = "dashed", linewidth = 0.7) +
  # add knee point based on ePF
  geom_point(data = kp2, aes(x = n_features, y = surv.cindex),
             shape = 4, size = 3, color = "red", stroke = 1.5) +
  labs(y = "Harrell's C-index", color = "Model") +
  theme(
    legend.position = c(0.55, 0.55),
    text = element_text(family = "Arial")
  )
supfig4 = plot_grid(p1, p2, labels = c("a", "b"))
ggsave("bench/img/fig_s4.png", plot = supfig4, width = 10, height = 5, dpi = 600, bg = "white")

# Figure 4 ----
p_nog3 = p_nog +
  labs(title = NULL, x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    #legend.text = element_text(size = 14),
    #legend.title = element_text(size = 16),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "inside",
    legend.position.inside = c(0.135, 0.8),
    text = element_text(family = "Arial"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(face = "bold", size = 14),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )

p_srp5_xicor3 = p_srp5_xicor +
  labs(title = NULL, x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none",
    text = element_text(family = "Arial"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(face = "bold", size = 14),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )
fig4 = plot_grid(p_nog3, p_srp5_xicor3, nrow = 2, labels = c("a", "b"))
ggsave("bench/img/fig4.png", plot = fig4, width = 12, height = 9, dpi = 600, bg = "white")
ggsave("bench/img/fig4.pdf", plot = fig4, width = 12, height = 9, dpi = 600, bg = "white", device = cairo_pdf)

# Figure 5 ----
top_row = plot_grid(p_bf2, labels = "a") # with measure as harrell_c

bottom_row = plot_grid(
  p_times + theme(legend.position = "none"),
  p_times2 + theme(
    legend.position = "right",
    legend.margin = margin(l = -10, unit = "pt"),
    plot.margin = margin(t = 9, b = 9, l = 10, r = 20, unit = "pt")
  ),
  rel_widths = c(1.5, 2),
  labels = c("b", "c"),
  ncol = 2)

fig5 = plot_grid(
  top_row,
  bottom_row,
  nrow = 2,
  rel_heights = c(1, 1),
  label_fontfamily = "Arial"
)
ggsave("bench/img/fig5.png", plot = fig5, width = 11, height = 8, dpi = 600, bg = "white")
ggsave("bench/img/fig5.pdf", plot = fig5, width = 11, height = 8, dpi = 600, bg = "white", device = cairo_pdf)

# Sup. Fig. 8 ----
## p_rsf, p_cox => both with meas = "harrell_c"
supfig8 = plot_grid(
  p_rsf,
  p_cox + theme(legend.position = "none"),
  labels = c("a", "b"),
  nrow = 2,
  rel_heights = c(1, 1),
  label_fontfamily = "Arial"
)
ggsave("bench/img/fig_s8.png", plot = supfig8, width = 11, height = 8, dpi = 600, bg = "white")
