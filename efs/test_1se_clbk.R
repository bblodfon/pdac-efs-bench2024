library(ggplot2)

# run with no 1se clbk for glmboost
efs_all = readRDS(file = "efs/gex/efs_all_100.rds")
# remove all other learners
efs_all$.__enclos_env__$private$.result = efs_all$.__enclos_env__$private$.result[learner_id == "glmb_cox" | learner_id == "glmb_loglog"]
efs_all

# run with glmboost learners only + 1se clbk
efs_glmb = readRDS(file = "efs/gex/efs_glmb_clbk.rds")
efs_glmb

id = "glmb_cox"
df = data.table(
  no_clbk   = efs_all$result[learner_id == id]$n_features,
  with_clbk = efs_glmb$result[learner_id == id]$n_features
)

df |> tidyr::pivot_longer(cols = tidyselect::everything(), names_to = "Group", values_to = "n_features") |>
ggplot(aes(x = Group, y = n_features, fill = Group)) +
  geom_boxplot(outlier.color = "red", outlier.size = 2) +
  labs(
    title = id,
    x = "Group",
    y = "Number of Features"
  ) +
  scale_fill_manual(values = c("#0073C2FF", "#EFC000FF")) + # Custom colors
  theme_minimal() +
  theme(
    legend.position = "none", # Remove legend as it's redundant
    plot.title = element_text(hjust = 0.5) # Center the title
  )
