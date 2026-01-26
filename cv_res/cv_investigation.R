#' Investigation for impact of #CV-folds in the stability of feature ranking [Review Question]
#' We executed `bench/efs.R` for some configurations, results stored in `cv_res/cv_{nfolds}`
#' and analyzed further in this script
library(mlr3)
library(mlr3fselect)
library(mlr3viz)
library(tidyverse)
source("bench/helpers.R")

datasets = c("wissel2023", "cao2021")
models = c("xgb", "glmb", "coxboost", "coxlasso", "all") # all => hEFS
nfolds = c(3, 5, 7, 10)

grid_df = expand.grid(dataset = datasets, model = models, nfolds = nfolds)

row_seq = seq_len(nrow(grid_df))

# read data
res = lapply(row_seq, function(i) {
  row = grid_df[i, ]
  dataset = row$dataset
  model   = row$model
  nfolds  = row$nfolds
  res_path = file.path("cv_res", paste0("cv_", nfolds), dataset)

  efs_file = file.path(res_path, paste0("efs_", model, ".rds"))

  if (file.exists(efs_file)) {
    # get efs result object
    efs = readRDS(efs_file)

    # get execution time
    if (model != "all") {
      times_file = file.path(res_path, paste0("times_", model, ".csv"))
      time = readr::read_csv(times_file, show_col_types = FALSE)$time
    } else {
      time = NA # just sum, can compute later
    }

    # return result
    tibble(
      dataset = row$dataset,
      model = row$model,
      nfolds = row$nfolds,
      time = time,
      efs = list(efs)
    )
  }
}) |> bind_rows()

# compute the total time for the hEFS
res = res |>
  group_by(dataset, nfolds) |>
  mutate(
    time = ifelse(
      model == "all",
      sum(time[model != "all"], na.rm = TRUE),
      time
    )) |>
  ungroup()

model_lvls = c("CoxLasso", "CoxBoost", "GLMBoost-Cox", "XGBoost-Cox", "hEFS")
res = res |>
  mutate(dataset = case_when(
    dataset == "wissel2023" ~ "TCGA",
    dataset == "cao2021" ~ "CPTAC",
    TRUE ~ dataset
  )) |>
  mutate(model = case_when(
    model == "all" ~ "hEFS",
    model == "xgb" ~ "XGBoost-Cox",
    model == "glmb" ~ "GLMBoost-Cox",
    model == "coxlasso" ~ "CoxLasso",
    model == "coxboost" ~ "CoxBoost",
    TRUE ~ model
  )) |>
  mutate(
    model = factor(model, levels = model_lvls),
    n_features = map_int(efs, get_nfeats), # see `helpers.R`
    feats = map2(efs, n_features, ~ .x$feature_ranking(
      method = "sav",
      use_weights = TRUE,
      committee_size = .y
    )[["feature"]])
  )

# nfeats vs #CV-folds
res |>
  ggplot(aes(x = factor(nfolds), y = n_features, fill = factor(nfolds))) +
  geom_col() +
  facet_grid(rows = vars(dataset), cols = vars(model)) +
  labs(
    x = "Number of CV folds",
    y = "Number of selected features (Pareto front)",
    fill = "CV folds"
  ) +
  theme_minimal(base_size = 16, base_family = "Arial") +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  )

# execution time table
res |>
  select(dataset, model, nfolds, time) |>
  arrange(dataset, model, nfolds) |>
  print(n = 60)

# feat similarity across CV-folds
heatmap_list = list()

for(ds in unique(res$dataset)) {
  for(mod in unique(res$model)) {
    tmp = res |> filter(dataset == ds, model == mod)

    # Get selected features per nfold
    feats_list = tmp |> pull(feats)
    nfolds = tmp |> pull(nfolds)
    names(feats_list) = nfolds

    # Compute pairwise metrics
    mat = matrix(NA, nrow = length(nfolds), ncol = length(nfolds),
                 dimnames = list(nfolds, nfolds))

    for(i in seq_along(feats_list)) {
      for(j in seq_along(feats_list)) {
        x = feats_list[[i]]
        y = feats_list[[j]]
        #mat[i,j] = length(intersect(x,y)) / length(union(x,y)) # jaccard index
        mat[i,j] = length(intersect(x,y)) # number of common elements
      }
    }

    # Convert to long format for ggplot
    hm_df = as.data.frame(mat) |>
      mutate(nfold_i = rownames(mat)) |>
      pivot_longer(cols = -nfold_i, names_to = "nfold_j", values_to = "metric")

    hm_df$dataset = ds
    hm_df$model = mod
    heatmap_list[[paste(ds, mod, sep = "_")]] = hm_df
  }
}

df = bind_rows(heatmap_list)

# Keep only upper triangle (including diagonal)
fold_levels = c("3", "5", "7", "10")
df_upper = df |>
  filter(as.numeric(nfold_j) >= as.numeric(nfold_i)) |>
  mutate(
    nfold_i = factor(nfold_i, levels = fold_levels),
    nfold_j = factor(nfold_j, levels = fold_levels),
    model = factor(model, levels = model_lvls)
  )

df_upper |>
  ggplot(aes(x = nfold_j, y = nfold_i, fill = metric)) +
    geom_tile(color = "white") + # heatmap tiles
    #scale_fill_viridis_c() +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      # midpoint = 0.5,
      # limits = c(0, 1),
      # name = "Jaccard"
      midpoint = 5,
      limits = c(0, 10),
      name = "No. Common Features"
    ) +
    geom_text(aes(label = round(metric, 2)), size = 5) +  # optional: show values
    facet_grid(rows = vars(dataset), cols = vars(model)) +
    labs(x = "#CV-Folds", y = "#CV-Folds", fill = "metric") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold")
    )
#
