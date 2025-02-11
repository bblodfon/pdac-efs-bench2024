#' Investigate how 1) number of resamples 2) number of learners influence
#' the number of selected features and exact features selected with the different voting theory methods (AV vs SAV)
library(dplyr)
library(tibble)
library(mlr3misc)
library(future)
library(future.apply)
library(progressr)
library(ggplot2)

# Data ----
## GEX or MUTATION data
efs_all = readRDS(file = "sensitivity_analysis/wissel2023/gex/efs_all.rds")
efs_all = readRDS(file = "sensitivity_analysis/wissel2023/mutation/efs_all.rds")
efs_all

autoplot(efs_all) +
  scale_color_brewer(palette = "Set1") +
  xlim(c(0, 100))
autoplot(efs_all, pareto_front = "estimated") +
  scale_color_brewer(palette = "Set1") +
  xlim(c(0, 100))

# Reference feature sets ----
# find REFERENCE (best) feature subset using all models and resamplings
max_nfeatures = min(max(efs_all$result$n_features), 100)
max_nfeatures
nfeats = efs_all$
  knee_points(type = "estimated", max_nfeatures = max_nfeatures)$
  n_features
print(nfeats)

av_all = efs_all$feature_ranking(method = "av", use_weights = TRUE,
                                 committee_size = nfeats)$feature
sav_all = efs_all$feature_ranking(method = "sav", use_weights = TRUE,
                                  committee_size = nfeats)$feature

# DICE/Jaccard ----
# DICE => % of the COMMON elements in the average size of the combined set
# less sensitive to differences in set sizes
dice = function(x1, x2) {
  2 * length(intersect(x1, x2)) / (length(x1) + length(x2))
}
dice(av_all, sav_all)

# JACCARD => % of the COMMON elements in the union of the two sets
# penalizes differences in set sizes
jaccard = function(x1, x2) {
  length(intersect(x1, x2)) / length(union(x1, x2))
}
jaccard(av_all, sav_all)

# example
x1 = 1:10
x2 = 2:15
x2 = 7:16
jaccard(x1, x2)
dice(x1, x2)

#' `lrn_ids` => learner ids to filter `efs` result
#' `group_id` => `id` for the learner group to put in the result
get_selfeats = function(efs, lrn_ids, rsmp_ids, methods = c("av", "sav")) {
  efs_copy = efs$clone()

  # subset learners and resampling iterations
  efs_copy$.__enclos_env__$private$.result =
    efs_copy$.__enclos_env__$private$.result[
      learner_id %in% lrn_ids & resampling_iteration %in% rsmp_ids
    ]

  # cap max features for estimated Pareto front to 100
  max_nfeatures = min(max(efs_copy$result$n_features), 100)

  # estimated Pareto front: select number of features
  nfeats = suppressWarnings(
    efs_copy$knee_points(type = "estimated", max_nfeatures = max_nfeatures)$n_features
  )
  if (is.na(nfeats)) stop("Too few points to calculate the pareto front")

  # subset resamplings
  data = list()
  i = 1
  for (method in methods) {
    # no weights
    # sel_feats = efs_copy$feature_ranking(method = method, use_weights = FALSE,
    #                                      committee_size = nfeats)$feature
    # data[[i]] = tibble::tibble(
    #   method = method,
    #   use_weights = FALSE,
    #   pf_nfeats = nfeats,
    #   sel_feats = list(sel_feats)
    # )
    # i = i + 1

    # with weights
    sel_feats = efs_copy$feature_ranking(method = method, use_weights = TRUE,
                                         committee_size = nfeats)$feature
    data[[i]] = tibble::tibble(
      method = method,
      #use_weights = TRUE,
      pf_nfeats = nfeats,
      sel_feats = list(sel_feats)
    )
    i = i + 1
  }

  dplyr::bind_rows(data)
}

# example
get_selfeats(efs_all, lrn_ids = "aorsf.fselector", rsmp_ids = c(1,10,100))

# RUN ANALYSIS ----
# get all learner ids
learner_ids = unique(efs_all$result$learner_id)

## Define learner groups ----
# Define the following 5 groups of learners:
# RSF group (3)
rsf_ids = learner_ids[grepl(pattern = "rsf", learner_ids)]
# XGBoost group (2)
xgb_ids = learner_ids[grepl(pattern = "xgb", learner_ids)]
# GLMBoost + Coxboost group (3) => component-wise boosting
cwb_ids = learner_ids[grepl(pattern = "glmb|coxboost", learner_ids)]
# CoxLasso (1) => id = "coxlasso"

# All possible combos: 2^4 - 1 = 15

# make (group_id => lrn_ids) mapping
id_map = list(rsf = rsf_ids, xgb = xgb_ids, cwb = cwb_ids, coxlasso = "coxlasso")

# make (all-group-id combinations - powerset => lrn_ids) mapping
group_names = names(id_map)
group_combos = unlist(
  lapply(1:length(group_names), function(k) combn(group_names, k, simplify = FALSE)),
  recursive = FALSE
)
length(group_combos) # 15

group_id_combos = unlist(lapply(group_combos, paste0, collapse = " + "))
group_id_combos

powerset_id_map = lapply(group_combos, function(groups) {
  unlist(unname(id_map[groups])) # Flatten IDs from selected groups
})
names(powerset_id_map) = group_id_combos
powerset_id_map

## Define each learner as a separate group ----
if (FALSE) {
  # It takes TOO much time...
  id_map = list(
    rsf_logrank = "rsf_logrank.fselector",
    rsf_maxstat = "rsf_logrank.maxstat",
    aorsf = "aorsf.fselector",
    xgb_cox = "xgb_cox.fselector",
    xgb_aft = "xgb_aft_log.fselector",
    glmb_cox = "glmb_cox",
    glmb_aft = "glmb_loglog",
    coxboost = "coxboost",
    coxlasso = "coxlasso"
  )
  group_names = names(id_map)
  group_combos = unlist(
    lapply(1:length(group_names), function(k) combn(group_names, k, simplify = FALSE)),
    recursive = FALSE
  )
  length(group_combos) # 511
  group_id_combos = unlist(lapply(group_combos, paste0, collapse = " + "))
  group_id_combos

  powerset_id_map = lapply(group_combos, function(groups) {
    unlist(unname(id_map[groups])) # Flatten IDs from selected groups
  })
  names(powerset_id_map) = group_id_combos
  tail(powerset_id_map)
}

n_resamples = seq(10, 100, 10)
n_resamples
total_rsmp_ids = unique(efs_all$result$resampling_iteration)

process_group = function(group_name, n_rsmps, index, total_rsmp_ids,
                         powerset_id_map, av_all, sav_all, p) {
  #cat("Lrn_grp: ", group_name, ", #Rsmps: ", n_rsmps,
  #    ", Iter: ", index, "\n", sep = "")

  #' Notify progress via `p = progressr::progressor()`
  p(sprintf("Lrn_grp: %s, #Rsmps: %i, #Iter: %i", group_name, n_rsmps, index))

  # Select a random sample of resamplings
  rsmp_ids = sample(x = total_rsmp_ids, size = n_rsmps)
  # Get the learner subset
  lrn_ids = powerset_id_map[[group_name]]
  n_learners = length(lrn_ids)

  # Get selected features via the estimated Pareto Front method
  res = tryCatch(
    {
      get_selfeats(efs_all, lrn_ids, rsmp_ids)
    }, error = function(e) {
      message("Error: ", e$message)
      return(NULL) # Return NULL to skip this iteration in case of errors
    }
  )

  if (is.null(res)) return(NULL)

  # Calculate similarity to reference feature subsets
  res = res |>
    add_column(
      group_name = group_name,
      n_learners = n_learners,
      n_rsmps = n_rsmps,
      .before = 1
    ) |>
    mutate(
      ref_feats = case_when(method == "av" ~ list(av_all), .default = list(sav_all))
    ) |>
    rowwise() |>
    mutate(
      jaccard = jaccard(sel_feats, ref_feats),
      dice = dice(sel_feats, ref_feats)
    ) |>
    ungroup() |>
    select(-ref_feats)

  return(res)
}

# Prepare grid for parallel processing
data_grid = expand.grid(
  group_name = group_id_combos,
  n_rsmps = n_resamples,
  index = 1:25, # For each `index` iteration, we select randomly `n_rsmps` resamples
  stringsAsFactors = FALSE
)
head(data_grid)
nrow(data_grid)

#' Run in parallel using `future.apply::future_lapply` + Progress bars
# Enable progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers("progress")

# Parallellize
plan("multisession", workers = 14)

# write it like this so that progress is reported correctly
execute = function() {
  row_seq = 1:nrow(data_grid)

  # Progress tracking
  p = progressr::progressor(along = row_seq)

  data = future.apply::future_lapply(row_seq, function(row_id) {
    combo = data_grid[row_id, ]
    process_group(
      group_name = combo$group_name,
      n_rsmps = combo$n_rsmps,
      index = combo$index,
      total_rsmp_ids = total_rsmp_ids,
      powerset_id_map = powerset_id_map,
      av_all = av_all,
      sav_all = sav_all,
      p
    )
  }, future.seed = TRUE)

  bind_rows(data)
}

# get Jaccard and DICE results
results = execute()

# save results
saveRDS(results, file = "sensitivity_analysis/results_wissel2023_gex.rds")
saveRDS(results, file = "sensitivity_analysis/results_wissel2023_mutation.rds")

#results = readRDS(file = "sensitivity_analysis/results_wissel2023_gex.rds")

# Summary stats ----
summary(results$pf_nfeats)

results |>
  filter(jaccard > 0.9) |> # or dice
  group_by(method) |>
  summarise(
    median_n_learners = median(n_learners),
    min_n_learners = min(n_learners),
    median_n_rsmps = median(n_rsmps),
    min_n_rsmps = min(n_rsmps)
  )

# Plot ----
# Choose particular `n_rsmps` for visualization purposes
sel_n_rsmps = c(10, 30, 50, 70, 90, 100)

## Number of selected features using the estimated Pareto Front
results |>
  filter(method == "sav") |>
  #filter(method == "sav", n_rsmps %in% sel_n_rsmps) |>
  mutate(lrn_grp = case_when(
    n_learners == 1 | n_learners == 2 ~ "1-2",
    n_learners == 3 | n_learners == 4 ~ "3-4",
    n_learners == 5 | n_learners == 6 ~ "5-6",
    n_learners == 7 | n_learners == 8 ~ "7-8",
    n_learners == 9 ~ "9"
  )) |>
  ggplot(aes(x = as.factor(lrn_grp), y = pf_nfeats)) +
  #ggplot(aes(x = as.factor(n_learners), y = pf_nfeats)) +
  geom_boxplot() +
  #facet_grid(rows = vars(n_rsmps)) +
  labs(
    title = "Dataset: TCGA (GEX) - Selecting Features via the Estimated Pareto Front",
    x = "Number of Learners",
    y = "Number of Selected Features"
  )

## Jaccard, x => group by #learners
results |>
  filter(method == "av", n_rsmps %in% sel_n_rsmps) |>
  mutate(lrn_grp = case_when(
    n_learners == 1 | n_learners == 2 ~ "1-2",
    n_learners == 3 | n_learners == 4 ~ "3-4",
    n_learners == 5 | n_learners == 6 ~ "5-6",
    n_learners == 7 | n_learners == 8 ~ "7-8",
    n_learners == 9 ~ "9"
  )) |>
  ggplot(aes(x = as.factor(lrn_grp), y = jaccard)) +
  #ggplot(aes(x = as.factor(n_learners), y = jaccard)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.9, colour = "red", linetype = "dashed") +
  facet_grid(rows = vars(n_rsmps)) +
  labs(
    #title = "Dataset: TCGA (GEX) - Sensitivity Analysis (SAV)",
    title = "Dataset: TCGA (Mutation) - Sensitivity Analysis (SAV)",
    x = "Number of Learners",
    y = "Jaccard Index"
  )

## Jaccard, x => group by the combined learner group name
results |>
  filter(method == "sav", n_rsmps %in% sel_n_rsmps) |>
  ggplot(aes(x = forcats::fct_reorder(group_name, n_learners), y = jaccard)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.9, colour = "red", linetype = "dashed") +
  facet_grid(rows = vars(n_rsmps)) +
  labs(
    #title = "Dataset: TCGA (GEX) - Sensitivity Analysis (SAV)",
    title = "Dataset: TCGA (Mutation) - Sensitivity Analysis (SAV)",
    x = "Combined Learner Group",
    y = "Jaccard Index"
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

## DICE, x => group by #learners
results |>
  filter(method == "sav", n_rsmps %in% sel_n_rsmps) |>
  mutate(lrn_grp = case_when(
    n_learners == 1 | n_learners == 2 ~ "1-2",
    n_learners == 3 | n_learners == 4 ~ "3-4",
    n_learners == 5 | n_learners == 6 ~ "5-6",
    n_learners == 7 | n_learners == 8 ~ "7-8",
    n_learners == 9 ~ "9"
  )) |>
  ggplot(aes(x = as.factor(lrn_grp), y = dice)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.9, colour = "red", linetype = "dashed") +
  facet_grid(rows = vars(n_rsmps)) +
  labs(
    title = "Dataset: TCGA (GEX) - Sensitivity Analysis (SAV)",
    x = "Number of Learners",
    y = "DICE score"
  )

## DICE, x => group by the combined learner group name
results |>
  filter(method == "sav", n_rsmps %in% sel_n_rsmps) |>
  ggplot(aes(x = forcats::fct_reorder(group_name, n_learners), y = dice)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.9, colour = "red", linetype = "dashed") +
  facet_grid(rows = vars(n_rsmps)) +
  labs(
    title = "Dataset: TCGA (GEX) - Sensitivity Analysis (SAV)",
    x = "Combined Learner Group",
    y = "DICE score"
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))


