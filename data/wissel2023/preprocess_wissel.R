#' Extract PDAC multi-omics data from Wissel et al. (2023)
#' "Systematic comparison of multi-omics survival models reveals a widespread lack of noise resistance"
#' Execute: `Rscript data/wissel2023/preprocess_wissel.R`
suppressPackageStartupMessages({
  library(mlr3proba)
  library(mlr3misc)
  library(caret)
  library(dplyr)
  library(readr)
  library(tibble)
  library(usefun)
  library(ComplexHeatmap)
})

#' Download from Zenodo (https://zenodo.org/records/7529459) the `preprocessed.zip` file
#' Get the `PAAD_data_preprocessed.csv` file
#' Original data is from: https://gdc.cancer.gov/about-data/publications/pancanatlas
#' Download script: https://github.com/BoevaLab/Multi-omics-noise-resistance/blob/main/scripts/sh/download_data.sh
#' Data preprocessing script: https://github.com/BoevaLab/Multi-omics-noise-resistance/blob/main/noise_resistance/R/prep/prepare_tcga_data.R

# DATA ----
dataset_path = file.path("data", "wissel2023")
csv_datafile = file.path(dataset_path, "PAAD_data_preprocessed.csv")
all_data = readr::read_csv(file = csv_datafile, show_col_types = FALSE)
dim(all_data) # 108 x 84720

# CLINICAL ----
clinical = all_data |>
  select(starts_with(c("patient", "clinical_", "OS")))

clinical = clinical |>
  # rename variables
  rename(age = clinical_age_at_initial_pathologic_diagnosis) |>
  rename(tumor_stage = clinical_ajcc_pathologic_tumor_stage) |>
  rename(histological_type = clinical_histological_type) |>
  rename(sex = clinical_gender) |>
  rename(race = clinical_race) |>
  rename(time = OS_days, status = OS) |>
  relocate(patient_id, time, status) |> # put the outcome columns in front
  # all NAs for clinical stage => remove
  select(-clinical_clinical_stage) |>
  # convert days => months
  mutate(time = ceiling(time/30.44))

# (time, status) are okay
stopifnot(all(clinical$time > 0)) # no NAs or negative numbers here
stopifnot(clinical$status == 1 | clinical$status == 0) # dead (1) or censored (0)

# 3 patients with Stage III, 1 with Stage IV, 2 with NA
table(clinical$tumor_stage, useNA = "ifany")
# most (but not all) are PDACs
table(clinical$histological_type, useNA = "ifany")
# 2 Asian, 3 Black, 3 `NA`, rest white
table(clinical$race, useNA = "ifany")

# keep only white patients with PDAC, stage I and II
# `NA`s are treated as FALSE by `which()`, so below we also exclude ids with
# `NA` for tumor_stages or histological type
ids_to_keep = which(
  clinical$race == "WHITE" &
  clinical$histological_type == "Pancreas-Adenocarcinoma Ductal Type" &
  clinical$tumor_stage != "Stage III" &
  clinical$tumor_stage != "Stage IV"
)
length(ids_to_keep) # 81 patients

clinical = clinical |>
  slice(ids_to_keep) |>
  select(-c(histological_type, race)) |> # remove constant features
  # `clinical_tumor_stage` => collapse categories, convert to numeric type
  mutate(tumor_stage = case_when(
    tumor_stage == "Stage IA" ~ 0,
    tumor_stage == "Stage IB" ~ 0,
    tumor_stage == "Stage IIA" ~ 1,
    tumor_stage == "Stage IIB" ~ 1,
  )) |>
  # convert types
  mutate(sex = case_when(
    sex == "MALE" ~ 0,
    sex == "FEMALE" ~ 1
  )) |>
  mutate(age = as.integer(age)) |> # years
  mutate(time = as.integer(time)) |> # months
  mutate(status = as.numeric(status))

# check missing clinical features => no missing features
mlr3misc::map(clinical, \(.x) sum(is.na(.x)))

# OMICS (6) ----
# see "SurvBoard: Standardised Benchmarking for Multi-omics Cancer Survival Models"
# paper for preprocessing details
## CNV ----
cnv = all_data |>
  select(starts_with("cnv_")) |>
  rename_with(make.names) |>
  slice(ids_to_keep)
cnv |> unlist() |> unique() # 0,1,2,-1,-2

## GEX ----
gex = all_data |>
  select(starts_with("gex_")) |>
  rename_with(make.names) |>
  slice(ids_to_keep)
counts = gex |> unlist()
sum(counts < 0) # normalized counts, always > 0

## RPPA ----
rppa = all_data |>
  select(starts_with("rppa_")) |>
  rename_with(make.names) |>
  slice(ids_to_keep) # normalized expression (not standardized)

## Mutation ----
mutation = all_data |>
  select(starts_with("mutation_")) |>
  rename_with(make.names) |>
  slice(ids_to_keep)
values = mutation |> unlist()
table(values) # number of non-silent mutations, always >= 0 (up to 50)

## Methylation ----
meth = all_data |>
  select(starts_with("meth_")) |>
  rename_with(make.names) |>
  slice(ids_to_keep)
values = meth |> unlist()
stopifnot(values >= 0, values <= 1) # methyl beta values from 0 to 1

## miRNA ----
mirna = all_data |>
  select(starts_with("mirna_")) |>
  rename_with(make.names) |>
  slice(ids_to_keep)
# Data distributions in MIRNA are mixed??? Don't include!
# eg compare:
#' `mirna$mirna_hsa.miR.548x.3p |> density() |> plot()` # seems perfect N(0,1)
#' `mirna$mirna_hsa.miR.10b.5p  |> density() |> plot()` # seems (normalized) counts

data_list = list(clinical = clinical, gex = gex, cnv = cnv, rppa = rppa,
                 mutation = mutation, meth = meth)

# nrows the same
stopifnot(all(mlr3misc::map(data_list, nrow) == length(ids_to_keep)))

# nfeatures is correct
#' `patient_id`, `race`, `histological_type` and miRNAs are not included
stopifnot(sum(mlr3misc::map_dbl(data_list, ncol)) == ncol(all_data) - 3 - ncol(mirna))

# no missing data in general
mlr3misc::map(data_list, function(.data) {
  .data |> as.matrix() |> as.vector() |> is.na() |> sum()
})

# Non-filtered OMICS ----
# Here we keep some omics with all features for further analyses
survival_outcome = clinical |> select(patient_id, time, status)

## GEX
gex_data = data_list$gex

# remove constant features
feats_to_exclude = names(which(apply(gex_data, 2, var, na.rm = TRUE) == 0))
#gex_data = gex_data |> select(-all_of(feats_to_exclude))

# keep the 10000 most variant genes
features_to_keep =
  apply(gex_data, 2, var) |>
  sort(decreasing = TRUE) |>
  names() |>
  head(n = 10000)

gex_data = gex_data |> select(all_of(features_to_keep))
data = bind_cols(survival_outcome, gex_data)
gex_task = as_task_surv(x = data, time = "time", event = "status", id = "gex")
gex_task$set_col_roles(cols = "patient_id", roles = "name")
saveRDS(gex_task, file = file.path(dataset_path, "gex_task10000.rds"))

## Mutation
mut_data = data_list$mutation
feats_to_exclude = names(which(apply(mut_data, 2, var, na.rm = TRUE) == 0)) # remove constant features
mut_data = mut_data |> select(-all_of(feats_to_exclude)) # < 10000 features remain
data = bind_cols(survival_outcome, mut_data)
mut_task = as_task_surv(x = data, time = "time", event = "status", id = "mutation")
mut_task$set_col_roles(cols = "patient_id", roles = "name")
saveRDS(mut_task, file = file.path(dataset_path, "mut_task.rds"))

# PRE-FILTER (variance) ----
# keep only the top 2000 features with the highest variance per omic
n_features = 2000
data_list_flt = mlr3misc::map(data_list, function(.data) {
  if (ncol(.data) > n_features) {
    features_to_keep =
      apply(.data, 2, var) |>
      sort(decreasing = TRUE) |>
      names() |>
      head(n = n_features)

    .data = .data |> select(all_of(features_to_keep))
  }

  .data
})

# TASKS ----
task_list = mlr3misc::map(names(data_list_flt), function(.data_name) {
  if (.data_name == "clinical") {
    data = data_list_flt[["clinical"]]
  } else {
    data = bind_cols(survival_outcome, data_list_flt[[.data_name]])
  }

  # create survival task
  task = as_task_surv(x = data, time = "time", event = "status", id = .data_name)
  # add "patient_id" as row names (not a feature)
  task$set_col_roles(cols = "patient_id", roles = "name")

  task
})
names(task_list) = names(data_list_flt)

# ADD EXTRA CLINICAL VAR ----
#' Why? => running a Cox model on the 3 clinical features alone results in C-index
#' ~0.47, which is worse than random (0.5) performance!
#' Here we add the number of lymph nodes from the GDC pathology sheet, which increases
#' that a bit. See `https://portal.gdc.cancer.gov/projects/TCGA-PAAD`, we downloaded
#' the `data/wissel2023/clinical.project-tcga-paad.2025-03-05.tar.gz` file and
#' extracted the pathological data sheet
path_tbl = read_tsv(file = file.path(dataset_path, "pathology_detail.tsv"))
data = data.table(
  id = path_tbl$case_submitter_id,
  lymph_nodes_positive = as.integer(path_tbl$lymph_nodes_positive)
)
data

clinical_task = task_list$clinical$clone()
ids = clinical_task$row_names$row_name # patient ids
all(ids %in% data$id) # check

data = data[id %in% ids] # Keep only matching IDs
data = data[order(match(id, ids))] # Reorder IDs
stopifnot(data$id == ids) # check order
data$id = NULL # drop id column

# add lymph node number to clinical data
clinical_task$cbind(data)
clinical_task$missings() # still no missing values
task_list$clinical = clinical_task # put it back
data_list_flt[["clinical"]] = cbind(data_list_flt[["clinical"]], data) # same

# RESAMPLING FOR BENCHMARK ----
# 100 times => train/test split, stratified by censoring status
clinical_task = task_list$clinical$clone()
clinical_task$set_col_roles(cols = "status", add_to = "stratum")
ss = rsmp("subsampling", ratio = 0.8, repeats = 100)
set.seed(42)
ss$instantiate(clinical_task)

# METADATA ----
metadata = tibble(
  n_patients = task_list$clinical$nrow,
  n_modalities = length(task_list),
  n_gex_features = task_list$gex$n_features,
  n_cnv_features = task_list$cnv$n_features,
  n_rppa_features = task_list$rppa$n_features,
  n_mutation_features = task_list$mutation$n_features,
  n_meth_features = task_list$meth$n_features,
  n_clinical_features = task_list$clinical$n_features,
  n_total_features = sum(mlr3misc::map_int(task_list, \(.t) {.t$n_features})),
  n_events = sum(task_list$clinical$status() == 1),
  cens_rate = round(task_list$clinical$cens_prop(), digits = 2)
)
print(metadata)

# SAVE ALL DATA TO FILES ----
saveRDS(task_list, file = file.path(dataset_path, "task_list.rds"))
saveRDS(ss, file = file.path(dataset_path, "subsampling.rds"))

omic_ids = data.frame(omic_id = setdiff(names(task_list), "clinical"))
write_csv(omic_ids, file = file.path(dataset_path, "omic_ids.csv"), col_names = FALSE)
write_csv(metadata, file = file.path(dataset_path, "metadata.csv"))
write_csv(bind_cols(data_list_flt), file = file.path(dataset_path, "all_data_preprocessed.csv"))
