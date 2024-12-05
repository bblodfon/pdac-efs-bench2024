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
csv_datafile = "data/wissel2023/PAAD_data_preprocessed.csv"
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

# keep only the top 2000 features with the highest variance per omic
n_features = 2000
data_list = mlr3misc::map(data_list, function(.data) {
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
survival_outcome = data_list$clinical |> select(patient_id, time, status)

task_list = mlr3misc::map(names(data_list), function(.data_name) {
  if (.data_name == "clinical") {
    data = data_list[["clinical"]]
  } else {
    data = bind_cols(survival_outcome, data_list[[.data_name]])
  }

  # create survival task
  task = as_task_surv(x = data, time = "time", event = "status", id = .data_name)
  # add "patient_id" as row names (not a feature)
  task$set_col_roles(cols = "patient_id", roles = "name")

  task
})
names(task_list) = names(data_list)
# task_list

# STRATIFIED TRAIN/TEST SPLIT ----
task_clinical = task_list$clinical$clone()

table(task_clinical$data(cols = "sex")) # 46 males (0), 35 females (1)
table(task_clinical$data(cols = "tumor_stage")) # 8 => stage I (0), 73 => stage II (1)
table(task_clinical$data(cols = "status")) # 33 => censored, 48 => events

# stratify on status and tumor_stage
task_clinical$set_col_roles(cols = c("status", "tumor_stage"), add_to = "stratum")
task_clinical$strata
# chose ratio to have 30 patients in the test set
set.seed(42)
part = partition(task_clinical, ratio = 0.63)
length(part$test)

# check stratification is done properly
task_clinical$cens_prop()
task_clinical$cens_prop(rows = part$train)
task_clinical$cens_prop(rows = part$test)

prop.table(table(task_clinical$data(cols = "tumor_stage")[[1L]]))
prop.table(table(task_clinical$data(rows = part$train, cols = "tumor_stage")[[1L]]))
prop.table(table(task_clinical$data(rows = part$test, cols = "tumor_stage")[[1L]]))

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
saveRDS(task_list, file = "data/wissel2023/task_list.rds")
saveRDS(part, file = "data/wissel2023/data_split.rds")
write_csv(metadata, file = "data/wissel2023/metadata.csv")

all_data_preprocessed = bind_cols(data_list)
write_csv(all_data_preprocessed, file = "data/wissel2023/all_data_preprocessed.csv")
write_csv(all_data_preprocessed[part$train, ],
          file = "data/wissel2023/all_data_preprocessed_train.csv")
write_csv(all_data_preprocessed[part$test, ],
          file = "data/wissel2023/all_data_preprocessed_test.csv")
