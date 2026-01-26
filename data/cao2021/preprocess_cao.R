#' Extract PDAC multi-omics data from Cao et al. (2021)
#' "Proteogenomic characterization of pancreatic ductal adenocarcinoma"
#' Execute: `Rscript data/cao2021/preprocess_cao.R`
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
source("data/helper.R")

#' Download from http://www.linkedomics.org/data_download/CPTAC-PDAC/ all the files,
#' put them in a directory called `PDAC_LinkedOmics_Data`.
data_dir = file.path("data", "cao2021")
dataset_path = file.path(data_dir, "PDAC_LinkedOmics_Data")
# extra - not used in this pre-processing
# mirna_file = file.path(dataset_path, "microRNA_TPM_log2_Tumor.cct")
# crna_file = file.path(dataset_path, "circRNA_RSEM_UQ_log2_Tumor.cct")
# meth_file = file.path(dataset_path, "methylation_betaValue_Tumor.cct")

# CLINICAL ----
clinical_file = file.path(dataset_path, "clinical_table_140.tsv")
clinical = readr::read_tsv(file = clinical_file, show_col_types = FALSE)

clinical = clinical |>
  # keep only PDACs
  filter(histology_diagnosis == "PDAC") |>
  # 5 `NA`s for followup time (can't do survival analysis with these)
  filter(!is.na(follow_up_days), follow_up_days > 0) |>
  # 4 `NA`s with tumor_stage
  filter(!is.na(tumor_stage_pathological)) |>
  # rename variables
  rename(
    patient_id = case_id,
    time = follow_up_days,
    status = vital_status,
    lost_to_follow_up = is_this_patient_lost_to_follow_up,
    tumor_stage = tumor_stage_pathological,
    lymph_nodes_positive = number_of_lymph_nodes_positive_for_tumor,
    lymph_nodes_examined = number_of_lymph_nodes_examined
  ) |>
  # convert days => months
  mutate(time = ceiling(time/30.44)) |>
  # put id and outcome columns first
  relocate(patient_id, time, status) |>
  # `status` = {Deceased, Living}
  # `lost_to_follow_up` = {Yes, No} => this one is confusing:
  # 1) If "No" => then `status` can be interpreted as usual (46:Deceased, 38:Living)
  # 2) If "Yes" and `status` == 'Living' => censored (16 patients)
  # 3) If "Yes" and `status` == 'Deceased' => ??? (21 patients are reported dead
  #    but marked lost) - we keep these as 'dead'/events
  #
  # Remove columns
  select(-c(
    # we discard `lost_to_follow_up` and keep the `status` as censoring indicator
    "lost_to_follow_up",
    # only PDACs left
    "histology_diagnosis",
    # for all it's the same ('yes')
    "tumor_included_for_the_study",
    # tumor_site: {body, head, tail, combos of these} => most of these are 'head's (99 patients)
    "tumor_site",
    # we don't care
    "normal_included_for_the_study",
    # remove ethnographic variables (most NA's for race either way)
    "race", "participant_country",
    # remove due to imbalance (2 classes): Multifocal (3) vs Unifocal  (135)
    "tumor_focality",
    # most values are "Not identified"
    "tumor_necrosis",
    # remove BMI due to it degrading baseline CoxPH C-index performance
    "bmi",
    # didn't boost baseline performance
    "perineural_invasion",
    # the pathological tumor stage is enough so we discard the following staging variables
    "pathologic_staging_primary_tumor_pt",
    "pathologic_staging_regional_lymph_nodes_pn",
    "pathologic_staging_distant_metastasis_pm",
    "clinical_staging_distant_metastasis_cm",
    # lots of categories, free-text most of the values
    "additional_pathologic_findings",
    "medical_condition",
    # 4 classes: ~half of the patients have no R0 (No residual tumor), rest have
    # R1/R2 and RX (cannot be assessed) in equal proportions => can't decide on
    # how to encode this optimally for modeling
    "residual_tumor",
    # these variables have values such as "X;Y;Z" so hard to encode these without
    # introducing more complexity
    "Acinar_fraction",
    "Fat_fraction",
    "Islet_fraction",
    "Stromal_fraction",
    "Inflammation_fraction",
    "Muscle_fraction",
    "Neoplastic_cellularity",
    "Non_neoplastic_duct",
    # exclude smoking and alcohol history (>=5 classes each)
    "alcohol_consumption",
    "tobacco_smoking_history",
    # 56 patients have `na` (these also have `status == Living`), 13 are labelled
    # as 'unknown' (and have `status == Deceased`), 49 have "pancreatic carcinoma".
    # 6 more causes have 1-3 patients. Exclude as we are interested only in
    # right-censoring (and maybe not enough data for competing risks either)
    "cause_of_death"
    )) |>
  # collapse tumor stage factor levels {I,II,II,IV} => {0,1,2,3}
  mutate(tumor_stage = case_when(
    tumor_stage == "Stage IA" ~ 0,
    tumor_stage == "Stage IB" ~ 0,
    tumor_stage == "Stage IIA" ~ 1,
    tumor_stage == "Stage IIB" ~ 1,
    tumor_stage == "Stage III" ~ 2,
    tumor_stage == "Stage IV" ~ 3
  )) |>
  # convert types
  mutate(status = case_when(
    status == "Living" ~ 0,
    status == "Deceased" ~ 1
  )) |>
  mutate(sex = case_when(
    sex == "Male" ~ 0,
    sex == "Female" ~ 1
  )) |>
  # the next + the `lympd_nodes_*` variables boost baseline performance
  mutate(lymph_vascular_invasion = case_when(
    lymph_vascular_invasion == "Not identified" ~ 0, # 34 (NO)
    lymph_vascular_invasion == "Indeterminate" ~ 0, # just 5 patients, so okay to collapse
    lymph_vascular_invasion == "Present" ~ 1 # 86 (YES)
  )) |>
  mutate(age = as.integer(age)) |> # years
  mutate(time = as.integer(time)) |> # months
  mutate(lymph_nodes_positive = as.integer(lymph_nodes_positive)) |>
  mutate(lymph_nodes_examined = as.integer(lymph_nodes_examined)) #|>
  # feat engineering => RSF doesn't benefit from it, CoxPH becomes bit worse
  #mutate(ln_ratio = lymph_nodes_positive / lymph_nodes_examined) |>
  #mutate(ln_diff = lymph_nodes_examined - lymph_nodes_positive) |>
  #select(-lymph_nodes_positive, -lymph_nodes_examined)

# I-IV tumor stages
table(clinical$tumor_stage)

# (time, status) are okay
stopifnot(all(clinical$time > 0)) # no NAs or negative numbers here
stopifnot(clinical$status == 1 | clinical$status == 0) # dead (1) or censored (0)

# check missing clinical features => no missing features
mlr3misc::map(clinical, \(.x) sum(is.na(.x)))

# keep patients ids for omic (row/observation) filtering
clin_pids = clinical$patient_id

# INVESTIGATE: Baseline performance (clinical data + CoxPH/RSF)
if (FALSE) {
  library(mlr3extralearners)
  library(mlr3viz)

  task = as_task_surv(x = clinical, time = "time", event = "status", id = "upto-IV")
  task$set_col_roles(cols = "patient_id", roles = "name")
  task$set_col_roles(cols = "status", add_to = "stratum")
  task$set_col_roles(cols = "tumor_stage", add_to = "stratum")

  clinical2 = clinical |> filter(tumor_stage < 3)
  task2 = as_task_surv(x = clinical2, time = "time", event = "status", id = "upto-III")
  task2$set_col_roles(cols = "patient_id", roles = "name")
  task2$set_col_roles(cols = "status", add_to = "stratum")
  task2$set_col_roles(cols = "tumor_stage", add_to = "stratum")

  clinical3 = clinical2 |> filter(tumor_stage < 2)
  task3 = as_task_surv(x = clinical3, time = "time", event = "status", id = "I-II")
  task3$set_col_roles(cols = "patient_id", roles = "name")
  task3$set_col_roles(cols = "status", add_to = "stratum")
  task3$set_col_roles(cols = "tumor_stage", add_to = "stratum")

  cox = lrn("surv.coxph")
  rsf = lrn("surv.ranger", num.trees = 250)

  grid = benchmark_grid(
    tasks = list(task, task2, task3),
    learners = list(cox, rsf),
    resamplings = list(rsmp("subsampling", repeats = 100, ratio = 0.8))
  )

  set.seed(42)
  bm = benchmark(design = grid)
  bm$aggregate() # bm2 => status stratified, bm => both, bm3 => none

  autoplot(bm) +
    geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", linewidth = 0.8)
}

# OMICS (5) ----
#' Note: In all omics datasets, we will use `Gene` level features and get only
#' the `Tumor` samples, unless otherwise specified

## GEX ----
#' Note: we have added a `names` string in the 1st line, 1st column of the file
#' for easier column name parsing
gex_file = file.path(dataset_path, "mRNA_RSEM_UQ_log2_Tumor.cct")
gex_data = readr::read_tsv(gex_file, show_col_types = FALSE)
gex_feats = make.names(paste0("gex_", gex_data$names))

gex_mat = gex_data |> select(!names) |> t()
gex_pids = rownames(gex_mat) # `pids` => patient ids
stopifnot(all(clin_pids %in% gex_pids)) # check
rownames(gex_mat) = NULL
colnames(gex_mat) = gex_feats
gex_tbl = gex_mat |> as_tibble()

dim(gex_tbl) # 140 28057
sum(gex_tbl < 0) # normalized counts, always > 0
any(is.na(gex_tbl)) # no NAs
sum(apply(gex_tbl, 2, var, na.rm = TRUE) == 0) # some constant features in here

## CNV ----
#' Note: we have added a `names` string in the 1st line, 1st column of the file
#' for easier column name parsing
cnv_file = file.path(dataset_path, "SCNA_log2_gene_level.cct")
cnv_data = readr::read_tsv(cnv_file, show_col_types = FALSE)
cnv_feats = make.names(paste0("cnv_", cnv_data$names))

cnv_mat = cnv_data |> select(!names) |> t()
cnv_pids = rownames(cnv_mat) # `pids` => patient ids
stopifnot(all(clinical$patient_id %in% cnv_pids)) # check
rownames(cnv_mat) = NULL
colnames(cnv_mat) = cnv_feats
cnv_tbl = cnv_mat |> as_tibble()

# Gene level copy number data, log2 ratio (compared to normal), centered, like in Osipov's dataset
summary(mlr3misc::map_dbl(cnv_tbl, \(.x) mean(.x, na.rm = TRUE))) # mean ~ 0
summary(mlr3misc::map_dbl(cnv_tbl, \(.x) sd(.x, na.rm = TRUE))) # mean ~ 0.08

any(is.na(cnv_tbl)) # some NAs - do imputation
cnv_tbl = impute(cnv_tbl)
# 18 features with more than 10% NAs across samples were removed.
# Imputing 45 features with <=10% missing values using median feature values.
any(is.na(cnv_tbl)) # no NAs anymore

sum(apply(cnv_tbl, 2, var, na.rm = TRUE) == 0) # no constant features
dim(cnv_tbl) # 140 19888

## Proteomics ----
#' Note: we have added a `names` string in the 1st line, 1st column of the file
#' for easier column name parsing
prot_file = file.path(dataset_path, "proteomics_gene_level_MD_abundance_tumor.cct")
prot_data = readr::read_tsv(prot_file, show_col_types = FALSE)
prot_feats = make.names(paste0("prot_", prot_data$names))

prot_mat = prot_data |> select(!names) |> t()
prot_pids = rownames(prot_mat) # `pids` => patient ids
stopifnot(all(clinical$patient_id %in% prot_pids)) # check
rownames(prot_mat) = NULL
colnames(prot_mat) = prot_feats
prot_tbl = prot_mat |> as_tibble()

any(is.na(prot_tbl)) # some NAs - do imputation
prot_tbl = impute(prot_tbl)
# 4798 features with more than 10% NAs across samples were removed.
# Imputing 1097 features with <=10% missing values using median feature values.
any(is.na(prot_tbl)) # no NAs anymore

dim(prot_tbl) # 140 6864
sum(prot_tbl < 0) # positive, median-normalized intensity values
sum(apply(prot_tbl, 2, var, na.rm = TRUE) == 0) # no constant features

# distributions in the proteomics data look very normal to me so no further
# pre-processing (scaling during modeling)

## Phosphoproteomics ----
#' Note: we have added a `names` string in the 1st line, 1st column of the file
#' for easier column name parsing
phos_file = file.path(dataset_path, "phosphoproteomics_gene_level_MD_abundance_tumor.cct")
phos_data = readr::read_tsv(phos_file, show_col_types = FALSE)
phos_feats = make.names(paste0("phos_", phos_data$names))

phos_mat = phos_data |> select(!names) |> t()
phos_pids = rownames(phos_mat) # `pids` => patient ids
stopifnot(all(clinical$patient_id %in% phos_pids)) # check
rownames(phos_mat) = NULL
colnames(phos_mat) = phos_feats
phos_tbl = phos_mat |> as_tibble()

any(is.na(phos_tbl)) # some NAs - do imputation
phos_tbl = impute(phos_tbl)
# 4169 features with more than 10% NAs across samples were removed.
# Imputing 816 features with <=10% missing values using median feature values.
any(is.na(phos_tbl)) # no NAs anymore

dim(phos_tbl) # 140 3835
sum(phos_tbl < 0) # positive, median-normalized intensity values
sum(apply(phos_tbl, 2, var, na.rm = TRUE) == 0) # no constant features

## N-Glycoproteomics ----
#' Note: we use the `Peptide` level here (and not `Site`) as our focus is to
#' reduce technical variability (less noise) and simplify preprocessing.
#' Peptide level is an aggregated version across one or more glyco-sites.
ngly_file = file.path(dataset_path, "N-glycoproteomics_peptide_level_ratio_tumor.cct")
ngly_data = readr::read_tsv(ngly_file, show_col_types = FALSE)
ngly_feats = make.names(paste0("ngly_", ngly_data$Sequence)) # peptide sequence

ngly_mat = ngly_data |> select(-Sequence, -Gene) |> t()
ngly_pids = rownames(ngly_mat) # `pids` => patient ids
stopifnot(all(clinical$patient_id %in% ngly_pids)) # check
rownames(ngly_mat) = NULL
colnames(ngly_mat) = ngly_feats
ngly_tbl = ngly_mat |> as_tibble()

any(is.na(ngly_tbl)) # some NAs - do imputation
ngly_tbl = impute(ngly_tbl)
# 27816 features with more than 10% NAs across samples were removed.
# Imputing 1425 features with <=10% missing values using median feature values.
any(is.na(ngly_tbl)) # no NAs anymore
sum(apply(ngly_tbl, 2, var, na.rm = TRUE) == 0) # no constant features

dim(ngly_tbl) # 140 2844
# Expression values (TMT, Log2ratio)
summary(mlr3misc::map_dbl(ngly_tbl, \(.x) mean(.x, na.rm = TRUE))) # mean ~ 0
summary(mlr3misc::map_dbl(ngly_tbl, \(.x) sd(.x, na.rm = TRUE))) # mean ~ 0.8
# distributions in the N-glycoproteomics data look very normal to me so no further
# pre-processing (scaling during modeling)

# COMBINE MODALITIES ----
# filter and reorder omics data by matching patient IDs
filter_omic = function(omic_tbl, omic_pids, clin_pids) {
  matched_idx = match(clin_pids, omic_pids)
  omic_tbl[matched_idx, , drop = FALSE]
}

data_list = list(
  clinical = clinical,
  gex = filter_omic(gex_tbl, gex_pids, clin_pids),
  cnv = filter_omic(cnv_tbl, gex_pids, clin_pids),
  prot = filter_omic(prot_tbl, gex_pids, clin_pids),
  phos = filter_omic(phos_tbl, gex_pids, clin_pids),
  ngly = filter_omic(ngly_tbl, gex_pids, clin_pids)
)

# nrows the same
stopifnot(all(mlr3misc::map(data_list, nrow) == length(clin_pids)))

# no missing data in general
mlr3misc::map(data_list, function(.data) {
  .data |> as.matrix() |> as.vector() |> is.na() |> sum()
})

# check #features
map(data_list, ncol) |> unlist()
# clinical      gex      cnv     prot     phos     ngly
# 10            28057    19888   6864     3835     2844

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
# check
map(data_list_flt, ncol) |> unlist()

# TASKS ----
task_list = mlr3misc::map(names(data_list_flt), function(.data_name) {
  if (.data_name == "clinical") {
    data = data_list_flt[["clinical"]]
  } else {
    data = bind_cols(
      patient_id = clinical$patient_id,
      time = clinical$time,
      status = clinical$status,
      data_list_flt[[.data_name]]
    )
  }

  # create survival task
  task = as_task_surv(x = data, time = "time", event = "status", id = .data_name)
  # add "patient_id" as row names (not a feature)
  task$set_col_roles(cols = "patient_id", roles = "name")

  task
})
names(task_list) = names(data_list_flt)

# RESAMPLING FOR BENCHMARK ----
# 100 times => train/test split, stratified by censoring status and tumor stage (as we have just 9 patients with stage IV tumor)
clinical_task = task_list$clinical$clone()
clinical_task$set_col_roles(cols = c("status", "tumor_stage"),
                            add_to = "stratum")
ss = rsmp("subsampling", ratio = 0.8, repeats = 100)
set.seed(42)
ss$instantiate(clinical_task)

# METADATA ----
metadata = tibble(
  n_patients = task_list$clinical$nrow,
  n_modalities = length(task_list),
  n_gex_features = task_list$gex$n_features,
  n_cnv_features = task_list$cnv$n_features,
  n_prot_features = task_list$prot$n_features,
  n_phos_features = task_list$phos$n_features,
  n_ngly_features = task_list$ngly$n_features,
  n_clinical_features = task_list$clinical$n_features,
  n_total_features = sum(mlr3misc::map_int(task_list, \(.t) {.t$n_features})),
  n_events = sum(task_list$clinical$status() == 1),
  cens_rate = round(task_list$clinical$cens_prop(), digits = 2)
)
print(metadata)

# SAVE ALL DATA TO FILES ----
saveRDS(task_list, file = file.path(data_dir, "task_list.rds"))
saveRDS(ss, file = file.path(data_dir, "subsampling.rds"))

omic_ids = data.frame(omic_id = setdiff(names(task_list), "clinical"))
write_csv(omic_ids, file = file.path(data_dir, "omic_ids.csv"), col_names = FALSE)
write_csv(metadata, file = file.path(data_dir, "metadata.csv"))
write_csv(bind_cols(data_list_flt), file = file.path(data_dir, "all_data_preprocessed.csv"))
