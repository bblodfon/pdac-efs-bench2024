#' Extract PDAC multi-omics data from Osipov et al. (2024)
#' "The Molecular Twin artificial-intelligence platform integrates multi-omic data
#' to predict outcomes for pancreatic adenocarcinoma patients"
#' Execute: `Rscript data/osipov2024/preprocess_osipov.R`
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

#' Download "Source Data Fig. 1" from the paper webpage:
#' https://www.nature.com/articles/s43018-023-00697-7#Sec38
#' we converted the `.xlsx` file to `.csv` format
#' Authors mention in the "Methods" section that
#' 1) all samples are PDACs and 2) all patients Stage III and IV have been excluded

# DATA ----
dataset_path = file.path("data", "osipov2024")
csv_datafile = file.path(dataset_path, "43018_2023_697_MOESM3_ESM.csv")
all_data = readr::read_csv(file = csv_datafile, show_col_types = FALSE)
dim(all_data) # 74 x 6361

# CLINICAL AND SURGICAL PATHOLOGY ----
clinical = all_data |>
  select(starts_with(c("Biobank", "clinical", "label")))

clinical = clinical |>
  # rename variables
  rename(patient_id = Biobank_Number) |>
  rename(time = label_days_to_death) |>
  # convert days => months
  mutate(time = ceiling(time/30.44)) |>
  rename(status = label_deceased) |>
  rename(age = clinical_Age_at_Diagnosis) |>
  rename(sex = clinical_Sex_ord) |>
  rename(height = clinical_Height) |>
  rename(weight = clinical_Weight) |>
  rename(BMI = clinical_BMI) |>
  # convert types
  mutate(age = as.integer(age)) |> # years
  mutate(time = as.integer(time)) |> # months
  mutate(status = as.numeric(status)) |>
  # 3 classes ({0,1,2} => last has only 3 patients so we 'collapse' it to 1)
  rename(lymph_invasion = clinical_Lymphovascular_Invasion_ord) |>
  mutate(lymph_invasion = if_else(lymph_invasion == 0, 0, 1)) |>
  # 2 stages (I and II is my guess :)
  rename(tumor_stage = clinical_TNM_Mixed_Stage_ord) |>
  # 2 classes from the International Classification of Diseases for Oncology, Third Edition (ICD-O-3)
  # 81403 => Adenocarcinoma (encode as 0), 85003 => invasive (encode as 1)
  rename(histology_behavior = `clinical_Histology_Behavior_ICD-O-3_ord`) |>
  mutate(histology_behavior = if_else(histology_behavior == 81403, 0, 1)) |>
  # 3 classes (0,1,2)
  rename(grade = clinical_Grade_Mixed_ord) |>
  # 5 classes => 51 (0 => Head?), 23 (>0) - this is most probably the primary tumor site
  # i.e. body, head, tail, and combinations of these
  rename(clinical_site = `clinical_Site_-_Primary_ICD-O-3_ord`) |>
  # collapse to 2 classes
  mutate(clinical_site =  if_else(clinical_site == 0, 0, 1)) |>
  relocate(patient_id, time, status) # put the outcome columns in front

# (time, status) are okay
stopifnot(all(clinical$time > 0)) # no NAs or negative numbers here
stopifnot(clinical$status == 1 | clinical$status == 0) # dead (1) or censored (0)

# exclude some clinical variables
clinical = clinical |>
  # remove one-hot encoded variables (from a factor with many levels)
  select(-starts_with("clinical_secondary_diagnosis_onehot")) |>
  # remove treatment variables (potential target data leakage, though some of
  # these might be neo-adjuvant, i.e. before surgery)
  select(-c(clinical_Chemotherapy_Binary, # 2 classes
            clinical_Chemotherapy_Summary_ord, # 4 classes
            clinical_Radiation_Summary_ord)) |> # 2 classes
  # 3 classes unbalanced (most are zeros)
  select(-clinical_Surgical_Margins_Summary_ord) |>
  # 4 classes (I guess IA, IB, IIA and IIB stages) - remove as binary tumor_stage is enough
  select(-clinical_TNM_Mixed_Substage) |>
  # 2 classes, unbalanced (most are 1s)
  select(-clinical_Perineural_Invasion_ord) |>
  # 5 classes and most are in classes 0 or 1
  select(-clinical_Patient_History_Alcohol_ord) |>
  # very cancer-specific, this is not lung cancer!
  select(-clinical_Patient_History_Tobacco_ord) |>
  # exclude any history related variables
  select(-starts_with("clinical_Family_history")) |>
  select(-starts_with("clinical_Patient_History")) |>
  # More than 10 ethnicity classes and unbalanced (most belong to 1 class)
  select(-clinical_merged_ethnicity_ord)

# check missing clinical features => no missing features
mlr3misc::map(clinical, \(.x) sum(is.na(.x)))
stopifnot(nrow(clinical) == nrow(all_data))
clinical_complete_ids = clinical |> complete.cases() |> which()
clinical

# some useful stats
table(clinical$sex)
table(clinical$tumor_stage) # 17 => stage I (0), 57 => stage II (1)
table(clinical$status) # 24 => censored/alive, 50 dead

# OMICS (9) ----

## GEX ----
#' Extra preprocessing from the paper ("Methods" section):
#' "Finally, we trained our classifiers using log10 estimated read counts for
#' these 2,000 genes as features"
gex = all_data |>
  select(starts_with("rna_expr")) |>
  rename_with(function(.x) gsub(x = .x, pattern = "rna_expr", replacement = "gex")) |>
  rename_with(make.names) |>
  mutate(across(everything(), \(.x) log(.x + 1, base = 2))) # log2 read counts
gex_ids = gex |> complete.cases() |> which()
length(gex_ids) # 57
ncol(gex) # 2000
stopifnot(nrow(gex) == nrow(all_data))
counts = gex |> unlist()
all(counts >= 0, na.rm = TRUE) # normalized count data, always > 0

## FUSION GENES ----
fgex = all_data |>
  select(starts_with("AF4")) |>
  rename_with(function(.x) gsub(x = .x, pattern = "AF4", replacement = "fgex")) |>
  rename_with(make.names)

# remove constant features
feats_to_exclude = names(which(apply(fgex, 2, var, na.rm = TRUE) == 0))
length(feats_to_exclude) # 4
fgex = fgex |> select(-all_of(feats_to_exclude))

fgex_ids = fgex |> complete.cases() |> which() # same as gex
length(fgex_ids) # 57
ncol(fgex) # 25 features
stopifnot(nrow(fgex) == nrow(all_data))
fgex |> unlist() |> table() # values in {0,1,...,8}, mostly zeros

#' `indel`, `cnv`, `snv` extra preprocessing from the paper ("Methods" section):
#' 1) univariate normalization (~ N(0,1))
#' 2) pruning of low-variance features (with variance threshold < 0.05)
#' 3) remove highly correlated features (Spearman coef < 0.95)
#' "Processed genomic features consisted of:
#' - `337` somatic SNV
#' - `219` CNV
#' - `72` INDEL gene-level features
## SNV ----
snv = all_data |>
  select(starts_with("freebayes")) |>
  rename_with(function(.x) gsub(x = .x, pattern = "freebayes_SNV", replacement = "snv")) |>
  rename_with(make.names)

# remove constant features
feats_to_exclude = names(which(apply(snv, 2, var, na.rm = TRUE) == 0))
length(feats_to_exclude) # 40
snv = snv |> select(-all_of(feats_to_exclude))

# remove highly correlated features
feats_to_exclude = caret::findCorrelation(cor(snv, use = "complete.obs", method = "spearman"), cutoff = 0.95, names = TRUE)
length(feats_to_exclude) # 297
snv = snv |> select(-all_of(feats_to_exclude))
snv_ids = snv |> complete.cases() |> which()
length(snv_ids) # 72
ncol(snv) # 274 features (we remove more than they did in the paper!)
stopifnot(nrow(snv) == nrow(all_data))
snv |> unlist() |> table() # values in {0,1,...,6}, mostly zeros

## CNV ----
cnv = all_data |>
  select(starts_with('CNV_')) |>
  rename_with(function(.x) gsub(x = .x, pattern = "CNV", replacement = "cnv")) |>
  rename_with(make.names)

# kinda continuous features, with lots of zeros and centered
summary(mlr3misc::map_dbl(cnv, \(.x) mean(.x, na.rm = TRUE))) # mean ~ 0
summary(mlr3misc::map_dbl(cnv, \(.x) sd(.x, na.rm = TRUE))) # mean ~ 0.1

# variance filtering leaves only ~16 features so we don't do it
feats_flt = cnv |>
  summarise(across(everything(), \(.x) var(.x, na.rm = TRUE))) |>
  select(where(~ . > 0.05)) |> # variance threshold
  names()
length(feats_flt) # 16!

# remove highly correlated features
feats_to_exclude = caret::findCorrelation(cor(cnv, use = "complete.obs", method = "spearman"), cutoff = 0.95, names = TRUE)
length(feats_to_exclude) # 470
cnv = cnv |> select(-all_of(feats_to_exclude))
cnv_ids = cnv |> complete.cases() |> which()
length(cnv_ids) # 72
ncol(cnv) # 178 features (we remove more than they did in the paper!)
stopifnot(nrow(cnv) == nrow(all_data))

## INDELS ----
indel = all_data |>
  select(starts_with("pindel")) |>
  rename_with(function(.x) gsub(x = .x, pattern = "pindel_INDEL", replacement = "indel")) |>
  rename_with(make.names)

# remove constant features
feats_to_exclude = names(which(apply(indel, 2, var, na.rm = TRUE) == 0))
length(feats_to_exclude) # 9
indel = indel |> select(-all_of(feats_to_exclude))

# remove highly correlated features
feats_to_exclude = caret::findCorrelation(cor(indel, use = "complete.obs", method = "spearman"), cutoff = 0.95, names = TRUE)
length(feats_to_exclude) # 54
indel = indel |> select(-all_of(feats_to_exclude))
indel_ids = indel |> complete.cases() |> which()
length(indel_ids) # 72
ncol(indel) # 63 features (we remove more than they did in the paper!)
stopifnot(nrow(indel) == nrow(all_data))
indel |> unlist() |> table() # values in {0,1,2,3}, mostly zeros

## PATHOLOGY ----
#' Nuclear features using different algos on the WSI images + aggregated order
#' statistics, see `source_data_fig2.xlsx` file in the "Source Data" section of
#' the paper and Fig. 2 for the computational pathology pipeline.
#'
#' Extra preprocessing from the paper:
#' "The `z-scored` case-level features were used to develop ML models for survival prediction"
path = all_data |>
  select(starts_with("pathology")) |>
  rename_with(function(.x) gsub(x = .x, pattern = "pathology", replacement = "path")) |>
  rename_with(function(.x) gsub(x = .x, pattern = "%", replacement = "")) |>
  rename_with(make.names)

# remove constant features
feats_to_exclude = names(which(apply(path, 2, var, na.rm = TRUE) == 0))
length(feats_to_exclude) # 26 features
path = path |> select(-all_of(feats_to_exclude))

# standardize (z-score)
path = path |>
  mutate(across(everything(), scale)) |>
  mutate(across(everything(), as.numeric))

path_ids = path |> complete.cases() |> which()
length(path_ids) # 71
ncol(path) # 794
stopifnot(nrow(path) == nrow(all_data))

#' `lipids` and `protein` data
#' Extra preprocessing steps mentioned in the paper ("Methods" section):
#' 1) filtering out proteins and lipids with more than 25% missing data not meeting
#' quality control criteria
#' 2) removing proteins with a low variance <0.1 threshold
#' 3) followed by imputation of remaining missing values using halved median values
#' for each column
#' 4) univariate normalization of each column
#'
#' We will not perform the above preprocessing as the number of samples is smaller
#' than the ones with also DNA omics.
#'
## PLASMA PROTEIN ----
plasma_protein = all_data |>
  select(starts_with("plasma_protein")) # 257
plasma_protein_ids = plasma_protein |> complete.cases() |> which()
length(plasma_protein_ids) # 51
stopifnot(nrow(plasma_protein) == nrow(all_data))

## PLASMA LIPID ----
plasma_lipid = all_data |>
  select(starts_with("plasma_lipid")) # 406
plasma_lipid_ids = plasma_lipid |> complete.cases() |> which()
length(plasma_lipid_ids) # 51
stopifnot(nrow(plasma_lipid) == nrow(all_data))

## TISSUE PROTEIN ----
tissue_protein = all_data |>
  select(starts_with("tissue_protein")) # 1130
tissue_protein_ids = tissue_protein |> complete.cases() |> which()
length(tissue_protein_ids) # 49
stopifnot(nrow(tissue_protein) == nrow(all_data))

# CHOOSE MODALITIES ----
#' We want:
#' - max number of modalities
#' - max number of patients
#' - with complete data (no missing values) in the associated modalities
ids_list = list(clinical = clinical_complete_ids, gex = gex_ids, fgex = fgex_ids,
                snv = snv_ids, cnv = cnv_ids,indel = indel_ids, path = path_ids,
                plasma_protein = plasma_protein_ids, plasma_lipid = plasma_lipid_ids,
                tissue_protein = tissue_protein_ids)
res = usefun::powerset_icounts(ids_list)
# at least 3 modalities, more than 60 patients with all data
res |>
  filter(num_subsets >= 3, count > 60) |>
  arrange(desc(count), desc(num_subsets))

# make UpSet plot
m = make_comb_mat(ids_list, mode = "intersect")
m2 = m[comb_size(m) > 60 & comb_degree(m) >= 3] # filtering same as above
UpSet(m2, comb_order = order(comb_size(m2)))

# best combo to keep: clinical-snv-cnv-indel-path
ids_to_keep = res |>
  filter(set_combo == "clinical-snv-cnv-indel-path") |>
  pull(common_ids) |>
  unlist()

data_list = list(
  clinical = clinical[ids_to_keep, ],
  snv = snv[ids_to_keep, ],
  cnv = cnv[ids_to_keep, ],
  indel = indel[ids_to_keep, ],
  path = path[ids_to_keep, ]
)

# nrows the same
stopifnot(all(mlr3misc::map(data_list, nrow) == length(ids_to_keep)))

# no missing data in general
mlr3misc::map(data_list, function(.data) {
  .data |> as.matrix() |> as.vector() |> is.na() |> sum()
})

# TASKS ----
survival_outcome = data_list$clinical |> select(patient_id, time, status)

task_list = mlr3misc::map(names(data_list), function(.data_name) {
  if (.data_name == "clinical") {
    data = data_list[[.data_name]]
  } else {
    data = bind_cols(survival_outcome, data_list[[.data_name]])
  }

  task = as_task_surv(x = data, time = "time", event = "status", id = .data_name)

  # not a feature, add as row names
  task$set_col_roles(cols = "patient_id", roles = "name")

  task
})
names(task_list) = names(data_list)

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
  n_snv = task_list$snv$n_features,
  n_cnv = task_list$cnv$n_features,
  n_indel = task_list$indel$n_features,
  n_path = task_list$path$n_features,
  n_clinical = task_list$clinical$n_features,
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
write_csv(bind_cols(data_list), file = file.path(dataset_path, "all_data_preprocessed.csv"))
