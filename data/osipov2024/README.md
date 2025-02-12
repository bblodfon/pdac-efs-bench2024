# Dataset description

## Publication

Osipov et al. (2024) The Molecular Twin artificial-intelligence platform integrates multi-omic data to predict outcomes for pancreatic adenocarcinoma patients

## Files

- `metadata.csv`: metadata info
- `preprocess_osipov.R`: the preprocessing script
- `all_data_preprocessed.csv`: preprocessed data from all patients and modalities
- `task_list.rds`: `mlr3` survival tasks, one per modality
- `subsampling.rds`: mlr3 subsampling object for the benchmark (100 subsamplings, 80/20 train/test set splits, stratified by status)

## Data summary

5 data modalities/types are included:

- `Clinical` + surgical pathology: `age`, `height`, `weight`, `BMI`, `sex`, `clinical_site` (2 classes: 0 or >0), `histology_behavior` (2 classes), `tumor_stage` (stage I and II), `grade` (3 classes: 0, 1 and 2), `lymph_invasion` (2 classes)
- `SNV`: values in {0,1,...,6}, mostly zeros
- `CNV`: continuous features, with lots of zeros and centered (mean ~ 0)
- `INDEL`: values in {0,1,2,3}, mostly zeros
- `Pathology`: standardized (z-scored) continuous features

## Notes

- The multi-modal dataset from the paper includes other modalities as well, but we chose the **maximum number of modalities (>=3)** for which we have **> 60 patients** with **complete data** across these modalities.
- The survival data is administratively censored and we use the dataset with a right-censored target outcome (`time`, `status`).
