# Dataset description

## Publication

Cao et al. (2021) Proteogenomic characterization of pancreatic ductal adenocarcinoma. [link](https://doi.org/10.1016/J.CELL.2021.08.023)

## Files

- `metadata.csv`: metadata info
- `preprocess_cao.R`: the preprocessing script
- `all_data_preprocessed.csv`: preprocessed data from all patients and data types
- `omic_ids.csv`: ids for all the omics preprocessed
- `task_list.rds`: `mlr3` survival tasks, one per data type
- `subsampling.rds`: `mlr3` subsampling object for the benchmark (100 subsamplings, 80/20 train/test set splits, stratified by status and tumor stage)

## Data summary

6 data modalities/types are included:

- `Clinical` + surgical pathology: `age`, `sex`, `tumor_stage` (stages I-IV), `tumor_size_cm`, `number of lymph nodes examined`, `number of positive lymph nodes`, `lymph_vascular_invasion` (2 classes)
- `GEX`: normalized counts, always > 0 (not standardized)
- `CNV`: Gene level copy number data, log2 ratio
- `Proteomics`: positive, median-normalized intensity values (not standardized)
- `Phosphoproteomics`: positive, median-normalized intensity values (not standardized)
- `N-Glycoproteomics`: expression values (TMT, Log2ratio, not standardized)

## Notes

- `miRNA`, `circRNA` and `methylation` data types were not used in this study.
We used all data types present in Fig 7A of the Cao et al. paper.
- We keep the highest variance 2000 features for modalities with more than 2000 features.
