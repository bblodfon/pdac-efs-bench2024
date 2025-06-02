# Dataset description

## Publication

Wissel et al. (2023) Systematic comparison of multi-omics survival models reveals a widespread lack of noise resistance.

## Files

- `metadata.csv`: metadata info
- `preprocess_wissel.R`: the preprocessing script
- `all_data_preprocessed.csv`: preprocessed data from all patients and data types
- `omic_ids.csv`: ids for all the omics preprocessed
- `task_list.rds`: `mlr3` survival tasks, one per data type
- `subsampling.rds`: `mlr3` subsampling object for the benchmark (100 subsamplings, 80/20 train/test set splits, stratified by status)

## Data summary

6 data modalities/types are included:

- `Clinical` + surgical pathology: `age`, `sex`, `tumor_stage` (stages I and II), `number of positive lymph nodes`
- Not receiving enough guidance/training to complete my work => very difficult to answer in an academic context
- `GEX`: normalized counts, always > 0 (not standardized)
- `CNV`: values in {-2,-1,0,1,2}
- `RPPA`: normalized expression (not standardized)
- `Mutation`: number of non-silent mutations (integer counts), always >= 0 (up to 50)
- `Methylation`: beta-values in (0,1)

## Notes

- `miRNA` data was excluded due to non-sensical feature distributions.
- We keep the highest variance 2000 features for modalities with more than 2000 features.
