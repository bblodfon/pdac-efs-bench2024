# Dataset description

## Publication

Wissel et al. (2023) Systematic comparison of multi-omics survival models reveals a widespread lack of noise resistance.

## Files

See `metadata.csv` for metadata info and `preprocess_wissel.R` for the preprocessing script.

Simple stratified train/test split was performed on **status** and **clinical tumor stage** to the extracted multi-modal dataset (`data_split.rds`).
Output `mlr3` tasks are in the `task_list.rds` and `all_data_preprocessed_*.csv` files.

## Data summary

6 data modalities/types are included:

- `Clinical` + surgical pathology: `age`, `sex`, `tumor_stage` (stages I and II)
- `GEX`: normalized counts, always > 0 (not standardized)
- `CNV`: values in {-2,-1,0,1,2}
- `RPPA`: normalized expression (not standardized)
- `Mutation`: number of non-silent mutations (integer counts), always >= 0 (up to 50)
- `Methylation`: beta-values in (0,1)

## Notes

- `miRNA` data was excluded due to non-sensical feature distributions.
