# Dataset description

## Publication

Wissel et al. (2023) Systematic comparison of multi-omics survival models reveals a widespread lack of noise resistance.

## Files

- `metadata.csv`: metadata info
- `preprocess_wissel.R`: the preprocessing script
- `task_list.rds`: `mlr3` survival tasks, one per modality
- `all_data_preprocessed.csv`: preprocessed data from all patients and modalities

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
- We keep the highest variance 2000 features for modalities with more than 2000 features.
