# pdac-efs-bench2024

This repository contains code and data for benchmarking multi-omics feature selection methods in **pancreatic ductal adenocarcinoma (PDAC)**, with a focus on survival prediction.

## Environment Setup

Restore R libraries required to run the benchmark using `renv`:
```
renv::restore(exclude = c("BiocGenerics", "BiocManager", "BiocVersion", "ComplexHeatmap", "IRanges", "S4Vectors", "MASS", "ipred", "class"))
```

## Data and Preprocessing

We use two multi-omics PDAC datasets. Preprocessing steps include:

- Patient and feature filtering
- Metadata curation
- Resampling definition for benchmark
- Creation of `mlr3` survival tasks per omic data type

All related code and metadata are found in the `data/` directory.

## Multi-omics Benchmark Overview

The benchmarking pipeline consists of three main steps, located in the `bench/` directory.

### Ensemble Feature Selection (EFS)

Script: [`bench/run_efs.sh`](https://github.com/bblodfon/pdac-efs-bench2024/blob/main/bench/run_efs.sh): Runs the Ensemble Feature Selection procedure ([efs.R](https://github.com/bblodfon/pdac-efs-bench2024/blob/main/bench/efs.R)) across:

- Both datasets
- All omics
- All resampling iterations

Stores results as [EnsembleFSResult](https://mlr3fselect.mlr-org.com/reference/ensemble_fs_result.html) objects in `bench/efs/`.
Contact the main author to share with you these intermediate results.

### Omic-wise Feature Selection

Script: `bench/run_fs.sh`: Performs feature selection per omic and per subsampling iteration (100 total).

Two methods available:
- Cox Lasso
- Pre-computed Ensemble Feature Selection (loaded from step 1). 
We automatically select the number of features via the Pareto front method.

Output: `bench/fs.rds` - a table with:

- `dataset_id`, `omic_id`, `rsmp_id`
- Selected features for each method

### Multi-omics Integration & Benchmarking

Script: `bench/run_mm_bench.R`: Combines selected features across omics (via late integration/fusion) per subsampling iteration, then trains and evaluates survival models on training/test splits.
Available models are Cox Proportional Hazards, Cox Lasso and Random Survival Forests (RSF).

Output: `bench/result.rds` - a table with:

- `dataset_id`, `fs_method_id`, `rsmp_id`, `model_data_config` (which model, which omics/clinical data config)
- Number of features for, task_feats
- Metrics such as Harrell's C-index, etc.
