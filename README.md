# pdac-efs-bench2024

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18375793.svg)](https://doi.org/10.5281/zenodo.18375793)

This repository contains code and data for benchmarking multi-omics feature selection methods in **pancreatic ductal adenocarcinoma (PDAC)**, with a focus on survival prediction.

## Citation

> Zobolas, J., George, A.-M., López, A., Fischer, S., Becker, M., & Aittokallio, T. (2025). Optimizing Prognostic Biomarker Discovery in Pancreatic Cancer Through Hybrid Ensemble Feature Selection and Multi-Omics Data. https://arxiv.org/pdf/2509.02648

## Environment Setup

Restore R libraries required to run the benchmark using `renv`:
```
renv::restore(exclude = c("BiocGenerics", "BiocManager", "BiocVersion", "ComplexHeatmap", "IRanges", "S4Vectors", "MASS", "ipred", "class"))
```

**`R` packages versions** are available [here](https://github.com/bblodfon/pdac-efs-bench2024/blob/main/bench/efs.R#L40), if executing more manually is desired.

## Data and Preprocessing

We use 3 multi-omics PDAC datasets. Preprocessing steps include:

- Patient and feature pre-filtering
- Metadata curation
- Resampling definition for benchmark
- Creation of `mlr3` survival tasks per omic data type

All related code and metadata are found in the [`data/`](https://github.com/bblodfon/pdac-efs-bench2024/tree/main/data) directory.

## Multi-omics Benchmark

The benchmarking pipeline consists of three main steps, see also the `bench/` directory.

### Hybrid Ensemble Feature Selection (hEFS)

Script: [`bench/run_efs.sh`](https://github.com/bblodfon/pdac-efs-bench2024/blob/main/bench/run_efs.sh): Runs the Ensemble Feature Selection procedure ([`bench/efs.R`](https://github.com/bblodfon/pdac-efs-bench2024/blob/main/bench/efs.R)) across:

- All datasets
- All omics
- All resampling iterations

Stores results as [EnsembleFSResult](https://mlr3fselect.mlr-org.com/reference/ensemble_fs_result.html) objects in `bench/efs/`.
Contact the main author to share with you these intermediate results (this is the most time consuming part of this benchmark).

### Omic-wise Feature Selection

Script: [`bench/run_fs.sh`](https://github.com/bblodfon/pdac-efs-bench2024/blob/main/bench/run_fs.R): Performs feature selection per omic and per subsampling iteration (100 total).

Two methods available:
1. Cox Lasso
2. Pre-computed Ensemble Feature Selection (loaded from step 1).
We automatically select the number of features via the Pareto front method.

Output: [`bench/fs.rds`](https://github.com/bblodfon/pdac-efs-bench2024/blob/main/bench/fs.rds) - a table with columns `dataset_id`, `omic_id`, `rsmp_id` and the selected features for each method.

### Multi-omics Integration & Benchmarking

Script: [`bench/run_mm_bench.R`](https://github.com/bblodfon/pdac-efs-bench2024/blob/main/bench/run_mm_bench.R): Combines selected features across omics (via **late integration/fusion**) per subsampling iteration, then trains and evaluates survival models on the respective training/test splits.
Available models are Cox Proportional Hazards, Cox Lasso, Random Survival Forests and BlockForest.

Output: [`bench/result.rds`](https://github.com/bblodfon/pdac-efs-bench2024/blob/main/bench/result.rds) - a table with:

- `dataset_id`: Identifier of the dataset used
- `fs_method_id`: Feature selection method applied (or `no FS` if otherwise)
- `rsmp_id`: Subsampling (resampling) iteration identifier
- `model_data_config`: Configuration used for model training, indicating the model type and which omics and/or clinical data were included (`ALL` means clinical + all omics, `Clinical + GEX` means configuration used only clinical and gene expression data, etc.)
- `task_nfeats`: Number of selected features used in the dataset
- `task_feats`: The specific features that were selected
- Performance scores for the test sets (Harrell's C-index, Uno's AUC, etc.)

## Number of CV folds investigation

See [`cv_res/`](https://github.com/bblodfon/pdac-efs-bench2024/tree/main/cv_res).

