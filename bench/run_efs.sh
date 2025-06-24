#!/bin/bash
#' This script runs the ensemble feature selection method on a specific
#' set of datasets, omics and resampling ids.

#' Execute: `bash bench/run_efs.sh` (from project root)

# Define paths
DATA_DIR="data"
EFS_SCRIPT="Rscript bench/efs.R"

# Define dataset IDs manually or leave empty to process all in data/
DATASET_IDS=()  # Example: ("wissel2023") or ("osipov2024") or ("cao2021") or leave empty to process all

# Define omic IDs manually or leave empty to read from omic_ids.csv
OMIC_IDS=()  # Example: ("gex" "cnv")

# Define rsmp_id range
RSMP_IDS=({1..100})  # Hard set to 1:100 or choose: `=42`, `=(1 42 99)`

# Get all dataset IDs if not specified
if [[ ${#DATASET_IDS[@]} -eq 0 ]]; then
    DATASET_IDS=($(basename -a "$DATA_DIR"/*/))
fi

# Calculate total number of iterations for progress tracking
TOTAL_RUNS=0
for dataset_id in "${DATASET_IDS[@]}"; do
    dataset_path="$DATA_DIR/$dataset_id"

    omic_file="$dataset_path/omic_ids.csv"
    if [[ -f "$omic_file" ]]; then
        if [[ ${#OMIC_IDS[@]} -eq 0 ]]; then
            mapfile -t omic_ids < "$omic_file"
        else
            omic_ids=("${OMIC_IDS[@]}")
        fi
        for omic_id in "${omic_ids[@]}"; do
            for rsmp_id in "${RSMP_IDS[@]}"; do
                ((TOTAL_RUNS++))
            done
        done
    fi
done

# Track progress
CURRENT_RUN=0

# Loop over datasets
for dataset_id in "${DATASET_IDS[@]}"; do
    dataset_path="$DATA_DIR/$dataset_id"

    # Ensure omic_ids.csv exists
    omic_file="$dataset_path/omic_ids.csv"
    if [[ ! -f "$omic_file" ]]; then
        echo "Warning: Missing $omic_file, skipping dataset $dataset_id"
        continue
    fi

    # Read available omic IDs from omic_ids.csv only if OMIC_IDS is not pre-defined
    if [[ ${#OMIC_IDS[@]} -eq 0 ]]; then
        mapfile -t omic_ids < "$omic_file"
    else
        omic_ids=("${OMIC_IDS[@]}")
    fi

    # Loop over omic IDs
    for omic_id in "${omic_ids[@]}"; do
        # Loop over rsmp_id values
        for rsmp_id in "${RSMP_IDS[@]}"; do
            ((CURRENT_RUN++))
            PERCENT_DONE=$((CURRENT_RUN * 100 / TOTAL_RUNS))
            echo "[$PERCENT_DONE%] Running: $EFS_SCRIPT $dataset_id $omic_id $rsmp_id"
            $EFS_SCRIPT "$dataset_id" "$omic_id" "$rsmp_id"
        done
    done
done
