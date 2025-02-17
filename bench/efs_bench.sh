#!/bin/bash

# Execute: bash run_benchmark.sh

# Define paths
DATA_DIR="data"
SCRIPT="Rscript bench/run_efs.R"

# Define dataset IDs manually or leave empty to process all in data/
DATASET_IDS=()  # Example: ("wissel2023" "another_dataset") or leave empty to process all

# Define omic IDs manually or leave empty to read from omic_ids.csv
OMIC_IDS=()  # Example: ("gex" "cnv")

# Define rsmp_id range
RSMP_IDS=({1..100})  # Hard set to 1:100 or choose: `=42`, `=(1 42 99)`

# Get all dataset IDs if not specified
if [[ ${#DATASET_IDS[@]} -eq 0 ]]; then
    DATASET_IDS=($(basename -a "$DATA_DIR"/*/))
fi

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
            echo "Running: $SCRIPT $dataset_id $omic_id $rsmp_id"
            $SCRIPT "$dataset_id" "$omic_id" "$rsmp_id"
        done
    done
done
