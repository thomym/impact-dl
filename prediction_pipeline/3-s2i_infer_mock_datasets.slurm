#!/bin/bash

#SBATCH --job-name=sei_mock_datasets
if [ "$#" -ne 1 ]; then
    echo "Usage: sbatch sei_run.sh <base_dir>"
    exit 1
fi

# Base directory (from user argument)
BASE_DIR="$1"

# --- Environment Setup ---
SEI_OUTPUT_DIR="$BASE_DIR/sei_outputs_mock_datasets"
VCF_DIR="$BASE_DIR/vcfs_mock_datasets"

# Create output directory if it doesn't exist
echo "Creating output directory: $SEI_OUTPUT_DIR"
mkdir -p "$SEI_OUTPUT_DIR" || echo "Failed to create directory: $SEI_OUTPUT_DIR"

# --- Loop Through Datasets ---
START_INDEX=1  # Change as needed
END_INDEX=300  # Adjust this based on the number of datasets

for TASK_ID in $(seq $START_INDEX $END_INDEX); do
    # Define variables for this iteration
    GWAS_NAME="mock_dataset${TASK_ID}"
    VCF_FILE="$VCF_DIR/${GWAS_NAME}/${GWAS_NAME}_snps_correct_ref_final.vcf"
    SEI_OUTPUT="$SEI_OUTPUT_DIR/${GWAS_NAME}"

    echo "Running Sei predictions for: $GWAS_NAME"
    echo "Using VCF: $VCF_FILE"
    echo "Output Directory: $SEI_OUTPUT"

    # --- Run Sei Predictions ---
    bash /path/to/your/sei/dir/sei-framework/1_variant_effect_prediction.sh $VCF_FILE hg19 $SEI_OUTPUT --cuda

    # Check for errors and exit if a VCF file is missing
    if [ $? -ne 0 ]; then
        echo "Error occurred while processing $GWAS_NAME."
        echo "VCF File: $VCF_FILE"
        echo "Output Directory: $SEI_OUTPUT"
        echo "Checking VCF file existence..."
        if [ ! -f "$VCF_FILE" ]; then
            echo "VCF file is missing: $VCF_FILE"
        fi
        echo "Exiting loop."
        exit 1
    fi
    

    echo "Sei predictions completed for $GWAS_NAME"
done

# --- Completion Message ---
echo "All Sei predictions completed for datasets ${START_INDEX} to ${END_INDEX}"
