#!/bin/bash

#SBATCH
# Check for required argument
if [ "$#" -ne 1 ]; then
    echo "Usage: sbatch sei_run.sh <base_dir>"
    exit 1
fi

# Base directory (from user argument)
base_dir="$1"

# Construct the VCF file path and output directory based on base_dir
vcf_file="${base_dir}/clumping_and_vcfs_outputs/main_gwas_snps_correct_ref_final.vcf"
output_dir="${base_dir}/sei_outputs_main_gwas"

# Ensure the output directory exists
mkdir -p "$output_dir"

# Print the paths for debugging
echo "Base Directory: $base_dir"
echo "VCF File: $vcf_file"
echo "Output Directory: $output_dir"

# Run Sei prediction
bash /path/to/your/sei/dir/sei-framework/1_variant_effect_prediction.sh "$vcf_file" hg19 "$output_dir" --cuda

if [ $? -eq 0 ]; then
    echo "Sei predictions completed successfully."
else
    echo "Sei predictions failed. Check logs for details." >&2
    exit 1
fi