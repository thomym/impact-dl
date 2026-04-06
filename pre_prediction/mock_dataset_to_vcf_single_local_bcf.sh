#!/bin/bash

# -------------------------------
# Single Mock Dataset to VCF Pipeline
# Description: Processes ONE mock dataset to create VCF ready for Sei predictions.
# -------------------------------

set -euo pipefail

usage() {
    cat <<EOF
Usage: $0 -k <100genomes_bfile> -b <base_dir> -i <mock_index> [options]

Mandatory arguments:
  -k  Path prefix for the 1kg plink file
  -b  Base directory for outputs
  -i  Mock dataset index (e.g., 1, 2, 3, ...)

Optional arguments:
  -f  Fasta file for analysis
  -r  r² threshold for LD clumping (default: 0.69)
  -p  pvalue threshold for clumping (default: 5e-8)

EOF
    exit 1
}

# --- Default Configuration ---
RENAME_CHR_FILE="rename_chrs.txt"

# --- Parse Command-Line Arguments ---
while getopts ":k:b:i:r:p:f:" opt; do
    case ${opt} in
        k) KG_BFILE="${OPTARG}" ;;
        b) BASE_DIR="${OPTARG}" ;;
        i) MOCK_INDEX="${OPTARG}" ;;
        f) REF_FASTA="${OPTARG}" ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done
shift $((OPTIND - 1))

# --- Validate Mandatory Parameters ---
if [[ -z "${KG_BFILE:-}" || -z "${BASE_DIR:-}" || -z "${MOCK_INDEX:-}" ]]; then
    echo "ERROR: Missing required arguments: -k, -b, -i" >&2
    usage
fi

# --- Configuration Variables ---
NAME="mock_dataset${MOCK_INDEX}"
FILE="${BASE_DIR}/mock_datasets/${NAME}.tsv"
OUTPUT_DIR="${BASE_DIR}/vcfs_mock_datasets"
VEP_OUTPUT_DIR="${BASE_DIR}/vep_results_mock_datasets"
PRE_CLUMPING="${BASE_DIR}/pre_clumping_and_ld"
EUR_1KG_IDS="${PRE_CLUMPING}/eur_1kg_ids.txt"

# Create necessary directories
mkdir -p $OUTPUT_DIR/$NAME
mkdir -p $VEP_OUTPUT_DIR/$NAME

echo "🚀 Processing dataset: $FILE"

# Check if file exists
if [[ ! -f "$FILE" ]]; then
    echo "ERROR: Mock dataset file not found: $FILE" >&2
    exit 1
fi

# Step 1: Extract SNP list from mock clump file
awk -F'\t' 'NR > 1 && $6 != "" {print $3; gsub(",", "\n", $6); print $6}' $FILE | sort | uniq > $OUTPUT_DIR/$NAME/${NAME}_snp_list.txt

# Step 2: Extract SNPs
plink --bfile $KG_BFILE --extract $OUTPUT_DIR/$NAME/${NAME}_snp_list.txt \
    --keep $EUR_1KG_IDS --make-bed \
    --out $OUTPUT_DIR/$NAME/${NAME}_snps_from_clumps

# Step 3: Create VCF
plink2 --bfile $OUTPUT_DIR/$NAME/${NAME}_snps_from_clumps \
    --export vcf --fa $REF_FASTA --out $OUTPUT_DIR/$NAME/${NAME}_snps_correct_ref --ref-from-fa force


#TODO: provide VEP container
# Step 4: VEP Annotation
udocker run --volume=/home vep_container vep \
    -i $OUTPUT_DIR/$NAME/${NAME}_snps_correct_ref.vcf \
    -o $VEP_OUTPUT_DIR/$NAME/${NAME}_vep_output.txt \
    --species homo_sapiens --assembly GRCh37 --cache --offline \
    --dir_cache /opt/vep/.vep --fields Uploaded_variation,Location,Consequence,Gene,Feature --force_overwrite

# Step 5: Identify Coding SNPs
echo "Identifying Coding SNPs..."
awk -F'\t' 'NR > 1 { 
    is_protein_changing = 0;
    split($7, consequences, ",");
    for (i in consequences) {
        if (consequences[i] ~ /transcript_ablation|stop_gained|frameshift_variant|stop_lost|start_lost|transcript_amplification|inframe_insertion|inframe_deletion|missense_variant|protein_altering_variant/) {
            is_protein_changing = 1;
            break;
        }
    }
    if (is_protein_changing == 1) print $1;
}' $VEP_OUTPUT_DIR/$NAME/${NAME}_vep_output.txt | sort | uniq > $OUTPUT_DIR/$NAME/${NAME}_coding_snps_ids.txt

# Step 5b: Generate non-coding SNP list
grep -vFxf $OUTPUT_DIR/$NAME/${NAME}_coding_snps_ids.txt $OUTPUT_DIR/$NAME/${NAME}_snp_list.txt > $OUTPUT_DIR/$NAME/${NAME}_non_coding_snps_ids.txt

# Add debugging info
echo "All SNPs: $(wc -l < $OUTPUT_DIR/$NAME/${NAME}_snp_list.txt)"
echo "Coding SNPs: $(wc -l < $OUTPUT_DIR/$NAME/${NAME}_coding_snps_ids.txt)"
echo "Non-coding SNPs: $(wc -l < $OUTPUT_DIR/$NAME/${NAME}_non_coding_snps_ids.txt)"

#TODO: change like in main script - but here it is inside the filter_mock_datasets_vcfs.sh
echo "Filtering and preparing VCF file..."

echo "Renaming Chromosomes..."
bcftools annotate --rename-chrs $RENAME_CHR_FILE \
    -o $OUTPUT_DIR/$NAME/${NAME}_snps_correct_ref_fixed.vcf $OUTPUT_DIR/$NAME/${NAME}_snps_correct_ref.vcf

# Step 7: Extract Only Non-Coding SNPs (filter individual SNPs, not clumps)
echo "Filtering out coding SNPs..."
bcftools view -i "ID=@$OUTPUT_DIR/$NAME/${NAME}_non_coding_snps_ids.txt" \
    -o $OUTPUT_DIR/$NAME/${NAME}_snps_correct_ref_filtered.vcf $OUTPUT_DIR/$NAME/${NAME}_snps_correct_ref_fixed.vcf

# Step 8: Final Normalization
echo "Normalizing VCF..."
bcftools norm --check-ref s --fasta-ref $REF_FASTA \
    -o $OUTPUT_DIR/$NAME/${NAME}_snps_correct_ref_final.vcf $OUTPUT_DIR/$NAME/${NAME}_snps_correct_ref_filtered.vcf


echo "Mock dataset $NAME processing completed!"