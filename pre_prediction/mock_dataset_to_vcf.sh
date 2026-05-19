#!/bin/bash

# -------------------------------
# Single Mock Dataset to VCF Pipeline
# Description: Processes ONE mock dataset to create VCF ready for Sei predictions.
# -------------------------------

set -euo pipefail

usage() {
    cat <<EOF
Usage: $0 -k <1kg_bfile_prefix> -b <base_dir> -i <mock_index> -f <ref_fasta> -C <vep_cache> [options]

Mandatory arguments:
  -k  Path prefix for the filtered 1KG PLINK bfile.
  -b  Base directory: \${results_root}/\${gwas_name}.
  -i  Mock dataset index (1..n_mock_datasets).
  -f  Reference FASTA (typically \${sei_framework}/resources/hg19_UCSC.fa).
  -C  Host path to VEP cache, mounted as /opt/vep/.vep in the container.

Optional arguments:
  -c  Container runtime: udocker (default) or docker.
  -I  VEP container image (default: ensemblorg/ensembl-vep:release_113.0).

EOF
    exit 1
}

# --- Default Configuration ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RENAME_CHR_FILE="${SCRIPT_DIR}/rename_chrs.txt"
CONTAINER_RUNTIME=udocker
VEP_IMAGE=ensemblorg/ensembl-vep:release_113.0

# --- Parse Command-Line Arguments ---
while getopts ":k:b:i:r:p:f:c:I:C:" opt; do
    case ${opt} in
        k) KG_BFILE="${OPTARG}" ;;
        b) BASE_DIR="${OPTARG}" ;;
        i) MOCK_INDEX="${OPTARG}" ;;
        f) REF_FASTA="${OPTARG}" ;;
        c) CONTAINER_RUNTIME="${OPTARG}" ;;
        I) VEP_IMAGE="${OPTARG}" ;;
        C) VEP_CACHE="${OPTARG}" ;;
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
if [[ -z "${KG_BFILE:-}" || -z "${BASE_DIR:-}" || -z "${MOCK_INDEX:-}" || -z "${REF_FASTA:-}" || -z "${VEP_CACHE:-}" ]]; then
    echo "ERROR: Missing required arguments: -k, -b, -i, -f, -C" >&2
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

echo "Processing dataset: $FILE"

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


# Step 4: VEP Annotation
${CONTAINER_RUNTIME} run \
    -v "${BASE_DIR}:${BASE_DIR}" \
    -v "${VEP_CACHE}:/opt/vep/.vep" \
    "${VEP_IMAGE}" vep \
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