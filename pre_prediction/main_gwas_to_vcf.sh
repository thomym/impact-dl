#!/bin/bash
set -euo pipefail

# -------------------------------
# Main GWAS + 1KG to VCF Pipeline
# Description: Processes a GWAS dataset, filters using 1KG, and prepares VCF for Sei predictions.
# -------------------------------

# --- Usage ---

usage() {
    cat <<EOF
Usage: $0 -n <gwas_name> -b <base_dir> -f <ref_fasta> [options]

Mandatory arguments:
  -n  GWAS name (used to locate inputs/outputs under -b).
  -b  Base directory: \${results_root}/\${gwas_name}. The pre_clumping_and_ld/
      directory must already exist here (output of gen_pre_clumping.sh).
  -f  Reference FASTA (typically \${sei_framework}/resources/hg19_UCSC.fa).

Optional arguments:
  -g  GWAS file for clumping (default: \${base_dir}/pre_clumping_and_ld/gwas_for_clumping.txt).
  -r  r² threshold for LD clumping (default: 0.69).
  -p  P-value threshold for clumping (default: 5e-8).
  -c  Container runtime: udocker (default) or docker.
  -I  VEP container image (default: ensemblorg/ensembl-vep:release_113.0).
  -C  Host path to VEP cache, mounted as /opt/vep/.vep in the container.
      Required for the VEP annotation step.

EOF
    exit 1
}



# --- Default Configuration ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RENAME_CHR_FILE="${SCRIPT_DIR}/rename_chrs.txt"
PVAL_THRESHOLD=5e-8
LD_R2_THRESHOLD=0.69
CONTAINER_RUNTIME=udocker
VEP_IMAGE=ensemblorg/ensembl-vep:release_113.0


# --- Parse Command-Line Arguments ---
while getopts ":g:b:n:r:p:f:c:I:C:" opt; do
    case ${opt} in
        g) GWAS_FILE="${OPTARG}" ;;
        b) BASE_DIR="${OPTARG}" ;;
        n) GWAS_NAME="${OPTARG}" ;;
        r) LD_R2_THRESHOLD="${OPTARG}" ;;
        p) PVAL_THRESHOLD="${OPTARG}" ;;
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
if [[ -z "${GWAS_NAME:-}" || -z "${BASE_DIR:-}" || -z "${REF_FASTA:-}" || -z "${VEP_CACHE:-}" ]]; then
    echo "ERROR: Missing required arguments: -n <gwas_name>, -b <base_dir>, -f <ref_fasta>, -C <vep_cache>." >&2
    usage
fi


KG_BFILE="$BASE_DIR/pre_clumping_and_ld/${GWAS_NAME}_FILTERED_1kg/ds"
VEP_OUTPUT_DIR="$BASE_DIR/vep_results_main_gwas"
PRE_CLUMPING="$BASE_DIR/pre_clumping_and_ld"
EUR_1KG_IDS="$PRE_CLUMPING/eur_1kg_ids.txt"
OUTPUT_DIR="$BASE_DIR/clumping_and_vcfs_outputs"
# Set GWAS_FILE to the provided path or default to the pre-clumping path
GWAS_FILE="${GWAS_FILE:-$PRE_CLUMPING/gwas_for_clumping.txt}"
mkdir -p $OUTPUT_DIR $VEP_OUTPUT_DIR

echo "Processing Main GWAS: $GWAS_FILE"

# --- Step 1: LD Clumping ---
plink --bfile $KG_BFILE \
 --keep $EUR_1KG_IDS \
 --clump $GWAS_FILE \
 --clump-p1 $PVAL_THRESHOLD --clump-r2 $LD_R2_THRESHOLD \
 --clump-kb 250 --out $OUTPUT_DIR/clumps_original_1kg

# Clean clumps file by:
# 1. Removing any occurrences of "(number)" patterns.
# 2. Replacing empty parentheses followed by whitespace with a comma.
# Clean the clumps to be in the plink2 format - no (1) after each snp in each clump
sed 's/([0-9]*)//g; s/()[[:space:]]/,/g; s/()//g' $OUTPUT_DIR/clumps_original_1kg.clumped > $OUTPUT_DIR/clumps_cleaned.clumped

# Create GWAS snp list - take snps from SNP column and from SP2 col, which contains the clumped snps
awk 'NR > 1 {print $3; gsub(",", "\n", $12); print $12}' $OUTPUT_DIR/clumps_cleaned.clumped | sort | uniq > $OUTPUT_DIR/main_gwas_snp_list.txt

# Create a SNP-to-clump mapping file
awk 'NR > 1 {
    clump_id = $3;  # Index SNP is the clump ID
    print clump_id "\t" clump_id;  # Index SNP belongs to its own clump
    
    # Only process SP2 if it contains actual SNPs (not NONE)
    if ($12 != "NONE") {
        split($12, snps, ",");
        for (i in snps) {
            if (snps[i] != "") print clump_id "\t" snps[i];  # Map each clumped SNP to its clump
        }
    }
}' $OUTPUT_DIR/clumps_cleaned.clumped > $OUTPUT_DIR/snp_to_clump_mapping.txt

# --- Step 2: Extract SNPs from 1KG ---
echo "Step 2: Extract SNPs from 1KG..."
plink --bfile $KG_BFILE --extract $OUTPUT_DIR/main_gwas_snp_list.txt \
    --keep $EUR_1KG_IDS --make-bed \
    --out $OUTPUT_DIR/main_gwas_snps_from_clumps

# --- Step 3: Create VCF with Correct REF Allele ---
echo "Step 3: Create VCF with Correct REF Allele..."
plink2 --bfile $OUTPUT_DIR/main_gwas_snps_from_clumps \
    --export vcf --fa $REF_FASTA --out $OUTPUT_DIR/main_gwas_snps_correct_ref --ref-from-fa force

# --- Step 4: VEP Annotation ---
echo "Step 4: Annotating with VEP..."
${CONTAINER_RUNTIME} run \
    -v "${BASE_DIR}:${BASE_DIR}" \
    -v "${VEP_CACHE}:/opt/vep/.vep" \
    "${VEP_IMAGE}" vep \
    -i $OUTPUT_DIR/main_gwas_snps_correct_ref.vcf \
    -o $VEP_OUTPUT_DIR/main_gwas_vep_output.txt \
    --species homo_sapiens --assembly GRCh37 --cache --offline \
    --dir_cache /opt/vep/.vep --fields Uploaded_variation,Location,Consequence,Gene,Feature --force_overwrite

# --- Step 5: Identify Coding (Protein-Changing) SNPs ---
echo "Step 5: Identifying Coding (Protein-Changing) SNPs..."
awk -F'\t' 'NR > 1 { 
    is_protein_changing = 0;
    split($7, consequences, ",");
    for (i in consequences) {
        # Check for protein-changing consequence types
        if (consequences[i] ~ /transcript_ablation|stop_gained|frameshift_variant|stop_lost|start_lost|transcript_amplification|inframe_insertion|inframe_deletion|missense_variant|protein_altering_variant/) {
            is_protein_changing = 1;
            break;
        }
    }
    if (is_protein_changing == 1) print $1;
}' $VEP_OUTPUT_DIR/main_gwas_vep_output.txt | sort | uniq > $OUTPUT_DIR/main_gwas_coding_snps_ids.txt

# --- Step 6: Identify Clumps Containing Coding SNPs ---
echo "Step 6: Identifying Clumps with Coding SNPs..."
awk 'NR==FNR {coding[$1]=1; next} 
    $2 in coding {print $1}' \
    $OUTPUT_DIR/main_gwas_coding_snps_ids.txt $OUTPUT_DIR/snp_to_clump_mapping.txt | sort | uniq > $OUTPUT_DIR/clumps_with_coding_snps.txt



# --- Step 6b: Create Filtered Clumps File ---
echo "Step 6b: Creating filtered clumps file (non-coding clumps only)..."
awk '
    # First file - build exclusion list
    NR==FNR {
        exclude_clump[$1]=1
        next
    }
    # Second file - header line
    FNR==1 {
        print
        next
    }
    # Second file - data lines, SNP is in column 3
    {
        if (!($3 in exclude_clump)) print
    }
' $OUTPUT_DIR/clumps_with_coding_snps.txt $OUTPUT_DIR/clumps_cleaned.clumped > $OUTPUT_DIR/clumps_cleaned_noncoding.clumped

# Count clumps before and after filtering
echo "Original clumps: $(grep -v '^CHR' $OUTPUT_DIR/clumps_cleaned.clumped | wc -l)"
echo "Non-coding clumps: $(grep -v '^CHR' $OUTPUT_DIR/clumps_cleaned_noncoding.clumped | wc -l)"


# --- Step 7: Extract SNPs from Non-Coding Clumps Only ---
echo "Step 7: Extracting SNPs from Non-Coding Clumps..."
awk 'NR==FNR {exclude_clump[$1]=1; next}
    !($1 in exclude_clump) {print $2}' \
    $OUTPUT_DIR/clumps_with_coding_snps.txt $OUTPUT_DIR/snp_to_clump_mapping.txt | sort | uniq > $OUTPUT_DIR/snps_from_noncoding_clumps.txt

# --- Step 8: Run the fix_vcf script on filtered VCF ---
echo "Step 8: Running filter_main_gwas_vcf.sh on filtered VCF..."
echo "Renaming Chromosomes for Reference Match..."
bcftools annotate --rename-chrs $RENAME_CHR_FILE -o $OUTPUT_DIR/main_gwas_snps_correct_ref_fixed.vcf $OUTPUT_DIR/main_gwas_snps_correct_ref.vcf

# --- Step 7: Filter Non-Coding SNPs ---
echo "Extracting Non-Coding SNPs..."
bcftools view -i ID=@$OUTPUT_DIR/snps_from_noncoding_clumps.txt -o $OUTPUT_DIR/main_gwas_snps_correct_ref_filtered.vcf $OUTPUT_DIR/main_gwas_snps_correct_ref_fixed.vcf


# --- Step 8: Validation ---
echo "Validating Final VCF File..."
bcftools norm --check-ref s --fasta-ref $REF_FASTA $OUTPUT_DIR/main_gwas_snps_correct_ref_filtered.vcf -o $OUTPUT_DIR/main_gwas_snps_correct_ref_final.vcf

echo "Pipeline Complete: Final VCF ready for Sei predictions."