#!/bin/bash
set -euo pipefail

# The usage() function prints a help message and then exits.
# It explains which parameters are mandatory and which are optional.
usage() {
    cat <<EOF
Usage: $0 -g <gwas_file> -p <pop_panel_file> -k <1kg_bfile_prefix> -o <output_dir> [options]

Mandatory arguments:
  -g  Path to the GWAS file. This file should be tab-separated with a header line.
      Expected column order (from left to right): 
         SNP, CHR, BP, A1, A2, MAF, SE, P, N, INFO, OR.
      For clumping, the script uses only the following columns:
         SNP (SNP identifier), CHR (chromosome), BP (base pair position), and P (p-value).
  -b  Base directory for ouyputs. in it a pre_clumping_and_ld directory will be created, and all outputs will be written to there
        other then genotyping files which are written to dec directories
  -n  GWAS name - for building directories 

  
Optional arguments:
  -k  Path prefix for the 1kg PLINK files (used with --bfile).

  -p  Path to the population panel file (e.g., pop.panel).

  -c  Comma-separated list of column numbers to extract from the GWAS file for clumping.
      Default is "1,2,3,8" corresponding to: SNP, CHR, BP, and P-value.
      Adjust if your GWAS file has a different column order.
  -m  Path to the gnomAD lookup file.

  -r  r² threshold for LD clumping.
      Default: 0.69

EOF
    exit 1
}

# Set default values for optional arguments.
CLUMP_COLS="1,2,3,8"
R2_THRESHOLD="0.69"


# Process command-line options.
while getopts ":g:p:k:b:c:m:r:n:" opt; do
    case ${opt} in
        g)
            GWAS_FILE="${OPTARG}"
            ;;
        p)
            POP_PANEL="${OPTARG}"
            ;;
        k)
            KG_BFILE="${OPTARG}"
            ;;
        b)
            BASE_DIR="${OPTARG}"
            ;;
        c)
            CLUMP_COLS="${OPTARG}"
            ;;
        m)
            GNOMAD_FILE="${OPTARG}"
            ;;
        r)
            R2_THRESHOLD="${OPTARG}"
            ;;
        n)
            GWAS_NAME="${OPTARG}"
            ;;
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

# Check that all mandatory parameters have been provided.
if [[ -z "${GWAS_FILE:-}" || -z "${BASE_DIR:-}" || -z "${GWAS_NAME:-}" ]]; then
    echo "ERROR: Missing one or more mandatory arguments." >&2
    usage
fi

# Create the output directory if it doesn't exist.
mkdir -p "$BASE_DIR"
OUTPUT_DIR="$BASE_DIR/pre_clumping_and_ld"
mkdir -p "$OUTPUT_DIR"


# --- Determine the location for the filtered 1kg file ---

# KG_BFILE is a PLINK prefix, e.g., "/path/to/default/1kg/original/ds"
# We want to create a new folder inside the parent of "original".
# 1. Get the directory that holds the KG_BFILE prefix.
#KG_DIR=$(dirname "$KG_BFILE")    # This would be .../1kg/original

# 2. Get the parent directory of KG_DIR.
#KG_PARENT=$(dirname "$KG_DIR")   # This should be .../1kg

# 3. Create a new folder named "<GWAS_NAME>_FILTERED_1kg" inside KG_PARENT.
#KG_NEW_FOLDER="${KG_PARENT}/${GWAS_NAME}_FILTERED_1kg"
#mkdir -p "$KG_NEW_FOLDER"

# 4. Define the new PLINK file prefix (we’ll keep the same file name as before, e.g., "ds").
#KG_OUTPUT_PREFIX="${KG_NEW_FOLDER}/ds"

KG_NEW_FOLDER="${OUTPUT_DIR}/${GWAS_NAME}_FILTERED_1kg"
mkdir -p "$KG_NEW_FOLDER"
KG_OUTPUT_PREFIX="${KG_NEW_FOLDER}/ds"


echo "Running pipeline with the following parameters:"
echo "GWAS name: $GWAS_NAME"
echo "GWAS file: $GWAS_FILE"
echo "Population panel: $POP_PANEL"
echo "1kg PLINK prefix: $KG_BFILE"
echo "Base Directory: $BASE_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Clump columns: $CLUMP_COLS (Expected order: SNP, CHR, BP, P)"
echo "gnomAD file: $GNOMAD_FILE"
echo "r² threshold: $R2_THRESHOLD"

#############################
# STEP 1: Create EUR IDs list from the population panel
#############################
grep EUR "$POP_PANEL" | awk '{print $1, $2}' > "$OUTPUT_DIR/eur_1kg_ids.txt"

#############################
# STEP 2: Extract list of GWAS SNPs
#############################
# Assumes the first column contains SNP IDs.
awk 'NR > 1 {print $1}' "$GWAS_FILE" > "$OUTPUT_DIR/full_gwas_snps.txt"

#############################
# STEP 3: Filter 1kg data for GWAS SNPs using PLINK
#############################

plink --bfile "$KG_BFILE" \
      --extract "$OUTPUT_DIR/full_gwas_snps.txt" \
      --make-bed \
      --out "$KG_OUTPUT_PREFIX"

#############################
# STEP 4: Filter GWAS file to include only SNPs present in the 1kg data
#############################
# Get the SNP list from the filtered 1kg data.
plink --bfile "$KG_OUTPUT_PREFIX" --write-snplist --out "$OUTPUT_DIR/kg_filtered_snplist"
# Write header to output file.
head -n 1 "$GWAS_FILE" > "$OUTPUT_DIR/gwas_filtered_by_1kg.txt"
# Append matching SNPs from the original GWAS file.
grep -wFf "$OUTPUT_DIR/kg_filtered_snplist.snplist" "$GWAS_FILE" >> "$OUTPUT_DIR/gwas_filtered_by_1kg.txt"

#############################
# STEP 5: Prepare GWAS file for clumping
#############################
# Here we extract the columns specified by CLUMP_COLS (default: 1,2,3,8)
# which correspond to SNP, CHR, BP, and P.
IFS=',' read -r col1 col2 col3 col4 <<< "$CLUMP_COLS"
awk_command="{print \$$col1\"\t\"\$$col2\"\t\"\$$col3\"\t\"\$$col4}"
echo "awk command creating gwas for clumping: awk $awk_command"
awk "$awk_command" "$OUTPUT_DIR/gwas_filtered_by_1kg.txt" > "$OUTPUT_DIR/gwas_for_clumping.txt"

# Create a list of SNPs (assumed to be in the first column) for LD calculation.
awk 'NR > 1 {print $1}' "$OUTPUT_DIR/gwas_for_clumping.txt" > "$OUTPUT_DIR/snp_list_for_r2.txt"

#############################
# STEP 6: Compute LD (r²) using PLINK
#############################
plink --bfile "$KG_OUTPUT_PREFIX" \
      --keep "$OUTPUT_DIR/eur_1kg_ids.txt" \
      --r2 yes-really \
      --ld-snp-list "$OUTPUT_DIR/snp_list_for_r2.txt" \
      --ld-window 99999 \
      --ld-window-kb 250 \
      --ld-window-r2 "$R2_THRESHOLD" \
      --out "$OUTPUT_DIR/r2_filtered_snps"

#############################
# STEP 7: Lookup MAF information from gnomAD
#############################
grep -wFf "$OUTPUT_DIR/snp_list_for_r2.txt" "$GNOMAD_FILE" > "$OUTPUT_DIR/maf_lookup.tsv"
header="CHROM\tPOS\tID\tREF\tALT\tMAF"
echo -e "$header" | cat - "$OUTPUT_DIR/maf_lookup.tsv" > "$OUTPUT_DIR/maf_lookup_with_header.tsv"

echo "Pipeline complete. The file for clumping is located at: $OUTPUT_DIR/gwas_for_clumping.txt"
