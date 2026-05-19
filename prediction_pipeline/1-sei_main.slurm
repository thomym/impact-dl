#!/bin/bash
#SBATCH --job-name=sei_main
#SBATCH --constraint="l40s"
#SBATCH --cpus-per-task=4
#SBATCH --mem=50000
#SBATCH --gres=gpu:1
#
# NOTE: --partition is intentionally NOT set here. Most users will want to pass
# their own (e.g. --partition=gpu-shamir). Submit as:
#   sbatch --partition=<your_partition> 1-sei_main.slurm --gwas_name <name>
# If you submit without --partition, your cluster's default partition will be used.

# Run Sei on the real (clumped) GWAS VCF.
#
# Usage (two equivalent forms):
#
#   # 1. Submit from inside the impact_dl repo — $SLURM_SUBMIT_DIR locates paths.py
#   cd /path/to/impact_dl
#   sbatch [--partition=...] prediction_pipeline/1-sei_main.slurm --gwas_name <name>
#
#   # 2. Submit from anywhere with --repo_root explicit
#   sbatch [--partition=...] /path/to/impact_dl/prediction_pipeline/1-sei_main.slurm \
#       --gwas_name <name> --repo_root /path/to/impact_dl
#
# Other overrides:
#   --paths_yaml /alt/paths.yaml      # override yaml file
#   --sei_dir /alt/sei-framework      # override sei_framework path

set -e

# --- Parse args ---
GWAS_NAME=""
PATHS_YAML=""
SEI_DIR=""
REPO_ROOT_ARG=""
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --gwas_name)  GWAS_NAME="$2"; shift 2 ;;
        --paths_yaml) PATHS_YAML="$2"; shift 2 ;;
        --sei_dir)    SEI_DIR="$2"; shift 2 ;;
        --repo_root)  REPO_ROOT_ARG="$2"; shift 2 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$GWAS_NAME" ]]; then
    echo "ERROR: --gwas_name is required" >&2
    exit 1
fi

# --- Locate the impact_dl repo (must contain paths.py) ---
# Priority: --repo_root flag → $SLURM_SUBMIT_DIR (or its parent) → script dir (non-SLURM only).
# BASH_SOURCE is unreliable under SLURM because the script is copied to /var/spool/slurmd.
REPO_ROOT=""
if [[ -n "$REPO_ROOT_ARG" ]]; then
    REPO_ROOT="$REPO_ROOT_ARG"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "$SLURM_SUBMIT_DIR/paths.py" ]]; then
    REPO_ROOT="$SLURM_SUBMIT_DIR"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "$SLURM_SUBMIT_DIR/../paths.py" ]]; then
    REPO_ROOT="$(cd "$SLURM_SUBMIT_DIR/.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)" || true
    if [[ -n "$SCRIPT_DIR" && -f "$SCRIPT_DIR/../paths.py" ]]; then
        REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
    fi
fi

if [[ -z "$REPO_ROOT" || ! -f "$REPO_ROOT/paths.py" ]]; then
    echo "ERROR: cannot locate impact_dl/paths.py. Either:" >&2
    echo "  - cd into impact_dl/ before sbatch (so SLURM_SUBMIT_DIR points at the repo), or" >&2
    echo "  - pass --repo_root /path/to/impact_dl as a script arg." >&2
    exit 1
fi

# --- Resolve paths from paths.yaml ---
paths_args=""
if [[ -n "$PATHS_YAML" ]]; then paths_args="--paths_yaml $PATHS_YAML"; fi
RESULTS_ROOT=$(python "$REPO_ROOT/paths.py" $paths_args --get results_root)
if [[ -z "$SEI_DIR" ]]; then
    SEI_DIR=$(python "$REPO_ROOT/paths.py" $paths_args --get sei_framework)
fi

BASE_DIR="${RESULTS_ROOT}/${GWAS_NAME}"
VCF_FILE="${BASE_DIR}/clumping_and_vcfs_outputs/main_gwas_snps_correct_ref_final.vcf"
OUTPUT_DIR="${BASE_DIR}/sei_outputs_main_gwas"
H5_DIR="${OUTPUT_DIR}/chromatin-profiles-hdf5"
DIFFS_FILE="${H5_DIR}/main_gwas_snps_correct_ref_final_diffs.h5"

# --- Redirect stdout/stderr into results dir ---
mkdir -p "${BASE_DIR}/slurm_logs"
exec > "${BASE_DIR}/slurm_logs/sei_main_${SLURM_JOB_ID:-local}.out" \
     2> "${BASE_DIR}/slurm_logs/sei_main_${SLURM_JOB_ID:-local}.err"

# Clean up SLURM's default output file in SLURM_SUBMIT_DIR on exit.
# Non-empty (= early errors happened): move to slurm_logs/ to preserve.
# Empty (= success / nothing wrote before exec): delete.
if [[ -n "${SLURM_SUBMIT_DIR:-}" && -n "${SLURM_JOB_ID:-}" ]]; then
    SLURM_DEFAULT_OUT="${SLURM_SUBMIT_DIR}/slurm-${SLURM_JOB_ID}.out"
    trap '
        if [[ -f "$SLURM_DEFAULT_OUT" ]]; then
            if [[ -s "$SLURM_DEFAULT_OUT" ]]; then
                mv -f "$SLURM_DEFAULT_OUT" "${BASE_DIR}/slurm_logs/" 2>/dev/null || true
            else
                rm -f "$SLURM_DEFAULT_OUT" 2>/dev/null || true
            fi
        fi
    ' EXIT
fi

echo "REPO_ROOT:     $REPO_ROOT"
echo "GWAS:          $GWAS_NAME"
echo "Partition:     ${SLURM_JOB_PARTITION:-<default>}"
echo "Base dir:      $BASE_DIR"
echo "VCF file:      $VCF_FILE"
echo "Sei dir:       $SEI_DIR"
echo "Output dir:    $OUTPUT_DIR"

if [[ -z "${SLURM_JOB_PARTITION:-}" ]]; then
    echo "WARNING: no SLURM partition set; falling back to cluster default."
fi

# --- Skip if already complete ---
if [[ -f "$DIFFS_FILE" ]]; then
    echo "Sei diffs already exist at ${DIFFS_FILE}; skipping."
    exit 0
fi

# --- Validate inputs ---
if [[ ! -f "$VCF_FILE" ]]; then
    echo "ERROR: VCF file not found: $VCF_FILE" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# --- Run Sei ---
cd "$SEI_DIR"
bash ./1_variant_effect_prediction.sh "$VCF_FILE" hg19 "$OUTPUT_DIR" --cuda

# --- Cleanup: keep diffs + row_labels, drop ref/alt prediction h5s ---
# Sei's actual filename pattern: <prefix>.ref_predictions.h5 / <prefix>.alt_predictions.h5
if [[ -d "$H5_DIR" ]]; then
    rm -f "$H5_DIR"/*.ref_predictions.h5 "$H5_DIR"/*.alt_predictions.h5
    echo "Removed ref/alt prediction h5s (kept diffs + row_labels)"
fi

echo "Sei predictions completed for $GWAS_NAME"
