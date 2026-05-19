#!/bin/bash
#SBATCH --job-name=sei_mocks
#SBATCH --constraint="l40s"
#SBATCH --cpus-per-task=4
#SBATCH --mem=50000
#SBATCH --gres=gpu:1
#
# NOTE: --partition is intentionally NOT set here — pass your own with sbatch.
# NOTE: --output / --error are NOT set here either; the script redirects logs
#       into ${results_root}/${gwas_name}/slurm_logs/ via exec below.

# Run Sei on ONE mock dataset (indexed by $SLURM_ARRAY_TASK_ID).
#
# Usage (submit as a SLURM array, two equivalent forms):
#
#   # 1. From inside impact_dl/ — $SLURM_SUBMIT_DIR locates paths.py
#   cd /path/to/impact_dl
#   sbatch [--partition=<p>] --array=1-300%20 \
#          prediction_pipeline/_sei_mocks_array.slurm --gwas_name <name>
#
#   # 2. From anywhere with --repo_root explicit
#   sbatch [...] --array=1-300%20 \
#          /path/to/impact_dl/prediction_pipeline/_sei_mocks_array.slurm \
#          --gwas_name <name> --repo_root /path/to/impact_dl
#
# - The array range (1-N) must match `n_mock_datasets` in paths.yaml.
# - The `%20` cap is concurrency — adjust to your cluster's QoS / GPU policy.
# - Failed tasks can be re-submitted with --array=<ids> — completed tasks are
#   skipped via the diffs.h5 existence check.
#
# Other overrides:
#   --paths_yaml /alt/paths.yaml
#   --sei_dir /alt/sei-framework

set -e

# --- Fail fast if not submitted as an array ---
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: this script must be submitted as a SLURM array job." >&2
    echo "       sbatch --array=1-N%C _sei_mocks_array.slurm --gwas_name <name>" >&2
    exit 1
fi
DATASET_ID="$SLURM_ARRAY_TASK_ID"

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
    echo "  - cd into impact_dl/ before sbatch, or" >&2
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
MOCK_NAME="mock_dataset${DATASET_ID}"
VCF_FILE="${BASE_DIR}/vcfs_mock_datasets/${MOCK_NAME}/${MOCK_NAME}_snps_correct_ref_final.vcf"
OUTPUT_DIR="${BASE_DIR}/sei_outputs_mock_datasets/${MOCK_NAME}"
H5_DIR="${OUTPUT_DIR}/chromatin-profiles-hdf5"
DIFFS_FILE="${H5_DIR}/${MOCK_NAME}_snps_correct_ref_final_diffs.h5"

# --- Redirect stdout/stderr into results dir ---
mkdir -p "${BASE_DIR}/slurm_logs"
exec > "${BASE_DIR}/slurm_logs/sei_mock_${DATASET_ID}.out" \
     2> "${BASE_DIR}/slurm_logs/sei_mock_${DATASET_ID}.err"

# Clean up SLURM's default array output file in SLURM_SUBMIT_DIR on exit.
if [[ -n "${SLURM_SUBMIT_DIR:-}" && -n "${SLURM_ARRAY_JOB_ID:-}" ]]; then
    SLURM_DEFAULT_OUT="${SLURM_SUBMIT_DIR}/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
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
echo "Mock ID:       $DATASET_ID"
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
    echo "Diffs already exist for $MOCK_NAME; skipping."
    exit 0
fi

# --- Validate input ---
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
    echo "Removed ref/alt prediction h5s for $MOCK_NAME (kept diffs + row_labels)"
fi

echo "Sei predictions completed for $MOCK_NAME"
