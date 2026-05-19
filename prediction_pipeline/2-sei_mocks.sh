#!/bin/bash
# Batched submission of _sei_mocks_array.slurm to avoid QOS submit limits.
#
# Submits the N mock datasets in chunks of BATCH_SIZE, waiting for each chunk
# to finish before submitting the next. Each per-task script is resumable
# (h5-existence check), so killing and re-running this helper is safe.
#
# RUN INSIDE tmux OR screen — this stays in the foreground until all mocks
# are done (typically several hours).
#
# Usage:
#   bash 2-sei_mocks.sh --gwas_name <name> [options]
#
# Options (defaults shown):
#   --gwas_name <name>          REQUIRED
#   --partition <p>             (no default; passed to sbatch if set)
#   --batch_size 80             how many array indices to submit at once
#   --array_concurrency 20      the "%N" cap on concurrent running tasks
#   --total <N>                 default: n_mock_datasets from paths.yaml
#   --paths_yaml <path>         override paths.yaml location
#   --repo_root <path>          override repo root (default: $SLURM_SUBMIT_DIR
#                               or this script's parent)

set -e

# --- Parse args ---
GWAS_NAME=""
PARTITION=""
BATCH_SIZE=80
ARRAY_CONCURRENCY=20
TOTAL=""
PATHS_YAML=""
REPO_ROOT_ARG=""
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --gwas_name)         GWAS_NAME="$2"; shift 2 ;;
        --partition)         PARTITION="$2"; shift 2 ;;
        --batch_size)        BATCH_SIZE="$2"; shift 2 ;;
        --array_concurrency) ARRAY_CONCURRENCY="$2"; shift 2 ;;
        --total)             TOTAL="$2"; shift 2 ;;
        --paths_yaml)        PATHS_YAML="$2"; shift 2 ;;
        --repo_root)         REPO_ROOT_ARG="$2"; shift 2 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$GWAS_NAME" ]]; then
    echo "ERROR: --gwas_name is required" >&2
    exit 1
fi

# --- Locate repo root (contains paths.py) ---
REPO_ROOT=""
if [[ -n "$REPO_ROOT_ARG" ]]; then
    REPO_ROOT="$REPO_ROOT_ARG"
elif [[ -f "$(pwd)/paths.py" ]]; then
    REPO_ROOT="$(pwd)"
elif [[ -f "$(pwd)/../paths.py" ]]; then
    REPO_ROOT="$(cd .. && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)" || true
    if [[ -n "$SCRIPT_DIR" && -f "$SCRIPT_DIR/../paths.py" ]]; then
        REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
    fi
fi

if [[ -z "$REPO_ROOT" || ! -f "$REPO_ROOT/paths.py" ]]; then
    echo "ERROR: cannot locate impact_dl/paths.py. Either:" >&2
    echo "  - cd into impact_dl/ before running this script, or" >&2
    echo "  - pass --repo_root /path/to/impact_dl." >&2
    exit 1
fi

# --- Resolve total mocks ---
paths_args=""
if [[ -n "$PATHS_YAML" ]]; then paths_args="--paths_yaml $PATHS_YAML"; fi
if [[ -z "$TOTAL" ]]; then
    TOTAL=$(python "$REPO_ROOT/paths.py" $paths_args --get n_mock_datasets)
fi

# --- Resolve base_dir (for completion counting) ---
RESULTS_ROOT=$(python "$REPO_ROOT/paths.py" $paths_args --get results_root)
BASE_DIR="${RESULTS_ROOT}/${GWAS_NAME}"
H5_PARENT="${BASE_DIR}/sei_outputs_mock_datasets"

# --- Build sbatch passthrough args for the array script ---
PARTITION_ARG=""
if [[ -n "$PARTITION" ]]; then PARTITION_ARG="--partition=$PARTITION"; fi

SCRIPT_ARGS=(--gwas_name "$GWAS_NAME" --repo_root "$REPO_ROOT")
if [[ -n "$PATHS_YAML" ]]; then SCRIPT_ARGS+=(--paths_yaml "$PATHS_YAML"); fi

# Count completed mock diffs.h5 files in [start..end] inclusive.
count_done() {
    local start="$1" end="$2" done=0
    for i in $(seq "$start" "$end"); do
        if [[ -f "${H5_PARENT}/mock_dataset${i}/chromatin-profiles-hdf5/mock_dataset${i}_snps_correct_ref_final_diffs.h5" ]]; then
            done=$((done + 1))
        fi
    done
    echo "$done"
}

# --- Submit in batches, polling for progress between them ---
NUM_BATCHES=$(( (TOTAL + BATCH_SIZE - 1) / BATCH_SIZE ))
echo "================================================================"
echo "GWAS:               $GWAS_NAME"
echo "Total mocks:        $TOTAL"
echo "Batch size:         $BATCH_SIZE"
echo "Array concurrency:  %$ARRAY_CONCURRENCY"
echo "Number of batches:  $NUM_BATCHES"
echo "Partition:          ${PARTITION:-<cluster default>}"
echo "Repo root:          $REPO_ROOT"
echo "Base dir:           $BASE_DIR"
echo "================================================================"
echo "NOTE: this will take several hours. Run inside tmux/screen."
echo "      Re-runs skip already-completed mocks via h5 existence check."
echo

POLL_INTERVAL=30  # seconds between progress polls within a batch

for b in $(seq 1 $NUM_BATCHES); do
    START=$(( (b - 1) * BATCH_SIZE + 1 ))
    END=$(( b * BATCH_SIZE ))
    if [[ $END -gt $TOTAL ]]; then END=$TOTAL; fi
    BATCH_N=$((END - START + 1))

    pre_done=$(count_done "$START" "$END")
    if [[ $pre_done -eq $BATCH_N ]]; then
        echo "[$(date +%T)] Batch $b/$NUM_BATCHES (mocks ${START}-${END}): all $BATCH_N already done — skipping."
        echo
        continue
    fi

    echo "[$(date +%T)] Batch $b/$NUM_BATCHES (mocks ${START}-${END}, $pre_done already done): submitting --array=${START}-${END}%${ARRAY_CONCURRENCY}"
    JOB_ID=$(sbatch --parsable $PARTITION_ARG \
                    --array="${START}-${END}%${ARRAY_CONCURRENCY}" \
                    "$REPO_ROOT/prediction_pipeline/_sei_mocks_array.slurm" \
                    "${SCRIPT_ARGS[@]}")
    echo "[$(date +%T)]   submitted as job $JOB_ID"

    # Poll until no tasks pending/running for this array.
    last_done=-1
    while true; do
        sleep "$POLL_INTERVAL"
        pending=$(squeue -j "$JOB_ID" -h -t PENDING 2>/dev/null | wc -l)
        running=$(squeue -j "$JOB_ID" -h -t RUNNING 2>/dev/null | wc -l)
        done=$(count_done "$START" "$END")
        cum_done=$(count_done 1 "$TOTAL")
        if [[ $done -ne $last_done ]]; then
            echo "[$(date +%T)]   batch $b: $done/$BATCH_N done, $running running, $pending queued | overall $cum_done/$TOTAL"
            last_done=$done
        fi
        if [[ $pending -eq 0 && $running -eq 0 ]]; then
            break
        fi
    done

    final_done=$(count_done "$START" "$END")
    echo "[$(date +%T)] Batch $b/$NUM_BATCHES done: $final_done/$BATCH_N mocks completed."
    echo
done

cum_done=$(count_done 1 "$TOTAL")
echo "================================================================"
echo "Done. $cum_done/$TOTAL mock diffs.h5 files present."
if [[ $cum_done -lt $TOTAL ]]; then
    echo "WARNING: some mocks failed to produce diffs.h5. Re-run this script"
    echo "         to retry (completed mocks are skipped automatically)."
fi
echo "Next: run step 3c (aggregate mocks) and 3d (empirical p-values) — locally."
echo "================================================================"
