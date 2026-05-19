import argparse
import os
import sys
from multiprocessing import Pool

import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _aggregation_helpers import (
    compute_max_per_clump,
    extract_clumps,
    load_snp_index_map,
    parse_list_column,
    print_env_info,
)

# ---------------------- Argument Parsing ---------------------- #

DEFAULT_MAX_WORKERS = 8  # Each worker loads one mock h5 (~600 MB); capping avoids OOM on big nodes.


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Aggregate mock-dataset Sei scores into a per-profile vector per mock.")
    parser.add_argument("--gwas_name", required=True,
                        help="GWAS name (subdir under results_root).")
    parser.add_argument("--base_dir", default=None,
                        help="Override base directory (default: {results_root}/{gwas_name}).")
    parser.add_argument("--paths_yaml", default=None,
                        help="Path to paths.yaml (default: <repo>/paths.yaml).")
    parser.add_argument("--n_maxes", type=int, default=50,
                        help="Top-K clump scores to average per annotation (default: 50).")
    parser.add_argument("--n_datasets", type=int, default=None,
                        help="Number of mock datasets (default: n_mock_datasets from paths.yaml).")
    parser.add_argument("--n_workers", type=int, default=None,
                        help=f"Parallel workers (default: min(cpu_count, {DEFAULT_MAX_WORKERS})). "
                             f"Each worker uses ~1 GB; raise consciously on large machines.")

    return parser.parse_args()


# ---------------------- Multiprocessing setup ---------------------- #

# Shared state: each worker gets a copy of this via the initializer, so
# mock_clumps_df is pickled once per worker process, not once per task.
_shared_clumps_df = None

def _init_worker(clumps_df):
    """Pool initializer — stashes the shared DataFrame in each worker."""
    global _shared_clumps_df
    _shared_clumps_df = clumps_df

def _process_single_dataset(task):
    """
    Worker function: runs the full pipeline for one mock dataset.

    task : tuple of (dataset_index, diff_file, label_file, n_maxes)
    Returns: (dataset_index, avg_top_n)   — index kept so we can reorder later.
    """
    global _shared_clumps_df
    i, diff_file, label_file, n_maxes = task

    snp_index_map = load_snp_index_map(label_file)

    with h5py.File(diff_file, 'r') as h5:
        differences = h5['data'][:]

    dataset_df = _shared_clumps_df[_shared_clumps_df['Mock_Dataset_ID'] == i]
    clumps_dict = extract_clumps(dataset_df)

    avg_top_n = compute_max_per_clump(differences, clumps_dict, snp_index_map, n=n_maxes)
    return (i, avg_top_n)


def main():
    args = parse_args()
    print_env_info()

    # Resolve base_dir + n_datasets from paths.yaml (with optional CLI overrides).
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from paths import load_paths
    paths = load_paths(args.paths_yaml)
    base_dir = args.base_dir or os.path.join(paths['results_root'], args.gwas_name)
    n_datasets = args.n_datasets if args.n_datasets is not None else paths['n_mock_datasets']

    output_dir = f"{base_dir}/sei_score_aggregations_all"
    pickle_dir = f"{base_dir}/pickle_dir"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # File paths
    mock_differences_files = [f"{base_dir}/sei_outputs_mock_datasets/mock_dataset{i}/chromatin-profiles-hdf5/mock_dataset{i}_snps_correct_ref_final_diffs.h5" for i in range(1, n_datasets + 1)]
    print(mock_differences_files[:5])
    mock_row_labels_files = [f"{base_dir}/sei_outputs_mock_datasets/mock_dataset{i}/chromatin-profiles-hdf5/mock_dataset{i}_snps_correct_ref_final_row_labels.txt" for i in range(1, n_datasets + 1)]
    print(mock_row_labels_files[:5])

    mock_clumps_df = pd.read_csv(
        f"{pickle_dir}/mock_datasets_df.csv",
        converters={'Clumped_SNPs': parse_list_column}
    )

    # --- Build the task list (one tuple per mock dataset) ---
    tasks = [
        (i, diff_file, label_file, args.n_maxes)
        for i, (diff_file, label_file) in enumerate(
            zip(mock_differences_files, mock_row_labels_files), start=1
        )
    ]

    cpu = os.cpu_count() or 1
    n_workers = args.n_workers or min(cpu, DEFAULT_MAX_WORKERS)
    print(f"Processing {len(tasks)} mock datasets with {n_workers} workers "
          f"(cpu_count={cpu}, default cap={DEFAULT_MAX_WORKERS}). "
          f"Override with --n_workers <N> if you have memory headroom.")

    # --- Parallel execution ---
    # imap_unordered returns results as they finish (faster than imap),
    # so we collect into a dict keyed by dataset index and sort afterwards.
    results_dict = {}
    with Pool(processes=n_workers,
              initializer=_init_worker,
              initargs=(mock_clumps_df,)) as pool:
        for idx, avg_top_n in tqdm(
            pool.imap_unordered(_process_single_dataset, tasks),
            total=len(tasks),
            desc="Mock datasets"
        ):
            results_dict[idx] = avg_top_n

    # Reassemble in dataset order
    results = [results_dict[i] for i in range(1, n_datasets + 1)]

    # Convert to DataFrame
    final_df = pd.DataFrame(results, index=[f'mock_{i}' for i in range(1, n_datasets + 1)])
    print(final_df.head())

    # Save results
    final_df.to_csv(f"{output_dir}/mock_scores_aggregation_{args.n_maxes}.csv")
    print("Saved full mock max matrix.")

    # Save SNP-to-clump mapping for each mock (optional, can be large)
    mock_snp_clump_map = []
    for dataset_id in range(1, n_datasets + 1):
        clumps_dict = extract_clumps(mock_clumps_df[mock_clumps_df['Mock_Dataset_ID'] == dataset_id])
        for clump, snps in clumps_dict.items():
            for snp in snps:
                mock_snp_clump_map.append({'Mock_Dataset': dataset_id, 'Clump': clump, 'SNP': snp})
    mock_snp_clump_df = pd.DataFrame(mock_snp_clump_map)
    mock_snp_clump_df.to_csv(f"{output_dir}/mock_snp_to_clump_map.csv", index=False)
    print("Saved mock SNP-to-clump mapping.")

if __name__ == "__main__":
    main()