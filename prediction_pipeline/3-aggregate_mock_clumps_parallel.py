import pickle 
import pandas as pd 
import h5py
import numpy as np
from tqdm import tqdm
import argparse
import os
from multiprocessing import Pool
import ast  

# ---------------------- Argument Parsing ---------------------- #

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Aggregate mock datasets Sei scores")
    parser.add_argument("base_dir", metavar="path", type=str, help="Path to the directory containing data files")
    parser.add_argument("--n_maxes", type=int, default=50, help="The number of maxes Sei scores to average per annotation, where each score is a clump's max  (default: 50)")
    parser.add_argument("--n_datasets", type=int, default=300, help="The number of mock datasets  (default: 300)")
    parser.add_argument("--n_workers", type=int, default=None,
            help="Number of parallel workers (default: all available CPUs)")

    return parser.parse_args()


def load_dataframe(filepath):
    """Load a DataFrame from a Pickle file."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"No file found at {filepath}")
    with open(filepath, "rb") as f:
        df = pickle.load(f)
    print(f"DataFrame loaded from {filepath}")
    return df

def load_hdf5(filepath, key="data", query=None):
    with pd.HDFStore(filepath, mode="r") as store:
        if query:
            return store.select(key, where=query)  # Fast querying
        else:
            return store[key]  # Load entire DataFrame

def parse_list_column(x):
    """Convert string representation of list back to actual list."""
    if pd.isna(x) or x == '' or x == 'nan':
        return []
    try:
        return ast.literal_eval(x)
    except (ValueError, SyntaxError):
        return []


def extract_clumps(clumps_df, use_real=False):
    """
    Extracts clumps from a dataset.
    Clumps are defined as Sampled_SNP + Clumped_SNPs.
    
    Parameters:
    - clumps_df (pd.DataFrame): DataFrame containing SNP and clumped SNP data.
    - use_real (bool): If True, uses 'SNP' as the clump key (for real data).
                       If False, uses 'Sampled_SNP' (for mock data).
    
    Returns:
    - dict: {clump_name: [snp1, snp2, ...]}
    """
    clump_mapping = {}
    
    key_column = "SNP" if use_real else "Sampled_SNP"  # Choose the key- Sampled is the leader of the mock clump, SNP is the leader for the real

    for _, row in clumps_df.iterrows():
        clump_name = row[key_column]  # Use the selected key
        clumped_snps = row['Clumped_SNPs']  # Get the list of clumped SNPs

        # Store in dictionary
        clump_mapping[clump_name] = [clump_name] + clumped_snps  # Include the key SNP

    return clump_mapping

def compute_max_per_clump(differences, clumps_dict, snp_index_map, n=10, stat="max"):
    """
    Computes the specified statistic (max or second max) per clump and averages the top `n` values per annotation.
    
    Parameters:
    - `differences`: NumPy array (rows = SNPs, cols = annotations)
    - `clumps_dict`: Dictionary {clump_name: [snp1, snp2, ...]}
    - `snp_index_map`: Dictionary {snp_name: row_index}
    - `n`: Number of top values to average (default: 10)
    - `stat`: Statistic to compute ("max" or "second_max")
    
    Returns:
    - 1D NumPy array: Averaged top `n` values per annotation.
    """
    num_annotations = differences.shape[1]
    num_clumps = len(clumps_dict)

    clump_stats = np.zeros((num_clumps, num_annotations))  # Store max per clump

    # Process each clump
    for i, (clump_name, snps) in enumerate(clumps_dict.items()):
        # Get row indices for all SNPs in this clump
        snp_indices = [snp_index_map[snp] for snp in snps if snp in snp_index_map]
        snp_indices = sorted(snp_indices)

        out_of_range_indices = [idx for idx in snp_indices if idx >= differences.shape[0] or idx < 0]
        if out_of_range_indices:
            print(f"  Warning: Clump '{clump_name}' has out-of-range indices: {out_of_range_indices}")
            print(f"  SNPs causing issue: {[snps[j] for j, idx in enumerate(snp_indices) if idx in out_of_range_indices]}")
            print(f"  Max allowed index: {differences.shape[0] - 1}")


        if len(snp_indices) > 0:
            clump_values = np.abs(differences[snp_indices])  # Take absolute values
            if stat == "max":
                clump_stats[i] = np.max(clump_values, axis=0)  # Max per annotation
            elif stat == "second_max" and len(snp_indices) > 1:
                clump_stats[i] = np.sort(clump_values, axis=0)[-2]  # Second largest per annotation


    # Compute top `n` values and average
    top_n_stats = np.sort(clump_stats, axis=0)[-n:]  # Last `n` values = Top `n` stats
    avg_top_n = np.mean(top_n_stats, axis=0)  # Average top `n` values per annotation

    return avg_top_n

def compute_mean_of_all_maxes_per_clump(differences, clumps_dict, snp_index_map):
    """
    Computes the max values per clump, all values per annotation, and averages them.
    - `differences`: NumPy array (rows = SNPs, cols = annotations)
    - `clumps_dict`: Dictionary {clump_name: [snp1, snp2, ...]}
    - `snp_index_map`: Dictionary {snp_name: row_index}
    """
    num_annotations = differences.shape[1]
    num_clumps = len(clumps_dict)

    clump_maxes = np.zeros((num_clumps, num_annotations))  # Store max per clump

    # Process each clump
    for i, (clump_name, snps) in enumerate(clumps_dict.items()):
        # Get row indices for all SNPs in this clump
        snp_indices = [snp_index_map[snp] for snp in snps if snp in snp_index_map]
        snp_indices = sorted(snp_indices)

        out_of_range_indices = [idx for idx in snp_indices if idx >= differences.shape[0] or idx < 0]
        if out_of_range_indices:
            print(f"  Warning: Clump '{clump_name}' has out-of-range indices: {out_of_range_indices}")
            print(f"  SNPs causing issue: {[snps[j] for j, idx in enumerate(snp_indices) if idx in out_of_range_indices]}")
            print(f"  Max allowed index: {differences.shape[0] - 1}")


        if len(snp_indices) > 0:
            clump_maxes[i] = np.max(np.abs(differences[snp_indices]), axis=0) # Max per annotation, taking absolute values for max CHANGE

    avg_of_all_maxes = np.mean(clump_maxes, axis=0)  # Mean of all maxes across clumps

    return avg_of_all_maxes

def load_snp_index_map(mock_row_labels_file):
    labels_df = pd.read_csv(mock_row_labels_file, sep='\t')
    labels_df['name'] = labels_df['name'].astype(str)
    return {snp: i for i, snp in enumerate(labels_df['name'])}


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
    base_dir = args.base_dir
    output_dir = f"{base_dir}/sei_score_aggregations_all"
    pickle_dir = f"{base_dir}/pickle_dir"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # File paths
    mock_differences_files = [f"{base_dir}/sei_outputs_mock_datasets/mock_dataset{i}/chromatin-profiles-hdf5/mock_dataset{i}_snps_correct_ref_final_diffs.h5" for i in range(1, args.n_datasets +1)]
    print(mock_differences_files[:5])
    mock_row_labels_files = [f"{base_dir}/sei_outputs_mock_datasets/mock_dataset{i}/chromatin-profiles-hdf5/mock_dataset{i}_snps_correct_ref_final_row_labels.txt" for i in range(1, args.n_datasets +1)]
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

    n_workers = args.n_workers or os.cpu_count()
    print(f"Processing {len(tasks)} mock datasets with {n_workers} workers...")

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
    results = [results_dict[i] for i in range(1, args.n_datasets + 1)]

    # Convert to DataFrame
    final_df = pd.DataFrame(results, index=[f'mock_{i}' for i in range(1, args.n_datasets +1 )])
    print(final_df.head())

    # Save results
    final_df.to_csv(f"{output_dir}/mock_scores_aggregation_{args.n_maxes}.csv")
    print("Saved full mock max matrix.")

    # Save SNP-to-clump mapping for each mock (optional, can be large)
    mock_snp_clump_map = []
    for dataset_id in range(1, args.n_datasets + 1):
        clumps_dict = extract_clumps(mock_clumps_df[mock_clumps_df['Mock_Dataset_ID'] == dataset_id])
        for clump, snps in clumps_dict.items():
            for snp in snps:
                mock_snp_clump_map.append({'Mock_Dataset': dataset_id, 'Clump': clump, 'SNP': snp})
    mock_snp_clump_df = pd.DataFrame(mock_snp_clump_map)
    mock_snp_clump_df.to_csv(f"{output_dir}/mock_snp_to_clump_map.csv", index=False)
    print("Saved mock SNP-to-clump mapping.")

if __name__ == "__main__":
    main()