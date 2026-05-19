"""
Shared helpers for Sei score aggregation.

Imported by both 3-aggregate_mocks.py (mock-side aggregation) and
4-empirical_pvalues.py (real-side aggregation + empirical p-values).
Logic identical to the inline versions that previously lived in
3-aggregate_mock_clumps_parallel.py.
"""

import ast
import os
import pickle

import numpy as np
import pandas as pd


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
            return store.select(key, where=query)
        return store[key]


def parse_list_column(x):
    """Convert string representation of list back to actual list."""
    if pd.isna(x) or x == '' or x == 'nan':
        return []
    try:
        return ast.literal_eval(x)
    except (ValueError, SyntaxError):
        return []


def load_snp_index_map(row_labels_file):
    """Read a Sei *_row_labels.txt and return {snp_id: row_index}."""
    labels_df = pd.read_csv(row_labels_file, sep='\t')
    labels_df['name'] = labels_df['name'].astype(str)
    return {snp: i for i, snp in enumerate(labels_df['name'])}


def extract_clumps(clumps_df, use_real=False):
    """
    Extracts clumps from a dataset. A clump is the index SNP plus its tag SNPs.

    use_real=True : 'SNP' is the leader column (real GWAS clumps).
    use_real=False: 'Sampled_SNP' is the leader column (mock clumps).
    """
    key_column = "SNP" if use_real else "Sampled_SNP"
    clump_mapping = {}
    for _, row in clumps_df.iterrows():
        leader = row[key_column]
        clump_mapping[leader] = [leader] + row['Clumped_SNPs']
    return clump_mapping


def compute_max_per_clump(differences, clumps_dict, snp_index_map, n=10, stat="max"):
    """
    For each clump, take max |effect| per Sei annotation. Then take the top-n
    clump scores per annotation and average them. Returns a 1D numpy array
    of length N_annotations.
    """
    num_annotations = differences.shape[1]
    num_clumps = len(clumps_dict)
    clump_stats = np.zeros((num_clumps, num_annotations))

    for i, (clump_name, snps) in enumerate(clumps_dict.items()):
        snp_indices = [snp_index_map[s] for s in snps if s in snp_index_map]
        snp_indices = sorted(snp_indices)

        out_of_range = [idx for idx in snp_indices if idx >= differences.shape[0] or idx < 0]
        if out_of_range:
            print(f"  Warning: Clump '{clump_name}' has out-of-range indices: {out_of_range}")
            print(f"  SNPs causing issue: {[snps[j] for j, idx in enumerate(snp_indices) if idx in out_of_range]}")
            print(f"  Max allowed index: {differences.shape[0] - 1}")

        if snp_indices:
            clump_values = np.abs(differences[snp_indices])
            if stat == "max":
                clump_stats[i] = np.max(clump_values, axis=0)
            elif stat == "second_max" and len(snp_indices) > 1:
                clump_stats[i] = np.sort(clump_values, axis=0)[-2]

    top_n_stats = np.sort(clump_stats, axis=0)[-n:]
    return np.mean(top_n_stats, axis=0)


def compute_mean_of_all_maxes_per_clump(differences, clumps_dict, snp_index_map):
    """
    Alternative aggregation: for each clump, max |effect| per annotation;
    then mean across ALL clumps (not top-n).
    """
    num_annotations = differences.shape[1]
    num_clumps = len(clumps_dict)
    clump_maxes = np.zeros((num_clumps, num_annotations))

    for i, (clump_name, snps) in enumerate(clumps_dict.items()):
        snp_indices = [snp_index_map[s] for s in snps if s in snp_index_map]
        snp_indices = sorted(snp_indices)

        out_of_range = [idx for idx in snp_indices if idx >= differences.shape[0] or idx < 0]
        if out_of_range:
            print(f"  Warning: Clump '{clump_name}' has out-of-range indices: {out_of_range}")
            print(f"  SNPs causing issue: {[snps[j] for j, idx in enumerate(snp_indices) if idx in out_of_range]}")
            print(f"  Max allowed index: {differences.shape[0] - 1}")

        if snp_indices:
            clump_maxes[i] = np.max(np.abs(differences[snp_indices]), axis=0)

    return np.mean(clump_maxes, axis=0)


def print_env_info():
    """Print numpy/scipy versions + active conda env. Numerical results can
    differ slightly between environments, so we record what we used."""
    import sys
    try:
        import scipy
        scipy_ver = scipy.__version__
    except ImportError:
        scipy_ver = "(not installed)"
    info = {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy_ver,
        "conda_env": os.environ.get("CONDA_DEFAULT_ENV", "<none>"),
    }
    print("Environment:")
    for k, v in info.items():
        print(f"  {k}: {v}")
    return info
