"""
Real-vs-mock comparison: computes empirical p-values per Sei profile.

Inputs (produced by earlier steps):
  - mock_scores_aggregation_{n_maxes}.csv   (from step 3)
  - main_gwas_snps_correct_ref_final_diffs.h5  (from step 1, Sei output)
  - main_gwas_snps_correct_ref_final_row_labels.txt
  - orig_clumps.pkl  (from pre_prediction/sample_null_blocks.py pickle dir)
  - mock_datasets_df.pkl

Outputs:
  - empirical_pvalues{n_maxes}.csv
  - top_{n_maxes}_snps_by_annotation_main.json
  - all_snp_scores_.json
  - empirical_pvalues{n_maxes}.versions.txt  (env info sidecar)

USAGE
  python 4-empirical_pvalues.py --gwas_name <gwas_name>

NOTE: run in the same Python environment as step 3 (m-env). Numerical results
can differ slightly between numpy/scipy versions; the versions file written
alongside the output captures what was used.
"""

import argparse
import json
import os
import sys

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _aggregation_helpers import (
    compute_max_per_clump,
    extract_clumps,
    load_dataframe,
    load_snp_index_map,
    print_env_info,
)


def parse_args():
    parser = argparse.ArgumentParser(description="Compute empirical p-values from real vs mock Sei scores.")
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
    return parser.parse_args()


def collect_all_snp_scores(differences, clumps_dict, snp_index_map):
    """
    Build {annotation_idx: [{'snp', 'clump', 'value'}, ...]} for downstream
    per-profile top-K SNP analysis.
    """
    num_annotations = differences.shape[1]
    all_snp_scores = {ann_idx: [] for ann_idx in range(num_annotations)}

    for clump_name, snps in clumps_dict.items():
        pairs = [(snp_index_map[s], s) for s in snps if s in snp_index_map]
        if not pairs:
            continue
        pairs.sort(key=lambda x: x[0])
        snp_indices = [idx for idx, _ in pairs]
        snps_sorted = [s for _, s in pairs]
        snp_values = np.abs(np.asarray(differences[snp_indices]))
        for row_i, snp in enumerate(snps_sorted):
            row = snp_values[row_i]
            for ann_idx in range(num_annotations):
                all_snp_scores[ann_idx].append({
                    "snp": snp,
                    "clump": clump_name,
                    "value": float(row[ann_idx]),
                })
    return all_snp_scores


def main():
    args = parse_args()
    env_info = print_env_info()

    # Resolve base_dir + n_datasets from paths.yaml
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from paths import load_paths
    paths = load_paths(args.paths_yaml)
    base_dir = args.base_dir or os.path.join(paths['results_root'], args.gwas_name)
    n_datasets = args.n_datasets if args.n_datasets is not None else paths['n_mock_datasets']

    output_dir = f"{base_dir}/sei_score_aggregations_all"
    pickle_dir = f"{base_dir}/pickle_dir"

    mock_values = pd.read_csv(f"{output_dir}/mock_scores_aggregation_{args.n_maxes}.csv")

    orig_clumps = load_dataframe(f"{pickle_dir}/orig_clumps.pkl")
    mock_clumps_df = load_dataframe(f'{pickle_dir}/mock_datasets_df.pkl')

    snp_leaders_big_clumps = list(mock_clumps_df[mock_clumps_df["Mock_Dataset_ID"] == 1]["Real_Clump_ID"])
    orig_clumps = orig_clumps[orig_clumps["SNP"].isin(snp_leaders_big_clumps)]

    with h5py.File(f'{base_dir}/sei_outputs_main_gwas/chromatin-profiles-hdf5/main_gwas_snps_correct_ref_final_diffs.h5', 'r') as f:
        differences = f['data'][:]

    snp_index_map = load_snp_index_map(f'{base_dir}/sei_outputs_main_gwas/chromatin-profiles-hdf5/main_gwas_snps_correct_ref_final_row_labels.txt')
    clumps_dict = extract_clumps(orig_clumps, use_real=True)
    real_values = compute_max_per_clump(differences, clumps_dict, snp_index_map, args.n_maxes)

    all_scores = collect_all_snp_scores(differences, clumps_dict, snp_index_map)
    with open(f"{output_dir}/all_snp_scores_.json", 'w') as f:
        json.dump(all_scores, f)

    top_n_by_annotation = {
        ann_idx: sorted(snp_list, key=lambda x: x['value'], reverse=True)[:args.n_maxes]
        for ann_idx, snp_list in all_scores.items()
    }
    with open(f"{output_dir}/top_{args.n_maxes}_snps_by_annotation_main.json", 'w') as f:
        json.dump(top_n_by_annotation, f, indent=2)

    if not np.issubdtype(mock_values.iloc[:, 0].dtype, np.number):
        print(f"Dropping extra column: {mock_values.columns[0]}")
        mock_values = mock_values.iloc[:, 1:]
    print(mock_values.columns)
    mock_values = mock_values.to_numpy()

    assert mock_values.shape == (n_datasets, 21907), f"Mock dataset shape mismatch! Found {mock_values.shape}"
    assert real_values.shape == (21907,), f"Real dataset shape mismatch! Found {real_values.shape}"

    mock_mean = np.mean(mock_values, axis=0)
    print(np.sum(mock_values >= real_values, axis=0))
    p_upper = np.sum(mock_values >= real_values, axis=0) / n_datasets
    p_lower = np.sum(mock_values <= real_values, axis=0) / n_datasets
    p_two_sided = 2 * np.minimum(p_upper, p_lower)
    effect_size = real_values - mock_mean

    p_values_df = pd.DataFrame({
        "Empirical_P_Upper": p_upper,
        "Empirical_P_Lower": p_lower,
        "Empirical_P_Two_Sided": p_two_sided,
        "Effect_Size": effect_size,
    })
    p_values_df.index.name = "Annotation_Index"

    print(p_values_df.head())
    out_csv = f"{output_dir}/empirical_pvalues{args.n_maxes}.csv"
    p_values_df.to_csv(out_csv)

    # Sidecar: record numpy/scipy/env used so a future reader can reproduce.
    versions_path = f"{output_dir}/empirical_pvalues{args.n_maxes}.versions.txt"
    with open(versions_path, "w") as f:
        for k, v in env_info.items():
            f.write(f"{k}: {v}\n")
    print(f"Saved env info to {versions_path}")


if __name__ == "__main__":
    main()
