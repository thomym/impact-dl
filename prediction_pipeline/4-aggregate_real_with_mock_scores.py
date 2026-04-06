import pandas as pd
import pickle
import numpy as np
import h5py
import argparse
import matplotlib.pyplot as plt
import json

from aggregate_mock_clumps_parallel import load_dataframe, load_snp_index_map, extract_clumps, compute_max_per_clump

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Aggregate mock datasets Sei scores")
    parser.add_argument("base_dir", metavar="path", type=str, help="Path to the directory containing data files")
    parser.add_argument("--n_maxes", type=int, default=50, help="The number of maxes Sei scores to average per annotation, where each score is a clump's max  (default: 50)")
    parser.add_argument("--n_datasets", type=int, default=300, help="The number of mock datasets  (default: 300)")

    return parser.parse_args()



def collect_all_snp_scores(differences, clumps_dict, snp_index_map):
    """
    Collects scores for all SNPs across all annotations, organized by annotation.
    
    Parameters:
    - `differences`: NumPy array (rows = SNPs, cols = annotations)
    - `clumps_dict`: Dictionary {clump_name: [snp1, snp2, ...]}
    - `snp_index_map`: Dictionary {snp_name: row_index}
    
    Returns:
    - Dictionary: {annotation_idx: [{'snp': snp_name, 'clump': clump_name, 'value': score}, ...]}
      Contains all SNPs and their scores for each annotation.
    """
    num_annotations = differences.shape[1]
    
    # Initialize dictionary for all annotations
    all_snp_scores = {ann_idx: [] for ann_idx in range(num_annotations)}
    
    # Process each clump
    for clump_name, snps in clumps_dict.items():
        # Build (row_index, snp_id) pairs, then sort by row_index
        pairs = [(snp_index_map[s], s) for s in snps if s in snp_index_map]
        if not pairs:
            continue

        pairs.sort(key=lambda x: x[0])
        snp_indices = [idx for idx, _ in pairs]
        snps_sorted = [s for _, s in pairs]

        # Pull the matching rows in the same order as snps_sorted
        snp_values = np.abs(np.asarray(differences[snp_indices]))  # shape: [k, num_annotations]

        # Record scores
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
    base_dir = args.base_dir
    output_dir = f"{base_dir}/sei_score_aggregations_all"
    pickle_dir = f"{base_dir}/pickle_dir"
    mock_values = pd.read_csv(f"{output_dir}/mock_scores_aggregation_{args.n_maxes}.csv")

    orig_clumps = load_dataframe(f"{pickle_dir}/orig_clumps.pkl")
    mock_clumps_df = load_dataframe(f'{pickle_dir}/mock_datasets_df.pkl')

    snp_leaders_big_clumps = list(mock_clumps_df[mock_clumps_df["Mock_Dataset_ID"] == 1]["Real_Clump_ID"])
    orig_clumps = orig_clumps[orig_clumps["SNP"].isin(snp_leaders_big_clumps)]
    f = h5py.File(f'{base_dir}/sei_outputs_main_gwas/chromatin-profiles-hdf5/main_gwas_snps_correct_ref_final_diffs.h5', 'r')
    differences = f['data']
    snp_index_map = load_snp_index_map(f'{base_dir}/sei_outputs_main_gwas/chromatin-profiles-hdf5/main_gwas_snps_correct_ref_final_row_labels.txt')  # Map SNP names to row indices
    clumps_dict = extract_clumps(orig_clumps, use_real=True)
    maxes = compute_max_per_clump(differences, clumps_dict, snp_index_map, args.n_maxes)#, stat="second_max")
    all_scores = collect_all_snp_scores(differences, clumps_dict, snp_index_map)
    with open(f"{output_dir}/all_snp_scores_.json", 'w') as f:
        json.dump(all_scores, f)
    top_n_by_annotation = {}
    for ann_idx, snp_list in all_scores.items():
        # Sort by value in descending order
        sorted_snps = sorted(snp_list, key=lambda x: x['value'], reverse=True)
        # Keep only the top 50
        top_n_by_annotation[ann_idx] = sorted_snps[:args.n_maxes]
    
    # Save top 50
    with open(f"{output_dir}/top_{args.n_maxes}_snps_by_annotation_main.json", 'w') as f:
        json.dump(top_n_by_annotation, f, indent=2)
        

    real_values = maxes
    if not np.issubdtype(mock_values.iloc[:, 0].dtype, np.number):
        print(f"Dropping extra column: {mock_values.columns[0]}")
        mock_values = mock_values.iloc[:, 1:]  # Drop the first column - Unnamed:0 index column
    print(mock_values.columns)
    # Convert to NumPy for fast calculations
    mock_values = mock_values.to_numpy()
    # Ensure shapes match: (args.n_datasets, 21907) for mocks and (21907,) for real
    assert mock_values.shape == (args.n_datasets, 21907), f"Mock dataset shape mismatch! Found {mock_values.shape}"
    assert real_values.shape == (21907,), f"Real dataset shape mismatch! Found {real_values.shape}"
    mock_mean = np.mean(mock_values, axis=0)  # Mean of mock distributions
    print(np.sum(mock_values >= real_values, axis=0))
    p_upper = np.sum(mock_values >= real_values, axis=0) / args.n_datasets  # Real is higher
    p_lower = np.sum(mock_values <= real_values, axis=0) / args.n_datasets  # Real is lower

    # Compute two-sided p-values
    p_two_sided = 2 * np.minimum(p_upper, p_lower)

    # Compute effect size
    effect_size = real_values - mock_mean  # Difference between real and null mean

    # Create DataFrame
    p_values_df = pd.DataFrame({
        "Empirical_P_Upper": p_upper,
        "Empirical_P_Lower": p_lower,
        "Empirical_P_Two_Sided": p_two_sided,
        "Effect_Size": effect_size  # Added Effect Size
    })

    p_values_df.index.name = "Annotation_Index"  # Annotation index is the row number

    # Display first few rows
    print(p_values_df.head())

    # Save to CSV
    p_values_df.to_csv(f"{output_dir}/empirical_pvalues{args.n_maxes}.csv")

if __name__ == "__main__":
    main()
