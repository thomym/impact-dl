import pandas as pd
import os
import pickle
import argparse


# ---------------------- Argument Parsing ---------------------- #
def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Process SNP clumping data")
    parser.add_argument("base_dir", metavar="path", type=str, help="Path to the directory containing data files")
    parser.add_argument("--n_mock_datasets", type=int, default=300, help="Number of mock datasets to generate (default: 300)")
    parser.add_argument("--use_pickles", action="store_true", help="Flag to use existing pickle files instead of generating them.")

    return parser.parse_args()


# ---------------------- File Handling ---------------------- #
def save_dataframe(df, filepath):
    """Save a DataFrame as a Pickle file."""
    with open(filepath, "wb") as f:
        pickle.dump(df, f)
    print(f"DataFrame saved to {filepath}")


def load_dataframe(filepath):
    """Load a DataFrame from a Pickle file."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"No file found at {filepath}")
    with open(filepath, "rb") as f:
        df = pickle.load(f)
    print(f"DataFrame loaded from {filepath}")
    return df


# ---------------------- Data Loading Functions ---------------------- #
def load_ld_data(pre_clumping_dir, pickle_dir, use_pickles):
    """Load LD r2-filtered SNPs, or load from pickle if specified."""
    pickle_path = f"{pickle_dir}/r2_filtered_snps.pkl"
    if use_pickles and os.path.exists(pickle_path):
        return load_dataframe(pickle_path)
    else:
        print("Loading raw LD data...")
        df = pd.read_csv(f"{pre_clumping_dir}/r2_filtered_snps.ld", sep=r"\s+")
        save_dataframe(df, pickle_path)
        return df


def load_clumps(clumping_dir, pickle_dir, use_pickles):
    """Load cleaned clump data, or load from pickle if specified."""
    pickle_path = f"{pickle_dir}/clumps_cleaned_noncoding.pkl"
    if use_pickles and os.path.exists(pickle_path):
        return load_dataframe(pickle_path)
    else:
        print("Loading raw clump data...")
        try:
            df = pd.read_csv(f"{clumping_dir}/clumps_cleaned_noncoding.clumped", sep=r"\s+")
            if df.empty:
                print(f"WARNING: No clumps found for this GWAS. Skipping.")
                return None
            save_dataframe(df, pickle_path)  
            return df
        except pd.errors.EmptyDataError:
            print(f"WARNING: Empty clumps file for this GWAS. Skipping.")
            return None


def load_maf(pre_clumping_dir, pickle_dir, use_pickles):
    """Load minor allele frequency (MAF) data, or load from pickle if specified."""
    pickle_path = f"{pickle_dir}/maf_lookup_with_header.pkl"
    if use_pickles and os.path.exists(pickle_path):
        return load_dataframe(pickle_path)
    else:
        print("Loading raw MAF data...")
        df = pd.read_csv(f"{pre_clumping_dir}/maf_lookup_with_header.tsv", sep=r"\s+")
        save_dataframe(df, pickle_path)
        return df


def load_gwas(pre_clumping_dir, pickle_dir, use_pickles):
    """Load GWAS data filtered by 1000 Genomes, or load from pickle if specified."""
    pickle_path = f"{pickle_dir}/gwas_filtered_by_1kg.pkl"
    if use_pickles and os.path.exists(pickle_path):
        return load_dataframe(pickle_path)
    else:
        print("Loading raw GWAS data...")
        df = pd.read_csv(f"{pre_clumping_dir}/gwas_filtered_by_1kg.txt", sep=r"\s+")
        save_dataframe(df, pickle_path)
        return df

# ---------------------- Data Processing ---------------------- #

def merge_maf_gwas(full_gwas, semi_filtered_mafs):
    """Check for discrepancies between GWAS alleles and MAF reference alleles with swapping and fallback to GWAS."""
        
    # Merging on direct match (gnomAD)
    maf_lookup = pd.merge(
        full_gwas,
        semi_filtered_mafs[["ID", "REF", "ALT", "MAF"]],
        how="left",
        left_on=["SNP", "A1", "A2"],
        right_on=["ID", "REF", "ALT"],
        suffixes=("_GWAS", "_gnomAD")
    )
    
    # Merging on swapped match (gnomAD)
    maf_lookup_swapped = pd.merge(
        full_gwas,
        semi_filtered_mafs[["ID", "REF", "ALT", "MAF"]],
        how="left",
        left_on=["SNP", "A2", "A1"],  # Swap A1 and A2
        right_on=["ID", "REF", "ALT"],
        suffixes=("_GWAS", "_gnomAD_swapped")
    )
    
    # Create the final MAF column
    maf_lookup["MAF"] = maf_lookup["MAF_gnomAD"]
    maf_lookup["MAF_Source"] = "gnomAD"  # Default to gnomAD
    
    # Fallback: if MAF is missing from gnomAD, take from swapped gnomAD
    maf_lookup.loc[maf_lookup["MAF_gnomAD"].isnull(), "MAF"] = maf_lookup_swapped["MAF_gnomAD_swapped"]
    maf_lookup.loc[maf_lookup["MAF_gnomAD"].isnull(), "MAF_Source"] = "gnomAD_swapped"
    
    # Fallback to GWAS only if MAF is missing from both gnomAD and swapped gnomAD
    maf_lookup.loc[maf_lookup["MAF"].isnull(), "MAF"] = maf_lookup["MAF_GWAS"]
    maf_lookup.loc[maf_lookup["MAF"].isnull(), "MAF_Source"] = "GWAS"

    # Count the sources
    maf_from_gnomad = (maf_lookup["MAF_Source"] == "gnomAD").sum()
    maf_from_swapped = (maf_lookup["MAF_Source"] == "gnomAD_swapped").sum()
    maf_from_gwas = (maf_lookup["MAF_Source"] == "GWAS").sum()
    
    print(f"Number of SNPs with MAF from gnomAD: {maf_from_gnomad}")
    print(f"Number of SNPs with MAF from swapped gnomAD: {maf_from_swapped}")
    print(f"Number of SNPs with MAF from GWAS: {maf_from_gwas}")
    
    # Identifying SNPs with missing MAF in all sources
    missing_maf_snps = maf_lookup[maf_lookup["MAF"].isnull()]
    print(f"Number of SNPs missing MAF in all sources: {len(missing_maf_snps)}")

    return maf_lookup


def compute_clumps(r2, maf_lookup):
    """Compute SNP clumps and merge relevant data."""
    print("creating all possible mock clumps from ld table")
    clumps = r2.groupby('SNP_A')['SNP_B'].apply(lambda snps: [snp for snp in snps if snp != snps.name])
    clumps = pd.DataFrame(clumps).reset_index()
    clumps.columns = ['SNP', 'Clumped_SNPs']
    print("merging clump df")
    clumps_merged = clumps.merge(
        r2[["SNP_A", "CHR_A", "BP_A"]].drop_duplicates().rename(columns={"CHR_A": "CHR", "BP_A": "BP", "SNP_A": "SNP"}),
        on="SNP", how="left"
    )
    clumps_merged = clumps_merged.merge(maf_lookup[["SNP", "MAF"]], on="SNP")

    clumps_merged["Clump_Size"] = clumps_merged["Clumped_SNPs"].apply(len)
    print(clumps_merged.head())
    return clumps_merged


def calculate_clump_kb(row, snp_to_bp):
    """Calculate maximum distance in kb between lead SNP and its clumped SNPs."""
    lead_snp_bp = row["BP"]
    clump_snps = row["Clumped_SNPs"]
    
    if not clump_snps:
        return 0

    clump_bps = [snp_to_bp[snp] for snp in clump_snps if snp in snp_to_bp]
    
    if not clump_bps:
        return 0
    
    max_distance = max(abs(bp - lead_snp_bp) for bp in clump_bps)
    return max_distance / 1000  # Convert to kb


def change_empty_clump(row):
    """Ensure that clumps with TOTAL=0 have an empty list."""
    return [] if row["TOTAL"] == 0 else row["Clumped_SNPs"]


def match_clumps(orig_clumps, clumps_merged, kb_plusminus_pctg = 0.2, snps_plusminus=15, MAF_tolerance = 0.1):
    """Match original clumps to merged clumps based on multiple criteria."""
    matches_list = []
    
    for _, real_clump in orig_clumps.iterrows():
        clump_id = real_clump["SNP"]
        size_tolerance_kb = max(1, real_clump["Clump_Size_kb"] * kb_plusminus_pctg)
        size_tolerance_count = max(1, real_clump["Clump_Size"] * kb_plusminus_pctg)
        
        if real_clump["Clump_Size_kb"] == 0:
            matching_clumps = clumps_merged[clumps_merged["Clump_Size_kb"] == 0]
        else:
            matching_clumps = clumps_merged[
                (clumps_merged["Clump_Size_kb"].between(real_clump["Clump_Size_kb"] - size_tolerance_kb, 
                                                        real_clump["Clump_Size_kb"] + size_tolerance_kb)) &
                (clumps_merged["Clump_Size"].between(real_clump["Clump_Size"] - snps_plusminus, real_clump["Clump_Size"] + snps_plusminus)) &
                (clumps_merged["CHR"] == real_clump["CHR"]) & 
                (clumps_merged["MAF"].between(real_clump["MAF"] - MAF_tolerance, real_clump["MAF"] + MAF_tolerance))
            ]
        
        matching_clumps = matching_clumps.assign(Real_Clump_ID=clump_id)
        matches_list.append(matching_clumps)

    matches_df = pd.concat(matches_list, ignore_index=True)
    matches_df = matches_df[~matches_df["SNP"].isin(matches_df["Real_Clump_ID"].unique())]
    
    return matches_df

def filter_large_clumps(matches_df, threshold):
    """Filter clumps exceeding a user-defined threshold of matches - name is misleading"""
    total_real_clumps = matches_df["Real_Clump_ID"].nunique()
    
    large_clumps = matches_df["Real_Clump_ID"].value_counts()
    large_clumps = large_clumps[large_clumps > threshold].index
    
    filtered_df = matches_df[matches_df["Real_Clump_ID"].isin(large_clumps)]
    filtered_real_clumps = filtered_df["Real_Clump_ID"].nunique()

    print(f"Kept {filtered_real_clumps} out of {total_real_clumps} real clumps.")

    return filtered_df


def generate_mock_datasets(filtered_matches_df, n_mock_datasets):
    """Generate mock datasets by sampling unique SNPs per Real_Clump_ID."""
    all_mock_datasets = []
    
    # Track globally sampled clumps
    sampled_clumps_global = {real_clump_id: set() for real_clump_id in filtered_matches_df["Real_Clump_ID"].unique()}
    
    for dataset_idx in range(n_mock_datasets):
        print(f"Creating mock dataset {dataset_idx + 1}/{n_mock_datasets}...")
        sampled_clumps_this_perm = set()
        mock_dataset = []
        
        # Group matches_df by Real_Clump_ID
        grouped_matches = filtered_matches_df.groupby("Real_Clump_ID")
        
        for real_clump_id, group in grouped_matches:
            eligible_matches = group[
                ~group["SNP"].isin(sampled_clumps_global[real_clump_id] | sampled_clumps_this_perm)
            ]
            
            if eligible_matches.empty:
                print(f"No eligible matches left for Real_Clump_ID {real_clump_id}. Skipping...")
                continue
            
            # Sample one SNP per clump
            sampled_clump = eligible_matches.sample(1, replace=False, random_state=dataset_idx)
            sampled_snp = sampled_clump["SNP"].values[0]
            
            mock_dataset.append({
                "Mock_Dataset_ID": dataset_idx + 1,
                "Real_Clump_ID": real_clump_id,
                "Sampled_SNP": sampled_snp,
                "Clump_Size": sampled_clump["Clump_Size"].values[0],
                "Clump_Size_kb": sampled_clump["Clump_Size_kb"].values[0],
                "Clumped_SNPs": sampled_clump["Clumped_SNPs"].values[0],
                "CHR": sampled_clump["CHR"].values[0],
            })
            
            sampled_clumps_global[real_clump_id].add(sampled_snp)
            sampled_clumps_this_perm.add(sampled_snp)
        
        all_mock_datasets.extend(mock_dataset)
    
    return pd.DataFrame(all_mock_datasets)


def save_pickles_datasets(mock_datasets_df, matches_df, filtered_matches_df, orig_clumps, pickle_dir):
    """Save mock datasets and filtered data."""
    save_dataframe(filtered_matches_df, f"{pickle_dir}/large_matches_df.pkl")
    save_dataframe(matches_df, f"{pickle_dir}/all_matches_df.pkl")
    save_dataframe(orig_clumps, f"{pickle_dir}/orig_clumps.pkl")
    save_dataframe(mock_datasets_df, f"{pickle_dir}/mock_datasets_df.pkl")
    mock_datasets_df.to_csv(f"{pickle_dir}/mock_datasets_df.csv")
    orig_clumps.to_csv(f"{pickle_dir}/orig_clumps.csv")

def export_mock_datasets(mock_datasets_df, base_dir):
    """Export mock datasets to TSV files."""
    output_dir = os.path.join(base_dir, "mock_datasets")
    
    # Check if directory exists before creating it
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for i, group in mock_datasets_df.groupby("Mock_Dataset_ID"):
        print(f"Processing dataset {i}")
        df = group.copy()
        df["Clumped_SNPs"] = df["Clumped_SNPs"].apply(lambda b: str(b)[1:-1].replace(" ", "").replace("'", ""))
        df.to_csv(f"{output_dir}/mock_dataset{i}.tsv", sep="\t", index=False)






# ---------------------- Main Execution ---------------------- #
def main():
    """Main function to execute SNP processing pipeline."""
    args = parse_args()
    base_dir = args.base_dir
    print("creating paths")
    # Define directories
    pre_clumping_and_ld_dir = f"{base_dir}/pre_clumping_and_ld"
    clumping_and_vcfs_outputs_dir = f"{base_dir}/clumping_and_vcfs_outputs"
    pickle_dir = f"{base_dir}/pickle_dir"
    print(pickle_dir)
    if not os.path.exists(pickle_dir):
        os.makedirs(pickle_dir, exist_ok=True)

    print("Loading data...")
    # Pass use_pickles flag to load functions. These are Datasets that dont change much - so we pickle them after the first run
    #Later on we can use these pickles if we modify the code \ for downstream analysis
    r2 = load_ld_data(pre_clumping_and_ld_dir, pickle_dir, args.use_pickles)
    orig_clumps = load_clumps(clumping_and_vcfs_outputs_dir, pickle_dir, args.use_pickles)
    semi_filtered_mafs = load_maf(pre_clumping_and_ld_dir, pickle_dir, args.use_pickles)
    full_gwas = load_gwas(pre_clumping_and_ld_dir, pickle_dir, args.use_pickles)
    
    # Process and save data only if not using pickles
    maf_lookup = merge_maf_gwas(full_gwas, semi_filtered_mafs)
    clumps_merged = compute_clumps(r2, maf_lookup)
    print("creating snp to bp dict")
    snp_to_bp = clumps_merged.set_index("SNP")["BP"].to_dict()
    clumps_merged["Clump_Size_kb"] = clumps_merged.apply(lambda row: calculate_clump_kb(row, snp_to_bp), axis=1)
    save_dataframe(clumps_merged, f"{pickle_dir}/clumps_merged.pkl")
    
    # Save intermediate data to pickle if not using pickles

    print("Raw data processed and saved to pickle files.")
    
    # Apply size to original clumps
    orig_clumps = orig_clumps.merge(maf_lookup[["SNP", "MAF"]], on = "SNP", how = "left")
    orig_clumps["Clumped_SNPs"] = orig_clumps["SP2"].apply(lambda a: a.split(","))
    orig_clumps["Clumped_SNPs"] = orig_clumps.apply(change_empty_clump, axis=1)
    orig_clumps["Clump_Size_kb"] = orig_clumps.apply(calculate_clump_kb, args=(snp_to_bp,), axis=1)
    orig_clumps["Clump_Size"] = orig_clumps["Clumped_SNPs"].apply(len)

    print("creating null matches to each clump")
    print(orig_clumps.head())
    matches_df = match_clumps(orig_clumps, clumps_merged)
    print(f"All matches: {len(matches_df)}")
    # Assuming `matches_per_clump` is a column or variable that tracks the number of matches per clump
    filtered_matches_df = filter_large_clumps(matches_df, args.n_mock_datasets)
    print("generating mock datasets")
    # Generate mock datasets
    mock_datasets_df = generate_mock_datasets(filtered_matches_df, args.n_mock_datasets)

    # Save and export data
    save_pickles_datasets(mock_datasets_df, matches_df, filtered_matches_df, orig_clumps, pickle_dir)
    export_mock_datasets(mock_datasets_df, base_dir)


if __name__ == "__main__":
    main()
