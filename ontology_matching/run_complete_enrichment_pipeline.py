import os
import re
import glob
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import mannwhitneyu

from functions_for_snp_gsealike import (
    compute_enrichment_score,
    permutation_test,
    statistical_tests,
)

# =====================================================================
# CONFIG
# =====================================================================

BASE_DIR = Path("") #config needed
RESULTS_DIR = BASE_DIR / "results"
CODE_DIR = BASE_DIR / "enrichment_pipeline_code"

# Collapsed (ancestor-expanded + filtered) mappings created by ancestor_and_collapse.py
MAPPING_DIR = CODE_DIR / "ontology_annotation_mappings"

COLLAPSED_MAPPING_FILES = {
    "EFO": MAPPING_DIR / "efo_term_to_annotations_collapsed.csv",
    "BTO": MAPPING_DIR / "bto_term_to_annotations_collapsed.csv",
    "CL":  MAPPING_DIR / "cl_term_to_annotations_collapsed.csv",
    "CLO": MAPPING_DIR / "clo_term_to_annotations_collapsed.csv",
}

# Minimum number of annotations per term to test
MIN_CELLS_PER_TERM = 10

# GWAS meta (for pretty labels)
GWAS_META_PATH = "/path/to/gwas_study_metainfo.xlsx" #provided in the paper


# =====================================================================
# GWAS META LABELS (for pretty trait names)
# =====================================================================

def _load_gwas_meta():
    meta = pd.read_excel(GWAS_META_PATH)
    meta.columns = [c.strip() for c in meta.columns]
    meta["study_str"] = meta["study"].astype(str)
    meta["study_lower"] = meta["study_str"].str.lower()
    return meta


def _parse_author_year(gwas_name: str):
    """
    From GWAS name like 'D2_ldlp_willer_2013' or 'D_asd_grove_2019'
    extract:
      author = surname token before year
      year   = 4-digit year
      code   = everything before the author+year (trait-ish code)
    """
    base = re.sub(r"^D2?_", "", gwas_name)  # drop D_ or D2_
    parts = base.split("_")
    if len(parts) < 3:
        return None, None, None

    year = parts[-1]
    author = parts[-2].lower()
    code = "_".join(parts[:-2]).lower()

    return author, year, code


def _find_meta_row_for_gwas(meta, author: str, year: str, code: str = None):
    if author is None or year is None:
        return None

    mask_year = meta["study_str"].str.contains(year)

    # author filter (avoid being too aggressive on very short names)
    if author:
        author_token = author + " "
        mask_author = meta["study_lower"].str.contains(author_token.lower())
    else:
        mask_author = True

    cand = meta[mask_year & mask_author]
    if len(cand) == 0:
        return None
    if len(cand) == 1:
        return cand.iloc[0]

    # tie-break using code vs abbreviation / trait
    if code:
        code = code.lower()
        cand = cand.copy()
        cand["abbr_str"] = cand["abbreviation"].astype(str).str.lower()
        cand["trait_str"] = cand["trait"].astype(str).str.lower()

        mask_abbr = (cand["abbr_str"] == code) | \
                    (cand["abbr_str"].str.contains(code))

        mask_trait = cand["trait_str"].str.contains(code)

        cand2 = cand[mask_abbr | mask_trait]
        if len(cand2) == 1:
            return cand2.iloc[0]
        elif len(cand2) > 1:
            cand = cand2

    return cand.iloc[0]


def build_gwas_label_map(gwas_list):
    meta = _load_gwas_meta()
    label_map = {}
    for g in gwas_list:
        if not (g.startswith("D") or g.startswith("D2")):
            continue

        author, year, code = _parse_author_year(g)
        row = _find_meta_row_for_gwas(meta, author, year, code=code)
        if row is None:
            print(f"[GWAS LABEL WARN] No meta match for {g} (author={author}, year={year}, code={code})")
            continue

        trait = str(row["trait"]).strip()
        abbr = str(row["abbreviation"]).strip()
        pretty = f"{abbr} ({trait})" if abbr and abbr.lower() != "nan" else trait
        label_map[g] = pretty

    print("\n[GWAS LABEL MAP] example entries:")
    for k in list(label_map.keys())[:8]:
        print(f"  {k} -> {label_map[k]}")
    return label_map


# =====================================================================
# GWAS FILTERING (mocks/clumps)
# =====================================================================

def count_vcf_rows(vcf_path: str) -> int:
    """Count non-header lines in a VCF file."""
    if not os.path.exists(vcf_path):
        return 0
    count = 0
    with open(vcf_path) as f:
        for line in f:
            if not line.startswith("#"):
                count += 1
    return count


def count_clumps(clumped_path: str) -> int:
    """Count number of clumps in PLINK .clumped file."""
    if not os.path.exists(clumped_path):
        return 0
    df = pd.read_csv(clumped_path, sep="\t", comment="#", engine="python")
    return df.shape[0]


def count_processed_mocks(stream_dir: str) -> int:
    """
    Count mock IDs that have BOTH:
      - mock_dataset_<ID>_contrib.npz
      - dataset_<ID>_complete.txt
    in streaming_results/.
    """
    if not os.path.exists(stream_dir):
        return 0

    mock_files = glob.glob(os.path.join(stream_dir, "mock_dataset_*_contrib.npz"))
    mock_ids = set()
    for path in mock_files:
        m = re.search(r"mock_dataset_(\d+)_contrib\.npz$", os.path.basename(path))
        if m:
            mock_ids.add(m.group(1))

    complete_files = glob.glob(os.path.join(stream_dir, "dataset_*_complete.txt"))
    complete_ids = set()
    for path in complete_files:
        m = re.search(r"dataset_(\d+)_complete\.txt$", os.path.basename(path))
        if m:
            complete_ids.add(m.group(1))

    processed_ids = mock_ids & complete_ids
    return len(processed_ids)


def build_filtered_gwas_list() -> list[str]:
    """
    Reproduce your original filtering:

    1. For each GWAS dir in RESULTS_DIR, compute:
         - 'VCF mock subdirs'
         - 'Processed mocks (streaming)'
         - 'Clumps'
         - 'SNPs in final VCF'
    2. For GWAS starting with 'D':
         keep only those with:
             Processed mocks (streaming) > 196
             Clumps > 65
    3. Add all non-'D' GWAS as-is.
    """
    summary = []

    for gwas in sorted(os.listdir(RESULTS_DIR)):
        gwas_path = os.path.join(RESULTS_DIR, gwas)
        if not os.path.isdir(gwas_path):
            continue

        # 1) mock VCF subdirs
        vcfs_mock_dir = os.path.join(gwas_path, "vcfs_mock_datasets")
        vcfs_mock_subdirs = 0
        if os.path.exists(vcfs_mock_dir):
            vcfs_mock_subdirs = len(
                [
                    d
                    for d in os.listdir(vcfs_mock_dir)
                    if os.path.isdir(os.path.join(vcfs_mock_dir, d))
                ]
            )

        # 2) processed mocks
        stream_dir = os.path.join(gwas_path, "streaming_results")
        processed_mocks = count_processed_mocks(stream_dir)

        # 3) clumps
        clumps_file = os.path.join(
            gwas_path,
            "clumping_and_vcfs_outputs",
            "clumps_cleaned.clumped",
        )
        clump_count = count_clumps(clumps_file)

        # 4) SNPs in final VCF
        final_vcf = os.path.join(
            gwas_path,
            "clumping_and_vcfs_outputs",
            "main_gwas_snps_correct_ref_final.vcf",
        )
        snp_count = count_vcf_rows(final_vcf)

        summary.append(
            {
                "GWAS": gwas,
                "VCF mock subdirs": vcfs_mock_subdirs,
                "Processed mocks (streaming)": processed_mocks,
                "Clumps": clump_count,
                "SNPs in final VCF": snp_count,
            }
        )

    df = pd.DataFrame(summary)

    print("\n=== GWAS SUMMARY (sorted by Clumps, Processed mocks) ===")
    if not df.empty:
        print(
            df.sort_values(
                by=["Clumps", "Processed mocks (streaming)"]
            ).to_string(index=False)
        )
    else:
        print("No GWAS directories found in RESULTS_DIR.")

    d_res = df[df["GWAS"].str.startswith("D")]
    notd_res = df[~df["GWAS"].str.startswith("D")]

    d_res_valid = d_res[
        (d_res["Processed mocks (streaming)"] > 196) #just empirical things that failed - FIX
        & (d_res["Clumps"] > 50)
    ].reset_index(drop=True)

    gwass_valid = pd.concat([d_res_valid, notd_res], ignore_index=True)
    gwass = list(gwass_valid["GWAS"])

    print("\n=== Filtered GWAS list (used for enrichment) ===")
    print(gwass)

    return gwass


# =====================================================================
# COLLAPSED TERM→ANNOTATION MAPPINGS
# =====================================================================

def load_collapsed_mapping(ontology_name: str):
    """
    Load collapsed mapping saved by ancestor_and_collapse.py:

      Ontology,Term_ID,Term_Name,Annotation_Index

    Returns:
      term_to_ann : dict term_id -> set(annotation_indices)
      term_to_name: dict term_id -> term_name
    """
    mapping_file = COLLAPSED_MAPPING_FILES.get(ontology_name)
    if mapping_file is None or not mapping_file.exists():
        raise FileNotFoundError(
            f"No collapsed mapping file for {ontology_name}: {mapping_file}"
        )

    df = pd.read_csv(mapping_file)
    if df.empty:
        print(f"[WARN] Collapsed mapping for {ontology_name} is empty: {mapping_file}")
        return {}, {}

    term_to_ann = defaultdict(set)
    term_to_name = {}

    for _, row in df.iterrows():
        term_id = str(row["Term_ID"])
        ann_idx = int(row["Annotation_Index"])
        term_to_ann[term_id].add(ann_idx)
        term_to_name[term_id] = str(row.get("Term_Name", term_id))

    print(
        f"[{ontology_name}] Loaded collapsed mapping: "
        f"{len(term_to_ann)} terms, "
        f"{df['Annotation_Index'].nunique()} unique annotations"
    )
    return term_to_ann, term_to_name


def build_ontology_masks_from_collapsed(
    term_to_ann_indices: dict,
    term_to_name: dict,
    annotation_indices: np.ndarray,
):
    """
    Turn term→annotation_indices into ontology_masks compatible with
    run_ontology_enrichment.

    annotation_indices: array of Annotation_Index that appear in the pvalues DF
                        (for a given GWAS).
    """
    annotation_indices = np.array(annotation_indices, dtype=int)
    all_ann_set = set(annotation_indices)

    ontology_masks = {}

    for term_id, inds in term_to_ann_indices.items():
        positive = sorted(all_ann_set.intersection(inds))
        positive_count = len(positive)
        negative_count = len(all_ann_set) - positive_count

        # boolean mask over annotation_indices
        mask = np.isin(annotation_indices, positive)

        ontology_masks[term_id] = {
            "mask": mask,
            "term_name": term_to_name.get(term_id, term_id),
            "positive_indices": positive,
            "positive_count": positive_count,
            "negative_count": negative_count,
        }

    return ontology_masks


# =====================================================================
# ENRICHMENT UTILITIES
# =====================================================================

def compute_es_nes_for_term(pvalues_df, term_indices, pval_col):
    """
    Use the existing GSEA-like functions to get ES, NES and ES_pvalue
    for a given ontology term (a set of Annotation_Index values).
    """
    df = pvalues_df.copy()

    if "Annotation_Index" in df.columns:
        df = df.set_index("Annotation_Index")

    if "Empirical_P_Upper" not in df.columns:
        df["Empirical_P_Upper"] = df[pval_col]

    enrichment_data = compute_enrichment_score(
        df, np.array(term_indices, dtype=int)
    )
    perm_res = permutation_test(enrichment_data)
    stats_res = statistical_tests(enrichment_data, df)

    ES = enrichment_data["max_ES"]
    NES = perm_res["NES"]
    ES_pval = perm_res["ES_pvalue"]

    return ES, NES, ES_pval, stats_res


def run_ontology_enrichment(
    ontology_masks,
    pvalues_df,
    gwas_name,
    ontology_name,
    min_cells=MIN_CELLS_PER_TERM,
    save_mappings=True,
    gwas_pretty=None,
):
    """
    Run enrichment analysis and (optionally) produce term→annotation mappings.

    ontology_masks: dict term_id -> {
        'mask', 'term_name', 'positive_indices',
        'positive_count', 'negative_count'
    }
    """
    pretty = gwas_pretty or gwas_name
    print(f"\n{'='*60}")
    print(f"Running enrichment for {gwas_name} [{pretty}] - {ontology_name}")
    print(f"{'='*60}")

    results_list = []
    mappings_list = []

    # Determine p-value column
    if "Empirical_P_Upper" in pvalues_df.columns:
        pval_col = "Empirical_P_Upper"
    elif "empirical_pvalue" in pvalues_df.columns:
        pval_col = "empirical_pvalue"
    elif "P_Value" in pvalues_df.columns:
        pval_col = "P_Value"
    else:
        print("ERROR: Could not find p-value column in pvalues_df.")
        return pd.DataFrame(), None

    print(f"Using p-value column: {pval_col}")

    total_unique_terms = len(ontology_masks)

    # filter by coverage
    testable_terms = {
        term_id: m_info
        for term_id, m_info in ontology_masks.items()
        if m_info["positive_count"] >= min_cells
    }

    tested_count = len(testable_terms)

    print("Terms Filtering (Statistical):")
    print(f"  N_min (annotations per term): {min_cells}")
    print(f"  Terms dropped due to low coverage: {total_unique_terms - tested_count}")
    print(f"  Terms to test: {tested_count}")

    for term_id, mask_info in tqdm(testable_terms.items(), desc="Testing terms"):
        try:
            term_indices = mask_info["positive_indices"]

            term_pvalues = pvalues_df[
                pvalues_df["Annotation_Index"].isin(term_indices)
            ]
            other_pvalues = pvalues_df[
                ~pvalues_df["Annotation_Index"].isin(term_indices)
            ]

            if len(term_pvalues) < min_cells:
                continue

            # 1) Mann-Whitney U (less = more enriched)
            stat, mwu_p = mannwhitneyu(
                term_pvalues[pval_col],
                other_pvalues[pval_col],
                alternative="less",
            )

            # 2) ES / NES using existing GSEA-like functions
            ES, NES, ES_pval, stats_res = compute_es_nes_for_term(
                pvalues_df=pvalues_df,
                term_indices=term_indices,
                pval_col=pval_col,
            )

            results_list.append(
                {
                    "GWAS": gwas_name,
                    "GWAS_Pretty": pretty,
                    "Ontology": ontology_name,
                    "Term_ID": term_id,
                    "Term_Name": mask_info["term_name"],
                    "N_Positive": mask_info["positive_count"],
                    "N_Negative": mask_info["negative_count"],
                    "MWU_Statistic": stat,
                    "MWU_PValue": mwu_p,
                    "ES": ES,
                    "NES": NES,
                    "ES_PValue": ES_pval,
                    "MWU_PValue_SNPDist": stats_res["mw_pvalue"],
                    "Fisher_PValue_SNPDist": stats_res["fisher_pvalue"],
                }
            )

            if save_mappings:
                for ann_idx in term_indices:
                    mappings_list.append(
                        {
                            "GWAS": gwas_name,
                            "GWAS_Pretty": pretty,
                            "Ontology": ontology_name,
                            "Term_ID": term_id,
                            "Term_Name": mask_info["term_name"],
                            "Annotation_Index": ann_idx,
                        }
                    )

        except Exception as e:
            print(f"[WARN] Error testing term {term_id}: {e}")
            continue

    results_df = pd.DataFrame(results_list)
    mappings_df = pd.DataFrame(mappings_list) if save_mappings else None

    if not results_df.empty:
        results_df = results_df.sort_values("MWU_PValue", ascending=True)
        print("\nResults summary:")
        print(f"  Terms tested: {len(results_df)}")
        print(f"  Significant (MWU p<0.05): {(results_df['MWU_PValue'] < 0.05).sum()}")

        print("\nTop 20 enriched terms (by MWU p):")
        top_20 = results_df.head(20)
        for rank, (_, row) in enumerate(top_20.iterrows(), start=1):
            print(
                f"  {rank:2d}. {row['Term_Name'][:50]:50s} "
                f"| p={row['MWU_PValue']:.2e} "
                f"| N={row['N_Positive']} "
                f"| NES={row['NES']:.3f}"
            )

    return results_df, mappings_df


# =====================================================================
# PER-GWAS PROCESSING
# =====================================================================

def process_single_gwas(
    gwas_name,
    ontology_name,
    term_to_ann_indices,
    term_to_name,
    gwas_label_map=None,
    min_cells=MIN_CELLS_PER_TERM,
    save_mappings=True,
):
    """
    Process a single GWAS with a single ontology using pre-collapsed mappings.
    """
    gwas_pretty = (
        gwas_label_map.get(gwas_name, gwas_name) if gwas_label_map else gwas_name
    )
    print(f"\n{'#' * 80}")
    print(f"# Processing: {gwas_name} [{gwas_pretty}] - {ontology_name}")
    print(f"{'#' * 80}")

    # Load p-values
    pval_file = (
        RESULTS_DIR
        / gwas_name
        / "sei_score_aggregations_all"
        / "empirical_pvalues50.csv"
    )
    if not pval_file.exists():
        print(f"[WARN] P-value file not found: {pval_file}")
        return None, None

    pvalues_df = pd.read_csv(pval_file)
    print(f"Loaded p-values: {len(pvalues_df)} annotations")

    # We assume pvalues_df has 'Annotation_Index'
    if "Annotation_Index" not in pvalues_df.columns:
        raise ValueError("pvalues_df must contain 'Annotation_Index' column")

    # Restrict to annotations that appear in this ontology mapping
    all_ann_with_terms = set()
    for inds in term_to_ann_indices.values():
        all_ann_with_terms.update(inds)

    pvalues_filtered = pvalues_df[
        pvalues_df["Annotation_Index"].isin(all_ann_with_terms)
    ].copy()

    print(
        f"Filtered p-values to annotations that have ontology terms: "
        f"{len(pvalues_filtered)} / {len(pvalues_df)}"
    )
    if pvalues_filtered.empty:
        print("[WARN] No overlapping annotations between GWAS and ontology mapping.")
        return None, None

    # Build ontology masks for this GWAS
    ontology_masks = build_ontology_masks_from_collapsed(
        term_to_ann_indices=term_to_ann_indices,
        term_to_name=term_to_name,
        annotation_indices=pvalues_filtered["Annotation_Index"].values,
    )

    results_df, mappings_df = run_ontology_enrichment(
        ontology_masks=ontology_masks,
        pvalues_df=pvalues_filtered,
        gwas_name=gwas_name,
        ontology_name=ontology_name,
        min_cells=min_cells,
        save_mappings=save_mappings,
        gwas_pretty=gwas_pretty,
    )

    return results_df, mappings_df


# =====================================================================
# PIPELINE DRIVER
# =====================================================================

def run_complete_pipeline(
    gwas_list=None,
    ontologies=None,
    min_cells=MIN_CELLS_PER_TERM,
):
    """
    Run enrichment for all GWAS × ontology combinations using pre-collapsed mappings.
    """
    if ontologies is None:
        ontologies = list(COLLAPSED_MAPPING_FILES.keys())

    print("\n" + "=" * 80)
    print("ONTOLOGY ENRICHMENT PIPELINE (using collapsed mappings)")
    print("=" * 80)
    print(f"Ontologies: {ontologies}")
    print(f"Min cells per term: {min_cells}")
    print("=" * 80 + "\n")

    # Build/filter GWAS list
    if gwas_list is None:
        gwas_list = build_filtered_gwas_list()

    print(f"\nTotal GWAS to process: {len(gwas_list)}\n")

    # GWAS pretty labels
    gwas_label_map = build_gwas_label_map(gwas_list)
    globals()["gwas_label_map"] = gwas_label_map


    # Load all ontology mappings once
    ontology_mappings = {}
    for ont_name in ontologies:
        term_to_ann, term_to_name = load_collapsed_mapping(ont_name)
        ontology_mappings[ont_name] = (term_to_ann, term_to_name)

    all_results = []
    all_mappings = []

    for gwas_name in gwas_list:
        for ont_name in ontologies:
            try:
                term_to_ann, term_to_name = ontology_mappings[ont_name]
                if not term_to_ann:
                    print(f"[WARN] No term mappings for ontology {ont_name}; skipping.")
                    continue

                results_df, mappings_df = process_single_gwas(
                    gwas_name=gwas_name,
                    ontology_name=ont_name,
                    term_to_ann_indices=term_to_ann,
                    term_to_name=term_to_name,
                    gwas_label_map=gwas_label_map,
                    min_cells=min_cells,
                    save_mappings=True,
                )

                if results_df is not None and not results_df.empty:
                    all_results.append(results_df)

                    out_res = (
                        RESULTS_DIR
                        / gwas_name
                        / f"ontology_enrichment_{ont_name.lower()}.csv"
                    )
                    out_res.parent.mkdir(parents=True, exist_ok=True)
                    results_df.to_csv(out_res, index=False)
                    print(f"Saved enrichment results to: {out_res}")

                if mappings_df is not None and not mappings_df.empty:
                    all_mappings.append(mappings_df)
                    out_map = (
                        RESULTS_DIR
                        / gwas_name
                        / f"term_to_annotation_mappings_{ont_name.lower()}.csv"
                    )
                    mappings_df.to_csv(out_map, index=False)

                    print(f"Saved term→annotation mappings to: {out_map}")

            except Exception as e:
                print(f"[ERROR] Processing {gwas_name} - {ont_name}: {e}")
                import traceback
                traceback.print_exc()
                continue

    # Combine enrichment results
    if all_results:
        combined_results = pd.concat(all_results, ignore_index=True)
        out_combined = CODE_DIR / "ontology_matching_res" /"ontology_enrichment_all_results.csv"
        combined_results.to_csv(out_combined, index=False)
        print(f"\nSaved combined enrichment results to: {out_combined}")
    else:
        combined_results = None
        print("\nNo enrichment results collected.")

    # Combine mappings (optional; can be handy)
    if all_mappings:
        combined_mappings = pd.concat(all_mappings, ignore_index=True)
        out_combined_map = CODE_DIR / "term_to_annotation_mappings_all.csv"
        combined_mappings.to_csv(out_combined_map, index=False)
        print(f"Saved combined term→annotation mappings to: {out_combined_map}")

    return combined_results


# =====================================================================
# MAIN
# =====================================================================

def main():
    ontologies_to_use = ["EFO", "CL", "BTO", "CLO"]
    _ = run_complete_pipeline(
        gwas_list=None,          # let it build + filter as before
        ontologies=ontologies_to_use,
        min_cells=MIN_CELLS_PER_TERM,
    )


if __name__ == "__main__":
    main()
