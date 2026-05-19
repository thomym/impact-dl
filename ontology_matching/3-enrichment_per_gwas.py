"""Step 3 of ontology_matching/ — per-GWAS ontology enrichment.

For one GWAS (via --gwas_name), test each collapsed ontology term for
enrichment among the most significant Sei annotations using:
  - one-sided Mann-Whitney U  on Empirical_P_Upper (mwu_pvalue)
  - GSEA-like running ES + NES + permutation p-value (seeded)
  - secondary MWU / Fisher reported by `functions_for_snp_gsealike`

Inputs:
    ${ontology_collapsed_dir}/<ont>_term_to_annotations_collapsed.csv   (from step 2)
    ${results_root}/${gwas_name}/sei_score_aggregations_all/empirical_pvalues50.csv
    ${results_root}/${gwas_name}/vcfs_mock_datasets/                    (sanity check)
    ${results_root}/${gwas_name}/clumping_and_vcfs_outputs/clumps_cleaned_noncoding.clumped (sanity check)

Outputs:
    ${results_root}/${gwas_name}/ontology_enrichment_<ont>.csv
    ${results_root}/${gwas_name}/term_to_annotation_mappings_<ont>.csv
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Set, Tuple

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from paths import load_paths  # noqa: E402

# Same-dir helper; uppercase identifier is fine here.
from functions_for_snp_gsealike import (  # noqa: E402
    compute_enrichment_score,
    permutation_test,
    statistical_tests,
)


ONTOLOGIES_DEFAULT = ["EFO", "BTO", "CL", "CLO"]


# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------

def assert_gwas_ready(gwas_dir: Path, n_mock_datasets: int) -> Path:
    """Strict pre-flight: refuse to run unless step 3 (prediction_pipeline) is complete."""
    pvalues_file = gwas_dir / "sei_score_aggregations_all" / "empirical_pvalues50.csv"
    if not pvalues_file.exists():
        raise FileNotFoundError(
            f"empirical_pvalues50.csv not found at {pvalues_file}. "
            "Run prediction_pipeline/4-empirical_pvalues.py first."
        )

    vcfs_mock_dir = gwas_dir / "vcfs_mock_datasets"
    mock_subdirs = (
        sum(1 for p in vcfs_mock_dir.iterdir() if p.is_dir())
        if vcfs_mock_dir.exists() else 0
    )
    if mock_subdirs != n_mock_datasets:
        raise AssertionError(
            f"Expected {n_mock_datasets} mock VCF subdirs under {vcfs_mock_dir}, "
            f"found {mock_subdirs}. Re-run pre_prediction / prediction_pipeline."
        )

    clumps_file = gwas_dir / "clumping_and_vcfs_outputs" / "clumps_cleaned_noncoding.clumped"
    n_clumps = 0
    if clumps_file.exists():
        n_clumps = pd.read_csv(clumps_file, sep="\t", comment="#", engine="python").shape[0]
    if n_clumps < 50:
        raise AssertionError(
            f"Only {n_clumps} non-coding clumps in {clumps_file} (need ≥50). "
            "GWAS likely lacks enough independent signals."
        )

    print(f"[OK] sanity: {mock_subdirs} mock VCF subdirs, {n_clumps} non-coding clumps")
    return pvalues_file


# ---------------------------------------------------------------------------
# Collapsed mapping loader
# ---------------------------------------------------------------------------

def load_collapsed_mapping(
    collapsed_dir: Path, ontology_name: str
) -> Tuple[Dict[str, Set[int]], Dict[str, str]]:
    mapping_file = collapsed_dir / f"{ontology_name.lower()}_term_to_annotations_collapsed.csv"
    if not mapping_file.exists():
        raise FileNotFoundError(f"No collapsed mapping for {ontology_name}: {mapping_file}")

    df = pd.read_csv(mapping_file)
    term_to_ann: Dict[str, Set[int]] = defaultdict(set)
    term_to_name: Dict[str, str] = {}
    for _, row in df.iterrows():
        term_id = str(row["Term_ID"])
        term_to_ann[term_id].add(int(row["Annotation_Index"]))
        term_to_name[term_id] = str(row.get("Term_Name", term_id))

    print(
        f"[{ontology_name}] Loaded collapsed mapping: "
        f"{len(term_to_ann)} terms, "
        f"{df['Annotation_Index'].nunique()} unique annotations"
    )
    return term_to_ann, term_to_name


def build_ontology_masks(
    term_to_ann: Dict[str, Set[int]],
    term_to_name: Dict[str, str],
    annotation_indices: np.ndarray,
) -> Dict[str, dict]:
    annotation_indices = np.asarray(annotation_indices, dtype=int)
    all_ann = set(annotation_indices.tolist())
    masks = {}
    for term_id, inds in term_to_ann.items():
        positive = sorted(all_ann & inds)
        masks[term_id] = {
            "mask": np.isin(annotation_indices, positive),
            "term_name": term_to_name.get(term_id, term_id),
            "positive_indices": positive,
            "positive_count": len(positive),
            "negative_count": len(all_ann) - len(positive),
        }
    return masks


# ---------------------------------------------------------------------------
# Enrichment
# ---------------------------------------------------------------------------

def compute_es_nes_for_term(pvalues_df, term_indices, pval_col, seed):
    df = pvalues_df.copy()
    if "Annotation_Index" in df.columns:
        df = df.set_index("Annotation_Index")
    if "Empirical_P_Upper" not in df.columns:
        df["Empirical_P_Upper"] = df[pval_col]
    enrichment = compute_enrichment_score(df, np.array(term_indices, dtype=int))
    perm = permutation_test(enrichment, seed=seed)
    stats_res = statistical_tests(enrichment, df)
    return enrichment["max_ES"], perm["NES"], perm["ES_pvalue"], stats_res


def run_ontology_enrichment(
    ontology_masks,
    pvalues_df,
    gwas_name: str,
    ontology_name: str,
    min_cells: int,
    seed: int,
    save_mappings: bool = True,
):
    print(f"\n{'='*60}\nRunning enrichment for {gwas_name} - {ontology_name}\n{'='*60}")

    if "Empirical_P_Upper" not in pvalues_df.columns:
        raise ValueError("pvalues_df missing required column 'Empirical_P_Upper'")
    pval_col = "Empirical_P_Upper"

    testable = {
        t: m for t, m in ontology_masks.items() if m["positive_count"] >= min_cells
    }
    print(
        f"  N_min annotations: {min_cells} | "
        f"terms dropped low-coverage: {len(ontology_masks) - len(testable)} | "
        f"terms to test: {len(testable)}"
    )

    results = []
    mappings = []

    for term_id, info in tqdm(testable.items(), desc="Testing terms"):
        try:
            term_idx = info["positive_indices"]
            term_p = pvalues_df[pvalues_df["Annotation_Index"].isin(term_idx)]
            other_p = pvalues_df[~pvalues_df["Annotation_Index"].isin(term_idx)]
            if len(term_p) < min_cells:
                continue

            stat, mwu_p = mannwhitneyu(
                term_p[pval_col], other_p[pval_col], alternative="less"
            )
            ES, NES, ES_pval, stats_res = compute_es_nes_for_term(
                pvalues_df=pvalues_df, term_indices=term_idx, pval_col=pval_col, seed=seed,
            )

            results.append({
                "GWAS": gwas_name,
                "Ontology": ontology_name,
                "Term_ID": term_id,
                "Term_Name": info["term_name"],
                "N_Positive": info["positive_count"],
                "N_Negative": info["negative_count"],
                "MWU_Statistic": stat,
                "MWU_PValue": mwu_p,
                "ES": ES,
                "NES": NES,
                "ES_PValue": ES_pval,
                "MWU_PValue_SNPDist": stats_res["mw_pvalue"],
                "Fisher_PValue_SNPDist": stats_res["fisher_pvalue"],
            })
            if save_mappings:
                for ann_idx in term_idx:
                    mappings.append({
                        "GWAS": gwas_name,
                        "Ontology": ontology_name,
                        "Term_ID": term_id,
                        "Term_Name": info["term_name"],
                        "Annotation_Index": ann_idx,
                    })
        except Exception as e:
            print(f"[WARN] Term {term_id} failed: {e}")
            continue

    results_df = pd.DataFrame(results)
    if not results_df.empty:
        results_df = results_df.sort_values("MWU_PValue", ascending=True)
        print(
            f"  -> {len(results_df)} terms tested, "
            f"{(results_df['MWU_PValue'] < 0.05).sum()} significant @ MWU p<0.05"
        )
    mappings_df = pd.DataFrame(mappings) if save_mappings else None
    return results_df, mappings_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    p.add_argument("--paths_yaml", default=None)
    p.add_argument("--gwas_name", required=True)
    p.add_argument("--ontologies", nargs="+", default=None,
                   help=f"Ontology acronyms to process. Default: {ONTOLOGIES_DEFAULT}")
    p.add_argument("--min_cells", type=int, default=None,
                   help="Minimum annotations per term to test. Default: paths.yaml/min_cells_per_term")
    p.add_argument("--seed", type=int, default=42,
                   help="RNG seed for permutation_test (NES/ES_PValue reproducibility)")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    paths = load_paths(args.paths_yaml)
    ontologies = args.ontologies or ONTOLOGIES_DEFAULT
    min_cells = args.min_cells if args.min_cells is not None else int(paths["min_cells_per_term"])

    gwas_dir = Path(paths["results_root"]) / args.gwas_name
    collapsed_dir = Path(paths["ontology_collapsed_dir"])
    n_mock_datasets = int(paths["n_mock_datasets"])

    print(f"3- config:")
    print(f"  gwas_name:              {args.gwas_name}")
    print(f"  gwas_dir:               {gwas_dir}")
    print(f"  collapsed_dir:          {collapsed_dir}")
    print(f"  ontologies:             {ontologies}")
    print(f"  min_cells_per_term:     {min_cells}")
    print(f"  seed:                   {args.seed}")

    pvalues_file = assert_gwas_ready(gwas_dir, n_mock_datasets)
    pvalues_df = pd.read_csv(pvalues_file)
    if "Annotation_Index" not in pvalues_df.columns:
        raise ValueError(f"{pvalues_file} missing required column 'Annotation_Index'")
    print(f"Loaded p-values: {len(pvalues_df)} annotations")

    for ont_name in ontologies:
        print(f"\n{'#' * 80}\n# {args.gwas_name} - {ont_name}\n{'#' * 80}")
        try:
            term_to_ann, term_to_name = load_collapsed_mapping(collapsed_dir, ont_name)
        except FileNotFoundError as e:
            print(f"[WARN] {e}; skipping")
            continue
        if not term_to_ann:
            print(f"[WARN] empty mapping for {ont_name}; skipping")
            continue

        # Restrict p-values to annotations covered by this ontology
        all_ann_with_terms = set().union(*term_to_ann.values())
        pvalues_filtered = pvalues_df[
            pvalues_df["Annotation_Index"].isin(all_ann_with_terms)
        ].copy()
        print(
            f"Filtered p-values to ontology-covered annotations: "
            f"{len(pvalues_filtered)} / {len(pvalues_df)}"
        )
        if pvalues_filtered.empty:
            print(f"[WARN] no overlap between p-values and {ont_name} mapping; skipping")
            continue

        masks = build_ontology_masks(
            term_to_ann=term_to_ann,
            term_to_name=term_to_name,
            annotation_indices=pvalues_filtered["Annotation_Index"].values,
        )

        results_df, mappings_df = run_ontology_enrichment(
            ontology_masks=masks,
            pvalues_df=pvalues_filtered,
            gwas_name=args.gwas_name,
            ontology_name=ont_name,
            min_cells=min_cells,
            seed=args.seed,
        )

        if results_df is not None and not results_df.empty:
            out = gwas_dir / f"ontology_enrichment_{ont_name.lower()}.csv"
            out.parent.mkdir(parents=True, exist_ok=True)
            results_df.to_csv(out, index=False)
            print(f"Saved enrichment results: {out}")

        if mappings_df is not None and not mappings_df.empty:
            out_map = gwas_dir / f"term_to_annotation_mappings_{ont_name.lower()}.csv"
            mappings_df.to_csv(out_map, index=False)
            print(f"Saved term→annotation mappings: {out_map}")


if __name__ == "__main__":
    main()
