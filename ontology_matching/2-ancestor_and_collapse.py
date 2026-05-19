"""Step 2 of ontology_matching/ — expand matched terms with ancestors, then
collapse nested terms with identical annotation coverage.

GWAS-independent: run once per refresh of the matched CSVs from step 1.

Inputs:
    ${ontology_matched_dir}/matched_sup2_<ont>.csv   (from step 1)
    ${{efo,bto,cl,clo}_owl}                          (downloaded; see README)
    ${target_names}

Outputs:
    ${ontology_collapsed_dir}/<ont>_term_to_annotations_collapsed.csv
    ${ontology_collapsed_dir}/<ont>_collapse_log.csv
"""
from __future__ import annotations

import argparse
import importlib.util
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional, Set

import pandas as pd
import pronto
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from paths import load_paths  # noqa: E402


# ---------------------------------------------------------------------------
# Load is_too_generic from `1-four_ontologies.py` (filename starts with a
# digit, so we can't `import` directly).
# ---------------------------------------------------------------------------

def _load_is_too_generic():
    sibling = Path(__file__).resolve().parent / "1-four_ontologies.py"
    spec = importlib.util.spec_from_file_location("_four_ontologies", sibling)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.is_too_generic


is_too_generic = _load_is_too_generic()


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ONTOLOGIES_DEFAULT = ["EFO", "BTO", "CL", "CLO"]
DEFAULT_RANKS_TO_USE = [1, 2, 3, 4]

DISEASE_KEYWORDS = [
    "disease", "disorder", "syndrome",
    "neoplasm", "tumor",
    "malignancy", "pathology", "abnormality", "malignant",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def normalize_term_id(term_id) -> Optional[str]:
    if term_id is None:
        return None
    term_id = str(term_id).strip()

    if term_id.startswith("http"):
        m = re.search(r"/obo/([A-Z]+)_(\d+)", term_id)
        if m:
            return f"{m.group(1)}:{m.group(2)}"
        m_efo = re.search(r"/efo/([A-Z]+)_(\d+)", term_id)
        if m_efo:
            return f"{m_efo.group(1)}:{m_efo.group(2)}"

    if ":" in term_id:
        return term_id
    if "_" in term_id:
        return term_id.replace("_", ":")
    return term_id


def clean_tissue(name):
    if pd.isna(name):
        return name
    return re.sub(r"\s+\(.*?\)", "", str(name)).strip()


def load_target_names(target_names_file: Path) -> pd.DataFrame:
    target_names = pd.read_csv(
        target_names_file, sep="|", header=None,
        names=["Tissue", "Feature", "ProjectID", "Else"],
    )
    target_names.reset_index(inplace=True)
    for col in target_names.columns:
        if target_names[col].dtype == "object":
            target_names[col] = (
                target_names[col].str.strip().str.replace(r"\s+", " ", regex=True)
            )
    target_names["Tissue_clean"] = target_names["Tissue"].apply(clean_tissue)
    return target_names


def load_ontology(ontology_name: str, owl_path: Path):
    if not owl_path.exists():
        print(f"[WARN] OWL not found for {ontology_name}: {owl_path}")
        return None
    print(f"Loading {ontology_name} ontology from {owl_path}...")
    ont = pronto.Ontology(str(owl_path))
    print(f"  Loaded {len(list(ont.terms()))} terms")
    return ont


def _candidate_keys(term_id: str):
    norm = normalize_term_id(term_id)
    if not norm or ":" not in norm:
        return [term_id]
    prefix, id_num = norm.split(":")
    keys = [
        f"{prefix}_{id_num}",
        norm,
        f"http://purl.obolibrary.org/obo/{prefix}_{id_num}",
    ]
    if prefix == "EFO":
        keys.append(f"http://www.ebi.ac.uk/efo/EFO_{id_num}")
    return keys


def get_term_name(term_id: str, ontology) -> str:
    if ontology is None:
        return term_id
    seen = set()
    for cid in _candidate_keys(term_id):
        if cid in seen:
            continue
        seen.add(cid)
        try:
            t = ontology.get(cid)
            if t and t.name:
                return t.name.strip()
        except Exception:
            continue
    return term_id


def get_all_ancestors(term_id: str, ontology) -> Set[str]:
    out: Set[str] = set()
    if ontology is None or term_id is None:
        return out
    norm = normalize_term_id(term_id)
    try:
        term_obj = None
        for cid in _candidate_keys(norm):
            if cid in ontology:
                term_obj = ontology[cid]
                break
        if term_obj is None:
            return {norm}
        out.add(norm)
        for anc in term_obj.superclasses():
            out.add(normalize_term_id(str(anc.id)))
    except KeyError:
        out.add(norm)
    return out


def is_disease_term_id(term_id: Optional[str]) -> bool:
    if term_id is None:
        return False
    norm = normalize_term_id(term_id)
    return bool(norm) and norm.startswith("MONDO:")


# ---------------------------------------------------------------------------
# Stage 1: matching -> ann_idx → {term_id}
# ---------------------------------------------------------------------------

def create_annotation_to_terms_mapping(
    ontology_name: str,
    target_names: pd.DataFrame,
    matched_df: pd.DataFrame,
    ranks_to_use,
):
    print("\n" + "=" * 60)
    print(f"Create annotation→terms mapping for {ontology_name}")
    print(f"Using ranks: {ranks_to_use}")
    print("=" * 60)

    matched_filtered = matched_df[
        matched_df["rank"].isin(ranks_to_use) & matched_df["term_id"].notna()
    ].copy()
    print(f"After rank filter: {len(matched_filtered)} rows")

    merged = matched_filtered.merge(
        target_names[["index", "Tissue_clean"]],
        left_on="Original Cell Type Name",
        right_on="Tissue_clean",
        how="left",
    )
    print(f"After merge with target_names: {len(merged)} rows")
    print(f"  Rows with annotation index: {merged['index'].notna().sum()}")

    ann_to_terms: Dict[int, Set[str]] = defaultdict(set)
    for _, row in merged.iterrows():
        ann_idx = row["index"]
        term_id = row["term_id"]
        if pd.notna(ann_idx) and pd.notna(term_id):
            ann_to_terms[int(ann_idx)].add(normalize_term_id(term_id))

    print(f"Annotations with ≥1 term: {len(ann_to_terms)} / {len(target_names)}")
    return ann_to_terms


# ---------------------------------------------------------------------------
# Stage 2: ancestor expand + filter diseases + generic terms
# ---------------------------------------------------------------------------

def build_term_universe_with_ancestors(
    ontology_name: str,
    ann_to_terms: Dict[int, Set[str]],
    ontology,
):
    print("\n" + "=" * 60)
    print(f"Build term universe (with ancestors) for {ontology_name}")
    print("=" * 60)

    term_to_ann_indices: Dict[str, Set[int]] = defaultdict(set)
    term_to_ancestors: Dict[str, Set[str]] = defaultdict(set)

    # Iterate sets in sorted order so dict-key insertion order downstream is
    # deterministic across runs (Python set iteration is unstable).
    for ann_idx, term_ids in tqdm(ann_to_terms.items(), desc="Expanding terms"):
        for t_id in sorted(term_ids):
            norm = normalize_term_id(t_id)
            term_to_ann_indices[norm].add(ann_idx)

            ancestors = get_all_ancestors(norm, ontology)
            term_to_ancestors[norm].update(ancestors)

            for anc in sorted(ancestors):
                anc_norm = normalize_term_id(anc)
                if anc_norm == norm:
                    continue
                if is_disease_term_id(anc_norm):
                    continue
                anc_name = get_term_name(anc_norm, ontology).lower()
                if any(k in anc_name for k in DISEASE_KEYWORDS):
                    continue
                term_to_ann_indices[anc_norm].add(ann_idx)

    print(f"\nTotal raw terms (incl. ancestors): {len(term_to_ann_indices)}")

    filtered_term_to_ann: Dict[str, Set[int]] = {}
    removed_generic = 0
    for term_id in sorted(term_to_ann_indices):
        inds = term_to_ann_indices[term_id]
        if not inds:
            continue
        if is_too_generic(get_term_name(term_id, ontology), term_id):
            removed_generic += 1
            continue
        filtered_term_to_ann[term_id] = inds

    print(f"  Terms removed as too generic: {removed_generic}")
    print(f"  Terms kept (post-filter): {len(filtered_term_to_ann)}")

    filtered_term_to_anc = {
        t_id: {normalize_term_id(a) for a in anc_set}
        for t_id, anc_set in term_to_ancestors.items()
        if t_id in filtered_term_to_ann
    }
    return filtered_term_to_ann, filtered_term_to_anc


# ---------------------------------------------------------------------------
# Stage 3: collapse identical-coverage terms (keep leaves)
# ---------------------------------------------------------------------------

def collapse_redundant_terms(term_to_ann_indices, term_to_ancestors, ontology, ontology_name):
    print("\n" + "=" * 60)
    print("Collapse redundant nested terms (same coverage, keep leaves only)")
    print("=" * 60)

    coverage_groups: Dict[frozenset, list] = defaultdict(list)
    for term_id in sorted(term_to_ann_indices):
        coverage_groups[frozenset(term_to_ann_indices[term_id])].append(term_id)

    print(f"Total distinct coverage signatures: {len(coverage_groups)}")

    # Sort groups deterministically: by the lexicographically-smallest term in
    # each group. This makes Coverage_Group_ID assignment stable across runs.
    sorted_groups = sorted(coverage_groups.items(), key=lambda kv: sorted(kv[1])[0])

    kept_terms: Set[str] = set()
    dropped_terms: Set[str] = set()
    collapse_log: list = []

    for group_idx, (_sig, terms) in enumerate(sorted_groups, start=1):
        group_id = f"group_{group_idx}"

        if len(terms) == 1:
            t = terms[0]
            kept_terms.add(t)
            collapse_log.append({
                "Ontology": ontology_name,
                "Coverage_Group_ID": group_id,
                "Term_ID": t,
                "Term_Name": get_term_name(t, ontology) if ontology is not None else None,
                "N_Annotations": len(term_to_ann_indices[t]),
                "Role": "kept_unique",
                "Within_Group_Ancestors": "",
                "Within_Group_Descendants": "",
                "Filtered_Group_Terms": "",
            })
            continue

        group_terms = set(terms)
        non_leaf: Set[str] = set()
        for t in sorted(group_terms):
            for a in term_to_ancestors.get(t, set()):
                if a in group_terms and a != t:
                    non_leaf.add(a)
        leaves = group_terms - non_leaf or group_terms

        kept_terms.update(leaves)
        dropped_terms.update(group_terms - leaves)
        filtered_terms = sorted(group_terms - leaves)

        for t in sorted(group_terms):
            anc_in_group = sorted(
                a for a in term_to_ancestors.get(t, set()) if a in group_terms and a != t
            )
            desc_in_group = sorted(
                d for d in group_terms if t in term_to_ancestors.get(d, set()) and d != t
            )
            collapse_log.append({
                "Ontology": ontology_name,
                "Coverage_Group_ID": group_id,
                "Term_ID": t,
                "Term_Name": get_term_name(t, ontology) if ontology is not None else None,
                "N_Annotations": len(term_to_ann_indices[t]),
                "Role": "kept_leaf" if t in leaves else "dropped_ancestor",
                "Within_Group_Ancestors": ";".join(anc_in_group),
                "Within_Group_Descendants": ";".join(desc_in_group),
                "Filtered_Group_Terms": ";".join(filtered_terms) if t in leaves else "",
            })

    print(f"Total terms before collapsing: {len(term_to_ann_indices)}")
    print(f"  Kept (leaf or unique): {len(kept_terms)}")
    print(f"  Dropped (redundant ancestors): {len(dropped_terms)}")

    collapsed = {t_id: term_to_ann_indices[t_id] for t_id in kept_terms}
    return collapsed, collapse_log


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def save_collapse_log(collapsed_dir: Path, ontology_name: str, collapse_log: list):
    if not collapse_log:
        print(f"No collapse log entries for {ontology_name}")
        return None
    df_log = pd.DataFrame(collapse_log)
    out = collapsed_dir / f"{ontology_name.lower()}_collapse_log.csv"
    df_log.to_csv(out, index=False)
    print(f"Saved collapse log for {ontology_name} to: {out}")
    return df_log


def save_collapsed_mapping(collapsed_dir: Path, ontology_name: str, collapsed, ontology):
    rows = []
    for term_id in sorted(collapsed):
        t_name = get_term_name(term_id, ontology)
        for ann_idx in sorted(collapsed[term_id]):
            rows.append({
                "Ontology": ontology_name,
                "Term_ID": term_id,
                "Term_Name": t_name,
                "Annotation_Index": int(ann_idx),
            })
    df = pd.DataFrame(rows)
    out = collapsed_dir / f"{ontology_name.lower()}_term_to_annotations_collapsed.csv"
    df.to_csv(out, index=False)
    print(f"Saved collapsed mapping for {ontology_name} to: {out}")
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    p.add_argument("--paths_yaml", default=None)
    p.add_argument("--ontologies", nargs="+", default=None,
                   help=f"Ontology acronyms to process. Default: {ONTOLOGIES_DEFAULT}")
    p.add_argument("--ranks_to_use", nargs="+", type=int, default=DEFAULT_RANKS_TO_USE,
                   help=f"Ranks (from step 1) to keep. Default: {DEFAULT_RANKS_TO_USE}")
    return p.parse_args()


def _owl_path_for(paths: dict, ontology: str) -> Path:
    return Path(paths[f"{ontology.lower()}_owl"])


def main() -> None:
    args = parse_args()
    paths = load_paths(args.paths_yaml)
    ontologies = args.ontologies or ONTOLOGIES_DEFAULT

    matched_dir = Path(paths["ontology_matched_dir"])
    collapsed_dir = Path(paths["ontology_collapsed_dir"])
    collapsed_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 80)
    print("STEP 2: ontology ancestor expansion + collapsing")
    print(f"Ontologies: {ontologies}")
    print(f"Ranks to use: {args.ranks_to_use}")
    print("=" * 80 + "\n")

    target_names = load_target_names(Path(paths["target_names"]))

    for ont_name in ontologies:
        print(f"\n\n######## {ont_name} ########")
        ontology = load_ontology(ont_name, _owl_path_for(paths, ont_name))
        if ontology is None:
            continue

        matched_file = matched_dir / f"matched_sup2_{ont_name.lower()}.csv"
        if not matched_file.exists():
            print(f"[WARN] No matched file for {ont_name}: {matched_file} — skip")
            continue
        matched_df = pd.read_csv(matched_file)

        ann_to_terms = create_annotation_to_terms_mapping(
            ont_name, target_names, matched_df, args.ranks_to_use
        )
        term_to_ann_raw, term_to_ancestors = build_term_universe_with_ancestors(
            ont_name, ann_to_terms, ontology
        )
        collapsed, collapse_log = collapse_redundant_terms(
            term_to_ann_indices=term_to_ann_raw,
            term_to_ancestors=term_to_ancestors,
            ontology=ontology,
            ontology_name=ont_name,
        )
        save_collapsed_mapping(collapsed_dir, ont_name, collapsed, ontology)
        save_collapse_log(collapsed_dir, ont_name, collapse_log)


if __name__ == "__main__":
    main()
