import pandas as pd
import numpy as np
from collections import defaultdict
from pathlib import Path
import pronto
import re
from tqdm import tqdm

from four_ontologies import is_too_generic

# ========================================================================
# CONFIG 
# ========================================================================

BASE_DIR = Path("" )  # adjust as needed
CODE_DIR = BASE_DIR / "enrichment_pipeline_code"

TARGET_NAMES_FILE = "/home/gaga/thomym/sei/sei-framework/model/target.names"
#SUP2_FILE = CODE_DIR / "ontology_matching_res/sup_table_2_sei_article.csv"

MATCHED_FILES = {
    'EFO': CODE_DIR / "ontology_matching_res/matched_sup2_efo.csv",
    'BTO': CODE_DIR / "ontology_matching_res/matched_sup2_bto.csv",
    'CL':  CODE_DIR / "ontology_matching_res/matched_sup2_cl.csv",
    'CLO': CODE_DIR / "ontology_matching_res/matched_sup2_clo.csv",
}

ONTOLOGY_FILES = {
    'EFO': BASE_DIR / "efo_ontology/efo.owl",
    'BTO': BASE_DIR / "brenda_ontology/bto.owl",
    'CL':  BASE_DIR / "cl_ontology/cl.owl",
    'CLO': BASE_DIR / "clo_ontology/clo.owl",
}

MAPPING_DIR = CODE_DIR / "ontology_annotation_mappings"
MAPPING_DIR.mkdir(parents=True, exist_ok=True)

# which ranks we keep from the matching output
DEFAULT_RANKS_TO_USE = [1, 2, 3, 4]

# ========================================================================
# BASIC HELPERS
# ========================================================================

def normalize_term_id(term_id: str):
    """Converts mixed IDs (URL, EFO_, UBERON_) to canonical PREFIX:ID."""
    if term_id is None:
        return None
    term_id = str(term_id).strip()

    # URLs like http://purl.obolibrary.org/obo/UBERON_0002530
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
    """Remove parentheses content from tissue names."""
    if pd.isna(name):
        return name
    return re.sub(r"\s+\(.*?\)", "", str(name)).strip()


def load_target_names():
    target_names = pd.read_csv(
        TARGET_NAMES_FILE,
        sep="|",
        header=None,
        names=["Tissue", "Feature", "ProjectID", "Else"],
    )
    target_names.reset_index(inplace=True)   # index column is annotation index

    for col in target_names.columns:
        if target_names[col].dtype == "object":
            target_names[col] = (
                target_names[col]
                .str.strip()
                .str.replace(r"\s+", " ", regex=True)
            )

    target_names["Tissue_clean"] = target_names["Tissue"].apply(clean_tissue)
    return target_names


def load_ontology(ontology_name: str):
    obo_file = ONTOLOGY_FILES.get(ontology_name)
    if not obo_file or not obo_file.exists():
        print(f"[WARN] Ontology file not found for {ontology_name}: {obo_file}")
        return None
    print(f"Loading {ontology_name} ontology from {obo_file}...")
    ont = pronto.Ontology(str(obo_file))
    print(f"  Loaded {len(list(ont.terms()))} terms")
    return ont


def get_term_name(term_id: str, ontology):
    """Try multiple key formats until we get a proper .name from pronto."""
    if ontology is None:
        return term_id

    normalized_id = normalize_term_id(term_id)
    if not normalized_id or ":" not in normalized_id:
        return term_id

    prefix, id_num = normalized_id.split(":")

    candidates = [
        f"{prefix}_{id_num}",  # OBO-style key
        normalized_id,         # PREFIX:ID
        f"http://purl.obolibrary.org/obo/{prefix}_{id_num}",
    ]
    if prefix == "EFO":
        candidates.append(f"http://www.ebi.ac.uk/efo/EFO_{id_num}")

    seen = []
    for cid in candidates:
        if cid in seen:
            continue
        seen.append(cid)
        try:
            t = ontology.get(cid)
            if t and t.name:
                return t.name.strip()
        except Exception:
            continue

    return term_id


def get_all_ancestors(term_id: str, ontology):
    """
    Return set of ancestor IDs (including the term itself) using pronto.
    All IDs are normalized.
    """
    out = set()
    if ontology is None or term_id is None:
        return out

    norm_id = normalize_term_id(term_id)
    try:
        # try several keys
        candidates = [norm_id.replace(":", "_"), norm_id,
                      f"http://purl.obolibrary.org/obo/{norm_id.replace(':', '_')}"]
        if norm_id.startswith("EFO:"):
            candidates.append(f"http://www.ebi.ac.uk/efo/{norm_id.replace(':', '_')}")

        term_obj = None
        for cid in candidates:
            if cid in ontology:
                term_obj = ontology[cid]
                break
        if term_obj is None:
            return {norm_id}

        out.add(norm_id)
        for anc in term_obj.superclasses():
            out.add(normalize_term_id(str(anc.id)))
    except KeyError:
        out.add(norm_id)
    return out



def is_disease_term_id(term_id: str):
    """Quick disease ID flag: MONDO or EFO disease branch (heuristic)."""
    if term_id is None:
        return False
    term_id = normalize_term_id(term_id)
    if term_id.startswith("MONDO:"):
        return True
    # for EFO we rely mainly on the name + extra filtering later
    return False

# ========================================================================
# 1. MATCHING: annotation_index -> list of matched term_ids
# ========================================================================

def create_annotation_to_terms_mapping(
    ontology_name: str,
    target_names: pd.DataFrame,
    matched_df: pd.DataFrame,
    ranks_to_use=None,
):
    """
    As in your original create_annotation_to_terms_mapping, but GWAS-agnostic.
    Returns: dict annotation_index -> set(term_id)
    """
    if ranks_to_use is None:
        ranks_to_use = matched_df["rank"].dropna().unique()

    print(f"\n{'='*60}")
    print(f"Create annotation→terms mapping for {ontology_name}")
    print(f"Using ranks: {ranks_to_use}")
    print(f"{'='*60}")

    matched_filtered = matched_df[
        matched_df["rank"].isin(ranks_to_use) & matched_df["term_id"].notna()
    ].copy()

    print(f"After rank filter: {len(matched_filtered)} rows")

    # Here we assume matched_df already has "Original Cell Type Name" and
    # "Standardized Celltype Name" from your sup2 pipeline.
    merged = matched_filtered.merge(
        target_names[["index", "Tissue_clean"]],
        left_on="Original Cell Type Name",
        right_on="Tissue_clean",
        how="left",
    )

    print(f"After merge with target_names: {len(merged)} rows")
    print(f"  Rows with annotation index: {merged['index'].notna().sum()}")

    ann_to_terms = defaultdict(set)
    for _, row in merged.iterrows():
        ann_idx = row["index"]
        term_id = row["term_id"]
        if pd.notna(ann_idx) and pd.notna(term_id):
            ann_to_terms[int(ann_idx)].add(normalize_term_id(term_id))

    print(f"Annotations with ≥1 term: {len(ann_to_terms)} / {len(target_names)}")
    return ann_to_terms

# ========================================================================
# 2. EXPANSION: add ancestors + filter diseases + generic terms
# ========================================================================

def build_term_universe_with_ancestors(
    ontology_name: str,
    ann_to_terms: dict,
    ontology,
):
    """
    - For each annotation, start from its direct matched terms.
    - Add ancestors (excluding diseases / MONDO).
    - Track for each term: which annotations it appears in.
    - Filter out 'too generic' terms with is_too_generic.
    Returns:
        term_to_ann_indices: dict term_id -> sorted list of ann_idx
        term_to_ancestors:   dict term_id -> set(ancestor_ids)
    """
    print(f"\n{'='*60}")
    print(f"Build term universe (with ancestors) for {ontology_name}")
    print(f"{'='*60}")

    term_to_ann_indices = defaultdict(set)
    term_to_ancestors = defaultdict(set)

    disease_keywords = [
        "disease", "disorder", "syndrome", 
        "neoplasm", "tumor",
        "malignancy", "pathology", "abnormality", "malignant",
    ]

    # 1. Expand matched terms with ancestors (excluding disease-y ones)
    for ann_idx, term_ids in tqdm(ann_to_terms.items(), desc="Expanding terms"):
        for t_id in term_ids:
            # always keep the matched term
            norm = normalize_term_id(t_id)
            term_to_ann_indices[norm].add(ann_idx)

            # get ancestors
            ancestors = get_all_ancestors(norm, ontology)
            term_to_ancestors[norm].update(ancestors)

            for anc in ancestors:
                anc_norm = normalize_term_id(anc)
                if anc_norm == norm:
                    continue
                if is_disease_term_id(anc_norm):
                    continue

                anc_name = get_term_name(anc_norm, ontology).lower()
                if any(k in anc_name for k in disease_keywords):
                    # skip disease-y ancestors
                    continue

                # keep this ancestor as a term to test
                term_to_ann_indices[anc_norm].add(ann_idx)

    print(f"\nTotal raw terms (incl. ancestors): {len(term_to_ann_indices)}")

    # 2. Generic-term re-filtering using is_too_generic
    filtered_term_to_ann = {}
    removed_generic = 0
    for term_id, inds in term_to_ann_indices.items():
        if not inds:
            continue
        t_name = get_term_name(term_id, ontology)
        if is_too_generic(t_name, term_id):
            removed_generic += 1
            continue
        filtered_term_to_ann[term_id] = inds

    print(f"  Terms removed as too generic: {removed_generic}")
    print(f"  Terms kept (post-filter): {len(filtered_term_to_ann)}")

    # 3. Keep ancestors only for the kept terms
    filtered_term_to_anc = {}
    for t_id, anc_set in term_to_ancestors.items():
        if t_id not in filtered_term_to_ann:
            continue
        filtered_term_to_anc[t_id] = {normalize_term_id(a) for a in anc_set}

    return filtered_term_to_ann, filtered_term_to_anc

# ========================================================================
# 3. COLLAPSING NESTED TERMS WITH IDENTICAL COVERAGE
# ========================================================================

def collapse_redundant_terms(
    term_to_ann_indices,
    term_to_ancestors,
    ontology=None,
    ontology_name=None,
):
    """
    term_to_ann_indices: dict term_id -> set of annotation indices
    term_to_ancestors:   dict term_id -> set of ancestor term_ids (normalized)

    For all terms that share the same exact annotation set:
      - If some are ancestors of others, keep only the LEAF terms
        (the most specific).

    Returns
    -------
    collapsed : dict
        term_id -> set(annotation_indices) for kept terms only
    collapse_log : list of dict
        Per-term log with group, role (kept/dropped), etc.
    """
    print(f"\n{'='*60}")
    print("Collapse redundant nested terms (same coverage, keep leaves only)")
    print(f"{'='*60}")

    # 1. Group terms by their coverage signature (frozenset of annotation indices)
    coverage_groups = defaultdict(list)
    for term_id, inds in term_to_ann_indices.items():
        sig = frozenset(inds)
        coverage_groups[sig].append(term_id)

    print(f"Total distinct coverage signatures: {len(coverage_groups)}")

    kept_terms = set()
    dropped_terms = set()
    collapse_log = []

    # 2. Within each coverage group, drop ancestors, keep leaves
    for group_idx, (sig, terms) in enumerate(coverage_groups.items(), start=1):
        group_id = f"group_{group_idx}"

        # ------------------------------------------------------------------
        # Case 1: only one term has this coverage → keep it as-is
        # ------------------------------------------------------------------
        if len(terms) == 1:
            t = terms[0]
            kept_terms.add(t)

            collapse_log.append({
                "Ontology": ontology_name,
                "Coverage_Group_ID": group_id,
                "Term_ID": t,
                "Term_Name": get_term_name(t, ontology) if ontology is not None else None,
                "N_Annotations": len(term_to_ann_indices[t]),
                "Role": "kept_unique",      # unique coverage signature
                "Within_Group_Ancestors": "",
                "Within_Group_Descendants": "",
                "Filtered_Group_Terms": "",  # nothing to filter within this group
            })
            continue

        # ------------------------------------------------------------------
        # Case 2: multiple terms with identical coverage → collapse ancestors
        # ------------------------------------------------------------------
        group_terms = set(terms)
        non_leaf = set()

        # Mark all terms that are ancestors of other terms in *this* group
        for t in group_terms:
            anc_set = term_to_ancestors.get(t, set())
            for a in anc_set:
                if a in group_terms and a != t:
                    non_leaf.add(a)

        # Leaves = those that are not ancestors of any other in this group
        leaves = group_terms - non_leaf
        if not leaves:
            # degenerate case: everything looks like an ancestor of something
            # → keep all of them
            leaves = group_terms

        # Bookkeeping
        kept_terms.update(leaves)
        dropped_terms.update(group_terms - leaves)

        # Terms that we filtered out in this group (for logging)
        filtered_terms = sorted(group_terms - leaves)

        # --- populate detailed log for this coverage group ---
        for t in group_terms:
            anc_in_group = sorted(
                a for a in term_to_ancestors.get(t, set())
                if a in group_terms and a != t
            )
            desc_in_group = sorted(
                d for d in group_terms
                if t in term_to_ancestors.get(d, set()) and d != t
            )

            # NEW COLUMN: terms filtered *within this group* —
            # only meaningful for kept leaf terms
            if t in leaves:
                filtered_list = filtered_terms
            else:
                filtered_list = []

            collapse_log.append({
                "Ontology": ontology_name,
                "Coverage_Group_ID": group_id,
                "Term_ID": t,
                "Term_Name": get_term_name(t, ontology) if ontology is not None else None,
                "N_Annotations": len(term_to_ann_indices[t]),
                "Role": "kept_leaf" if t in leaves else "dropped_ancestor",
                "Within_Group_Ancestors": ";".join(anc_in_group),
                "Within_Group_Descendants": ";".join(desc_in_group),
                "Filtered_Group_Terms": ";".join(filtered_list),
            })

    print(f"Total terms before collapsing: {len(term_to_ann_indices)}")
    print(f"  Kept (leaf or unique): {len(kept_terms)}")
    print(f"  Dropped (redundant ancestors): {len(dropped_terms)}")

    collapsed = {
        t_id: term_to_ann_indices[t_id]
        for t_id in kept_terms
    }
    return collapsed, collapse_log

def save_collapse_log(ontology_name: str, collapse_log: list):
    """
    Save a per-term log of what happened during collapsing:
    kept vs dropped, which group, etc.
    """
    if not collapse_log:
        print(f"No collapse log entries for {ontology_name}")
        return None

    df_log = pd.DataFrame(collapse_log)
    out_file = MAPPING_DIR / f"{ontology_name.lower()}_collapse_log.csv"
    df_log.to_csv(out_file, index=False)
    print(f"Saved collapse log for {ontology_name} to: {out_file}")
    return df_log

def save_collapsed_mapping(
    ontology_name: str,
    collapsed_term_to_ann_indices: dict,
    ontology,
):
    """
    Flatten to rows:
        Ontology, Term_ID, Term_Name, Annotation_Index
    and save as one CSV per ontology.
    """
    rows = []
    for term_id, inds in collapsed_term_to_ann_indices.items():
        t_name = get_term_name(term_id, ontology)
        for ann_idx in sorted(inds):
            rows.append(
                {
                    "Ontology": ontology_name,
                    "Term_ID": term_id,
                    "Term_Name": t_name,
                    "Annotation_Index": int(ann_idx),
                }
            )

    df = pd.DataFrame(rows)
    out_file = MAPPING_DIR / f"{ontology_name.lower()}_term_to_annotations_collapsed.csv"
    df.to_csv(out_file, index=False)
    print(f"Saved collapsed mapping for {ontology_name} to: {out_file}")
    return df

# ========================================================================
# MAIN
# ========================================================================


def main(ontologies=None, ranks_to_use=None):
    if ontologies is None:
        ontologies = list(MATCHED_FILES.keys())
    if ranks_to_use is None:
        ranks_to_use = DEFAULT_RANKS_TO_USE

    print(f"\n{'='*80}")
    print("STEP 1: Ontology matching + ancestor expansion + collapsing")
    print(f"Ontologies: {ontologies}")
    print(f"Ranks to use: {ranks_to_use}")
    print(f"{'='*80}\n")

    target_names = load_target_names()

    for ont_name in ontologies:
        print(f"\n\n######## {ont_name} ########")
        ontology = load_ontology(ont_name)
        if ontology is None:
            continue

        matched_df = pd.read_csv(MATCHED_FILES[ont_name])
        ann_to_terms = create_annotation_to_terms_mapping(
            ont_name, target_names, matched_df, ranks_to_use
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

        save_collapsed_mapping(ont_name, collapsed, ontology)
        save_collapse_log(ont_name, collapse_log)


if __name__ == "__main__":
    main()
