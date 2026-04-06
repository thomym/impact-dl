import pandas as pd
import numpy as np
import pickle
import re
import requests
import os, time, urllib.parse, requests
import logging
from tqdm import tqdm
import json
from typing import Optional, Dict, List, Any

# ---- LOGGING SETUP ----
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('efo_annotation.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# ---- CONFIG ----
BIOPORTAL_API_KEY = os.getenv("BIOPORTAL_API_KEY", "API_KEY")#your own API key here)
"""
BIOPORTAL_API_KEY = os.getenv("BIOPORTAL_API_KEY")
if not BIOPORTAL_API_KEY:
    raise RuntimeError(
        "BIOPORTAL_API_KEY is not set. "
        "Set it in your environment, e.g. export BIOPORTAL_API_KEY=..."
    )

"""

BASE = "https://data.bioontology.org"
SLEEP = 0.2            # gentle rate limit between Annotator calls
TARGET_COL = "Standardized Celltype Name"  # column in `table` to annotate

# Cache files for persistence
CACHE_FILE = "bioportal_class_cache.json"               # Class details cache
ANNOTATION_CACHE_PREFIX = "annotation_cache"            # Prefix for annotation caches (one per ontology)

# Ontologies to analyze separately
ONTOLOGIES = ["EFO", "CLO", "CL", "BTO"]

# Filtering mode: 
# "strict" = only native terms (EFO_ for EFO, CLO_ for CLO, etc.)
# "permissive" = native + legitimate imports (UBERON, CL) but filter diseases/proteins
FILTERING_MODE = "permissive"  # Change to "strict" if you want native-only

TARGET_NAMES_FILE = "/your/sei/path/sei-framework/model/target.names" #your own path here
SUP2_FIL = "path/to/sup_table_2_sei_article" #your own path here

# ---- Helpers ----
def variants(s: str):
    """Generate minimally different variants:
       raw, _→space, -→space, remove -, remove _, both→spaces, remove both."""
    s = (s or "").strip()
    candidates = {
        s,
        s.replace("_", " "),
        s.replace("-", " "),
        s.replace("-", ""),
        s.replace("_", ""),
        s.replace("_", " ").replace("-", " "),
        s.replace("_", "").replace("-", ""),
    }
    out, seen = [], set()
    for v in (" ".join(x.split()) for x in candidates):
        if v and v not in seen:
            seen.add(v); out.append(v)
    return out

def compact_efo(iri: str) -> Optional[str]:
    """Convert IRI to compact form (e.g., EFO:0000814)"""
    if not iri:
        return iri
    # Handle different ontology formats
    for ont in ["EFO", "CL", "BTO", "CLO", "MONDO", "PR"]:
        if f"{ont}_" in iri:
            return f"{ont}:" + iri.rsplit(f"{ont}_", 1)[-1]
    return iri

def score(query: str, pref: Optional[str], match_type: Optional[str]) -> float:
    q = (query or "").lower()
    p = (pref or "").lower()
    base = 1.0 if q == p else (1.0 if (p and (p in q or q in p)) else 0.0)
    mt = (match_type or "").upper()
    bonus = 0.15 if mt == "PREF" else 0.05 if mt == "SYN" else 0.0
    return base + bonus

def extract_core_celltype(name: str) -> str:
    """
    Extract the core cell type name, removing specific identifiers.
    
    Examples:
        "501-Mel_melanoma" → "melanoma"
        "501-Mel_Melanoma_Cell" → "melanoma cell"
        "A549_lung_cancer" → "lung cancer"
    """
    # Remove leading cell line identifiers (e.g., "501-Mel", "A549")
    name = re.sub(r'^\d+[-_]?\w+[-_]', '', name)
    
    # Handle trailing _Cell
    name = re.sub(r'_[Cc]ell$', ' cell', name)
    
    # Clean up underscores and dashes
    name = name.replace('_', ' ').replace('-', ' ')
    
    # Collapse multiple spaces
    name = ' '.join(name.split())
    
    return name.lower().strip()

def enhanced_variants(s: str) -> List[str]:
    """
    Generate variants including core cell type extraction.
    
    For "501-Mel_melanoma":
    - Original variants: ["501-Mel_melanoma", "501 Mel melanoma", ...]
    - Core type variants: ["melanoma", "melanoma cell", "melanoma cell line"]
    """
    # Get original variants
    original_variants = variants(s)
    
    # Get core cell type
    core = extract_core_celltype(s)
    
    # Generate core variants if different from original
    core_variants = []
    if core and core not in [v.lower() for v in original_variants]:
        core_variants = [
            core,
            f"{core} cell",
            f"{core} cell line",
        ]
    
    # Combine and deduplicate
    all_variants = original_variants + core_variants
    seen = set()
    result = []
    for v in all_variants:
        v_clean = ' '.join(v.split())
        if v_clean and v_clean.lower() not in seen:
            seen.add(v_clean.lower())
            result.append(v_clean)
    
    return result

def is_human_relevant(cls_data: Dict, term_label: str, synonyms: str) -> bool:
    """
    Check if the ontology term is relevant for human cell lines/tissues.
    
    Filters out:
    - Non-human organisms (mouse, rat, etc.)
    - Disease terms (MONDO ontology)
    """
    # Combine all text to search
    all_text = f"{term_label} {synonyms}".lower()
    
    # Exclude non-human organisms
    non_human_keywords = [
        'mouse', 'murine', 'mus musculus',
        'rat', 'rattus',
        'bovine', 'cow',
        'porcine', 'pig', 'swine',
        'canine', 'dog',
        'feline', 'cat',
        'chicken', 'gallus',
        'xenopus', 'frog',
        'zebrafish', 'danio',
        'drosophila', 'fly',
        'c. elegans', 'c elegans', 'worm',
        'yeast', 'saccharomyces',
        'hamster',
    ]
    
    for keyword in non_human_keywords:
        # Create a regex pattern to match the keyword as a whole word
        # re.escape is used to handle keywords with special regex characters (e.g., 'c. elegans')
        pattern = r'\b' + re.escape(keyword) + r'\b' 
        
        if re.search(pattern, all_text):
            logger.info(f"Filtered out (non-human: {keyword}): {term_label}")
            return False
            
    # Check if it's a disease term (from MONDO or disease-related)
    if cls_data:
        iri = cls_data.get('@id', '')
        if 'MONDO_' in iri or 'mondo' in iri.lower():
            logger.info(f"Filtered out (MONDO disease): {term_label}")
            return False
    
    # Additional disease filtering (only for non-cell-line terms)
    if 'cell line' not in all_text and 'cell' not in all_text and 'tissue' not in all_text:
        disease_keywords = ['disease', 'disorder', 'syndrome']
        for keyword in disease_keywords:
            if keyword in all_text:
                logger.info(f"Filtered out (disease term): {term_label}")
                return False
    
    return True

# --- HTTP helpers ---
def _get(url: str, params: Optional[Dict] = None, timeout: int = 30):
    params = dict(params or {})
    params.setdefault("apikey", BIOPORTAL_API_KEY)
    
    logger.debug(f"GET request to: {url}")
    try:
        r = requests.get(url, params=params, timeout=timeout)
        if r.status_code == 429:
            logger.warning("Rate limit hit (429), waiting 2 seconds...")
            time.sleep(2.0)
            r = requests.get(url, params=params, timeout=timeout)
        r.raise_for_status()
        logger.debug(f"Request successful (status {r.status_code})")
        return r
    except requests.exceptions.RequestException as e:
        logger.error(f"Request failed: {e}")
        raise

# Tiny cache for class JSON - now with persistence
_cls_cache = {}

def load_cache():
    """Load cache from disk if it exists."""
    global _cls_cache
    if os.path.exists(CACHE_FILE):
        try:
            with open(CACHE_FILE, 'r') as f:
                _cls_cache = json.load(f)
            logger.info(f"Loaded {len(_cls_cache)} entries from cache file")
        except Exception as e:
            logger.warning(f"Could not load cache file: {e}")
            _cls_cache = {}
    else:
        logger.info("No cache file found, starting fresh")

def save_cache():
    """Save cache to disk."""
    try:
        with open(CACHE_FILE, 'w') as f:
            json.dump(_cls_cache, f, indent=2)
        logger.info(f"Saved {len(_cls_cache)} entries to cache file")
    except Exception as e:
        logger.error(f"Could not save cache file: {e}")

def load_annotation_cache(ontology: str) -> Dict:
    """Load previously completed annotations from disk for specific ontology."""
    cache_file = f"{ANNOTATION_CACHE_PREFIX}_{ontology.lower()}.pkl"
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'rb') as f:
                cache = pickle.load(f)
            logger.info(f"Loaded {len(cache)} pre-annotated cell types from {ontology} cache")
            return cache
        except Exception as e:
            logger.warning(f"Could not load annotation cache for {ontology}: {e}")
            return {}
    else:
        logger.info(f"No annotation cache found for {ontology}, starting fresh")
        return {}

def save_annotation_cache(label_to_cands: Dict, ontology: str):
    """Save completed annotations to disk for resume capability."""
    cache_file = f"{ANNOTATION_CACHE_PREFIX}_{ontology.lower()}.pkl"
    try:
        with open(cache_file, 'wb') as f:
            pickle.dump(label_to_cands, f)
        logger.info(f"Saved {len(label_to_cands)} annotations to {ontology} cache")
    except Exception as e:
        logger.error(f"Could not save annotation cache for {ontology}: {e}")

def fetch_class_from_self(self_url: str):
    """Fetch class JSON using the 'self' URL from BioPortal."""
    if not self_url:
        return None
    if self_url in _cls_cache:
        logger.debug(f"Cache hit for: {self_url}")
        return _cls_cache[self_url]
    
    logger.debug(f"Cache miss, fetching: {self_url}")
    try:
        data = _get(self_url).json()
        _cls_cache[self_url] = data
        return data
    except Exception as e:
        logger.error(f"Failed to fetch class from self URL: {e}")
        return None

def fetch_class_fallbacks(iri: str, ont_acronym: Optional[str]):
    """Try ontology-scoped then global classes endpoints."""
    key_iri = iri or ""
    if key_iri in _cls_cache:
        logger.info(f"Cache hit for IRI: {key_iri}")
        return _cls_cache[key_iri]

    encoded = urllib.parse.quote_plus(key_iri)
    
    # 1) /ontologies/{ACRONYM}/classes/{encoded_IRI}
    if ont_acronym:
        try:
            url = f"{BASE}/ontologies/{ont_acronym}/classes/{encoded}"
            logger.info(f"Trying ontology-scoped fetch: {ont_acronym}")
            data = _get(url).json()
            _cls_cache[key_iri] = data
            return data
        except requests.HTTPError as e:
            logger.debug(f"Ontology-scoped fetch failed: {e}")
    
    # 2) /classes/{encoded_IRI}
    try:
        url = f"{BASE}/classes/{encoded}"
        logger.debug("Trying global classes fetch")
        data = _get(url).json()
        _cls_cache[key_iri] = data
        return data
    except requests.HTTPError as e:
        logger.debug(f"Global fetch failed: {e}")
        return None

def fetch_class(iri: str, self_url: Optional[str], ont_acronym: Optional[str]):
    """Robust class fetch:
       1) use 'self' URL if present; 2) fallback to ontology-scoped; 3) fallback to global."""
    # prefer 'self' URL
    if self_url:
        data = fetch_class_from_self(self_url)
        if data:
            return data
    # fallbacks
    return fetch_class_fallbacks(iri, ont_acronym)
def is_too_generic(term_label: str, term_id: str) -> bool:
    """
    Filter out overly generic terms that don't add information for enrichment.
    """
    if not term_label:
        return True
    
    label_lower = term_label.lower().strip()
    
    # List of overly generic terms
    generic_terms = [
        "cell",
        "cell line",
        "cell type", 
        "cells",
        "cellular_component",
        "tissue",
        "tissues",
        "sample",
        "culture",
        "specimen",
        "material",
        "anatomical structure",
        "anatomical entity",
        "experimental factor",
        "material entity",
        "biological sample",
        "organ",
        "organism part",
        "structure",
        "body part",
        "organism",
        "nucleate cell", 
        "single nucleate cell",
        "material anatomical entity",
        "immortal cell line cell", 
        "processed material", 
        "cell in vitro", 
        "cultured cell", 
        "experimentally modified cell in vitro", 
        "secondary cultured cell", 
        "cell line cell", 
        "animal cell",
        "quality",
        "transparent",
        "continuant",
        "entity",
        "independent continuant", 
        "obsolete_class", 
        "multicellular anatomical structure",
        "material property",
        "disposition", 
        "hapmap cell line",
        "encode cell line",
        "homo sapiens cell line",
        "anatomy basic component",
        "anatomical modifier",
        "geometric modifier"



    ]
    
    # Exact match to generic terms
    if label_lower in generic_terms:
        logger.info(f"Filtered out (too generic): {term_label}")
        return True
    
    # Very short non-specific terms (1-3 characters)
    if len(label_lower) <= 2:
        logger.info(f"Filtered out (too short/generic): {term_label}")
        return True
    
    return False

def annotate_single_ontology(label: str, ontology: str, top_k: int = 5):
    """
    Annotate a single label using one specific ontology.
    
    Args:
        label: Cell type name to annotate
        ontology: Ontology acronym (EFO, CLO, CL, BTO)
        top_k: Number of top results to return
    
    Returns:
        List of candidate dictionaries (deduplicated by term ID)
    """
    logger.info(f"Annotating '{label}' with {ontology}")
    url = f"{BASE}/annotator"
    
    # Use enhanced variants (includes core cell type)
    variant_list = enhanced_variants(label)
    logger.debug(f"  Testing {len(variant_list)} variants: {variant_list}")
    
    # Track unique term IDs to avoid duplicates across queries
    # Key: term_id, Value: best scoring candidate for that term
    seen_term_ids = {}
    
    for i, q in enumerate(variant_list, 1):
        logger.debug(f"  Variant {i}/{len(variant_list)}: '{q}'")
        try:
            r = _get(url, params={
                "text": q,
                "ontologies": ontology,
                "longest_only": "true",
                "exclude_numbers": "true",
            }, timeout=45)
            anns = r.json()
            logger.debug(f"    Found {len(anns)} annotations")
            time.sleep(SLEEP)

            for ann in anns:
                ac = ann.get("annotatedClass", {}) or {}
                iri = ac.get("@id")
                links = ac.get("links", {}) or {}
                ont_acronym = (links.get("ontology") or "").split("/")[-1] or ontology
                self_url = links.get("self")
                mt = (ann.get("annotations") or [{}])[0].get("matchType")

                # SMART FILTERING: Remove garbage imports (diseases, proteins)
                # but KEEP legitimate imports (anatomy, cell types) in permissive mode
                
                # Always filter these (garbage for cell line annotation):
                garbage_prefixes = ['MONDO_', '/MONDO_', 'DOID_', '/DOID_', '/PR_', 'PR_', 'HP_', '/HP_']
                is_garbage = any(prefix in iri for prefix in garbage_prefixes)
                
                if is_garbage:
                    logger.info(f"Filtered out (disease/protein): {iri}")
                    continue
                
                # Strict mode: only keep native terms for each ontology
                if FILTERING_MODE == "strict":
                    is_native = False
                    if ontology == "EFO" and 'EFO_' in iri:
                        is_native = True
                    elif ontology == "CLO" and 'CLO_' in iri:
                        is_native = True
                    elif ontology == "CL" and 'CL_' in iri:
                        is_native = True
                    elif ontology == "BTO" and ('BTO:' in iri or 'BTO_' in iri):
                        is_native = True
                    
                    if not is_native:
                        logger.info(f"  Filtered out (strict mode, not native {ontology}): {iri}")
                        continue
                
                # Permissive mode: keep native + legitimate imports (UBERON, CL)
                # No additional filtering needed - garbage already filtered above

                cls = fetch_class(iri=iri, self_url=self_url, ont_acronym=ont_acronym) or {}
                pref = cls.get("prefLabel")
                syns = cls.get("synonym", [])
                syns_str = "; ".join(syns) if isinstance(syns, list) else syns

                # FILTER 1: Remove overly generic terms (not informative for enrichment)
                if is_too_generic(pref, iri):
                    continue  # ← SKIP THIS TERM!

                # FILTER 2: Check if human-relevant (removes non-human organisms and diseases)
                if not is_human_relevant(cls, pref or "", syns_str or ""):
                    continue
                # Calculate score with boost for cell line matches
                base_score = score(q, pref, mt)
                
                # Boost score if it matches cell line/tissue terms
                pref_lower = (pref or "").lower()
                if 'cell line' in pref_lower:
                    base_score += 0.1
                elif 'cell' in pref_lower or 'tissue' in pref_lower:
                    base_score += 0.05

                term_id = compact_efo(iri)
                
                # Deduplicate: Keep only the best scoring instance of each term ID
                if term_id not in seen_term_ids or base_score > seen_term_ids[term_id]["score"]:
                    seen_term_ids[term_id] = {
                        "standardized": label,
                        "query_used": q,
                        "term_iri": iri,
                        "term_id": term_id,
                        "term_label": pref,
                        "ontology": ont_acronym,
                        "matchType": mt,
                        "synonyms": syns_str,
                        "score": round(base_score, 3),
                    }
        except Exception as e:
            logger.error(f"  Error processing variant '{q}': {e}")
            continue
    
    # Convert dict to list
    unique_cands = list(seen_term_ids.values())
    
    # Sort by score
    unique_cands.sort(key=lambda x: x["score"], reverse=True)
    
    # Add rank
    for i, c in enumerate(unique_cands[:top_k], start=1):
        c["rank"] = i
    
    logger.info(f"  {ontology}: Found {len(unique_cands[:top_k])} unique candidates (score range: {unique_cands[0]['score'] if unique_cands else 'N/A'} - {unique_cands[min(top_k-1, len(unique_cands)-1)]['score'] if unique_cands else 'N/A'})")
    return unique_cands[:top_k]

def build_ontology_review(df: pd.DataFrame, ontology: str, col: str = TARGET_COL, head_n: Optional[int] = None) -> pd.DataFrame:
    """Build annotation review for a single ontology."""
    assert col in df.columns, f"Column '{col}' missing in DataFrame"
    
    logger.info(f"\n{'='*80}")
    logger.info(f"Processing {ontology} Ontology")
    logger.info(f"{'='*80}")
    
    # If head_n is None -> all rows, else first head_n
    series = df[col].astype(str)
    if head_n is not None:
        series = series.head(head_n)
        logger.info(f"Processing first {head_n} rows only")
    else:
        logger.info(f"Processing all {len(series)} rows")
    
    # Check for duplicates to optimize
    unique_labels = series.unique()
    logger.info(f"Found {len(unique_labels)} unique cell types out of {len(series)} total rows")
    
    # Load previously completed annotations (RESUME CAPABILITY!)
    label_to_cands = load_annotation_cache(ontology)
    already_done = set(label_to_cands.keys())
    remaining = [label for label in unique_labels if label not in already_done]
    
    if already_done:
        logger.info(f"RESUMING: {len(already_done)} cell types already annotated")
        logger.info(f"{len(remaining)} cell types remaining to process")
    else:
        logger.info(f"Starting fresh annotation of {len(unique_labels)} unique cell types")
    
    # Only process cell types not in cache
    if remaining:
        logger.info("Starting annotation process...")
        for i, label in enumerate(tqdm(remaining, desc=f"Annotating with {ontology}"), 1):
            cands = annotate_single_ontology(label, ontology)
            label_to_cands[label] = cands
            
            # Save BOTH caches periodically (every 10 unique labels)
            if i % 10 == 0:
                save_cache()  # Class details cache
                save_annotation_cache(label_to_cands, ontology)  # Annotation results cache
        
        # Save both caches after all annotations complete
        save_cache()
        save_annotation_cache(label_to_cands, ontology)
    else:
        logger.info("All cell types already annotated! Using cached results.")
    
    # Now build rows for all original indices
    logger.info("Building result dataframe...")
    rows = []
    for idx, label in tqdm(series.items(), total=len(series), desc=f"Creating {ontology} output rows"):
        cands = label_to_cands.get(label, [])
        if not cands:
            rows.append({
                "row_idx": idx,
                "standardized": label,
                "rank": None,
                "term_id": None,
                "term_label": None,
                "ontology": ontology,
                "matchType": None,
                "score": None,
                "query_used": None,
                "synonyms": None,
            })
        else:
            for c in cands:
                c_row = {
                    "row_idx": idx,
                    "standardized": c.get("standardized", label),
                    "rank": c.get("rank"),
                    "term_id": c.get("term_id"),
                    "term_label": c.get("term_label"),
                    "ontology": c.get("ontology", ontology),
                    "matchType": c.get("matchType"),
                    "score": c.get("score"),
                    "query_used": c.get("query_used"),
                    "synonyms": c.get("synonyms"),
                }
                rows.append(c_row)

    return pd.DataFrame(rows, columns=[
        "row_idx", "standardized", "rank", "term_id", "term_label",
        "ontology", "matchType", "score", "query_used", "synonyms"
    ])

def clean_tissue(name):
    if pd.isna(name):
        return name
    # Remove everything in parentheses and strip spaces
    return re.sub(r"\s+\(.*?\)", "", name).strip()

# ---- MAIN EXECUTION ----
if __name__ == "__main__":
    logger.info("="*80)
    logger.info("Starting Multi-Ontology Annotation Pipeline")
    logger.info(f"Ontologies: {', '.join(ONTOLOGIES)}")
    logger.info(f"Filtering mode: {FILTERING_MODE}")
    logger.info("="*80)
    
    # Load cache
    load_cache()
    
    # File paths
    
    # Load data with error handling
    logger.info("Loading input files...")
    try:
        target_names = pd.read_csv(TARGET_NAMES_FILE, sep='|', header=None, 
                                   names=['Tissue', 'Feature', 'ProjectID', 'Else'])
        logger.info(f"Loaded target_names: {len(target_names)} rows")
    except Exception as e:
        logger.error(f"Failed to load target_names: {e}")
        raise
    
    try:
        sup2 = pd.read_csv(SUP2_FIL)
        logger.info(f"Loaded sup2: {len(sup2)} rows, columns: {list(sup2.columns)}")
    except Exception as e:
        logger.error(f"Failed to load sup2: {e}")
        raise
    
    # Validate column exists
    if TARGET_COL not in sup2.columns:
        logger.error(f"Column '{TARGET_COL}' not found in sup2! Available columns: {list(sup2.columns)}")
        raise ValueError(f"Column '{TARGET_COL}' not found")
    
    # Clean target_names
    logger.info("Cleaning target_names data...")
    target_names.reset_index(inplace=True)
    for col in target_names.columns:
        if target_names[col].dtype == "object":
            target_names[col] = target_names[col].str.strip()
            target_names[col] = target_names[col].str.replace(r"\s+", " ", regex=True)
    
    target_names["Tissue_clean"] = target_names["Tissue"].apply(clean_tissue)
    logger.info("Cleaned target_names")
    
    # Process each ontology separately
    for ontology in ONTOLOGIES:
        logger.info(f"\n\n{'#'*80}")
        logger.info(f"# Processing Ontology: {ontology}")
        logger.info(f"{'#'*80}\n")
        
        # Build review for this ontology
        review = build_ontology_review(sup2, ontology, col=TARGET_COL, head_n=None)
        logger.info(f"Completed {ontology} annotation: {len(review)} result rows")
        
        # Merge with original data
        logger.info(f"Merging {ontology} results with original data...")
        merged = sup2.merge(
            review,
            left_index=True,
            right_on="row_idx",
            how="left"
        ).drop(columns=["row_idx"])
        
        # Save output for this ontology
        output_file = f"ontology_matching_res/matched_sup2_{ontology.lower()}.csv"
        logger.info(f"Saving {ontology} results to {output_file}...")
        merged.to_csv(output_file, index=False)
        logger.info(f"Saved {len(merged)} rows to {output_file}")
        
        # Summary statistics for this ontology
        logger.info(f"\n{'-'*80}")
        logger.info(f"{ontology} SUMMARY STATISTICS")
        logger.info(f"{'-'*80}")
        logger.info(f"Total rows processed: {len(sup2)}")
        logger.info(f"Unique cell types: {sup2[TARGET_COL].nunique()}")
        logger.info(f"Rows with {ontology} matches: {merged['term_id'].notna().sum()}")
        logger.info(f"Rows without matches: {merged['term_id'].isna().sum()}")
        logger.info(f"Output saved to: {output_file}")
        logger.info(f"{'-'*80}\n")
    
    # Final summary
    logger.info("\n" + "="*80)
    logger.info("PIPELINE COMPLETED SUCCESSFULLY!")
    logger.info("="*80)
    logger.info(f"Created {len(ONTOLOGIES)} output files:")
    for ontology in ONTOLOGIES:
        logger.info(f"  - ontology_matching_res/matched_sup2_{ontology.lower()}.csv")
    logger.info(f"\nCache entries: {len(_cls_cache)}")
