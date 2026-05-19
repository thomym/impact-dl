"""Step 1 of ontology_matching/ — sup2 cell-type names -> ontology terms via BioPortal Annotator.

GWAS-independent: run once for the whole pipeline (or skip entirely if the
bundled caches under `${bioportal_cache_dir}/` already cover every cell type
in your sup_table_2).

Inputs:
    ${target_names}  — Sei target.names (used to clean tissue strings)
    ${sup2_file}     — Standardized celltype names (one row per Sei profile)
    ${bioportal_cache_dir}/  — annotation_cache_<ont>.pkl + bioportal_class_cache.json

Outputs:
    ${ontology_matched_dir}/matched_sup2_<ont>.csv   (one per ontology)

BioPortal API:
    Set BIOPORTAL_API_KEY in your environment. The key is only needed when the
    cache misses; with the bundled caches a fresh run hits no API at all.
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import pickle
import re
import sys
import time
import urllib.parse
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
import requests
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from paths import load_paths  # noqa: E402


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BASE = "https://data.bioontology.org"
SLEEP = 0.2
TARGET_COL = "Standardized Celltype Name"
ONTOLOGIES_DEFAULT = ["EFO", "CLO", "CL", "BTO"]
FILTERING_MODE = "permissive"   # "permissive" (native + UBERON/CL imports) | "strict"


# ---------------------------------------------------------------------------
# Variant generation
# ---------------------------------------------------------------------------

def variants(s: str) -> List[str]:
    """Minor textual variants: raw, _→space, -→space, drops, etc."""
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
            seen.add(v)
            out.append(v)
    return out


def extract_core_celltype(name: str) -> str:
    """Strip leading cell-line tokens ("501-Mel_melanoma" -> "melanoma")."""
    name = re.sub(r"^\d+[-_]?\w+[-_]", "", name)
    name = re.sub(r"_[Cc]ell$", " cell", name)
    name = name.replace("_", " ").replace("-", " ")
    return " ".join(name.split()).lower().strip()


def enhanced_variants(s: str) -> List[str]:
    original = variants(s)
    core = extract_core_celltype(s)
    extras: List[str] = []
    if core and core not in {v.lower() for v in original}:
        extras = [core, f"{core} cell", f"{core} cell line"]
    seen: set[str] = set()
    out: List[str] = []
    for v in original + extras:
        clean = " ".join(v.split())
        if clean and clean.lower() not in seen:
            seen.add(clean.lower())
            out.append(clean)
    return out


def compact_id(iri: str) -> Optional[str]:
    if not iri:
        return iri
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


# ---------------------------------------------------------------------------
# Term filtering
# ---------------------------------------------------------------------------

_NON_HUMAN = [
    "mouse", "murine", "mus musculus",
    "rat", "rattus",
    "bovine", "cow",
    "porcine", "pig", "swine",
    "canine", "dog",
    "feline", "cat",
    "chicken", "gallus",
    "xenopus", "frog",
    "zebrafish", "danio",
    "drosophila", "fly",
    "c. elegans", "c elegans", "worm",
    "yeast", "saccharomyces",
    "hamster",
]


def is_human_relevant(cls_data: Dict, term_label: str, synonyms: str) -> bool:
    all_text = f"{term_label} {synonyms}".lower()
    for keyword in _NON_HUMAN:
        if re.search(r"\b" + re.escape(keyword) + r"\b", all_text):
            logger.info("Filtered out (non-human: %s): %s", keyword, term_label)
            return False
    if cls_data:
        iri = cls_data.get("@id", "")
        if "MONDO_" in iri or "mondo" in iri.lower():
            logger.info("Filtered out (MONDO disease): %s", term_label)
            return False
    if "cell line" not in all_text and "cell" not in all_text and "tissue" not in all_text:
        for keyword in ("disease", "disorder", "syndrome"):
            if keyword in all_text:
                logger.info("Filtered out (disease term): %s", term_label)
                return False
    return True


_GENERIC_TERMS = {
    "cell", "cell line", "cell type", "cells",
    "cellular_component", "tissue", "tissues", "sample", "culture",
    "specimen", "material", "anatomical structure", "anatomical entity",
    "experimental factor", "material entity", "biological sample",
    "organ", "organism part", "structure", "body part", "organism",
    "nucleate cell", "single nucleate cell", "material anatomical entity",
    "immortal cell line cell", "processed material", "cell in vitro",
    "cultured cell", "experimentally modified cell in vitro",
    "secondary cultured cell", "cell line cell", "animal cell",
    "quality", "transparent", "continuant", "entity",
    "independent continuant", "obsolete_class",
    "multicellular anatomical structure", "material property",
    "disposition", "hapmap cell line", "encode cell line",
    "homo sapiens cell line", "anatomy basic component",
    "anatomical modifier", "geometric modifier",
}


def is_too_generic(term_label: Optional[str], _term_id: str) -> bool:
    if not term_label:
        return True
    label = term_label.lower().strip()
    if label in _GENERIC_TERMS:
        logger.info("Filtered out (too generic): %s", term_label)
        return True
    if len(label) <= 2:
        logger.info("Filtered out (too short/generic): %s", term_label)
        return True
    return False


# ---------------------------------------------------------------------------
# BioPortal HTTP layer + caches
# ---------------------------------------------------------------------------

class BioPortalClient:
    """Per-run client owning the API key + class-detail cache."""

    def __init__(self, api_key: Optional[str], cache_dir: Path):
        self.api_key = api_key
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.class_cache_file = cache_dir / "bioportal_class_cache.json"
        self._cls_cache: Dict[str, Any] = self._load_class_cache()

    # ----- class-detail cache -----
    def _load_class_cache(self) -> Dict[str, Any]:
        if self.class_cache_file.exists():
            try:
                with open(self.class_cache_file) as f:
                    cache = json.load(f)
                logger.info("Loaded %d entries from class cache", len(cache))
                return cache
            except Exception as e:
                logger.warning("Could not load class cache: %s", e)
        else:
            logger.info("No class cache found; starting fresh")
        return {}

    def save_class_cache(self) -> None:
        try:
            with open(self.class_cache_file, "w") as f:
                json.dump(self._cls_cache, f, indent=2)
            logger.info("Saved %d entries to class cache", len(self._cls_cache))
        except Exception as e:
            logger.error("Could not save class cache: %s", e)

    # ----- per-ontology annotation cache -----
    def annotation_cache_path(self, ontology: str) -> Path:
        return self.cache_dir / f"annotation_cache_{ontology.lower()}.pkl"

    def load_annotation_cache(self, ontology: str) -> Dict[str, List[Dict]]:
        path = self.annotation_cache_path(ontology)
        if not path.exists():
            logger.info("No annotation cache for %s; starting fresh", ontology)
            return {}
        try:
            with open(path, "rb") as f:
                cache = pickle.load(f)
            logger.info("Loaded %d cached annotations for %s", len(cache), ontology)
            return cache
        except Exception as e:
            logger.warning("Could not load annotation cache for %s: %s", ontology, e)
            return {}

    def save_annotation_cache(self, ontology: str, label_to_cands: Dict[str, List[Dict]]) -> None:
        path = self.annotation_cache_path(ontology)
        try:
            with open(path, "wb") as f:
                pickle.dump(label_to_cands, f)
            logger.info("Saved %d annotations to %s cache", len(label_to_cands), ontology)
        except Exception as e:
            logger.error("Could not save annotation cache for %s: %s", ontology, e)

    # ----- HTTP -----
    def _ensure_api_key(self) -> str:
        if not self.api_key:
            raise RuntimeError(
                "BIOPORTAL_API_KEY is not set. Either set it in your environment "
                "(`export BIOPORTAL_API_KEY=...`) or use the bundled caches that "
                "cover all cell types in sup_table_2."
            )
        return self.api_key

    def _get(self, url: str, params: Optional[Dict] = None, timeout: int = 30):
        params = dict(params or {})
        params.setdefault("apikey", self._ensure_api_key())
        r = requests.get(url, params=params, timeout=timeout)
        if r.status_code == 429:
            logger.warning("Rate limit hit (429), backing off 2s...")
            time.sleep(2.0)
            r = requests.get(url, params=params, timeout=timeout)
        r.raise_for_status()
        return r

    # ----- class fetch -----
    def fetch_class(self, iri: str, self_url: Optional[str], ont_acronym: Optional[str]):
        if self_url:
            if self_url in self._cls_cache:
                return self._cls_cache[self_url]
            try:
                data = self._get(self_url).json()
                self._cls_cache[self_url] = data
                return data
            except Exception as e:
                logger.debug("self-URL fetch failed: %s", e)

        key_iri = iri or ""
        if key_iri in self._cls_cache:
            return self._cls_cache[key_iri]

        encoded = urllib.parse.quote_plus(key_iri)
        for url in filter(None, [
            f"{BASE}/ontologies/{ont_acronym}/classes/{encoded}" if ont_acronym else None,
            f"{BASE}/classes/{encoded}",
        ]):
            try:
                data = self._get(url).json()
                self._cls_cache[key_iri] = data
                return data
            except requests.HTTPError as e:
                logger.debug("class fetch %s failed: %s", url, e)
        return None

    def annotator(self, text: str, ontology: str):
        return self._get(
            f"{BASE}/annotator",
            params={
                "text": text,
                "ontologies": ontology,
                "longest_only": "true",
                "exclude_numbers": "true",
            },
            timeout=45,
        )


# ---------------------------------------------------------------------------
# Per-label annotation
# ---------------------------------------------------------------------------

def annotate_single_ontology(
    client: BioPortalClient,
    label: str,
    ontology: str,
    top_k: int = 5,
) -> List[Dict]:
    """Annotate one cell-type label against one ontology; returns up to top_k candidates."""
    logger.info("Annotating '%s' with %s", label, ontology)
    variant_list = enhanced_variants(label)
    seen_term_ids: Dict[str, Dict] = {}

    garbage_prefixes = ("MONDO_", "/MONDO_", "DOID_", "/DOID_", "/PR_", "PR_", "HP_", "/HP_")

    for q in variant_list:
        try:
            anns = client.annotator(q, ontology).json()
            time.sleep(SLEEP)

            for ann in anns:
                ac = ann.get("annotatedClass", {}) or {}
                iri = ac.get("@id") or ""
                links = ac.get("links", {}) or {}
                ont_acronym = (links.get("ontology") or "").split("/")[-1] or ontology
                self_url = links.get("self")
                mt = (ann.get("annotations") or [{}])[0].get("matchType")

                if any(prefix in iri for prefix in garbage_prefixes):
                    logger.info("Filtered out (disease/protein): %s", iri)
                    continue

                if FILTERING_MODE == "strict":
                    native_marker = {"EFO": "EFO_", "CLO": "CLO_", "CL": "CL_", "BTO": "BTO_"}.get(ontology)
                    if native_marker and native_marker not in iri and not (ontology == "BTO" and "BTO:" in iri):
                        logger.info("Filtered out (strict mode, not native %s): %s", ontology, iri)
                        continue

                cls = client.fetch_class(iri=iri, self_url=self_url, ont_acronym=ont_acronym) or {}
                pref = cls.get("prefLabel")
                syns = cls.get("synonym", [])
                syns_str = "; ".join(syns) if isinstance(syns, list) else syns

                if is_too_generic(pref, iri):
                    continue
                if not is_human_relevant(cls, pref or "", syns_str or ""):
                    continue

                base_score = score(q, pref, mt)
                pref_lower = (pref or "").lower()
                if "cell line" in pref_lower:
                    base_score += 0.1
                elif "cell" in pref_lower or "tissue" in pref_lower:
                    base_score += 0.05

                term_id = compact_id(iri)
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
            logger.error("  Error processing variant '%s': %s", q, e)
            continue

    cands = sorted(seen_term_ids.values(), key=lambda x: x["score"], reverse=True)
    for i, c in enumerate(cands[:top_k], start=1):
        c["rank"] = i
    return cands[:top_k]


# ---------------------------------------------------------------------------
# Per-ontology review + I/O
# ---------------------------------------------------------------------------

def build_ontology_review(
    client: BioPortalClient,
    df: pd.DataFrame,
    ontology: str,
    col: str = TARGET_COL,
    head_n: Optional[int] = None,
) -> pd.DataFrame:
    assert col in df.columns, f"Column '{col}' missing in sup_table_2"

    logger.info("=" * 80)
    logger.info("Processing %s ontology", ontology)
    logger.info("=" * 80)

    series = df[col].astype(str)
    if head_n is not None:
        series = series.head(head_n)
    unique_labels = list(series.unique())
    logger.info("%d unique labels (of %d total rows)", len(unique_labels), len(series))

    label_to_cands = client.load_annotation_cache(ontology)
    remaining = [lbl for lbl in unique_labels if lbl not in label_to_cands]
    logger.info("Resuming: %d cached, %d to process", len(label_to_cands), len(remaining))

    for i, label in enumerate(tqdm(remaining, desc=f"Annotating with {ontology}"), 1):
        label_to_cands[label] = annotate_single_ontology(client, label, ontology)
        if i % 10 == 0:
            client.save_class_cache()
            client.save_annotation_cache(ontology, label_to_cands)
    client.save_class_cache()
    client.save_annotation_cache(ontology, label_to_cands)

    logger.info("Building result dataframe...")
    rows: List[Dict] = []
    empty_row_template = {
        "rank": None, "term_id": None, "term_label": None, "ontology": ontology,
        "matchType": None, "score": None, "query_used": None, "synonyms": None,
    }
    for idx, label in series.items():
        cands = label_to_cands.get(label, [])
        if not cands:
            rows.append({"row_idx": idx, "standardized": label, **empty_row_template})
        else:
            for c in cands:
                rows.append({
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
                })

    return pd.DataFrame(rows, columns=[
        "row_idx", "standardized", "rank", "term_id", "term_label",
        "ontology", "matchType", "score", "query_used", "synonyms",
    ])


def clean_tissue(name: Optional[str]):
    if pd.isna(name):
        return name
    return re.sub(r"\s+\(.*?\)", "", str(name)).strip()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    p.add_argument("--paths_yaml", default=None)
    p.add_argument("--ontologies", nargs="+", default=None,
                   help=f"Ontology acronyms to process. Default: {ONTOLOGIES_DEFAULT}")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    paths = load_paths(args.paths_yaml)
    ontologies = args.ontologies or ONTOLOGIES_DEFAULT

    cache_dir = Path(paths["bioportal_cache_dir"])
    matched_dir = Path(paths["ontology_matched_dir"])
    matched_dir.mkdir(parents=True, exist_ok=True)

    target_names_file = Path(paths["target_names"])
    sup2_file = Path(paths["sup2_file"])

    logger.info("4a config:")
    logger.info("  target_names:           %s", target_names_file)
    logger.info("  sup2_file:              %s", sup2_file)
    logger.info("  bioportal_cache_dir:    %s", cache_dir)
    logger.info("  ontology_matched_dir:   %s", matched_dir)
    logger.info("  ontologies:             %s", ontologies)

    client = BioPortalClient(api_key=os.getenv("BIOPORTAL_API_KEY"), cache_dir=cache_dir)

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

    sup2 = pd.read_csv(sup2_file)
    if TARGET_COL not in sup2.columns:
        raise ValueError(
            f"Column '{TARGET_COL}' not found in sup2 ({sup2_file}). "
            f"Available: {list(sup2.columns)}"
        )

    for ontology in ontologies:
        logger.info("\n%s\n# Processing ontology: %s\n%s\n", "#" * 80, ontology, "#" * 80)

        review = build_ontology_review(client, sup2, ontology, col=TARGET_COL)
        merged = sup2.merge(review, left_index=True, right_on="row_idx", how="left").drop(
            columns=["row_idx"]
        )

        out_file = matched_dir / f"matched_sup2_{ontology.lower()}.csv"
        merged.to_csv(out_file, index=False)
        logger.info("Saved %s -> %s (%d rows; %d with matches)",
                    ontology, out_file, len(merged), int(merged["term_id"].notna().sum()))

    logger.info("Pipeline complete. Cached %d class entries.", len(client._cls_cache))


if __name__ == "__main__":
    main()
