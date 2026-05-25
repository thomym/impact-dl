"""GSEA-like enrichment helpers used by 3-enrichment_per_gwas.py.

Provides:
- compute_enrichment_score: one-sided GSEA-like enrichment score for a
  foreground annotation set.
- permutation_test: NES and empirical ES p-value (seeded for reproducible
  output).
- statistical_tests: Mann-Whitney U and Fisher's exact test on the SNP-score
  distribution.
"""
import numpy as np
from scipy.stats import mannwhitneyu, fisher_exact


def compute_enrichment_score(pvalues_df, tissue_indices):
    """One-sided (positive) GSEA-like enrichment score.

    `pvalues_df` is indexed by Annotation_Index and contains a column
    `Empirical_P_Upper`. `tissue_indices` is an array of Annotation_Index
    values marking the foreground set.
    """
    pvalues_df["Is_Tissue"] = pvalues_df.index.isin(tissue_indices)

    p_values_df_sorted = pvalues_df.sort_values("Empirical_P_Upper")
    sorted_pvals = p_values_df_sorted["Empirical_P_Upper"].values
    is_tissue_sorted = p_values_df_sorted["Is_Tissue"].values

    N_tissue = int(np.sum(is_tissue_sorted))
    N_nontissue = len(is_tissue_sorted) - N_tissue

    step_up = 1 / N_tissue if N_tissue > 0 else 0
    step_down = 1 / N_nontissue if N_nontissue > 0 else 0
    ES = np.cumsum(is_tissue_sorted * step_up - (~is_tissue_sorted) * step_down)

    max_ES_index = int(np.argmax(ES))
    max_ES = float(ES[max_ES_index])

    return {
        "ES": ES,
        "max_ES": max_ES,
        "max_ES_index": max_ES_index,
        "sorted_pvals": sorted_pvals,
        "is_tissue_sorted": is_tissue_sorted,
        "N_tissue": N_tissue,
        "N_nontissue": N_nontissue,
    }


def permutation_test(enrichment_data, n_permutations=1000, seed=42):
    """Permutation null for the (positive) ES.

    `seed` controls the per-call RNG so NES and ES_pvalue are reproducible.
    Pass `seed=None` to fall back to fresh randomness each call.
    """
    is_tissue_sorted = enrichment_data["is_tissue_sorted"]
    max_ES = enrichment_data["max_ES"]
    step_up = 1 / enrichment_data["N_tissue"] if enrichment_data["N_tissue"] > 0 else 0
    step_down = 1 / enrichment_data["N_nontissue"] if enrichment_data["N_nontissue"] > 0 else 0

    rng = np.random.default_rng(seed)
    permuted_ES = np.zeros(n_permutations)
    for i in range(n_permutations):
        shuffled_labels = rng.permutation(is_tissue_sorted)
        perm_ES = np.cumsum(shuffled_labels * step_up - (~shuffled_labels) * step_down)
        permuted_ES[i] = np.max(perm_ES)

    ES_pvalue = float(np.mean(permuted_ES >= max_ES))

    std = np.std(permuted_ES)
    NES = float((max_ES - np.mean(permuted_ES)) / std) if std > 0 else 0.0

    return {"ES_pvalue": ES_pvalue, "NES": NES}


def statistical_tests(enrichment_data, pvalues_df):
    """Mann-Whitney U + Fisher's exact on the foreground vs background p-values."""
    is_tissue_sorted = enrichment_data["is_tissue_sorted"]
    sorted_pvals = enrichment_data["sorted_pvals"]

    tissue_pvals = sorted_pvals[is_tissue_sorted]
    non_tissue_pvals = sorted_pvals[~is_tissue_sorted]

    _, mw_pvalue = mannwhitneyu(tissue_pvals, non_tissue_pvals, alternative="less")

    tissue_sig = int(np.sum(tissue_pvals < 0.01))
    tissue_not_sig = int(np.sum(tissue_pvals >= 0.01))
    non_tissue_sig = int(np.sum(non_tissue_pvals < 0.01))
    non_tissue_not_sig = int(np.sum(non_tissue_pvals >= 0.01))

    fisher_table = np.array(
        [[tissue_sig, tissue_not_sig], [non_tissue_sig, non_tissue_not_sig]]
    )
    _, fisher_pvalue = fisher_exact(fisher_table, alternative="greater")

    return {"mw_pvalue": float(mw_pvalue), "fisher_pvalue": float(fisher_pvalue)}
