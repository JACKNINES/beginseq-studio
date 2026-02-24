"""
classification.py — Expression-based sample classification.

Provides lightweight CPM normalization and expression-based sample
classification to create sub-conditions (e.g., tumor_TNBC vs
tumor_nonTNBC) BEFORE running DESeq2.

Workflow
--------
1. CPM normalization:  cpm = counts / lib_size * 1e6
2. Log-transform:      log_cpm = log2(cpm + 1)
3. Extract marker genes from log_cpm matrix.
4. Classify each sample by expression thresholds.
5. Merge classification into metadata as compound condition.

This is a lightweight pre-processing step — NOT a replacement for
DESeq2's own normalization.  CPM is only used for classification
purposes so the user can define biologically meaningful sub-groups.

Functions
---------
compute_log2_cpm(counts_df)
    → Compute log2(CPM + 1) from raw counts.

classify_samples(log_cpm_df, rules)
    → Classify each sample based on marker gene expression rules.

merge_classification(metadata_df, classification, label, condition_col)
    → Merge classification result into metadata as compound condition.

build_compound_condition(metadata_df, class_col_name, …)
    → Build compound condition labels (e.g. tumor_TNBC) from condition + classification.

Example
-------
    from classification import compute_log2_cpm, classify_samples, merge_classification

    log_cpm = compute_log2_cpm(counts_df)
    classification = classify_samples(log_cpm, rules)
    metadata_df = merge_classification(metadata_df, classification, "TNBC", "condition")
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def compute_log2_cpm(counts_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute log2(CPM + 1) from a raw count matrix.

    This is a lightweight normalization used ONLY for sample
    classification — NOT for differential expression (DESeq2 does
    its own normalization).

    CPM (Counts Per Million) normalizes each sample's counts by
    its library size (total counts), making expression levels
    comparable across samples with different sequencing depths.

    The log2(x + 1) transformation compresses the dynamic range,
    making it easier to set meaningful thresholds for classification.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Raw count matrix (genes × samples).  Values must be
        non-negative integers.

    Returns
    -------
    pd.DataFrame
        log2(CPM + 1) matrix with the same shape and indices.

    Notes
    -----
    Memory: creates one float64 matrix.  For 25K genes × 800 samples
    this is ~150 MB.  Acceptable for a transient pre-processing step.
    """
    # Library sizes: sum of counts per sample (column sums)
    lib_sizes = counts_df.sum(axis=0).astype(np.float64)

    # Avoid division by zero for empty samples
    lib_sizes = lib_sizes.replace(0, 1)

    # CPM = (count / library_size) × 1,000,000
    cpm = counts_df.div(lib_sizes, axis=1) * 1e6

    # log2(CPM + 1) — the +1 handles zeros gracefully
    log_cpm = np.log2(cpm + 1)

    return log_cpm


def classify_samples(
    log_cpm_df: pd.DataFrame,
    rules: list[dict],
    positive_label: str = "positive",
    negative_label: str = "negative",
) -> pd.Series:
    """
    Classify each sample based on marker gene expression rules.

    A sample is classified as ``positive_label`` if ALL rules are
    satisfied, otherwise as ``negative_label``.

    Each rule specifies:
    - gene: gene identifier (must exist in log_cpm_df index)
    - threshold: numeric threshold for log2(CPM+1)
    - direction: "below" (gene < threshold → positive) or
                 "above" (gene > threshold → positive)

    For TNBC example (ESR1 low AND PGR low AND ERBB2 low):
    rules = [
        {"gene": "ESR1",  "threshold": 1.0, "direction": "below"},
        {"gene": "PGR",   "threshold": 1.0, "direction": "below"},
        {"gene": "ERBB2", "threshold": 1.0, "direction": "below"},
    ]

    Parameters
    ----------
    log_cpm_df : pd.DataFrame
        log2(CPM+1) matrix (genes × samples).
    rules : list[dict]
        List of classification rules.  Each dict must have keys:
        'gene', 'threshold', 'direction'.
    positive_label : str
        Label for samples that satisfy ALL rules.
    negative_label : str
        Label for samples that do NOT satisfy all rules.

    Returns
    -------
    pd.Series
        Classification for each sample (index = sample names).

    Raises
    ------
    ValueError
        If a gene in the rules is not found in the matrix.
    """
    if not rules:
        raise ValueError("At least one classification rule is required.")

    # Start with all True, then AND each rule
    mask = pd.Series(True, index=log_cpm_df.columns)

    for rule in rules:
        gene = rule["gene"]
        threshold = float(rule["threshold"])
        direction = rule["direction"]

        if gene not in log_cpm_df.index:
            raise ValueError(
                f"Gene '{gene}' not found in the count matrix. "
                f"Available genes are indexed by: {log_cpm_df.index.name or 'row index'}."
            )

        gene_expr = log_cpm_df.loc[gene]  # Series (samples)

        if direction == "below":
            rule_mask = gene_expr < threshold
        else:  # "above"
            rule_mask = gene_expr > threshold

        mask = mask & rule_mask

    result = pd.Series(negative_label, index=log_cpm_df.columns)
    result[mask] = positive_label

    return result


def merge_classification(
    metadata_df: pd.DataFrame,
    classification: pd.Series,
    class_col_name: str,
    condition_col: str = "condition",
) -> pd.DataFrame:
    """
    Merge expression-based classification into metadata as a compound
    condition column.

    Creates compound labels by combining the original condition with
    the classification result.  For example:

    - Original conditions: "tumor", "control"
    - Classification: "TNBC", "nonTNBC"
    - Compound: "tumor_TNBC", "tumor_nonTNBC", "control"

    Control/reference samples keep their original label (they are not
    reclassified).  Only non-reference samples get the compound label.

    Parameters
    ----------
    metadata_df : pd.DataFrame
        Metadata with the condition column.  Index = sample names.
    classification : pd.Series
        Classification result (index = sample names, values = labels).
    class_col_name : str
        Name for the new classification column added to metadata
        (e.g., "tnbc_status").
    condition_col : str
        Name of the existing condition column.

    Returns
    -------
    pd.DataFrame
        Updated metadata with:
        - New column ``class_col_name`` with the classification.
        - The ``condition_col`` updated with compound labels.
        - Original condition saved as ``{condition_col}_original``.
    """
    meta = metadata_df.copy()

    # Add classification as a new column
    meta[class_col_name] = classification.reindex(meta.index)

    # Save original condition
    original_col = f"{condition_col}_original"
    if original_col not in meta.columns:
        meta[original_col] = meta[condition_col]

    return meta


def build_compound_condition(
    metadata_df: pd.DataFrame,
    class_col_name: str,
    condition_col: str = "condition",
    reference_level: str | None = None,
) -> pd.DataFrame:
    """
    Build compound condition labels from original condition + classification.

    Samples whose original condition matches ``reference_level`` keep
    their original label.  All other samples get a compound label:
    ``{original_condition}_{classification}``.

    Parameters
    ----------
    metadata_df : pd.DataFrame
        Metadata with condition_col and class_col_name columns.
    class_col_name : str
        Name of the classification column.
    condition_col : str
        Name of the condition column.
    reference_level : str or None
        The reference/control condition.  Samples with this condition
        keep their label unchanged.  If None, all samples get compound
        labels.

    Returns
    -------
    pd.DataFrame
        Metadata with updated condition_col containing compound labels.
    """
    meta = metadata_df.copy()

    # Save original if not already saved
    original_col = f"{condition_col}_original"
    if original_col not in meta.columns:
        meta[original_col] = meta[condition_col]

    # Build compound labels
    def _make_compound(row):
        original = str(row[original_col])
        classification = str(row[class_col_name])
        if reference_level and original == reference_level:
            return original
        return f"{original}_{classification}"

    meta[condition_col] = meta.apply(_make_compound, axis=1)

    return meta
