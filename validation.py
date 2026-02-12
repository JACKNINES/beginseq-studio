"""
validation.py — Input data validation and normalization.

This module is responsible for verifying that user-provided data is correct
BEFORE executing any statistical analysis. Catching errors here saves time
and avoids cryptic error messages from pydeseq2.

Functions
---------
validate_counts_df(counts_df)
    → Validates the structure of the count matrix.

validate_metadata_df(metadata_df, condition_col)
    → Validates the structure of the metadata DataFrame.

align_samples(counts_df, metadata_df)
    → Verifies that samples match and aligns them.

normalize_indices(counts_df, metadata_df)
    → Converts all indices and columns to strings.

validate_condition_levels(metadata_df, condition_col, reference_level)
    → Verifies that conditions and the reference level are valid.

check_counts_are_raw(counts_df)
    → Heuristic to detect if the matrix appears normalized (FPKM/TPM/etc).

Usage example
--------------
    from validation import normalize_indices, align_samples

    counts_df, metadata_df = normalize_indices(counts_df, metadata_df)
    counts_df, metadata_df = align_samples(counts_df, metadata_df)
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from config import DESEQ2_DEFAULTS, GENE_FILTER_DEFAULTS, MEMORY_CONFIG


def validate_counts_df(counts_df: pd.DataFrame) -> None:
    """
    Validate that the count matrix has the expected structure.

    Checks performed:
    1. It is not empty.
    2. All values are numeric.
    3. It does not contain negative values (raw counts must be >= 0).

    Parameters
    ----------
    counts_df : pd.DataFrame
        Raw count matrix. Rows = genes, Columns = samples.

    Raises
    ------
    ValueError
        If any check fails, with a descriptive message.
    """
    if counts_df.empty:
        raise ValueError(
            "The count matrix is empty. "
            "Verify that the file was loaded correctly "
            "and that the separator (TSV=tab, CSV=comma) is appropriate."
        )

    # Verify that all values are numeric
    non_numeric = counts_df.select_dtypes(exclude=["number"])
    if not non_numeric.empty:
        bad_cols = list(non_numeric.columns[:5])
        raise ValueError(
            f"The count matrix contains non-numeric columns: "
            f"{bad_cols}. Make sure that the first column (genes) "
            f"is the index and not a data column."
        )

    # Verify that there are no negative values
    # Use numpy min() — O(1) extra memory instead of building a full
    # boolean matrix with (counts_df < 0).any().any()
    if counts_df.values.min() < 0:
        raise ValueError(
            "The count matrix contains negative values. "
            "DESeq2 requires raw counts (integers >= 0)."
        )


def check_counts_are_raw(counts_df: pd.DataFrame) -> dict:
    """
    Heuristic to detect if the count matrix appears to contain
    normalized values (FPKM, TPM, CPM, log-transformed) instead of
    raw counts.

    DESeq2 REQUIRES raw integer counts. Normalized matrices produce
    statistically invalid results.

    Warning signs:
    - High proportion of values with decimals (raw counts are integers).
    - Very low maximum values (suggests log-normalization).
    - Per-sample values summing to ~1e6 (CPM) or ~1e6 (TPM).

    Parameters
    ----------
    counts_df : pd.DataFrame
        Count matrix (genes x samples).

    Returns
    -------
    dict with keys:
        - "is_suspect" : bool
            True if the matrix appears to be normalized data.
        - "reason" : str | None
            Description of the detected problem, or None if everything is OK.
        - "pct_decimal" : float
            Percentage of values that have decimals (0-100).
    """
    # Sample a subset of values for heuristic check instead of
    # materialising the entire matrix as a flat array.  For a
    # 60K × 800 matrix, ravel() alone creates a 384 MB float64 copy.
    vals = counts_df.values
    n_total = vals.size

    # If matrix is small, check all.  Otherwise, sample ~100K values.
    SAMPLE_LIMIT = 100_000
    if n_total <= SAMPLE_LIMIT:
        flat = vals.ravel()
    else:
        rng = np.random.default_rng(42)
        idx_rows = rng.integers(0, vals.shape[0], size=SAMPLE_LIMIT)
        idx_cols = rng.integers(0, vals.shape[1], size=SAMPLE_LIMIT)
        flat = vals[idx_rows, idx_cols]

    nonzero = flat[flat != 0]

    if len(nonzero) == 0:
        return {"is_suspect": False, "reason": None, "pct_decimal": 0.0}

    # 1. What proportion of non-zero values have a decimal part?
    has_decimal = nonzero != np.floor(nonzero)
    pct_decimal = float(has_decimal.sum() / len(nonzero) * 100)

    # 2. Is the maximum suspiciously low? (log-normalized is usually < 20)
    max_val = float(np.max(nonzero))

    # 3. Detect if it appears log-transformed (majority of values < 20 with decimals)
    if pct_decimal > 50 and max_val < 25:
        return {
            "is_suspect": True,
            "reason": (
                f"The values appear to be **log-normalized** "
                f"({pct_decimal:.0f}% have decimals, max value = {max_val:.1f}). "
                f"DESeq2 requires raw integer counts."
            ),
            "pct_decimal": pct_decimal,
        }

    # 4. Detect FPKM/TPM/CPM: many decimals but larger values
    if pct_decimal > 50:
        return {
            "is_suspect": True,
            "reason": (
                f"The values appear to be **normalized** (FPKM, TPM, or CPM) — "
                f"{pct_decimal:.0f}% of non-zero values have decimals. "
                f"DESeq2 requires raw integer counts, not normalized values."
            ),
            "pct_decimal": pct_decimal,
        }

    # 5. Moderate percentage of decimals: soft warning
    if pct_decimal > 10:
        return {
            "is_suspect": True,
            "reason": (
                f"**{pct_decimal:.0f}%** of non-zero values have decimals. "
                f"Raw counts should be integers. This matrix may contain "
                f"partially normalized data. Verify your input."
            ),
            "pct_decimal": pct_decimal,
        }

    return {"is_suspect": False, "reason": None, "pct_decimal": pct_decimal}


def _normalize_columns(metadata_df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize metadata column names:
    1. Remove leading and trailing whitespace (strip).
    2. Convert everything to lowercase.

    This allows files with columns like "Condition",
    "CONDITION", "  condition  " or "Sample " to work without errors.

    Parameters
    ----------
    metadata_df : pd.DataFrame
        DataFrame with unnormalized columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with normalized columns (lowercase, no whitespace).
    """
    metadata_df = metadata_df.copy()
    metadata_df.columns = metadata_df.columns.str.strip().str.lower()
    return metadata_df


def validate_metadata_df(
    metadata_df: pd.DataFrame,
    condition_col: str = DESEQ2_DEFAULTS["condition_col"],
) -> pd.DataFrame:
    """
    Validate the structure of the metadata DataFrame and set the index.

    Checks performed:
    1. It is not empty.
    2. Normalizes column names (strip + lowercase) so that
       "Condition", "CONDITION", " condition " are equivalent.
    3. If it has a "sample" column, uses it as the index.
    4. Contains the specified condition column.
    5. Normalizes the values of the condition column (strip)
       to avoid issues with invisible whitespace.

    Parameters
    ----------
    metadata_df : pd.DataFrame
        Metadata DataFrame with sample information.
    condition_col : str
        Name of the column defining experimental conditions.
        The search is case-insensitive thanks to normalization.

    Returns
    -------
    pd.DataFrame
        Validated metadata, with normalized columns and the index
        set correctly.

    Raises
    ------
    ValueError
        If any check fails.
    """
    if metadata_df.empty:
        raise ValueError(
            "The metadata file is empty. "
            "Verify that the file contains at least one row of data."
        )

    # Normalize column names: strip + lowercase
    metadata_df = _normalize_columns(metadata_df)

    # Set "sample" column as the index if it exists
    if "sample" in metadata_df.columns:
        metadata_df = metadata_df.set_index("sample")

    # Verify that the condition column exists
    if condition_col not in metadata_df.columns:
        available = list(metadata_df.columns)
        raise ValueError(
            f"Column '{condition_col}' does not exist in the metadata. "
            f"Available columns: {available}. "
            f"Verify that your file has a column named "
            f"'{condition_col}' or specify the correct name."
        )

    # Normalize values of the condition column (strip whitespace)
    metadata_df[condition_col] = (
        metadata_df[condition_col].astype(str).str.strip()
    )

    return metadata_df


def normalize_indices(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Convert all indices and column names to strings and remove
    leading/trailing whitespace.

    This prevents comparison errors when:
    - Sample names are numbers that pandas interprets as integers
      in one file but as strings in another.
    - There is invisible trailing whitespace in names (e.g. "GSM123 ").

    **Memory note:** We mutate the index/columns in-place on copies of
    the *metadata* only (small). For the counts DataFrame we update
    the index/column objects without copying the underlying data.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Raw count matrix.
    metadata_df : pd.DataFrame
        Metadata DataFrame.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        Tuple with (counts_df, metadata_df) with normalized indices.
    """
    # Index / column objects are lightweight — replacing them does NOT
    # copy the entire underlying data array.  However we must avoid
    # mutating the *caller's* DataFrame.  pd.DataFrame._set_axis is
    # cheap (no data copy), but we first wrap in a new DataFrame
    # header that shares the same block manager (zero-copy).
    new_idx = counts_df.index.astype(str).str.strip()
    new_cols = counts_df.columns.astype(str).str.strip()

    # Only create a thin wrapper if the index/columns actually changed
    if not (counts_df.index.equals(new_idx) and counts_df.columns.equals(new_cols)):
        counts_df = counts_df.copy(deep=False)   # shallow: shares data blocks
        counts_df.index = new_idx
        counts_df.columns = new_cols

    metadata_df = metadata_df.copy()
    metadata_df.index = metadata_df.index.astype(str).str.strip()

    return counts_df, metadata_df


def check_sample_overlap(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
) -> dict:
    """
    Analyze the sample overlap between the count matrix and the
    metadata. Does NOT raise errors; returns a dictionary with the
    full diagnostic so the UI can decide what to do.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Count matrix (columns = samples).
    metadata_df : pd.DataFrame
        Metadata (index = samples).

    Returns
    -------
    dict with the following keys:
        - "match" : bool
            True if the samples match exactly.
        - "common" : set[str]
            Samples present in BOTH files (the intersection).
        - "only_in_counts" : set[str]
            Samples that are in counts but NOT in metadata.
        - "only_in_metadata" : set[str]
            Samples that are in metadata but NOT in counts.
        - "n_common" : int
            Number of samples in common.
        - "n_counts" : int
            Total number of samples in the count matrix.
        - "n_metadata" : int
            Total number of samples in the metadata.
    """
    counts_samples = set(counts_df.columns)
    metadata_samples = set(metadata_df.index)

    common = counts_samples & metadata_samples
    only_in_counts = counts_samples - metadata_samples
    only_in_metadata = metadata_samples - counts_samples

    return {
        "match": counts_samples == metadata_samples,
        "common": common,
        "only_in_counts": only_in_counts,
        "only_in_metadata": only_in_metadata,
        "n_common": len(common),
        "n_counts": len(counts_samples),
        "n_metadata": len(metadata_samples),
    }


def align_samples(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    allow_subset: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Align the samples between counts and metadata.

    If allow_subset=False (default), requires an exact match.
    If allow_subset=True, filters both DataFrames to keep only
    the common samples (intersection). This is the standard
    behavior in most real pipelines.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Count matrix (columns = samples).
    metadata_df : pd.DataFrame
        Metadata (index = samples).
    allow_subset : bool, default False
        If True, allows samples to not match exactly and filters
        both DataFrames to the intersection.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        Tuple (counts_df, metadata_df) with aligned samples.

    Raises
    ------
    ValueError
        If samples do not match and allow_subset=False.
        If there are NO common samples at all (regardless of
        allow_subset).
    """
    overlap = check_sample_overlap(counts_df, metadata_df)

    # Regardless of the mode, if there are no common samples it is a fatal error
    if overlap["n_common"] == 0:
        raise ValueError(
            "There are NO common samples between the count matrix "
            "and the metadata file. Verify that the sample names "
            "match between both files."
        )

    if not overlap["match"]:
        if not allow_subset:
            # Strict mode: raise error (should not reach here
            # if app.py handles the flow correctly)
            msg_parts = [
                "Samples do not match between the count matrix "
                "and the metadata file."
            ]
            if overlap["only_in_counts"]:
                msg_parts.append(
                    f"  • Samples in counts but NOT in metadata: "
                    f"{sorted(overlap['only_in_counts'])}"
                )
            if overlap["only_in_metadata"]:
                msg_parts.append(
                    f"  • Samples in metadata but NOT in counts: "
                    f"{sorted(overlap['only_in_metadata'])}"
                )
            raise ValueError("\n".join(msg_parts))

        # Subset mode: filter both DataFrames to the intersection
        common_samples = sorted(overlap["common"])
        counts_df = counts_df[common_samples]
        metadata_df = metadata_df.loc[common_samples]
    else:
        # Exact match: just reorder
        metadata_df = metadata_df.loc[counts_df.columns]

    return counts_df, metadata_df


def get_condition_levels(
    metadata_df: pd.DataFrame,
    condition_col: str = DESEQ2_DEFAULTS["condition_col"],
) -> list[str]:
    """
    Extract and validate the unique levels of the condition column.

    Only verifies that there are at least 2 levels (required for a
    contrast). Does NOT require a level with a specific name like
    "control". The user chooses the reference level in the UI.

    Parameters
    ----------
    metadata_df : pd.DataFrame
        Metadata with the condition column.
    condition_col : str
        Name of the condition column.

    Returns
    -------
    list[str]
        Sorted list of all unique levels found.

    Raises
    ------
    ValueError
        If there are fewer than 2 levels (cannot perform a contrast).
    """
    levels = sorted(metadata_df[condition_col].unique().tolist())

    if len(levels) < 2:
        raise ValueError(
            f"At least 2 levels are needed in the column "
            f"'{condition_col}' to perform a contrast. "
            f"Only found: {levels}."
        )

    return levels


def validate_condition_levels(
    metadata_df: pd.DataFrame,
    condition_col: str = DESEQ2_DEFAULTS["condition_col"],
    reference_level: str = DESEQ2_DEFAULTS["reference_level"],
) -> list[str]:
    """
    Validate that the reference level chosen by the user exists
    in the data and return the test levels.

    Parameters
    ----------
    metadata_df : pd.DataFrame
        Metadata with the condition column.
    condition_col : str
        Name of the condition column.
    reference_level : str
        Reference level chosen by the user.

    Returns
    -------
    list[str]
        List of condition levels that are NOT the reference.

    Raises
    ------
    ValueError
        If the reference level does not exist or there are fewer than 2 levels.
    """
    levels = get_condition_levels(metadata_df, condition_col)

    if reference_level not in levels:
        raise ValueError(
            f"The reference level '{reference_level}' was not found "
            f"in the column '{condition_col}'. "
            f"Available levels: {levels}."
        )

    test_levels = [lvl for lvl in levels if lvl != reference_level]

    return test_levels


# ──────────────────────────────────────────────────────────────────────
# Gene pre-filtering
# ──────────────────────────────────────────────────────────────────────

def filter_low_expression_genes(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame | None = None,
    condition_col: str = DESEQ2_DEFAULTS["condition_col"],
    min_total_count: int = GENE_FILTER_DEFAULTS["min_total_count"],
    min_samples_expressing: int = GENE_FILTER_DEFAULTS["min_samples_expressing"],
    min_count_per_sample: int = GENE_FILTER_DEFAULTS["min_count_per_sample"],
) -> tuple[pd.DataFrame, dict]:
    """
    Filter out lowly-expressed genes before running DESeq2.

    This is **standard practice** in RNA-seq differential expression
    analysis (recommended by Love, Huber & Anders 2014). Genes with
    very low counts across all samples:

      - Add noise to dispersion estimation.
      - Waste statistical power through unnecessary multiple-testing
        corrections.
      - Can cause numerical instability in the GLM fitting.

    Filtering criteria
    ------------------
    A gene is **kept** if BOTH conditions are met:

    1. **Total count ≥ min_total_count**
       Sum of counts across all samples must reach a minimum.
       Default: 10 (equivalent to roughly 1 CPM in a typical library).

    2. **Expressed in ≥ min_samples_expressing samples**
       The gene must have ≥ ``min_count_per_sample`` count in at least
       this many samples.

       If ``min_samples_expressing == 0`` (auto mode), the threshold
       is set to ``max(smallest_group, n_samples × auto_fraction)``.
       The ``auto_fraction`` (default 0.5) ensures that on large TCGA
       datasets (~800+ samples) the filter is aggressive enough to
       reduce ~60K genes to ~20-25K, instead of barely filtering
       anything when using only the smallest group size.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Raw count matrix (genes × samples).
    metadata_df : pd.DataFrame or None
        Metadata with condition column.  Only needed when
        ``min_samples_expressing == 0`` (auto mode) so we can
        determine the smallest group size.
    condition_col : str
        Name of the condition column in ``metadata_df``.
    min_total_count : int
        Minimum total count threshold.
    min_samples_expressing : int
        Minimum number of samples where the gene must be "expressed"
        (count ≥ ``min_count_per_sample``).
        0 = auto: ``max(smallest_group, n_samples × auto_fraction)``.
    min_count_per_sample : int
        Minimum count in a sample for the gene to be considered
        "expressed" there.  Default 1 (any non-zero count).

    Returns
    -------
    tuple[pd.DataFrame, dict]
        - Filtered counts DataFrame (genes × samples).
        - Stats dict with keys:
          - ``n_before``: genes before filtering.
          - ``n_after``: genes after filtering.
          - ``n_removed``: genes removed.
          - ``min_total_count``: threshold used.
          - ``min_samples_expressing``: threshold used (resolved).
          - ``min_count_per_sample``: threshold used.
    """
    n_before = counts_df.shape[0]
    n_samples = counts_df.shape[1]

    # Resolve auto min_samples_expressing
    if min_samples_expressing == 0:
        # Start with the smallest condition group
        smallest_group = 2  # fallback
        if metadata_df is not None and condition_col in metadata_df.columns:
            group_sizes = metadata_df[condition_col].value_counts()
            smallest_group = int(group_sizes.min())

        # For large datasets, also consider a fraction of total samples.
        # On TCGA datasets (800+ samples), the smallest group (e.g. 135)
        # is too low — most genes are expressed in that many samples.
        # Using a higher threshold ensures meaningful gene reduction.
        auto_fraction = GENE_FILTER_DEFAULTS.get("auto_samples_fraction", 0.5)
        fraction_based = int(n_samples * auto_fraction)

        min_samples_expressing = max(smallest_group, fraction_based)

        # Clamp to not exceed total samples
        min_samples_expressing = min(min_samples_expressing, n_samples)

    # ── Criterion 1: total count per gene ─────────────────────────
    # Use the underlying numpy array directly — avoids pandas overhead
    # and prevents creation of intermediate Series objects.
    arr = np.asarray(counts_df)            # zero-copy view when numeric
    gene_totals = arr.sum(axis=1)          # 1-D ndarray, ~480 KB for 60K genes
    mask_total = gene_totals >= min_total_count

    # ── Criterion 2: number of samples with count ≥ threshold ─────
    # When min_count_per_sample == 1, np.count_nonzero is optimal
    # (no boolean matrix).  For higher thresholds, we compare directly.
    if min_count_per_sample <= 1:
        samples_expressing = np.count_nonzero(arr, axis=1)
    else:
        # Count samples where count >= min_count_per_sample.
        # Use sum of boolean — creates a temporary bool array but it's
        # the only way to threshold per-element without a loop.
        samples_expressing = (arr >= min_count_per_sample).sum(axis=1)
    mask_samples = samples_expressing >= min_samples_expressing
    del arr  # release reference

    # ── Apply both ────────────────────────────────────────────────
    keep_mask = mask_total & mask_samples
    filtered_df = counts_df.loc[keep_mask]

    n_after = filtered_df.shape[0]

    stats = {
        "n_before": n_before,
        "n_after": n_after,
        "n_removed": n_before - n_after,
        "min_total_count": min_total_count,
        "min_samples_expressing": min_samples_expressing,
        "min_count_per_sample": min_count_per_sample,
    }

    return filtered_df, stats


def downcast_counts(counts_df: pd.DataFrame) -> pd.DataFrame:
    """
    Downcast count matrix to the most compact integer type.

    RNA-seq raw counts are non-negative integers.  The default
    ``int64`` (8 bytes per value) is massive overkill.  ``int32``
    supports values up to 2,147,483,647 — more than enough for any
    realistic RNA-seq experiment.  This halves memory usage.

    For a TCGA dataset (25 000 genes × 800 samples = 20 million values),
    this saves ~150 MB of RAM.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Count matrix with integer values.

    Returns
    -------
    pd.DataFrame
        Same matrix with int32 dtype.
    """
    target_dtype = MEMORY_CONFIG.get("counts_dtype", "int32")
    if counts_df.dtypes.eq(target_dtype).all():
        return counts_df  # already optimal
    return counts_df.astype(target_dtype, copy=False)


# ──────────────────────────────────────────────────────────────────────
# Pre-analysis checks (warnings, not errors)
# ──────────────────────────────────────────────────────────────────────

def check_class_imbalance(
    metadata_df: pd.DataFrame,
    condition_col: str = "condition",
    threshold_ratio: float = 5.0,
) -> dict | None:
    """
    Check if condition groups are severely imbalanced.

    Parameters
    ----------
    metadata_df : pd.DataFrame
        Metadata with condition column.
    condition_col : str
        Name of condition column.
    threshold_ratio : float
        Ratio above which imbalance is flagged (default 5:1).

    Returns
    -------
    dict or None
        If imbalanced: dict with keys 'ratio', 'large_group', 'large_n',
        'small_group', 'small_n'.  Otherwise None.
    """
    if condition_col not in metadata_df.columns:
        return None
    group_sizes = metadata_df[condition_col].value_counts()
    if len(group_sizes) < 2:
        return None
    max_n = int(group_sizes.max())
    min_n = int(group_sizes.min())
    if min_n == 0:
        return None
    ratio = max_n / min_n
    if ratio > threshold_ratio:
        return {
            "ratio": round(ratio, 1),
            "large_group": str(group_sizes.idxmax()),
            "large_n": max_n,
            "small_group": str(group_sizes.idxmin()),
            "small_n": min_n,
        }
    return None


def check_sparse_data(counts_df: pd.DataFrame) -> bool:
    """
    Check if the count matrix is sparse (many genes have zeros).

    A dataset is considered sparse if every gene has at least one zero
    across all samples.  This is common in large TCGA datasets and
    causes the default 'ratio' size-factor method in PyDESeq2 to fail,
    falling back to an extremely slow 'iterative' method.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Count matrix (genes × samples).

    Returns
    -------
    bool
        True if every gene has at least one zero (all columns have a zero
        somewhere), which means 'poscounts' normalisation should be used.
    """
    arr = np.asarray(counts_df)
    # Check: does every gene (row) have at least one zero?
    return bool((arr == 0).any(axis=1).all())
