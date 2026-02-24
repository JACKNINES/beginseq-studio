"""
deseq_runner.py — DESeq2 statistical model execution.

This module encapsulates ALL interaction with the pydeseq2 library.
It is the core of the analysis: it takes validated data and produces
the differential expression results table.

What does DESeq2 do internally?
-------------------------------
1. Count normalization (median-of-ratios / size factors).
2. Dispersion estimation (gene-by-gene variability).
3. Fitting a generalized linear model (GLM) with a negative
   binomial distribution.
4. Wald test to identify differentially expressed genes.
5. P-value correction for multiple comparisons (Benjamini-Hochberg).

Functions
---------
build_deseq_dataset(counts_df, metadata_df, condition_col, reference_level)
    -> Creates the DeseqDataSet object from the data.

run_deseq2(dds, progress_callback)
    -> Runs the DESeq2 analysis with step-by-step progress.

compute_contrast(dds, metadata_df, condition_col, reference_level, alpha)
    -> Computes the statistical contrast and returns the results table.

Usage example
--------------
    from deseq_runner import build_deseq_dataset, run_deseq2, compute_contrast

    dds = build_deseq_dataset(counts_df, metadata_df)
    dds = run_deseq2(dds)
    results_df = compute_contrast(dds, metadata_df, "condition", "control")
"""

from __future__ import annotations

import gc
import logging
import time
from typing import Callable

import numpy as np
import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from engine.config import DESEQ2_DEFAULTS, MEMORY_CONFIG

logger = logging.getLogger(__name__)

# Number of progress steps reported by run_deseq2 (used by analysis.py)
N_DESEQ2_STEPS = 6


def build_deseq_dataset(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    condition_col: str = DESEQ2_DEFAULTS["condition_col"],
    reference_level: str = DESEQ2_DEFAULTS["reference_level"],
) -> DeseqDataSet:
    """
    Builds the DeseqDataSet object that pydeseq2 needs to
    run the analysis.

    IMPORTANT: pydeseq2 expects the counts matrix TRANSPOSED
    (samples x genes), but our input is (genes x samples).
    This function handles the transposition.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Raw counts matrix (genes x samples).
        Must be previously validated with validation.py.
    metadata_df : pd.DataFrame
        Metadata aligned with the samples in the counts.
        Must be previously validated and aligned.
    condition_col : str
        Name of the condition column in metadata.
    reference_level : str
        Reference level for the contrast (e.g. "control").

    Returns
    -------
    DeseqDataSet
        DESeq2 object ready to run the analysis.

    Notes
    -----
    The `design` parameter defines the statistical model formula.
    "~ condition" means that the only factor explaining the variation
    in the counts is the experimental condition.

    The `ref_level` parameter tells DESeq2 which is the baseline level
    for computing fold-changes. The log2FC values are interpreted as
    "test_level vs reference_level".
    """
    # ── Memory-safe transpose ──────────────────────────────────────
    # counts_df.T creates a full copy (~190 MB for 25K genes × 800
    # samples at int32, much more if still float64).
    #
    # Strategy: build the transposed DataFrame directly from the numpy
    # array, then explicitly delete the intermediate to let the GC
    # reclaim memory before DeseqDataSet copies it internally.
    counts_t = pd.DataFrame(
        counts_df.values.T,
        index=counts_df.columns,
        columns=counts_df.index,
    )

    # ── Limit joblib parallelism ─────────────────────────────────
    # PyDESeq2 uses joblib for per-gene dispersion fitting.  Each
    # worker gets a copy of the full counts matrix (~163 MB for
    # 855×25K int64).  On an 8 GB M2, using all 8 cores would
    # allocate 8×163 MB = 1.3 GB *just* for worker copies.  Limit
    # to a sensible number to keep total memory under control.
    n_cpus = MEMORY_CONFIG.get("deseq2_n_cpus", 4)

    # For very large matrices, reduce parallelism further to avoid OOM.
    # Each worker gets a full copy of the counts matrix (int64 internally).
    # Matrix size in bytes: n_samples × n_genes × 8 (int64)
    n_genes = counts_df.shape[0]
    n_samples = counts_df.shape[1]
    matrix_bytes = n_samples * n_genes * 8  # int64

    if matrix_bytes > 300_000_000:  # > 300 MB per worker copy
        n_cpus = min(n_cpus, 2)
    elif matrix_bytes > 150_000_000:  # > 150 MB per worker copy
        n_cpus = min(n_cpus, 3)

    # ── Safety Catch for "no types given" (pandas/pydeseq2 issue) ──
    try:
        # Force metadata string conversion to avoid mixed-type issues in formulaic
        metadata_df = metadata_df.astype(str)

        dds = DeseqDataSet(
            counts=counts_t,
            metadata=metadata_df,
            design=f"~ {condition_col}",
            n_cpus=n_cpus,
            quiet=True,
        )
    except ValueError as e:
        # Check for the specific "no types given" error from pandas/formulaic
        if "no types given" in str(e):
            import sys
            print(f"\n!!! DESeq2 Crash Debug Info !!!", file=sys.stderr)
            print(f"Error: {e}", file=sys.stderr)
            print(f"Counts shape: {counts_t.shape}, dtypes: {counts_t.dtypes.iloc[0]}", file=sys.stderr)
            print(f"Metadata shape: {metadata_df.shape}", file=sys.stderr)
            print(f"Metadata dtypes:\n{metadata_df.dtypes}", file=sys.stderr)
            print(f"Condition col: {condition_col}", file=sys.stderr)
            print(f"First 5 rows of metadata:\n{metadata_df.head()}", file=sys.stderr)
            
            # Re-raise with a user-friendly message
            raise ValueError(
                f"Internal compatibility error (pandas/pydeseq2): {e}. "
                "Details have been logged for debugging. "
                "Try simplifying your sample names or verifying "
                "that there are no unusual characters in the metadata."
            ) from e
        raise e
    del counts_t
    gc.collect()

    return dds


def run_deseq2(
    dds: DeseqDataSet,
    progress_callback: Callable[[int, int, str], None] | None = None,
) -> tuple[DeseqDataSet, dict[str, float]]:
    """
    Runs the full DESeq2 pipeline on the dataset,
    reporting granular progress at each sub-step.

    Replicates exactly the internal logic of ``dds.deseq2()`` but
    reporting progress between each sub-step:

    1. fit_size_factors — Normalization factors (median-of-ratios).
    2. fit_genewise_dispersions — Gene-by-gene dispersions.
    3. fit_dispersion_trend + fit_dispersion_prior — Trend curve.
    4. fit_MAP_dispersions — MAP posterior dispersions.
    5. fit_LFC — Log-fold changes.
    6. calculate_cooks + refit + cooks_outlier — Outlier detection.

    Parameters
    ----------
    dds : DeseqDataSet
        DESeq2 object built with build_deseq_dataset().
    progress_callback : callable or None
        Called as ``progress_callback(current_step, total_steps, i18n_key)``
        before each sub-step starts. The i18n_key can be used to look up
        a translated message.

    Returns
    -------
    tuple[DeseqDataSet, dict[str, float]]
        - The fitted DeseqDataSet.
        - Dict mapping step names to elapsed seconds.
    """
    total = N_DESEQ2_STEPS
    step_timings: dict[str, float] = {}

    def _report(step: int, key: str):
        if progress_callback:
            progress_callback(step, total, key)

    # Step 1: Size factors (must pass the dds's own params)
    _report(0, "progress.size_factors")
    t0 = time.monotonic()

    # ── Automatic "poscounts" for sparse data ──
    # If every gene has at least one zero (common in large TCGA datasets),
    # the default "ratio" method fails and pydeseq2 falls back to
    # "iterative" which is extremely slow.  Force "poscounts" instead.
    fit_type = dds.size_factors_fit_type
    if (dds.X == 0).any(axis=0).all():
        fit_type = "poscounts"

    dds.fit_size_factors(
        fit_type=fit_type,
        control_genes=dds.control_genes,
    )
    step_timings["size_factors"] = time.monotonic() - t0
    gc.collect()

    # Step 2: Gene-wise dispersions (heaviest step — per-gene optimization)
    _report(1, "progress.genewise_disp")
    t0 = time.monotonic()
    dds.fit_genewise_dispersions()
    step_timings["genewise_disp"] = time.monotonic() - t0
    gc.collect()

    # Step 3: Dispersion trend + prior
    _report(2, "progress.disp_trend")
    t0 = time.monotonic()
    dds.fit_dispersion_trend()
    dds.fit_dispersion_prior()
    step_timings["disp_trend"] = time.monotonic() - t0

    # Step 4: MAP dispersions (second heaviest — per-gene optimization)
    _report(3, "progress.map_disp")
    t0 = time.monotonic()
    dds.fit_MAP_dispersions()
    step_timings["map_disp"] = time.monotonic() - t0
    gc.collect()

    # Step 5: Log-fold changes (per-gene IRLS)
    _report(4, "progress.fit_lfc")
    t0 = time.monotonic()
    dds.fit_LFC()
    step_timings["fit_lfc"] = time.monotonic() - t0
    gc.collect()

    # Step 6: Cook's distances + optional refit + outlier mask
    _report(5, "progress.cooks")
    t0 = time.monotonic()
    dds.calculate_cooks()
    if dds.refit_cooks:
        dds.refit()
    dds.cooks_outlier()
    step_timings["cooks"] = time.monotonic() - t0

    # Signal completion
    _report(total, "progress.deseq2_done")

    return dds, step_timings


def compute_contrast(
    dds: DeseqDataSet,
    metadata_df: pd.DataFrame,
    condition_col: str = DESEQ2_DEFAULTS["condition_col"],
    reference_level: str = DESEQ2_DEFAULTS["reference_level"],
    alpha: float = DESEQ2_DEFAULTS["alpha"],
    test_level: str | None = None,
    progress_callback: Callable[[int, int, str], None] | None = None,
) -> tuple[pd.DataFrame, str]:
    """
    Computes the statistical contrast between conditions and returns
    the results table.

    A "contrast" in DESeq2 compares two levels of the condition
    variable. For example: treated vs control.

    The results table contains for each gene:
    - baseMean: normalized mean of counts.
    - log2FoldChange: log2 of the expression change.
    - lfcSE: standard error of the log2FC.
    - stat: Wald statistic.
    - pvalue: raw p-value.
    - padj: adjusted p-value (Benjamini-Hochberg).

    Parameters
    ----------
    dds : DeseqDataSet
        DESeq2 object ALREADY fitted (after run_deseq2).
    metadata_df : pd.DataFrame
        Metadata with the condition column.
    condition_col : str
        Name of the condition column.
    reference_level : str
        Reference level (baseline).
    alpha : float
        Significance threshold for p-value adjustment.
    test_level : str, optional
        Condition level to compare against the reference.
        If not specified, the first available level is used.
    progress_callback : callable or None
        Called as ``progress_callback(current, total, i18n_key)``.

    Returns
    -------
    tuple[pd.DataFrame, str]
        - results_df: Full DESeq2 results table.
        - test_level: The condition level compared against the reference.
          Needed for the volcano plot labels.

    Raises
    ------
    ValueError
        If the condition levels cannot be identified.
    """
    # Identify the test level
    if test_level is None:
        # Use the first level that is not the reference
        condition_levels = metadata_df[condition_col].unique().tolist()
        condition_levels.remove(reference_level)
        test_level = condition_levels[0]

    if progress_callback:
        progress_callback(0, 1, "progress.wald_test")

    # Build and run the contrast
    stats = DeseqStats(
        dds,
        contrast=[condition_col, test_level, reference_level],
        alpha=alpha,
    )
    stats.summary()

    # Save MLE (unshrunk) log2FoldChange before attempting shrinkage
    results_df = stats.results_df.copy()
    results_df["log2FoldChange_MLE"] = results_df["log2FoldChange"].copy()

    # ── LFC shrinkage (apeGLM prior) ──────────────────────────────
    # Shrinks noisy fold-changes toward zero, especially for low-count
    # genes.  Improves volcano plot interpretability and reduces false
    # extremes.  P-values are left unchanged.
    shrinkage_applied = False
    try:
        coeff_name = f"{condition_col}[T.{test_level}]"
        stats.lfc_shrink(coeff=coeff_name)
        # lfc_shrink overwrites stats.results_df["log2FoldChange"] in-place
        results_df["log2FoldChange"] = stats.results_df["log2FoldChange"]
        results_df["lfcSE"] = stats.results_df["lfcSE"]
        shrinkage_applied = True
        logger.info("LFC shrinkage (apeGLM) applied successfully.")
    except Exception as e:
        # Shrinkage can fail with very few genes or singular designs.
        # Fall back to MLE estimates silently.
        logger.warning("LFC shrinkage failed, using MLE estimates: %s", e)

    # Store shrinkage status as DataFrame attribute
    results_df.attrs["shrinkage_applied"] = shrinkage_applied

    if progress_callback:
        progress_callback(1, 1, "progress.wald_done")

    return results_df, test_level
