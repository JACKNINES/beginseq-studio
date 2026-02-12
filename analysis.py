"""
analysis.py — Orchestrator for the complete analysis pipeline.

This module acts as a FACADE (Facade Pattern) that coordinates all
specialized modules to execute the differential expression pipeline
from start to finish.

Before refactoring, this file contained ALL the logic (~140 lines).
Now it delegates each responsibility to the corresponding module:

    analysis.py (orchestrator)
        ├── validation.py       → Validate data
        ├── deseq_runner.py     → Run DESeq2
        └── visualization.py    → Generate plots

NOTE: Data is expected to arrive pre-processed from app.py:
- Normalized columns (lowercase, stripped).
- Normalized indices (str, stripped).
- Samples already filtered/aligned.
- reference_level already validated by the UI.

The Step 1 validations are DEFENSIVE (idempotent): if the data is
already clean, they do not break anything. This allows the pipeline
to also be used standalone without the UI.

Functions
---------
run_deseq2_pipeline(counts_df, metadata_df, ...)
    → Runs the complete pipeline and returns results + volcano plot.

Usage Example
-------------
    from analysis import run_deseq2_pipeline

    results_df, volcano_fig = run_deseq2_pipeline(counts_df, metadata_df)
"""

from __future__ import annotations

import gc
import time
from typing import Callable

import pandas as pd
import matplotlib.figure

from config import DESEQ2_DEFAULTS
from config import GENE_FILTER_DEFAULTS
from config import TIME_ESTIMATION
from validation import (
    validate_counts_df,
    validate_metadata_df,
    normalize_indices,
    align_samples,
    validate_condition_levels,
    filter_low_expression_genes,
    downcast_counts,
)
from deseq_runner import (
    build_deseq_dataset, run_deseq2, compute_contrast, N_DESEQ2_STEPS,
)
from visualization import (
    prepare_volcano_data, create_volcano_plot,
    prepare_pca_data, create_pca_plot,
    prepare_ma_data, create_ma_plot,
    prepare_heatmap_data, create_heatmap,
)


def _format_seconds_human(total: float) -> str:
    """Format seconds as a human-readable '~X min Y s' string."""
    if total < 60:
        return f"~{int(total)} s"
    elif total < 3600:
        mins = int(total // 60)
        secs = int(total % 60)
        return f"~{mins} min {secs} s" if secs > 0 else f"~{mins} min"
    else:
        hrs = int(total // 3600)
        mins = int((total % 3600) // 60)
        return f"~{hrs} h {mins} min"


def estimate_pipeline_time(
    n_samples: int, n_genes: int,
) -> dict:
    """
    Estimate total and per-step wall-clock time for the DESeq2 pipeline.

    Uses empirical coefficients from ``config.TIME_ESTIMATION`` calibrated
    on Apple M2.  For large datasets (> 20M elements), applies a
    super-linear scaling factor to account for memory pressure and
    swap activity on 8 GB machines.

    Parameters
    ----------
    n_samples : int
        Number of samples in the dataset.
    n_genes : int
        Number of genes (after filtering, if applied).

    Returns
    -------
    dict
        - "total_seconds": float — estimated total wall time.
        - "total_human": str — human-readable string, e.g. "~2 min 30 s".
        - "per_step": dict[str, float] — estimated seconds per step.
        - "n_samples": int
        - "n_genes": int
    """
    M = (n_samples * n_genes) / 1_000_000  # millions of elements

    # For large datasets, apply super-linear scaling to account for
    # memory pressure: PyDESeq2 creates multiple float64 copies of the
    # full matrix, and on 8 GB machines this causes swap → slowdown.
    threshold = TIME_ESTIMATION.get("superlinear_threshold_M", 20.0)
    exponent = TIME_ESTIMATION.get("superlinear_exponent", 1.3)
    if M > threshold:
        # Scale the excess portion super-linearly
        effective_M = threshold + (M - threshold) ** exponent
    else:
        effective_M = M

    per_step = {}
    for key in TIME_ESTIMATION["step_keys"]:
        coeff = TIME_ESTIMATION.get(f"{key}_per_M", 0.0)
        per_step[key] = effective_M * coeff

    total = TIME_ESTIMATION["fixed_overhead_s"] + sum(per_step.values())

    return {
        "total_seconds": total,
        "total_human": _format_seconds_human(total),
        "per_step": per_step,
        "n_samples": n_samples,
        "n_genes": n_genes,
    }


def estimate_pipeline_time_with_filter(
    n_samples: int, n_genes_raw: int, filter_genes: bool,
) -> dict:
    """
    Estimate pipeline time, accounting for expected gene filtering.

    When gene filtering is enabled and we don't yet know the post-filter
    count, we use the typical_filter_ratio from config to estimate the
    effective gene count that DESeq2 will process.

    Parameters
    ----------
    n_samples : int
    n_genes_raw : int
        Raw gene count before filtering.
    filter_genes : bool
        Whether gene filtering will be applied.

    Returns
    -------
    dict
        Same as estimate_pipeline_time, plus:
        - "n_genes_estimated": int — estimated post-filter gene count.
        - "filter_applied": bool
    """
    if filter_genes and n_genes_raw > 30_000:
        ratio = TIME_ESTIMATION.get("typical_filter_ratio", 0.42)
        n_genes_est = int(n_genes_raw * ratio)
    else:
        n_genes_est = n_genes_raw

    est = estimate_pipeline_time(n_samples, n_genes_est)
    est["n_genes_estimated"] = n_genes_est
    est["n_genes_raw"] = n_genes_raw
    est["filter_applied"] = filter_genes and n_genes_raw > 30_000
    return est


def run_deseq2_pipeline(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    condition_col: str = DESEQ2_DEFAULTS["condition_col"],
    reference_level: str = DESEQ2_DEFAULTS["reference_level"],
    test_level: str | None = None,
    alpha: float = DESEQ2_DEFAULTS["alpha"],
    log2fc_threshold: float = DESEQ2_DEFAULTS["log2fc_threshold"],
    progress_callback: Callable[[int, int, str], None] | None = None,
    estimate_callback: Callable[[dict], None] | None = None,
    # Gene pre-filtering parameters
    filter_genes: bool = GENE_FILTER_DEFAULTS["enabled"],
    min_total_count: int = GENE_FILTER_DEFAULTS["min_total_count"],
    min_samples_expressing: int = GENE_FILTER_DEFAULTS["min_samples_expressing"],
    min_count_per_sample: int = GENE_FILTER_DEFAULTS["min_count_per_sample"],
    # Visualization parameters
    min_base_mean: float = DESEQ2_DEFAULTS.get("min_base_mean", 0),
    batch_col: str | None = None,
    legend_loc: str = "upper right",
) -> tuple[pd.DataFrame, dict]:
    """
    Run the complete differential expression analysis pipeline
    with DESeq2.

    This is the main entry point for the analysis. It coordinates
    the following steps in order:

    Step 1 — Defensive validation (validation.py)
    ├── Validate counts matrix (structure, types, values)
    ├── Validate metadata (structure, condition column)
    ├── Normalize indices (idempotent if already normalized)
    ├── Align samples (allow_subset=True to be tolerant)
    ├── Validate condition levels (reference exists, ≥2 levels)
    ├── **Filter low-expression genes** (reduces RAM + improves power)
    └── **Downcast to int32** (reduces RAM 50%)

    Step 2 — Statistical analysis (deseq_runner.py)
    ├── Build DeseqDataSet object
    ├── Run DESeq2 (normalization + model)
    └── Compute contrast (test vs reference)

    Step 3 — Visualization (visualization.py)
    ├── PCA (with VST if available, top-variance genes)
    ├── Volcano plot
    ├── MA plot
    └── Heatmap (auto-capped genes, adaptive clustering)

    Parameters
    ----------
    counts_df : pd.DataFrame
        Raw RNA-seq counts matrix.
        - Rows: genes (gene_id as index).
        - Columns: samples.
        - Values: integers >= 0 (raw counts, NOT normalized).
        Can arrive pre-processed from app.py or raw from a script.
    metadata_df : pd.DataFrame
        Sample metadata.
        - Must have a "sample" column (or samples as index).
        - Must have a column with the experimental conditions.
    condition_col : str, default "condition"
        Name of the metadata column that defines the experimental
        conditions (e.g. "treated" vs "control").
    reference_level : str, default "control"
        Reference level for the contrast. Fold-changes are computed
        RELATIVE to this level.
    test_level : str, optional
        Condition level to compare against the reference. If not
        specified, the first level that is not the reference is used.
        Useful when there are more than 2 conditions.
    alpha : float, default 0.05
        Adjusted p-value threshold (Benjamini-Hochberg) to consider
        a gene as statistically significant.
    log2fc_threshold : float, default 1.0
        Absolute log2 fold-change threshold to consider a gene as
        biologically relevant.
    progress_callback : callable or None
        Called as ``progress_callback(current_step, total_steps, i18n_key)``
        to report progress. Total steps = 11 (validation + filtering +
        6 DESeq2 + Wald test + visualizations + done).
    estimate_callback : callable or None
        Called once after gene filtering with the actual time estimate dict
        (from ``estimate_pipeline_time``) based on the real post-filter
        dimensions. The UI can use this to recalculate the ETA.
    filter_genes : bool
        Whether to apply gene pre-filtering before DESeq2.
    min_total_count : int
        Minimum total counts per gene (sum across samples).
    min_samples_expressing : int
        Minimum samples with ≥ min_count_per_sample count.
        0 = auto: max(smallest_group, n_samples × auto_fraction).
    min_count_per_sample : int
        Minimum count per sample for a gene to be "expressed" there.

    Returns
    -------
    tuple[pd.DataFrame, dict]
        - results_df: DESeq2 results table with columns:
            - baseMean: normalized mean of counts.
            - log2FoldChange: log2 of expression change.
            - lfcSE: standard error of log2FC.
            - stat: Wald statistic.
            - pvalue: raw p-value.
            - padj: adjusted p-value (BH).
        - figures: Dictionary with the generated figures:
            - "volcano": Volcano plot.
            - "pca": Sample PCA plot.
            - "ma": MA plot.
            - "heatmap": Top genes heatmap.
            - "filter_stats": dict with gene filtering statistics
              (n_before, n_after, n_removed, thresholds).

    Raises
    ------
    ValueError
        If the input data does not pass validation.

    Example
    -------
        import pandas as pd
        from analysis import run_deseq2_pipeline

        counts = pd.read_csv("counts.tsv", sep="\\t", index_col=0)
        metadata = pd.read_csv("metadata.csv")

        results, figures = run_deseq2_pipeline(
            counts, metadata,
            condition_col="treatment",
            reference_level="placebo",
            alpha=0.01,
            log2fc_threshold=1.5,
        )

        print(results.head())
        figures["volcano"].savefig("volcano.png")
        figures["pca"].savefig("pca.png")
        figures["ma"].savefig("ma.png")
        figures["heatmap"].savefig("heatmap.png")
    """

    # Total pipeline steps:
    # 1 validation + 1 gene_filter + N DESeq2 sub-steps + 1 Wald test + 1 visualizations + 1 done
    TOTAL_STEPS = 1 + 1 + N_DESEQ2_STEPS + 1 + 1 + 1  # = 11
    step = 0
    filter_stats = {}
    pipeline_t0 = time.monotonic()

    def _report(i18n_key: str):
        """Report progress for the current step, then advance."""
        nonlocal step
        if progress_callback:
            progress_callback(step, TOTAL_STEPS, i18n_key)
        step += 1

    # ── Step 1: Defensive validation ──────────────────────────────
    _report("progress.validating")

    validate_counts_df(counts_df)
    metadata_df = validate_metadata_df(metadata_df, condition_col)
    counts_df, metadata_df = normalize_indices(counts_df, metadata_df)

    counts_df, metadata_df = align_samples(
        counts_df, metadata_df, allow_subset=True
    )

    available_test_levels = validate_condition_levels(
        metadata_df, condition_col, reference_level
    )

    if test_level is None:
        test_level = available_test_levels[0]
    elif test_level not in available_test_levels:
        raise ValueError(
            f"The test level '{test_level}' is not valid. "
            f"Available options: {available_test_levels}"
        )

    # ── Step 1b: Low-expression gene filtering ────────────────────
    _report("progress.filtering_genes")

    if filter_genes:
        counts_df, filter_stats = filter_low_expression_genes(
            counts_df,
            metadata_df=metadata_df,
            condition_col=condition_col,
            min_total_count=min_total_count,
            min_samples_expressing=min_samples_expressing,
            min_count_per_sample=min_count_per_sample,
        )

    # Downcast to int32 to halve memory footprint
    counts_df = downcast_counts(counts_df)
    gc.collect()  # Free pre-filter / pre-downcast arrays

    # ── Compute time estimate AFTER gene filtering ─────────────────
    # Now we know the actual dataset dimensions that DESeq2 will see.
    n_samples_final = counts_df.shape[1]
    n_genes_final = counts_df.shape[0]
    time_est = estimate_pipeline_time(n_samples_final, n_genes_final)

    # Notify the UI with accurate post-filter estimate so it can
    # recalculate the ETA display.
    if estimate_callback:
        estimate_callback(time_est)

    # ── Prepare visualization subset BEFORE building DDS ─────────────
    # We need counts_df for visualizations later, but only for the two
    # selected conditions.  Extract that subset now (lightweight column
    # selection, no data copy) so we can free counts_df right after
    # building the DeseqDataSet — saving ~80-160 MB of RAM during the
    # expensive dispersion-fitting steps.
    selected_conditions = [reference_level, test_level]
    mask = metadata_df[condition_col].isin(selected_conditions)
    metadata_subset = metadata_df[mask].copy()
    counts_subset = counts_df[metadata_subset.index].copy()

    # ── Step 2: Statistical analysis ────────────────────────────────
    dds = build_deseq_dataset(
        counts_df, metadata_df, condition_col, reference_level
    )
    # Free the original genes×samples matrix — DDS owns a transposed
    # copy internally, and we've already extracted counts_subset above.
    del counts_df
    gc.collect()

    # DESeq2 sub-steps with granular progress
    # Now shifted by 2 (validation + gene_filter)
    def _deseq_progress(current: int, total: int, i18n_key: str):
        """Translate DESeq2 sub-step progress into pipeline-level progress."""
        nonlocal step
        # DESeq2 steps occupy positions 2..7 in the pipeline
        step = 2 + current
        if progress_callback:
            progress_callback(step, TOTAL_STEPS, i18n_key)

    dds, step_timings = run_deseq2(dds, progress_callback=_deseq_progress)

    # Wald test (step 8)
    step = 8
    def _contrast_progress(current: int, total: int, i18n_key: str):
        nonlocal step
        if progress_callback:
            progress_callback(step, TOTAL_STEPS, i18n_key)

    t0_wald = time.monotonic()
    results_df, _ = compute_contrast(
        dds, metadata_df, condition_col, reference_level, alpha, test_level,
        progress_callback=_contrast_progress,
    )
    step_timings["wald_test"] = time.monotonic() - t0_wald

    # ── Step 3: Visualization ─────────────────────────────────────────
    step = 9
    _report("progress.visualizations")

    t0_viz = time.monotonic()
    figures = {}

    # Volcano plot
    shrinkage = results_df.attrs.get("shrinkage_applied", False)
    volcano_df = prepare_volcano_data(
        results_df, alpha, log2fc_threshold, min_base_mean=min_base_mean,
    )
    figures["volcano"] = create_volcano_plot(
        volcano_df, alpha, log2fc_threshold, test_level, reference_level,
        shrinkage_applied=shrinkage, legend_loc=legend_loc,
    )

    # PCA plot — pass dds for VST if available, top-variance gene selection
    try:
        pca_df = prepare_pca_data(
            counts_subset, metadata_subset, condition_col, dds=dds,
            batch_col=batch_col,
        )
    except Exception:
        # Fallback without dds if VST extraction fails
        pca_df = prepare_pca_data(
            counts_subset, metadata_subset, condition_col, dds=None,
            batch_col=batch_col,
        )
    figures["pca"] = create_pca_plot(pca_df, condition_col, legend_loc=legend_loc)

    # Free the dds object — no longer needed after PCA (which uses VST)
    del dds
    gc.collect()

    # MA plot
    ma_df = prepare_ma_data(results_df, alpha, log2fc_threshold)
    figures["ma"] = create_ma_plot(ma_df, test_level, reference_level, legend_loc=legend_loc)

    # Heatmap (only samples from the two conditions)
    heatmap_df = prepare_heatmap_data(
        counts_subset, metadata_subset, results_df, condition_col
    )
    figures["heatmap"] = create_heatmap(heatmap_df)
    step_timings["visualizations"] = time.monotonic() - t0_viz

    # Free visualization intermediates
    del counts_subset, metadata_subset
    gc.collect()

    # Attach filter stats and timing info to the figures dict for the UI
    figures["filter_stats"] = filter_stats
    total_elapsed = time.monotonic() - pipeline_t0
    figures["timing"] = {
        "step_timings": step_timings,
        "total_elapsed": total_elapsed,
        "time_estimate": time_est,
    }

    # Done
    step = 10
    _report("progress.done")

    return results_df, figures
