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

from typing import Callable

import pandas as pd

from engine.config import DESEQ2_DEFAULTS
from engine.config import GENE_FILTER_DEFAULTS
from engine.config import TIME_ESTIMATION


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

    # ── Delegate to the class-based pipeline ───────────────────────
    # The DifferentialExpressionPipeline class encapsulates the same
    # logic that previously lived here as procedural code.  Delegating
    # preserves backward compatibility (same signature, same return)
    # while gaining reproducibility, step-by-step testing, and pickle.
    from engine.pipeline import DifferentialExpressionPipeline

    pipeline = DifferentialExpressionPipeline(counts_df, metadata_df)
    pipeline.configure(
        condition_col=condition_col,
        reference_level=reference_level,
        test_level=test_level,
        alpha=alpha,
        log2fc_threshold=log2fc_threshold,
        filter_genes=filter_genes,
        min_total_count=min_total_count,
        min_samples_expressing=min_samples_expressing,
        min_count_per_sample=min_count_per_sample,
        min_base_mean=min_base_mean,
        batch_col=batch_col,
        legend_loc=legend_loc,
    )
    pipeline.progress_callback = progress_callback
    pipeline.estimate_callback = estimate_callback

    pipeline.run()

    return pipeline.get_results()
