"""
engine/pipeline.py — Class-based orchestrator for DESeq2 differential expression.

Encapsulates the complete pipeline state, allowing step-by-step execution,
fluent chaining, reproducibility (pickle), and isolated unit testing.

Usage
-----
Fluent chaining (full pipeline)::

    pipeline = (
        DifferentialExpressionPipeline(counts_df, metadata_df)
        .configure(reference_level="control", test_level="treated")
        .validate()
        .filter_genes()
        .fit()
        .compute_results()
        .generate_plots()
    )
    results_df = pipeline.results_df
    figures = pipeline.figures

Step-by-step (for testing / inspection)::

    pipeline = DifferentialExpressionPipeline(counts_df, metadata_df)
    pipeline.configure(reference_level="control")
    pipeline.validate()
    # inspect pipeline.counts_df, pipeline.metadata_df
    pipeline.filter_genes()
    # inspect pipeline.filter_stats
    pipeline.fit()
    pipeline.compute_results()
    pipeline.generate_plots()

Backward compatible — ``run_deseq2_pipeline()`` delegates here internally.
"""

from __future__ import annotations

import gc
import time
from dataclasses import dataclass, field
from typing import Callable

import pandas as pd

from engine.config import (
    DESEQ2_DEFAULTS,
    GENE_FILTER_DEFAULTS,
)
from engine.validation import (
    validate_counts_df,
    validate_metadata_df,
    normalize_indices,
    align_samples,
    validate_condition_levels,
    filter_low_expression_genes,
    downcast_counts,
    # Strict validation
    ValidationReport,
    Severity,
    check_duplicate_samples,
    check_duplicate_gene_ids,
    check_condition_nan,
    check_metadata_types,
    check_counts_are_raw,
    check_class_imbalance,
    check_sparse_data,
    check_sample_overlap,
    check_alignment_quality,
    check_post_filter_sanity,
    check_filter_params,
)
from engine.deseq_runner import (
    build_deseq_dataset,
    run_deseq2,
    compute_contrast,
    N_DESEQ2_STEPS,
)
from engine.visualization import (
    prepare_volcano_data,
    create_volcano_plot,
    prepare_pca_data,
    create_pca_plot,
    prepare_ma_data,
    create_ma_plot,
    prepare_heatmap_data,
    create_heatmap,
)
from engine.analysis import estimate_pipeline_time


# ─────────────────────────────────────────────────────────────────────
# Parameter snapshot
# ─────────────────────────────────────────────────────────────────────

@dataclass
class PipelineParams:
    """Immutable parameter snapshot for reproducibility.

    Stores every tuneable knob of the DESeq2 pipeline so the exact
    configuration can be inspected, compared, or serialised alongside
    the results.
    """

    condition_col: str = DESEQ2_DEFAULTS["condition_col"]
    reference_level: str = DESEQ2_DEFAULTS["reference_level"]
    test_level: str | None = None
    alpha: float = DESEQ2_DEFAULTS["alpha"]
    log2fc_threshold: float = DESEQ2_DEFAULTS["log2fc_threshold"]

    # Gene pre-filtering
    filter_genes: bool = GENE_FILTER_DEFAULTS["enabled"]
    min_total_count: int = GENE_FILTER_DEFAULTS["min_total_count"]
    min_samples_expressing: int = GENE_FILTER_DEFAULTS["min_samples_expressing"]
    min_count_per_sample: int = GENE_FILTER_DEFAULTS["min_count_per_sample"]

    # Visualisation
    min_base_mean: float = DESEQ2_DEFAULTS.get("min_base_mean", 0)
    batch_col: str | None = None
    legend_loc: str = "upper right"


# ─────────────────────────────────────────────────────────────────────
# Pipeline class
# ─────────────────────────────────────────────────────────────────────

class DifferentialExpressionPipeline:
    """Class-based orchestrator for the DESeq2 differential expression pipeline.

    Encapsulates the complete pipeline state so that each stage can be
    executed independently, intermediate results inspected, and the
    full analysis serialised for reproducibility.

    Every mutating method returns ``self`` to allow fluent chaining::

        pipeline.configure(...).validate().filter_genes().fit()
            .compute_results().generate_plots()

    Pipeline stages
    ~~~~~~~~~~~~~~~~
    1. ``validate()``        — input validation + quality report
    2. ``filter_genes()``    — low-expression gene filtering + downcast
    3. ``fit()``             — build DDS + DESeq2 model (expensive)
    4. ``compute_results()`` — Wald test + LFC shrinkage
    5. ``generate_plots()``  — volcano, PCA, MA, heatmap
    6. ``run()``             — executes all stages sequentially

    Key attributes
    ~~~~~~~~~~~~~~~
    params : PipelineParams
        Immutable snapshot of all tuneable pipeline parameters.
    validation_report : ValidationReport | None
        Populated by ``validate()``; contains quality findings.
    results_df : pd.DataFrame | None
        DESeq2 results table (after ``compute_results()``).
    figures : dict
        Matplotlib figures + intermediate DataFrames (after ``generate_plots()``).
    filter_stats : dict
        Gene filtering statistics (after ``filter_genes()``).
    step_timings : dict[str, float]
        Per-step wall-clock times.

    Constants
    ~~~~~~~~~
    TOTAL_STEPS : int
        Number of progress steps reported (= 11).

    Serialisation
    ~~~~~~~~~~~~~~
    The pipeline is picklable.  Non-serialisable objects (callbacks, DDS,
    matplotlib Figures) are excluded via ``__getstate__``/``__setstate__``.
    """

    # Total pipeline steps (must match the legacy function):
    # 1 validation + 1 gene_filter + 6 DESeq2 sub-steps + 1 Wald + 1 viz + 1 done
    TOTAL_STEPS: int = 1 + 1 + N_DESEQ2_STEPS + 1 + 1 + 1  # = 11

    # ── Constructor ────────────────────────────────────────────────

    def __init__(
        self,
        counts_df: pd.DataFrame,
        metadata_df: pd.DataFrame,
    ) -> None:
        """Initialise the pipeline with raw input data.

        Parameters
        ----------
        counts_df : pd.DataFrame
            Gene-by-sample raw count matrix (genes as rows, samples as
            columns).  Stored by reference — not copied.
        metadata_df : pd.DataFrame
            Sample metadata with at least a condition column.  A defensive
            copy is made internally.

        Notes
        -----
        After construction, call ``configure()`` to set analysis
        parameters, then either ``run()`` for the full pipeline or
        individual steps (``validate()``, ``filter_genes()``, etc.).

        Callbacks (``progress_callback``, ``estimate_callback``) are
        set as attributes after construction and are NOT serialised.
        """
        # Store original inputs (defensive copy of metadata only — counts is large)
        self._counts_input: pd.DataFrame = counts_df
        self._metadata_input: pd.DataFrame = metadata_df.copy()

        # Parameters (set via configure())
        self.params: PipelineParams = PipelineParams()

        # ── Intermediate state (populated by pipeline steps) ───────
        self.counts_df: pd.DataFrame | None = None
        self.metadata_df: pd.DataFrame | None = None
        self.test_level: str | None = None
        self.available_test_levels: list[str] | None = None

        self.validation_report: ValidationReport | None = None
        self.filter_stats: dict = {}
        self.time_estimate: dict | None = None

        # Visualisation subsets (2-condition slice for PCA / heatmap)
        self.counts_subset: pd.DataFrame | None = None
        self.metadata_subset: pd.DataFrame | None = None

        self._dds = None  # DeseqDataSet (not picklable)
        self.step_timings: dict[str, float] = {}

        # Final outputs
        self.results_df: pd.DataFrame | None = None
        self.figures: dict = {}

        # Pipeline tracking
        self._step_log: list[str] = []
        self._pipeline_t0: float | None = None
        self._total_elapsed: float | None = None

        # Progress step counter (mirrors legacy `step` local var)
        self._step: int = 0

        # Callbacks (set externally, NOT serialised)
        self.progress_callback: Callable[[int, int, str], None] | None = None
        self.estimate_callback: Callable[[dict], None] | None = None

    # ── Configuration ──────────────────────────────────────────────

    def configure(self, **kwargs) -> DifferentialExpressionPipeline:
        """Set pipeline parameters.  Returns *self* for chaining.

        Parameters are validated against the ``PipelineParams`` dataclass
        fields.  Unknown keys raise ``ValueError``.
        """
        for key, value in kwargs.items():
            if hasattr(self.params, key):
                setattr(self.params, key, value)
            else:
                raise ValueError(
                    f"Unknown parameter: '{key}'.  "
                    f"Valid keys: {[f.name for f in self.params.__dataclass_fields__.values()]}"
                )
        return self

    # ── Progress helpers ───────────────────────────────────────────

    def _report(self, i18n_key: str) -> None:
        """Report progress for the current step, then advance."""
        if self.progress_callback:
            self.progress_callback(self._step, self.TOTAL_STEPS, i18n_key)
        self._step += 1

    def _make_deseq_progress(self) -> Callable:
        """Return a callback that maps DESeq2 sub-steps → pipeline steps 2..7."""
        def callback(current: int, total: int, i18n_key: str) -> None:
            self._step = 2 + current
            if self.progress_callback:
                self.progress_callback(self._step, self.TOTAL_STEPS, i18n_key)
        return callback

    def _make_contrast_progress(self) -> Callable:
        """Return a callback for the Wald-test step (position 8)."""
        def callback(current: int, total: int, i18n_key: str) -> None:
            if self.progress_callback:
                self.progress_callback(self._step, self.TOTAL_STEPS, i18n_key)
        return callback

    # ── Step 1: Validate ───────────────────────────────────────────

    def validate(self) -> DifferentialExpressionPipeline:
        """Defensive validation of inputs.

        Validates counts matrix and metadata structure, normalises
        indices, aligns samples, and resolves condition levels.
        Populates ``self.validation_report`` with findings.

        After this step
        ~~~~~~~~~~~~~~~~
        - ``self.counts_df``  — validated, normalised, aligned counts
        - ``self.metadata_df`` — validated, normalised, aligned metadata
        - ``self.test_level``  — resolved test-condition level
        - ``self.available_test_levels`` — all non-reference levels
        - ``self.validation_report`` — populated ValidationReport
        """
        self._report("progress.validating")

        p = self.params
        counts = self._counts_input
        meta = self._metadata_input

        # ── Initialise validation report ────────────────────────────
        report = ValidationReport()
        self.validation_report = report

        # ── Pre-structural checks (before core validation) ──────────
        dup_samples = check_duplicate_samples(meta)
        if dup_samples:
            report.add(
                "duplicate_samples", Severity.ERROR,
                "strict.duplicate_samples",
                i18n_params={"count": len(dup_samples),
                             "samples": ", ".join(dup_samples[:10])},
                details=dup_samples,
            )

        dup_genes = check_duplicate_gene_ids(counts)
        if dup_genes:
            report.add(
                "duplicate_gene_ids", Severity.ERROR,
                "strict.duplicate_gene_ids",
                i18n_params={"count": len(dup_genes),
                             "genes": ", ".join(dup_genes[:10])},
                details=dup_genes,
            )

        # Fatal errors → bail out early
        if report.has_errors():
            summary = "; ".join(
                f.i18n_key for f in report.findings
                if f.severity == Severity.ERROR
            )
            raise ValueError(
                f"Validation failed with {report.n_errors} error(s): {summary}"
            )

        # ── Core validation (existing logic) ────────────────────────
        validate_counts_df(counts)

        # Condition NaN check (BEFORE validate_metadata_df, which
        # converts NaN → the string "nan" via .astype(str))
        n_nan = check_condition_nan(meta, p.condition_col)
        if n_nan:
            report.add(
                "condition_nan", Severity.ERROR,
                "strict.condition_nan",
                i18n_params={"count": n_nan},
            )
            raise ValueError(
                f"Validation failed: {n_nan} sample(s) have missing "
                f"condition values in column '{p.condition_col}'."
            )

        meta = validate_metadata_df(meta, p.condition_col)

        # Normalize and align
        counts, meta = normalize_indices(counts, meta)

        # ── Alignment quality analysis ──────────────────────────────
        overlap_info = check_sample_overlap(counts, meta)
        check_alignment_quality(overlap_info, report)

        counts, meta = align_samples(counts, meta, allow_subset=True)

        self.available_test_levels = validate_condition_levels(
            meta, p.condition_col, p.reference_level,
        )

        # ── Post-alignment advisory checks ──────────────────────────
        # Check if data looks like raw counts
        raw_check = check_counts_are_raw(counts)
        if raw_check.get("is_suspect", False):
            report.add(
                "normalized_data", Severity.WARNING,
                "strict.normalized_data",
                i18n_params={"reason": raw_check.get("reason", "unknown")},
                details=raw_check,
            )

        # Check class imbalance
        imbalance = check_class_imbalance(meta, p.condition_col)
        if imbalance is not None:
            report.add(
                "class_imbalance", Severity.WARNING,
                "strict.class_imbalance",
                i18n_params={"ratio": imbalance.get("ratio", "?")},
                details=imbalance,
            )

        # Check sparse data
        if check_sparse_data(counts):
            report.add(
                "sparse_data", Severity.INFO,
                "strict.sparse_data",
            )

        # Check metadata types (numeric conditions)
        type_info = check_metadata_types(meta, p.condition_col)
        if type_info.get("is_numeric", False):
            report.add(
                "condition_numeric", Severity.WARNING,
                "strict.condition_numeric",
                i18n_params={"sample_values": type_info.get("sample_values", "")},
                details=type_info,
            )

        # Check filter params if filtering is enabled
        if p.filter_genes:
            check_filter_params(
                p.min_samples_expressing, counts.shape[1], report,
            )

        # Resolve test level
        if p.test_level is None:
            self.test_level = self.available_test_levels[0]
        elif p.test_level not in self.available_test_levels:
            raise ValueError(
                f"The test level '{p.test_level}' is not valid. "
                f"Available options: {self.available_test_levels}"
            )
        else:
            self.test_level = p.test_level

        self.counts_df = counts
        self.metadata_df = meta

        self._step_log.append("validate")
        return self

    # ── Step 2: Filter genes ───────────────────────────────────────

    def filter_genes(self) -> DifferentialExpressionPipeline:
        """Apply low-expression gene filtering and downcast to int32.

        After this step
        ~~~~~~~~~~~~~~~~
        - ``self.counts_df``       — filtered and downcasted
        - ``self.filter_stats``    — dict with n_before, n_after, …
        - ``self.time_estimate``   — updated with actual post-filter dims
        - ``self.counts_subset``   — 2-condition slice for visualisation
        - ``self.metadata_subset`` — matching metadata slice
        """
        self._report("progress.filtering_genes")

        p = self.params
        counts = self.counts_df
        meta = self.metadata_df

        if p.filter_genes:
            n_before = counts.shape[0]
            counts, self.filter_stats = filter_low_expression_genes(
                counts,
                metadata_df=meta,
                condition_col=p.condition_col,
                min_total_count=p.min_total_count,
                min_samples_expressing=p.min_samples_expressing,
                min_count_per_sample=p.min_count_per_sample,
            )
            n_after = counts.shape[0]

            # Post-filter sanity check
            if self.validation_report:
                check_post_filter_sanity(n_before, n_after, self.validation_report)

        # Downcast to int32 to halve memory footprint
        counts = downcast_counts(counts)
        gc.collect()

        # Compute time estimate AFTER gene filtering — now we know
        # the actual dimensions DESeq2 will process.
        n_samples_final = counts.shape[1]
        n_genes_final = counts.shape[0]
        self.time_estimate = estimate_pipeline_time(n_samples_final, n_genes_final)

        if self.estimate_callback:
            self.estimate_callback(self.time_estimate)

        # Extract visualisation subset BEFORE building DDS (so we can
        # free counts_df right after DDS construction).
        selected = [p.reference_level, self.test_level]
        mask = meta[p.condition_col].isin(selected)
        self.metadata_subset = meta[mask].copy()
        self.counts_subset = counts[self.metadata_subset.index].copy()

        self.counts_df = counts

        self._step_log.append("filter_genes")
        return self

    # ── Step 3: Fit model ──────────────────────────────────────────

    def fit(self) -> DifferentialExpressionPipeline:
        """Build DeseqDataSet and run the full DESeq2 model.

        This is the computationally expensive step (size factors,
        dispersion estimation, MAP, LFC, Cook's).

        After this step
        ~~~~~~~~~~~~~~~~
        - ``self._dds``          — fitted DeseqDataSet
        - ``self.step_timings``  — per-step timing dict
        """
        p = self.params

        self._dds = build_deseq_dataset(
            self.counts_df, self.metadata_df,
            p.condition_col, p.reference_level,
        )
        # Free the original genes×samples matrix — DDS owns a transposed
        # copy internally, and counts_subset was already extracted.
        self.counts_df = None
        gc.collect()

        self._dds, self.step_timings = run_deseq2(
            self._dds, progress_callback=self._make_deseq_progress(),
        )

        self._step_log.append("fit")
        return self

    # ── Step 4: Compute results ────────────────────────────────────

    def compute_results(self) -> DifferentialExpressionPipeline:
        """Run the Wald test contrast and (optionally) LFC shrinkage.

        After this step
        ~~~~~~~~~~~~~~~~
        - ``self.results_df`` — full DESeq2 results table
        """
        p = self.params
        self._step = 8

        t0_wald = time.monotonic()
        self.results_df, _ = compute_contrast(
            self._dds,
            self.metadata_df,
            p.condition_col,
            p.reference_level,
            p.alpha,
            self.test_level,
            progress_callback=self._make_contrast_progress(),
        )
        self.step_timings["wald_test"] = time.monotonic() - t0_wald

        self._step_log.append("compute_results")
        return self

    # ── Step 5: Generate plots ─────────────────────────────────────

    def generate_plots(self) -> DifferentialExpressionPipeline:
        """Generate all visualisation figures.

        After this step
        ~~~~~~~~~~~~~~~~
        - ``self.figures`` — dict with ``"volcano"``, ``"pca"``,
          ``"ma"``, ``"heatmap"`` (matplotlib Figures) plus intermediate
          DataFrames (``"pca_df"``, ``"ma_df"``, ``"heatmap_df"``) for
          theme-reactive re-rendering.
        """
        self._step = 9
        self._report("progress.visualizations")

        p = self.params
        t0_viz = time.monotonic()
        figures: dict = {}

        # ── Volcano plot ───────────────────────────────────────────
        shrinkage = self.results_df.attrs.get("shrinkage_applied", False)
        volcano_df = prepare_volcano_data(
            self.results_df, p.alpha, p.log2fc_threshold,
            min_base_mean=p.min_base_mean,
        )
        figures["volcano"] = create_volcano_plot(
            volcano_df, p.alpha, p.log2fc_threshold,
            self.test_level, p.reference_level,
            shrinkage_applied=shrinkage, legend_loc=p.legend_loc,
        )

        # ── PCA plot ───────────────────────────────────────────────
        try:
            pca_df = prepare_pca_data(
                self.counts_subset, self.metadata_subset,
                p.condition_col, dds=self._dds, batch_col=p.batch_col,
            )
        except Exception:
            # Fallback without dds if VST extraction fails
            pca_df = prepare_pca_data(
                self.counts_subset, self.metadata_subset,
                p.condition_col, dds=None, batch_col=p.batch_col,
            )
        figures["pca"] = create_pca_plot(
            pca_df, p.condition_col, legend_loc=p.legend_loc,
        )
        figures["pca_df"] = pca_df

        # Free the dds — no longer needed after PCA (which uses VST)
        self._dds = None
        gc.collect()

        # ── MA plot ────────────────────────────────────────────────
        ma_df = prepare_ma_data(self.results_df, p.alpha, p.log2fc_threshold)
        figures["ma"] = create_ma_plot(
            ma_df, self.test_level, p.reference_level,
            legend_loc=p.legend_loc,
        )
        figures["ma_df"] = ma_df

        # ── Heatmap ───────────────────────────────────────────────
        heatmap_df = prepare_heatmap_data(
            self.counts_subset, self.metadata_subset,
            self.results_df, p.condition_col,
        )
        figures["heatmap"] = create_heatmap(heatmap_df)
        figures["heatmap_df"] = heatmap_df
        self.step_timings["visualizations"] = time.monotonic() - t0_viz

        # Free visualisation intermediates
        self.counts_subset = None
        self.metadata_subset = None
        gc.collect()

        self.figures = figures

        self._step_log.append("generate_plots")
        return self

    # ── Convenience: run all ───────────────────────────────────────

    def run(self) -> DifferentialExpressionPipeline:
        """Run the complete pipeline from validation through plots.

        Equivalent to calling each step in sequence::

            self.validate().filter_genes().fit()
                .compute_results().generate_plots()

        Returns *self* for chaining.
        """
        self._pipeline_t0 = time.monotonic()
        self._step = 0

        (
            self
            .validate()
            .filter_genes()
            .fit()
            .compute_results()
            .generate_plots()
        )

        self._total_elapsed = time.monotonic() - self._pipeline_t0

        # Final "done" report
        self._step = 10
        self._report("progress.done")

        self._step_log.append("run_complete")
        return self

    # ── Output accessors ───────────────────────────────────────────

    def build_audit(self) -> dict:
        """Build a scientific audit log from the completed pipeline state.

        Returns a JSON-serializable dict capturing parameters, library
        versions, step timings, and results summary.  Requires the
        pipeline to have been run (``results_df`` must be populated).
        """
        from engine.audit import build_bulk_audit
        return build_bulk_audit(self)

    def get_results(self) -> tuple[pd.DataFrame, dict]:
        """Return ``(results_df, figures)`` matching the legacy function signature.

        This is the bridge method that ``run_deseq2_pipeline()`` uses.
        The ``figures`` dict is augmented with filter stats and timing
        metadata, exactly as the original function returned.

        Raises
        ------
        RuntimeError
            If the pipeline has not been run yet.
        """
        if self.results_df is None:
            raise RuntimeError(
                "Pipeline has not been run yet. Call .run() or "
                "execute steps individually first."
            )

        figures = dict(self.figures)
        figures["filter_stats"] = self.filter_stats
        figures["timing"] = {
            "step_timings": self.step_timings,
            "total_elapsed": self._total_elapsed,
            "time_estimate": self.time_estimate,
        }
        figures["validation_report"] = (
            self.validation_report.to_dict() if self.validation_report else None
        )
        figures["audit"] = self.build_audit()
        return self.results_df, figures

    # ── Serialisation ──────────────────────────────────────────────

    def __getstate__(self) -> dict:
        """Exclude non-picklable objects for serialisation.

        Removed from pickle:
        - ``progress_callback`` / ``estimate_callback`` (closures)
        - ``_dds`` (DeseqDataSet — contains C extensions)
        - matplotlib Figure objects in ``figures``

        Preserved:
        - ``results_df``, ``filter_stats``, ``step_timings``, ``params``
        - ``_step_log`` (audit trail)
        - Intermediate DataFrames in ``figures`` (``pca_df``, etc.)
        """
        state = self.__dict__.copy()

        # Callbacks are closures — not picklable
        state["progress_callback"] = None
        state["estimate_callback"] = None

        # DeseqDataSet contains C extensions
        state["_dds"] = None

        # Remove matplotlib figures (not reliably picklable)
        if "figures" in state and state["figures"]:
            state["figures"] = {
                k: v for k, v in state["figures"].items()
                if isinstance(v, (pd.DataFrame, dict, list, str, int, float, bool, type(None)))
            }

        return state

    def __setstate__(self, state: dict) -> None:
        """Restore state from pickle, setting non-picklable fields to None."""
        self.__dict__.update(state)
        # Ensure callbacks are always initialised
        if not hasattr(self, "progress_callback"):
            self.progress_callback = None
        if not hasattr(self, "estimate_callback"):
            self.estimate_callback = None
        if not hasattr(self, "_dds"):
            self._dds = None

    # ── Repr ───────────────────────────────────────────────────────

    def __repr__(self) -> str:
        status = self._step_log[-1] if self._step_log else "not started"
        n_genes = (
            self.results_df.shape[0] if self.results_df is not None else "?"
        )
        return (
            f"<DifferentialExpressionPipeline "
            f"status={status!r} "
            f"genes={n_genes} "
            f"params={self.params!r}>"
        )
