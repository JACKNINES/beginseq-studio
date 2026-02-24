"""
engine/audit.py -- Scientific audit log for reproducibility.

Captures parameters, library versions, data dimensions, step timings,
seeds, and results summaries for both Bulk RNA-seq and scRNA-seq
pipelines.  No Streamlit dependency.

Functions
---------
get_library_versions()
    Return a dict of library name -> version string.

build_bulk_audit(pipeline)
    Build a JSON-serializable audit log from a completed DESeq2 pipeline.

build_scrna_audit(adata, params, elapsed_seconds)
    Build a JSON-serializable audit log from a completed scRNA-seq pipeline.

format_audit_text(audit_dict)
    Format an audit dict as human-readable text for lab notebooks.
"""

from __future__ import annotations

import datetime
import platform
import sys
from dataclasses import asdict
from typing import Any, TYPE_CHECKING

if TYPE_CHECKING:
    pass


# ══════════════════════════════════════════════════════════════════════
# Library versions
# ══════════════════════════════════════════════════════════════════════

def get_library_versions() -> dict[str, str]:
    """Return ``{library: version}`` for all relevant scientific packages.

    Libraries that are not installed return ``"not installed"``.
    """
    libs: dict[str, str] = {}
    for name in (
        "pydeseq2",
        "scanpy",
        "anndata",
        "pandas",
        "numpy",
        "scipy",
        "matplotlib",
        "sklearn",
        "harmonypy",
    ):
        try:
            mod = __import__(name)
            libs[name] = getattr(mod, "__version__", "unknown")
        except ImportError:
            libs[name] = "not installed"
    return libs


# ══════════════════════════════════════════════════════════════════════
# JSON-safe serialiser
# ══════════════════════════════════════════════════════════════════════

def _safe_serialize(obj: Any) -> Any:
    """Recursively convert *obj* to JSON-safe Python primitives.

    Handles numpy scalars, numpy arrays, dataclass dicts, pandas
    Timestamps, and other common scientific-Python types.
    """
    import numpy as np

    if obj is None or isinstance(obj, (str, bool)):
        return obj
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, (int, float)):
        return obj
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, dict):
        return {str(k): _safe_serialize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_safe_serialize(v) for v in obj]
    if isinstance(obj, datetime.datetime):
        return obj.isoformat()
    # Fallback: stringify
    return str(obj)


# ══════════════════════════════════════════════════════════════════════
# Bulk RNA-seq (DESeq2) audit builder
# ══════════════════════════════════════════════════════════════════════

def build_bulk_audit(pipeline) -> dict:
    """Build a JSON-serializable audit log from a completed DESeq2 pipeline.

    Parameters
    ----------
    pipeline : DifferentialExpressionPipeline
        A pipeline instance whose ``.run()`` method has completed.

    Returns
    -------
    dict
        Complete audit log, safe for ``json.dumps()``.
    """
    p = pipeline.params
    results_df = pipeline.results_df

    # ── Input data dimensions ─────────────────────────────────────
    filter_stats = pipeline.filter_stats or {}
    n_genes_before = filter_stats.get("n_before", "unknown")
    n_genes_after = filter_stats.get(
        "n_after",
        results_df.shape[0] if results_df is not None else "unknown",
    )
    time_est = pipeline.time_estimate or {}
    n_samples = time_est.get("n_samples", "unknown")

    # ── Results summary ───────────────────────────────────────────
    results_summary: dict[str, Any] = {}
    if results_df is not None:
        sig_mask = results_df["padj"] < p.alpha
        fc_mask = results_df["log2FoldChange"].abs() > p.log2fc_threshold
        both = sig_mask & fc_mask
        results_summary = {
            "n_genes_tested": int(results_df.shape[0]),
            "n_significant_padj": int(sig_mask.sum()),
            "n_significant_padj_and_lfc": int(both.sum()),
            "n_up": int((both & (results_df["log2FoldChange"] > 0)).sum()),
            "n_down": int((both & (results_df["log2FoldChange"] < 0)).sum()),
            "shrinkage_applied": bool(
                results_df.attrs.get("shrinkage_applied", False)
            ),
        }

    # ── Assemble ──────────────────────────────────────────────────
    audit = {
        "beginseq_studio": {
            "pipeline": "bulk_rnaseq",
            "timestamp": datetime.datetime.now(
                datetime.timezone.utc
            ).isoformat(),
            "pipeline_class": "DifferentialExpressionPipeline",
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "libraries": get_library_versions(),
        },
        "input_data": {
            "n_samples": n_samples,
            "n_genes_raw": n_genes_before,
            "n_genes_after_filter": n_genes_after,
        },
        "parameters": asdict(p),
        "gene_filtering": filter_stats,
        "execution": {
            "steps_completed": list(pipeline._step_log),
            "step_timings_seconds": dict(pipeline.step_timings),
            "total_seconds": pipeline._total_elapsed,
            "time_estimate": time_est,
        },
        "validation": (
            pipeline.validation_report.to_dict()
            if pipeline.validation_report
            else None
        ),
        "results_summary": results_summary,
    }

    return _safe_serialize(audit)


# ══════════════════════════════════════════════════════════════════════
# scRNA-seq audit builder
# ══════════════════════════════════════════════════════════════════════

def build_scrna_audit(
    adata,
    params: dict,
    elapsed_seconds: float | None = None,
) -> dict:
    """Build a JSON-serializable audit log from a completed scRNA-seq pipeline.

    Parameters
    ----------
    adata : AnnData
        The AnnData object after ``run_scrna_pipeline()``.
    params : dict
        The ``params`` dict that was passed to ``run_scrna_pipeline()``.
    elapsed_seconds : float, optional
        Total wall-clock time of the pipeline run.

    Returns
    -------
    dict
        Complete audit log, safe for ``json.dumps()``.
    """
    uns = adata.uns

    # ── Collect per-step stats already stored in adata.uns ────────
    step_stats: dict[str, Any] = {}
    for key in (
        "filtering_stats",
        "doublet_stats",
        "hvg_stats",
        "harmony_stats",
        "leiden_stats",
    ):
        if key in uns:
            step_stats[key] = dict(uns[key])

    # ── Results summary ───────────────────────────────────────────
    n_clusters = 0
    if "leiden_stats" in uns:
        n_clusters = uns["leiden_stats"].get("n_clusters", 0)

    results_summary = {
        "n_cells_final": int(adata.n_obs),
        "n_genes_final": int(adata.n_vars),
        "n_clusters": int(n_clusters),
        "has_marker_genes": "rank_genes_groups" in uns,
    }

    # ── Seeds ─────────────────────────────────────────────────────
    seeds = {
        "umap_random_state": 42,
    }

    # ── Assemble ──────────────────────────────────────────────────
    audit = {
        "beginseq_studio": {
            "pipeline": "scrna_seq",
            "timestamp": datetime.datetime.now(
                datetime.timezone.utc
            ).isoformat(),
            "pipeline_function": "run_scrna_pipeline",
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "libraries": get_library_versions(),
        },
        "input_data": {
            "n_cells_before_filter": (
                uns.get("filtering_stats", {}).get("cells_before", "unknown")
            ),
            "n_genes_before_filter": (
                uns.get("filtering_stats", {}).get("genes_before", "unknown")
            ),
        },
        "parameters": params,
        "seeds": seeds,
        "execution": {
            "total_seconds": elapsed_seconds,
            "step_stats": step_stats,
        },
        "results_summary": results_summary,
    }

    return _safe_serialize(audit)


# ══════════════════════════════════════════════════════════════════════
# Human-readable text formatter
# ══════════════════════════════════════════════════════════════════════

def format_audit_text(audit: dict) -> str:
    """Format an audit dict as a human-readable text report.

    Suitable for pasting into a methods section, lab notebook, or
    supplementary materials.

    Parameters
    ----------
    audit : dict
        Output of ``build_bulk_audit()`` or ``build_scrna_audit()``.

    Returns
    -------
    str
        Multi-line plain-text summary.
    """
    lines: list[str] = []
    sep = "=" * 60

    lines.append(sep)
    lines.append("BeginSeq Studio -- Scientific Audit Log")
    lines.append(sep)
    lines.append("")

    # ── Header ────────────────────────────────────────────────────
    meta = audit.get("beginseq_studio", {})
    lines.append(f"Pipeline : {meta.get('pipeline', 'unknown')}")
    lines.append(f"Timestamp: {meta.get('timestamp', 'unknown')}")
    lines.append("")

    # ── Environment ───────────────────────────────────────────────
    lines.append("--- Environment ---")
    env = audit.get("environment", {})
    for k, v in env.items():
        if isinstance(v, dict):
            lines.append(f"  {k}:")
            for sk, sv in v.items():
                lines.append(f"    {sk}: {sv}")
        else:
            lines.append(f"  {k}: {v}")
    lines.append("")

    # ── Input data ────────────────────────────────────────────────
    lines.append("--- Input Data ---")
    for k, v in audit.get("input_data", {}).items():
        lines.append(f"  {k}: {v}")
    lines.append("")

    # ── Parameters ────────────────────────────────────────────────
    lines.append("--- Parameters ---")
    params = audit.get("parameters", {})
    for k, v in params.items():
        if isinstance(v, dict):
            lines.append(f"  {k}:")
            for sk, sv in v.items():
                lines.append(f"    {sk}: {sv}")
        else:
            lines.append(f"  {k}: {v}")
    lines.append("")

    # ── Seeds (scRNA only) ────────────────────────────────────────
    if "seeds" in audit:
        lines.append("--- Seeds ---")
        for k, v in audit["seeds"].items():
            lines.append(f"  {k}: {v}")
        lines.append("")

    # ── Execution ─────────────────────────────────────────────────
    lines.append("--- Execution ---")
    exe = audit.get("execution", {})
    total = exe.get("total_seconds")
    if total is not None:
        lines.append(f"  Total time: {total:.1f}s")
    timings = exe.get("step_timings_seconds", {})
    if timings:
        lines.append("  Step timings:")
        for step, secs in timings.items():
            lines.append(f"    {step}: {secs:.2f}s")
    step_stats = exe.get("step_stats", {})
    if step_stats:
        lines.append("  Step statistics:")
        for step_name, stats in step_stats.items():
            lines.append(f"    {step_name}:")
            if isinstance(stats, dict):
                for sk, sv in stats.items():
                    lines.append(f"      {sk}: {sv}")
    lines.append("")

    # ── Gene filtering (bulk only) ────────────────────────────────
    gf = audit.get("gene_filtering", {})
    if gf:
        lines.append("--- Gene Filtering ---")
        for k, v in gf.items():
            lines.append(f"  {k}: {v}")
        lines.append("")

    # ── Results summary ───────────────────────────────────────────
    lines.append("--- Results Summary ---")
    for k, v in audit.get("results_summary", {}).items():
        lines.append(f"  {k}: {v}")
    lines.append("")

    lines.append(sep)
    return "\n".join(lines)
