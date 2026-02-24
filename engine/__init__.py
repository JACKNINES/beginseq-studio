"""
engine -- Pure analysis engine for BeginSeq Studio.

No Streamlit dependency. Importable for headless use, testing,
notebooks, CLI tools, or future web frameworks.

Usage:
    from engine import run_deseq2_pipeline
    from engine.data_io import read_counts_file
    from engine.scrna_pipeline import run_scrna_pipeline
"""

from engine.analysis import (
    run_deseq2_pipeline,
    estimate_pipeline_time,
    estimate_pipeline_time_with_filter,
)
from engine.pipeline import (
    DifferentialExpressionPipeline,
    PipelineParams,
)
from engine.data_io import (
    read_counts_file,
    read_metadata_file,
    results_to_csv,
)
from engine.gdc_client import (
    fetch_tcga_projects,
    fetch_rnaseq_files,
    download_and_extract_batch,
    build_count_matrix,
)
from engine.protocols import ProgressCallback, EstimateCallback, FileData
from engine.config import set_theme, THEME, THEME_PRESETS
from engine.audit import (
    build_bulk_audit,
    build_scrna_audit,
    get_library_versions,
    format_audit_text,
)

# --- Lazy imports for optional heavy dependencies (scanpy) --------
_SCRNA_NAMES = {"run_scrna_pipeline", "load_h5ad", "integrate_10x_files"}

def __getattr__(name: str):
    if name in _SCRNA_NAMES:
        from engine import scrna_pipeline as _sp
        return getattr(_sp, name)
    raise AttributeError(f"module 'engine' has no attribute {name!r}")

__all__ = [
    "run_deseq2_pipeline",
    "estimate_pipeline_time",
    "estimate_pipeline_time_with_filter",
    "DifferentialExpressionPipeline",
    "PipelineParams",
    "read_counts_file",
    "read_metadata_file",
    "results_to_csv",
    "run_scrna_pipeline",
    "load_h5ad",
    "integrate_10x_files",
    "fetch_tcga_projects",
    "fetch_rnaseq_files",
    "download_and_extract_batch",
    "build_count_matrix",
    "ProgressCallback",
    "EstimateCallback",
    "FileData",
    "set_theme",
    "THEME",
    "THEME_PRESETS",
    "build_bulk_audit",
    "build_scrna_audit",
    "get_library_versions",
    "format_audit_text",
]
