"""
pages/2_🔬_scRNA-seq.py — scRNA-seq Analysis dashboard.

Complete single-cell RNA-seq analysis pipeline powered by Scanpy.
Workflow: Load -> QC -> Filter -> Doublets -> Normalize -> HVG ->
          Scale -> PCA -> Neighbors -> UMAP -> Leiden -> Markers.
"""

import base64
import gc
import re
import tempfile
from pathlib import Path

import streamlit as st
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from i18n import t, tf
from config import SCRNA_CONFIG, UPLOAD_LIMITS_LOCAL
from runtime_utils import is_running_locally, apply_local_upload_limit
from visualization import LEGEND_LOCATIONS

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plt_close_all():
    """Close all matplotlib figures to free memory."""
    plt.close("all")

# ─────────────────────────────────────────────────────────────────────
# Page configuration
# ─────────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="scRNA-seq Analysis",
    layout="wide",
    page_icon="🔬",
)

# ─────────────────────────────────────────────────────────────────────
# Sidebar: language selector
# ─────────────────────────────────────────────────────────────────────
if "language" not in st.session_state:
    st.session_state.language = "en"

with st.sidebar:
    _lang_options = {"English": "en", "Español": "es"}
    _current_label = "English" if st.session_state.language == "en" else "Español"
    _selected_label = st.selectbox(
        "🌐",
        options=list(_lang_options.keys()),
        index=list(_lang_options.keys()).index(_current_label),
        key="lang_selector_scrna",
        label_visibility="collapsed",
    )
    if _lang_options[_selected_label] != st.session_state.language:
        st.session_state.language = _lang_options[_selected_label]
        st.rerun()

# ─────────────────────────────────────────────────────────────────────
# Localhost gate + dependency check
# ─────────────────────────────────────────────────────────────────────
apply_local_upload_limit()

# Check dependencies FIRST — on cloud, scanpy is not installed, so this
# catches the cloud case even if is_running_locally() has a false positive.
_scrna_deps_available = True
try:
    import scanpy  # noqa: F401 — quick availability check
except ImportError:
    _scrna_deps_available = False

if not is_running_locally() or not _scrna_deps_available:
    st.title(f"🔬 {t('scrna.title')}")
    st.markdown("---")
    st.warning(t("scrna.remote_blocked"))
    st.stop()

# ─────────────────────────────────────────────────────────────────────
# Title with animation
# ─────────────────────────────────────────────────────────────────────
_cell_path = Path(__file__).resolve().parent.parent / "assets" / "Cell.MP4"
if _cell_path.exists():
    _cell_b64 = base64.b64encode(_cell_path.read_bytes()).decode()
    st.markdown(
        f"""
        <div style="display:flex; align-items:center; gap:12px; margin-bottom:0.25em;">
            <video autoplay loop muted playsinline
                   style="height:60px; width:auto; border-radius:6px;">
                <source src="data:video/mp4;base64,{_cell_b64}" type="video/mp4">
            </video>
            <h1 style="margin:0; padding:0; font-size:2em;">{t('scrna.title')}</h1>
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.title(f"🔬 {t('scrna.title')}")

st.markdown(t("scrna.subtitle"))
st.markdown("---")

# ─────────────────────────────────────────────────────────────────────
# AnnData Glossary
# ─────────────────────────────────────────────────────────────────────
with st.expander(f"📖 {t('scrna.glossary_title')}", expanded=False):
    st.markdown(t("scrna.glossary_desc"))
    st.markdown("")

    _glossary_data = {
        t("scrna.glossary_col_element"): ["`adata.X`", "`adata.obs`", "`adata.var`", "`adata.obsm`", "`adata.uns`", "`adata.layers`"],
        t("scrna.glossary_col_stands_for"): [
            t("scrna.glossary_x_name"),
            t("scrna.glossary_obs_name"),
            t("scrna.glossary_var_name"),
            t("scrna.glossary_obsm_name"),
            t("scrna.glossary_uns_name"),
            t("scrna.glossary_layers_name"),
        ],
        t("scrna.glossary_col_contains"): [
            t("scrna.glossary_x_desc"),
            t("scrna.glossary_obs_desc"),
            t("scrna.glossary_var_desc"),
            t("scrna.glossary_obsm_desc"),
            t("scrna.glossary_uns_desc"),
            t("scrna.glossary_layers_desc"),
        ],
    }

    _glossary_df = pd.DataFrame(_glossary_data)
    st.table(_glossary_df)

st.markdown("---")

# ─────────────────────────────────────────────────────────────────────
# 10x File Integrator (matrix + features + barcodes → H5AD)
# ─────────────────────────────────────────────────────────────────────
with st.expander(f"🧩 {t('scrna.integrator_title')}", expanded=False):
    st.markdown(t("scrna.integrator_desc"))

    int_col1, int_col2, int_col3 = st.columns(3)

    with int_col1:
        mtx_file = st.file_uploader(
            t("scrna.integrator_matrix_label"),
            type=["mtx", "gz"],
            help=t("scrna.integrator_matrix_help"),
            key="scrna_mtx_upload",
        )
    with int_col2:
        feat_file = st.file_uploader(
            t("scrna.integrator_features_label"),
            type=["tsv", "gz"],
            help=t("scrna.integrator_features_help"),
            key="scrna_feat_upload",
        )
    with int_col3:
        bc_file = st.file_uploader(
            t("scrna.integrator_barcodes_label"),
            type=["tsv", "gz"],
            help=t("scrna.integrator_barcodes_help"),
            key="scrna_bc_upload",
        )

    # ── Optional metadata file ──────────────────────────────────────
    st.markdown("---")
    st.markdown(f"**{t('scrna.integrator_meta_title')}**")
    st.caption(t("scrna.integrator_meta_desc"))

    meta_file = st.file_uploader(
        t("scrna.integrator_meta_label"),
        type=["csv", "tsv", "gz"],
        help=t("scrna.integrator_meta_help"),
        key="scrna_meta_upload",
    )

    # If metadata is uploaded, let the user preview the auto-detection
    # or override it — but only after the 3 core files are also present.
    _meta_force_target = None
    if meta_file is not None:
        st.caption(t("scrna.integrator_meta_auto_note"))
        _meta_override = st.radio(
            t("scrna.integrator_meta_target_label"),
            options=["auto", "obs", "var"],
            format_func=lambda x: {
                "auto": t("scrna.integrator_meta_target_auto"),
                "obs": t("scrna.integrator_meta_target_obs"),
                "var": t("scrna.integrator_meta_target_var"),
            }[x],
            index=0,
            horizontal=True,
            key="scrna_meta_target",
        )
        if _meta_override != "auto":
            _meta_force_target = _meta_override

    # ── Build H5AD button ───────────────────────────────────────────
    if mtx_file and feat_file and bc_file:
        integrate_btn = st.button(
            t("scrna.integrator_run"),
            type="primary",
            key="scrna_integrate_btn",
        )
        if integrate_btn:
            from scrna_pipeline import integrate_10x_files

            with st.spinner(t("scrna.integrator_running")):
                try:
                    h5ad_bytes = integrate_10x_files(
                        mtx_file, feat_file, bc_file,
                        metadata_file=meta_file,
                        metadata_force_target=_meta_force_target,
                    )
                    st.success(t("scrna.integrator_success"))

                    # If metadata was included, show what happened
                    if meta_file is not None:
                        from scrna_pipeline import (
                            _read_metadata_file,
                            detect_metadata_target,
                        )
                        import anndata as _ad
                        import tempfile as _tmpf
                        # Quick re-read to show detection info
                        _tmp = _tmpf.NamedTemporaryFile(
                            suffix=".h5ad", delete=False
                        )
                        _tmp.write(h5ad_bytes)
                        _tmp.close()
                        _adata_check = _ad.read_h5ad(_tmp.name)
                        Path(_tmp.name).unlink(missing_ok=True)
                        meta_file.seek(0)
                        _meta_df = _read_metadata_file(meta_file)
                        _det = detect_metadata_target(_adata_check, _meta_df)
                        _target_used = (
                            _meta_force_target
                            if _meta_force_target
                            else _det["target"]
                        )
                        if _target_used == "obs":
                            st.info(tf(
                                "scrna.integrator_meta_merged_obs",
                                overlap=_det["obs_overlap"],
                                pct=_det["obs_pct"],
                                total=_det["meta_ids"],
                            ))
                        elif _target_used == "var":
                            st.info(tf(
                                "scrna.integrator_meta_merged_var",
                                overlap=_det["var_overlap"],
                                pct=_det["var_pct"],
                                total=_det["meta_ids"],
                            ))
                        else:
                            st.warning(t("scrna.integrator_meta_no_match"))

                    _raw_name = st.text_input(
                        t("scrna.integrator_filename_label"),
                        value="integrated",
                        help=t("scrna.integrator_filename_help"),
                        key="integrator_filename",
                    )
                    _clean = re.sub(r'[^\w\-. ]+', '_', _raw_name.strip())
                    _fname = (_clean if _clean else "integrated") + ".h5ad"

                    st.download_button(
                        t("scrna.integrator_download"),
                        data=h5ad_bytes,
                        file_name=_fname,
                        mime="application/octet-stream",
                        key="dl_integrated_h5ad",
                    )
                except Exception as e:
                    st.error(tf("scrna.integrator_error", error=str(e)))

# ─────────────────────────────────────────────────────────────────────
# Ambient RNA Cleaner — SoupX via rpy2
# ─────────────────────────────────────────────────────────────────────

# Cache dependency check in session_state so re-runs don't re-import R
if "_soupx_ok" not in st.session_state:
    from scrna_pipeline import check_soupx_available
    st.session_state["_soupx_status"] = check_soupx_available()
    st.session_state["_soupx_ok"] = st.session_state["_soupx_status"]["soupx"]

with st.expander(f"🧹 {t('scrna.soupx_title')}", expanded=False):
    # Disclaimer: CellBender is better but too heavy
    st.info(t("scrna.soupx_disclaimer"))
    st.markdown(t("scrna.soupx_desc"))

    _soupx_status = st.session_state["_soupx_status"]

    if not _soupx_status["rpy2"]:
        st.warning(t("scrna.soupx_missing_rpy2"))
    elif not _soupx_status["r"]:
        st.warning(t("scrna.soupx_missing_r"))
    elif not _soupx_status["soupx"]:
        st.warning(t("scrna.soupx_missing_pkg"))
    else:
        st.success(t("scrna.soupx_ready"))

    # Always render widgets (even if deps missing — they stay disabled-looking)
    # This prevents Streamlit widget-tree mismatches on re-runs.
    if st.session_state["_soupx_ok"]:
        sx_col1, sx_col2 = st.columns(2)

        with sx_col1:
            sx_raw_file = st.file_uploader(
                t("scrna.soupx_raw_label"),
                type=["h5ad"],
                help=t("scrna.soupx_raw_help"),
                key="scrna_soupx_raw",
            )
        with sx_col2:
            sx_filt_file = st.file_uploader(
                t("scrna.soupx_filt_label"),
                type=["h5ad"],
                help=t("scrna.soupx_filt_help"),
                key="scrna_soupx_filt",
            )

        # Optional: manual contamination fraction
        sx_auto = st.checkbox(
            t("scrna.soupx_auto_label"),
            value=True,
            help=t("scrna.soupx_auto_help"),
            key="scrna_soupx_auto",
        )
        sx_contam = None
        if not sx_auto:
            sx_contam = st.slider(
                t("scrna.soupx_contam_label"),
                min_value=0.01,
                max_value=0.50,
                value=0.10,
                step=0.01,
                help=t("scrna.soupx_contam_help"),
                key="scrna_soupx_contam",
            )

        if sx_raw_file and sx_filt_file:
            soupx_btn = st.button(
                t("scrna.soupx_run"),
                type="primary",
                key="scrna_soupx_btn",
            )
            if soupx_btn:
                from scrna_pipeline import run_soupx

                with st.spinner(t("scrna.soupx_running")):
                    try:
                        cleaned_bytes = run_soupx(
                            raw_h5ad_bytes=sx_raw_file.getvalue(),
                            filtered_h5ad_bytes=sx_filt_file.getvalue(),
                            contamination_fraction=sx_contam,
                        )
                        st.success(t("scrna.soupx_success"))
                        st.download_button(
                            t("scrna.soupx_download"),
                            data=cleaned_bytes,
                            file_name="soupx_cleaned.h5ad",
                            mime="application/octet-stream",
                            key="dl_soupx_h5ad",
                        )
                    except Exception as e:
                        st.error(tf("scrna.soupx_error", error=str(e)))

# ─────────────────────────────────────────────────────────────────────
# H5AD File Upload (analysis input)
# ─────────────────────────────────────────────────────────────────────
st.subheader(t("scrna.upload_header"))

uploaded_file = st.file_uploader(
    t("scrna.upload_label"),
    type=SCRNA_CONFIG["file_extensions"],
    help=t("scrna.upload_help"),
    key="scrna_file_upload",
)

if uploaded_file is None:
    st.stop()

# ── Per-page soft file-size check (scRNA: 5 GB) ─────────────────────
if is_running_locally():
    _scrna_limit_mb = UPLOAD_LIMITS_LOCAL["scrna_mb"]
    _size_mb = len(uploaded_file.getvalue()) / (1024 * 1024)
    if _size_mb > _scrna_limit_mb:
        st.error(
            tf("upload.file_too_large",
               file=uploaded_file.name,
               size_mb=f"{_size_mb:,.0f}",
               limit_mb=f"{_scrna_limit_mb:,}")
        )
        st.stop()

# ─────────────────────────────────────────────────────────────────────
# Load Data (H5AD only)  — backed (disk) mode for RAM-efficient preview
# ─────────────────────────────────────────────────────────────────────

class _AnnDataPreview:
    """Lightweight proxy that mimics AnnData for metadata-only preview.

    Holds obs, var, and shape so the sidebar / preview expander can work
    without loading the full expression matrix into RAM.  The actual
    expression matrix is only loaded when the user clicks "Run Pipeline"
    (via ``_load_h5ad``).
    """
    def __init__(self, obs, var, shape, _unused=None):
        self.obs = obs
        self.var = var
        self._shape = shape

    @property
    def n_obs(self):
        return self._shape[0]

    @property
    def n_vars(self):
        return self._shape[1]


@st.cache_data(show_spinner=False)
def _load_h5ad(file_bytes: bytes):
    """Load H5AD data with caching based on file content."""
    from scrna_pipeline import load_h5ad as _load

    class _FakeUpload:
        def __init__(self, data):
            self._data = data
        def getvalue(self):
            return self._data

    return _load(_FakeUpload(file_bytes))


@st.cache_data(show_spinner=False)
def _load_h5ad_backed(file_bytes: bytes):
    """Load only metadata (backed/disk mode) for RAM-efficient preview.

    Returns (obs, var, shape) — the expression matrix is NOT loaded,
    saving several GB of RAM for large datasets during parameter tuning.
    """
    from scrna_pipeline import load_h5ad_backed as _load_backed
    obs, var, shape, tmp_path = _load_backed(file_bytes)
    # Clean up the temporary file immediately — we only needed metadata
    Path(tmp_path).unlink(missing_ok=True)
    return obs, var, shape


# Use backed mode for the initial load — only reads obs/var metadata,
# saving several GB of RAM for large datasets during parameter tuning.
# The full matrix is loaded only when the user clicks "Run Pipeline".
try:
    with st.spinner(t("scrna.running")):
        _file_bytes = uploaded_file.getvalue()
        _obs, _var, _shape = _load_h5ad_backed(_file_bytes)
        adata_raw = _AnnDataPreview(_obs, _var, _shape, None)
except Exception as e:
    st.error(tf("scrna.error_loading", error=str(e)))
    st.stop()

st.success(
    tf("scrna.data_loaded", n_cells=adata_raw.n_obs, n_genes=adata_raw.n_vars)
)

with st.expander(t("scrna.preview_title"), expanded=False):
    st.caption(tf("scrna.cells_x_genes", n_cells=adata_raw.n_obs, n_genes=adata_raw.n_vars))

    # ── Cell metadata (obs) ──
    st.markdown(f"**{t('scrna.preview_obs_title')}**")
    if adata_raw.obs.shape[1] > 0:
        st.dataframe(adata_raw.obs.head(10), width="stretch")
        st.caption(tf("scrna.preview_obs_caption",
                       rows=adata_raw.obs.shape[0],
                       cols=adata_raw.obs.shape[1]))
    else:
        st.info(t("scrna.preview_obs_empty"))

    # ── Gene metadata (var) ──
    st.markdown(f"**{t('scrna.preview_var_title')}**")
    if adata_raw.var.shape[1] > 0:
        st.dataframe(adata_raw.var.head(10), width="stretch")
        st.caption(tf("scrna.preview_var_caption",
                       rows=adata_raw.var.shape[0],
                       cols=adata_raw.var.shape[1]))
    else:
        st.info(t("scrna.preview_var_empty"))

st.markdown("---")

# ─────────────────────────────────────────────────────────────────────
# QC Parameters (sidebar)
# ─────────────────────────────────────────────────────────────────────
with st.sidebar:
    st.subheader(t("scrna.qc_header"))
    st.caption(t("scrna.qc_description"))

    _qc_defaults = SCRNA_CONFIG["qc_defaults"]

    min_genes = st.number_input(
        t("scrna.min_genes_label"),
        min_value=0, max_value=10000,
        value=_qc_defaults["min_genes"],
        step=50,
        help=t("scrna.min_genes_help"),
        key="scrna_min_genes",
    )
    max_genes = st.number_input(
        t("scrna.max_genes_label"),
        min_value=500, max_value=50000,
        value=_qc_defaults["max_genes"],
        step=500,
        help=t("scrna.max_genes_help"),
        key="scrna_max_genes",
    )
    min_counts = st.number_input(
        t("scrna.min_counts_label"),
        min_value=0, max_value=100000,
        value=_qc_defaults["min_counts"],
        step=100,
        help=t("scrna.min_counts_help"),
        key="scrna_min_counts",
    )
    max_counts = st.number_input(
        t("scrna.max_counts_label"),
        min_value=1000, max_value=500000,
        value=_qc_defaults["max_counts"],
        step=5000,
        help=t("scrna.max_counts_help"),
        key="scrna_max_counts",
    )
    max_pct_mt = st.slider(
        t("scrna.max_pct_mt_label"),
        min_value=0.0, max_value=100.0,
        value=_qc_defaults["max_pct_mt"],
        step=1.0,
        help=t("scrna.max_pct_mt_help"),
        key="scrna_max_pct_mt",
    )
    min_cells = st.number_input(
        t("scrna.min_cells_label"),
        min_value=1, max_value=100,
        value=_qc_defaults["min_cells"],
        step=1,
        help=t("scrna.min_cells_help"),
        key="scrna_min_cells",
    )

    st.markdown("---")
    st.subheader(t("scrna.params_header"))

    n_hvg = st.number_input(
        t("scrna.n_hvg_label"),
        min_value=500, max_value=10000,
        value=SCRNA_CONFIG["hvg"]["n_top_genes"],
        step=500,
        help=t("scrna.n_hvg_help"),
        key="scrna_n_hvg",
    )
    n_pcs = st.slider(
        t("scrna.n_pcs_label"),
        min_value=5, max_value=50,
        value=SCRNA_CONFIG["neighbors"]["n_pcs"],
        step=5,
        help=t("scrna.n_pcs_help"),
        key="scrna_n_pcs",
    )
    n_neighbors = st.slider(
        t("scrna.n_neighbors_label"),
        min_value=5, max_value=100,
        value=SCRNA_CONFIG["neighbors"]["n_neighbors"],
        step=5,
        help=t("scrna.n_neighbors_help"),
        key="scrna_n_neighbors",
    )
    umap_min_dist = st.slider(
        t("scrna.umap_min_dist_label"),
        min_value=0.0, max_value=1.0,
        value=SCRNA_CONFIG["umap"]["min_dist"],
        step=0.05,
        help=t("scrna.umap_min_dist_help"),
        key="scrna_umap_min_dist",
    )
    leiden_res = st.slider(
        t("scrna.leiden_res_label"),
        min_value=0.05, max_value=3.0,
        value=SCRNA_CONFIG["leiden"]["resolution"],
        step=0.05,
        help=t("scrna.leiden_res_help"),
        key="scrna_leiden_res",
    )
    de_method = st.selectbox(
        t("scrna.de_method_label"),
        options=["wilcoxon", "t-test", "t-test_overestim_var", "logreg"],
        index=0,
        help=t("scrna.de_method_help"),
        key="scrna_de_method",
    )
    n_marker_genes = st.number_input(
        t("scrna.n_marker_genes_label"),
        min_value=5, max_value=100,
        value=SCRNA_CONFIG["rank_genes"]["n_genes"],
        step=5,
        help=t("scrna.n_marker_genes_help"),
        key="scrna_n_marker_genes",
    )
    enable_doublet = st.checkbox(
        t("scrna.doublet_checkbox"),
        value=True,
        help=t("scrna.doublet_help"),
        key="scrna_enable_doublet",
    )

    # ── Batch Effect Correction ───────────────────────────────────────
    st.markdown("---")
    st.subheader(t("scrna.batch_header"))

    enable_batch_correction = st.checkbox(
        t("scrna.batch_checkbox"),
        value=False,
        help=t("scrna.batch_help"),
        key="scrna_enable_batch_correction",
    )

    batch_col = None
    if enable_batch_correction:
        _batch_candidates = []
        for _col in adata_raw.obs.columns:
            n_unique = adata_raw.obs[_col].nunique()
            if 2 <= n_unique < adata_raw.n_obs:
                _batch_candidates.append(_col)

        if _batch_candidates:
            batch_col = st.selectbox(
                t("scrna.batch_col_label"),
                options=_batch_candidates,
                index=0,
                help=t("scrna.batch_col_help"),
                key="scrna_batch_col",
            )
            if batch_col:
                _bc = adata_raw.obs[batch_col].value_counts()
                st.caption(
                    tf("scrna.batch_preview",
                       n_batches=len(_bc), col=batch_col)
                )
        else:
            st.warning(t("scrna.batch_no_columns"))

    # ── Annotation Columns ─────────────────────────────────────────────
    st.markdown("---")
    st.subheader(t("scrna.annotation_header"))

    # Build list of categorical obs columns from the uploaded data
    _cat_cols = ["leiden"]
    for _col in adata_raw.obs.columns:
        if hasattr(adata_raw.obs[_col], "cat") or adata_raw.obs[_col].dtype == "object":
            if _col not in _cat_cols:
                _cat_cols.append(_col)

    cell_type_col = st.selectbox(
        f"🏷️ {t('scrna.celltype_col_label')}",
        options=[t("scrna.option_none")] + _cat_cols,
        index=0,
        help=t("scrna.celltype_col_help"),
        key="scrna_celltype_col",
    )

    marker_groupby = st.selectbox(
        f"📊 {t('scrna.marker_groupby_label')}",
        options=_cat_cols,
        index=0,
        help=t("scrna.marker_groupby_help"),
        key="scrna_marker_groupby",
    )

    # ── Visualization Settings ─────────────────────────────────────────
    st.markdown("---")
    st.subheader(t("scrna.viz_header"))

    scrna_legend_loc = st.selectbox(
        f"📍 {t('scrna.legend_position_label')}",
        options=LEGEND_LOCATIONS,
        index=0,
        help=t("scrna.legend_position_help"),
        key="scrna_legend_loc",
    )

# ─────────────────────────────────────────────────────────────────────
# Run Pipeline
# ─────────────────────────────────────────────────────────────────────
run_button = st.button(
    t("scrna.run_button"),
    type="primary",
    width="stretch",
    key="scrna_run_btn",
)

if run_button:
    from scrna_pipeline import run_scrna_pipeline, PIPELINE_STEPS, annotate_qc

    # Free previous results from session_state to avoid accumulating
    # multiple AnnData copies in memory across repeated runs.
    for _sk in ("scrna_adata", "scrna_adata_qc_preview", "scrna_params"):
        if _sk in st.session_state:
            del st.session_state[_sk]
    gc.collect()

    # Build params dict from sidebar inputs
    params = {
        "qc": {
            "min_genes": min_genes,
            "max_genes": max_genes,
            "min_counts": min_counts,
            "max_counts": max_counts,
            "max_pct_mt": max_pct_mt,
            "min_cells": min_cells,
        },
        "n_top_genes": n_hvg,
        "n_pcs": n_pcs,
        "n_neighbors": n_neighbors,
        "umap_min_dist": umap_min_dist,
        "leiden_resolution": leiden_res,
        "de_method": de_method,
        "n_marker_genes": n_marker_genes,
        "remove_doublets": enable_doublet,
        "enable_batch_correction": enable_batch_correction,
        "batch_key": batch_col,
        "marker_groupby": marker_groupby,
    }

    # ── Full data load (backed mode → full AnnData) ─────────────────
    # Until now only metadata (obs/var) was in RAM.  Now we need the
    # expression matrix for the actual pipeline run.
    progress_bar = st.progress(0.0, text=t("scrna.running"))
    status_text = st.empty()

    if isinstance(adata_raw, _AnnDataPreview):
        _adata_full = _load_h5ad(uploaded_file.getvalue())
    else:
        _adata_full = adata_raw  # already fully loaded (fallback)

    # Pre-annotate QC on raw data for violin plots BEFORE filtering.
    _adata_for_qc = _adata_full.copy()
    _adata_for_qc = annotate_qc(_adata_for_qc)
    st.session_state["scrna_adata_qc_preview"] = _adata_for_qc
    del _adata_for_qc
    gc.collect()

    def _progress_callback(step_idx, total, step_name):
        frac = (step_idx + 1) / total
        step_key = f"scrna.step.{step_name}"
        msg = t(step_key)
        progress_bar.progress(frac, text=msg)
        status_text.caption(tf("scrna.step_caption", step=step_idx + 1, total=total, msg=msg))

    try:
        adata_result = run_scrna_pipeline(
            _adata_full.copy(),
            params=params,
            progress_callback=_progress_callback,
        )

        if adata_result.n_obs == 0:
            st.error(t("scrna.error_no_cells"))
            st.stop()

        progress_bar.progress(1.0, text=t("scrna.step.done"))
        status_text.empty()
        st.session_state["scrna_adata"] = adata_result
        st.session_state["scrna_params"] = params
        # Free the full raw data — no longer needed after pipeline
        del _adata_full
        gc.collect()

    except Exception as e:
        _err_msg = str(e)
        # LAPACK / scipy raises LinAlgError with a *bytes* message like
        # b'reciprocal condition number ...' when a matrix is nearly singular.
        # Convert to a human-readable message instead of showing raw bytes.
        if "reciprocal condition number" in _err_msg:
            st.error(t("scrna.error_lapack"))
            st.info(t("scrna.error_lapack_hint"))
        else:
            st.error(tf("scrna.error_pipeline", error=_err_msg))
        st.stop()

# ─────────────────────────────────────────────────────────────────────
# Results Dashboard
# ─────────────────────────────────────────────────────────────────────
if "scrna_adata" not in st.session_state:
    st.info(t("scrna.qc_description"))
    st.stop()

adata = st.session_state["scrna_adata"]

st.markdown("---")
st.subheader(t("scrna.results_header"))

# ── Summary Metrics ─────────────────────────────────────────────────
col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric(t("scrna.stats_cells"), f"{adata.n_obs:,}")
with col2:
    st.metric(t("scrna.stats_genes"), f"{adata.n_vars:,}")
with col3:
    n_clusters = adata.uns.get("leiden_stats", {}).get("n_clusters", "—")
    st.metric(t("scrna.stats_clusters"), n_clusters)
with col4:
    n_hvg_val = adata.uns.get("hvg_stats", {}).get("n_hvg", "—")
    st.metric(t("scrna.stats_hvg"), f"{n_hvg_val:,}" if isinstance(n_hvg_val, int) else n_hvg_val)

# ── Pipeline Stats ──────────────────────────────────────────────────
with st.expander(t("scrna.stats_title"), expanded=False):
    if "filtering_stats" in adata.uns:
        fs = adata.uns["filtering_stats"]
        st.markdown(tf("scrna.filter_stats", **fs))

    if "doublet_stats" in adata.uns:
        ds = adata.uns["doublet_stats"]
        st.markdown(tf("scrna.doublet_stats", **ds))

    if "hvg_stats" in adata.uns:
        hs = adata.uns["hvg_stats"]
        st.markdown(tf("scrna.hvg_stats", **hs))

    if "harmony_stats" in adata.uns:
        hms = adata.uns["harmony_stats"]
        st.markdown(tf("scrna.harmony_stats", **hms))

    if "leiden_stats" in adata.uns:
        ls = adata.uns["leiden_stats"]
        st.markdown(tf("scrna.leiden_stats", **ls))

st.markdown("---")

# ─────────────────────────────────────────────────────────────────────
# Visualization Tabs
# ─────────────────────────────────────────────────────────────────────
try:
    from scrna_visualization import (
        plot_qc_violin,
        plot_qc_scatter,
        plot_hvg,
        plot_pca_variance,
        plot_pca_embedding,
        plot_umap,
        plot_marker_genes_ranking,
        plot_dotplot,
        plot_marker_heatmap,
        _fig_to_bytes,
    )
except ImportError:
    st.error("Visualization dependencies are not available.")
    st.stop()

tab_qc, tab_hvg, tab_pca, tab_umap, tab_markers = st.tabs([
    t("scrna.tab_qc"),
    t("scrna.tab_hvg"),
    t("scrna.tab_pca"),
    t("scrna.tab_umap"),
    t("scrna.tab_markers"),
])

# ── QC Tab ──────────────────────────────────────────────────────────
with tab_qc:
    # Use pre-filter QC data if available
    adata_qc = st.session_state.get("scrna_adata_qc_preview", adata)

    fig_violin = plot_qc_violin(adata_qc)
    st.pyplot(fig_violin, width="stretch")

    dl_col1, dl_col2 = st.columns(2)
    with dl_col1:
        st.download_button(
            t("scrna.download_png"),
            data=_fig_to_bytes(fig_violin, "png"),
            file_name="qc_violin.png",
            mime="image/png",
            key="dl_qc_violin_png",
        )
    with dl_col2:
        st.download_button(
            t("scrna.download_svg"),
            data=_fig_to_bytes(fig_violin, "svg"),
            file_name="qc_violin.svg",
            mime="image/svg+xml",
            key="dl_qc_violin_svg",
        )

    st.markdown("---")

    fig_scatter = plot_qc_scatter(adata_qc)
    st.pyplot(fig_scatter, width="stretch")

    dl_sc1, dl_sc2 = st.columns(2)
    with dl_sc1:
        st.download_button(
            f"{t('scrna.download_png')} {t('scrna.dl_suffix_scatter')}",
            data=_fig_to_bytes(fig_scatter, "png"),
            file_name="qc_scatter.png",
            mime="image/png",
            key="dl_qc_scatter_png",
        )
    with dl_sc2:
        st.download_button(
            f"{t('scrna.download_svg')} {t('scrna.dl_suffix_scatter')}",
            data=_fig_to_bytes(fig_scatter, "svg"),
            file_name="qc_scatter.svg",
            mime="image/svg+xml",
            key="dl_qc_scatter_svg",
        )

    plt_close_all()

# ── HVG Tab ─────────────────────────────────────────────────────────
with tab_hvg:
    fig_hvg = plot_hvg(adata)
    st.pyplot(fig_hvg, width="stretch")

    dl_col1, dl_col2 = st.columns(2)
    with dl_col1:
        st.download_button(
            t("scrna.download_png"),
            data=_fig_to_bytes(fig_hvg, "png"),
            file_name="hvg.png",
            mime="image/png",
            key="dl_hvg_png",
        )
    with dl_col2:
        st.download_button(
            t("scrna.download_svg"),
            data=_fig_to_bytes(fig_hvg, "svg"),
            file_name="hvg.svg",
            mime="image/svg+xml",
            key="dl_hvg_svg",
        )

    plt_close_all()

# ── PCA Tab ─────────────────────────────────────────────────────────
with tab_pca:
    col_pca1, col_pca2 = st.columns(2)

    with col_pca1:
        fig_pca_var = plot_pca_variance(adata)
        st.pyplot(fig_pca_var, width="stretch")

    with col_pca2:
        _legend_loc = st.session_state.get("scrna_legend_loc", "best")
        _ct = st.session_state.get("scrna_celltype_col", t("scrna.option_none"))
        _pca_color = _ct if _ct != t("scrna.option_none") and _ct in adata.obs.columns else "leiden"
        fig_pca_emb = plot_pca_embedding(adata, color=_pca_color, legend_loc=_legend_loc)
        st.pyplot(fig_pca_emb, width="stretch")

    dl_col1, dl_col2, dl_col3, dl_col4 = st.columns(4)
    with dl_col1:
        st.download_button(
            f"{t('scrna.download_png')} {t('scrna.dl_suffix_elbow')}",
            data=_fig_to_bytes(fig_pca_var, "png"),
            file_name="pca_elbow.png",
            mime="image/png",
            key="dl_pca_elbow_png",
        )
    with dl_col2:
        st.download_button(
            f"{t('scrna.download_svg')} {t('scrna.dl_suffix_elbow')}",
            data=_fig_to_bytes(fig_pca_var, "svg"),
            file_name="pca_elbow.svg",
            mime="image/svg+xml",
            key="dl_pca_elbow_svg",
        )
    with dl_col3:
        st.download_button(
            f"{t('scrna.download_png')} {t('scrna.dl_suffix_pca')}",
            data=_fig_to_bytes(fig_pca_emb, "png"),
            file_name="pca_embedding.png",
            mime="image/png",
            key="dl_pca_emb_png",
        )
    with dl_col4:
        st.download_button(
            f"{t('scrna.download_svg')} {t('scrna.dl_suffix_pca')}",
            data=_fig_to_bytes(fig_pca_emb, "svg"),
            file_name="pca_embedding.svg",
            mime="image/svg+xml",
            key="dl_pca_emb_svg",
        )

    plt_close_all()

# ── UMAP Tab ────────────────────────────────────────────────────────
with tab_umap:
    # Color selector — prioritize cell type column if set
    _ct = st.session_state.get("scrna_celltype_col", t("scrna.option_none"))
    if _ct != t("scrna.option_none") and _ct in adata.obs.columns:
        obs_options = [_ct, "leiden"]
    else:
        obs_options = ["leiden"]
    for col in adata.obs.columns:
        if col not in obs_options and col not in ("predicted_doublet", "doublet_score"):
            obs_options.append(col)

    umap_color = st.selectbox(
        t("scrna.umap_color_label"),
        options=obs_options,
        index=0,
        help=t("scrna.umap_color_help"),
        key="scrna_umap_color",
    )

    _legend_loc = st.session_state.get("scrna_legend_loc", "best")
    fig_umap = plot_umap(adata, color=umap_color, legend_loc=_legend_loc)
    st.pyplot(fig_umap, width="stretch")

    dl_col1, dl_col2 = st.columns(2)
    with dl_col1:
        st.download_button(
            t("scrna.download_png"),
            data=_fig_to_bytes(fig_umap, "png"),
            file_name=f"umap_{umap_color}.png",
            mime="image/png",
            key="dl_umap_png",
        )
    with dl_col2:
        st.download_button(
            t("scrna.download_svg"),
            data=_fig_to_bytes(fig_umap, "svg"),
            file_name=f"umap_{umap_color}.svg",
            mime="image/svg+xml",
            key="dl_umap_svg",
        )

    st.markdown("---")

    # Gene expression on UMAP
    gene_input = st.text_input(
        t("scrna.gene_search_label"),
        placeholder=t("scrna.gene_search_placeholder"),
        help=t("scrna.gene_search_help"),
        key="scrna_gene_search",
    )

    if gene_input:
        gene_name = gene_input.strip()
        # Search in adata.raw (all genes) first, then fall back to adata (HVGs only)
        _gene_found = (
            (adata.raw is not None and gene_name in adata.raw.var_names)
            or gene_name in adata.var_names
        )
        if _gene_found:
            fig_gene = plot_umap(adata, color=gene_name, title=tf("plot.umap_title_prefix", color=gene_name), legend_loc=_legend_loc)
            st.pyplot(fig_gene, width="stretch")
            dl_g1, dl_g2 = st.columns(2)
            with dl_g1:
                st.download_button(
                    f"{t('scrna.download_png')} ({gene_name})",
                    data=_fig_to_bytes(fig_gene, "png"),
                    file_name=f"umap_{gene_name}.png",
                    mime="image/png",
                    key="dl_umap_gene_png",
                )
            with dl_g2:
                st.download_button(
                    f"{t('scrna.download_svg')} ({gene_name})",
                    data=_fig_to_bytes(fig_gene, "svg"),
                    file_name=f"umap_{gene_name}.svg",
                    mime="image/svg+xml",
                    key="dl_umap_gene_svg",
                )
        else:
            st.warning(tf("scrna.gene_not_found", gene=gene_name))

    plt_close_all()

# ── Marker Genes Tab ────────────────────────────────────────────────
with tab_markers:
    if "rank_genes_groups" in adata.uns:
        # Recover the groupby used during marker gene analysis
        _groupby = adata.uns["rank_genes_groups"]["params"]["groupby"]

        fig_ranking = plot_marker_genes_ranking(adata, n_genes=5, groupby=_groupby)
        st.pyplot(fig_ranking, width="stretch")

        dl_m1, dl_m2 = st.columns(2)
        with dl_m1:
            st.download_button(
                f"{t('scrna.download_png')} {t('scrna.dl_suffix_markers')}",
                data=_fig_to_bytes(fig_ranking, "png"),
                file_name="marker_genes_ranking.png",
                mime="image/png",
                key="dl_markers_png",
            )
        with dl_m2:
            st.download_button(
                f"{t('scrna.download_svg')} {t('scrna.dl_suffix_markers')}",
                data=_fig_to_bytes(fig_ranking, "svg"),
                file_name="marker_genes_ranking.svg",
                mime="image/svg+xml",
                key="dl_markers_svg",
            )

        st.markdown("---")

        # Dotplot — top marker genes per cluster
        fig_dotplot = plot_dotplot(adata, n_genes=5, groupby=_groupby)
        st.pyplot(fig_dotplot, width="stretch")

        dl_d1, dl_d2 = st.columns(2)
        with dl_d1:
            st.download_button(
                f"{t('scrna.download_png')} {t('scrna.dl_suffix_dotplot')}",
                data=_fig_to_bytes(fig_dotplot, "png"),
                file_name="marker_dotplot.png",
                mime="image/png",
                key="dl_dotplot_png",
            )
        with dl_d2:
            st.download_button(
                f"{t('scrna.download_svg')} {t('scrna.dl_suffix_dotplot')}",
                data=_fig_to_bytes(fig_dotplot, "svg"),
                file_name="marker_dotplot.svg",
                mime="image/svg+xml",
                key="dl_dotplot_svg",
            )

        plt_close_all()

        st.markdown("---")

        # Heatmap
        fig_heatmap = plot_marker_heatmap(adata, n_genes=5, groupby=_groupby)
        st.pyplot(fig_heatmap, width="stretch")

        dl_h1, dl_h2 = st.columns(2)
        with dl_h1:
            st.download_button(
                f"{t('scrna.download_png')} {t('scrna.dl_suffix_heatmap')}",
                data=_fig_to_bytes(fig_heatmap, "png"),
                file_name="marker_heatmap.png",
                mime="image/png",
                key="dl_heatmap_png",
            )
        with dl_h2:
            st.download_button(
                f"{t('scrna.download_svg')} {t('scrna.dl_suffix_heatmap')}",
                data=_fig_to_bytes(fig_heatmap, "svg"),
                file_name="marker_heatmap.svg",
                mime="image/svg+xml",
                key="dl_heatmap_svg",
            )

        st.markdown("---")

        # Marker genes table per group
        from scrna_pipeline import get_marker_genes_df

        cluster_col1, cluster_col2 = st.columns([1, 3])
        with cluster_col1:
            selected_cluster = st.selectbox(
                t("scrna.marker_cluster_label"),
                options=list(adata.obs[_groupby].cat.categories),
                key="scrna_marker_cluster",
            )
            n_show = st.number_input(
                t("scrna.marker_n_genes_label"),
                min_value=5, max_value=100, value=25, step=5,
                key="scrna_marker_n_show",
            )
        with cluster_col2:
            df_markers = get_marker_genes_df(adata, group=selected_cluster, n_genes=n_show)
            st.dataframe(df_markers, width="stretch", height=400)

        # Download all markers
        df_all_markers = get_marker_genes_df(adata, group=None, n_genes=n_marker_genes)
        csv_markers = df_all_markers.to_csv(index=False)
        st.download_button(
            t("scrna.download_markers_csv"),
            data=csv_markers,
            file_name="marker_genes_all_clusters.csv",
            mime="text/csv",
            key="dl_all_markers_csv",
        )
    else:
        st.info(t("scrna.no_markers_info"))

    plt_close_all()

# ─────────────────────────────────────────────────────────────────────
# Downloads Section
# ─────────────────────────────────────────────────────────────────────
st.markdown("---")
st.subheader(t("scrna.download_header"))

dl_col1, dl_col2, dl_col3 = st.columns(3)

# Download H5AD
with dl_col1:
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp_h5ad:
        adata.write_h5ad(tmp_h5ad.name)
        h5ad_bytes = Path(tmp_h5ad.name).read_bytes()
        Path(tmp_h5ad.name).unlink(missing_ok=True)

    st.download_button(
        t("scrna.download_h5ad"),
        data=h5ad_bytes,
        file_name="scrna_analysis.h5ad",
        mime="application/octet-stream",
        help=t("scrna.download_h5ad_help"),
        key="dl_h5ad",
    )

# Download cell metadata
with dl_col2:
    csv_obs = adata.obs.to_csv()
    st.download_button(
        t("scrna.download_obs_csv"),
        data=csv_obs,
        file_name="cell_metadata.csv",
        mime="text/csv",
        key="dl_obs_csv",
    )

# Download marker genes
with dl_col3:
    if "rank_genes_groups" in adata.uns:
        from scrna_pipeline import get_marker_genes_df
        df_all = get_marker_genes_df(adata, group=None, n_genes=25)
        csv_all = df_all.to_csv(index=False)
        st.download_button(
            t("scrna.download_markers_csv"),
            data=csv_all,
            file_name="marker_genes.csv",
            mime="text/csv",
            key="dl_markers_csv_bottom",
        )
