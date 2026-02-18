"""
pages/1_🧬_Bulk_RNA-seq.py — Bulk RNA-seq Classic analysis page.

Differential expression analysis from raw count matrices using DESeq2.
This is the main analysis tool, migrated from the original app.py.
"""

import streamlit as st
import sys
import numpy as np
from pathlib import Path

# Ensure project root is in path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import time

from config import STREAMLIT_CONFIG, FILE_CONFIG, DESEQ2_DEFAULTS, GENE_FILTER_DEFAULTS, UPLOAD_LIMITS_LOCAL
from runtime_utils import is_running_locally, apply_local_upload_limit
from data_io import read_counts_file, read_metadata_file, results_to_csv
from validation import (
    normalize_indices,
    validate_metadata_df,
    check_sample_overlap,
    get_condition_levels,
    check_counts_are_raw,
    check_class_imbalance,
    check_sparse_data,
)
from classification import (
    compute_log2_cpm,
    classify_samples,
    merge_classification,
    build_compound_condition,
)
from analysis import run_deseq2_pipeline, estimate_pipeline_time, estimate_pipeline_time_with_filter
from visualization import prepare_volcano_data, create_volcano_plot, create_volcano_highlight, LEGEND_LOCATIONS
from i18n import t, tf

import io


# ─────────────────────────────────────────────────────────────────────
# Cached wrappers — avoid redundant work on Streamlit reruns
# ─────────────────────────────────────────────────────────────────────

@st.cache_data(show_spinner=False)
def _cached_read_counts(file_data: bytes, file_name: str):
    """Read counts matrix (cached by file content hash)."""
    buf = io.BytesIO(file_data)
    buf.name = file_name
    return read_counts_file(buf)


@st.cache_data(show_spinner=False)
def _cached_read_metadata(file_data: bytes, file_name: str):
    """Read metadata file (cached by file content hash)."""
    buf = io.BytesIO(file_data)
    buf.name = file_name
    return read_metadata_file(buf)


@st.cache_data(show_spinner=False)
def _cached_prepare_volcano(results_df, alpha, log2fc_threshold, min_base_mean):
    """Classify genes for volcano plot (cached by parameters)."""
    return prepare_volcano_data(results_df, alpha, log2fc_threshold, min_base_mean)


@st.cache_data(show_spinner=False)
def _cached_create_volcano(volcano_df, alpha, log2fc_threshold,
                           test_level, reference_level, label_genes,
                           title, shrinkage_applied, legend_loc):
    """Render volcano plot (cached by parameters)."""
    return create_volcano_plot(
        volcano_df, alpha, log2fc_threshold, test_level, reference_level,
        label_genes=label_genes, title=title,
        shrinkage_applied=shrinkage_applied, legend_loc=legend_loc,
    )


@st.cache_data(show_spinner=False)
def _cached_create_volcano_highlight(volcano_df, highlight_genes,
                                     alpha, log2fc_threshold,
                                     test_level, reference_level,
                                     legend_loc):
    """Render volcano highlight plot (cached by parameters)."""
    return create_volcano_highlight(
        volcano_df, highlight_genes, alpha, log2fc_threshold,
        test_level, reference_level, legend_loc=legend_loc,
    )


def _format_seconds(seconds: float) -> str:
    """Format seconds into a human-readable string."""
    if seconds < 60:
        return f"{int(seconds)} s"
    elif seconds < 3600:
        mins = int(seconds // 60)
        secs = int(seconds % 60)
        return f"{mins} min {secs} s" if secs > 0 else f"{mins} min"
    else:
        hrs = int(seconds // 3600)
        mins = int((seconds % 3600) // 60)
        return f"{hrs} h {mins} min"


def _fig_download_buttons(fig, base_name: str, key_prefix: str):
    """Render PNG + SVG download buttons side-by-side for a matplotlib figure."""
    import io

    col_png, col_svg = st.columns(2)

    buf_png = io.BytesIO()
    fig.savefig(buf_png, format="png", dpi=300, bbox_inches="tight")
    buf_png.seek(0)
    with col_png:
        st.download_button(
            label=f"📥 {t('viz.download_png')}",
            data=buf_png,
            file_name=f"{base_name}.png",
            mime="image/png",
            key=f"{key_prefix}_png",
        )

    buf_svg = io.BytesIO()
    fig.savefig(buf_svg, format="svg", bbox_inches="tight")
    buf_svg.seek(0)
    with col_svg:
        st.download_button(
            label=f"📥 {t('viz.download_svg')}",
            data=buf_svg,
            file_name=f"{base_name}.svg",
            mime="image/svg+xml",
            key=f"{key_prefix}_svg",
        )


# ─────────────────────────────────────────────────────────────────────
# Page configuration
# ─────────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title=f"{STREAMLIT_CONFIG['page_title']} — Bulk RNA-seq",
    layout=STREAMLIT_CONFIG["layout"],
    page_icon="🧬",
)

# ─────────────────────────────────────────────────────────────────────
# Local-only: raise upload limit
# ─────────────────────────────────────────────────────────────────────
apply_local_upload_limit()

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
        key="lang_selector_bulk",
        label_visibility="collapsed",
    )
    if _lang_options[_selected_label] != st.session_state.language:
        st.session_state.language = _lang_options[_selected_label]
        st.rerun()

# ─────────────────────────────────────────────────────────────────────
# Page title (with helix animation)
# ─────────────────────────────────────────────────────────────────────
import base64, pathlib

_helix_path = pathlib.Path(__file__).resolve().parent.parent / "helix.MP4"
if _helix_path.exists():
    _helix_b64 = base64.b64encode(_helix_path.read_bytes()).decode()
    st.markdown(
        f"""
        <div style="display:flex; align-items:center; gap:12px; margin-bottom:0.25em;">
            <video autoplay loop muted playsinline
                   style="height:60px; width:auto; border-radius:6px;">
                <source src="data:video/mp4;base64,{_helix_b64}" type="video/mp4">
            </video>
            <h1 style="margin:0; padding:0; font-size:2em;">{t('app.title')}</h1>
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.title(f"🧬 {t('app.title')}")
st.markdown(t("app.welcome"))

# ─────────────────────────────────────────────────────────────────────
# File upload
# ─────────────────────────────────────────────────────────────────────
st.header(f"📁 {t('upload.header')}")

counts_file = st.file_uploader(
    t("upload.counts_label"),
    type=FILE_CONFIG["counts_extensions"],
    help=t("upload.counts_help"),
)

with st.expander(f"ℹ️ {t('upload.expander_title')}"):
    st.markdown(t("upload.counts_info"))

metadata_file = st.file_uploader(
    t("upload.metadata_label"),
    type=FILE_CONFIG["metadata_extensions"],
    help=t("upload.metadata_help"),
)

# ── Per-page soft file-size check (Bulk: 2 GB) ──────────────────────
if is_running_locally():
    _bulk_limit_mb = UPLOAD_LIMITS_LOCAL["bulk_rna_mb"]
    for _uf, _label in [(counts_file, t("upload.counts_label")),
                         (metadata_file, t("upload.metadata_label"))]:
        if _uf is not None:
            _size_mb = len(_uf.getvalue()) / (1024 * 1024)
            if _size_mb > _bulk_limit_mb:
                st.error(
                    tf("upload.file_too_large",
                       file=_uf.name,
                       size_mb=f"{_size_mb:,.0f}",
                       limit_mb=f"{_bulk_limit_mb:,}")
                )
                st.stop()

# ─────────────────────────────────────────────────────────────────────
# Analysis execution
# ─────────────────────────────────────────────────────────────────────
if counts_file and metadata_file:

    # Read files using data_io (auto-detects separator)
    # Cached by file content — changing a widget won't re-read the files.
    try:
        counts_df = _cached_read_counts(counts_file.getvalue(), counts_file.name)
        metadata_df_raw = _cached_read_metadata(metadata_file.getvalue(), metadata_file.name)
    except ValueError as e:
        st.error(f"❌ {tf('error.reading_files', error=e)}")
        st.stop()

    # ─────────────────────────────────────────────────────────────────
    # Drop non-numeric columns (e.g. gene_name, gene_type from STAR files)
    # Memory note: .drop(copy=False) avoids copying the underlying data
    # blocks when possible (pandas ≥ 2.0 CoW).  We also del the original
    # reference immediately to free any old block.
    # ─────────────────────────────────────────────────────────────────
    non_numeric_cols = counts_df.select_dtypes(exclude="number").columns.tolist()
    if non_numeric_cols:
        counts_df = counts_df[counts_df.columns.difference(non_numeric_cols)]

    # ─────────────────────────────────────────────────────────────────
    # Filter genes with zero counts across all samples.
    # Memory note: use numpy to compute row sums without pandas overhead.
    # np.asarray avoids copying when the dtype is already numeric.
    # ─────────────────────────────────────────────────────────────────
    n_genes_original = counts_df.shape[0]
    row_sums = np.asarray(counts_df).sum(axis=1)
    keep_mask = row_sums > 0
    if not keep_mask.all():
        counts_df = counts_df.iloc[keep_mask]
    n_genes_filtered = n_genes_original - counts_df.shape[0]

    # Show preview of loaded data
    with st.expander(f"👀 {t('preview.counts_title')}", expanded=False):
        st.dataframe(counts_df.head())
        st.caption(
            tf("preview.counts_caption",
               genes=counts_df.shape[0],
               samples=counts_df.shape[1])
        )
        if n_genes_filtered > 0:
            st.info(
                f"ℹ️ {tf('preview.genes_removed', count=n_genes_filtered)}"
            )

    # ─────────────────────────────────────────────────────────────────
    # Heuristic: Does the matrix appear to be normalized?
    # ─────────────────────────────────────────────────────────────────
    raw_check = check_counts_are_raw(counts_df)
    if raw_check["is_suspect"]:
        st.warning(
            f"⚠️ {tf('warning.normalized_matrix', reason=raw_check['reason'])}"
        )

    with st.expander(f"👀 {t('preview.metadata_title')}", expanded=False):
        st.dataframe(metadata_df_raw.head(10))
        st.caption(
            tf("preview.metadata_caption",
               rows=metadata_df_raw.shape[0],
               cols=metadata_df_raw.shape[1])
        )

    # ─────────────────────────────────────────────────────────────────
    # Metadata column selection
    # ─────────────────────────────────────────────────────────────────
    st.header(f"🔧 {t('config.header')}")

    # Normalize column names for display
    metadata_df_raw.columns = metadata_df_raw.columns.str.strip()
    all_columns = metadata_df_raw.columns.tolist()

    st.markdown(t("config.instructions"))

    col1, col2 = st.columns(2)

    with col1:
        sample_col_default = 0
        for i, col in enumerate(all_columns):
            if col.lower() in ["sample", "samples", "sample_id", "sampleid", "id"]:
                sample_col_default = i
                break

        sample_col = st.selectbox(
            f"📋 {t('config.sample_col_label')}",
            options=all_columns,
            index=sample_col_default,
            help=t("config.sample_col_help"),
        )

    with col2:
        condition_options = [c for c in all_columns if c != sample_col]

        condition_col_default = 0
        for i, col in enumerate(condition_options):
            if col.lower() in ["condition", "group", "treatment", "disease", "status", "type"]:
                condition_col_default = i
                break

        condition_col = st.selectbox(
            f"🏷️ {t('config.condition_col_label')}",
            options=condition_options,
            index=condition_col_default,
            help=t("config.condition_col_help"),
        )

    # Filter columns (optional) — only if there are more than 2 columns
    filter_options = [c for c in all_columns if c not in [sample_col, condition_col]]

    if filter_options:
        st.markdown("---")
        st.subheader(f"🔍 {t('filter.subheader')}")

        use_filters = st.checkbox(
            t("filter.checkbox_label"),
            value=False,
            help=t("filter.checkbox_help"),
        )

        if use_filters:
            # ── Initialize multiple filter state ───────────
            if "sample_filters" not in st.session_state:
                st.session_state.sample_filters = [{"id": 0, "column": None, "values": []}]
            if "filter_counter" not in st.session_state:
                st.session_state.filter_counter = 0

            filters = st.session_state.sample_filters
            working_metadata = metadata_df_raw.copy()
            filters_valid = True

            for idx, filt in enumerate(filters):
                used_cols = [
                    f["column"] for j, f in enumerate(filters)
                    if j != idx and f["column"] is not None
                ]
                available_cols = [c for c in filter_options if c not in used_cols]

                if not available_cols:
                    break

                filt_id = filt["id"]
                st.markdown(f"**{tf('filter.title', n=idx + 1)}**")
                col_sel, col_vals, col_del = st.columns([2, 4, 1])

                with col_sel:
                    default_idx = 0
                    if filt["column"] and filt["column"] in available_cols:
                        default_idx = available_cols.index(filt["column"])

                    selected_col = st.selectbox(
                        t("filter.column_label"),
                        options=available_cols,
                        index=default_idx,
                        key=f"filter_col_{filt_id}",
                        label_visibility="collapsed" if idx > 0 else "visible",
                    )
                    filters[idx]["column"] = selected_col

                with col_vals:
                    unique_vals = (
                        working_metadata[selected_col].dropna().unique().tolist()
                    )
                    selected_vals = st.multiselect(
                        tf("filter.values_label", column=selected_col),
                        options=unique_vals,
                        default=unique_vals,
                        key=f"filter_vals_{filt_id}",
                        help=t("filter.values_help"),
                    )
                    filters[idx]["values"] = selected_vals

                with col_del:
                    st.markdown("<br>", unsafe_allow_html=True)
                    if len(filters) > 1:
                        if st.button("🗑️", key=f"remove_filter_{filt_id}"):
                            filters.pop(idx)
                            st.session_state.sample_filters = filters
                            st.rerun()

                # Apply this filter to working_metadata
                if selected_vals:
                    working_metadata = working_metadata[
                        working_metadata[selected_col].isin(selected_vals)
                    ].copy()
                else:
                    st.warning(
                        f"⚠️ {tf('filter.no_values_warning', n=idx + 1)}"
                    )
                    filters_valid = False
                    break

            # ── Button to add another filter ────────────────────
            remaining_cols = [
                c for c in filter_options
                if c not in [f["column"] for f in filters if f["column"]]
            ]
            if remaining_cols and filters_valid:
                if st.button(f"➕ {t('filter.add_button')}"):
                    st.session_state.filter_counter += 1
                    new_id = st.session_state.filter_counter
                    filters.append({"id": new_id, "column": None, "values": []})
                    st.session_state.sample_filters = filters
                    st.rerun()

            # ── Apply filters and show summary ─────────────────
            if filters_valid and len(working_metadata) > 0:
                metadata_df_raw = working_metadata
                st.info(
                    f"ℹ️ {tf('filter.summary', count=len(metadata_df_raw), n_filters=len(filters))}"
                )
            elif filters_valid:
                st.warning(f"⚠️ {t('filter.all_removed_warning')}")
                st.stop()

            st.session_state.sample_filters = filters

    # ─────────────────────────────────────────────────────────────────
    # Prepare metadata with selected columns
    # ─────────────────────────────────────────────────────────────────
    metadata_df = metadata_df_raw[[sample_col, condition_col]].copy()
    metadata_df.columns = ["sample", "condition"]

    # ─────────────────────────────────────────────────────────────────
    # Pre-validate and normalize data BEFORE showing options.
    # ─────────────────────────────────────────────────────────────────
    try:
        metadata_validated = validate_metadata_df(metadata_df)
        counts_norm, meta_norm = normalize_indices(
            counts_df, metadata_validated
        )
        overlap = check_sample_overlap(counts_norm, meta_norm)
    except ValueError as e:
        st.error(f"❌ {tf('error.validation', error=e)}")
        st.stop()

    allow_subset = False

    if not overlap["match"]:
        st.warning(f"⚠️ {t('mismatch.warning')}")

        col1, col2, col3 = st.columns(3)
        col1.metric(t("mismatch.metric_counts"), overlap["n_counts"])
        col2.metric(t("mismatch.metric_metadata"), overlap["n_metadata"])
        col3.metric(t("mismatch.metric_common"), overlap["n_common"])

        if overlap["n_common"] == 0:
            st.error(f"❌ {t('mismatch.no_common')}")
            st.stop()

        with st.expander(f"🔍 {t('mismatch.details_title')}", expanded=True):
            if overlap["only_in_counts"]:
                st.markdown(
                    tf("mismatch.only_in_counts",
                       count=len(overlap["only_in_counts"]))
                )
                st.code(
                    ", ".join(sorted(overlap["only_in_counts"])),
                    language=None,
                )

            if overlap["only_in_metadata"]:
                st.markdown(
                    tf("mismatch.only_in_metadata",
                       count=len(overlap["only_in_metadata"]))
                )
                st.code(
                    ", ".join(sorted(overlap["only_in_metadata"])),
                    language=None,
                )

        st.info(
            f"ℹ️ {tf('mismatch.proceed_info', n_common=overlap['n_common'])}"
        )

        allow_subset = st.checkbox(
            f"✅ {t('mismatch.proceed_checkbox')}",
            value=False,
        )

        if not allow_subset:
            st.stop()

        common_samples = sorted(overlap["common"])
        counts_norm = counts_norm[common_samples]
        meta_norm = meta_norm.loc[common_samples]

    # ─────────────────────────────────────────────────────────────────
    # Detect condition levels and allow reference selection
    # ─────────────────────────────────────────────────────────────────
    try:
        condition_levels = get_condition_levels(meta_norm)
    except ValueError:
        unique_levels = meta_norm["condition"].unique().tolist()
        st.error(
            f"❌ {tf('error.single_condition', n_levels=len(unique_levels), levels=unique_levels, n_samples=meta_norm.shape[0], first_level=unique_levels[0] if unique_levels else '?')}"
        )
        st.stop()

    st.header(f"⚙️ {t('params.header')}")

    col_ref, col_test = st.columns(2)

    with col_ref:
        reference_level = st.selectbox(
            f"🔵 {t('params.reference_label')}",
            options=condition_levels,
            index=0,
            help=t("params.reference_help"),
        )

    with col_test:
        test_options = [lvl for lvl in condition_levels if lvl != reference_level]

        if len(test_options) == 1:
            test_level = test_options[0]
            st.selectbox(
                f"🔴 {t('params.test_label')}",
                options=test_options,
                index=0,
                disabled=True,
                help=t("params.test_help_single"),
            )
        else:
            test_level = st.selectbox(
                f"🔴 {t('params.test_label')}",
                options=test_options,
                index=0,
                help=t("params.test_help_multi"),
            )

    # Show contrast summary
    st.success(
        f"📐 {tf('params.contrast_summary', test=test_level, ref=reference_level)}"
    )

    # ─────────────────────────────────────────────────────────────────
    # Expression-based sample classification (optional)
    # ─────────────────────────────────────────────────────────────────
    st.markdown("---")
    st.subheader(f"🧪 {t('classify.subheader')}")
    st.markdown(t("classify.description"))

    use_classification = st.checkbox(
        f"✅ {t('classify.checkbox')}",
        value=False,
        help=t("classify.checkbox_help"),
    )

    if use_classification:
        # ── Gene input: text area for pasting gene lists ──────────
        genes_input = st.text_area(
            f"🧬 {t('classify.gene_input_label')}",
            height=100,
            placeholder=t("classify.gene_input_placeholder"),
            help=t("classify.gene_input_help"),
            key="classify_genes_text",
        )

        # Parse gene names from text (comma, newline, space, tab separated)
        # Also strip quotes that are common when copying from code/papers
        import re
        raw_genes = re.split(r'[,\n\t;]+', genes_input) if genes_input.strip() else []
        parsed_genes = []
        for g in raw_genes:
            g = g.strip().strip('"').strip("'").strip()
            if g:
                parsed_genes.append(g)
        # Remove duplicates preserving order
        seen = set()
        parsed_genes = [g for g in parsed_genes if not (g in seen or seen.add(g))]

        # Match against available genes in the count matrix
        available_gene_set = set(counts_norm.index.tolist())
        found_genes = [g for g in parsed_genes if g in available_gene_set]
        not_found_genes = [g for g in parsed_genes if g not in available_gene_set]

        if not parsed_genes:
            st.warning(f"⚠️ {t('classify.no_genes_warning')}")
        else:
            if not_found_genes:
                _nf_display = ", ".join(not_found_genes[:15])
                _nf_msg = tf("classify.genes_not_found_summary",
                             n_found=len(found_genes),
                             n_not_found=len(not_found_genes),
                             not_found_list=_nf_display)
                if len(not_found_genes) > 15:
                    _nf_msg += tf("bulk.more_genes_suffix", extra=len(not_found_genes) - 15)
                st.warning(f"⚠️ {_nf_msg}")

            if found_genes:
                st.success(
                    f"✅ {tf('classify.genes_found_summary', n=len(found_genes))}"
                )

                # Compute log2(CPM+1) once — cache in session state
                if "classify_log_cpm" not in st.session_state or \
                   st.session_state.get("classify_log_cpm_hash") != id(counts_norm):
                    st.session_state.classify_log_cpm = compute_log2_cpm(counts_norm)
                    st.session_state.classify_log_cpm_hash = id(counts_norm)

                log_cpm = st.session_state.classify_log_cpm

                # ── Global threshold & direction (applied to ALL markers) ──
                st.markdown(f"**{t('classify.global_params_title')}**")
                col_thr, col_dir = st.columns(2)
                with col_thr:
                    global_threshold = st.number_input(
                        t("classify.global_threshold_label"),
                        min_value=0.0,
                        max_value=25.0,
                        value=1.0,
                        step=0.1,
                        format="%.2f",
                        help=t("classify.threshold_help"),
                        key="classify_global_thr",
                    )
                with col_dir:
                    dir_options = [
                        t("classify.direction_below"),
                        t("classify.direction_above"),
                    ]
                    direction_label = st.selectbox(
                        t("classify.global_direction_label"),
                        options=dir_options,
                        index=0,
                        help=t("classify.direction_help"),
                        key="classify_global_dir",
                    )
                    global_direction = "below" if direction_label == dir_options[0] else "above"

                # Show expression stats (with z-score of the selected threshold)
                with st.expander(
                    f"📊 {tf('classify.stats_expander', n=len(found_genes))}",
                    expanded=False,
                ):
                    import pandas as _pd
                    import numpy as _np
                    stats_rows = []
                    for gene in found_genes:
                        gv = log_cpm.loc[gene]
                        _mean = float(gv.mean())
                        _sd = float(gv.std())
                        _z = round((global_threshold - _mean) / _sd, 2) if _sd > 0 else _np.nan
                        stats_rows.append({
                            "gene": gene,
                            "min": round(float(gv.min()), 2),
                            "median": round(float(gv.median()), 2),
                            "mean": round(_mean, 2),
                            "sd": round(_sd, 2),
                            "max": round(float(gv.max()), 2),
                            f"z(thr={global_threshold})": _z,
                        })
                    st.dataframe(_pd.DataFrame(stats_rows).set_index("gene"))
                    st.caption(t("classify.zscore_caption"))

                # Build rules: same threshold & direction for all genes
                rules = [
                    {"gene": g, "threshold": global_threshold, "direction": global_direction}
                    for g in found_genes
                ]

                # Labels configuration
                col_pos, col_neg = st.columns(2)
                with col_pos:
                    positive_label = st.text_input(
                        f"✅ {t('classify.positive_label')}",
                        value="positive",
                        help=t("classify.positive_help"),
                        key="classify_pos_label",
                    )
                with col_neg:
                    negative_label = st.text_input(
                        f"❌ {t('classify.negative_label')}",
                        value="negative",
                        help=t("classify.negative_help"),
                        key="classify_neg_label",
                    )

                # Option to keep reference samples unchanged
                keep_reference = st.checkbox(
                    f"🔒 {t('classify.reference_keep_label')}",
                    value=True,
                    help=t("classify.reference_keep_help"),
                    key="classify_keep_ref",
                )

                # Run classification
                try:
                    classification = classify_samples(
                        log_cpm, rules,
                        positive_label=positive_label,
                        negative_label=negative_label,
                    )

                    # Merge classification and build compound conditions
                    class_col_name = f"{positive_label}_status"
                    meta_classified = merge_classification(
                        meta_norm, classification,
                        class_col_name=class_col_name,
                        condition_col="condition",
                    )
                    meta_classified = build_compound_condition(
                        meta_classified,
                        class_col_name=class_col_name,
                        condition_col="condition",
                        reference_level=reference_level if keep_reference else None,
                    )

                    # Show preview with compound condition counts
                    compound_counts = meta_classified["condition"].value_counts()
                    cc_df = compound_counts.reset_index()
                    cc_df.columns = ["condition", "count"]

                    with st.expander(
                        f"👀 {t('classify.preview_title')}", expanded=True
                    ):
                        st.markdown(f"**{t('classify.compound_preview')}**")
                        st.dataframe(cc_df)

                    # Apply classification to meta_norm
                    meta_norm = meta_classified

                    # Re-detect condition levels with the new compound conditions
                    try:
                        condition_levels = get_condition_levels(meta_norm)
                    except ValueError:
                        unique_levels = meta_norm["condition"].unique().tolist()
                        st.error(
                            f"❌ {tf('error.single_condition', n_levels=len(unique_levels), levels=unique_levels, n_samples=meta_norm.shape[0], first_level=unique_levels[0] if unique_levels else '?')}"
                        )
                        st.stop()

                    # Let user re-select reference and test levels
                    st.markdown("---")
                    st.markdown(f"**⚙️ {t('params.header')}** {t('bulk.params_updated_classification')}")

                    col_ref2, col_test2 = st.columns(2)

                    with col_ref2:
                        # Default to the same reference if it still exists
                        ref_idx = 0
                        if reference_level in condition_levels:
                            ref_idx = condition_levels.index(reference_level)

                        reference_level = st.selectbox(
                            f"🔵 {t('params.reference_label')}",
                            options=condition_levels,
                            index=ref_idx,
                            help=t("params.reference_help"),
                            key="ref_level_classified",
                        )

                    with col_test2:
                        test_options_2 = [
                            lvl for lvl in condition_levels
                            if lvl != reference_level
                        ]

                        if len(test_options_2) == 1:
                            test_level = test_options_2[0]
                            st.selectbox(
                                f"🔴 {t('params.test_label')}",
                                options=test_options_2,
                                index=0,
                                disabled=True,
                                help=t("params.test_help_single"),
                                key="test_level_classified",
                            )
                        else:
                            test_level = st.selectbox(
                                f"🔴 {t('params.test_label')}",
                                options=test_options_2,
                                index=0,
                                help=t("params.test_help_multi"),
                                key="test_level_classified",
                            )

                    st.success(
                        f"📐 {tf('params.contrast_summary', test=test_level, ref=reference_level)}"
                    )

                    # Download button for classified metadata (sample + new condition only)
                    _dl_meta = meta_classified[["condition"]].copy()
                    _dl_meta.index.name = "sample"
                    classified_csv = _dl_meta.to_csv()
                    st.download_button(
                        label=f"📥 {tf('classify.download_metadata', fmt='CSV')}",
                        data=classified_csv,
                        file_name="classified_metadata.csv",
                        mime="text/csv",
                        key="download_classified_metadata",
                    )

                except ValueError as e:
                    st.error(f"❌ {tf('classify.error', error=str(e))}")

    # ─────────────────────────────────────────────────────────────────
    # Gene pre-filtering options
    # ─────────────────────────────────────────────────────────────────
    st.markdown("---")
    st.subheader(f"🧹 {t('gene_filter.subheader')}")
    st.markdown(t("gene_filter.description"))

    filter_genes = st.checkbox(
        f"✅ {t('gene_filter.checkbox')}",
        value=GENE_FILTER_DEFAULTS["enabled"],
        help=t("gene_filter.checkbox_help"),
    )

    if filter_genes:
        col_mtc, col_mse, col_mcs = st.columns(3)
        with col_mtc:
            min_total_count = st.number_input(
                t("gene_filter.min_total_label"),
                min_value=0,
                max_value=1000,
                value=GENE_FILTER_DEFAULTS["min_total_count"],
                step=5,
                help=t("gene_filter.min_total_help"),
            )
        with col_mse:
            min_samples_expressing = st.number_input(
                t("gene_filter.min_samples_label"),
                min_value=0,
                max_value=counts_norm.shape[1],
                value=GENE_FILTER_DEFAULTS["min_samples_expressing"],
                step=1,
                help=t("gene_filter.min_samples_help"),
            )
        with col_mcs:
            min_count_per_sample = st.number_input(
                t("gene_filter.min_count_label"),
                min_value=1,
                max_value=100,
                value=GENE_FILTER_DEFAULTS["min_count_per_sample"],
                step=1,
                help=t("gene_filter.min_count_help"),
            )
    else:
        min_total_count = 0
        min_samples_expressing = 0
        min_count_per_sample = 1

    # ─────────────────────────────────────────────────────────────────
    # Subset to selected conditions (reference + test) for estimates
    # and pre-analysis warnings.  The full pipeline also does this
    # internally, but we need accurate numbers for the UI display.
    # ─────────────────────────────────────────────────────────────────
    _selected_mask = meta_norm["condition"].isin([reference_level, test_level])
    _meta_selected = meta_norm[_selected_mask]
    _counts_selected = counts_norm[_meta_selected.index]

    # ─────────────────────────────────────────────────────────────────
    # Show time estimate BEFORE the run button
    # ─────────────────────────────────────────────────────────────────
    n_samples_for_est = _counts_selected.shape[1]
    n_genes_for_est = _counts_selected.shape[0]
    pre_estimate = estimate_pipeline_time_with_filter(
        n_samples_for_est, n_genes_for_est, filter_genes,
    )

    if pre_estimate.get("filter_applied"):
        est_genes_display = pre_estimate["n_genes_estimated"]
        st.info(
            f"⏱️ {tf('time.estimate_detail_filtered', n_samples=n_samples_for_est, n_genes_raw=n_genes_for_est, n_genes_est=est_genes_display, estimate=pre_estimate['total_human'])}"
        )
    else:
        st.info(
            f"⏱️ {tf('time.estimate_detail', n_samples=n_samples_for_est, n_genes=n_genes_for_est, estimate=pre_estimate['total_human'])}"
        )

    # ─────────────────────────────────────────────────────────────────
    # Pre-analysis warnings (shown BEFORE the run button)
    # Only check the two selected conditions (reference + test),
    # not all conditions in the metadata.
    # ─────────────────────────────────────────────────────────────────
    imbalance = check_class_imbalance(_meta_selected, "condition")
    if imbalance:
        st.warning(
            f"⚠️ {tf('warning.class_imbalance', **imbalance)}"
        )

    is_sparse = check_sparse_data(_counts_selected)
    if is_sparse:
        st.warning(
            f"⚠️ {tf('warning.sparse_data')}"
        )

    # ─────────────────────────────────────────────────────────────────
    # Batch correction for PCA (optional)
    # ─────────────────────────────────────────────────────────────────
    st.markdown("---")
    batch_col = None
    _meta_cols = [
        c for c in meta_norm.columns
        if c not in ("condition", "sample", "condition_original")
    ]
    if _meta_cols:
        use_batch = st.checkbox(
            f"🧪 {t('pca.batch_checkbox')}",
            value=False,
            help=t("pca.batch_checkbox_help"),
            key="pca_batch_checkbox",
        )
        if use_batch:
            batch_col = st.selectbox(
                f"📦 {t('pca.batch_select_label')}",
                options=_meta_cols,
                help=t("pca.batch_select_help"),
                key="pca_batch_col",
            )

    # Button to run the analysis
    if st.button(f"🚀 {t('params.run_button')}"):
        progress_bar = st.progress(0, text=t("progress.validating"))
        status_text = st.empty()
        pipeline_start = time.monotonic()

        # ── Mutable ETA state ──────────────────────────────────────
        # These get recalculated when the pipeline reports actual
        # post-filter dimensions via "progress.time_estimate".
        _eta_state = {
            "est_total": pre_estimate["total_seconds"],
            "per_step": dict(pre_estimate["per_step"]),
            "cumulative": {},  # step_index -> cumulative seconds
        }

        def _build_cumulative(per_step: dict, est_total: float) -> dict:
            """Build cumulative time map: step_index -> cumulative s."""
            cum = {0: 2.0, 1: 5.0}  # validation + gene filtering
            running = 5.0
            deseq_keys = ["size_factors", "genewise_disp", "disp_trend",
                          "map_disp", "fit_lfc", "cooks"]
            for i, key in enumerate(deseq_keys):
                running += per_step.get(key, 0.0)
                cum[2 + i] = running
            running += per_step.get("wald_test", 0.0)
            cum[8] = running
            running += per_step.get("visualizations", 0.0)
            cum[9] = running
            cum[10] = est_total
            return cum

        _eta_state["cumulative"] = _build_cumulative(
            _eta_state["per_step"], _eta_state["est_total"],
        )

        def _estimate_callback(real_estimate: dict):
            """Called by the pipeline after gene filtering with actual dimensions."""
            _eta_state["est_total"] = real_estimate["total_seconds"]
            _eta_state["per_step"] = dict(real_estimate["per_step"])
            _eta_state["cumulative"] = _build_cumulative(
                real_estimate["per_step"], real_estimate["total_seconds"],
            )

        def _pipeline_progress(current: int, total: int, i18n_key: str):
            """Update Streamlit progress bar with ETA from pipeline callback."""
            pct = current / total if total > 0 else 0
            step_msg = t(i18n_key) if i18n_key else ""

            # Calculate ETA based on elapsed time and estimated progress
            elapsed = time.monotonic() - pipeline_start
            est_total = _eta_state["est_total"]
            est_at_step = _eta_state["cumulative"].get(current, 0.0)
            remaining = max(0.0, est_total - est_at_step)

            # If we have real elapsed data, adjust the remaining estimate
            if elapsed > 0 and est_at_step > 0:
                # Use ratio of actual vs estimated time to scale remaining
                speed_factor = elapsed / est_at_step
                remaining = max(0.0, (est_total - est_at_step) * speed_factor)

            if remaining > 2 and current < total:
                remaining_str = _format_seconds(remaining)
                label = tf("time.step_with_eta",
                           step_msg=step_msg, remaining=remaining_str)
            else:
                label = step_msg

            progress_bar.progress(min(pct, 1.0), text=label)

        try:
            # Subset to only the two selected conditions BEFORE the
            # pipeline.  This ensures DESeq2 models only the relevant
            # samples — critical for memory on 8 GB machines and
            # consistent with the time estimate shown to the user.
            _run_mask = meta_norm["condition"].isin(
                [reference_level, test_level]
            )
            _meta_run = meta_norm[_run_mask]
            _counts_run = counts_norm[_meta_run.index]

            results_df, figures = run_deseq2_pipeline(
                _counts_run, _meta_run,
                reference_level=reference_level,
                test_level=test_level,
                progress_callback=_pipeline_progress,
                estimate_callback=_estimate_callback,
                filter_genes=filter_genes,
                min_total_count=int(min_total_count),
                min_samples_expressing=int(min_samples_expressing),
                min_count_per_sample=int(min_count_per_sample),
                min_base_mean=DESEQ2_DEFAULTS.get("min_base_mean", 0),
                batch_col=batch_col,
                legend_loc=st.session_state.get("viz_legend_loc", "upper right"),
            )

            total_elapsed = time.monotonic() - pipeline_start
            elapsed_str = _format_seconds(total_elapsed)
            progress_bar.progress(
                1.0,
                text=tf("time.completed_in", elapsed=elapsed_str),
            )
            status_text.empty()

            st.session_state.results_df = results_df
            st.session_state.figures = figures
            st.session_state.analysis_done = True
            st.session_state.reference_level = reference_level
            st.session_state.test_level = test_level

            # Show filter statistics if gene filtering was applied
            fstats = figures.get("filter_stats", {})
            if fstats and fstats.get("n_removed", 0) > 0:
                st.info(
                    f"🧹 {tf('gene_filter.stats_summary', n_before=fstats['n_before'], n_after=fstats['n_after'], n_removed=fstats['n_removed'], min_total=fstats['min_total_count'], min_samples=fstats['min_samples_expressing'], min_count=fstats.get('min_count_per_sample', 1))}"
                )



            st.success(
                f"✅ {tf('params.success', n_genes=len(results_df))}"
            )

            # Show step timing breakdown in an expander
            timing = figures.get("timing", {})
            step_timings = timing.get("step_timings", {})
            if step_timings:
                with st.expander(f"⏱️ {t('time.step_breakdown')}", expanded=False):
                    _step_name_keys = {
                        "size_factors": "time.step_name_size_factors",
                        "genewise_disp": "time.step_name_genewise_disp",
                        "disp_trend": "time.step_name_disp_trend",
                        "map_disp": "time.step_name_map_disp",
                        "fit_lfc": "time.step_name_fit_lfc",
                        "cooks": "time.step_name_cooks",
                        "wald_test": "time.step_name_wald_test",
                        "visualizations": "time.step_name_visualizations",
                    }
                    for step_key, secs in step_timings.items():
                        name_key = _step_name_keys.get(step_key, step_key)
                        name = t(name_key) if name_key.startswith("time.") else step_key
                        st.text(f"  {name}: {secs:.1f} s")
                    st.text(f"  {'─' * 30}")
                    st.text(f"  {tf('bulk.timing_total', elapsed=_format_seconds(total_elapsed))}")

                    # Compare with estimate
                    est_total = timing.get("time_estimate", {}).get("total_seconds", 0)
                    if est_total > 0:
                        diff = total_elapsed - est_total
                        if diff < -5:
                            st.caption(f"🚀 {t('time.faster_than_expected')}")

        except ValueError as e:
            progress_bar.empty()
            st.error(f"❌ {tf('error.validation', error=e)}")
            st.exception(e)
        except Exception as e:
            progress_bar.empty()
            st.error(f"❌ {tf('error.unexpected', error=e)}")
            st.exception(e)

    # ─────────────────────────────────────────────────────────────────
    # Show results (only if the analysis has been run)
    # ─────────────────────────────────────────────────────────────────
    if st.session_state.get("analysis_done", False):
        results_df = st.session_state.results_df
        figures = st.session_state.figures

        # ─────────────────────────────────────────────────────────
        # Post-analysis filtering options
        # ─────────────────────────────────────────────────────────
        st.header(f"🔬 {t('results_filter.header')}")

        min_basemean = results_df["baseMean"].min()
        max_basemean = results_df["baseMean"].max()
        median_basemean = results_df["baseMean"].median()

        st.markdown(
            tf("results_filter.stats",
               min_bm=min_basemean,
               med_bm=median_basemean,
               max_bm=max_basemean)
        )

        _radio_default = t("results_filter.radio_default")
        _radio_custom = t("results_filter.radio_custom")
        _radio_none = t("results_filter.radio_none")

        filter_option = st.radio(
            t("results_filter.radio_label"),
            options=[_radio_default, _radio_custom, _radio_none],
            index=0,
            help=t("results_filter.radio_help"),
        )

        if filter_option == _radio_default:
            basemean_threshold = 10
            padj_threshold = 0.05
            log2fc_threshold = 0.5

            st.info(
                f"📋 {tf('results_filter.default_info', basemean=basemean_threshold, padj=padj_threshold, log2fc=log2fc_threshold)}"
            )

            filtered_df = results_df[
                (results_df["baseMean"] >= basemean_threshold) &
                (results_df["padj"] < padj_threshold) &
                (results_df["log2FoldChange"].abs() > log2fc_threshold)
            ].copy()

        elif filter_option == _radio_custom:
            st.markdown("---")
            col_f1, col_f2, col_f3 = st.columns(3)

            with col_f1:
                basemean_threshold = st.number_input(
                    t("bulk.basemean_filter_label"), min_value=0.0, max_value=float(max_basemean),
                    value=10.0, step=1.0, help=t("results_filter.basemean_help"),
                )
            with col_f2:
                padj_threshold = st.number_input(
                    t("bulk.padj_filter_label"), min_value=0.0001, max_value=1.0,
                    value=0.05, step=0.01, format="%.4f",
                    help=t("results_filter.padj_help"),
                )
            with col_f3:
                log2fc_threshold = st.number_input(
                    t("bulk.log2fc_filter_label"), min_value=0.0, max_value=10.0,
                    value=0.5, step=0.1, help=t("results_filter.log2fc_help"),
                )

            mask = results_df["baseMean"] >= basemean_threshold
            if padj_threshold < 1.0:
                mask = mask & (results_df["padj"] < padj_threshold)
            if log2fc_threshold > 0:
                mask = mask & (results_df["log2FoldChange"].abs() > log2fc_threshold)
            filtered_df = results_df[mask].copy()

        else:  # No filtering
            filtered_df = results_df.copy()
            padj_threshold = 1.0
            log2fc_threshold = 0.0
            st.warning(f"⚠️ {t('results_filter.no_filter_warning')}")

        n_original = len(results_df)
        n_filtered = len(filtered_df)
        n_removed = n_original - n_filtered

        if filter_option != _radio_none:
            st.success(
                f"📊 {tf('results_filter.summary', n_filtered=n_filtered, n_removed=n_removed, n_total=n_original)}"
            )

        # ─────────────────────────────────────────────────────────
        # Results table
        # ─────────────────────────────────────────────────────────
        st.subheader(f"📊 {t('table.subheader')}")

        if len(filtered_df) == 0:
            st.warning(f"⚠️ {t('table.no_genes')}")
        else:
            filtered_df_sorted = filtered_df.sort_values("padj")
            st.dataframe(filtered_df_sorted.head(50))
            st.caption(tf("table.caption", n=len(filtered_df)))

        # ─────────────────────────────────────────────────────────
        # Ranking table
        # ─────────────────────────────────────────────────────────
        st.subheader(f"📈 {t('ranking.subheader')}")
        st.markdown(t("ranking.description"))

        first_col_name = results_df.index.name or "gene"
        col_rank1, col_rank2 = st.columns([2, 1])

        with col_rank1:
            ranking_id_col = st.text_input(
                t("ranking.id_label"), value=first_col_name,
                help=t("ranking.id_help"),
            )
        with col_rank2:
            ranking_top_n = st.number_input(
                t("ranking.top_n_label"), min_value=10, max_value=500,
                value=50, step=10, help=t("ranking.top_n_help"),
            )

        ranking_df = results_df[["stat"]].dropna().sort_values("stat", ascending=False)
        ranking_df = ranking_df.reset_index()
        ranking_df.columns = [ranking_id_col, "stat"]

        st.dataframe(ranking_df.head(int(ranking_top_n)))
        st.caption(tf("ranking.caption", n=int(ranking_top_n), total=len(ranking_df)))

        # ─────────────────────────────────────────────────────────
        # Download
        # ─────────────────────────────────────────────────────────
        st.subheader(f"💾 {t('download.subheader')}")
        col_name, col_format = st.columns([2, 1])

        with col_name:
            file_basename = st.text_input(
                t("download.filename_label"), value="deseq2_results",
                help=t("download.filename_help"),
            )
        with col_format:
            file_format = st.selectbox(
                t("download.format_label"), options=["CSV", "TSV"],
                index=0, help=t("download.format_help"),
            )

        if file_format == "TSV":
            separator, extension, mime_type = "\t", ".tsv", "text/tab-separated-values"
        else:
            separator, extension, mime_type = ",", ".csv", "text/csv"

        filtered_csv = filtered_df.to_csv(sep=separator)
        all_csv = results_df.to_csv(sep=separator)
        ranking_csv = ranking_df.head(50).to_csv(sep=separator, index=False)

        col_dl1, col_dl2, col_dl3 = st.columns(3)
        with col_dl1:
            st.download_button(
                label=f"📥 {tf('download.btn_filtered', fmt=file_format)}",
                data=filtered_csv,
                file_name=f"{file_basename}_filtered{extension}",
                mime=mime_type,
            )
        with col_dl2:
            st.download_button(
                label=f"📥 {tf('download.btn_all', fmt=file_format)}",
                data=all_csv,
                file_name=f"{file_basename}_all{extension}",
                mime=mime_type,
            )
        with col_dl3:
            st.download_button(
                label=f"📥 {tf('download.btn_ranking', fmt=file_format)}",
                data=ranking_csv,
                file_name=f"{file_basename}_ranking{extension}",
                mime=mime_type,
            )

        # ─────────────────────────────────────────────────────────
        # Visualizations
        # ─────────────────────────────────────────────────────────
        st.subheader(f"📈 {t('viz.subheader')}")

        col_viz1, col_viz2, col_viz3 = st.columns(3)
        with col_viz1:
            label_top_genes = st.checkbox(
                f"🏷️ {t('viz.label_checkbox')}",
                value=False,
                help=t("viz.label_help"),
            )
        with col_viz2:
            viz_min_base_mean = st.number_input(
                f"🔬 {t('viz.basemean_filter_label')}",
                min_value=0,
                value=int(DESEQ2_DEFAULTS.get("min_base_mean", 10)),
                step=5,
                help=t("viz.basemean_filter_help"),
                key="viz_min_base_mean",
            )
        with col_viz3:
            legend_loc = st.selectbox(
                f"📍 {t('viz.legend_position_label')}",
                options=LEGEND_LOCATIONS,
                index=LEGEND_LOCATIONS.index("upper right"),
                help=t("viz.legend_position_help"),
                key="viz_legend_loc",
            )

        top_10_genes = ranking_df.head(10)[ranking_id_col].tolist() if label_top_genes else None
        _shrinkage = results_df.attrs.get("shrinkage_applied", False)

        tab_volcano, tab_volcano_filt, tab_pca, tab_ma, tab_heatmap = st.tabs([
            f"🌋 {t('viz.tab_volcano_all')}",
            f"🌋 {t('viz.tab_volcano_filt')}",
            f"📊 {t('viz.tab_pca')}",
            f"🎯 {t('viz.tab_ma')}",
            f"🔥 {t('viz.tab_heatmap')}",
        ])

        volcano_title = t("viz.volcano_title")

        with tab_volcano:
            # Always re-render volcano with current baseMean filter
            volcano_all_df = _cached_prepare_volcano(
                results_df,
                alpha=DESEQ2_DEFAULTS["alpha"],
                log2fc_threshold=DESEQ2_DEFAULTS["log2fc_threshold"],
                min_base_mean=viz_min_base_mean,
            )
            fig_volcano_all = _cached_create_volcano(
                volcano_all_df,
                alpha=DESEQ2_DEFAULTS["alpha"],
                log2fc_threshold=DESEQ2_DEFAULTS["log2fc_threshold"],
                test_level=st.session_state.get("test_level", "test"),
                reference_level=st.session_state.get("reference_level", "reference"),
                label_genes=tuple(top_10_genes) if label_top_genes and top_10_genes else None,
                title=volcano_title,
                shrinkage_applied=_shrinkage,
                legend_loc=legend_loc,
            )
            st.pyplot(fig_volcano_all)
            _fig_download_buttons(fig_volcano_all, "volcano_all", "dl_volc_all")

        with tab_volcano_filt:
            if len(filtered_df) == 0:
                st.warning(t("viz.no_genes_volcano"))
            else:
                filt_padj = padj_threshold if filter_option != _radio_none else DESEQ2_DEFAULTS["alpha"]
                filt_log2fc = log2fc_threshold if filter_option != _radio_none else DESEQ2_DEFAULTS["log2fc_threshold"]
                top_10_filtered = [g for g in (top_10_genes or []) if g in filtered_df.index]

                volcano_filt_df = _cached_prepare_volcano(
                    filtered_df, alpha=filt_padj, log2fc_threshold=filt_log2fc,
                    min_base_mean=viz_min_base_mean,
                )
                fig_volcano_filt = _cached_create_volcano(
                    volcano_filt_df,
                    alpha=filt_padj, log2fc_threshold=filt_log2fc,
                    test_level=st.session_state.get("test_level", "test"),
                    reference_level=st.session_state.get("reference_level", "reference"),
                    label_genes=tuple(top_10_filtered) if label_top_genes and top_10_filtered else None,
                    title=volcano_title,
                    shrinkage_applied=_shrinkage,
                    legend_loc=legend_loc,
                )
                st.pyplot(fig_volcano_filt)
                st.caption(
                    tf("viz.filtered_caption",
                       n=len(filtered_df), padj=filt_padj, log2fc=filt_log2fc)
                )
                _fig_download_buttons(fig_volcano_filt, "volcano_filtered", "dl_volc_filt")

        with tab_pca:
            st.pyplot(figures["pca"])
            _fig_download_buttons(figures["pca"], "pca", "dl_pca")

        with tab_ma:
            st.pyplot(figures["ma"])
            _fig_download_buttons(figures["ma"], "ma_plot", "dl_ma")

        with tab_heatmap:
            st.pyplot(figures["heatmap"])
            _fig_download_buttons(figures["heatmap"], "heatmap", "dl_heatmap")

        # ─────────────────────────────────────────────────────────
        # Highlight genes of interest
        # ─────────────────────────────────────────────────────────
        st.header(f"🔎 {t('highlight.header')}")
        st.markdown(t("highlight.instructions"))

        highlight_input = st.text_area(
            t("highlight.textarea_label"),
            height=100,
            placeholder=t("highlight.textarea_placeholder"),
            help=t("highlight.textarea_help"),
        )

        if highlight_input.strip():
            raw_genes = highlight_input.replace(",", "\n").split("\n")
            highlight_genes = [g.strip() for g in raw_genes if g.strip()]

            if highlight_genes:
                volcano_hl_df = _cached_prepare_volcano(
                    results_df,
                    alpha=DESEQ2_DEFAULTS["alpha"],
                    log2fc_threshold=DESEQ2_DEFAULTS["log2fc_threshold"],
                    min_base_mean=viz_min_base_mean,
                )

                found = [g for g in highlight_genes if g in volcano_hl_df.index]
                not_found = [g for g in highlight_genes if g not in volcano_hl_df.index]

                if not_found:
                    _nf_genes = ", ".join(not_found[:10])
                    _nf_msg = tf("highlight.not_found",
                                 count=len(not_found), genes=_nf_genes)
                    if len(not_found) > 10:
                        _nf_msg += tf("highlight.not_found_more",
                                      extra=len(not_found) - 10)
                    st.warning(f"⚠️ {_nf_msg}")

                if found:
                    st.info(
                        f"✅ {tf('highlight.found', count=len(found), genes=', '.join(found))}"
                    )

                    with st.expander(f"📋 {t('highlight.details_title')}", expanded=True):
                        hl_details = results_df.loc[found].sort_values("padj")
                        st.dataframe(hl_details)

                    fig_highlight = _cached_create_volcano_highlight(
                        volcano_hl_df,
                        highlight_genes=tuple(found),
                        alpha=DESEQ2_DEFAULTS["alpha"],
                        log2fc_threshold=DESEQ2_DEFAULTS["log2fc_threshold"],
                        test_level=st.session_state.get("test_level", "test"),
                        reference_level=st.session_state.get("reference_level", "reference"),
                        legend_loc=legend_loc,
                    )
                    st.pyplot(fig_highlight)
                    _fig_download_buttons(fig_highlight, "volcano_highlighted", "dl_highlight")
                else:
                    st.error(f"❌ {t('highlight.none_found')}")
