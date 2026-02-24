"""
pages/3_📦_Dataset_Creator.py — GDC/TCGA Dataset Creator.

Downloads RNA-seq STAR-Counts datasets from the Genomic Data Commons,
assembles DESeq2-ready count matrices with auto-generated metadata.
"""

import base64
import pandas as pd
import streamlit as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from i18n import t, tf
from engine.config import set_theme, THEME_PRESETS
from engine.gdc_client import (
    fetch_tcga_projects,
    fetch_rnaseq_files,
    download_and_extract_batch,
    build_count_matrix,
)


# ─────────────────────────────────────────────────────────────────────
# Cached wrappers — avoid redundant GDC API calls on Streamlit reruns
# ─────────────────────────────────────────────────────────────────────

@st.cache_data(ttl=86400, show_spinner=False)
def _cached_fetch_projects():
    """Fetch TCGA projects (cached for 24 h)."""
    return fetch_tcga_projects()


@st.cache_data(ttl=3600, show_spinner=False)
def _cached_fetch_files(project_id: str):
    """Fetch RNA-seq files for a project (cached for 1 h)."""
    return fetch_rnaseq_files(project_id)

# ─────────────────────────────────────────────────────────────────────
# Page config
# ─────────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Dataset Creator",
    layout="wide",
    page_icon="📦",
)

# ─────────────────────────────────────────────────────────────────────
# Sidebar: language selector
# ─────────────────────────────────────────────────────────────────────
if "language" not in st.session_state:
    st.session_state.language = "en"
if "plot_theme" not in st.session_state:
    st.session_state.plot_theme = "dark"

with st.sidebar:
    _lang_options = {"English": "en", "Español": "es"}
    _current_label = "English" if st.session_state.language == "en" else "Español"
    _selected_label = st.selectbox(
        "🌐",
        options=list(_lang_options.keys()),
        index=list(_lang_options.keys()).index(_current_label),
        key="lang_selector_dataset",
        label_visibility="collapsed",
    )
    if _lang_options[_selected_label] != st.session_state.language:
        st.session_state.language = _lang_options[_selected_label]
        st.rerun()

    # ── Theme selector ────────────────────────────────────────────────
    _theme_labels = {"🌑 Dark": "dark", "☀️ Light": "light", "⚡ Cyberpunk": "cyberpunk"}
    _current_theme_label = {v: k for k, v in _theme_labels.items()}[st.session_state.plot_theme]
    _selected_theme = st.selectbox(
        "🎨",
        options=list(_theme_labels.keys()),
        index=list(_theme_labels.keys()).index(_current_theme_label),
        key="theme_selector_dataset",
        label_visibility="collapsed",
    )
    if _theme_labels[_selected_theme] != st.session_state.plot_theme:
        st.session_state.plot_theme = _theme_labels[_selected_theme]
        st.rerun()

# Apply the selected theme to the engine + UI
set_theme(st.session_state.plot_theme)
from theme_css import inject_theme_css
inject_theme_css()

# ─────────────────────────────────────────────────────────────────────
# Title
# ─────────────────────────────────────────────────────────────────────
_box_path = Path(__file__).resolve().parent.parent / "assets" / "Box.MP4"
if _box_path.exists():
    _box_b64 = base64.b64encode(_box_path.read_bytes()).decode()
    st.markdown(
        f"""
        <div style="display:flex; align-items:center; gap:12px; margin-bottom:0.25em;">
            <video autoplay loop muted playsinline
                   style="height:60px; width:auto; border-radius:6px;">
                <source src="data:video/mp4;base64,{_box_b64}" type="video/mp4">
            </video>
            <h1 style="margin:0; padding:0; font-size:2em;">{t('dc.title')}</h1>
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.title(f"📦 {t('dc.title')}")
st.markdown(t("dc.subtitle"))
st.markdown("---")

# ─────────────────────────────────────────────────────────────────────
# Initialize session state
# ─────────────────────────────────────────────────────────────────────
for key in ("dc_projects", "dc_files", "dc_counts", "dc_metadata"):
    if key not in st.session_state:
        st.session_state[key] = None

# ═════════════════════════════════════════════════════════════════════
# STEP 1: Load TCGA projects
# ═════════════════════════════════════════════════════════════════════
st.subheader(f"1️⃣ {t('dc.step1_title')}")
st.markdown(t("dc.step1_desc"))

if st.button(t("dc.load_projects"), key="btn_load_projects"):
    with st.spinner(t("dc.loading_projects")):
        try:
            projects = _cached_fetch_projects()
            st.session_state.dc_projects = projects
            st.session_state.dc_files = None
            st.session_state.dc_counts = None
            st.session_state.dc_metadata = None
        except Exception as exc:
            st.error(tf("dc.error_api", error=str(exc)))

if st.session_state.dc_projects is not None:
    projects_df = st.session_state.dc_projects

    if projects_df.empty:
        st.warning(t("dc.no_projects"))
    else:
        st.success(tf("dc.projects_loaded", count=len(projects_df)))

        with st.expander(t("dc.view_projects"), expanded=False):
            st.dataframe(projects_df, use_container_width=True, hide_index=True)

        # ═════════════════════════════════════════════════════════════
        # STEP 2: Select a project and fetch files
        # ═════════════════════════════════════════════════════════════
        st.markdown("---")
        st.subheader(f"2️⃣ {t('dc.step2_title')}")

        proj_options = [
            f"{row['project_id']} — {row['name']}"
            for _, row in projects_df.iterrows()
        ]
        selected_proj_label = st.selectbox(
            t("dc.select_project"),
            options=proj_options,
            key="dc_project_selector",
        )
        selected_project_id = selected_proj_label.split(" — ")[0]

        if st.button(t("dc.fetch_files"), key="btn_fetch_files"):
            with st.spinner(tf("dc.fetching_files", project=selected_project_id)):
                try:
                    files_df = _cached_fetch_files(selected_project_id)
                    st.session_state.dc_files = files_df
                    st.session_state.dc_counts = None
                    st.session_state.dc_metadata = None
                except Exception as exc:
                    st.error(tf("dc.error_api", error=str(exc)))

        if st.session_state.dc_files is not None:
            files_df = st.session_state.dc_files

            if files_df.empty:
                st.warning(tf("dc.no_files", project=selected_project_id))
            else:
                # ─────────────────────────────────────────────────────
                # Summary metrics
                # ─────────────────────────────────────────────────────
                n_total = len(files_df)
                n_tumor = len(files_df[files_df["condition"] == "Tumor"])
                n_normal = len(files_df[files_df["condition"] == "Normal"])
                total_size_mb = files_df["file_size"].sum() / (1024 * 1024)

                st.success(tf("dc.files_found", count=n_total, project=selected_project_id))

                mc1, mc2, mc3, mc4 = st.columns(4)
                mc1.metric(t("dc.metric_total"), n_total)
                mc2.metric(t("dc.metric_tumor"), n_tumor)
                mc3.metric(t("dc.metric_normal"), n_normal)
                mc4.metric(t("dc.metric_size"), f"{total_size_mb:.1f} MB")

                with st.expander(t("dc.view_files"), expanded=False):
                    display_cols = ["file_id", "case_id", "sample_type", "condition", "file_size"]
                    st.dataframe(
                        files_df[display_cols],
                        use_container_width=True,
                        hide_index=True,
                    )

                # ═════════════════════════════════════════════════════
                # STEP 3: Filter and configure
                # ═════════════════════════════════════════════════════
                st.markdown("---")
                st.subheader(f"3️⃣ {t('dc.step3_title')}")

                # ── Filter by condition ───────────────────────────
                st.markdown(f"**{t('dc.filter_samples')}**")
                available_conditions = sorted(files_df["condition"].unique().tolist())
                selected_conditions = st.multiselect(
                    t("dc.select_conditions"),
                    options=available_conditions,
                    default=available_conditions,
                    key="dc_condition_filter",
                )

                if not selected_conditions:
                    st.warning(t("dc.no_conditions_selected"))
                    filtered_df = files_df.head(0)
                else:
                    filtered_df = files_df[files_df["condition"].isin(selected_conditions)]

                if not filtered_df.empty:
                    n_filt_tumor = len(filtered_df[filtered_df["condition"] == "Tumor"])
                    n_filt_normal = len(filtered_df[filtered_df["condition"] == "Normal"])

                    st.info(
                        tf("dc.filtered_summary",
                           count=len(filtered_df),
                           tumor=n_filt_tumor,
                           normal=n_filt_normal)
                    )

                    # ── Sample limit ──────────────────────────────
                    st.markdown(f"**{t('dc.sample_limit_title')}**")
                    st.caption(t("dc.sample_limit_desc"))

                    use_limit = st.checkbox(
                        t("dc.use_sample_limit"),
                        value=False,
                        key="dc_use_limit",
                    )

                    if use_limit:
                        max_available = len(filtered_df)
                        max_per_condition = max(
                            filtered_df["condition"].value_counts().max(), 1
                        )

                        limit_mode = st.radio(
                            t("dc.limit_mode_label"),
                            options=[t("dc.limit_mode_total"), t("dc.limit_mode_per_cond")],
                            key="dc_limit_mode",
                            horizontal=True,
                        )

                        if limit_mode == t("dc.limit_mode_total"):
                            sample_limit = st.slider(
                                t("dc.limit_total_slider"),
                                min_value=2,
                                max_value=max_available,
                                value=min(100, max_available),
                                step=1,
                                key="dc_limit_total",
                            )
                            # Proportional sampling across conditions
                            sampled_parts = []
                            for cond in filtered_df["condition"].unique():
                                cond_df = filtered_df[filtered_df["condition"] == cond]
                                proportion = len(cond_df) / len(filtered_df)
                                n_for_cond = max(1, round(sample_limit * proportion))
                                sampled_parts.append(
                                    cond_df.sample(n=min(n_for_cond, len(cond_df)), random_state=42)
                                )
                            filtered_df = pd.concat(sampled_parts).head(sample_limit)

                        else:  # per condition
                            per_cond_limit = st.slider(
                                t("dc.limit_per_cond_slider"),
                                min_value=1,
                                max_value=int(max_per_condition),
                                value=min(50, int(max_per_condition)),
                                step=1,
                                key="dc_limit_per_cond",
                            )
                            sampled_parts = []
                            for cond in filtered_df["condition"].unique():
                                cond_df = filtered_df[filtered_df["condition"] == cond]
                                sampled_parts.append(
                                    cond_df.sample(n=min(per_cond_limit, len(cond_df)), random_state=42)
                                )
                            filtered_df = pd.concat(sampled_parts)

                        # Show updated counts after limiting
                        n_lim_tumor = len(filtered_df[filtered_df["condition"] == "Tumor"])
                        n_lim_normal = len(filtered_df[filtered_df["condition"] == "Normal"])
                        st.success(
                            tf("dc.limited_summary",
                               count=len(filtered_df),
                               tumor=n_lim_tumor,
                               normal=n_lim_normal)
                        )

                    # ── Configuration options ─────────────────────
                    cfg_col1, cfg_col2 = st.columns(2)
                    with cfg_col1:
                        gene_id_type = st.selectbox(
                            t("dc.gene_id_label"),
                            options=["gene_name", "gene_id"],
                            index=0,
                            help=t("dc.gene_id_help"),
                            key="dc_gene_id_type",
                        )
                    with cfg_col2:
                        count_column = st.selectbox(
                            t("dc.count_col_label"),
                            options=["unstranded", "stranded_first", "stranded_second"],
                            index=0,
                            help=t("dc.count_col_help"),
                            key="dc_count_column",
                        )

                    # ═════════════════════════════════════════════════
                    # STEP 4: Download & Build
                    # ═════════════════════════════════════════════════
                    st.markdown("---")
                    st.subheader(f"4️⃣ {t('dc.step4_title')}")

                    est_size = filtered_df["file_size"].sum() / (1024 * 1024)
                    st.markdown(
                        tf("dc.download_info",
                           n_files=len(filtered_df),
                           size=f"{est_size:.1f}")
                    )

                    # ── Download directory selector ───────────────
                    default_dir = str(Path.home() / "Desktop" / "gdc_downloads")
                    download_dir_input = st.text_input(
                        t("dc.download_dir_label"),
                        value=default_dir,
                        help=t("dc.download_dir_help"),
                        key="dc_download_dir",
                    )

                    # Resolve the directory path
                    user_dest_dir: Path | None = None
                    if download_dir_input.strip():
                        user_dest_dir = Path(download_dir_input.strip()).expanduser()
                        st.caption(
                            tf("dc.download_dir_persistent_info",
                               path=str(user_dest_dir))
                        )

                    # Warn about large downloads
                    if len(filtered_df) > 200:
                        st.warning(t("dc.large_download_warning"))

                    if st.button(t("dc.start_download"), key="btn_start_download", type="primary"):
                        file_ids = filtered_df["file_id"].tolist()

                        # ── Phase 1: Download & Extract ───────────
                        progress_bar = st.progress(0, text=t("dc.downloading"))
                        status_text = st.empty()

                        def download_progress(current: int, total: int, msg: str):
                            """Update progress bar during file download (0–50%)."""
                            pct = current / total * 0.5 if total > 0 else 0
                            progress_bar.progress(
                                min(pct, 0.5),
                                text=tf("dc.download_progress",
                                        current=current, total=total),
                            )
                            status_text.text(msg)

                        try:
                            extracted_files, tmp_dir = download_and_extract_batch(
                                file_ids,
                                dest_dir=user_dest_dir,
                                progress_callback=download_progress,
                            )
                        except Exception as exc:
                            st.error(tf("dc.error_download", error=str(exc)))
                            st.stop()

                        if not extracted_files:
                            st.error(t("dc.no_files_downloaded"))
                            st.stop()

                        status_text.text(
                            tf("dc.downloaded_count", count=len(extracted_files))
                        )

                        # ── Phase 2: Build matrix ─────────────────
                        progress_bar.progress(0.5, text=t("dc.building_matrix"))

                        def parse_progress(current: int, total: int, msg: str):
                            """Update progress bar during matrix assembly (50–100%)."""
                            pct = 0.5 + (current / total * 0.5) if total > 0 else 0.5
                            progress_bar.progress(
                                min(pct, 1.0),
                                text=tf("dc.parse_progress",
                                        current=current, total=total),
                            )

                        try:
                            counts_df, metadata_df = build_count_matrix(
                                extracted_files,
                                filtered_df,
                                gene_id_type=gene_id_type,
                                count_column=count_column,
                                progress_callback=parse_progress,
                            )
                        except Exception as exc:
                            st.error(tf("dc.error_building", error=str(exc)))
                            st.stop()

                        progress_bar.progress(1.0, text=t("dc.done"))
                        status_text.empty()

                        st.session_state.dc_counts = counts_df
                        st.session_state.dc_metadata = metadata_df

                        st.success(
                            tf("dc.build_success",
                               genes=counts_df.shape[0],
                               samples=counts_df.shape[1])
                        )
                        st.rerun()

                    # ═════════════════════════════════════════════════
                    # STEP 5: Preview & Download results
                    # ═════════════════════════════════════════════════
                    if st.session_state.dc_counts is not None:
                        counts_df = st.session_state.dc_counts
                        metadata_df = st.session_state.dc_metadata

                        st.markdown("---")
                        st.subheader(f"5️⃣ {t('dc.step5_title')}")

                        st.success(
                            tf("dc.build_success",
                               genes=counts_df.shape[0],
                               samples=counts_df.shape[1])
                        )

                        # Preview counts matrix
                        st.markdown(f"**{t('dc.preview_counts')}**")
                        st.caption(
                            tf("dc.preview_counts_caption",
                               genes=counts_df.shape[0],
                               samples=counts_df.shape[1])
                        )
                        st.dataframe(counts_df.head(20), use_container_width=True)

                        # Preview metadata
                        st.markdown(f"**{t('dc.preview_metadata')}**")
                        st.dataframe(metadata_df, use_container_width=True, hide_index=True)

                        # ── Condition summary ─────────────────────
                        cond_counts = metadata_df["condition"].value_counts()
                        st.markdown(f"**{t('dc.condition_summary')}**")
                        for cond, cnt in cond_counts.items():
                            st.markdown(f"- **{cond}**: {cnt}")

                        # ── Download buttons ──────────────────────
                        st.markdown("---")
                        st.markdown(f"**{t('dc.download_results')}**")

                        dl_col1, dl_col2 = st.columns(2)

                        with dl_col1:
                            counts_csv = counts_df.to_csv()
                            st.download_button(
                                label=t("dc.download_counts"),
                                data=counts_csv,
                                file_name="counts_matrix.csv",
                                mime="text/csv",
                                key="btn_dl_counts",
                            )

                        with dl_col2:
                            meta_csv = metadata_df.to_csv(index=False)
                            st.download_button(
                                label=t("dc.download_metadata"),
                                data=meta_csv,
                                file_name="metadata.csv",
                                mime="text/csv",
                                key="btn_dl_metadata",
                            )

                        # ── Tip for Bulk RNA-seq ──────────────────
                        st.info(t("dc.tip_bulk"))
