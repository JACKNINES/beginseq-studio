"""
app.py — Landing page for BeginSeq Studio.

This is the main entry point of the multipage Streamlit application.
It displays a front page with language selector and navigation cards
for the three analysis tools.

To run:
    streamlit run app.py
"""

import streamlit as st
from i18n import t
from config import AUTO_SHUTDOWN_CONFIG
from auto_shutdown import start_shutdown_watcher

# ─────────────────────────────────────────────────────────────────────
# Auto-shutdown: stop the server when all browser tabs are closed
# ─────────────────────────────────────────────────────────────────────
if AUTO_SHUTDOWN_CONFIG["enabled"]:
    start_shutdown_watcher(
        grace_period=AUTO_SHUTDOWN_CONFIG["grace_period_seconds"],
        poll_interval=AUTO_SHUTDOWN_CONFIG["poll_interval_seconds"],
    )

# ─────────────────────────────────────────────────────────────────────
# Page configuration
# ─────────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="BeginSeq Studio",
    layout="centered",
    page_icon="🧬",
)

# ─────────────────────────────────────────────────────────────────────
# Sidebar: language selector (persists across all pages)
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
        key="lang_selector_home",
        label_visibility="collapsed",
    )
    if _lang_options[_selected_label] != st.session_state.language:
        st.session_state.language = _lang_options[_selected_label]
        st.rerun()

# ─────────────────────────────────────────────────────────────────────
# Main title (with helix animation)
# ─────────────────────────────────────────────────────────────────────
import base64, pathlib

_helix_path = pathlib.Path(__file__).resolve().parent / "helix.MP4"
if _helix_path.exists():
    _helix_b64 = base64.b64encode(_helix_path.read_bytes()).decode()
    st.markdown(
        f"""
        <div style="display:flex; align-items:center; gap:12px; margin-bottom:0.25em;">
            <video autoplay loop muted playsinline
                   style="height:60px; width:auto; border-radius:6px;">
                <source src="data:video/mp4;base64,{_helix_b64}" type="video/mp4">
            </video>
            <h1 style="margin:0; padding:0; font-size:2em;">{t('landing.title')}</h1>
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.title(f"🧬 {t('landing.title')}")
st.markdown(t("landing.subtitle"))
st.markdown("---")

# ─────────────────────────────────────────────────────────────────────
# Tool cards
# ─────────────────────────────────────────────────────────────────────
col1, col2, col3 = st.columns(3)

with col1:
    st.markdown(f"### 🧬 {t('landing.bulk_title')}")
    st.markdown(t("landing.bulk_desc"))
    st.page_link(
        "pages/1_🧬_Bulk_RNA-seq.py",
        label=f"→ {t('landing.open_tool')}",
        icon="🚀",
    )

with col2:
    st.markdown(f"### 🔬 {t('landing.scrna_title')}")
    st.markdown(t("landing.scrna_desc"))
    st.markdown(f"*🚧 {t('landing.coming_soon')}*")
    st.page_link(
        "pages/2_🔬_scRNA-seq.py",
        label=f"→ {t('landing.open_tool')}",
        icon="🔬",
        disabled=True,
    )

with col3:
    st.markdown(f"### 📦 {t('landing.dataset_title')}")
    st.markdown(t("landing.dataset_desc"))
    st.page_link(
        "pages/3_📦_Dataset_Creator.py",
        label=f"→ {t('landing.open_tool')}",
        icon="📦",
    )

# ─────────────────────────────────────────────────────────────────────
# Footer: disclaimer
# ─────────────────────────────────────────────────────────────────────
st.markdown("---")
st.caption(t("landing.disclaimer"))
