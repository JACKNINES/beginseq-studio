"""
pages/2_🔬_scRNA-seq.py — scRNA-seq Analysis placeholder page.

This tool is under development.
"""

import base64
import streamlit as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from i18n import t

# ─────────────────────────────────────────────────────────────────────
# Page configuration
# ─────────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="scRNA-seq Analysis",
    layout="centered",
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
# Content
# ─────────────────────────────────────────────────────────────────────
_cell_path = Path(__file__).resolve().parent.parent / "Cell.MP4"
if _cell_path.exists():
    _cell_b64 = base64.b64encode(_cell_path.read_bytes()).decode()
    st.markdown(
        f"""
        <div style="display:flex; align-items:center; gap:12px; margin-bottom:0.25em;">
            <video autoplay loop muted playsinline
                   style="height:60px; width:auto; border-radius:6px;">
                <source src="data:video/mp4;base64,{_cell_b64}" type="video/mp4">
            </video>
            <h1 style="margin:0; padding:0; font-size:2em;">{t('placeholder.title_scrna')}</h1>
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.title(f"🔬 {t('placeholder.title_scrna')}")

st.markdown("---")

st.info(f"🚧 {t('placeholder.message')}")

st.markdown(t("placeholder.scrna_preview"))
