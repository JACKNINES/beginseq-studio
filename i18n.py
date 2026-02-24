"""
i18n.py -- Internationalization for the RNA-seq Analyzer.

Thin Streamlit adapter over engine.i18n_labels.
Reads st.session_state.language to determine which language to use.

Usage:
    from i18n import t, tf

    st.title(t("app.title"))
    st.caption(tf("preview.counts_caption", genes=100, samples=12))
"""

import streamlit as st

from engine.i18n_labels import (
    TRANSLATIONS,
    get_label,
    get_label_formatted,
)


def t(key: str) -> str:
    """
    Return the translated string for the given key.

    Reads st.session_state.language to determine which language to use.
    Falls back to English if key is missing in the selected language.
    Falls back to the key itself if not found in any language.
    """
    lang = st.session_state.get("language", "en")
    return get_label(key, lang)


def tf(key: str, **kwargs) -> str:
    """
    Return the translated string with format placeholders replaced.

    Usage:
        tf("preview.counts_caption", genes=100, samples=12)
        # English: "100 genes x 12 samples"
        # Spanish: "100 genes x 12 muestras"
    """
    lang = st.session_state.get("language", "en")
    return get_label_formatted(key, lang, **kwargs)
