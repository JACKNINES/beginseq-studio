"""
runtime_utils.py -- Shared runtime utilities for BeginSeq Studio.

Provides helpers used across multiple pages:
    - Localhost detection for gating cloud vs local features.
    - Upload size limit enforcement per page.
"""

import os

_LOCAL_HOSTS = {"localhost", "127.0.0.1", "::1", ""}

# Streamlit Cloud sets these environment variables on their managed instances.
_CLOUD_ENV_MARKERS = (
    "STREAMLIT_SHARING_MODE",   # legacy Streamlit sharing
    "IS_STREAMLIT_CLOUD",       # newer Streamlit Cloud
    "HOME",                     # /home/adminuser or /home/appuser on Streamlit Cloud
)


def _is_streamlit_cloud() -> bool:
    """Detect Streamlit Cloud by environment markers."""
    # Streamlit Cloud runs under /home/adminuser OR /home/appuser
    # (the username changed in recent Cloud platform updates).
    home = os.environ.get("HOME", "")
    if home in ("/home/adminuser", "/home/appuser"):
        return True
    # Code is mounted under /mount/src/ on Streamlit Cloud
    if os.getcwd().startswith("/mount/src"):
        return True
    # Explicit cloud flags
    if os.environ.get("STREAMLIT_SHARING_MODE"):
        return True
    if os.environ.get("IS_STREAMLIT_CLOUD"):
        return True
    return False


def is_running_locally() -> bool:
    """
    Return True if Streamlit is served from localhost.

    Detection strategy (in order):
    1. Check for Streamlit Cloud environment markers.
       If any are present, we are definitely NOT local.
    2. Check ``streamlit.config.get_option("server.address")``.
       ``None`` or empty string means the default (localhost).
       Note: ``0.0.0.0`` is excluded because cloud deployments
       also bind to ``0.0.0.0``.
    3. Conservative fallback: assume local.
    """
    # Method 1: cloud environment detection (most reliable)
    if _is_streamlit_cloud():
        return False

    # Method 2: configured server address
    try:
        from streamlit import config
        server_addr = config.get_option("server.address")
        # None or empty → default → localhost
        if server_addr is None or server_addr in _LOCAL_HOSTS:
            return True
        # 0.0.0.0 means "listen on all interfaces" — common in both local
        # Docker setups and cloud.  Since we already ruled out cloud above,
        # treat it as local.
        if server_addr == "0.0.0.0":
            return True
        return False
    except Exception:
        pass

    # Fallback: allow
    return True


def apply_local_upload_limit() -> None:
    """
    When running locally, raise the Streamlit upload limit to 5 GB
    (5120 MB) so that large scRNA-seq datasets can be uploaded.

    This MUST be called early — before any ``st.file_uploader`` is
    rendered — because Streamlit reads the option once at startup.

    On cloud deployments the default Streamlit limit (200 MB) is kept.
    """
    if not is_running_locally():
        return

    try:
        from streamlit import config
        config.set_option("server.maxUploadSize", 5120)
    except Exception:
        pass
