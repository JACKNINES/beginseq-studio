"""
runtime_utils.py -- Shared runtime utilities for BeginSeq Studio.

Provides helpers used across multiple pages:
    - Localhost detection for gating cloud vs local features.
    - Upload size limit enforcement per page.
"""

_LOCAL_HOSTS = {"localhost", "127.0.0.1", "0.0.0.0", "::1", ""}


def is_running_locally() -> bool:
    """
    Return True if Streamlit is served from localhost.

    Detection strategy (in order):
    1. Check ``streamlit.config.get_option("server.address")``.
       ``None`` or empty string means the default (localhost).
    2. Parse the URL built by the Streamlit server utility.
    3. Conservative fallback: assume local.
    """
    # Method 1: configured server address
    try:
        from streamlit import config
        server_addr = config.get_option("server.address")
        if server_addr is None or server_addr in _LOCAL_HOSTS:
            return True
        return False
    except Exception:
        pass

    # Method 2: parse the URL Streamlit builds
    try:
        from urllib.parse import urlparse
        from streamlit.web.server.server_util import get_url
        host = urlparse(get_url("")).hostname or ""
        if host in _LOCAL_HOSTS:
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
