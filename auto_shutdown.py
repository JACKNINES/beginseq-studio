"""
auto_shutdown.py — Automatically stop the Streamlit server when all
browser tabs have been closed.

A single daemon thread polls the Streamlit Runtime's active session count.
When no sessions remain for longer than a configurable grace period, it
calls Runtime.stop(), which triggers a clean shutdown of the Tornado
server and all associated resources.

Usage
-----
Call once from the main entry point (app.py):

    from auto_shutdown import start_shutdown_watcher
    start_shutdown_watcher()

The function is idempotent — subsequent calls (from Streamlit reruns)
are silently ignored.
"""

import logging
import os
import signal
import threading
import time

logger = logging.getLogger(__name__)

# ── Defaults (can be overridden via function arguments) ──────────────
_DEFAULT_GRACE_PERIOD = 30   # seconds after last session disconnects
_DEFAULT_POLL_INTERVAL = 5   # seconds between session-count checks

# ── Module-level guard ───────────────────────────────────────────────
_watcher_started = False
_lock = threading.Lock()


def start_shutdown_watcher(
    grace_period: float = _DEFAULT_GRACE_PERIOD,
    poll_interval: float = _DEFAULT_POLL_INTERVAL,
) -> None:
    """Start the shutdown watcher daemon thread.

    This function is idempotent — calling it multiple times (from
    multiple Streamlit reruns) will only start one watcher thread.

    Parameters
    ----------
    grace_period : float
        Seconds to wait after the last session disconnects before
        shutting down.  30 s is enough to survive page refreshes and
        multipage navigation.
    poll_interval : float
        Seconds between active-session checks.
    """
    global _watcher_started

    with _lock:
        if _watcher_started:
            return
        _watcher_started = True

    thread = threading.Thread(
        target=_watcher_loop,
        args=(grace_period, poll_interval),
        daemon=True,
        name="auto-shutdown-watcher",
    )
    thread.start()
    logger.info(
        "Auto-shutdown watcher started (grace=%ds, poll=%ds)",
        grace_period,
        poll_interval,
    )


# ── Private helpers ──────────────────────────────────────────────────

def _watcher_loop(grace_period: float, poll_interval: float) -> None:
    """Background loop that monitors active sessions."""
    from streamlit.runtime import Runtime

    # Wait for the Runtime to be fully initialised
    while not Runtime.exists():
        time.sleep(1)

    runtime = Runtime.instance()

    # Wait for at least one session to connect so we don't shut down
    # before the browser has even opened
    while runtime._session_mgr.num_active_sessions() == 0:
        time.sleep(poll_interval)

    # ── Monitor for disconnections ───────────────────────────────
    no_session_since: float | None = None

    while True:
        time.sleep(poll_interval)

        active = runtime._session_mgr.num_active_sessions()

        if active > 0:
            no_session_since = None
            continue

        # No active sessions — start / continue the countdown
        if no_session_since is None:
            no_session_since = time.monotonic()
            logger.info("All sessions disconnected. Grace period started.")

        elapsed = time.monotonic() - no_session_since

        if elapsed >= grace_period:
            logger.info(
                "No active sessions for %.0f s. Shutting down.", elapsed
            )
            _shutdown(runtime)
            return


def _shutdown(runtime) -> None:
    """Perform a clean shutdown of the Streamlit server."""
    try:
        runtime.stop()
    except Exception:
        logger.exception("Runtime.stop() failed — sending SIGTERM")
        os.kill(os.getpid(), signal.SIGTERM)
