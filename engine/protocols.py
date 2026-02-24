"""
engine/protocols.py -- Abstract protocols and shared types for the engine.

Defines callback protocols and lightweight data containers used across
engine modules.  No Streamlit dependency.

Types
-----
ProgressCallback
    Protocol — ``(current, total, message_key) -> None``.

EstimateCallback
    Protocol — ``(estimate_dict) -> None``.  The dict contains keys
    like ``estimated_seconds``, ``n_genes``, ``n_samples``,
    ``category`` (see :func:`engine.analysis.estimate_pipeline_time`).

FileData
    NamedTuple — ``(content: bytes, name: str)``.  Wraps uploaded
    files for engine I/O without depending on Streamlit's UploadedFile.
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable, NamedTuple


@runtime_checkable
class ProgressCallback(Protocol):
    """Progress reporting callback.

    Parameters
    ----------
    current : int
        Current step index (0-based).
    total : int
        Total number of steps.
    message_key : str
        An i18n key or plain message describing the current step.
        The frontend is responsible for translating this.
    """

    def __call__(self, current: int, total: int, message_key: str) -> None: ...


@runtime_checkable
class EstimateCallback(Protocol):
    """Called with a time-estimate dict after gene filtering."""

    def __call__(self, estimate: dict) -> None: ...


class FileData(NamedTuple):
    """Minimal file representation for engine I/O.

    Attributes
    ----------
    content : bytes
        Raw file bytes.
    name : str
        Original filename (used for extension-based separator detection).
    """

    content: bytes
    name: str
