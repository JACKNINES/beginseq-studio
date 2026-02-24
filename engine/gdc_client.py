"""
gdc_client.py — GDC REST API client for downloading TCGA RNA-seq datasets.

This module provides functions to:
1. List available TCGA projects from the GDC.
2. Search for STAR-Counts RNA-seq files for a given project.
3. Download gene expression files in batch (single tar.gz request).
4. Extract and decompress STAR-Counts files from GDC archives.
5. Parse STAR-Counts TSV files into count vectors.
6. Assemble a DESeq2-ready count matrix + metadata.

All communication uses the official GDC REST API (https://api.gdc.cancer.gov).
No R, no TCGAbiolinks, no external wrappers.

Functions
---------
fetch_tcga_projects()
    → List available TCGA projects from the GDC.

fetch_rnaseq_files(project_id)
    → Search for STAR-Counts RNA-seq files for a given project.

download_and_extract_batch(file_ids, …)
    → Download and decompress gene expression files in a single tar.gz request.

parse_star_counts(text)
    → Parse a single STAR-Counts TSV into a ``{gene_id: count}`` dict.

build_count_matrix(files, …)
    → Assemble a DESeq2-ready count matrix + metadata from downloaded files.

Usage:
    from gdc_client import (
        fetch_tcga_projects,
        fetch_rnaseq_files,
        download_and_extract_batch,
        build_count_matrix,
    )
"""

from __future__ import annotations

import gzip
import io
import json
import shutil
import tarfile
import tempfile
import time
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
import requests

from engine.config import MEMORY_CONFIG

# ══════════════════════════════════════════════════════════════════════
# Constants
# ══════════════════════════════════════════════════════════════════════

GDC_BASE = "https://api.gdc.cancer.gov"
GDC_PROJECTS_ENDPOINT = f"{GDC_BASE}/projects"
GDC_FILES_ENDPOINT = f"{GDC_BASE}/files"
GDC_DATA_ENDPOINT = f"{GDC_BASE}/data"

# Rows in STAR-Counts TSV that are QC summaries, not gene counts
SPECIAL_ROWS = {"N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"}

# HTTP settings
REQUEST_TIMEOUT = 60  # seconds per request
DOWNLOAD_TIMEOUT = 600  # seconds for batch download (large archives)
MAX_RETRIES = 3
BACKOFF_BASE = 2  # exponential backoff base (seconds)

# Pagination
PAGE_SIZE = 2000  # max records per API page

# Batch download: GDC caps at ~600 files per POST to /data
BATCH_CHUNK_SIZE = 500


# ══════════════════════════════════════════════════════════════════════
# Internal helpers
# ══════════════════════════════════════════════════════════════════════

def _request_with_retry(
    method: str,
    url: str,
    retries: int = MAX_RETRIES,
    **kwargs,
) -> requests.Response:
    """
    Make an HTTP request with exponential-backoff retry.

    Parameters
    ----------
    method : str
        HTTP method ("GET" or "POST").
    url : str
        Full URL to request.
    retries : int
        Number of retry attempts.
    **kwargs
        Extra arguments forwarded to ``requests.request()``.

    Returns
    -------
    requests.Response

    Raises
    ------
    requests.HTTPError
        After all retries are exhausted.
    """
    kwargs.setdefault("timeout", REQUEST_TIMEOUT)
    last_exc: Exception | None = None

    for attempt in range(retries):
        try:
            resp = requests.request(method, url, **kwargs)
            resp.raise_for_status()
            return resp
        except (requests.RequestException, requests.HTTPError) as exc:
            last_exc = exc
            if attempt < retries - 1:
                wait = BACKOFF_BASE ** attempt
                time.sleep(wait)

    raise last_exc  # type: ignore[misc]


def _build_filter(filters: dict) -> dict:
    """
    Build a GDC API filter payload from a flat dict of {field: value}.

    If value is a list, uses the ``in`` operator; otherwise ``=``.
    """
    content = []
    for field, value in filters.items():
        if isinstance(value, list):
            content.append({
                "op": "in",
                "content": {"field": field, "value": value},
            })
        else:
            content.append({
                "op": "=",
                "content": {"field": field, "value": value},
            })
    return {"op": "and", "content": content}


def _classify_sample_type(sample_type: str) -> str:
    """Classify a TCGA sample_type string into Tumor/Normal/Other."""
    if not sample_type:
        return "Unknown"
    st_lower = sample_type.lower()
    if "normal" in st_lower:
        return "Normal"
    if "tumor" in st_lower or "cancer" in st_lower or "metastatic" in st_lower:
        return "Tumor"
    return "Other"


def _decompress_gz(gz_path: Path) -> Path:
    """
    Decompress a .gz file in-place, returning the path to the
    decompressed file. The .gz file is removed after extraction.
    """
    out_path = gz_path.with_suffix("")  # remove .gz extension
    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    gz_path.unlink()
    return out_path


# ══════════════════════════════════════════════════════════════════════
# 1. List TCGA projects
# ══════════════════════════════════════════════════════════════════════

def fetch_tcga_projects() -> pd.DataFrame:
    """
    Fetch the list of TCGA projects from GDC.

    Returns
    -------
    pd.DataFrame
        Columns: project_id, name, primary_site, disease_type, file_count.
        Sorted by project_id.
    """
    filters = _build_filter({"program.name": "TCGA"})
    params = {
        "filters": json.dumps(filters),
        "fields": (
            "project_id,"
            "name,"
            "primary_site,"
            "disease_type,"
            "summary.file_count"
        ),
        "size": "100",
        "format": "json",
    }

    resp = _request_with_retry("GET", GDC_PROJECTS_ENDPOINT, params=params)
    data = resp.json()

    hits = data.get("data", {}).get("hits", [])
    rows = []
    for h in hits:
        rows.append({
            "project_id": h.get("project_id", ""),
            "name": h.get("name", ""),
            "primary_site": (
                ", ".join(h["primary_site"])
                if isinstance(h.get("primary_site"), list)
                else h.get("primary_site", "")
            ),
            "disease_type": (
                ", ".join(h["disease_type"])
                if isinstance(h.get("disease_type"), list)
                else h.get("disease_type", "")
            ),
            "file_count": (
                h.get("summary", {}).get("file_count", 0)
                if isinstance(h.get("summary"), dict)
                else 0
            ),
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(
            columns=["project_id", "name", "primary_site", "disease_type", "file_count"]
        )
    return df.sort_values("project_id").reset_index(drop=True)


# ══════════════════════════════════════════════════════════════════════
# 2. Search RNA-seq files for a project
# ══════════════════════════════════════════════════════════════════════

def fetch_rnaseq_files(project_id: str) -> pd.DataFrame:
    """
    Find open-access STAR-Counts RNA-seq files for a TCGA project.

    Parameters
    ----------
    project_id : str
        TCGA project ID, e.g. ``"TCGA-BRCA"``.

    Returns
    -------
    pd.DataFrame
        Columns: file_id, file_name, file_size, case_id,
                 sample_type, sample_type_id, condition.
        One row per file.
    """
    filters = _build_filter({
        "cases.project.project_id": project_id,
        "data_category": "Transcriptome Profiling",
        "data_type": "Gene Expression Quantification",
        "experimental_strategy": "RNA-Seq",
        "analysis.workflow_type": "STAR - Counts",
        "access": "open",
    })

    fields = ",".join([
        "file_id",
        "file_name",
        "file_size",
        "cases.submitter_id",
        "cases.samples.sample_type",
        "cases.samples.sample_type_id",
    ])

    all_hits: list[dict] = []
    offset = 0

    while True:
        payload = {
            "filters": json.dumps(filters),
            "fields": fields,
            "format": "json",
            "size": str(PAGE_SIZE),
            "from": str(offset),
        }
        resp = _request_with_retry("POST", GDC_FILES_ENDPOINT, data=payload)
        data = resp.json()

        hits = data.get("data", {}).get("hits", [])
        if not hits:
            break
        all_hits.extend(hits)

        pagination = data.get("data", {}).get("pagination", {})
        total = pagination.get("total", 0)
        offset += PAGE_SIZE
        if offset >= total:
            break

    rows = []
    for h in all_hits:
        file_id = h.get("file_id", "")
        file_name = h.get("file_name", "")
        file_size = h.get("file_size", 0)

        # Extract case/sample info (nested structure)
        cases = h.get("cases", [])
        case_id = ""
        sample_type = ""
        sample_type_id = ""

        if cases:
            case = cases[0]
            case_id = case.get("submitter_id", "")
            samples = case.get("samples", [])
            if samples:
                sample = samples[0]
                sample_type = sample.get("sample_type", "")
                sample_type_id = sample.get("sample_type_id", "")

        rows.append({
            "file_id": file_id,
            "file_name": file_name,
            "file_size": file_size,
            "case_id": case_id,
            "sample_type": sample_type,
            "sample_type_id": sample_type_id,
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(
            columns=["file_id", "file_name", "file_size",
                      "case_id", "sample_type", "sample_type_id", "condition"]
        )

    df["condition"] = df["sample_type"].apply(_classify_sample_type)
    return df.reset_index(drop=True)


# ══════════════════════════════════════════════════════════════════════
# 3. Batch download & extraction
# ══════════════════════════════════════════════════════════════════════

def download_and_extract_batch(
    file_ids: list[str],
    dest_dir: Path | None = None,
    progress_callback: Callable[[int, int, str], None] | None = None,
) -> tuple[dict[str, Path], Path]:
    """
    Download multiple GDC files using the batch endpoint (POST /data)
    which returns a tar.gz archive, then extract and decompress all
    STAR-Counts TSV files.

    The GDC batch endpoint packs files into:
        archive.tar.gz
        └── <uuid>/
            └── <filename>.rna_seq.augmented_star_gene_counts.tsv  (may be .gz)

    Parameters
    ----------
    file_ids : list[str]
        List of GDC file UUIDs to download.
    dest_dir : Path or None
        Directory to save extracted files. If None, a temp directory is created.
    progress_callback : callable or None
        Called as ``progress_callback(current, total, message)``
        to report progress.

    Returns
    -------
    tuple[dict[str, Path], Path]
        (mapping of file_id -> extracted TSV path, directory used).
    """
    if dest_dir is None:
        dest_dir = Path(tempfile.mkdtemp(prefix="gdc_downloads_"))
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    total = len(file_ids)
    extracted: dict[str, Path] = {}

    # Split into chunks to respect GDC limits
    chunks = [
        file_ids[i : i + BATCH_CHUNK_SIZE]
        for i in range(0, total, BATCH_CHUNK_SIZE)
    ]

    files_processed = 0

    for chunk_idx, chunk in enumerate(chunks):
        if progress_callback:
            progress_callback(
                files_processed, total,
                f"Downloading batch {chunk_idx + 1}/{len(chunks)} "
                f"({len(chunk)} files)..."
            )

        # POST to /data with JSON body containing file IDs
        payload = json.dumps({"ids": chunk})
        try:
            resp = requests.post(
                GDC_DATA_ENDPOINT,
                data=payload,
                headers={"Content-Type": "application/json"},
                timeout=DOWNLOAD_TIMEOUT,
                stream=True,
            )
            resp.raise_for_status()
        except Exception as exc:
            if progress_callback:
                progress_callback(
                    files_processed, total,
                    f"Error downloading batch {chunk_idx + 1}: {exc}"
                )
            continue

        # ── Stream response to disk with progress ─────────────────
        content_length = int(resp.headers.get("Content-Length", 0))
        archive_path = dest_dir / f"_batch_{chunk_idx}.tar.gz"

        downloaded_bytes = 0
        _STREAM_CHUNK = 1024 * 256  # 256 KB chunks

        with open(archive_path, "wb") as f:
            for data in resp.iter_content(chunk_size=_STREAM_CHUNK):
                f.write(data)
                downloaded_bytes += len(data)
                if progress_callback and content_length > 0:
                    mb_done = downloaded_bytes / (1024 * 1024)
                    mb_total = content_length / (1024 * 1024)
                    progress_callback(
                        files_processed, total,
                        f"Downloading batch {chunk_idx + 1}/{len(chunks)}: "
                        f"{mb_done:.1f} / {mb_total:.1f} MB"
                    )

        resp.close()

        if progress_callback:
            progress_callback(
                files_processed, total,
                f"Extracting batch {chunk_idx + 1}/{len(chunks)}..."
            )

        # ── Extract archive ───────────────────────────────────────
        content_type = resp.headers.get("Content-Type", "")

        if "octet-stream" in content_type or "gzip" in content_type or len(chunk) > 1:
            # Batch response: tar.gz archive
            try:
                with tarfile.open(archive_path, mode="r:gz") as tar:
                    tar.extractall(path=dest_dir, filter="data")

                # Walk extracted directories: each UUID folder contains the TSV
                for uuid in chunk:
                    uuid_dir = dest_dir / uuid
                    if not uuid_dir.is_dir():
                        continue

                    for fpath in uuid_dir.iterdir():
                        if fpath.suffix == ".gz":
                            fpath = _decompress_gz(fpath)
                        if fpath.suffix in (".tsv", ".txt"):
                            extracted[uuid] = fpath
                            break

                    files_processed += 1
                    if progress_callback:
                        progress_callback(
                            files_processed, total,
                            f"Extracted {files_processed}/{total}"
                        )

            except (tarfile.TarError, gzip.BadGzipFile, OSError) as exc:
                if progress_callback:
                    progress_callback(
                        files_processed, total,
                        f"Error extracting archive: {exc}"
                    )
        else:
            # Single-file response (when only 1 file requested)
            uuid = chunk[0]
            cd = resp.headers.get("Content-Disposition", "")
            if "filename=" in cd:
                fname = cd.split("filename=")[-1].strip('" ')
            else:
                fname = f"{uuid}.tsv"

            fpath = dest_dir / fname
            archive_path.rename(fpath)

            # Decompress if gzipped
            if fpath.suffix == ".gz":
                fpath = _decompress_gz(fpath)

            extracted[uuid] = fpath
            files_processed += 1
            if progress_callback:
                progress_callback(files_processed, total, fname)

        # Clean up archive file after extraction
        if archive_path.exists():
            archive_path.unlink()

    return extracted, dest_dir


# ══════════════════════════════════════════════════════════════════════
# 4. Parse STAR-Counts TSV
# ══════════════════════════════════════════════════════════════════════

def parse_star_counts(
    filepath: Path,
    count_column: str = "unstranded",
    gene_id_type: str = "gene_name",
) -> pd.Series:
    """
    Parse a STAR-Counts TSV file into a Series of gene counts.

    The file format is:
      - Comment line starting with ``#``
      - Header: gene_id, gene_name, gene_type, unstranded,
                stranded_first, stranded_second, tpm_unstranded, ...
      - 4 special QC rows (N_unmapped, etc.)
      - Gene rows

    Memory optimization: reads only the two columns needed (gene ID +
    count column) instead of all ~8 columns, and reads counts directly
    as int32.  For ~60K genes this reduces per-file memory from ~3.8 MB
    to ~0.5 MB (75% savings × 800 files = significant).

    Parameters
    ----------
    filepath : Path
        Path to the downloaded TSV file.
    count_column : str
        Which count column to extract. Default ``"unstranded"``.
    gene_id_type : str
        Which column to use as the gene index: ``"gene_name"``
        or ``"gene_id"``. Default ``"gene_name"``.

    Returns
    -------
    pd.Series
        Index = gene identifiers, values = raw integer counts.
    """
    target_dtype = MEMORY_CONFIG.get("counts_dtype", "int32")

    # ── Determine which columns to read ──────────────────────────
    # gene_id is always needed to filter SPECIAL_ROWS even if
    # gene_id_type is "gene_name" (SPECIAL_ROWS are in gene_id col).
    cols_needed = {gene_id_type, count_column, "gene_id"}

    try:
        df = pd.read_csv(
            filepath,
            sep="\t",
            comment="#",
            usecols=list(cols_needed),
            dtype={count_column: target_dtype},
        )
    except (ValueError, TypeError):
        # Fallback: usecols may fail if columns don't exist or TSV
        # has unexpected format.  Read everything and validate below.
        df = pd.read_csv(filepath, sep="\t", comment="#")

    # Remove special QC rows
    if "gene_id" in df.columns:
        df = df[~df["gene_id"].isin(SPECIAL_ROWS)]

    # Validate requested columns exist
    if count_column not in df.columns:
        available = [c for c in df.columns if c not in ("gene_id", "gene_name", "gene_type")]
        raise ValueError(
            f"Count column '{count_column}' not found. "
            f"Available: {available}"
        )
    if gene_id_type not in df.columns:
        raise ValueError(
            f"Gene ID column '{gene_id_type}' not found. "
            f"Available: {list(df.columns)}"
        )

    # Handle duplicate gene names by summing (common for gene_name)
    series = df.groupby(gene_id_type)[count_column].sum()
    series.name = filepath.stem  # Will be replaced by sample ID later

    return series


# ══════════════════════════════════════════════════════════════════════
# 5. Build count matrix + metadata
# ══════════════════════════════════════════════════════════════════════

def build_count_matrix(
    extracted_files: dict[str, Path],
    file_metadata: pd.DataFrame,
    gene_id_type: str = "gene_name",
    count_column: str = "unstranded",
    progress_callback: Callable[[int, int, str], None] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Assemble a DESeq2-ready count matrix and metadata table.

    Parameters
    ----------
    extracted_files : dict[str, Path]
        Mapping of file_id (UUID) -> path to extracted TSV file.
    file_metadata : pd.DataFrame
        Must have columns: file_id, case_id, sample_type, condition.
        One row per file.
    gene_id_type : str
        ``"gene_name"`` or ``"gene_id"``.
    count_column : str
        Count column to use from the TSV files.
    progress_callback : callable or None
        Called as ``progress_callback(current, total, message)``.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        (counts_df, metadata_df)
        - counts_df: genes (rows) x samples (columns), integer counts.
        - metadata_df: columns = [sample, condition, case_id, sample_type].
    """
    target_dtype = MEMORY_CONFIG.get("counts_dtype", "int32")

    # Build lookup: file_id -> metadata row
    meta_lookup: dict[str, dict] = {}
    for _, row in file_metadata.iterrows():
        meta_lookup[row.get("file_id", "")] = row.to_dict()

    # ── Pass 1: Parse all files into Series + build metadata ──────
    # Each Series is ~240 KB (60K genes × 4 bytes int32).
    # We keep them in a list temporarily to discover the master gene
    # index (union of all gene sets — usually identical across files).
    parsed: list[tuple[str, pd.Series]] = []  # (sample_name, counts)
    meta_rows: list[dict] = []
    total = len(extracted_files)
    existing_names: set[str] = set()

    for i, (file_id, fpath) in enumerate(extracted_files.items(), start=1):
        if progress_callback:
            progress_callback(i, total, f"Parsing {fpath.name}")

        meta = meta_lookup.get(file_id, {})
        case_id = meta.get("case_id", file_id[:8])
        condition = meta.get("condition", "Unknown")
        sample_type = meta.get("sample_type", "Unknown")

        # Build a unique sample name: case_id + condition initial
        sample_name = f"{case_id}_{condition[0] if condition else 'X'}"

        # Ensure uniqueness by appending suffix if needed
        if sample_name in existing_names:
            suffix = 2
            while f"{sample_name}_{suffix}" in existing_names:
                suffix += 1
            sample_name = f"{sample_name}_{suffix}"

        try:
            counts = parse_star_counts(fpath, count_column, gene_id_type)
            parsed.append((sample_name, counts))
            existing_names.add(sample_name)

            meta_rows.append({
                "sample": sample_name,
                "condition": condition,
                "case_id": case_id,
                "sample_type": sample_type,
            })
        except Exception:
            # Skip unparseable files
            continue

    if not parsed:
        raise ValueError("No files could be parsed into count vectors.")

    # ── Pass 2: Build master gene index ───────────────────────────
    # All GDC files from the same project use the same GTF, so gene
    # lists are almost always identical.  We take the union to be
    # defensive (handles rare edge cases like different GTF versions).
    # Cost: a few set operations on ~60K strings — negligible.
    master_genes = parsed[0][1].index
    if len(parsed) > 1:
        # Quick check: if first and last have the same index, skip
        # the expensive union (common case).
        if not parsed[0][1].index.equals(parsed[-1][1].index):
            gene_set: set[str] = set()
            for _, s in parsed:
                gene_set.update(s.index)
            master_genes = pd.Index(sorted(gene_set))

    n_genes = len(master_genes)
    n_samples = len(parsed)
    gene_to_row = {g: i for i, g in enumerate(master_genes)}

    # ── Pass 3: Pre-allocate numpy matrix + fill column by column ─
    # Single allocation: n_genes × n_samples × 4 bytes (int32).
    # For 60K × 800 = ~190 MB — one contiguous block, zero fragmentation.
    # np.zeros pre-fills with 0, so missing genes are already handled.
    matrix = np.zeros((n_genes, n_samples), dtype=target_dtype)
    sample_names: list[str] = []

    for col_idx, (sample_name, series) in enumerate(parsed):
        sample_names.append(sample_name)

        # Fast path: if gene index matches master exactly, direct copy
        if series.index.equals(master_genes):
            matrix[:, col_idx] = series.values.astype(target_dtype)
        else:
            # Slow path: align by gene name (rare for same-project data)
            for gene, val in series.items():
                row_idx = gene_to_row.get(gene)
                if row_idx is not None:
                    matrix[row_idx, col_idx] = val

    # Free all parsed Series — the data is now in the numpy matrix
    del parsed

    # ── Pass 4: Build final DataFrame (zero-copy from numpy) ──────
    counts_df = pd.DataFrame(matrix, index=master_genes, columns=sample_names)
    counts_df.index.name = gene_id_type

    # Build metadata DataFrame
    metadata_df = pd.DataFrame(meta_rows)

    return counts_df, metadata_df
