"""
scrna_pipeline.py -- Backend pipeline for scRNA-seq analysis.

Wraps scanpy functions into a step-by-step pipeline following
the standard Sanbomics/scanpy workflow:

    1. Load data (h5ad, 10x h5, CSV/TSV)
    2. QC annotation (mitochondrial, ribosomal, hemoglobin genes)
    3. Cell & gene filtering
    4. Doublet detection (Scrublet)
    5. Normalization + log1p
    6. Highly variable gene selection
    7. Scaling
    8. PCA
    9. Neighborhood graph
   10. UMAP embedding
   11. Leiden clustering
   12. Marker gene identification

Each function operates on an AnnData object in-place and returns it.
"""

import gc
import gzip
import io
import shutil
import tempfile
from pathlib import Path
from typing import Optional

import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io

from config import SCRNA_CONFIG


# ── LAPACK / scipy condition-number guard ────────────────────────────
# Some LAPACK routines (eigh, inv, solve) raise scipy.linalg.LinAlgError
# with a *bytes* message like b'reciprocal condition number ...'.
# We intercept these and convert them to a clean Python-string warning
# so they never crash the pipeline.

def _is_lapack_condition_error(exc: Exception) -> bool:
    """Return True if *exc* is a LAPACK reciprocal-condition-number error."""
    msg = str(exc)
    return "reciprocal condition number" in msg or b"reciprocal" in str(exc).encode(errors="replace")


# ══════════════════════════════════════════════════════════════════════
# 1. Data Loading (H5AD only)
# ══════════════════════════════════════════════════════════════════════

def load_h5ad(uploaded_file) -> ad.AnnData:
    """
    Load an AnnData object from an uploaded .h5ad file.

    Parameters
    ----------
    uploaded_file : UploadedFile
        Streamlit uploaded file object (.h5ad).

    Returns
    -------
    AnnData
    """
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        tmp.write(uploaded_file.getvalue())
        tmp_path = tmp.name

    try:
        adata = sc.read_h5ad(tmp_path)
    finally:
        Path(tmp_path).unlink(missing_ok=True)

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    return adata


def load_h5ad_backed(file_bytes: bytes) -> tuple:
    """
    Load only metadata (obs, var, shape) from an H5AD file using backed
    (disk) mode.  This avoids reading the full expression matrix into RAM,
    which can save several GB for large datasets.

    Parameters
    ----------
    file_bytes : bytes
        Raw bytes of the .h5ad file.

    Returns
    -------
    tuple
        (obs: pd.DataFrame, var: pd.DataFrame, shape: tuple[int, int],
         tmp_path: str)
        The caller is responsible for deleting tmp_path when no longer needed.
    """
    tmp = tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False)
    tmp.write(file_bytes)
    tmp.close()
    tmp_path = tmp.name

    try:
        adata_b = sc.read_h5ad(tmp_path, backed="r")
        obs = adata_b.obs.copy()
        var = adata_b.var.copy()
        shape = adata_b.shape
        adata_b.file.close()
        del adata_b
        gc.collect()
    except Exception:
        # Fallback: if backed mode fails (e.g. non-chunked HDF5), load fully
        adata = sc.read_h5ad(tmp_path)
        obs = adata.obs.copy()
        var = adata.var.copy()
        shape = adata.shape
        del adata
        gc.collect()

    return obs, var, shape, tmp_path


def load_h5ad_from_path(tmp_path: str) -> ad.AnnData:
    """
    Load a full AnnData object from a temporary file path (created by
    ``load_h5ad_backed``).  The file is deleted after loading.

    Parameters
    ----------
    tmp_path : str
        Path to the .h5ad temporary file.

    Returns
    -------
    AnnData
    """
    try:
        adata = sc.read_h5ad(tmp_path)
    finally:
        Path(tmp_path).unlink(missing_ok=True)

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    return adata


# ══════════════════════════════════════════════════════════════════════
# 1b. 10x File Integrator → H5AD
# ══════════════════════════════════════════════════════════════════════

def _open_maybe_gzipped(file_bytes: bytes, name: str):
    """Return a text-mode file handle, decompressing .gz if needed."""
    if name.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(io.BytesIO(file_bytes), "rb"), encoding="utf-8")
    return io.StringIO(file_bytes.decode("utf-8"))


def _fix_features_columns(
    raw_text: str,
    gene_id_map: Optional[dict] = None,
) -> str:
    """
    Ensure a features/genes file has exactly 3 tab-separated columns
    (gene_id, gene_name, feature_type) as expected by
    ``scanpy.read_10x_mtx()``.

    Handles:
    - **1 column** (gene_name only): looks up gene_id from
      *gene_id_map* if available, otherwise duplicates the name.
      Appends 'Gene Expression'.
    - **2 columns** (gene_id + gene_name): appends 'Gene Expression'.
    - **3+ columns**: returned as-is.

    Parameters
    ----------
    raw_text : str
        Raw text content of features.tsv / genes.tsv.
    gene_id_map : dict, optional
        Mapping ``{gene_name: gene_id}`` extracted from an external
        source (e.g. metadata with Ensembl IDs).  Used only when
        features has 1 column.
    """
    lines = raw_text.rstrip("\n").split("\n")
    if not lines:
        return raw_text

    # Detect number of columns from the first non-empty line
    first_line = lines[0]
    n_cols = len(first_line.split("\t"))

    if n_cols >= 3:
        return raw_text  # already has 3 columns, nothing to do

    fixed = []
    if n_cols == 1:
        # Only gene names — try to resolve gene_id from map
        _map = gene_id_map or {}
        for line in lines:
            name = line.rstrip()
            if name:
                gid = _map.get(name, name)  # fallback: duplicate
                fixed.append(f"{gid}\t{name}\tGene Expression")
    else:
        # 2 columns — just append feature type
        for line in lines:
            fixed.append(line.rstrip() + "\tGene Expression")

    return "\n".join(fixed) + "\n"


def _try_extract_gene_id_map(metadata_file) -> Optional[dict]:
    """
    Try to extract a ``{gene_name: gene_id}`` mapping from a metadata
    file.  Returns ``None`` if no suitable columns are found.

    Detection heuristics:
    - Looks for columns whose values start with Ensembl-like prefixes
      (ENSG, ENSMUSG, ENSDARG, …).
    - If the index looks like gene symbols and a column looks like
      Ensembl IDs (or vice-versa), build the map.
    """
    if metadata_file is None:
        return None

    try:
        raw = metadata_file.getvalue()
        name = metadata_file.name.lower()

        sep = "\t" if (".tsv" in name) else ","

        if name.endswith(".gz"):
            text = gzip.decompress(raw).decode("utf-8")
        else:
            text = raw.decode("utf-8")

        df = pd.read_csv(io.StringIO(text), sep=sep, index_col=0, nrows=50)

        # Ensembl-like pattern: ENS + optional species code + digits
        import re
        _ens_re = re.compile(r"^ENS[A-Z]*G\d{5,}")

        def _looks_ensembl(series: pd.Series) -> bool:
            """Return True if ≥50 % of non-null values match Ensembl."""
            sample = series.dropna().astype(str).head(30)
            if len(sample) == 0:
                return False
            return (sample.str.match(_ens_re).sum() / len(sample)) >= 0.5

        # Check if the index itself is Ensembl
        idx_is_ens = _looks_ensembl(pd.Series(df.index))

        # Find columns that look like Ensembl IDs
        ens_cols = [c for c in df.columns if _looks_ensembl(df[c])]

        # Find columns that look like gene symbols (short, no dots-only)
        def _looks_symbols(series: pd.Series) -> bool:
            sample = series.dropna().astype(str).head(30)
            if len(sample) == 0:
                return False
            # Symbols: typically 1-20 chars, letters/digits/dashes
            sym_pat = re.compile(r"^[A-Za-z][A-Za-z0-9._\-]{0,30}$")
            return (sample.str.match(sym_pat).sum() / len(sample)) >= 0.5

        sym_cols = [c for c in df.columns if _looks_symbols(df[c])]

        # Re-read full file only if we found a promising pair
        if idx_is_ens and sym_cols:
            # index = Ensembl ID, column = gene symbol → map symbol→id
            df_full = pd.read_csv(
                io.StringIO(text), sep=sep, index_col=0,
            )
            metadata_file.seek(0)
            sym_col = sym_cols[0]
            return dict(
                zip(
                    df_full[sym_col].astype(str),
                    df_full.index.astype(str),
                )
            )

        if ens_cols and not idx_is_ens:
            # index = gene symbol, column = Ensembl ID → map symbol→id
            idx_looks_sym = _looks_symbols(pd.Series(df.index))
            if idx_looks_sym:
                df_full = pd.read_csv(
                    io.StringIO(text), sep=sep, index_col=0,
                )
                metadata_file.seek(0)
                ens_col = ens_cols[0]
                return dict(
                    zip(
                        df_full.index.astype(str),
                        df_full[ens_col].astype(str),
                    )
                )

        metadata_file.seek(0)
        return None

    except Exception:
        try:
            metadata_file.seek(0)
        except Exception:
            pass
        return None


def integrate_10x_files(
    matrix_file,
    features_file,
    barcodes_file,
    metadata_file=None,
    metadata_force_target: Optional[str] = None,
) -> bytes:
    """
    Integrate 10x Genomics files (matrix.mtx, features/genes.tsv,
    barcodes.tsv) into a single .h5ad file.

    Handles both plain and .gz compressed files.
    Silently fixes features/genes files with fewer than 3 columns:
    - 1 column (gene names only): tries to extract gene IDs from
      the metadata file (Ensembl IDs), falls back to duplicating.
    - 2 columns (gene_id + gene_name): appends 'Gene Expression'.

    If a metadata CSV/TSV is provided, it is automatically assigned
    to ``adata.obs`` or ``adata.var`` based on index overlap.

    Parameters
    ----------
    matrix_file : UploadedFile
        matrix.mtx or matrix.mtx.gz
    features_file : UploadedFile
        features.tsv(.gz) or genes.tsv(.gz)
    barcodes_file : UploadedFile
        barcodes.tsv(.gz)
    metadata_file : UploadedFile, optional
        CSV or TSV with extra annotations.

    Returns
    -------
    bytes
        The .h5ad file content ready for download or direct loading.
    """
    tmpdir = tempfile.mkdtemp(prefix="beginseq_10x_")
    tmpdir_path = Path(tmpdir)

    try:
        # ── Write matrix.mtx.gz ─────────────────────────────────────
        mtx_bytes = matrix_file.getvalue()
        mtx_name = matrix_file.name

        if mtx_name.endswith(".gz"):
            (tmpdir_path / "matrix.mtx.gz").write_bytes(mtx_bytes)
        else:
            # Compress to .gz so scanpy's read_10x_mtx finds it
            with gzip.open(tmpdir_path / "matrix.mtx.gz", "wb") as gz:
                gz.write(mtx_bytes)

        # ── Write barcodes.tsv.gz ───────────────────────────────────
        bc_bytes = barcodes_file.getvalue()
        bc_name = barcodes_file.name

        if bc_name.endswith(".gz"):
            (tmpdir_path / "barcodes.tsv.gz").write_bytes(bc_bytes)
        else:
            with gzip.open(tmpdir_path / "barcodes.tsv.gz", "wb") as gz:
                gz.write(bc_bytes)

        # ── Write features.tsv.gz (with silent 2→3 column fix) ─────
        feat_bytes = features_file.getvalue()
        feat_name = features_file.name

        # Decompress if needed to inspect columns
        if feat_name.endswith(".gz"):
            with gzip.open(io.BytesIO(feat_bytes), "rt", encoding="utf-8") as f:
                raw_text = f.read()
        else:
            raw_text = feat_bytes.decode("utf-8")

        # If features has 1 column, try to get gene IDs from metadata
        _first_line = raw_text.split("\n", 1)[0]
        _n_feat_cols = len(_first_line.split("\t"))
        gene_id_map = None
        if _n_feat_cols == 1 and metadata_file is not None:
            gene_id_map = _try_extract_gene_id_map(metadata_file)

        # Fix features columns (1→3 or 2→3) silently
        fixed_text = _fix_features_columns(raw_text, gene_id_map=gene_id_map)

        # Write as features.tsv.gz (scanpy expects this name)
        with gzip.open(tmpdir_path / "features.tsv.gz", "wt", encoding="utf-8") as gz:
            gz.write(fixed_text)

        # ── Validate files exist before reading ─────────────────────
        expected = ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"]
        for fname in expected:
            fpath = tmpdir_path / fname
            if not fpath.exists():
                raise FileNotFoundError(
                    f"Expected file '{fname}' not found in temp directory."
                )
            if fpath.stat().st_size == 0:
                raise ValueError(f"File '{fname}' is empty (0 bytes).")

        # ── Read with scanpy ────────────────────────────────────────
        try:
            adata = sc.read_10x_mtx(
                tmpdir_path,
                var_names="gene_symbols",
                cache=False,
            )
        except Exception as read_err:
            import traceback
            tb = traceback.format_exc()
            raise RuntimeError(
                f"scanpy.read_10x_mtx failed: {read_err}\n\n"
                f"Details:\n"
                f"  matrix.mtx.gz  = {(tmpdir_path / 'matrix.mtx.gz').stat().st_size:,} bytes\n"
                f"  features.tsv.gz = {(tmpdir_path / 'features.tsv.gz').stat().st_size:,} bytes\n"
                f"  barcodes.tsv.gz = {(tmpdir_path / 'barcodes.tsv.gz').stat().st_size:,} bytes\n"
                f"\nTraceback:\n{tb}"
            ) from read_err

        adata.var_names_make_unique()
        adata.obs_names_make_unique()

        # ── Merge optional metadata ─────────────────────────────────
        if metadata_file is not None:
            meta_result = detect_and_merge_metadata(
                adata, metadata_file, force_target=metadata_force_target,
            )
            adata = meta_result["adata"]

        # ── Serialise to h5ad bytes ─────────────────────────────────
        h5ad_path = tmpdir_path / "integrated.h5ad"
        adata.write_h5ad(h5ad_path)
        h5ad_bytes = h5ad_path.read_bytes()

        return h5ad_bytes

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ══════════════════════════════════════════════════════════════════════
# 1b-extra. Metadata Auto-Detection & Merge
# ══════════════════════════════════════════════════════════════════════

def _read_metadata_file(uploaded_file) -> pd.DataFrame:
    """Read a CSV or TSV metadata file into a DataFrame, using the
    first column as index."""
    raw = uploaded_file.getvalue()
    name = uploaded_file.name.lower()

    # Detect separator
    if name.endswith(".tsv") or name.endswith(".tsv.gz"):
        sep = "\t"
    else:
        sep = ","

    # Decompress if needed
    if name.endswith(".gz"):
        import gzip as _gz
        text = _gz.decompress(raw).decode("utf-8")
    else:
        text = raw.decode("utf-8")

    df = pd.read_csv(io.StringIO(text), sep=sep, index_col=0)
    return df


def detect_metadata_target(
    adata: ad.AnnData,
    metadata_df: pd.DataFrame,
) -> dict:
    """
    Automatically detect whether a metadata DataFrame should be merged
    into ``adata.obs`` (cell annotations) or ``adata.var`` (gene
    annotations) based on index overlap.

    Parameters
    ----------
    adata : AnnData
        The AnnData object with obs_names (barcodes) and var_names (genes).
    metadata_df : pd.DataFrame
        Metadata with its index set to the first column.

    Returns
    -------
    dict
        {
            "target": "obs" | "var" | "unknown",
            "obs_overlap": int,
            "var_overlap": int,
            "obs_pct": float,   # 0-100
            "var_pct": float,   # 0-100
            "meta_ids": int,    # total IDs in metadata
        }
    """
    meta_ids = set(metadata_df.index.astype(str))
    obs_ids = set(adata.obs_names.astype(str))
    var_ids = set(adata.var_names.astype(str))

    obs_overlap = len(meta_ids & obs_ids)
    var_overlap = len(meta_ids & var_ids)

    n_meta = len(meta_ids)
    obs_pct = (obs_overlap / n_meta * 100) if n_meta > 0 else 0.0
    var_pct = (var_overlap / n_meta * 100) if n_meta > 0 else 0.0

    # Decision: whichever has the higher overlap wins.
    # Require at least 10% overlap to be considered a match.
    if obs_overlap > var_overlap and obs_pct >= 10:
        target = "obs"
    elif var_overlap > obs_overlap and var_pct >= 10:
        target = "var"
    elif obs_overlap > 0 and obs_overlap == var_overlap:
        # Tie-break: if same overlap, prefer obs (more common)
        target = "obs"
    else:
        target = "unknown"

    return {
        "target": target,
        "obs_overlap": obs_overlap,
        "var_overlap": var_overlap,
        "obs_pct": round(obs_pct, 1),
        "var_pct": round(var_pct, 1),
        "meta_ids": n_meta,
    }


def detect_and_merge_metadata(
    adata: ad.AnnData,
    metadata_file,
    force_target: Optional[str] = None,
) -> dict:
    """
    Read a metadata file, auto-detect its target (obs or var),
    and merge it into the AnnData object.

    Parameters
    ----------
    adata : AnnData
        The AnnData to annotate.
    metadata_file : UploadedFile
        Streamlit uploaded CSV/TSV file.
    force_target : str, optional
        If ``"obs"`` or ``"var"``, skip auto-detection and merge
        directly into the specified slot.

    Returns
    -------
    dict
        {
            "adata": AnnData (modified in-place and returned),
            "target": "obs" | "var" | "unknown",
            "detection": dict from detect_metadata_target,
            "columns_added": list[str],
        }
    """
    meta_df = _read_metadata_file(metadata_file)
    meta_df.index = meta_df.index.astype(str)

    detection = detect_metadata_target(adata, meta_df)
    target = force_target if force_target in ("obs", "var") else detection["target"]

    columns_added = []

    if target == "obs":
        # Align to adata.obs_names — left join so all cells are kept
        for col in meta_df.columns:
            safe_col = col if col not in adata.obs.columns else f"meta_{col}"
            adata.obs[safe_col] = (
                adata.obs_names.astype(str)
                .to_series()
                .map(meta_df[col])
                .values
            )
            columns_added.append(safe_col)

    elif target == "var":
        # Align to adata.var_names — left join so all genes are kept
        for col in meta_df.columns:
            safe_col = col if col not in adata.var.columns else f"meta_{col}"
            adata.var[safe_col] = (
                adata.var_names.astype(str)
                .to_series()
                .map(meta_df[col])
                .values
            )
            columns_added.append(safe_col)

    return {
        "adata": adata,
        "target": target,
        "detection": detection,
        "columns_added": columns_added,
    }


# ══════════════════════════════════════════════════════════════════════
# 1c. Ambient RNA Removal — SoupX via rpy2
# ══════════════════════════════════════════════════════════════════════

def check_soupx_available() -> dict:
    """
    Check whether R, rpy2, and SoupX are available.

    Returns
    -------
    dict  {"r": bool, "rpy2": bool, "soupx": bool, "error": str|None}
    """
    import warnings

    result = {"r": False, "rpy2": False, "soupx": False, "error": None}

    # 1. Check rpy2 import
    try:
        import rpy2  # noqa: F401
        result["rpy2"] = True
    except ImportError:
        result["error"] = "rpy2 is not installed. Run: pip install rpy2"
        return result

    # 2. Check R is reachable
    # Suppress rpy2 warnings about R not being initialized from the main
    # thread and the reticulate segfault notice — they are harmless and
    # would pollute the terminal on every Streamlit rerun.
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="R is not initialized")
            warnings.filterwarnings("ignore", message=".*reticulate.*")
            import rpy2.robjects  # noqa: F401
        result["r"] = True
    except Exception as exc:
        result["error"] = f"R not found or rpy2 cannot connect to R: {exc}"
        return result

    # 3. Check SoupX package
    try:
        from rpy2.robjects.packages import importr
        importr("SoupX")
        result["soupx"] = True
    except Exception:
        result["error"] = (
            "SoupX R package not installed. "
            "Open R and run: install.packages('SoupX')"
        )

    return result


def run_soupx(
    raw_h5ad_bytes: bytes,
    filtered_h5ad_bytes: bytes,
    contamination_fraction: Optional[float] = None,
) -> bytes:
    """
    Remove ambient RNA from a filtered count matrix using SoupX (via rpy2).

    SoupX estimates the contamination profile from the raw (unfiltered)
    droplets and removes it from the filtered (cell-containing) droplets.

    Parameters
    ----------
    raw_h5ad_bytes : bytes
        Unfiltered (raw) count matrix as .h5ad bytes.
    filtered_h5ad_bytes : bytes
        Filtered count matrix as .h5ad bytes.
    contamination_fraction : float, optional
        If provided, override SoupX's automatic contamination estimation.
        Typical range: 0.01 – 0.20.

    Returns
    -------
    bytes
        Cleaned .h5ad file content.
    """
    # ── Lazy imports (rpy2 is optional) ──────────────────────────────
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.packages import importr

    numpy2ri.activate()
    pandas2ri.activate()

    SoupX = importr("SoupX")
    Matrix = importr("Matrix")
    base = importr("base")

    tmpdir = tempfile.mkdtemp(prefix="beginseq_soupx_")
    tmpdir_path = Path(tmpdir)

    try:
        # ── Load both AnnData objects ────────────────────────────────
        raw_path = tmpdir_path / "raw.h5ad"
        filt_path = tmpdir_path / "filtered.h5ad"
        raw_path.write_bytes(raw_h5ad_bytes)
        filt_path.write_bytes(filtered_h5ad_bytes)

        adata_raw = sc.read_h5ad(str(raw_path))
        adata_filt = sc.read_h5ad(str(filt_path))

        adata_raw.var_names_make_unique()
        adata_raw.obs_names_make_unique()
        adata_filt.var_names_make_unique()
        adata_filt.obs_names_make_unique()

        # ── Quick clustering on filtered data (SoupX needs clusters) ─
        adata_tmp = adata_filt.copy()
        sc.pp.normalize_total(adata_tmp)
        sc.pp.log1p(adata_tmp)
        sc.pp.highly_variable_genes(adata_tmp, n_top_genes=2000, flavor="seurat_v3", layer=None)
        sc.pp.pca(adata_tmp, n_comps=min(30, adata_tmp.n_vars - 1, adata_tmp.n_obs - 1))
        sc.pp.neighbors(adata_tmp, n_neighbors=15)
        sc.tl.leiden(adata_tmp, resolution=0.5, flavor="igraph", n_iterations=2)

        clusters = adata_tmp.obs["leiden"].values.astype(str)

        # ── Align gene space between raw & filtered ──────────────────
        common_genes = list(
            set(adata_raw.var_names).intersection(set(adata_filt.var_names))
        )
        common_genes.sort()

        adata_raw = adata_raw[:, common_genes].copy()
        adata_filt = adata_filt[:, common_genes].copy()

        # ── Convert to dense numpy arrays for R ──────────────────────
        import scipy.sparse as sp

        raw_X = adata_raw.X
        if sp.issparse(raw_X):
            raw_X = raw_X.toarray()
        raw_X = np.asarray(raw_X, dtype=np.float64)

        filt_X = adata_filt.X
        if sp.issparse(filt_X):
            filt_X = filt_X.toarray()
        filt_X = np.asarray(filt_X, dtype=np.float64)

        # ── Prepare R objects ────────────────────────────────────────
        # SoupX expects genes-as-rows, cells-as-columns  (dgCMatrix)
        # We have cells-as-rows, genes-as-columns (AnnData convention)
        # So we transpose.
        r_tod = ro.r("function(m) as(t(m), 'dgCMatrix')")

        r_raw_mat = r_tod(ro.r.matrix(
            ro.FloatVector(raw_X.flatten()),
            nrow=raw_X.shape[0],
            ncol=raw_X.shape[1],
        ))
        r_filt_mat = r_tod(ro.r.matrix(
            ro.FloatVector(filt_X.flatten()),
            nrow=filt_X.shape[0],
            ncol=filt_X.shape[1],
        ))

        # Set dimnames  (genes = rownames, cells = colnames)
        gene_names = ro.StrVector(common_genes)
        raw_barcodes = ro.StrVector(list(adata_raw.obs_names))
        filt_barcodes = ro.StrVector(list(adata_filt.obs_names))

        ro.r.assign("raw_mat", r_raw_mat)
        ro.r.assign("filt_mat", r_filt_mat)
        ro.r.assign("gene_names", gene_names)
        ro.r.assign("raw_bc", raw_barcodes)
        ro.r.assign("filt_bc", filt_barcodes)

        ro.r("""
            rownames(raw_mat) <- gene_names
            colnames(raw_mat) <- raw_bc
            rownames(filt_mat) <- gene_names
            colnames(filt_mat) <- filt_bc
        """)

        # ── Build SoupChannel ────────────────────────────────────────
        r_clusters = ro.StrVector(list(clusters))
        ro.r.assign("clusters", r_clusters)
        ro.r.assign("filt_bc_vec", filt_barcodes)

        ro.r("""
            cluster_df <- data.frame(
                row.names = filt_bc_vec,
                clusters  = clusters
            )
        """)

        ro.r("sc_obj <- SoupX::SoupChannel(raw_mat, filt_mat)")
        ro.r("sc_obj <- SoupX::setClusters(sc_obj, cluster_df$clusters)")

        # ── Estimate / set contamination & adjust counts ─────────────
        if contamination_fraction is not None:
            ro.r.assign("contam", ro.FloatVector([contamination_fraction]))
            ro.r("sc_obj <- SoupX::setContaminationFraction(sc_obj, contam)")
            ro.r("adjusted <- SoupX::adjustCounts(sc_obj)")
        else:
            ro.r("sc_obj <- SoupX::autoEstCont(sc_obj)")
            ro.r("adjusted <- SoupX::adjustCounts(sc_obj)")

        # ── Pull corrected matrix back into Python ───────────────────
        ro.r("adj_dense <- as.matrix(adjusted)")
        adj_matrix = np.array(ro.r("adj_dense"))  # genes x cells
        adj_matrix = adj_matrix.T  # back to cells x genes

        # ── Build cleaned AnnData ────────────────────────────────────
        adata_clean = ad.AnnData(
            X=adj_matrix.astype(np.float32),
            obs=adata_filt.obs.copy(),
            var=adata_filt[:, common_genes].var.copy(),
        )
        adata_clean.var_names = common_genes
        adata_clean.obs_names = list(adata_filt.obs_names)

        # Store metadata about the SoupX run
        adata_clean.uns["soupx"] = {
            "method": "SoupX",
            "contamination_fraction": (
                contamination_fraction
                if contamination_fraction is not None
                else "auto"
            ),
            "n_raw_droplets": adata_raw.n_obs,
            "n_filtered_cells": adata_filt.n_obs,
            "n_genes_used": len(common_genes),
        }

        # ── Serialise to h5ad bytes ──────────────────────────────────
        out_path = tmpdir_path / "soupx_cleaned.h5ad"
        adata_clean.write_h5ad(out_path)
        return out_path.read_bytes()

    finally:
        # Deactivate rpy2 converters
        try:
            numpy2ri.deactivate()
            pandas2ri.deactivate()
        except Exception:
            pass
        shutil.rmtree(tmpdir, ignore_errors=True)
        gc.collect()


# ══════════════════════════════════════════════════════════════════════
# 2. Quality Control Annotation
# ══════════════════════════════════════════════════════════════════════

def annotate_qc(adata: ad.AnnData) -> ad.AnnData:
    """
    Annotate mitochondrial, ribosomal, and hemoglobin genes,
    then compute QC metrics.
    """
    # Gene category flags
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.upper().str.contains("^HB[^(P)]", regex=True)

    # Compute QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        inplace=True,
        log1p=True,
    )

    return adata


# ══════════════════════════════════════════════════════════════════════
# 3. Cell & Gene Filtering
# ══════════════════════════════════════════════════════════════════════

def filter_cells_and_genes(
    adata: ad.AnnData,
    min_genes: int = 200,
    max_genes: int = 5000,
    min_counts: int = 500,
    max_counts: int = 50000,
    max_pct_mt: float = 20.0,
    min_cells: int = 3,
) -> ad.AnnData:
    """
    Filter cells by QC metrics and genes by minimum cell count.

    Returns the filtered AnnData (subset, not in-place).
    """
    n_before_cells = adata.n_obs
    n_before_genes = adata.n_vars

    # Filter genes first
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # Filter cells by gene count
    sc.pp.filter_cells(adata, min_genes=min_genes)

    # Build a single boolean mask for all cell-level filters, then apply
    # once.  This avoids creating 4 intermediate AnnData copies which
    # spike peak memory (each .copy() duplicates the full X matrix).
    mask = adata.obs["n_genes_by_counts"] <= max_genes
    mask &= adata.obs["total_counts"] >= min_counts
    mask &= adata.obs["total_counts"] <= max_counts
    if "pct_counts_mt" in adata.obs.columns:
        mask &= adata.obs["pct_counts_mt"] < max_pct_mt
    adata = adata[mask, :].copy()

    n_after_cells = adata.n_obs
    n_after_genes = adata.n_vars

    adata.uns["filtering_stats"] = {
        "cells_before": n_before_cells,
        "cells_after": n_after_cells,
        "cells_removed": n_before_cells - n_after_cells,
        "genes_before": n_before_genes,
        "genes_after": n_after_genes,
        "genes_removed": n_before_genes - n_after_genes,
    }

    return adata


# ══════════════════════════════════════════════════════════════════════
# 4. Doublet Detection (Scrublet)
# ══════════════════════════════════════════════════════════════════════

def detect_doublets(
    adata: ad.AnnData,
    batch_key: Optional[str] = None,
) -> ad.AnnData:
    """
    Run Scrublet for doublet detection.
    Adds 'predicted_doublet' and 'doublet_score' to adata.obs.
    """
    try:
        sc.pp.scrublet(adata, batch_key=batch_key)
        adata.uns["doublet_stats"] = {
            "n_doublets": int(adata.obs["predicted_doublet"].sum()),
            "n_singlets": int((~adata.obs["predicted_doublet"]).sum()),
            "doublet_rate": float(adata.obs["predicted_doublet"].mean()),
        }
    except Exception:
        # Scrublet can fail on very small datasets
        adata.obs["predicted_doublet"] = False
        adata.obs["doublet_score"] = 0.0
        adata.uns["doublet_stats"] = {
            "n_doublets": 0,
            "n_singlets": adata.n_obs,
            "doublet_rate": 0.0,
        }

    return adata


def remove_doublets(adata: ad.AnnData) -> ad.AnnData:
    """Remove cells predicted as doublets."""
    if "predicted_doublet" in adata.obs.columns:
        adata = adata[~adata.obs["predicted_doublet"]].copy()
    return adata


# ══════════════════════════════════════════════════════════════════════
# 5. Normalization
# ══════════════════════════════════════════════════════════════════════

def normalize_data(
    adata: ad.AnnData,
    target_sum: Optional[float] = None,
) -> ad.AnnData:
    """
    Save raw counts to a layer, then normalize and log-transform.

    For sparse matrices the copy is cheap (only the data array is
    duplicated, not the full dense expansion).
    """
    from scipy import sparse

    if sparse.issparse(adata.X):
        # Sparse copy: only duplicates the .data array (~5% of dense size)
        adata.layers["counts"] = adata.X.copy()
    else:
        adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=target_sum, inplace=True)
    sc.pp.log1p(adata)

    return adata


# ══════════════════════════════════════════════════════════════════════
# 6. Highly Variable Gene Selection
# ══════════════════════════════════════════════════════════════════════

def select_hvg(
    adata: ad.AnnData,
    n_top_genes: int = 2000,
    flavor: str = "seurat_v3",
    batch_key: Optional[str] = None,
) -> ad.AnnData:
    """
    Select highly variable genes.

    For flavor='seurat_v3', uses the raw counts layer.
    Otherwise uses the current (log-normalized) X.
    """
    _kw = dict(n_top_genes=n_top_genes, flavor=flavor, batch_key=batch_key)
    if flavor == "seurat_v3":
        _kw["layer"] = "counts"

    try:
        sc.pp.highly_variable_genes(adata, **_kw)
    except Exception as exc:
        if _is_lapack_condition_error(exc) and batch_key is not None:
            # Retry without batch correction (avoids batch-wise regression)
            warnings.warn("HVG selection hit LAPACK error; retrying without batch_key.")
            _kw.pop("batch_key", None)
            sc.pp.highly_variable_genes(adata, **_kw)
        else:
            raise

    adata.uns["hvg_stats"] = {
        "n_hvg": int(adata.var["highly_variable"].sum()),
        "n_total": adata.n_vars,
    }

    return adata


# ══════════════════════════════════════════════════════════════════════
# 7. Scaling
# ══════════════════════════════════════════════════════════════════════

def scale_data(adata: ad.AnnData, max_value: float = 10.0) -> ad.AnnData:
    """
    Subset to highly variable genes and scale to unit variance.

    Subsetting before scaling avoids densifying the full gene matrix,
    which would consume too much RAM on large datasets.
    """
    # Keep all genes accessible in adata.raw before subsetting
    adata.raw = adata

    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=max_value, zero_center=False)
    return adata


# ══════════════════════════════════════════════════════════════════════
# 8. PCA
# ══════════════════════════════════════════════════════════════════════

def run_pca(adata: ad.AnnData, n_comps: int = 50) -> ad.AnnData:
    """Run PCA on the scaled data (uses HVGs by default)."""
    n_comps = min(n_comps, adata.n_vars - 1, adata.n_obs - 1)
    try:
        sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack")
    except Exception as exc:
        if _is_lapack_condition_error(exc):
            # Fallback: default solver (randomized SVD)
            warnings.warn("PCA arpack hit LAPACK condition error; retrying with default solver.")
            sc.tl.pca(adata, n_comps=n_comps)
        else:
            raise
    return adata


# ══════════════════════════════════════════════════════════════════════
# 8b. Batch Effect Correction (Harmony)
# ══════════════════════════════════════════════════════════════════════


def _harmony_worker(pca_matrix, batch_series, conn):
    """Run Harmony inside an isolated subprocess.

    PyTorch + Harmony allocate large temporary buffers.  By running in a
    child process all that memory is returned to the OS when the process
    exits, preventing the parent from accumulating a multi-GB RSS.

    The corrected matrix is sent back through a ``multiprocessing.Connection``.
    """
    try:
        import numpy as _np
        import pandas as _pd
        import harmonypy as hm

        ho = hm.run_harmony(
            pca_matrix,
            _pd.DataFrame({"batch": batch_series}),
            "batch",
            verbose=False,
        )
        # v0.2: Z_corr is (n_cells, n_pcs)
        # v0.0.x: Z_corr was (n_pcs, n_cells)
        result = _np.asarray(ho.Z_corr, dtype=_np.float32)
        if result.shape[0] != pca_matrix.shape[0]:
            result = result.T  # transpose if old API
        conn.send(("ok", result))
    except Exception as e:
        conn.send(("error", str(e)))
    finally:
        conn.close()


def run_harmony(
    adata: ad.AnnData,
    batch_key: str,
    n_pcs: int = 50,
) -> ad.AnnData:
    """
    Correct batch effects in PCA space using Harmony.

    Runs Harmony in an **isolated child process** so that the PyTorch
    runtime and all temporary allocations are fully released when the
    subprocess exits.  Only the corrected PCA matrix (float32) is sent
    back to the parent, keeping memory usage bounded.

    Parameters
    ----------
    adata : AnnData
        Must have ``X_pca`` computed.
    batch_key : str
        Column in ``adata.obs`` identifying the batch variable.
    n_pcs : int
        Number of PCs to use from ``X_pca``.

    Returns
    -------
    AnnData with ``obsm["X_pca_harmony"]`` and ``uns["harmony_stats"]``.
    """
    import multiprocessing as mp

    if batch_key not in adata.obs.columns:
        raise ValueError(
            f"Batch key '{batch_key}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    n_pcs = min(n_pcs, adata.obsm["X_pca"].shape[1])
    pca_matrix = np.asarray(adata.obsm["X_pca"][:, :n_pcs], dtype=np.float32)
    batch_series = adata.obs[batch_key].values

    # Run in a child process — all torch/harmony RAM is freed on exit
    parent_conn, child_conn = mp.Pipe(duplex=False)
    ctx = mp.get_context("spawn")       # clean process, no fork-bomb
    proc = ctx.Process(
        target=_harmony_worker,
        args=(pca_matrix, batch_series, child_conn),
        daemon=True,
    )
    proc.start()
    child_conn.close()                  # parent doesn't write

    # Wait for result (timeout 10 min for very large datasets)
    if parent_conn.poll(timeout=600):
        status, payload = parent_conn.recv()
    else:
        proc.kill()
        raise RuntimeError("Harmony subprocess timed out after 10 minutes.")
    parent_conn.close()
    proc.join(timeout=10)

    if status == "error":
        raise RuntimeError(f"Harmony failed: {payload}")

    adata.obsm["X_pca_harmony"] = payload  # (n_cells, n_pcs) float32

    n_batches = adata.obs[batch_key].nunique()
    adata.uns["harmony_stats"] = {
        "batch_key": batch_key,
        "n_batches": int(n_batches),
        "n_pcs": n_pcs,
    }
    return adata


# ══════════════════════════════════════════════════════════════════════
# 9. Neighborhood Graph
# ══════════════════════════════════════════════════════════════════════

def compute_neighbors(
    adata: ad.AnnData,
    n_neighbors: int = 15,
    n_pcs: int = 40,
    use_rep: Optional[str] = None,
) -> ad.AnnData:
    """Build the k-nearest-neighbor graph from PCA (or corrected) coordinates."""
    try:
        if use_rep is not None:
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep)
        else:
            n_pcs = min(n_pcs, adata.obsm["X_pca"].shape[1])
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    except Exception as exc:
        if _is_lapack_condition_error(exc):
            # Retry with fewer PCs which often avoids the degenerate subspace
            _safe_pcs = min(n_pcs, 20)
            warnings.warn(f"Neighbors hit LAPACK condition error; retrying with {_safe_pcs} PCs.")
            if use_rep is not None:
                sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=_safe_pcs, use_rep=use_rep)
            else:
                sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=_safe_pcs)
        else:
            raise
    return adata


# ══════════════════════════════════════════════════════════════════════
# 10. UMAP
# ══════════════════════════════════════════════════════════════════════

def run_umap(adata: ad.AnnData, min_dist: float = 0.5, spread: float = 1.0) -> ad.AnnData:
    """Compute UMAP embedding (with fixed seed for reproducibility)."""
    sc.tl.umap(adata, min_dist=min_dist, spread=spread, random_state=42)
    return adata


# ══════════════════════════════════════════════════════════════════════
# 11. Clustering
# ══════════════════════════════════════════════════════════════════════

def run_leiden(
    adata: ad.AnnData,
    resolution: float = 0.5,
    key_added: str = "leiden",
) -> ad.AnnData:
    """Run Leiden clustering."""
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=key_added,
        flavor="igraph",
        n_iterations=2,
    )
    adata.uns["leiden_stats"] = {
        "n_clusters": int(adata.obs[key_added].nunique()),
        "resolution": resolution,
    }
    return adata


# ══════════════════════════════════════════════════════════════════════
# 12. Marker Gene Identification
# ══════════════════════════════════════════════════════════════════════

def find_marker_genes(
    adata: ad.AnnData,
    groupby: str = "leiden",
    method: str = "wilcoxon",
    n_genes: int = 25,
) -> ad.AnnData:
    """
    Rank genes for characterizing each cluster.
    Uses the log-normalized expression (not scaled).
    """
    try:
        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            method=method,
            n_genes=n_genes,
            use_raw=True,
        )
    except Exception as exc:
        if _is_lapack_condition_error(exc):
            # Fallback: try t-test which doesn't require matrix inversion
            _fallback = "t-test" if method != "t-test" else "t-test_overestim_var"
            warnings.warn(f"Marker genes ({method}) hit LAPACK error; retrying with {_fallback}.")
            sc.tl.rank_genes_groups(
                adata,
                groupby=groupby,
                method=_fallback,
                n_genes=n_genes,
                use_raw=True,
            )
        else:
            raise
    return adata


def get_marker_genes_df(
    adata: ad.AnnData,
    group: Optional[str] = None,
    n_genes: int = 25,
) -> pd.DataFrame:
    """
    Extract marker genes as a DataFrame.

    If group is None, returns markers for all groups.
    """
    if group is not None:
        return sc.get.rank_genes_groups_df(adata, group=group).head(n_genes)

    # All groups — recover groupby from rank_genes_groups params
    _groupby = adata.uns["rank_genes_groups"]["params"]["groupby"]
    groups = list(adata.obs[_groupby].cat.categories)
    dfs = []
    for g in groups:
        df = sc.get.rank_genes_groups_df(adata, group=g).head(n_genes)
        df["cluster"] = g
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


# ══════════════════════════════════════════════════════════════════════
# Full Pipeline Orchestrator
# ══════════════════════════════════════════════════════════════════════

PIPELINE_STEPS = [
    "qc_annotation",
    "filtering",
    "doublet_detection",
    "doublet_removal",
    "normalization",
    "hvg_selection",
    "scaling",
    "pca",
    "batch_correction",
    "neighbors",
    "umap",
    "clustering",
    "marker_genes",
]


def run_scrna_pipeline(
    adata: ad.AnnData,
    params: dict,
    progress_callback=None,
) -> ad.AnnData:
    """
    Run the full scRNA-seq analysis pipeline.

    Parameters
    ----------
    adata : AnnData
        Raw loaded data.
    params : dict
        Pipeline parameters (merged with SCRNA_CONFIG defaults).
    progress_callback : callable, optional
        Called with (step_index, total_steps, step_name) for progress tracking.

    Returns
    -------
    AnnData with all analysis results.
    """
    total = len(PIPELINE_STEPS)
    cfg = SCRNA_CONFIG

    def _progress(i, name):
        if progress_callback:
            progress_callback(i, total, name)

    # Step 1: QC annotation
    _progress(0, "qc_annotation")
    adata = annotate_qc(adata)

    # Step 2: Filtering
    _progress(1, "filtering")
    qc = params.get("qc", cfg["qc_defaults"])
    adata = filter_cells_and_genes(
        adata,
        min_genes=qc.get("min_genes", cfg["qc_defaults"]["min_genes"]),
        max_genes=qc.get("max_genes", cfg["qc_defaults"]["max_genes"]),
        min_counts=qc.get("min_counts", cfg["qc_defaults"]["min_counts"]),
        max_counts=qc.get("max_counts", cfg["qc_defaults"]["max_counts"]),
        max_pct_mt=qc.get("max_pct_mt", cfg["qc_defaults"]["max_pct_mt"]),
        min_cells=qc.get("min_cells", cfg["qc_defaults"]["min_cells"]),
    )
    gc.collect()

    # Step 3: Doublet detection
    _progress(2, "doublet_detection")
    adata = detect_doublets(adata, batch_key=params.get("batch_key"))

    # Step 4: Remove doublets
    _progress(3, "doublet_removal")
    if params.get("remove_doublets", True):
        adata = remove_doublets(adata)
    gc.collect()

    # Step 5: Normalization
    _progress(4, "normalization")
    adata = normalize_data(
        adata,
        target_sum=params.get("target_sum", cfg["normalization"]["target_sum"]),
    )
    gc.collect()

    # Step 6: HVG selection
    _progress(5, "hvg_selection")
    hvg_cfg = cfg["hvg"]
    adata = select_hvg(
        adata,
        n_top_genes=params.get("n_top_genes", hvg_cfg["n_top_genes"]),
        flavor=params.get("hvg_flavor", hvg_cfg["flavor"]),
        batch_key=params.get("batch_key"),
    )

    # Step 7: Scaling (subsets to HVGs → drops non-HVG genes from X)
    _progress(6, "scaling")
    adata = scale_data(adata)
    gc.collect()

    # Step 8: PCA
    _progress(7, "pca")
    adata = run_pca(
        adata,
        n_comps=params.get("n_pcs_compute", cfg["pca"]["n_comps"]),
    )
    gc.collect()

    # Step 9: Batch correction (optional — Harmony)
    _progress(8, "batch_correction")
    _use_rep = None
    _batch_key = params.get("batch_key")
    _enable_batch = params.get("enable_batch_correction", False)
    if _enable_batch and _batch_key:
        adata = run_harmony(
            adata,
            batch_key=_batch_key,
            n_pcs=params.get("n_pcs", cfg["neighbors"]["n_pcs"]),
        )
        _use_rep = "X_pca_harmony"
        gc.collect()

    # Step 10: Neighbors
    _progress(9, "neighbors")
    nb_cfg = cfg["neighbors"]
    adata = compute_neighbors(
        adata,
        n_neighbors=params.get("n_neighbors", nb_cfg["n_neighbors"]),
        n_pcs=params.get("n_pcs", nb_cfg["n_pcs"]),
        use_rep=_use_rep,
    )

    # Step 11: UMAP
    _progress(10, "umap")
    adata = run_umap(
        adata,
        min_dist=params.get("umap_min_dist", cfg["umap"]["min_dist"]),
    )
    gc.collect()

    # Step 12: Clustering
    _progress(11, "clustering")
    adata = run_leiden(
        adata,
        resolution=params.get("leiden_resolution", cfg["leiden"]["resolution"]),
    )

    # Step 13: Marker genes
    _progress(12, "marker_genes")
    rg_cfg = cfg["rank_genes"]
    _groupby = params.get("marker_groupby", "leiden")
    adata = find_marker_genes(
        adata,
        groupby=_groupby,
        method=params.get("de_method", rg_cfg["method"]),
        n_genes=params.get("n_marker_genes", rg_cfg["n_genes"]),
    )

    gc.collect()
    return adata
