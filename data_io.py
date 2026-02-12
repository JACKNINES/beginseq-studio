"""
data_io.py — Reading and writing data files.

This module centralizes ALL of the application's I/O (Input/Output) logic.
It fixes the original bug where the separator was hardcoded as tab ("\\t"),
causing CSV files with commas to not be read correctly.

It now automatically detects the separator based on the file extension,
and supports ZIP-compressed files (e.g. metadata.csv.zip).

Functions
---------
detect_separator(filename)
    → Detects the correct separator based on the file extension.

extract_from_zip(uploaded_file)
    → If the file is a ZIP, extracts the first CSV/TSV contained within
      and returns it as a BytesIO with the real name of the internal file.

read_counts_file(uploaded_file)
    → Reads the counts matrix from a Streamlit UploadedFile.

read_metadata_file(uploaded_file)
    → Reads the metadata file from a Streamlit UploadedFile.

results_to_csv(results_df)
    → Converts the results DataFrame to a CSV string for download.

Usage example
--------------
    from data_io import read_counts_file, read_metadata_file

    counts_df = read_counts_file(uploaded_counts_file)
    metadata_df = read_metadata_file(uploaded_metadata_file)
"""

import io
import zipfile

import numpy as np
import pandas as pd
from streamlit.runtime.uploaded_file_manager import UploadedFile

from config import FILE_CONFIG, MEMORY_CONFIG


def detect_separator(filename: str) -> str:
    """
    Detects the column separator based on the file extension.

    Logic:
    - .tsv → tab ("\\t")
    - .csv → comma (",")
    - Other extension → raises ValueError

    Parameters
    ----------
    filename : str
        File name (e.g. "counts.tsv", "data.csv").

    Returns
    -------
    str
        The separator character to use with pd.read_csv().

    Raises
    ------
    ValueError
        If the extension is not recognized (not in FILE_CONFIG).

    Example
    -------
        >>> detect_separator("my_counts.tsv")
        '\\t'
        >>> detect_separator("metadata.csv")
        ','
    """
    extension = filename.rsplit(".", maxsplit=1)[-1].lower()
    separators = FILE_CONFIG["separators"]

    if extension in separators:
        return separators[extension]

    raise ValueError(
        f"Extension '.{extension}' not recognized. "
        f"Supported extensions: {list(separators.keys())}."
    )


def extract_from_zip(uploaded_file: UploadedFile) -> tuple[io.BytesIO, str]:
    """
    If the uploaded file is a ZIP, extracts the first CSV or TSV file
    found inside and returns it as a BytesIO ready to read.

    Supports names such as:
    - metadata.csv.zip  → contains metadata.csv
    - counts.tsv.zip    → contains counts.tsv
    - datos.zip         → contains any .csv or .tsv inside

    Parameters
    ----------
    uploaded_file : UploadedFile
        Streamlit UploadedFile object.

    Returns
    -------
    tuple[io.BytesIO, str]
        - BytesIO with the content of the extracted file.
        - Name of the internal file (e.g. "metadata.csv").

    Raises
    ------
    ValueError
        If the ZIP is corrupt, empty, or does not contain CSV/TSV files.
    """
    uploaded_file.seek(0)

    try:
        with zipfile.ZipFile(uploaded_file, "r") as zf:
            # List files inside the ZIP
            names = zf.namelist()

            if not names:
                raise ValueError(
                    "The ZIP file is empty. "
                    "It must contain at least one CSV or TSV file."
                )

            # Find the first file with a CSV or TSV extension
            valid_extensions = tuple(FILE_CONFIG["separators"].keys())
            target = None
            for name in names:
                if name.lower().endswith(valid_extensions):
                    target = name
                    break

            if target is None:
                raise ValueError(
                    f"The ZIP file does not contain files with valid extensions "
                    f"({valid_extensions}). "
                    f"Files found: {names}."
                )

            # Extract the file into memory
            data = zf.read(target)
            buffer = io.BytesIO(data)
            return buffer, target

    except zipfile.BadZipFile:
        raise ValueError(
            f"The file '{uploaded_file.name}' is not a valid ZIP or "
            f"is corrupt."
        )


def _prepare_file(
    uploaded_file: UploadedFile,
) -> tuple[io.BytesIO | UploadedFile, str]:
    """
    Prepares the file for reading: if it is a ZIP it extracts it,
    otherwise it simply resets the cursor.

    Parameters
    ----------
    uploaded_file : UploadedFile
        File uploaded by the user.

    Returns
    -------
    tuple
        - The file-like object ready to read (BytesIO or UploadedFile).
        - The real file name (the internal one if it came from a ZIP).
    """
    filename = uploaded_file.name

    if filename.lower().endswith(".zip"):
        buffer, inner_name = extract_from_zip(uploaded_file)
        return buffer, inner_name
    else:
        uploaded_file.seek(0)
        return uploaded_file, filename


def read_counts_file(uploaded_file: UploadedFile) -> pd.DataFrame:
    """
    Reads a raw counts matrix from a file uploaded by the user.
    Supports CSV, TSV, or ZIP files containing a CSV/TSV.

    The file must have:
    - First column: gene identifiers (used as the index).
    - Remaining columns: samples with raw integer counts.

    Expected format:
        gene_id     sample_1    sample_2    sample_3
        GENE_A      120         340         256
        GENE_B      0           15          8
        ...

    Parameters
    ----------
    uploaded_file : UploadedFile
        Streamlit UploadedFile object (returned by st.file_uploader).

    Returns
    -------
    pd.DataFrame
        DataFrame with genes as rows and samples as columns.

    Raises
    ------
    ValueError
        If the file cannot be read or the extension is not valid.
    """
    file_obj, real_name = _prepare_file(uploaded_file)
    sep = detect_separator(real_name)
    target_dtype = MEMORY_CONFIG.get("counts_dtype", "int32")

    try:
        # ── Memory-critical read ───────────────────────────────────
        # Default pd.read_csv uses float64 (8 bytes/value).  For a
        # typical TCGA matrix (60K genes × 800 samples = 48M values)
        # that's ~3.8 GB.  Reading directly as int32 cuts this to
        # ~190 MB — a 20× reduction.
        #
        # comment="#" skips STAR-Counts header lines (e.g. "# ...")
        # that are common in GDC/TCGA augmented STAR gene count TSVs.
        #
        # IMPORTANT: We cannot pass dtype=int32 as a scalar because
        # that also applies to the index column (gene IDs = strings)
        # and pd.read_csv will raise.  Instead we:
        #   1. Peek at the header to discover column names.
        #   2. Build a dtype dict mapping only data columns to int32.
        #   3. Read the full file with that dict + index_col=0.
        # If the data has floats (FPKM/TPM), the int32 read raises
        # and we fall back to float32 (still 50% less than float64).
        #
        # Peek at header line (skip # comment lines)
        start_pos = file_obj.tell()
        header_line = ""
        for raw_line in file_obj:
            line = raw_line.decode("utf-8") if isinstance(raw_line, bytes) else raw_line
            if not line.startswith("#"):
                header_line = line.strip()
                break
        file_obj.seek(start_pos)  # rewind for pd.read_csv

        if header_line:
            col_names = header_line.split(sep)
            # First column is the index (gene IDs) — skip it in dtype dict
            dtype_dict = {col: target_dtype for col in col_names[1:]}
        else:
            dtype_dict = None  # empty file — will fail below

        try:
            df = pd.read_csv(
                file_obj,
                sep=sep,
                index_col=0,
                comment="#",
                dtype=dtype_dict,
                compression=None,
            )
        except (ValueError, TypeError, OverflowError):
            # Fallback: data has floats, mixed types, or values
            # exceeding int32 range.  Re-read as float32 (50% less
            # than the default float64).
            file_obj.seek(start_pos)
            fallback_dict = (
                {col: np.float32 for col in col_names[1:]}
                if header_line else None
            )
            df = pd.read_csv(
                file_obj,
                sep=sep,
                index_col=0,
                comment="#",
                dtype=fallback_dict,
                compression=None,
            )
    except Exception as e:
        raise ValueError(
            f"Error leyendo la matriz de conteos '{real_name}': "
            f"{str(e)}. Verifica el formato del archivo."
        ) from e

    return df


def read_metadata_file(uploaded_file: UploadedFile) -> pd.DataFrame:
    """
    Reads a metadata file from a file uploaded by the user.
    Supports CSV or ZIP files containing a CSV.

    The file must have at least:
    - A "sample" column with the sample names, OR the sample
      names as the first column.
    - A "condition" column (or the configured name) with the
      experimental conditions.

    Expected format:
        sample      condition
        sample_1    control
        sample_2    control
        sample_3    treated
        ...

    Parameters
    ----------
    uploaded_file : UploadedFile
        Streamlit UploadedFile object.

    Returns
    -------
    pd.DataFrame
        Metadata DataFrame (without "sample" index yet; that is
        done by validation.validate_metadata_df).

    Raises
    ------
    ValueError
        If the file cannot be read.
    """
    file_obj, real_name = _prepare_file(uploaded_file)
    sep = detect_separator(real_name)

    try:
        df = pd.read_csv(
            file_obj,
            sep=sep,
            compression=None,
        )
    except Exception as e:
        raise ValueError(
            f"Error leyendo el metadata '{real_name}': "
            f"{str(e)}. Verifica el formato del archivo."
        ) from e

    return df


def results_to_csv(results_df: pd.DataFrame) -> str:
    """
    Converts the DESeq2 results DataFrame to a CSV string
    ready for download.

    Parameters
    ----------
    results_df : pd.DataFrame
        DataFrame with the DESeq2 results (log2FC, padj, etc.).

    Returns
    -------
    str
        String with the CSV content (includes index = gene names).

    Example
    -------
        csv_string = results_to_csv(results_df)
        st.download_button("Descargar", csv_string, "results.csv")
    """
    return results_df.to_csv(index=True)
