<p align="center">
  <img src="assets/beginseq_logo.png" alt="BeginSeq Studio" width="300">
</p>

<h1 align="center">BeginSeq Studio</h1>

<p align="center">
  <strong>Interactive RNA-seq analysis — no code required.</strong>
</p>

**BeginSeq Studio** is an interactive web application for RNA-seq analysis built with [Streamlit](https://streamlit.io/), created by **Elliot Ridout-Buhl**. It is designed for **beginners and researchers with no programming or command-line experience**, providing a fully code-free interface to perform bulk and single-cell differential expression analysis, build datasets from public repositories, and visualize results — all from the browser.

---

## What does the app do?

BeginSeq Studio is a multipage application with three integrated modules:

### 1. Bulk RNA-seq — Differential Expression Analysis

Upload a raw count matrix and a metadata file, and the app runs a complete **DESeq2** pipeline (Python implementation via [pydeseq2](https://github.com/owkin/PyDESeq2)):

- Median-of-ratios normalization
- Per-gene dispersion estimation
- Negative binomial GLM fitting
- Wald test with Benjamini-Hochberg correction
- Interactive volcano plot, MA plot, PCA, heatmap, and gene highlighting
- Downloadable CSV with full results (log2FoldChange, padj, baseMean, etc.)

### 2. scRNA-seq — Single-Cell RNA-seq Analysis

A complete single-cell analysis pipeline powered by [Scanpy](https://scanpy.readthedocs.io/), following the standard Sanbomics/scanpy workflow. The module includes three data-ingestion tools and a full analysis pipeline:

#### 10x File Integrator

Combine 10x Genomics output files (`matrix.mtx`, `features.tsv`/`genes.tsv`, `barcodes.tsv`) into a single `.h5ad` AnnData object:

- Accepts both plain-text and `.gz` compressed files
- Automatically fixes malformed `features.tsv` files:
  - **1-column** files (gene names only) — attempts to recover Ensembl gene IDs from an accompanying metadata file; falls back to duplicating the name
  - **2-column** files (gene_id + gene_name) — appends the missing `Gene Expression` feature-type column
- **Optional metadata integration** — upload a CSV/TSV with additional annotations:
  - Auto-detects whether the metadata belongs to **cells** (`adata.obs`) or **genes** (`adata.var`) by comparing index overlap
  - Manual override available (force obs or var)
  - Left-join merge preserves all cells/genes

#### SoupX Ambient RNA Removal

Remove ambient RNA contamination from raw single-cell data using the [SoupX](https://github.com/constantAmateur/SoupX) R package (called via [rpy2](https://rpy2.github.io/)):

- Requires **R** and the **SoupX** R package to be installed (see [Optional Dependencies](#optional-dependencies) below)
- Upload raw (unfiltered) and filtered `.h5ad` files
- Automatic contamination fraction estimation, or set a manual value
- Download the corrected `.h5ad` ready for downstream analysis
- *Note:* For production-grade ambient RNA removal, [CellBender](https://github.com/broadinstitute/CellBender) is recommended but requires a GPU. SoupX provides a lightweight CPU-only alternative.

#### H5AD Upload

Directly upload a pre-built `.h5ad` AnnData file to skip integration and proceed straight to analysis.

#### Analysis Pipeline

Once data is loaded (from any of the three sources above), the full scanpy pipeline runs interactively:

| Step | Description |
|------|-------------|
| QC Annotation | Flag mitochondrial, ribosomal, and hemoglobin genes |
| Cell & Gene Filtering | Filter by min counts, min genes, max genes, % mitochondrial |
| Doublet Detection | Identify and remove doublets via Scrublet |
| Normalization | Library-size normalization + log1p transform |
| HVG Selection | Select highly variable genes (flavors: seurat, seurat_v3, cell_ranger) |
| Scaling | Z-score scaling with optional max clipping |
| PCA | Principal component analysis (configurable number of components) |
| Batch Correction | Optional Harmony integration for multi-sample/multi-batch datasets |
| Neighbors | k-nearest-neighbor graph construction |
| UMAP | 2D UMAP embedding |
| Leiden Clustering | Community detection with adjustable resolution |
| Marker Genes | Rank genes per cluster (Wilcoxon, t-test, or logreg) |

All parameters are configurable through the sidebar. Results can be downloaded as `.h5ad` at the end.

#### Memory Considerations for the scRNA-seq Module

Single-cell datasets are inherently large. The scRNA-seq module processes everything in-memory, which means **RAM is the primary bottleneck**. Below are practical guidelines for users with limited RAM (8 GB machines).

##### Approximate RAM requirements

| Dataset size | Peak RAM (no batch correction) | Peak RAM (with Harmony) |
|:-------------|:-------------------------------|:------------------------|
| 5,000 cells  | ~2 GB                          | ~2.5 GB                |
| 10,000 cells | ~3 GB                          | ~4 GB                  |
| 30,000 cells | ~6 GB                          | ~7 GB                  |
| 50,000+ cells | ~10+ GB                       | ~12+ GB                |

> Harmony batch correction runs in an **isolated subprocess**: PyTorch and all Harmony buffers are loaded in a child process and fully freed when it exits. This prevents the torch runtime (~400 MB) from accumulating in the main application.

##### Strategies for low-RAM machines

**1. Aggressive early filtering (most effective)**

The biggest memory savings come from removing cells and genes early. Use stricter QC thresholds in the sidebar:

| Parameter | Default | Aggressive |
|:----------|:--------|:-----------|
| Min genes per cell | 200 | 500 |
| Max genes per cell | 5,000 | 3,000 |
| Min counts per cell | 500 | 1,000 |
| Max % mitochondrial | 20% | 10% |
| Min cells per gene | 3 | 10 |

Filtering a 50k-cell dataset down to 15k cells can reduce peak RAM by **60-70%**.

**2. Subsample for exploratory runs**

If you are tuning parameters (resolution, n_neighbors, etc.), consider subsampling your `.h5ad` before uploading:

```python
import scanpy as sc

adata = sc.read_h5ad("full_dataset.h5ad")

# Random subsample to 10k cells for fast exploration
sc.pp.subsample(adata, n_obs=10000)

adata.write_h5ad("subset_10k.h5ad")
```

Run the full dataset only after you have found satisfactory parameters on the subset.

**3. Reduce HVG count**

The default 2,000 highly variable genes works well for most datasets, but reducing to 1,000-1,500 can save memory during PCA and downstream steps with minimal impact on biological signal.

**4. Fewer principal components**

Reducing from 50 PCs to 30 PCs saves memory during PCA, Harmony, and the neighbor graph. Most biological variation is captured in the first 20-30 PCs.

**5. Disable optional steps**

- **Doublet detection** (Scrublet) requires temporary dense matrix operations. Disable it if your dataset was already filtered externally.
- **Batch correction** (Harmony) adds a subprocess with its own memory footprint. Only enable it if your data contains multiple batches.

**6. Close other applications**

Browsers, IDEs, and other applications compete for RAM. On an 8 GB machine, close unnecessary programs before running the pipeline.

##### What the pipeline already does to save memory

- **Single-mask filtering**: all cell QC filters are combined into one boolean mask and applied in a single `.copy()`, avoiding 4 intermediate copies.
- **Sparse matrix preservation**: the expression matrix stays in sparse (CSR) format throughout QC, filtering, and normalization. Densification only happens at the scaling step, and only for the HVG subset (~2k genes, not 20k+).
- **Early HVG subsetting**: genes are subset to HVGs *before* scaling and PCA, so the dense matrix is ~10x smaller than the full gene set.
- **Subprocess isolation for Harmony**: PyTorch + Harmony run in a `spawn` child process. When the subprocess exits, the OS reclaims all of its memory — nothing leaks into the main application.
- **Session state cleanup**: previous analysis results are freed from memory before each new pipeline run.
- **Aggressive garbage collection**: `gc.collect()` is called after every major pipeline step (filtering, doublets, normalization, scaling, PCA, batch correction, UMAP).

### 3. Dataset Creator — GDC/TCGA Downloader

Build analysis-ready datasets directly from the **NCI Genomic Data Commons (GDC)**:

- Browse and select TCGA projects
- Download STAR-Counts RNA-seq files via the GDC REST API
- Auto-assemble a DESeq2-ready count matrix + metadata
- No R, no TCGAbiolinks, no external wrappers needed

### Additional features

- Bilingual interface (English / Spanish)
- Configurable significance thresholds and plot aesthetics via `config.py`
- Input validation with clear, actionable error messages
- Local file-upload limit of **5 GB** (automatically configured when running on localhost)

---

## Installation

### Prerequisites

- Python >= 3.10

### Steps

```bash
# 1. Clone or download the project
cd "BeginSeq Studio"

# 2. Create a virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # Linux/Mac
# venv\Scripts\activate   # Windows

# 3. Install core dependencies (Bulk RNA-seq + Dataset Creator)
pip install -r requirements.txt

# 4. (Optional) Install scRNA-seq dependencies — local use only
pip install -r requirements-scrna.txt

# 5. Run the application
streamlit run app.py
```

The application will open at `http://localhost:8501`.

> **Tip:** By default the app runs in headless mode (no browser pop-up). To have the browser open automatically on launch, edit `.streamlit/config.toml` and set `headless = false`.

> **Note on Streamlit Cloud:** Only `requirements.txt` is needed for cloud deployments. The scRNA-seq module is automatically disabled on cloud (it requires more RAM than cloud instances provide), so `requirements-scrna.txt` is not needed there.

### Optional Dependencies

#### scRNA-seq module

The scRNA-seq module requires additional packages listed in `requirements-scrna.txt` (scanpy, anndata, harmonypy, etc.). These are **not included in the base `requirements.txt`** because they pull in heavy dependencies (PyTorch, LLVM) that are unnecessary for cloud deployments where the module is disabled.

#### SoupX Ambient RNA Removal

The **SoupX** feature in the scRNA-seq module additionally requires:

1. **R** (>= 4.0) installed and available on `PATH`
2. The **SoupX** R package — install it from R:
   ```r
   install.packages("SoupX")
   ```
3. The **rpy2** Python package (included in `requirements-scrna.txt`)

If R or SoupX is not available, the rest of BeginSeq Studio works normally — only the SoupX section will show a status message explaining what is missing.

---

## Project Structure

```
BeginSeq Studio/
├── app.py                      # Main Streamlit entry point
├── config.py                   # Central configuration (thresholds, plot styles)
├── i18n.py                     # Internationalisation (EN / ES translations)
├── runtime_utils.py            # Localhost detection & upload-limit helpers
├── analysis.py                 # DESeq2 pipeline orchestrator (facade)
├── deseq_runner.py             # PyDESeq2 wrapper
├── validation.py               # Input data validation & normalization
├── data_io.py                  # File reading/writing (CSV, TSV, ZIP)
├── visualization.py            # Bulk RNA-seq plot generation
├── classification.py           # Expression-based sample classification
├── gdc_client.py               # GDC REST API client (Dataset Creator)
├── scrna_pipeline.py           # scRNA-seq backend (scanpy, SoupX, 10x integrator)
├── scrna_visualization.py      # scRNA-seq plot generation
├── auto_shutdown.py            # Auto-shutdown on browser disconnect
├── assets/
│   ├── beginseq_logo.png       # Project logo
│   ├── helix.MP4               # DNA helix animation (landing page)
│   ├── Cell.MP4                # Cell animation (scRNA-seq page)
│   └── Box.MP4                 # Box animation (Dataset Creator page)
├── requirements.txt            # Core dependencies (Bulk + Dataset Creator)
├── requirements-scrna.txt      # scRNA-seq dependencies (local-only)
├── .streamlit/
│   └── config.toml             # Streamlit server settings (upload limits)
├── pages/
│   ├── 1_🧬_Bulk_RNA-seq.py    # Bulk RNA-seq DE analysis page
│   ├── 2_🔬_scRNA-seq.py       # Single-cell analysis page
│   └── 3_📦_Dataset_Creator.py # GDC/TCGA dataset downloader
├── tutorial_bulk_rnaseq.md     # Bulk RNA-seq tutorial
├── tutorial_scrna_seq.md       # scRNA-seq tutorial
├── tutorial_dataset_creator.md # Dataset Creator tutorial
├── LICENSE                     # MIT License
├── THIRD_PARTY_NOTICES.txt     # Third-party dependency attributions
└── README.md                   # This file
```

---

## License

This project is released under the **MIT License**.

Copyright © 2026 Elliot Ridout-Buhl.

You are free to use, modify, and distribute this software for academic, educational, or commercial purposes, provided that the original copyright notice and this permission notice are included in all copies or substantial portions of the software.

See the [LICENSE](LICENSE) file for the full text.

---

## Disclaimer

### Intended use

This software is intended for **research and educational purposes only**. It is not designed for clinical, diagnostic, therapeutic, or medical decision-making use.

BeginSeq Studio is primarily designed for **beginners and users without computational or programming experience**. While the interface simplifies complex bioinformatics workflows, users should be aware that simplified tools do not substitute for a thorough understanding of the underlying statistical methods and biological context.

### No warranty

This software is provided **"as is"**, without warranty of any kind, express or implied, including but not limited to merchantability, fitness for a particular purpose, or non-infringement.

Results generated by this software depend on input data quality, preprocessing choices, and parameter selection. No warranty is made regarding accuracy, completeness, or fitness for a particular scientific purpose.

### User responsibility

- Interpretation of datasets and biological conclusions derived from this software remain the **sole responsibility of the user**.
- Computational results should be **independently validated experimentally** before drawing biological or clinical conclusions.
- Visualizations are exploratory tools and should not be considered definitive evidence of biological effects.
- This software **does not replace** expert bioinformatics or statistical consultation.
- Automated analysis pipelines simplify workflows but do not eliminate the need for methodological understanding.

### Data and compliance

- Users are responsible for complying with applicable data usage agreements, privacy regulations, and ethical guidelines when analyzing human or clinical datasets.
- Users are responsible for securing their data and ensuring compliance with institutional data protection policies.

### Reproducibility

Results may vary depending on software versions, computational environment, and analysis parameters. Reproducibility is not guaranteed without controlled environments.

### Third-party dependencies

This software relies on third-party libraries. Their respective licenses and limitations apply independently. See [THIRD_PARTY_NOTICES.txt](THIRD_PARTY_NOTICES.txt) for details.

### Limitation of liability

The authors shall not be liable for any damages arising from the use or inability to use this software, including data loss, incorrect analysis results, or research delays.

---

## Dataset Policy

### User-uploaded data

- All data uploaded to BeginSeq Studio is processed **locally** on the machine where the application is running.
- No user data is transmitted to external servers (except when using the Dataset Creator, which queries the GDC API).
- Uploaded files are held in memory during the session and are not persisted to disk unless the user explicitly downloads results.

### GDC / TCGA data

The Dataset Creator module retrieves open-access data from the [NCI Genomic Data Commons (GDC)](https://gdc.cancer.gov/) via its public REST API.

- Only **open-access** data is retrieved. BeginSeq Studio does not access controlled-access data and does not handle dbGaP authorization tokens.
- Users must comply with the [GDC Data Use Agreement](https://gdc.cancer.gov/about-gdc/gdc-policies) and the [NIH Genomic Data Sharing Policy](https://sharing.nih.gov/genomic-data-sharing-policy) when using downloaded datasets.
- TCGA data is subject to the policies described by the [TCGA publication guidelines](https://www.cancer.gov/ccg/research/genome-sequencing/tcga/using-tcga-data). If you use TCGA data in a publication, cite the original TCGA project and the GDC.
- Downloaded data may contain information derived from human samples. Users are responsible for ensuring that their use of this data complies with all applicable ethical guidelines and institutional policies.

---

## Suggested Citation

If you use BeginSeq Studio in your research, please cite:

```
Ridout-Buhl, E. (2025). BeginSeq Studio: Interactive RNA-seq differential expression analysis tool.
https://github.com/JACKNINES/beginseq-studio
```

If your analysis uses the DESeq2 statistical method, please also cite:

> Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* **15**, 550 (2014). https://doi.org/10.1186/s13059-014-0550-8

And the Python implementation:

> Muzellec, B., Telenczuk, M., Cabeli, V. & Andreux, M. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. *Bioinformatics* **39**, btad547 (2023). https://doi.org/10.1093/bioinformatics/btad547

If your analysis uses the scRNA-seq module, please also cite Scanpy:

> Wolf, F.A., Angerer, P. & Theis, F.J. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology* **19**, 15 (2018). https://doi.org/10.1186/s13059-017-1382-0

If you use the Harmony batch correction feature, cite:

> Korsunsky, I. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods* **16**, 1289–1296 (2019). https://doi.org/10.1038/s41592-019-0619-0

If you use the SoupX ambient RNA removal feature, cite:

> Young, M.D. & Behjati, S. SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. *GigaScience* **9**, giaa151 (2020). https://doi.org/10.1093/gigascience/giaa151

If you use TCGA data obtained through the Dataset Creator, cite:

> The Cancer Genome Atlas Research Network. Comprehensive genomic characterization defines human glioblastoma genes and core pathways. *Nature* **455**, 1061–1068 (2008).

> Grossman, R.L. et al. Toward a Shared Vision for Cancer Genomic Data. *N Engl J Med* **375**, 1109–1112 (2016). https://doi.org/10.1056/NEJMp1607591
