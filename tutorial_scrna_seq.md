# Tutorial: scRNA-seq — Single-Cell RNA-seq Analysis

## Table of Contents

1. [Overview](#overview)
2. [RAM Considerations — Read This First](#ram-considerations--read-this-first)
3. [What You Need Before Starting](#what-you-need-before-starting)
4. [Data Ingestion — Three Methods](#data-ingestion--three-methods)
   - [Method A: 10x File Integrator](#method-a-10x-file-integrator)
   - [Method B: SoupX Ambient RNA Removal](#method-b-soupx-ambient-rna-removal)
   - [Method C: Direct H5AD Upload](#method-c-direct-h5ad-upload)
5. [Understanding the Sidebar Parameters](#understanding-the-sidebar-parameters)
6. [The Pipeline — Step by Step](#the-pipeline--step-by-step)
   - [Step 1: QC Annotation](#step-1-qc-annotation)
   - [Step 2: Cell & Gene Filtering](#step-2-cell--gene-filtering)
   - [Step 3: Doublet Detection](#step-3-doublet-detection)
   - [Step 4: Doublet Removal](#step-4-doublet-removal)
   - [Step 5: Normalization](#step-5-normalization)
   - [Step 6: Highly Variable Gene (HVG) Selection](#step-6-highly-variable-gene-hvg-selection)
   - [Step 7: Scaling](#step-7-scaling)
   - [Step 8: PCA](#step-8-pca)
   - [Step 9: Batch Effect Correction (Harmony)](#step-9-batch-effect-correction-harmony)
   - [Step 10: Neighborhood Graph](#step-10-neighborhood-graph)
   - [Step 11: UMAP Embedding](#step-11-umap-embedding)
   - [Step 12: Leiden Clustering](#step-12-leiden-clustering)
   - [Step 13: Marker Gene Identification](#step-13-marker-gene-identification)
7. [Interpreting Results](#interpreting-results)
8. [Visualization Tabs](#visualization-tabs)
9. [Downloading Results](#downloading-results)
10. [RAM Optimization Guide — In Depth](#ram-optimization-guide--in-depth)
11. [Frequently Asked Questions](#frequently-asked-questions)
12. [Glossary](#glossary)

---

## Overview

The scRNA-seq module provides a **complete single-cell RNA-seq analysis pipeline** powered by [Scanpy](https://scanpy.readthedocs.io/). It takes raw or pre-processed single-cell data and produces:

- Quality control metrics and visualizations
- Filtered, normalized, and scaled expression data
- Dimensionality reduction (PCA + UMAP)
- Optional batch effect correction (Harmony)
- Cell clustering (Leiden algorithm)
- Marker gene identification per cluster
- Interactive plots and downloadable results

The entire pipeline follows the standard Scanpy workflow and is designed to be run **without any programming experience**.

---

## RAM Considerations — Read This First

> **This is the most important section of this tutorial.** Single-cell datasets are large. Unlike bulk RNA-seq (tens of samples), scRNA-seq can involve **thousands to hundreds of thousands of cells**, each measured across 20,000+ genes. This means **RAM (memory) is the primary bottleneck**.

### How much RAM do I need?

| Dataset Size | Peak RAM (no batch correction) | Peak RAM (with Harmony) |
|:-------------|:-------------------------------|:------------------------|
| 5,000 cells  | ~2 GB                          | ~2.5 GB                 |
| 10,000 cells | ~3 GB                          | ~4 GB                   |
| 30,000 cells | ~6 GB                          | ~7 GB                   |
| 50,000+ cells | ~10+ GB                       | ~12+ GB                 |

### What happens if I run out of RAM?

The pipeline **will crash** without warning. On macOS, the process gets killed by the system. On Linux, the OOM (Out of Memory) killer terminates it. You will see the Streamlit progress bar freeze and the pipeline will report an error.

### Key takeaway

If you have an **8 GB machine**:
- Datasets up to **~30,000 cells** should work with default settings.
- Datasets of **50,000+ cells** will require aggressive filtering or subsampling.
- **Always close other applications** (browsers, IDEs, Slack) before running the pipeline.

If you have a **16 GB machine**, you can comfortably handle most single-cell datasets up to ~100,000 cells.

**See [RAM Optimization Guide](#ram-optimization-guide--in-depth) at the end of this tutorial for detailed strategies.**

---

## What You Need Before Starting

### Input Data

You need single-cell expression data in one of these formats:

| Format | Description | How to get it |
|--------|-------------|---------------|
| **H5AD** (.h5ad) | AnnData file — the standard format for Scanpy | Output from Cell Ranger (`filtered_feature_bc_matrix.h5ad`), Scanpy, or any AnnData-compatible tool |
| **10x Genomics files** | Three files: `matrix.mtx`, `features.tsv`, `barcodes.tsv` | Output from Cell Ranger (`outs/filtered_feature_bc_matrix/`) |

### Optional Data

- **Metadata CSV/TSV** — additional cell or gene annotations (sample IDs, donor IDs, batch labels).
- **Raw (unfiltered) H5AD** — for SoupX ambient RNA removal (from Cell Ranger `raw_feature_bc_matrix`).

---

## Data Ingestion — Three Methods

The scRNA-seq module provides three ways to load data. Each is accessible from the top of the page.

---

### Method A: 10x File Integrator

**When to use:** You have the three output files from Cell Ranger (`matrix.mtx`, `features.tsv`, `barcodes.tsv`) and want to combine them into a single `.h5ad` file.

#### Files required

| File | Description | Typical name |
|------|-------------|--------------|
| **Matrix** | Sparse count matrix | `matrix.mtx` or `matrix.mtx.gz` |
| **Features** | Gene identifiers | `features.tsv`, `genes.tsv`, or `.gz` variants |
| **Barcodes** | Cell barcodes | `barcodes.tsv` or `barcodes.tsv.gz` |

#### Smart fixes for malformed files

The integrator automatically fixes common problems with `features.tsv` files:

- **1-column files** (gene names only): The integrator duplicates the column and adds `Gene Expression` as the feature type. If you provide a metadata file with Ensembl gene IDs, it attempts to recover them.
- **2-column files** (gene_id + gene_name): The integrator adds the missing `Gene Expression` column.
- **3-column files**: Used as-is.

#### Optional metadata integration

You can upload an additional CSV/TSV metadata file. The integrator:

1. **Auto-detects** whether it contains cell metadata (`obs`) or gene metadata (`var`) by checking which index overlaps more.
2. **Merges** it via left-join (all original cells/genes are preserved).
3. You can **override** the auto-detection if needed (force obs or force var).

#### How to use

1. Upload the three core files in the respective uploaders.
2. (Optional) Upload a metadata file and choose the target (auto, obs, or var).
3. Click **"Build H5AD"**.
4. Choose a filename and click **"Download"**.

The downloaded `.h5ad` file can then be uploaded in the H5AD uploader below to run the analysis pipeline.

#### RAM impact

The 10x File Integrator is lightweight. It creates the AnnData object in memory and writes it to a file. Even for datasets with 100,000+ cells, this step typically uses less than 1 GB of RAM.

---

### Method B: SoupX Ambient RNA Removal

**When to use:** You want to remove ambient RNA contamination from your single-cell data before analysis. Ambient RNA is mRNA released from lysed cells during droplet-based capture — it contaminates every cell's profile with a "background soup" of irrelevant transcripts.

#### Requirements

This feature requires:
- **R** (>= 4.0) installed and available on PATH
- The **SoupX** R package installed: `install.packages("SoupX")` from an R console
- **rpy2** Python package (included in `requirements.txt`)

If any of these is missing, the app will show a status message explaining what's needed.

#### Files required

| File | Description |
|------|-------------|
| **Raw (unfiltered) H5AD** | Contains ALL droplets (cells + empty droplets). From Cell Ranger: `raw_feature_bc_matrix/` converted to H5AD |
| **Filtered H5AD** | Contains only cell-containing droplets. From Cell Ranger: `filtered_feature_bc_matrix/` converted to H5AD |

#### How to use

1. Upload both the raw and filtered `.h5ad` files.
2. Choose contamination estimation:
   - **Automatic** (recommended): SoupX estimates the contamination fraction from the data.
   - **Manual**: Set a fixed contamination fraction (0.01 to 0.50). Typical values: 0.05–0.15.
3. Click **"Run SoupX"**.
4. Download the corrected `.h5ad` file.

#### Why this matters

Ambient RNA can:
- Inflate expression of highly-expressed genes (e.g., hemoglobin in blood samples).
- Create false marker genes that are actually ambient contamination.
- Distort clustering by making cell types look more similar than they are.

#### When to skip this

- If your data was already processed with **CellBender** (GPU-based, more thorough).
- If your data has very low contamination (< 2%).
- If you don't have access to the raw (unfiltered) matrix.

#### RAM impact

SoupX runs through R via rpy2. It loads both matrices into memory simultaneously, so the peak RAM is roughly **2x the size of your raw matrix**. For a raw matrix with 500,000 droplets × 30,000 genes, expect ~4-6 GB peak RAM.

---

### Method C: Direct H5AD Upload

**When to use:** You already have a `.h5ad` file ready for analysis (from Cell Ranger, the 10x Integrator, SoupX, or any other preprocessing).

Simply upload the file and the pipeline will start from Step 1 (QC Annotation).

#### File size limits

When running locally, the maximum upload size is **5 GB**.

#### What the app shows after loading

- **Cell count** and **gene count**
- Preview of **cell metadata** (`obs`) — sample IDs, pre-existing annotations
- Preview of **gene metadata** (`var`) — gene names, pre-existing flags

> **Tip:** If your H5AD already has QC annotations (e.g., from a previous Scanpy run), the pipeline will **re-annotate** from scratch. This ensures consistency with the pipeline's own mitochondrial/ribosomal/hemoglobin gene detection patterns.

#### RAM impact

Loading the H5AD file into memory requires approximately **1-2x the file size**. A 500 MB `.h5ad` file will use ~500 MB to 1 GB of RAM just for loading.

---

## Understanding the Sidebar Parameters

All pipeline parameters are configured in the **sidebar** (left panel). Here's a complete reference:

### Quality Control Parameters

| Parameter | Default | Range | Purpose |
|-----------|---------|-------|---------|
| **Min genes per cell** | 200 | 0–10,000 | Remove cells expressing fewer than this many genes. These are likely empty droplets or debris. |
| **Max genes per cell** | 5,000 | 500–50,000 | Remove cells expressing more than this many genes. These are likely doublets (two cells in one droplet). |
| **Min counts per cell** | 500 | 0–100,000 | Remove cells with fewer total UMI counts. Low counts = poor capture quality. |
| **Max counts per cell** | 50,000 | 1,000–500,000 | Remove cells with too many UMI counts. Very high counts = likely doublets. |
| **Max % mitochondrial** | 20% | 0–100% | Remove cells where more than this % of reads map to mitochondrial genes. High mt% = dying/stressed cells. |
| **Min cells per gene** | 3 | 1–100 | Remove genes detected in fewer than this many cells. Rarely-detected genes add noise. |

### Analysis Parameters

| Parameter | Default | Range | Purpose |
|-----------|---------|-------|---------|
| **Number of HVGs** | 2,000 | 500–10,000 | How many highly variable genes to select. More HVGs capture more biology but use more RAM. |
| **Number of PCs** | 40 | 5–50 | How many principal components to compute and use. More PCs capture more variance but add noise. |
| **Number of neighbors** | 15 | 5–100 | k for the k-nearest-neighbor graph. Higher = smoother/larger clusters. Lower = finer detail. |
| **UMAP min_dist** | 0.5 | 0.0–1.0 | How tightly UMAP packs points. Lower = tighter clusters. Higher = more spread out. |
| **Leiden resolution** | 0.5 | 0.05–3.0 | Controls number of clusters. Higher = more clusters. Lower = fewer, larger clusters. |
| **DE method** | wilcoxon | wilcoxon, t-test, t-test_overestim_var, logreg | Statistical test for marker gene identification. |
| **Marker genes per cluster** | 25 | 5–100 | How many top marker genes to compute per cluster. |
| **Enable doublet detection** | Checked | On/Off | Whether to run Scrublet for doublet detection. Disable if your data was already filtered. |

### Batch Effect Correction

| Parameter | Default | Purpose |
|-----------|---------|---------|
| **Enable batch correction** | Unchecked | Turn on Harmony batch correction. Only enable if your data has multiple batches. |
| **Batch variable column** | (auto) | Which column in `obs` identifies the batch (e.g., `sample`, `donor`, `experiment`). |

### Annotation Columns

| Parameter | Default | Purpose |
|-----------|---------|---------|
| **Cell type column** | (none) | If your H5AD has pre-annotated cell types, select that column to color plots by cell type instead of cluster. |
| **Marker gene groupby** | leiden | Which grouping variable to use for marker gene analysis. Usually `leiden`, but can be a cell type column. |

### Visualization Settings

| Parameter | Default | Purpose |
|-----------|---------|---------|
| **Legend position** | best | Where to place the legend in plots. |

---

## The Pipeline — Step by Step

When you click **"Run Pipeline"**, the app executes 13 steps in sequence. Here's what each step does and why it matters.

---

### Step 1: QC Annotation

**What it does:** Flags special gene categories in your data:
- **Mitochondrial genes** (`MT-` prefix) — markers of cell stress/death.
- **Ribosomal genes** (`RPS`, `RPL` prefixes) — ubiquitously expressed, can dominate signal.
- **Hemoglobin genes** (`HB` prefix, e.g., `HBA1`, `HBB`) — markers of red blood cell contamination.

**Why it matters:** These annotations are used to calculate quality metrics:
- `pct_counts_mt` — percentage of a cell's counts mapping to mitochondrial genes.
- `pct_counts_ribo` — percentage from ribosomal genes.
- `pct_counts_hb` — percentage from hemoglobin genes.

These metrics drive the filtering in Step 2.

**RAM impact:** Minimal. This step adds a few columns to the metadata without copying the expression matrix.

---

### Step 2: Cell & Gene Filtering

**What it does:** Removes low-quality cells and rarely-detected genes based on the QC thresholds you set in the sidebar.

**Cells removed if:**
- Too few genes detected (< `min_genes`)
- Too many genes detected (> `max_genes`) — likely doublets
- Too few total counts (< `min_counts`)
- Too many total counts (> `max_counts`) — likely doublets
- Too much mitochondrial expression (> `max_pct_mt`) — dying cells

**Genes removed if:**
- Detected in fewer than `min_cells` cells

**Why it matters:** Filtering is the **most critical step for data quality AND memory usage**. Removing poor-quality cells:
- Eliminates noise from empty droplets, debris, and dying cells.
- Reduces the matrix size (fewer cells × fewer genes = less RAM).
- Improves downstream analysis by focusing on real, healthy cells.

**RAM impact: HIGH.** This is where the biggest memory savings happen.

> **For 8 GB machines:** This is your most powerful lever. Using aggressive filtering (see [RAM Optimization Guide](#ram-optimization-guide--in-depth)) can reduce a 50,000-cell dataset to 15,000 cells, cutting peak RAM by **60-70%**.

**Pipeline optimization:** The app applies all cell filters as a **single boolean mask** in one operation, avoiding 4 separate `.copy()` calls that would each duplicate the entire expression matrix.

---

### Step 3: Doublet Detection

**What it does:** Uses **Scrublet** to identify potential doublets — droplets that captured two cells instead of one.

**How it works:**
1. Scrublet simulates doublets by combining random pairs of real cells.
2. It scores each real cell on how similar it is to the simulated doublets.
3. Cells scoring above the threshold are flagged as `predicted_doublet = True`.

**Why it matters:** Doublets appear as cells with unusually high gene counts that fall between two clusters. They create:
- False intermediate cell states
- Incorrect cluster assignments
- Spurious marker genes

**When to disable:** If your data was already filtered for doublets externally (e.g., with CellRanger's built-in doublet filter, or Demuxlet for multiplexed samples).

**RAM impact: MODERATE.** Scrublet temporarily creates a dense matrix for the simulation step. On a 30,000-cell dataset, this can add ~1-2 GB of peak RAM.

---

### Step 4: Doublet Removal

**What it does:** Removes cells flagged as doublets by Scrublet (if enabled).

**Why it's a separate step:** The detection and removal are separated so you can:
1. See how many doublets were detected in the statistics.
2. Optionally keep them (by disabling the "remove doublets" checkbox).

**RAM impact:** Minimal. Only removes rows from the matrix.

---

### Step 5: Normalization

**What it does:** Two-step normalization:
1. **Library-size normalization** — scales each cell so that its total counts equal the median total counts across all cells. This corrects for differences in sequencing depth between cells.
2. **Log1p transformation** — applies `log(x + 1)` to all values. This reduces the skewness of count data and makes the expression distribution more symmetric.

**Why it matters:** Raw counts are not directly comparable between cells because:
- Some cells were sequenced more deeply than others.
- Count data is heavily right-skewed (a few genes dominate).

After normalization, expression values are on a comparable scale across cells.

**RAM impact:** Moderate. The normalization is performed on the sparse matrix in-place, but a `.copy()` is needed to store the raw counts in `adata.raw` for later use by the marker gene analysis.

---

### Step 6: Highly Variable Gene (HVG) Selection

**What it does:** Identifies the genes with the highest cell-to-cell variability, which are most likely to drive biological differences between cell types.

**Available flavors:**
- **seurat_v3** (default): Uses variance-stabilizing transformation. Best for raw count data.
- **seurat**: Original Seurat method. Works on log-normalized data.
- **cell_ranger**: Cell Ranger's HVG selection method.

**Why it matters:** A typical scRNA-seq dataset measures ~20,000-30,000 genes, but most are:
- Not expressed at all (near-zero in most cells).
- Expressed at constant levels (housekeeping genes — informative but not discriminative).

Only ~2,000-5,000 genes vary enough to distinguish cell types. By selecting HVGs:
- **PCA and clustering focus on signal, not noise.**
- **RAM usage drops dramatically** (downstream steps work with ~2K genes instead of 20K+).
- **Computation is much faster.**

**RAM impact:** The HVG selection itself is lightweight. The major savings come in the next step (scaling), because only HVGs are densified.

> **For 8 GB machines:** Reducing HVGs from 2,000 to 1,000-1,500 saves memory during PCA and downstream steps with minimal impact on biological signal.

---

### Step 7: Scaling

**What it does:** Z-score standardizes each gene: mean = 0, variance = 1. Optionally clips extreme values (default: max = 10) to prevent a single outlier cell from dominating.

**Critical implementation detail:** Before scaling, the pipeline **subsets the data to only HVGs**. This means the dense (non-sparse) scaled matrix is ~10x smaller than the full gene set.

**Why it matters:**
- PCA requires features on the same scale. Without scaling, highly-expressed genes would dominate the first principal components.
- Clipping extreme values prevents outliers from distorting the embedding.

**What happens to non-HVG genes?** They are preserved in `adata.raw`, which stores the pre-scaling, pre-subsetting normalized data. Marker gene analysis (Step 13) uses `adata.raw` to rank all genes, not just HVGs.

**RAM impact: MODERATE TO HIGH.** This is the step where the sparse matrix becomes dense (for the HVG subset). With 2,000 HVGs and 30,000 cells, the dense matrix is ~230 MB (float32). With 20,000 genes (no HVG subsetting), it would be ~2.3 GB.

---

### Step 8: PCA

**What it does:** Reduces the ~2,000 HVG dimensions to 50 principal components (or fewer, depending on your setting).

**Why it matters:**
- 2,000 dimensions is too many for neighbor graph construction and clustering.
- PCA captures the main axes of variation in a compact representation.
- Most biological signal is captured in the first 20-30 PCs.

**How to choose the number of PCs:**
- After running the pipeline, check the **Elbow Plot** in the PCA tab.
- Look for the "elbow" — the point where additional PCs add diminishing variance.
- Typical values: 20-40 PCs.

**RAM impact:** Moderate. PCA computes on the dense scaled matrix. The output is a small matrix (n_cells × n_PCs).

> **For 8 GB machines:** Reducing from 50 PCs to 30 PCs saves memory during PCA, Harmony, and the neighbor graph. Most biological variation is captured in the first 20-30 PCs anyway.

---

### Step 9: Batch Effect Correction (Harmony)

**What it does:** If enabled, corrects batch effects in the PCA embedding space using [Harmony](https://github.com/immunogenomics/harmony).

**When to use:**
- Your data comes from **multiple samples**, **multiple donors**, or **multiple experiments**.
- Cells cluster primarily by batch (sample/donor) instead of by biology (cell type).
- You can see batch-driven separation in the PCA plot.

**When NOT to use:**
- Your data is from a single sample/donor/experiment.
- Batch and biology are completely confounded (all cell types come from one batch).
- The data is already batch-corrected (e.g., processed with Scanorama or BBKNN upstream).

**How it works:**
1. Takes the PCA embedding matrix (n_cells × n_PCs).
2. Iteratively adjusts the embedding so that cells from different batches overlap while preserving biological structure.
3. Stores the corrected embedding in `adata.obsm["X_pca_harmony"]`.
4. The corrected embedding is used for all downstream steps (neighbors, UMAP, clustering).

**How to enable:**
1. In the sidebar, check **"Enable batch correction (Harmony)"**.
2. Select the column that identifies the batch (e.g., `sample`, `donor`, `experiment`).
3. The app shows a preview: `Found X batches in column Y`.

#### RAM Impact: IMPORTANT

Harmony uses **PyTorch** internally (harmonypy >= 0.2), which adds ~400 MB just for the torch runtime. To prevent this from exploding RAM:

> **BeginSeq Studio runs Harmony in an isolated subprocess.** The entire Harmony computation (including PyTorch) happens in a separate child process. When it finishes, the operating system reclaims ALL of that memory. Only the corrected PCA matrix (a small float32 array) is sent back to the main process.

This means:
- **Peak RAM increases temporarily** during Harmony (~1-2 GB extra for PyTorch + Harmony buffers).
- **After Harmony finishes**, the RAM returns to pre-Harmony levels.
- The main application process never accumulates PyTorch's memory footprint.

> **For 8 GB machines:** Only enable Harmony if your data genuinely has multiple batches. The temporary RAM spike can push you close to the limit with 30,000+ cells.

---

### Step 10: Neighborhood Graph

**What it does:** Builds a k-nearest-neighbor (kNN) graph where each cell is connected to its closest neighbors in PCA space (or Harmony-corrected PCA space).

**Why it matters:**
- The neighbor graph defines the "similarity structure" of the data.
- UMAP and Leiden clustering both operate on this graph.
- The number of neighbors (k) controls the balance between local and global structure.

**How to interpret the `n_neighbors` parameter:**
- **Lower (5-10):** Captures fine-grained local structure. Results in smaller, more specific clusters.
- **Default (15):** Good balance for most datasets.
- **Higher (30-100):** Captures broader global structure. Results in larger, more general clusters.

**RAM impact:** The neighbor graph is sparse and memory-efficient. This step is rarely a bottleneck.

---

### Step 11: UMAP Embedding

**What it does:** Computes a 2D embedding (UMAP) for visualization. UMAP preserves local neighborhood relationships while arranging cells in a 2D space.

**How to interpret the `min_dist` parameter:**
- **Lower (0.0-0.2):** Points are packed tightly. Clusters appear as dense blobs.
- **Default (0.5):** Moderate spacing. Good for most datasets.
- **Higher (0.8-1.0):** Points are more spread out. Useful when clusters overlap.

**Important caveats:**
- UMAP is a **visualization tool**, not a clustering tool. The positions are for visual interpretation only.
- Distances between clusters on the UMAP are **not meaningful** — two clusters that look far apart may actually be similar.
- The UMAP embedding can change with different random seeds or parameters.

**RAM impact:** Moderate. UMAP stores a small matrix (n_cells × 2).

---

### Step 12: Leiden Clustering

**What it does:** Identifies communities (clusters) in the neighbor graph using the Leiden algorithm — a graph-based community detection method.

**How to interpret the `resolution` parameter:**
- **Lower (0.1-0.3):** Fewer, larger clusters. Broad cell type groupings.
- **Default (0.5):** Good starting point for most datasets.
- **Higher (1.0-3.0):** More, smaller clusters. Finer sub-populations.

**Practical advice:**
- Start with the default (0.5) and adjust based on the UMAP plot.
- If you see clusters that obviously should be split, increase resolution.
- If you see too many clusters that look biologically similar, decrease resolution.
- There is no "correct" resolution — it depends on the biological question.

**RAM impact:** Minimal. Clustering operates on the sparse neighbor graph.

---

### Step 13: Marker Gene Identification

**What it does:** For each cluster, identifies the genes that are most specifically expressed in that cluster compared to all other clusters.

**Available statistical tests:**
- **Wilcoxon rank-sum** (default): Non-parametric, robust, recommended for most cases.
- **t-test**: Parametric, faster, assumes roughly normal expression.
- **t-test_overestim_var**: t-test with overestimated variance — more conservative.
- **logreg**: Logistic regression — identifies genes that best classify each cluster.

**Why it matters:** Marker genes are how you assign **biological identity** to clusters:
- Cluster 0 has high CD3D, CD3E, CD2 → T cells
- Cluster 1 has high CD19, MS4A1, CD79A → B cells
- Cluster 2 has high CD14, LYZ, CST3 → Monocytes

**Output:** For each cluster, the top N genes (default: 25) with their:
- **names** — gene symbol
- **scores** — test statistic (higher = more specific)
- **logfoldchanges** — log2 fold-change vs. other clusters
- **pvals** — raw p-value
- **pvals_adj** — adjusted p-value (Benjamini-Hochberg)

**RAM impact:** Moderate. Uses `adata.raw` (the pre-scaled, pre-subsetted normalized data) for the statistical test.

---

## Interpreting Results

### Summary Metrics

After the pipeline completes, the top of the results section shows four key metrics:

| Metric | What it tells you |
|--------|------------------|
| **Cells** | How many cells survived all filtering steps |
| **Genes** | How many genes are in the HVG-subsetted data |
| **Clusters** | How many Leiden clusters were found |
| **HVGs** | How many highly variable genes were selected |

### Pipeline Statistics

Expand the **"Pipeline Statistics"** section to see detailed information:

- **Filtering stats:** How many cells/genes were removed and why.
- **Doublet stats:** How many doublets were detected and removed.
- **HVG stats:** Number of highly variable genes selected.
- **Harmony stats** (if enabled): Number of batches corrected and PCs used.
- **Leiden stats:** Number of clusters and resolution used.

---

## Visualization Tabs

### QC Tab

Shows quality control distributions **before filtering**:

- **Violin plots:** Distribution of `n_genes`, `total_counts`, and `pct_counts_mt` across all cells.
- **Scatter plots:** Relationship between total counts and number of genes, colored by mitochondrial percentage.

**What to look for:**
- A clear separation between high-quality cells and low-quality cells.
- Outlier cells with extremely high counts or gene numbers (potential doublets).
- Cells with high mitochondrial percentage (dying cells — should be filtered).

### HVG Tab

Shows the highly variable gene selection:

- **Scatter plot:** Mean expression vs. variance for all genes.
- **Highlighted points:** HVGs are colored, non-HVGs are grey.

**What to look for:**
- HVGs should span a range of expression levels, not concentrate at one end.
- Known marker genes (CD3D, MS4A1, etc.) should be among the HVGs.

### PCA Tab

Two plots:

- **Elbow plot:** Variance explained by each PC. Look for the "elbow" point.
- **PCA embedding:** First two PCs, colored by cluster (or cell type if annotated).

**What to look for:**
- The elbow plot helps you decide how many PCs to use. If most variance is captured in the first 15 PCs, you could safely reduce `n_pcs` to 20-25.
- In the PCA embedding, well-separated cell populations should appear as distinct groups.

### UMAP Tab

The main visualization for exploring your data:

- **UMAP colored by cluster (or cell type):** See the overall structure.
- **Color selector:** Change the coloring to any metadata column (e.g., sample, donor, condition).
- **Gene search:** Type a gene name to visualize its expression on the UMAP.

**What to look for:**
- Distinct, well-separated clusters suggest clear cell populations.
- If cells cluster by sample/batch instead of biology, consider enabling Harmony.
- Gene expression overlays help validate cluster identities (e.g., CD3D should light up in T cell clusters).

### Marker Genes Tab

Three visualizations:

- **Ranking plot:** Top 5 marker genes per cluster, ranked by score.
- **Heatmap:** Expression of top marker genes across all cells, grouped by cluster.
- **Interactive table:** Select a cluster and view its marker genes with statistics.

**What to look for:**
- Each cluster should have distinct, non-overlapping marker genes.
- Known marker genes should appear in the expected clusters.
- Very high fold-changes (> 2-3) with low adjusted p-values indicate strong markers.

---

## Downloading Results

### Available Downloads

| Download | Format | Content |
|----------|--------|---------|
| **H5AD** | .h5ad | Complete AnnData object with all analysis results (embeddings, clusters, markers) |
| **Cell metadata** | .csv | Full `obs` table (cluster assignments, QC metrics, sample info) |
| **Marker genes** | .csv | Top 25 marker genes per cluster with statistics |

### Plot Downloads

Every plot offers PNG (300 DPI) and SVG download buttons.

### Using the H5AD downstream

The downloaded `.h5ad` file contains everything from the pipeline and can be opened in:
- **Python:** `import scanpy as sc; adata = sc.read_h5ad("scrna_analysis.h5ad")`
- **R:** Using the `anndata` R package or converting to Seurat format
- **Cell type annotation tools:** CellTypist, SingleR, ScType

---

## RAM Optimization Guide — In Depth

This section provides practical strategies for running the scRNA-seq pipeline on machines with limited RAM (8 GB). These strategies are ordered by effectiveness — **the first one is the most impactful**.

### Strategy 1: Aggressive Early Filtering (Most Effective)

The biggest memory savings come from removing cells and genes **early** in the pipeline. The expression matrix (cells × genes) is the dominant consumer of RAM, so reducing either dimension has multiplicative effects.

**Aggressive QC thresholds:**

| Parameter | Default | Aggressive (8 GB) | Why |
|-----------|---------|-------------------|-----|
| Min genes per cell | 200 | 500 | Removes more empty droplets and debris |
| Max genes per cell | 5,000 | 3,000 | Catches more doublets |
| Min counts per cell | 500 | 1,000 | Removes more low-quality cells |
| Max % mitochondrial | 20% | 10% | Removes more dying cells |
| Min cells per gene | 3 | 10 | Removes more noise genes |

**Impact:** Filtering a 50,000-cell dataset down to 15,000 cells can reduce peak RAM by **60-70%**.

### Strategy 2: Subsample for Exploratory Runs

If you are tuning parameters (resolution, n_neighbors, etc.), consider subsampling your `.h5ad` before uploading:

```python
import scanpy as sc

adata = sc.read_h5ad("full_dataset.h5ad")

# Random subsample to 10k cells for fast exploration
sc.pp.subsample(adata, n_obs=10000)

adata.write_h5ad("subset_10k.h5ad")
```

Run the full dataset only after you have found satisfactory parameters on the subset.

### Strategy 3: Reduce HVG Count

The default 2,000 HVGs works well for most datasets. Reducing to 1,000-1,500:
- Saves memory during scaling (the dense matrix is smaller).
- Saves memory during PCA and Harmony.
- Minimal impact on biological signal — most cell type differences are captured in the top 1,000 genes.

### Strategy 4: Fewer Principal Components

Reducing from 50 PCs to 30 PCs saves memory during:
- PCA computation
- Harmony batch correction
- Neighbor graph construction

Most biological variation is captured in the first 20-30 PCs. Use the elbow plot to decide.

### Strategy 5: Disable Optional Steps

- **Doublet detection (Scrublet):** Temporarily creates a dense matrix for simulation. Disable it if your data was already filtered externally.
- **Batch correction (Harmony):** Adds a subprocess with its own memory footprint (~1-2 GB extra temporarily). Only enable it if your data genuinely contains multiple batches.

### Strategy 6: Close Other Applications

Browsers (especially Chrome with multiple tabs), IDEs (VS Code, PyCharm), and communication apps (Slack, Discord) all compete for RAM. On an 8 GB machine, close unnecessary programs before running the pipeline.

### What the Pipeline Already Does to Save RAM

You don't need to do anything special for these — they're built into the code:

1. **Single-mask filtering:** All cell QC filters are combined into one boolean mask and applied in a single `.copy()`, avoiding 4 intermediate copies.
2. **Sparse matrix preservation:** The expression matrix stays in sparse (CSR) format throughout QC, filtering, and normalization. Densification only happens at scaling, and only for the HVG subset (~2K genes, not 20K+).
3. **Early HVG subsetting:** Genes are subset to HVGs *before* scaling and PCA, so the dense matrix is ~10x smaller than the full gene set.
4. **Subprocess isolation for Harmony:** PyTorch + Harmony run in a `spawn` child process. When the subprocess exits, the OS reclaims all of its memory — nothing leaks into the main application.
5. **Session state cleanup:** Previous analysis results are freed from memory before each new pipeline run.
6. **Aggressive garbage collection:** `gc.collect()` is called after every major pipeline step (filtering, doublets, normalization, scaling, PCA, batch correction, UMAP).

---

## Frequently Asked Questions

### Q: The pipeline crashed without an error message. What happened?

Almost certainly an **out-of-memory (OOM) event**. The operating system killed the process. See [RAM Optimization Guide](#ram-optimization-guide--in-depth) for solutions.

### Q: My data clusters by sample/donor instead of by cell type. What should I do?

Enable **Harmony batch correction** (Step 9):
1. In the sidebar, check "Enable batch correction (Harmony)".
2. Select the column that identifies the sample/donor.
3. Re-run the pipeline.

Harmony will adjust the PCA embedding so that cells from different samples/donors integrate properly, while preserving biological differences.

### Q: How do I know if my filtering thresholds are too aggressive or too lenient?

Check the **QC violin plots** (QC tab):
- If the violins show a clear bimodal distribution, your thresholds should separate the two peaks.
- If you're removing > 50% of cells, your thresholds may be too aggressive.
- If the UMAP shows a large cluster of low-quality cells, your thresholds may be too lenient.

### Q: What resolution should I use for Leiden clustering?

There's no universal answer. Guidelines:
- **0.1-0.3:** Broad cell types (T cells, B cells, myeloid cells)
- **0.5 (default):** Good starting point for most datasets
- **0.8-1.5:** Sub-populations (CD4+ vs. CD8+ T cells, naive vs. memory)
- **2.0+:** Very fine resolution — risk of over-clustering

Start with 0.5, look at the UMAP, and adjust up or down.

### Q: Can I re-run the pipeline with different parameters without re-uploading?

**Yes.** Change the sidebar parameters and click "Run Pipeline" again. The H5AD data stays cached — only the pipeline re-runs. Each run clears the previous results from memory.

### Q: Why does the pipeline re-annotate QC even if my H5AD already has annotations?

The pipeline annotates from scratch to ensure consistency with its own mitochondrial/ribosomal/hemoglobin gene detection patterns (prefix-based: `MT-`, `RPS`, `RPL`, `HB`). Different tools use different detection patterns, so re-annotating ensures the filtering thresholds work correctly.

### Q: How do I identify cell types from the clusters?

After running the pipeline:
1. Go to the **Marker Genes** tab.
2. For each cluster, look at the top marker genes.
3. Compare with known marker gene databases:
   - [CellMarker](http://bio-bigdata.hrbmu.edu.cn/CellMarker/)
   - [PanglaoDB](https://panglaodb.se/)
   - Literature from your tissue/organism of interest.

Example:
- CD3D, CD3E → T cells
- CD79A, MS4A1, CD19 → B cells
- CD14, LYZ, CST3 → Monocytes
- NKG7, GNLY → NK cells
- COL1A1, COL1A2 → Fibroblasts

### Q: Can I use the output H5AD in R (Seurat)?

Yes. You can convert the H5AD to a Seurat object using the `anndata` R package or the `SeuratDisk` package:

```r
library(SeuratDisk)
Convert("scrna_analysis.h5ad", dest = "h5seurat")
sobj <- LoadH5Seurat("scrna_analysis.h5seurat")
```

### Q: What's the difference between Wilcoxon and t-test for marker genes?

- **Wilcoxon rank-sum** (default): Non-parametric. Does not assume normal distribution. More robust for scRNA-seq data, which is often zero-inflated.
- **t-test**: Parametric. Assumes normal distribution. Faster but can be less accurate for sparse single-cell data.
- **logreg**: Logistic regression. Identifies genes that best *discriminate* each cluster from all others. Can be more specific but slower.

**Recommendation:** Use Wilcoxon (default) unless you have a specific reason to prefer another method.

---

## Glossary

| Term | Definition |
|------|-----------|
| **AnnData** | Annotated data matrix — the standard data structure for Scanpy. Contains the expression matrix (X), cell metadata (obs), gene metadata (var), embeddings (obsm), and analysis results (uns) |
| **Barcode** | A unique DNA sequence identifying each cell in a droplet-based experiment |
| **Batch effect** | Technical variation between experimental groups (samples, donors, sequencing runs) that can confound biological signal |
| **Clustering** | Grouping cells with similar expression profiles into clusters (cell populations) |
| **Doublet** | A droplet that captured two cells — appears as a single "cell" with abnormally high gene count and mixed expression profile |
| **H5AD** | HDF5-based file format for storing AnnData objects |
| **Harmony** | A batch correction algorithm that aligns cells across batches in PCA space |
| **HVG** | Highly Variable Gene — a gene with high cell-to-cell expression variability, likely driving biological differences |
| **kNN graph** | k-Nearest-Neighbor graph — a graph where each cell is connected to its k most similar cells |
| **Leiden** | Graph-based community detection algorithm for clustering cells |
| **Library size** | Total number of UMI counts in a cell — a proxy for sequencing depth |
| **Log1p** | Log(x + 1) transformation — reduces skewness of count data |
| **Marker gene** | A gene specifically expressed in one cluster/cell type compared to others |
| **Mitochondrial genes** | Genes encoded in the mitochondrial genome (MT- prefix). High expression indicates cell stress or death |
| **Normalization** | Adjusting expression values to account for technical differences (library size, sequencing depth) |
| **PCA** | Principal Component Analysis — dimensionality reduction that captures the main axes of variation |
| **Resolution** | Parameter controlling the number of Leiden clusters. Higher = more clusters |
| **Scanpy** | Python library for single-cell analysis |
| **Scrublet** | Tool for detecting doublets by simulating artificial doublets and comparing real cells to them |
| **SoupX** | Tool for removing ambient RNA contamination from droplet-based scRNA-seq data |
| **Sparse matrix** | Memory-efficient representation of a matrix with many zeros (typical for scRNA-seq) |
| **UMAP** | Uniform Manifold Approximation and Projection — 2D visualization preserving local neighborhood structure |
| **UMI** | Unique Molecular Identifier — a short barcode attached to each mRNA molecule to count unique transcripts |

---

*This tutorial is part of the [BeginSeq Studio](https://github.com/JACKNINES/beginseq-studio) documentation.*
