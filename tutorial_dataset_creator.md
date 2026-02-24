# Tutorial: Dataset Creator — GDC/TCGA Downloader

## Table of Contents

1. [Overview](#overview)
2. [What is GDC/TCGA?](#what-is-gdctcga)
3. [What You Need Before Starting](#what-you-need-before-starting)
4. [Sidebar — Theme and Language](#sidebar--theme-and-language)
5. [Step 1 — Load TCGA Projects](#step-1--load-tcga-projects)
6. [Step 2 — Select a Project and Fetch Files](#step-2--select-a-project-and-fetch-files)
7. [Step 3 — Filter and Configure](#step-3--filter-and-configure)
8. [Step 4 — Download and Build](#step-4--download-and-build)
9. [Step 5 — Preview and Download Results](#step-5--preview-and-download-results)
10. [Using Your Dataset in the Bulk RNA-seq Module](#using-your-dataset-in-the-bulk-rna-seq-module)
11. [Understanding the Output Files](#understanding-the-output-files)
12. [Practical Examples](#practical-examples)
13. [Frequently Asked Questions](#frequently-asked-questions)
14. [Glossary](#glossary)

---

## Overview

The Dataset Creator module allows you to **build analysis-ready RNA-seq datasets** directly from the [NCI Genomic Data Commons (GDC)](https://gdc.cancer.gov/) — without writing code, without using R, and without installing TCGAbiolinks or any command-line tools.

### What it does

1. Connects to the GDC REST API to list all available TCGA projects.
2. Lets you select a cancer type (e.g., TCGA-BRCA for breast cancer).
3. Downloads STAR-Counts RNA-seq files from the GDC.
4. Assembles a **DESeq2-ready count matrix** and **metadata file**.
5. Saves both files locally for immediate use in the Bulk RNA-seq module.

### Why this is useful

Traditionally, downloading and assembling TCGA data requires:
- Installing R and TCGAbiolinks (or GDCquery)
- Writing scripts to query the API, download files, and parse them
- Handling complex file formats (STAR-Counts TSV files with multiple count columns)
- Manual metadata assembly from GDC case information

This module does all of that automatically through a visual interface.

---

## What is GDC/TCGA?

### The Cancer Genome Atlas (TCGA)

TCGA is the largest public cancer genomics dataset in the world. It contains:
- **33 cancer types** (breast, lung, colon, brain, etc.)
- **~11,000 patients** with matched tumor and normal tissue samples
- **Multiple data types:** RNA-seq, DNA sequencing, methylation, copy number, etc.

### The Genomic Data Commons (GDC)

The GDC is the NIH data portal that hosts TCGA (and other) datasets. It provides:
- A REST API for programmatic access
- Open-access data (available to anyone, no special permissions needed)
- Standardized processing pipelines (all RNA-seq data processed with STAR)

### What data does this module access?

This module downloads **open-access RNA-seq STAR-Counts files** from TCGA projects. These are:
- Raw gene-level read counts (integers)
- Aligned and quantified using the STAR aligner
- Available for both tumor and normal tissue samples
- Suitable for differential expression analysis with DESeq2

> **Important:** This module only accesses **open-access** data. It does NOT access controlled-access data, and does NOT require dbGaP authorization or any special permissions.

---

## What You Need Before Starting

### Requirements

- **Internet connection** — the module downloads data from the GDC API.
- **Disk space** — downloads are saved to a directory on your machine. See estimates below.
- **Time** — downloading depends on your internet speed and the number of files.

### No special software needed

Unlike traditional TCGA access methods, you do NOT need:
- ❌ R or Bioconductor
- ❌ TCGAbiolinks
- ❌ Command-line tools
- ❌ GDC Data Transfer Tool
- ❌ dbGaP authorization (only open-access data is used)

### Disk space estimates

| Project size | Approximate download size |
|:-------------|:--------------------------|
| 50 samples | ~50-100 MB |
| 200 samples | ~200-400 MB |
| 500 samples | ~500 MB–1 GB |
| 1,000+ samples | ~1-2 GB |

Files are downloaded as compressed archives and extracted locally. The extracted TSV files are kept in the download directory for reference.

---

## Sidebar — Theme and Language

Before starting, you can customize the interface from the sidebar:

### Language (🌐)

Switch between **English** and **Spanish**. All labels and messages update instantly.

### Visual Theme (🎨)

Choose a visual theme for the interface:

| Theme | Description |
|-------|-------------|
| 🌑 **Dark** (default) | Dark background, comfortable for extended use |
| ☀️ **Light** | White background, best for readability |
| ⚡ **Cyberpunk** | Deep black with neon accents |

The Dataset Creator does not generate analysis plots, but the theme applies to the entire interface (buttons, tables, progress indicators, etc.) and carries over when you navigate to the Bulk RNA-seq module.

---

## Step 1 — Load TCGA Projects

### What to do

1. Open the **Dataset Creator** page from the sidebar (📦 Dataset Creator).
2. Click **"Load TCGA Projects"**.

### What happens

The app queries the GDC API for all available TCGA projects with RNA-seq data. This typically takes 2-5 seconds depending on your internet connection.

### What you see

A success message showing the number of available projects (usually ~33), and an expandable table with:

| Column | Description |
|--------|-------------|
| **project_id** | TCGA project identifier (e.g., `TCGA-BRCA`) |
| **name** | Full project name (e.g., `Breast Invasive Carcinoma`) |
| **primary_site** | Anatomical site (e.g., `Breast`) |
| **file_count** | Number of RNA-seq files available |

### Why this step exists

The GDC hosts hundreds of projects, but not all have RNA-seq data. This step filters to only TCGA projects with STAR-Counts files available, so you can focus on relevant datasets.

---

## Step 2 — Select a Project and Fetch Files

### What to do

1. Use the dropdown to select a TCGA project (e.g., `TCGA-BRCA — Breast Invasive Carcinoma`).
2. Click **"Fetch Files"**.

### What happens

The app queries the GDC API for all RNA-seq STAR-Counts files in the selected project. This returns detailed information about each file, including:
- File ID (unique identifier for downloading)
- Case ID (patient identifier)
- Sample type (e.g., "Primary Tumor", "Solid Tissue Normal")
- Condition (derived from sample type: "Tumor" or "Normal")
- File size

### What you see

After fetching:

- **Summary metrics:**
  - Total files found
  - Number of Tumor samples
  - Number of Normal samples
  - Total download size (MB)

- **Expandable file table** with detailed information per file.

### Understanding sample types

The GDC classifies tissue samples using standardized sample type codes:

| Sample Type | Condition Label | Description |
|------------|----------------|-------------|
| Primary Tumor | Tumor | Main tumor tissue |
| Recurrent Tumor | Tumor | Recurrent tumor |
| Metastatic | Tumor | Metastatic tissue |
| Solid Tissue Normal | Normal | Adjacent normal tissue |
| Blood Derived Normal | Normal | Normal blood sample |

> **Note:** Not all projects have both tumor and normal samples. Some projects (e.g., rare cancers) may have very few normal samples or none at all.

### What to check

- **Sample balance:** Ideally, you want a reasonable number of both Tumor and Normal samples. Extreme imbalances (e.g., 800 Tumor vs. 5 Normal) can reduce statistical power.
- **Total file count:** Very large projects (500+ files) will take longer to download and process.

---

## Step 3 — Filter and Configure

This step lets you customize which files to download and how to build the count matrix.

### Filter by Condition

By default, both Tumor and Normal samples are selected. You can:
- Deselect "Normal" to download only tumor samples (e.g., for tumor subtype analysis).
- Deselect "Tumor" to download only normal samples (e.g., for a normal tissue reference).

**Why this is useful:** If you only need tumor data (e.g., for survival analysis or subtype classification), skipping normal samples saves download time and disk space.

### Sample Limit (Optional)

**When to use:** When the project has hundreds of files and you want to:
- **Test the workflow** with a small subset first.
- **Reduce download time** for a quick exploratory analysis.
- **Balance conditions** to have equal numbers of Tumor and Normal samples.

#### Two limiting modes:

**1. Total limit:** Set a maximum total number of files. Files are sampled proportionally across conditions. For example, if the project has 800 Tumor and 100 Normal, a total limit of 100 would give ~89 Tumor and ~11 Normal (proportional).

**2. Per-condition limit:** Set a maximum per condition. For example, 50 per condition = 50 Tumor + 50 Normal = 100 total. This creates a perfectly balanced dataset.

> **Tip:** For your first analysis of a large project, try **per-condition = 20-30**. This gives you a manageable dataset (40-60 files) that downloads quickly and runs fast in the Bulk RNA-seq module.

### Gene ID Type

Choose how genes are identified in the output count matrix:

| Option | Example | When to use |
|--------|---------|-------------|
| **gene_name** (default) | TP53, BRCA1, MYC | Most intuitive. Use this unless you have a specific need for Ensembl IDs. |
| **gene_id** | ENSG00000141510, ENSG00000012048 | Needed for tools that require Ensembl gene IDs. More precise (no name ambiguity). |

**Why it matters:** Gene names are human-readable but can have duplicates (e.g., some genes share names). Ensembl gene IDs are unique and unambiguous. For most users, gene names are the better choice.

### Count Column

STAR produces multiple count columns for each gene. Choose which one to use:

| Column | Description | When to use |
|--------|-------------|-------------|
| **unstranded** (default) | Total counts regardless of strand | Use this for most RNA-seq library preparations |
| **stranded_first** | Counts from the first read strand | Use for strand-specific libraries (e.g., dUTP method) |
| **stranded_second** | Counts from the second read strand | Use for strand-specific libraries (reverse strand) |

**If you're not sure:** Use **unstranded**. It works for both stranded and unstranded libraries, though stranded-specific counting is slightly more accurate for stranded data.

### Why these options matter

The gene ID type and count column affect the downstream analysis:
- **Gene ID type** determines how you'll search for and interpret genes in the Bulk RNA-seq module.
- **Count column** must match your library preparation — using the wrong strand can halve your apparent expression levels.

---

## Step 4 — Download and Build

### Download Directory

The app shows a default download directory (usually `~/Desktop/gdc_downloads`). You can change this to any directory on your machine.

**What gets saved there:**
- The raw GDC tar.gz archives (downloaded from the API).
- Extracted STAR-Counts TSV files (one per sample).
- These files persist after the app closes — they are NOT temporary.

> **Why files are kept:** If you need to re-run the matrix building step with different settings (e.g., different gene ID type), the app can reuse the already-downloaded files without re-downloading them from the GDC.

### Starting the Download

Click **"Start Download"** to begin. The process has two phases:

#### Phase 1: Download and Extract (~50% of progress bar)

- The app batches files into groups and downloads them as tar.gz archives from the GDC REST API.
- Each archive is extracted to the download directory.
- The progress bar shows: "Downloading file X of Y".

**Download speed** depends on your internet connection and GDC server load:
- 50 files: ~1-3 minutes
- 200 files: ~5-15 minutes
- 800+ files: ~20-60 minutes

#### Phase 2: Build Matrix (~50% of progress bar)

- The app parses each STAR-Counts TSV file.
- Extracts the selected gene ID and count column.
- Assembles all samples into a single count matrix.
- Generates metadata (sample ID, case ID, condition).

**Build speed** depends on the number of files and genes:
- 50 files: ~10-30 seconds
- 200 files: ~30-90 seconds
- 800+ files: ~2-5 minutes

### Large download warning

If you selected more than 200 files, the app shows a warning that the download may take a while. Consider using a **sample limit** (Step 3) if this is an exploratory run.

---

## Step 5 — Preview and Download Results

After the build completes, you see:

### Success Summary

> ✅ Built count matrix: **X genes × Y samples**

### Count Matrix Preview

The first 20 rows of the count matrix. Check that:
- Gene names (or IDs) look correct.
- Sample column names are meaningful (case IDs).
- Values are integers (raw counts).

### Metadata Preview

The full metadata table. Each row is a sample with:

| Column | Content |
|--------|---------|
| **sample** | Sample identifier (used as column name in the count matrix) |
| **case_id** | Patient identifier |
| **condition** | "Tumor" or "Normal" |
| **sample_type** | Original GDC sample type (e.g., "Primary Tumor") |

### Condition Summary

Shows how many samples are in each condition. For example:
- **Tumor**: 120
- **Normal**: 15

### Download Buttons

| Button | File | Format |
|--------|------|--------|
| **📥 Download Count Matrix** | `counts_matrix.csv` | CSV (genes as rows, samples as columns) |
| **📥 Download Metadata** | `metadata.csv` | CSV (one row per sample) |

### Tip for Bulk RNA-seq

At the bottom, a tip reminds you that these files can be directly uploaded to the Bulk RNA-seq module for differential expression analysis.

---

## Using Your Dataset in the Bulk RNA-seq Module

The Dataset Creator produces files in **exactly the format** the Bulk RNA-seq module expects. Here's a quick walkthrough:

### 1. Download both files

- `counts_matrix.csv` — the count matrix
- `metadata.csv` — the sample metadata

### 2. Open the Bulk RNA-seq module

Navigate to 🧬 **Bulk RNA-seq** from the sidebar.

### 3. Upload the files

- Upload `counts_matrix.csv` as the **Count Matrix**.
- Upload `metadata.csv` as the **Metadata**.

### 4. Configure columns

- **Sample column:** Select `sample` (auto-detected).
- **Condition column:** Select `condition` (auto-detected).

### 5. Select the contrast

- **Reference level:** `Normal`
- **Test level:** `Tumor`

### 6. Enable gene pre-filtering (recommended)

TCGA count matrices typically have ~60,000 genes. Most are lowly expressed. Gene pre-filtering is **highly recommended**:
- It reduces the gene count from ~60K to ~20-25K.
- This cuts runtime by 60-70% and improves statistical power.
- Default filter settings work well for TCGA data.

### 7. Run the analysis

Click "Run Analysis" and wait. With 200 samples and gene filtering enabled, expect ~3-10 minutes on a modern laptop.

---

## Understanding the Output Files

### Count Matrix (`counts_matrix.csv`)

```
,TCGA-A7-A13E-01A,TCGA-A7-A13F-01A,TCGA-BH-A0BJ-11A,...
TP53,1205,982,1547,...
BRCA1,342,410,289,...
MYC,5621,4890,6102,...
...
```

- **First column (index):** Gene names or Ensembl IDs.
- **Column headers:** Sample identifiers (derived from TCGA barcode).
- **Values:** Raw integer counts from STAR alignment.

**Important:** The first 4 rows of STAR-Counts files contain summary statistics (`N_unmapped`, `N_multimapping`, `N_noFeature`, `N_ambiguous`). The Dataset Creator **automatically removes** these rows — only gene counts are included.

### Metadata (`metadata.csv`)

```
sample,case_id,condition,sample_type
TCGA-A7-A13E-01A,TCGA-A7-A13E,Tumor,Primary Tumor
TCGA-A7-A13F-01A,TCGA-A7-A13F,Tumor,Primary Tumor
TCGA-BH-A0BJ-11A,TCGA-BH-A0BJ,Normal,Solid Tissue Normal
...
```

- **sample:** Matches column names in the count matrix.
- **case_id:** Patient identifier (same patient may have both tumor and normal samples).
- **condition:** Simplified to "Tumor" or "Normal" for DESeq2.
- **sample_type:** Original GDC sample type classification.

---

## Practical Examples

### Example 1: Breast Cancer (TCGA-BRCA)

One of the largest TCGA projects with excellent tumor/normal representation.

1. Load projects → Select `TCGA-BRCA — Breast Invasive Carcinoma`.
2. Fetch files → ~1,200 files, ~1,100 Tumor, ~112 Normal.
3. For a first analysis: Set per-condition limit to 30 (60 files total).
4. Gene ID: `gene_name`, Count column: `unstranded`.
5. Download → ~60 MB, ~2 minutes.
6. Build → ~30 seconds.
7. Upload to Bulk RNA-seq → Enable gene filtering → Run.

### Example 2: Lung Adenocarcinoma (TCGA-LUAD)

1. Load projects → Select `TCGA-LUAD — Lung Adenocarcinoma`.
2. Fetch files → ~600 files, ~530 Tumor, ~59 Normal.
3. Keep all samples for a comprehensive analysis.
4. Download → ~600 MB, ~15 minutes.
5. Build → ~90 seconds.
6. Upload to Bulk RNA-seq → **Strongly recommend gene filtering** (60K genes × 600 samples is large).

### Example 3: Quick Exploratory Analysis

If you want to quickly explore any TCGA project:

1. Select any project.
2. Fetch files.
3. Set total limit to **20 samples** (fastest possible).
4. Download → ~20 MB, ~30 seconds.
5. Build → ~10 seconds.
6. Upload to Bulk RNA-seq → Runs in ~30-60 seconds.

> **Note:** With only 20 samples, statistical power is limited. Use this for workflow testing, not for drawing biological conclusions.

---

## Frequently Asked Questions

### Q: How long does the download take?

It depends on your internet speed and the number of files:
- 50 files: 1-3 minutes
- 200 files: 5-15 minutes
- 800+ files: 20-60 minutes

The GDC API can sometimes be slow during peak hours (US business hours). If downloads are very slow, try again later.

### Q: What if the download fails partway through?

The app downloads files in batches. If a batch fails, the error message will tell you which batch had a problem. You can:
1. Check your internet connection.
2. Retry the download — already-extracted files won't be re-downloaded.

### Q: Can I use controlled-access GDC data?

**No.** This module only accesses **open-access** data via the public GDC REST API. Controlled-access data requires dbGaP authorization and the GDC Data Transfer Tool, which are beyond the scope of this app.

### Q: What if a project has no Normal samples?

Some TCGA projects (especially rare cancers) have very few or no matched normal samples. In this case:
- You can still download the tumor data.
- For DE analysis, you would need to find normal reference samples elsewhere (e.g., from GTEx or another TCGA project).
- The Dataset Creator will show 0 Normal samples in the metrics.

### Q: Why does the count matrix have ~60,000 genes?

STAR quantifies against the full GENCODE annotation, which includes ~60,000 entries:
- ~20,000 protein-coding genes
- ~15,000 lncRNAs
- ~10,000 pseudogenes
- ~15,000 other features (rRNAs, miRNAs, etc.)

Most of these are lowly expressed or not expressed at all. Use **gene pre-filtering** in the Bulk RNA-seq module to focus on the ~20-25K genes that matter.

### Q: Can I download data from non-TCGA projects (e.g., TARGET, CGCI)?

Currently, the module is designed for **TCGA projects** only. The API queries filter for TCGA project IDs. Other GDC programs (TARGET, CGCI, etc.) use different data structures and would require modifications.

### Q: Are the downloaded files kept on my machine?

**Yes.** The extracted STAR-Counts files are saved in the download directory you specified. They are NOT deleted when the app closes. This is intentional — it lets you:
- Re-build the matrix with different settings.
- Access the raw data for other analyses.
- Share the files with collaborators.

If you no longer need them, you can manually delete the download directory.

### Q: What happens if the same patient has multiple tumor samples?

Some TCGA patients have multiple samples (e.g., primary tumor + metastatic tumor). Each sample gets its own row in the metadata and its own column in the count matrix. The `case_id` column lets you identify samples from the same patient.

For standard DE analysis (Tumor vs. Normal), this is usually fine. For more sophisticated analyses (e.g., paired tumor/normal from the same patient), you would need to handle the pairing in a custom DESeq2 design formula (not currently supported in the Bulk RNA-seq module).

### Q: How do I cite TCGA data?

If you use TCGA data in a publication, cite:

1. The **original TCGA project** for your cancer type.
2. The **GDC**: Grossman, R.L. et al. *N Engl J Med* 375, 1109–1112 (2016).
3. **BeginSeq Studio** if you used this tool for download/analysis.

See the main [README](README.md) for full citation information.

---

## Glossary

| Term | Definition |
|------|-----------|
| **Case ID** | Patient identifier in TCGA (e.g., `TCGA-A7-A13E`) |
| **Condition** | Simplified biological grouping: "Tumor" or "Normal" |
| **Count matrix** | Table of raw integer read counts (genes × samples) |
| **DESeq2-ready** | A count matrix + metadata file in the format required by DESeq2 for differential expression analysis |
| **GDC** | Genomic Data Commons — NIH data portal hosting TCGA and other cancer genomics data |
| **Gene ID (Ensembl)** | Unique identifier for a gene (e.g., `ENSG00000141510` for TP53) |
| **Gene name** | Human-readable gene symbol (e.g., `TP53`, `BRCA1`) |
| **GENCODE** | The gene annotation catalog used by STAR — includes ~60,000 features |
| **Metadata** | Sample-level information: which samples belong to which conditions |
| **Open-access** | Data available to anyone without special permissions |
| **REST API** | Programming interface for accessing GDC data over the internet |
| **Sample type** | GDC classification of tissue type (e.g., "Primary Tumor", "Solid Tissue Normal") |
| **STAR** | Spliced Transcripts Alignment to a Reference — the RNA-seq aligner used by GDC/TCGA |
| **STAR-Counts** | Gene-level count files produced by STAR alignment |
| **TCGA** | The Cancer Genome Atlas — the largest public cancer genomics dataset |
| **Unstranded** | Read counts regardless of which DNA strand the read maps to |

---

*This tutorial is part of the [BeginSeq Studio](https://github.com/JACKNINES/beginseq-studio) documentation.*
