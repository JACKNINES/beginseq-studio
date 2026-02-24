# Tutorial: Bulk RNA-seq — Differential Expression Analysis

## Table of Contents

1. [Overview](#overview)
2. [What You Need Before Starting](#what-you-need-before-starting)
3. [Sidebar — Theme and Language](#sidebar--theme-and-language)
4. [Step 1 — Upload Your Files](#step-1--upload-your-files)
5. [Step 2 — Preview and Validate Data](#step-2--preview-and-validate-data)
6. [Step 3 — Configure Metadata Columns](#step-3--configure-metadata-columns)
7. [Step 4 — Filter Samples (Optional)](#step-4--filter-samples-optional)
8. [Step 5 — Select the Contrast](#step-5--select-the-contrast)
9. [Step 6 — Expression-Based Classification (Optional)](#step-6--expression-based-classification-optional)
10. [Step 7 — Gene Pre-Filtering (Optional but Recommended)](#step-7--gene-pre-filtering-optional-but-recommended)
11. [Step 8 — Batch Correction for PCA (Optional)](#step-8--batch-correction-for-pca-optional)
12. [Step 9 — Run the Analysis](#step-9--run-the-analysis)
13. [Step 10 — Interpret Results](#step-10--interpret-results)
14. [Step 11 — Visualize Results](#step-11--visualize-results)
15. [Step 12 — Highlight Genes of Interest](#step-12--highlight-genes-of-interest)
16. [Step 13 — Download Results](#step-13--download-results)
17. [Step 14 — Audit Log (Reproducibility)](#step-14--audit-log-reproducibility)
18. [Frequently Asked Questions](#frequently-asked-questions)
19. [Glossary](#glossary)

---

## Overview

The Bulk RNA-seq module performs **differential expression (DE) analysis** on raw count data using the DESeq2 algorithm (Python implementation via [PyDESeq2](https://github.com/owkin/PyDESeq2)).

### What does it do?

Given two experimental conditions (e.g., "Tumor" vs. "Normal", "Treated" vs. "Control"), the pipeline identifies which genes are **significantly up-regulated or down-regulated** between them. It performs:

1. **Median-of-ratios normalization** — adjusts for library size differences between samples.
2. **Per-gene dispersion estimation** — models the biological variability of each gene.
3. **Negative binomial GLM fitting** — fits a generalized linear model to each gene.
4. **Wald test** — tests whether the fold-change is significantly different from zero.
5. **Benjamini-Hochberg correction** — controls the false discovery rate across all genes.
6. **Visualization** — generates volcano plots, MA plots, PCA, and heatmaps.

### Why DESeq2? (I personally prefer to work with python over R)

DESeq2 is one of the most widely used and well-validated tools for bulk RNA-seq DE analysis. It correctly handles:

- Count data (integers, not normalized values)
- Small sample sizes (as few as 3 replicates per group)
- Overdispersed data (biological variability beyond Poisson noise)
- Multiple testing correction

---

## What You Need Before Starting

### Required Files

You need exactly **two files**:

#### 1. Count Matrix

A table where:
- **Rows** = genes (gene IDs or gene names)
- **Columns** = samples
- **Values** = raw, unnormalized integer counts

| gene | Sample_1 | Sample_2 | Sample_3 | Sample_4 |
|------|----------|----------|----------|----------|
| TP53 | 1205 | 982 | 1547 | 1103 |
| BRCA1 | 342 | 410 | 289 | 376 |
| MYC | 5621 | 4890 | 6102 | 5234 |

**Accepted formats:** `.csv`, `.tsv`, `.zip`

> **Important:** The values must be **raw counts** (integers like 0, 15, 1205), NOT normalized values (like FPKM, TPM, or RPKM). DESeq2 performs its own normalization internally. If you upload normalized data, the app will show a warning.

**Where to get count data:**
- **STAR aligner** output (`ReadsPerGene.out.tab`)
- **featureCounts** output
- **HTSeq** count files
- **The Dataset Creator module** in this app (downloads from TCGA/GDC)

#### 2. Metadata File

A table that maps each sample to its experimental condition:

| sample | condition |
|--------|-----------|
| Sample_1 | Control |
| Sample_2 | Control |
| Sample_3 | Treated |
| Sample_4 | Treated |

**Accepted formats:** `.csv`, `.zip`

**Requirements:**
- Must contain at least two columns: one for sample IDs and one for conditions.
- Sample IDs must match the column names in the count matrix.
- There must be at least 2 different conditions.
- Each condition should have at least 2 samples (ideally 3+).

> **Tip:** You can include additional columns (e.g., batch, sex, age). These can be used for sample filtering or batch correction in PCA.

### File Size Limits

When running locally, each file can be up to **2 GB**. This is more than enough for any bulk RNA-seq experiment.

---

## Sidebar — Theme and Language

Before starting your analysis, you can customize the interface from the sidebar:

### Language (🌐)

Switch between **English** and **Spanish**. All labels, instructions, and messages update instantly.

### Visual Theme (🎨)

Choose a visual theme that applies to both the interface AND all generated plots:

| Theme | Description |
|-------|-------------|
| 🌑 **Dark** (default) | Dark background with soft-contrast colors. Comfortable for extended use. |
| ☀️ **Light** | White background with high-contrast colors. Best for print-ready figures. |
| ⚡ **Cyberpunk** | Deep black with neon accents. High-saturation for presentations on dark backgrounds. |

The theme persists across all pages during your session. Choose the theme **before** running the analysis — all plots will be generated using the selected palette.

> **Tip:** If you plan to include plots in a publication or poster with a white background, use the **Light** theme. For presentations on a projector, **Dark** or **Cyberpunk** work best.

---

## Step 1 — Upload Your Files

1. Open the **Bulk RNA-seq** page from the sidebar (🧬 Bulk RNA-seq).
2. Use the first file uploader to upload your **count matrix** (CSV or TSV).
3. Use the second file uploader to upload your **metadata** file (CSV).

### What happens behind the scenes

- The app reads both files and auto-detects the delimiter (comma for CSV, tab for TSV).
- Non-numeric columns in the count matrix (e.g., `gene_name`, `gene_type` from STAR files) are automatically dropped — only numeric count columns are kept.
- Genes with **zero counts across all samples** are removed (they carry no information and waste computational resources).

### Why this matters

DESeq2 requires raw integer counts. By automatically removing non-numeric columns and zero-count genes, the app ensures the data is in the correct format before analysis begins.

---

## Step 2 — Preview and Validate Data

After uploading, expand the preview sections to verify your data:

### Count Matrix Preview (👀)

- Shows the first 5 rows of your count matrix.
- Displays: **X genes × Y samples**.
- Reports how many zero-count genes were removed.

### Normalized Data Warning (⚠️)

If the data looks like it might already be normalized (fractional values, suspiciously uniform distributions), the app will show a warning. This is a **heuristic** — if you are confident your data is raw counts, you can ignore it.

### Metadata Preview (👀)

- Shows the first 10 rows of your metadata file.
- Displays the number of rows and columns.

### What to check

- **Sample counts match:** The number of columns in the count matrix should roughly match the number of rows in the metadata.
- **Gene names make sense:** Look for familiar gene names (TP53, GAPDH, ACTB, etc.).
- **Conditions are correct:** The metadata should show the expected groups (e.g., "Tumor" and "Normal").

---

## Step 3 — Configure Metadata Columns

### Sample ID Column (📋)

Select which column in your metadata contains the **sample identifiers**. These must match the column names in your count matrix.

The app auto-detects common names: `sample`, `Sample`, `sample_id`, `SampleID`, `id`.

### Condition Column (🏷️)

Select which column defines the **experimental conditions** you want to compare (e.g., `condition`, `group`, `treatment`, `disease`).

### Why this matters

The app needs to know:
1. Which metadata rows correspond to which count matrix columns (sample IDs).
2. Which grouping variable defines the comparison for differential expression.

### What if samples don't match?

If some samples appear in the count matrix but not in the metadata (or vice versa), the app shows a detailed mismatch report:
- How many samples are **only in the count matrix**
- How many samples are **only in the metadata**
- How many **match** (these will be used for analysis)

You can choose to proceed with only the matching samples by checking the "Proceed with common samples" checkbox.

---

## Step 4 — Filter Samples (Optional)

**When to use this:** If your metadata contains additional columns (e.g., `batch`, `sex`, `age_group`, `tissue_type`) and you want to restrict the analysis to a subset of samples.

### How to filter

1. Check ✅ **"Filter samples by metadata columns"**.
2. A filter interface appears where you can:
   - Select a column to filter by (e.g., `tissue_type`).
   - Choose which values to keep (e.g., keep only `"Lung"` and `"Liver"`, exclude `"Brain"`).
3. Click ➕ to add additional filters (e.g., also filter by `sex`).
4. Filters are applied **sequentially** — each filter narrows down the result of the previous one.

### Why this is useful

- **Focus on specific subsets:** Analyze only tumor samples from a specific tissue type.
- **Remove confounders:** Exclude samples from an unexpected batch or technical artifact.
- **Balanced comparisons:** Include only samples that are comparable.

> **Note:** Filtering happens before the analysis, reducing the number of samples DESeq2 needs to process.

---

## Step 5 — Select the Contrast

DESeq2 calculates fold-changes **relative to a reference level**. You need to choose:

### Reference Level (🔵)

The **baseline** condition. Fold-changes are calculated "relative to this." Think of it as the denominator.

- In a treatment study: this is usually `"Control"` or `"Untreated"`.
- In a cancer study: this is usually `"Normal"`.

### Test Level (🔴)

The condition being **compared**. Think of it as the numerator.

- In a treatment study: `"Treated"` or `"Drug_X"`.
- In a cancer study: `"Tumor"`.

### How to read the result

If reference = `Control` and test = `Treated`:
- A gene with **log2FC = 2** is **4x more expressed** in Treated vs. Control (up-regulated).
- A gene with **log2FC = -3** is **8x less expressed** in Treated vs. Control (down-regulated).
- A gene with **log2FC ≈ 0** shows no significant change.

### Contrast summary

The app shows a green box confirming your contrast:

> 📐 Contrast: **Treated** vs. **Control** (reference)

---

## Step 6 — Expression-Based Classification (Optional)

**When to use this:** When you want to **reclassify samples** based on the expression level of specific marker genes, creating compound conditions (e.g., `"Tumor_TP53_positive"` vs. `"Tumor_TP53_negative"`).

### How it works

1. Check ✅ **"Enable expression-based classification"**.
2. Enter one or more gene names (comma-separated or one per line).
3. Set a **threshold** in log2(CPM+1) units:
   - The app calculates log2(Counts-Per-Million + 1) for each gene in each sample.
   - Samples below (or above) the threshold are classified as "positive" or "negative".
4. Choose the **direction**:
   - "Below threshold" = classify as positive if expression is below the threshold.
   - "Above threshold" = classify as positive if expression is above the threshold.
5. Name the positive and negative labels (e.g., `"high"`, `"low"`).
6. Optionally keep reference samples unchanged (recommended).

### Expression statistics

Expand the 📊 **"Expression stats"** expander to see:
- Min, median, mean, standard deviation, and max expression for each gene.
- A z-score for your chosen threshold — helps you gauge how stringent your cutoff is.

### Why this is useful

- **Subgroup analysis:** Compare TP53-high tumors vs. TP53-low tumors, regardless of the original metadata grouping.
- **Biomarker-driven stratification:** Split samples by expression of a known biomarker.

### After classification

The conditions are updated to compound labels (e.g., `"Tumor_positive"`, `"Tumor_negative"`, `"Control"`). You can re-select the reference and test levels from these new conditions.

---

## Step 7 — Gene Pre-Filtering (Optional but Recommended)

### What it does

Before running DESeq2, you can remove lowly-expressed genes that are unlikely to be biologically meaningful. This is **standard practice** in RNA-seq analysis.

### Why filter genes?

1. **Reduces runtime:** Fewer genes = faster analysis (sometimes dramatically, from 60,000 genes to ~20,000).
2. **Improves statistical power:** Fewer genes = fewer multiple-testing corrections = more power to detect real differences.
3. **Stabilizes dispersion estimates:** Near-zero-expression genes have unreliable dispersion estimates that can distort the overall model.
4. **Reduces memory usage:** Especially important for large TCGA datasets (~800 samples × 60,000 genes).

### Parameters

| Parameter | Default | What it does |
|-----------|---------|-------------|
| **Min total count** | 10 | Gene must have ≥ 10 total counts across all samples |
| **Min samples expressing** | Auto | Gene must be expressed in ≥ this many samples. Auto mode uses max(smallest_group, n_samples × 0.5) |
| **Min count per sample** | 1 | A gene is "expressed" in a sample if it has ≥ this count |

### Recommendations

- For **small datasets** (< 30 samples): Use defaults. They are conservative enough.
- For **large datasets** (100+ samples, e.g., TCGA): Enable filtering. It can reduce 60K genes to ~20K, cutting runtime by 60-70%.
- For **exploratory analysis**: Use more aggressive filtering (min total = 50, min samples = 10%) to focus on robustly expressed genes.

> **Tip:** Filtered genes are completely removed before DESeq2 runs. They will NOT appear in the results table, so choose your thresholds thoughtfully.

---

## Step 8 — Batch Correction for PCA (Optional)

### When to use this

If your samples were processed in different **batches** (e.g., different sequencing runs, different labs, different dates) and you have a metadata column that identifies the batch, you can correct for batch effects in the PCA visualization.

### How it works

1. Check 🧪 **"Apply batch correction (ComBat) to PCA"**.
2. Select the batch column from your metadata.

### What this does

- Applies **ComBat** batch correction to the expression matrix **only for the PCA plot**.
- The DESeq2 analysis itself is NOT affected — only the visualization.
- This helps you see whether samples cluster by biology (condition) rather than by technical batch.

### When NOT to use this

- If you only have one batch (no batch information in metadata).
- If batch and condition are completely confounded (all Control samples in Batch 1, all Treated in Batch 2) — batch correction cannot separate the two effects.

---

## Step 9 — Run the Analysis

Click 🚀 **"Run Analysis"** to start the DESeq2 pipeline.

### What happens during the run

A progress bar shows each step with an estimated time remaining:

1. **Validation** — Checks data integrity.
2. **Gene filtering** — Removes low-expression genes (if enabled).
3. **Size factors** — Estimates library size normalization factors.
4. **Genewise dispersion** — Estimates per-gene variability (the slowest step).
5. **Dispersion trend** — Fits a parametric trend to stabilize estimates.
6. **MAP dispersion** — Shrinks per-gene dispersions toward the trend.
7. **Fit LFC** — Estimates log2 fold-changes via iterative reweighted least squares.
8. **Cook's distance** — Identifies outlier samples for each gene.
9. **Wald test** — Tests each gene for differential expression + BH correction.
10. **Visualizations** — Generates PCA, volcano, MA, and heatmap plots.

### Time estimates

| Dataset size | Approximate time |
|:-------------|:----------------|
| 20 samples × 15K genes | 30–60 seconds |
| 100 samples × 25K genes | 3–5 minutes |
| 800 samples × 25K genes | 15–30 minutes |
| 800 samples × 60K genes (no filter) | 45–90 minutes |

> **Tip:** Gene pre-filtering (Step 7) is the single most effective way to reduce runtime on large datasets.

### After completion

The progress bar reaches 100% and shows the total elapsed time. A step-by-step timing breakdown is available in an expander.

---

## Step 10 — Interpret Results

### Results Table (📊)

The main results table shows one row per gene, sorted by adjusted p-value:

| Column | Description |
|--------|-------------|
| **gene** (index) | Gene identifier |
| **baseMean** | Average normalized expression across all samples |
| **log2FoldChange** | log2 ratio of expression: test vs. reference |
| **lfcSE** | Standard error of the log2FC estimate |
| **stat** | Wald test statistic (log2FC / lfcSE) |
| **pvalue** | Raw p-value from the Wald test |
| **padj** | Adjusted p-value (Benjamini-Hochberg) |

### How to read the columns

- **padj < 0.05** → The gene is statistically significant (with 95% confidence after multiple testing correction).
- **log2FC > 1** → The gene is at least 2-fold up-regulated in the test condition.
- **log2FC < -1** → The gene is at least 2-fold down-regulated.
- **baseMean** → Higher values indicate more robustly expressed genes. Very low baseMean genes are unreliable.

### Post-Analysis Filtering

After the analysis, you can filter the results table with three modes:

1. **Default filters** (recommended):
   - baseMean ≥ 10
   - padj < 0.05
   - |log2FC| > 0.5

2. **Custom filters**: Set your own thresholds for each parameter.

3. **No filter**: Show all genes (useful for gene-specific lookups, but the table may be very large).

### Ranking Table (📈)

The ranking table shows genes sorted by the **Wald test statistic** (strongest signal first). This is useful for:
- Gene Set Enrichment Analysis (GSEA) — use this ranked list as input.
- Identifying the most confidently changed genes.

---

## Step 11 — Visualize Results

### Volcano Plot (🌋)

The most commonly used visualization for DE analysis.

- **X-axis:** log2 fold-change (effect size)
- **Y-axis:** -log10(adjusted p-value) (statistical significance)
- **Green dots:** Significantly up-regulated genes
- **Red dots:** Significantly down-regulated genes
- **Grey dots:** Not significant

**How to read it:**
- Genes in the **upper-left corner** are significantly down-regulated (large negative FC, high significance).
- Genes in the **upper-right corner** are significantly up-regulated.
- The **dashed lines** mark the significance thresholds.

Two volcano plots are available:
- **All genes:** Shows every gene in the dataset.
- **Filtered genes:** Shows only genes passing your post-analysis filters.

**Options:**
- 🏷️ **Label top genes:** Annotate the top 10 most significant genes with their names.
- 🔬 **baseMean filter:** Hide low-expression genes (reduces clutter).
- 📍 **Legend position:** Move the legend to avoid obscuring data points.

### PCA Plot (📊)

Principal Component Analysis shows how samples cluster in 2D.

- Each point is a sample, colored by condition.
- **Confidence ellipses** (2 standard deviations) show the spread of each group.
- **Good result:** Conditions form distinct, non-overlapping clusters.
- **Concerning result:** Conditions overlap heavily, suggesting weak differences or strong confounders.

### MA Plot (🎯)

Alternative view of DE results.

- **X-axis:** Mean expression (log10 scale)
- **Y-axis:** log2 fold-change
- Shows how fold-change relates to expression level.
- Well-expressed genes (right side) tend to have more reliable fold-changes.

### Heatmap (🔥)

Shows the top differentially expressed genes across all samples.

- **Rows:** Top 30 DE genes (by adjusted p-value)
- **Columns:** Samples (with condition color bar)
- **Colors:** Z-scored expression (red = high, blue = low)
- Hierarchical clustering reveals sample and gene groupings.

---

## Step 12 — Highlight Genes of Interest

After viewing the global results, you can highlight specific genes of interest on the volcano plot.

### How to use it

1. Scroll to the 🔎 **"Highlight genes"** section.
2. Enter gene names (comma-separated or one per line).
3. The app generates a special volcano plot where:
   - Your genes are shown as **orange diamonds** with labels.
   - All other genes are shown in grey.
   - A table shows the full DE statistics for your genes.

### Why this is useful

- **Validate candidate genes:** Check if your known genes of interest are DE.
- **Pathway genes:** Enter a list of genes from a specific pathway to see their collective behavior.
- **Publication figures:** The highlighted volcano plot is often used in publications.

> **Note:** If a gene is not found, the app will tell you which genes were not found and how many were successfully matched.

---

## Step 13 — Download Results

### Available Downloads

| Download | Content | Format |
|----------|---------|--------|
| **Filtered results** | Only genes passing your post-analysis filters | CSV or TSV |
| **All results** | Every gene with DE statistics | CSV or TSV |
| **Ranking** | Genes sorted by Wald statistic (for GSEA) | CSV or TSV |
| **Classified metadata** | Sample metadata with classification labels (if used) | CSV |

### Plot Downloads

Every plot has two download buttons:
- 📥 **PNG** (300 DPI) — for presentations and documents.
- 📥 **SVG** — for publication-quality vector graphics (editable in Illustrator/Inkscape).

---

## Step 14 — Audit Log (Reproducibility)

After the analysis completes, an expandable **Audit Log** section appears at the bottom of the results.

### What it contains

The audit log captures everything needed to reproduce your analysis:

| Section | Content |
|---------|---------|
| **Timestamp** | When the analysis was run (UTC) |
| **Platform** | Operating system, Python version, CPU architecture |
| **Library versions** | Exact versions of pydeseq2, pandas, numpy, scipy, matplotlib, etc. |
| **Input data** | Gene count, sample count, conditions, reference/test levels |
| **Parameters** | All DESeq2 parameters, gene filter settings, significance thresholds |
| **Filter statistics** | How many genes were removed by pre-filtering and why |
| **Results summary** | Total significant genes, up/down-regulated counts, top genes by padj |
| **Step timings** | How long each pipeline step took (in seconds) |

### Downloads

Two download buttons are available:

| Button | Format | Best for |
|--------|--------|----------|
| **📥 Download audit log (JSON)** | `.json` | Programmatic access, automated pipelines, archiving |
| **📥 Download audit log (TXT)** | `.txt` | Pasting into lab notebooks, methods sections, supplementary materials |

### Why this matters

Reproducibility is a cornerstone of scientific research. The audit log ensures that:
- You (or a reviewer) can verify exactly which parameters and software versions produced your results.
- Re-running the analysis with the same data + parameters should give identical results.
- The TXT format is ready to paste into a Methods section or lab notebook entry.

---

## Frequently Asked Questions

### Q: Can I use normalized data (FPKM, TPM)?

**No.** DESeq2 requires raw, unnormalized integer counts. It performs its own normalization internally (median-of-ratios). Using pre-normalized data will give incorrect results.

### Q: How many replicates do I need?

A minimum of **2 per condition**, but **3 or more** is strongly recommended. With only 2 replicates, the statistical power is very low and few genes will reach significance.

### Q: Why are some genes showing padj = NaN?

This can happen when:
- A gene has **zero counts in all samples of one condition** — the model cannot estimate a fold-change.
- A gene was flagged by **Cook's distance** as having an outlier sample — the p-value is set to NaN as a safety measure.
- A gene had an **extreme dispersion estimate** that the model could not stabilize.

### Q: My conditions have very different sample sizes (e.g., 400 Tumor vs. 60 Normal). Is that a problem?

DESeq2 can handle imbalanced designs, but extreme imbalances reduce statistical power. The app will show a ⚠️ class imbalance warning. Consider:
- Whether the imbalance reflects the real population.
- Whether you need to subsample the larger group for balanced analysis.

### Q: What if my samples are from different batches?

Use the **batch correction for PCA** option (Step 8) to visually assess batch effects. For the DE analysis itself, DESeq2 models the design formula without explicit batch terms in this app — if batch is a strong confounder, consider more advanced tools like DESeq2 in R with a multi-factor design formula.

### Q: How long will the analysis take?

See the time estimates in [Step 9](#step-9--run-the-analysis). The single most effective speedup is enabling **gene pre-filtering** (Step 7).

### Q: Can I use data from the Dataset Creator module?

**Yes!** The Dataset Creator produces a count matrix and metadata file in exactly the format this module expects. Just download them and upload here.

---

## Glossary

| Term | Definition |
|------|-----------|
| **Audit log** | A record of all parameters, library versions, data dimensions, and step timings from an analysis run — used for reproducibility |
| **baseMean** | Average normalized expression across all samples |
| **Benjamini-Hochberg (BH)** | Multiple testing correction method that controls the false discovery rate |
| **CPM** | Counts Per Million — a simple normalization for comparing across samples |
| **DESeq2** | Statistical method for differential expression analysis of count data |
| **Dispersion** | A measure of the variance of a gene's expression across samples, beyond what is expected by Poisson noise |
| **FDR** | False Discovery Rate — the expected proportion of false positives among all declared significant genes |
| **log2FC** | log2 Fold-Change — the log2 ratio of expression between two conditions |
| **MA plot** | A plot of log fold-change vs. mean expression |
| **padj** | Adjusted p-value (after multiple testing correction) |
| **PCA** | Principal Component Analysis — dimensionality reduction for visualizing sample relationships |
| **Raw counts** | Unnormalized integer read counts from RNA-seq alignment |
| **Volcano plot** | A plot of -log10(p-value) vs. log2 fold-change |
| **Wald test** | Statistical test used by DESeq2 to assess differential expression |

---

*This tutorial is part of the [BeginSeq Studio](https://github.com/JACKNINES/beginseq-studio) documentation.*
