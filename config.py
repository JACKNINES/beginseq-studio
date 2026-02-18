"""
config.py — Central configuration for the RNA-seq Analyzer application.

This module centralises ALL default parameters, thresholds and
configurations used throughout the application. Keeping them in a
single place avoids duplication of magic numbers scattered across the
code and makes modification easier without touching the logic.

Sections
--------
1. STREAMLIT_CONFIG : dict
   → Streamlit web page configuration.

2. DESEQ2_DEFAULTS : dict
   → Default parameters for the DESeq2 analysis.

3. VOLCANO_PLOT_CONFIG : dict
   → Visual configuration for the volcano plot.

4. PCA_PLOT_CONFIG : dict
   → Visual configuration for the PCA plot.

5. MA_PLOT_CONFIG : dict
   → Visual configuration for the MA plot.

6. HEATMAP_CONFIG : dict
   → Visual configuration for the top-genes heatmap.

7. FILE_CONFIG : dict
   → File read/write configuration.

Usage example
-------------
    from config import DESEQ2_DEFAULTS, VOLCANO_PLOT_CONFIG

    alpha = DESEQ2_DEFAULTS["alpha"]
    figsize = VOLCANO_PLOT_CONFIG["figsize"]
"""

# ──────────────────────────────────────────────────────────────────────
# 1. Streamlit configuration (web page)
# ──────────────────────────────────────────────────────────────────────
STREAMLIT_CONFIG: dict = {
    "page_title": "RNA-seq Analyzer",
    "layout": "centered",
    "page_icon": "🧬",
}

# ──────────────────────────────────────────────────────────────────────
# 2. Default parameters for DESeq2
# ──────────────────────────────────────────────────────────────────────
DESEQ2_DEFAULTS: dict = {
    # Name of the metadata column that defines the experimental
    # conditions (e.g. "treated" vs "control").
    "condition_col": "condition",

    # Reference level for the contrast.
    # Fold-changes are computed RELATIVE to this level.
    "reference_level": "control",

    # Adjusted p-value threshold (Benjamini-Hochberg).
    # Genes with padj < alpha are considered statistically significant.
    "alpha": 0.05,

    # Absolute log2 fold-change threshold.
    # Genes with |log2FC| > this value are considered biologically relevant.
    "log2fc_threshold": 1.0,

    # Minimum baseMean to include in volcano plot.
    # Genes with baseMean < this are low-expression noise and clutter the plot.
    "min_base_mean": 10,
}

# ──────────────────────────────────────────────────────────────────────
# 2b. Gene pre-filtering defaults (applied BEFORE DESeq2)
# ──────────────────────────────────────────────────────────────────────
# Filtering low-expressed genes is standard practice in RNA-seq:
#   - Reduces memory ~40-70% (removes noise genes).
#   - Improves statistical power (fewer multiple-testing corrections).
#   - Stabilises dispersion estimation (fewer near-zero genes).
#
# Default criterion: keep a gene only if it has ≥ min_total_count
# total counts AND is expressed (≥ min_count_per_sample count) in
# ≥ min_samples_expressing samples.
#
# IMPORTANT: For large TCGA datasets (~800+ samples), the old
# "expressed in >= smallest_group samples" criterion was far too weak
# because most genes (85%+) are detected in the smallest group.
# The auto mode now uses max(smallest_group, n_samples × auto_fraction)
# to ensure meaningful reduction (typically 60K → 20-25K genes).
GENE_FILTER_DEFAULTS: dict = {
    # Minimum sum of counts across ALL samples for a gene to be kept.
    "min_total_count": 10,

    # Minimum raw count to consider a gene "expressed" in a sample.
    # The gene must have ≥ this value in ≥ min_samples_expressing samples.
    # Default 1 = any non-zero count.  Higher values (e.g. 5, 10) give
    # more aggressive filtering for very large datasets.
    "min_count_per_sample": 1,

    # Minimum number of samples in which the gene must be "expressed".
    # Special value 0 means "auto": max(smallest_group, n × auto_fraction).
    "min_samples_expressing": 0,

    # Fraction of total samples used in auto mode for min_samples_expressing.
    # For 835 samples with 0.5 → min_samples = max(135, 417) = 417.
    # This ensures ~40-60% gene reduction on TCGA-scale datasets.
    "auto_samples_fraction": 0.5,

    # Whether to apply gene filtering by default.
    "enabled": True,
}

# ──────────────────────────────────────────────────────────────────────
# 3. Visual configuration for the volcano plot
# ──────────────────────────────────────────────────────────────────────
VOLCANO_PLOT_CONFIG: dict = {
    # Figure size in inches (width, height).
    "figsize": (10, 7),

    # ── Colours by category ──────────────────────────────────────────
    # Not significant: soft grey.
    "color_ns": "#D5D5D5",
    "alpha_ns": 0.5,
    "size_ns": 18,

    # Up-regulated: pastel green.
    "color_up": "#81C784",          # Material Green 300
    "alpha_up": 0.75,
    "size_up": 28,

    # Down-regulated: pastel red.
    "color_down": "#E57373",        # Material Red 300
    "alpha_down": 0.75,
    "size_down": 28,

    # ── Threshold lines ──────────────────────────────────────────────
    "threshold_linestyle": "--",
    "threshold_color": "#BDBDBD",   # light grey for minimalism
    "threshold_linewidth": 0.9,

    # ── Labels ───────────────────────────────────────────────────────
    "ylabel": r"$-\log_{10}$(adjusted p-value)",
    "xlabel_template": r"$\log_{{2}}$ Fold Change ({test} vs {ref})",
    "title": "Volcano Plot — Differential Expression",

    # ── Y axis: lower limit ────────────────────────────────────────
    "y_bottom": -1,

    # ── Fonts ────────────────────────────────────────────────────────
    "font_title": 14,
    "font_axes": 12,
    "font_legend": 9,
    "font_annotation": 9,

    # ── Density shading (gaussian KDE for NS points) ─────────────
    "density_shading": True,
    "density_cmap_ns": "Greys",

    # ── Highlight (genes of interest) ──────────────────────────────
    "highlight": {
        "color": "#FF6F00",         # intense orange for highlighted genes
        "edge_color": "#BF360C",    # dark border
        "size": 90,                 # point size
        "marker": "D",              # diamond
        "label_fontsize": 8.5,
        "label_color": "#BF360C",
        "bg_color": "#E0E0E0",      # background colour (non-highlighted genes)
        "bg_alpha": 0.3,
        "bg_size": 12,
    },
}

# ──────────────────────────────────────────────────────────────────────
# 4. Visual configuration for the PCA plot
# ──────────────────────────────────────────────────────────────────────
PCA_PLOT_CONFIG: dict = {
    "figsize": (9, 7),

    # Point size and transparency
    "point_size": 120,
    "point_alpha": 0.8,

    # Colour palette for conditions (up to 10)
    "color_palette": [
        "#4CAF50",  # green
        "#E53935",  # red
        "#1E88E5",  # blue
        "#FF9800",  # orange
        "#9C27B0",  # purple
        "#00BCD4",  # cyan
        "#795548",  # brown
        "#607D8B",  # blue-grey
        "#FFEB3B",  # yellow
        "#E91E63",  # pink
    ],

    # Sample labels
    "show_labels": True,
    "label_fontsize": 8,
    "label_offset": (5, 5),

    # Confidence ellipses per group
    "show_ellipses": True,
    "ellipse_alpha": 0.15,
    "ellipse_std": 2.0,  # standard deviations

    # Titles and axes
    "title": "PCA — Sample Clustering",
    "font_title": 14,
    "font_axes": 12,
    "font_legend": 9,
}

# ──────────────────────────────────────────────────────────────────────
# 5. Visual configuration for the MA plot
# ──────────────────────────────────────────────────────────────────────
MA_PLOT_CONFIG: dict = {
    "figsize": (10, 7),

    # Colours (reuses the volcano palette for consistency)
    "color_ns": "#D5D5D5",
    "alpha_ns": 0.5,
    "size_ns": 18,

    "color_up": "#81C784",
    "alpha_up": 0.75,
    "size_up": 28,

    "color_down": "#E57373",
    "alpha_down": 0.75,
    "size_down": 28,

    # Reference line at y=0
    "hline_color": "#555555",
    "hline_style": "--",
    "hline_width": 1.0,

    # Titles and axes
    "title": "MA Plot — Differential Expression",
    "xlabel": r"$\log_{10}$(mean expression)",
    "ylabel": r"$\log_{2}$ Fold Change",
    "font_title": 14,
    "font_axes": 12,
    "font_legend": 9,
    "font_annotation": 9,
}

# ──────────────────────────────────────────────────────────────────────
# 6. Visual configuration for the Heatmap
# ──────────────────────────────────────────────────────────────────────
HEATMAP_CONFIG: dict = {
    "figsize": (12, 10),

    # Number of top genes to display (hard-capped at 50 in code)
    "n_top_genes": 30,

    # Colormap (diverging, centred at zero for z-scores)
    "cmap": "RdBu_r",
    "center": 0,

    # Clustering
    "cluster_rows": True,
    "cluster_cols": True,
    "linkage_method": "ward",
    "distance_metric": "euclidean",

    # Dendrograms
    "show_row_dendrogram": True,
    "show_col_dendrogram": True,
    "dendrogram_ratio": 0.15,

    # Colour bar for conditions
    "show_condition_bar": True,
    "condition_bar_height": 0.03,

    # Titles and fonts
    "title": "Heatmap — Top Differentially Expressed Genes",
    "font_title": 14,
    "font_axes": 10,
    "font_colorbar": 9,

    # Gene labels
    "show_gene_labels": True,
    "gene_label_fontsize": 8,
}

# ──────────────────────────────────────────────────────────────────────
# 7. Memory optimisation
# ──────────────────────────────────────────────────────────────────────
MEMORY_CONFIG: dict = {
    # Counts are non-negative integers.  int32 supports up to ~2.1 billion
    # which covers any realistic RNA-seq count.  Saves 50% RAM vs int64.
    "counts_dtype": "int32",

    # Maximum samples before PCA disables sample labels (too cluttered).
    "pca_label_max_samples": 60,

    # Maximum samples before heatmap disables column clustering
    # (O(n²) distance matrix).
    "heatmap_cluster_cols_max_samples": 200,

    # PCA: keep only top-variance genes for dimensionality reduction.
    # 0 means "use all genes" (small datasets). Recommended: 2000-5000.
    "pca_top_var_genes": 2000,

    # NOTE: matrix assembly in gdc_client.build_count_matrix() now uses
    # numpy pre-allocation (zero-merge), so no chunk size is needed.

    # PyDESeq2: limit parallel workers to avoid memory duplication.
    # Each joblib worker gets a copy of the counts matrix.  With 8 cores
    # and a 163 MB matrix (855×25K int64), that's 8×163 = 1.3 GB just
    # for worker copies.  Limiting to 4 halves this to ~650 MB while
    # still being significantly faster than serial.
    "deseq2_n_cpus": 4,
}

# ──────────────────────────────────────────────────────────────────────
# 7b. Time estimation for DESeq2 pipeline
# ──────────────────────────────────────────────────────────────────────
# Empirical coefficients for estimating DESeq2 runtime.
#
# Model: T_step ≈ M × coeff_per_M × (1 + scale_factor_above_20M)
# where M = (n_samples × n_genes) / 1_000_000.
#
# For datasets > 20M elements (e.g. 800×25K), PyDESeq2 operations
# scale super-linearly due to:
#   - Memory pressure → swapping on 8 GB machines
#   - joblib worker copies grow with matrix size
#   - fit_size_factors median-of-ratios creates 2-3 float64 copies
#
# Calibrated on Apple M2 8 GB, PyDESeq2 v0.5.x, n_cpus=4.
# Measured with 200×15K, 855×25K, and 835×59K datasets.
TIME_ESTIMATION: dict = {
    # Per-step coefficients (seconds per million sample×gene elements)
    # e.g. 855 samples × 25K genes = 21.375M elements
    "size_factors_per_M": 0.4,        # median-of-ratios: heavy float64 allocs
    "genewise_disp_per_M": 2.5,       # heaviest: per-gene NB optimization
    "disp_trend_per_M": 0.08,         # parametric curve fit
    "map_disp_per_M": 1.8,            # second heaviest: per-gene MAP
    "fit_lfc_per_M": 1.2,             # per-gene IRLS
    "cooks_per_M": 0.15,              # distance calculation
    "wald_test_per_M": 0.4,           # Wald + BH correction
    "visualizations_per_M": 0.05,     # PCA, volcano, MA, heatmap

    # Above this threshold (in M elements), apply a super-linear
    # penalty to account for memory pressure and swapping.
    "superlinear_threshold_M": 20.0,
    # Exponent for the penalty: effective_M = M^exponent for M > threshold
    "superlinear_exponent": 1.3,

    # Fixed overhead (validation, gene filter, data prep) in seconds
    "fixed_overhead_s": 8.0,

    # Typical gene filtering reduction ratio (60K → ~25K)
    # Used to estimate post-filter dimensions when showing pre-run estimate.
    "typical_filter_ratio": 0.42,

    # Step names for per-step estimates (ordered)
    "step_keys": [
        "size_factors",
        "genewise_disp",
        "disp_trend",
        "map_disp",
        "fit_lfc",
        "cooks",
        "wald_test",
        "visualizations",
    ],
}

# ──────────────────────────────────────────────────────────────────────
# 8. scRNA-seq Analysis configuration
# ──────────────────────────────────────────────────────────────────────
SCRNA_CONFIG: dict = {
    # ── Quality Control ───────────────────────────────────────────────
    "qc_defaults": {
        "min_genes": 200,
        "max_genes": 5000,
        "min_counts": 500,
        "max_counts": 50000,
        "max_pct_mt": 20.0,
        "min_cells": 3,
    },

    # ── Normalization ─────────────────────────────────────────────────
    "normalization": {
        "target_sum": None,  # None = median of total counts
    },

    # ── Highly Variable Genes ────────────────────────────────────────
    "hvg": {
        "n_top_genes": 2000,
        "flavor": "seurat_v3",
        "batch_key": None,
    },

    # ── PCA ────────────────────────────────────────────────────────────
    "pca": {
        "n_comps": 50,
    },

    # ── Batch Effect Correction ──────────────────────────────────────
    "batch_correction": {
        "enabled": False,       # disabled by default
        "method": "harmony",    # only harmony for now
    },

    # ── Neighbors & Embedding ────────────────────────────────────────
    "neighbors": {
        "n_neighbors": 15,
        "n_pcs": 40,
    },
    "umap": {
        "min_dist": 0.5,
    },

    # ── Clustering ────────────────────────────────────────────────────
    "leiden": {
        "resolution": 0.5,
        "flavor": "igraph",
        "n_iterations": 2,
    },

    # ── Marker Genes ──────────────────────────────────────────────────
    "rank_genes": {
        "method": "wilcoxon",
        "n_genes": 25,
    },

    # ── Supported file formats (analysis accepts H5AD only) ──────────
    "file_extensions": ["h5ad"],

    # ── 10x integrator accepted formats ───────────────────────────────
    "integrator_matrix_ext": ["mtx", "gz"],
    "integrator_features_ext": ["tsv", "gz"],
    "integrator_barcodes_ext": ["tsv", "gz"],

    # ── Plot styling (consistent with bulk RNA-seq palette) ──────────
    "plot": {
        "figsize_qc": (12, 4),
        "figsize_pca": (9, 7),
        "figsize_umap": (9, 7),
        "figsize_heatmap": (12, 8),
        "figsize_dotplot": (12, 6),
        "figsize_violin": (12, 4),
        "color_palette": [
            "#4CAF50", "#E53935", "#1E88E5", "#FF9800", "#9C27B0",
            "#00BCD4", "#795548", "#607D8B", "#FFEB3B", "#E91E63",
            "#8BC34A", "#3F51B5", "#FF5722", "#009688", "#673AB7",
            "#CDDC39", "#2196F3", "#F44336", "#4DB6AC", "#BA68C8",
        ],
        "cmap_expression": "viridis",
        "cmap_heatmap": "RdBu_r",
        "point_size": 20,
        "point_alpha": 0.7,
        "dpi": 300,
        "font_title": 14,
        "font_axes": 12,
        "font_legend": 9,
    },
}

# ──────────────────────────────────────────────────────────────────────
# 9. File configuration (I/O)
# ──────────────────────────────────────────────────────────────────────
FILE_CONFIG: dict = {
    # Allowed file extensions for the counts matrix.
    "counts_extensions": ["tsv", "csv", "zip"],

    # Allowed file extensions for the metadata.
    "metadata_extensions": ["csv", "zip"],

    # Default filename for the downloadable results file.
    "default_output_filename": "deseq2_results.csv",

    # Separators automatically recognised by file extension.
    "separators": {
        "tsv": "\t",
        "csv": ",",
    },
}

# ──────────────────────────────────────────────────────────────────────
# 9b. Local-only upload size limits (MB)
# ──────────────────────────────────────────────────────────────────────
# When running locally the Streamlit server.maxUploadSize is raised to
# the maximum of these two values (5120 MB).  Each page then applies
# its own *soft* limit so the user gets a clear warning.
UPLOAD_LIMITS_LOCAL: dict = {
    "bulk_rna_mb": 2048,     # 2 GB for Bulk RNA-seq
    "scrna_mb": 5120,        # 5 GB for scRNA-seq
}

# ──────────────────────────────────────────────────────────────────────
# 10. Auto-shutdown configuration
# ──────────────────────────────────────────────────────────────────────
AUTO_SHUTDOWN_CONFIG: dict = {
    # Seconds to wait after the last browser tab closes before
    # shutting down the Streamlit server.  30 s is enough to survive
    # page refreshes and multipage navigation.
    "grace_period_seconds": 30,

    # How often (in seconds) to check for active sessions.
    "poll_interval_seconds": 5,

    # Set to False to keep the server running indefinitely.
    "enabled": True,
}
