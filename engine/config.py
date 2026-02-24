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
# 0. Cohesive visual theme — THREE presets: Dark · Light · Cyberpunk
# ──────────────────────────────────────────────────────────────────────
# A unified visual identity shared by ALL plots (bulk + scRNA-seq).
# Every visualisation module reads from the mutable ``THEME`` dict so
# changing it at runtime (via ``set_theme()``) propagates everywhere.
#
# Available presets: "dark", "light", "cyberpunk"

THEME_DARK: dict = {
    # ── Surfaces ──────────────────────────────────────────────────────
    "bg":          "#0e1117",
    "surface":     "#161b22",
    "surface2":    "#1c2333",
    "border":      "#2a3444",
    # ── Typography ────────────────────────────────────────────────────
    "text":        "#e6edf3",
    "text_muted":  "#8b949e",
    "text_subtle": "#6e7681",
    "font_family": "monospace",
    # ── Accents ───────────────────────────────────────────────────────
    "cyan":   "#58d5c1", "green": "#3fb950", "red":   "#f85149",
    "amber":  "#d29922", "blue":  "#58a6ff", "violet":"#bc8cff",
    "pink":   "#f778ba",
    # ── Semantic ──────────────────────────────────────────────────────
    "color_up": "#3fb950", "color_down": "#f85149", "color_ns": "#484f58",
    # ── 20-colour palette ─────────────────────────────────────────────
    "palette": [
        "#58d5c1","#f85149","#58a6ff","#d29922","#bc8cff",
        "#f778ba","#3fb950","#db6d28","#79c0ff","#d2a8ff",
        "#56d364","#ff7b72","#a5d6ff","#e3b341","#7ee787",
        "#ffa657","#b392f0","#f0883e","#39d353","#da3633",
    ],
    # ── Grid / Points ─────────────────────────────────────────────────
    "grid_color": "#21262d", "grid_width": 0.6, "grid_style": ":",
    "point_edge": "#0e1117", "point_edge_width": 0.4,
    # ── Colormaps ─────────────────────────────────────────────────────
    "cmap_diverging": "RdBu_r", "cmap_sequential": "viridis",
    "cmap_density":   "Greys",
    # ── Font sizes ────────────────────────────────────────────────────
    "font_title": 14, "font_axes": 11, "font_legend": 9,
    "font_annotation": 9, "font_ticks": 9,
    # ── Legend / Annotation ───────────────────────────────────────────
    "legend_bg": "#161b22", "legend_edge": "#2a3444", "legend_alpha": 0.92,
    "annot_bg":  "#1c2333", "annot_edge":  "#2a3444", "annot_alpha":  0.92,
}

THEME_LIGHT: dict = {
    # ── Surfaces ──────────────────────────────────────────────────────
    "bg":          "#ffffff",
    "surface":     "#f8f9fa",
    "surface2":    "#f0f2f5",
    "border":      "#d0d7de",
    # ── Typography ────────────────────────────────────────────────────
    "text":        "#1f2328",
    "text_muted":  "#656d76",
    "text_subtle": "#8b949e",
    "font_family": "sans-serif",
    # ── Accents (Material-inspired, high-contrast on white) ───────────
    "cyan":   "#00897b", "green": "#2e7d32", "red":   "#c62828",
    "amber":  "#f57f17", "blue":  "#1565c0", "violet":"#6a1b9a",
    "pink":   "#ad1457",
    # ── Semantic ──────────────────────────────────────────────────────
    "color_up": "#2e7d32", "color_down": "#c62828", "color_ns": "#bdbdbd",
    # ── 20-colour palette (Material Design 700-series) ────────────────
    "palette": [
        "#00897b","#c62828","#1565c0","#f57f17","#6a1b9a",
        "#ad1457","#2e7d32","#e65100","#0277bd","#8e24aa",
        "#00695c","#d32f2f","#1976d2","#f9a825","#7b1fa2",
        "#c2185b","#388e3c","#ef6c00","#0288d1","#9c27b0",
    ],
    # ── Grid / Points ─────────────────────────────────────────────────
    "grid_color": "#e1e4e8", "grid_width": 0.6, "grid_style": ":",
    "point_edge": "#ffffff", "point_edge_width": 0.4,
    # ── Colormaps ─────────────────────────────────────────────────────
    "cmap_diverging": "RdBu_r", "cmap_sequential": "viridis",
    "cmap_density":   "Greys",
    # ── Font sizes ────────────────────────────────────────────────────
    "font_title": 14, "font_axes": 11, "font_legend": 9,
    "font_annotation": 9, "font_ticks": 9,
    # ── Legend / Annotation ───────────────────────────────────────────
    "legend_bg": "#ffffff", "legend_edge": "#d0d7de", "legend_alpha": 0.95,
    "annot_bg":  "#f0f2f5", "annot_edge":  "#d0d7de", "annot_alpha":  0.95,
}

THEME_CYBERPUNK: dict = {
    # ── Surfaces (deeper blacks) ──────────────────────────────────────
    "bg":          "#0a0a12",
    "surface":     "#0d1420",
    "surface2":    "#111a2e",
    "border":      "#1a2744",
    # ── Typography ────────────────────────────────────────────────────
    "text":        "#e8edf5",
    "text_muted":  "#7a8ba8",
    "text_subtle": "#4a5d7a",
    "font_family": "monospace",
    # ── Accents (neon-glow, high-saturation) ──────────────────────────
    "cyan":   "#00f0ff", "green": "#00ff88", "red":   "#ff2e5f",
    "amber":  "#ffb300", "blue":  "#448aff", "violet":"#d500f9",
    "pink":   "#ff4081",
    # ── Semantic ──────────────────────────────────────────────────────
    "color_up": "#00ff88", "color_down": "#ff2e5f", "color_ns": "#2a3450",
    # ── 20-colour palette (neon-glow) ─────────────────────────────────
    "palette": [
        "#00f0ff","#ff2e5f","#448aff","#ffb300","#d500f9",
        "#ff4081","#00ff88","#ff6e40","#40c4ff","#ea80fc",
        "#69f0ae","#ff5252","#82b1ff","#ffd740","#b388ff",
        "#ff80ab","#00e676","#ff9100","#18ffff","#e040fb",
    ],
    # ── Grid / Points ─────────────────────────────────────────────────
    "grid_color": "#152038", "grid_width": 0.5, "grid_style": ":",
    "point_edge": "#0a0a12", "point_edge_width": 0.3,
    # ── Colormaps ─────────────────────────────────────────────────────
    "cmap_diverging": "RdBu_r", "cmap_sequential": "plasma",
    "cmap_density":   "Greys",
    # ── Font sizes ────────────────────────────────────────────────────
    "font_title": 14, "font_axes": 11, "font_legend": 9,
    "font_annotation": 9, "font_ticks": 9,
    # ── Legend / Annotation ───────────────────────────────────────────
    "legend_bg": "#0d1420", "legend_edge": "#1a2744", "legend_alpha": 0.94,
    "annot_bg":  "#111a2e", "annot_edge":  "#1a2744", "annot_alpha":  0.94,
}

# ── Preset registry ──────────────────────────────────────────────────
THEME_PRESETS: dict = {
    "dark":      THEME_DARK,
    "light":     THEME_LIGHT,
    "cyberpunk": THEME_CYBERPUNK,
}

# ── Active theme (mutable — swapped in place by set_theme()) ─────────
THEME: dict = dict(THEME_DARK)   # default is Dark


def set_theme(name: str) -> None:
    """Switch the active theme in place.

    Parameters
    ----------
    name : str
        One of ``"dark"``, ``"light"``, ``"cyberpunk"``.

    After calling this, all subsequent reads of ``THEME[key]`` will
    return values from the chosen preset.  Plot configs that were
    computed at import time (``VOLCANO_PLOT_CONFIG``, etc.) are also
    refreshed so they stay in sync.
    """
    preset = THEME_PRESETS.get(name)
    if preset is None:
        raise ValueError(f"Unknown theme '{name}'. Choose from: {list(THEME_PRESETS)}")
    THEME.clear()
    THEME.update(preset)
    # Refresh import-time plot configs that cached THEME values
    _refresh_plot_configs()


def _refresh_plot_configs() -> None:
    """Re-sync every plot config dict with the current THEME values.

    Called automatically by ``set_theme()``.  The dicts are updated
    **in place** so that any module that already imported them via
    ``from engine.config import VOLCANO_PLOT_CONFIG`` sees the new
    colours immediately.
    """
    # ── Volcano ───────────────────────────────────────────────────────
    VOLCANO_PLOT_CONFIG.update({
        "color_ns":         THEME["color_ns"],
        "color_up":         THEME["color_up"],
        "color_down":       THEME["color_down"],
        "threshold_color":  THEME["border"],
        "font_title":       THEME["font_title"],
        "font_axes":        THEME["font_axes"],
        "font_legend":      THEME["font_legend"],
        "font_annotation":  THEME["font_annotation"],
        "density_cmap_ns":  THEME["cmap_density"],
    })
    VOLCANO_PLOT_CONFIG["highlight"].update({
        "color":       THEME["amber"],
        "label_color": THEME["amber"],
        "bg_color":    THEME["color_ns"],
    })
    # ── PCA ───────────────────────────────────────────────────────────
    PCA_PLOT_CONFIG.update({
        "color_palette": THEME["palette"],
        "font_title":    THEME["font_title"],
        "font_axes":     THEME["font_axes"],
        "font_legend":   THEME["font_legend"],
    })
    # ── MA ────────────────────────────────────────────────────────────
    MA_PLOT_CONFIG.update({
        "color_ns":        THEME["color_ns"],
        "color_up":        THEME["color_up"],
        "color_down":      THEME["color_down"],
        "hline_color":     THEME["text_subtle"],
        "font_title":      THEME["font_title"],
        "font_axes":       THEME["font_axes"],
        "font_legend":     THEME["font_legend"],
        "font_annotation": THEME["font_annotation"],
    })
    # ── Heatmap ───────────────────────────────────────────────────────
    HEATMAP_CONFIG.update({
        "cmap":       THEME["cmap_diverging"],
        "font_title": THEME["font_title"],
    })
    # ── scRNA plot sub-dict ───────────────────────────────────────────
    SCRNA_CONFIG["plot"].update({
        "color_palette":    THEME["palette"],
        "cmap_expression":  THEME["cmap_sequential"],
        "cmap_heatmap":     THEME["cmap_diverging"],
        "font_title":       THEME["font_title"],
        "font_axes":        THEME["font_axes"],
        "font_legend":      THEME["font_legend"],
    })


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
    "figsize": (10, 7),

    # ── Colours by category (from THEME) ─────────────────────────────
    "color_ns":   THEME["color_ns"],
    "alpha_ns":   0.6,
    "size_ns":    18,

    "color_up":   THEME["color_up"],
    "alpha_up":   0.85,
    "size_up":    30,

    "color_down": THEME["color_down"],
    "alpha_down": 0.85,
    "size_down":  30,

    # ── Threshold lines ──────────────────────────────────────────────
    "threshold_linestyle": "--",
    "threshold_color": THEME["border"],
    "threshold_linewidth": 0.9,

    # ── Labels ───────────────────────────────────────────────────────
    "ylabel": r"$-\log_{10}$(adjusted p-value)",
    "xlabel_template": r"$\log_{{2}}$ Fold Change ({test} vs {ref})",
    "title": "Volcano Plot — Differential Expression",

    # ── Y axis: lower limit ────────────────────────────────────────
    "y_bottom": -1,

    # ── Fonts ────────────────────────────────────────────────────────
    "font_title":      THEME["font_title"],
    "font_axes":       THEME["font_axes"],
    "font_legend":     THEME["font_legend"],
    "font_annotation": THEME["font_annotation"],

    # ── Density shading (gaussian KDE for NS points) ─────────────
    "density_shading": True,
    "density_cmap_ns": THEME["cmap_density"],

    # ── Highlight (genes of interest) ──────────────────────────────
    "highlight": {
        "color":       THEME["amber"],
        "edge_color":  "#a67615",
        "size":        90,
        "marker":      "D",
        "label_fontsize": 8.5,
        "label_color": THEME["amber"],
        "bg_color":    THEME["color_ns"],
        "bg_alpha":    0.3,
        "bg_size":     12,
    },
}

# ──────────────────────────────────────────────────────────────────────
# 4. Visual configuration for the PCA plot
# ──────────────────────────────────────────────────────────────────────
PCA_PLOT_CONFIG: dict = {
    "figsize": (9, 7),

    "point_size":  120,
    "point_alpha": 0.85,

    # Uses the shared THEME palette
    "color_palette": THEME["palette"],

    # Sample labels
    "show_labels":    True,
    "label_fontsize": 8,
    "label_offset":   (5, 5),

    # Confidence ellipses per group
    "show_ellipses":  True,
    "ellipse_alpha":  0.12,
    "ellipse_std":    2.0,

    # Titles and axes
    "title":       "PCA — Sample Clustering",
    "font_title":  THEME["font_title"],
    "font_axes":   THEME["font_axes"],
    "font_legend": THEME["font_legend"],
}

# ──────────────────────────────────────────────────────────────────────
# 5. Visual configuration for the MA plot
# ──────────────────────────────────────────────────────────────────────
MA_PLOT_CONFIG: dict = {
    "figsize": (10, 7),

    "color_ns":   THEME["color_ns"],
    "alpha_ns":   0.6,
    "size_ns":    18,

    "color_up":   THEME["color_up"],
    "alpha_up":   0.85,
    "size_up":    30,

    "color_down": THEME["color_down"],
    "alpha_down": 0.85,
    "size_down":  30,

    # Reference line at y=0
    "hline_color": THEME["text_subtle"],
    "hline_style": "--",
    "hline_width": 1.0,

    # Titles and axes
    "title":           "MA Plot — Differential Expression",
    "xlabel":          r"$\log_{10}$(mean expression)",
    "ylabel":          r"$\log_{2}$ Fold Change",
    "font_title":      THEME["font_title"],
    "font_axes":       THEME["font_axes"],
    "font_legend":     THEME["font_legend"],
    "font_annotation": THEME["font_annotation"],
}

# ──────────────────────────────────────────────────────────────────────
# 6. Visual configuration for the Heatmap
# ──────────────────────────────────────────────────────────────────────
HEATMAP_CONFIG: dict = {
    "figsize": (12, 10),

    "n_top_genes": 30,

    "cmap":   THEME["cmap_diverging"],
    "center": 0,

    # Clustering
    "cluster_rows":    True,
    "cluster_cols":    True,
    "linkage_method":  "ward",
    "distance_metric": "euclidean",

    # Dendrograms
    "show_row_dendrogram": True,
    "show_col_dendrogram": True,
    "dendrogram_ratio":    0.15,

    # Colour bar for conditions
    "show_condition_bar":    True,
    "condition_bar_height":  0.03,

    # Titles and fonts
    "title":              "Heatmap — Top Differentially Expressed Genes",
    "font_title":         THEME["font_title"],
    "font_axes":          10,
    "font_colorbar":      9,

    # Gene labels
    "show_gene_labels":    True,
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

    # ── Plot styling (uses shared THEME palette) ─────────────────────
    "plot": {
        "figsize_qc": (12, 4),
        "figsize_pca": (9, 7),
        "figsize_umap": (9, 7),
        "figsize_heatmap": (12, 8),
        "figsize_dotplot": (12, 6),
        "figsize_violin": (12, 4),
        "color_palette": THEME["palette"],
        "cmap_expression": THEME["cmap_sequential"],
        "cmap_heatmap":    THEME["cmap_diverging"],
        "point_size":  20,
        "point_alpha": 0.8,
        "dpi": 300,
        "font_title":  THEME["font_title"],
        "font_axes":   THEME["font_axes"],
        "font_legend": THEME["font_legend"],
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
