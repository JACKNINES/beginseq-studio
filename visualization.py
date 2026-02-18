"""
visualization.py — Visualization generation for RNA-seq.

This module contains all the visualization functions for the
application. It currently implements the volcano plot, but is
designed to be extended with other common RNA-seq plots
(MA plot, heatmaps, PCA, etc.).

What is a Volcano Plot?
-----------------------
It is a scatter plot that shows:
- X-axis: log2 Fold Change (magnitude of expression change).
  → Positive values = UP-regulated genes in the test condition.
  → Negative values = DOWN-regulated genes in the test condition.
- Y-axis: -log10(adjusted p-value) (statistical significance).
  → Higher values = more statistically significant genes.

Visualized gene categories:
1. Not Significant  — light grey (#D5D5D5)
2. Up-regulated      — pastel green (#81C784)
3. Down-regulated    — pastel red (#E57373)
4. Highly Significant — orange (#FF9800)
   The top 5 significant genes with the lowest padj.

Functions
---------
prepare_volcano_data(results_df, alpha, log2fc_threshold)
    → Prepares and classifies genes into 4 categories.

create_volcano_plot(volcano_df, alpha, log2fc_threshold, test_level,
                    reference_level)
    → Generates the professional matplotlib volcano plot figure.

Usage example
-------------
    from visualization import prepare_volcano_data, create_volcano_plot

    volcano_df = prepare_volcano_data(results_df)
    fig = create_volcano_plot(volcano_df, test_level="treated",
                              reference_level="control")
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.figure

from config import (
    DESEQ2_DEFAULTS,
    VOLCANO_PLOT_CONFIG,
    PCA_PLOT_CONFIG,
    MA_PLOT_CONFIG,
    HEATMAP_CONFIG,
    MEMORY_CONFIG,
)
from i18n import t, tf


# ── Classification categories ────────────────────────────────────────
# Deliberate order: least prominent first (NS), then the significant ones.
CATEGORY_ORDER = ["NS", "Down", "Up"]

# ── Valid matplotlib legend locations ────────────────────────────────
LEGEND_LOCATIONS = [
    "best",
    "none",
    "upper right",
    "upper left",
    "lower left",
    "lower right",
    "right",
    "center left",
    "center right",
    "lower center",
    "upper center",
    "center",
]


def _apply_legend(
    ax,
    loc: str = "upper right",
    fontsize: int = 9,
    draggable: bool = True,
    **extra_kw,
) -> None:
    """
    Create and style a legend on *ax* with consistent scientific aesthetics.

    If *draggable* is True **and** the matplotlib backend supports
    interactive mode, the legend is made draggable.  On non-interactive
    backends (e.g. Agg, used by Streamlit) the call is silently skipped
    so no errors are raised.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    loc : str
        Any valid matplotlib legend *loc* string.
    fontsize : int
        Legend font size.
    draggable : bool
        Attempt to make the legend draggable (interactive backends only).
    **extra_kw
        Forwarded to ``ax.legend()``.
    """
    legend_kw = dict(
        loc=loc,
        fontsize=fontsize,
        frameon=True,
        framealpha=0.9,
        facecolor="white",
        edgecolor="#E0E0E0",
    )
    legend_kw.update(extra_kw)
    legend = ax.legend(**legend_kw)
    legend.set_zorder(10)

    if draggable:
        try:
            backend = matplotlib.get_backend().lower()
            if backend != "agg":
                legend.set_draggable(True)
        except Exception:
            pass  # Non-interactive backend — ignore silently


def prepare_volcano_data(
    results_df: pd.DataFrame,
    alpha: float = DESEQ2_DEFAULTS["alpha"],
    log2fc_threshold: float = DESEQ2_DEFAULTS["log2fc_threshold"],
    min_base_mean: float = 0,
) -> pd.DataFrame:
    """
    Prepare the DESeq2 results DataFrame for volcano plot
    visualization, classifying each gene into one of four categories.

    Categories assigned in the "category" column:
    ──────────────────────────────────────────────
    1. "NS"                 → Not significant (does not pass base thresholds).
    2. "Up"                 → Significant and up-regulated (log2FC > 0).
    3. "Down"               → Significant and down-regulated (log2FC < 0).
    4. "Highly significant" → Top N significant genes with the lowest
                              padj (default N=5, configurable in
                              VOLCANO_PLOT_CONFIG["highly_sig_top_n"]).

    Parameters
    ----------
    results_df : pd.DataFrame
        DESeq2 results table (output from compute_contrast).
        Must contain columns: "padj", "log2FoldChange".
    alpha : float
        Adjusted p-value threshold for base significance.
    log2fc_threshold : float
        |log2FC| threshold for base biological relevance.
    min_base_mean : float
        Minimum baseMean threshold.  Genes with baseMean below this
        value are excluded (removes low-expression noise).

    Returns
    -------
    pd.DataFrame
        DataFrame with additional columns:
        - "neg_log10_padj" : float — Y-axis value.
        - "significant"    : bool  — True if it passes base thresholds.
        - "category"       : str   — one of the 4 categories.
    """
    cfg = VOLCANO_PLOT_CONFIG

    volcano_df = results_df.dropna(subset=["padj"]).copy()

    # ── Filter low baseMean genes (noise removal) ─────────────────
    if min_base_mean > 0 and "baseMean" in volcano_df.columns:
        volcano_df = volcano_df[volcano_df["baseMean"] >= min_base_mean]

    # ── Y-axis ─────────────────────────────────────────────────────────
    volcano_df["neg_log10_padj"] = -np.log10(volcano_df["padj"])

    # ── Base significance (padj < alpha AND |log2FC| > threshold) ───
    volcano_df["significant"] = (
        (volcano_df["padj"] < alpha)
        & (volcano_df["log2FoldChange"].abs() > log2fc_threshold)
    )

    # ── Assign category ────────────────────────────────────────────────
    volcano_df["category"] = "NS"

    up_mask = volcano_df["significant"] & (volcano_df["log2FoldChange"] > 0)
    down_mask = volcano_df["significant"] & (volcano_df["log2FoldChange"] < 0)

    volcano_df.loc[up_mask, "category"] = "Up"
    volcano_df.loc[down_mask, "category"] = "Down"

    return volcano_df


def create_volcano_plot(
    volcano_df: pd.DataFrame,
    alpha: float = DESEQ2_DEFAULTS["alpha"],
    log2fc_threshold: float = DESEQ2_DEFAULTS["log2fc_threshold"],
    test_level: str = "treated",
    reference_level: str = DESEQ2_DEFAULTS["reference_level"],
    label_genes: list = None,
    title: str = None,
    shrinkage_applied: bool = False,
    legend_loc: str = "upper right",
) -> matplotlib.figure.Figure:
    """
    Generate a professional volcano plot with minimalist style.

    Visual features:
    ────────────────
    • 3 categories with distinct colors (NS, Up, Down).
    • Y-axis starts at -1 to provide lower visual margin.
    • Clean theme: white background, no top/right borders.
    • Threshold lines in light grey so they don't compete with the data.
    • Organized legend with gene counts per category.
    • Summary annotation in the upper left corner.
    • Optionally labels specific genes.

    Parameters
    ----------
    volcano_df : pd.DataFrame
        DataFrame prepared with prepare_volcano_data().
        Must contain: "log2FoldChange", "neg_log10_padj", "category".
    alpha : float
        Base significance threshold (for horizontal line).
    log2fc_threshold : float
        Base fold-change threshold (for vertical lines).
    test_level : str
        Name of the test condition (X-axis label).
    reference_level : str
        Name of the reference condition (X-axis label).
    label_genes : list, optional
        List of gene names to label on the plot.
    title : str, optional
        Custom title for the plot. If None, uses the config default.

    Returns
    -------
    matplotlib.figure.Figure
        Matplotlib Figure object ready for st.pyplot().
    """
    cfg = VOLCANO_PLOT_CONFIG

    # ── Minimalist style ────────────────────────────────────────────
    with plt.style.context("seaborn-v0_8-whitegrid"):
        fig, ax = plt.subplots(figsize=cfg["figsize"])

    # Clean background
    ax.set_facecolor("#FAFAFA")
    fig.patch.set_facecolor("white")

    # Remove top and right borders (spines)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#CCCCCC")
    ax.spines["bottom"].set_color("#CCCCCC")

    # Soft grid on Y-axis only
    ax.yaxis.grid(True, linestyle=":", linewidth=0.5, color="#E0E0E0")
    ax.xaxis.grid(False)
    ax.set_axisbelow(True)

    # Tick styling
    ax.tick_params(
        axis="both", which="major",
        labelsize=10, colors="#555555",
        direction="out", length=4, width=0.8,
    )

    # ── Visual properties map per category ──────────────────────────
    cat_style = {
        "NS": {
            "color": cfg["color_ns"],
            "alpha": cfg["alpha_ns"],
            "size": cfg["size_ns"],
            "label": t("plot.not_significant"),
            "zorder": 1,
            "marker": "o",
        },
        "Down": {
            "color": cfg["color_down"],
            "alpha": cfg["alpha_down"],
            "size": cfg["size_down"],
            "label": t("plot.down_regulated"),
            "zorder": 2,
            "marker": "v",       # downward triangle
        },
        "Up": {
            "color": cfg["color_up"],
            "alpha": cfg["alpha_up"],
            "size": cfg["size_up"],
            "label": t("plot.up_regulated"),
            "zorder": 2,
            "marker": "^",       # upward triangle
        },
    }

    # ── Draw each category in order (NS first, Highly last) ───────────
    counts_by_cat = {}

    for cat in CATEGORY_ORDER:
        mask = volcano_df["category"] == cat
        n = mask.sum()
        counts_by_cat[cat] = n

        if n == 0:
            continue

        style = cat_style[cat]

        # ── Density shading for NS points ─────────────────────────
        if cat == "NS" and cfg.get("density_shading", False) and n > 30:
            try:
                from scipy.stats import gaussian_kde

                x_ns = volcano_df.loc[mask, "log2FoldChange"].values
                y_ns = volcano_df.loc[mask, "neg_log10_padj"].values
                xy = np.vstack([x_ns, y_ns])
                kde = gaussian_kde(xy)
                density = kde(xy)
                # Normalise density to [0.3, 0.85] so the lightest
                # points are a visible light-grey and the densest are
                # dark-grey (not black).
                d_min, d_max = density.min(), density.max()
                if d_max > d_min:
                    density_norm = (density - d_min) / (d_max - d_min)
                else:
                    density_norm = np.ones_like(density)
                # Map to [0.3, 0.85] range of the Greys colormap
                density_mapped = 0.3 + density_norm * 0.55

                sc = ax.scatter(
                    x_ns, y_ns,
                    c=density_mapped,
                    cmap=cfg.get("density_cmap_ns", "Greys"),
                    alpha=0.65,
                    s=style["size"],
                    marker=style["marker"],
                    edgecolors="none",
                    linewidths=0,
                    zorder=style["zorder"],
                    vmin=0, vmax=1,
                )
                # Manual legend entry (scatter with cmap can't auto-legend)
                ax.scatter([], [], c=style["color"], alpha=style["alpha"],
                           s=style["size"], marker=style["marker"],
                           label=f'{style["label"]} ({n:,})')
                continue
            except Exception:
                pass  # Fall back to normal scatter below

        ax.scatter(
            volcano_df.loc[mask, "log2FoldChange"],
            volcano_df.loc[mask, "neg_log10_padj"],
            c=style["color"],
            alpha=style["alpha"],
            s=style["size"],
            marker=style["marker"],
            label=f'{style["label"]} ({n:,})',
            edgecolors="white",
            linewidths=0.3,
            zorder=style["zorder"],
        )

    # ── Threshold lines (subtle, don't compete with data) ──────────
    threshold_kw = dict(
        linestyle=cfg["threshold_linestyle"],
        color=cfg["threshold_color"],
        linewidth=cfg["threshold_linewidth"],
        zorder=0,
    )

    # Horizontal line: -log10(alpha)
    ax.axhline(-np.log10(alpha), **threshold_kw)

    # Vertical lines: ±log2fc_threshold
    ax.axvline(-log2fc_threshold, **threshold_kw)
    ax.axvline(log2fc_threshold, **threshold_kw)

    # ── Threshold line labels ────────────────────────────────────────
    y_top = ax.get_ylim()[1]
    ax.text(
        log2fc_threshold, y_top * 0.02,
        f"  |log₂FC| = {log2fc_threshold}",
        fontsize=7.5, color="#999999", va="bottom",
    )
    ax.text(
        ax.get_xlim()[1] * 0.55, -np.log10(alpha) + y_top * 0.01,
        f"padj = {alpha}",
        fontsize=7.5, color="#999999", va="bottom",
    )

    # ── Y-axis: lower limit at -1 ───────────────────────────────────
    y_bottom = cfg.get("y_bottom", -1)
    ax.set_ylim(bottom=y_bottom)

    # ── Axis labels and title ────────────────────────────────────────
    xlabel = cfg["xlabel_template"].format(test=test_level, ref=reference_level)
    if shrinkage_applied:
        xlabel += t("plot.volcano_shrunk_suffix")
    ax.set_xlabel(xlabel, fontsize=cfg["font_axes"], color="#333333")
    ax.set_ylabel(cfg["ylabel"], fontsize=cfg["font_axes"], color="#333333")

    plot_title = title if title is not None else cfg["title"]
    ax.set_title(
        plot_title,
        fontsize=cfg["font_title"],
        fontweight="bold",
        color="#222222",
        pad=15,
    )

    # ── Organized legend ────────────────────────────────────────────
    _apply_legend(
        ax, loc=legend_loc, fontsize=cfg["font_legend"],
        borderpad=0.8, labelspacing=0.6, handletextpad=0.5,
    )

    # ── Summary annotation ───────────────────────────────────────────
    n_total = len(volcano_df)
    n_up = counts_by_cat.get("Up", 0)
    n_down = counts_by_cat.get("Down", 0)
    n_sig_total = n_up + n_down

    summary_text = tf(
        "plot.volcano_summary",
        n_total=n_total, n_sig=n_sig_total, n_up=n_up, n_down=n_down,
    )
    ax.annotate(
        summary_text,
        xy=(0.02, 0.97),
        xycoords="axes fraction",
        verticalalignment="top",
        fontsize=cfg["font_annotation"],
        fontstyle="italic",
        color="#555555",
        bbox=dict(
            boxstyle="round,pad=0.4",
            facecolor="white",
            edgecolor="#E0E0E0",
            alpha=0.9,
        ),
    )

    # ── Specific gene labels ────────────────────────────────────────
    if label_genes:
        from adjustText import adjust_text

        texts = []
        genes_to_label = [g for g in label_genes if g in volcano_df.index]

        for gene_name in genes_to_label:
            if gene_name in volcano_df.index:
                row = volcano_df.loc[gene_name]
                txt = ax.text(
                    row["log2FoldChange"],
                    row["neg_log10_padj"],
                    f"  {gene_name}",
                    fontsize=7.5,
                    color="#E65100",
                    fontweight="bold",
                    ha="left",
                    va="center",
                    zorder=10,
                )
                texts.append(txt)

        # Adjust label positions to avoid overlap
        if len(texts) > 1:
            try:
                adjust_text(
                    texts,
                    ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#E65100", lw=0.5, alpha=0.6),
                    expand_points=(1.5, 1.5),
                    force_text=(0.8, 1.0),
                )
            except Exception:
                pass

    # ── Wide margins ──────────────────────────────────────────────────
    plt.tight_layout(pad=1.5)

    return fig


def create_volcano_highlight(
    volcano_df: pd.DataFrame,
    highlight_genes: list,
    alpha: float = DESEQ2_DEFAULTS["alpha"],
    log2fc_threshold: float = DESEQ2_DEFAULTS["log2fc_threshold"],
    test_level: str = "treated",
    reference_level: str = DESEQ2_DEFAULTS["reference_level"],
    legend_loc: str = "upper right",
) -> matplotlib.figure.Figure:
    """
    Generate a volcano plot where the specified genes are visually
    highlighted and the rest of the points are dimmed in the background.

    Parameters
    ----------
    volcano_df : pd.DataFrame
        DataFrame prepared with prepare_volcano_data().
    highlight_genes : list
        List of gene names/IDs to highlight.
    alpha, log2fc_threshold, test_level, reference_level :
        Same parameters as create_volcano_plot.

    Returns
    -------
    matplotlib.figure.Figure
    """
    from adjustText import adjust_text

    cfg = VOLCANO_PLOT_CONFIG
    hl_cfg = cfg.get("highlight", {})

    # ── Separate highlighted genes vs background ─────────────────────
    highlight_set = set(highlight_genes)
    found = [g for g in highlight_genes if g in volcano_df.index]
    mask_hl = volcano_df.index.isin(highlight_set)
    df_bg = volcano_df[~mask_hl]
    df_hl = volcano_df[mask_hl]

    # ── Figure ────────────────────────────────────────────────────────
    with plt.style.context("seaborn-v0_8-whitegrid"):
        fig, ax = plt.subplots(figsize=cfg["figsize"])

    ax.set_facecolor("#FAFAFA")
    fig.patch.set_facecolor("white")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#CCCCCC")
    ax.spines["bottom"].set_color("#CCCCCC")
    ax.yaxis.grid(True, linestyle=":", linewidth=0.5, color="#E0E0E0")
    ax.xaxis.grid(False)
    ax.set_axisbelow(True)
    ax.tick_params(
        axis="both", which="major",
        labelsize=10, colors="#555555",
        direction="out", length=4, width=0.8,
    )

    # ── Background: all genes dimmed ─────────────────────────────────
    ax.scatter(
        df_bg["log2FoldChange"],
        df_bg["neg_log10_padj"],
        c=hl_cfg.get("bg_color", "#E0E0E0"),
        alpha=hl_cfg.get("bg_alpha", 0.3),
        s=hl_cfg.get("bg_size", 12),
        marker="o",
        edgecolors="none",
        zorder=1,
        label=tf("plot.highlight_other_genes", count=len(df_bg)),
    )

    # ── Highlight: genes of interest ──────────────────────────────────
    hl_color = hl_cfg.get("color", "#FF6F00")
    hl_edge = hl_cfg.get("edge_color", "#BF360C")

    if len(df_hl) > 0:
        ax.scatter(
            df_hl["log2FoldChange"],
            df_hl["neg_log10_padj"],
            c=hl_color,
            alpha=0.95,
            s=hl_cfg.get("size", 90),
            marker=hl_cfg.get("marker", "D"),
            edgecolors=hl_edge,
            linewidths=1.2,
            zorder=5,
            label=tf("plot.highlight_highlighted", count=len(df_hl)),
        )

        # ── Labels for highlighted genes ──────────────────────────────
        texts = []
        for gene_name in found:
            row = volcano_df.loc[gene_name]
            txt = ax.text(
                row["log2FoldChange"],
                row["neg_log10_padj"],
                f"  {gene_name}",
                fontsize=hl_cfg.get("label_fontsize", 8.5),
                color=hl_cfg.get("label_color", "#BF360C"),
                fontweight="bold",
                ha="left",
                va="center",
                zorder=10,
            )
            texts.append(txt)

        if len(texts) > 1:
            try:
                adjust_text(
                    texts,
                    ax=ax,
                    arrowprops=dict(
                        arrowstyle="-",
                        color=hl_cfg.get("label_color", "#BF360C"),
                        lw=0.6,
                        alpha=0.7,
                    ),
                    expand_points=(1.5, 1.5),
                    force_text=(0.8, 1.0),
                )
            except Exception:
                pass

    # ── Threshold lines ─────────────────────────────────────────────
    threshold_kw = dict(
        linestyle=cfg["threshold_linestyle"],
        color=cfg["threshold_color"],
        linewidth=cfg["threshold_linewidth"],
        zorder=0,
    )
    ax.axhline(-np.log10(alpha), **threshold_kw)
    ax.axvline(-log2fc_threshold, **threshold_kw)
    ax.axvline(log2fc_threshold, **threshold_kw)

    # ── Threshold labels ──────────────────────────────────────────────
    y_top = ax.get_ylim()[1]
    ax.text(
        log2fc_threshold, y_top * 0.02,
        f"  |log₂FC| = {log2fc_threshold}",
        fontsize=7.5, color="#999999", va="bottom",
    )
    ax.text(
        ax.get_xlim()[1] * 0.55, -np.log10(alpha) + y_top * 0.01,
        f"padj = {alpha}",
        fontsize=7.5, color="#999999", va="bottom",
    )

    # ── Y-axis ────────────────────────────────────────────────────────
    y_bottom = cfg.get("y_bottom", -1)
    ax.set_ylim(bottom=y_bottom)

    # ── Axes and title ────────────────────────────────────────────────
    xlabel = cfg["xlabel_template"].format(test=test_level, ref=reference_level)
    ax.set_xlabel(xlabel, fontsize=cfg["font_axes"], color="#333333")
    ax.set_ylabel(cfg["ylabel"], fontsize=cfg["font_axes"], color="#333333")
    ax.set_title(
        t("plot.highlight_title"),
        fontsize=cfg["font_title"],
        fontweight="bold",
        color="#222222",
        pad=15,
    )

    # ── Legend ────────────────────────────────────────────────────────
    _apply_legend(ax, loc=legend_loc, fontsize=cfg["font_legend"])

    # ── Info box with genes not found ─────────────────────────────────
    not_found = [g for g in highlight_genes if g not in volcano_df.index]
    info_lines = [tf("plot.highlight_info", n_found=len(found))]
    if not_found:
        info_lines.append(tf("plot.highlight_not_found", n_not_found=len(not_found)))
    info_text = "\n".join(info_lines)

    ax.annotate(
        info_text,
        xy=(0.02, 0.97),
        xycoords="axes fraction",
        verticalalignment="top",
        fontsize=cfg["font_annotation"],
        fontstyle="italic",
        color="#555555",
        bbox=dict(
            boxstyle="round,pad=0.4",
            facecolor="white",
            edgecolor="#E0E0E0",
            alpha=0.9,
        ),
    )

    plt.tight_layout(pad=1.5)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# PCA PLOT
# ═══════════════════════════════════════════════════════════════════════

def prepare_pca_data(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    condition_col: str = "condition",
    dds=None,
    batch_col: str | None = None,
) -> pd.DataFrame:
    """
    Compute PCA of samples using a variance-stabilising transformation.

    **Why not raw log2(counts + 1)?**
    The naïve log2 transform is biased: genes with low mean counts have
    artificially inflated variance, which dominates PCA axes.  DESeq2's
    VST (via ``dds.layers["vst_counts"]``) corrects for the mean-variance
    relationship.  If the fitted ``dds`` object is provided we use its
    VST counts; otherwise we fall back to log2(counts + 1) as before.

    **Memory optimisation:**
    PCA on ~60K genes × 800 samples (48M floats) is wasteful because
    most genes contribute only noise.  We select the ``top_n`` genes
    by variance *before* running PCA, reducing dimensions to ~2K.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Raw count matrix (genes × samples).
    metadata_df : pd.DataFrame
        Metadata with sample index and condition column.
    condition_col : str
        Condition column name in metadata.
    dds : DeseqDataSet or None
        If provided and it has ``vst_counts`` in its layers, those
        are used instead of the naïve log transform.
    batch_col : str or None
        If provided, apply ComBat batch correction on the transformed
        matrix BEFORE PCA.  The column must exist in metadata_df.

    Returns
    -------
    pd.DataFrame
        Columns: PC1, PC2, condition.  Index = sample names.
        attrs["var_explained"] = array of explained variance ratios.
        attrs["transform"] = "vst" | "log2(n+1)" (for display).
    """
    from sklearn.decomposition import PCA

    top_n = MEMORY_CONFIG.get("pca_top_var_genes", 2000)
    transform_name = "log2(n+1)"

    # ── Choose transformation ─────────────────────────────────────
    if dds is not None and hasattr(dds, "layers") and "vst_counts" in dds.layers:
        # VST counts from pydeseq2: already (samples × genes) in AnnData
        vst = pd.DataFrame(
            dds.layers["vst_counts"],
            index=dds.obs_names,
            columns=dds.var_names,
        )
        # Transpose to (genes × samples) for variance selection
        transformed = vst.T
        transform_name = "VST"
        del vst
    else:
        # Fallback: log2(counts + 1) — reasonable approximation
        transformed = np.log2(counts_df.astype("float32") + 1)

    # ── Batch correction (ComBat) ──────────────────────────────────
    if batch_col and batch_col in metadata_df.columns:
        try:
            from combat.pycombat import pycombat

            # pycombat expects (genes × samples) DataFrame and batch list
            batch_labels = metadata_df.loc[
                transformed.columns, batch_col
            ].tolist()
            transformed = pycombat(transformed, batch_labels)
            transform_name += " + ComBat"
        except Exception:
            pass  # If ComBat fails, proceed without correction

    # ── Variance-based gene selection ─────────────────────────────
    n_genes = transformed.shape[0]
    if top_n > 0 and n_genes > top_n:
        gene_var = transformed.var(axis=1)
        top_genes = gene_var.nlargest(top_n).index
        transformed = transformed.loc[top_genes]

    # ── PCA (samples × selected_genes) ────────────────────────────
    X = transformed.T.values  # (n_samples, n_genes)
    sample_names = transformed.columns.tolist()

    pca = PCA(n_components=2)
    pca_coords = pca.fit_transform(X)

    # Free the large matrix before building the small result
    del X, transformed

    pca_df = pd.DataFrame(
        {"PC1": pca_coords[:, 0], "PC2": pca_coords[:, 1]},
        index=sample_names,
    )
    pca_df.index.name = "sample"

    # Add condition from metadata
    pca_df[condition_col] = metadata_df.loc[pca_df.index, condition_col].values

    # Store variance explained + transform used
    pca_df.attrs["var_explained"] = pca.explained_variance_ratio_
    pca_df.attrs["transform"] = transform_name

    return pca_df


def create_pca_plot(
    pca_df: pd.DataFrame,
    condition_col: str = "condition",
    legend_loc: str = "best",
) -> matplotlib.figure.Figure:
    """
    Generate a professional PCA plot with colors per condition.

    Parameters
    ----------
    pca_df : pd.DataFrame
        DataFrame with PC1, PC2, and condition column.
        Must have the attrs["var_explained"] attribute.
    condition_col : str
        Name of the condition column.

    Returns
    -------
    matplotlib.figure.Figure
    """
    from matplotlib.patches import Ellipse

    cfg = PCA_PLOT_CONFIG
    var_explained = pca_df.attrs.get("var_explained", [0, 0])

    # ── Dynamic size based on number of samples ─────────────────────
    n_samples = len(pca_df)
    base_width, base_height = cfg["figsize"]

    # Expand if there are many samples
    if n_samples > 20:
        width = min(base_width + (n_samples - 20) * 0.15, 14)
        height = min(base_height + (n_samples - 20) * 0.1, 10)
    else:
        width, height = base_width, base_height

    with plt.style.context("seaborn-v0_8-whitegrid"):
        fig, ax = plt.subplots(figsize=(width, height))

    ax.set_facecolor("#FAFAFA")
    fig.patch.set_facecolor("white")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#CCCCCC")
    ax.spines["bottom"].set_color("#CCCCCC")

    # Color map per condition
    conditions = pca_df[condition_col].unique().tolist()
    color_map = {
        cond: cfg["color_palette"][i % len(cfg["color_palette"])]
        for i, cond in enumerate(conditions)
    }

    # Draw confidence ellipses if enabled
    if cfg.get("show_ellipses", True):
        for cond in conditions:
            subset = pca_df[pca_df[condition_col] == cond]
            if len(subset) < 3:
                continue
            mean_x, mean_y = subset["PC1"].mean(), subset["PC2"].mean()
            std_x, std_y = subset["PC1"].std(), subset["PC2"].std()
            ellipse = Ellipse(
                (mean_x, mean_y),
                width=std_x * cfg["ellipse_std"] * 2,
                height=std_y * cfg["ellipse_std"] * 2,
                facecolor=color_map[cond],
                alpha=cfg["ellipse_alpha"],
                edgecolor=color_map[cond],
                linewidth=1,
                zorder=1,
            )
            ax.add_patch(ellipse)

    # Draw points
    for cond in conditions:
        subset = pca_df[pca_df[condition_col] == cond]
        ax.scatter(
            subset["PC1"],
            subset["PC2"],
            c=color_map[cond],
            s=cfg["point_size"],
            alpha=cfg["point_alpha"],
            label=f"{cond} (n={len(subset)})",
            edgecolors="white",
            linewidths=0.8,
            zorder=2,
        )

    # ── Sample labels with automatic position adjustment ─────────────
    max_label_samples = MEMORY_CONFIG.get("pca_label_max_samples", 60)
    show_labels = cfg.get("show_labels", True) and n_samples <= max_label_samples

    if show_labels:
        # Adjust font size based on number of samples
        label_fontsize = cfg.get("label_fontsize", 8)
        if n_samples > 30:
            label_fontsize = max(6, label_fontsize - 2)
        elif n_samples > 50:
            label_fontsize = max(5, label_fontsize - 3)

        texts = []
        for sample_name, row in pca_df.iterrows():
            txt = ax.text(
                row["PC1"],
                row["PC2"],
                f"  {sample_name}",
                fontsize=label_fontsize,
                color="#555555",
                zorder=3,
                ha="left",
                va="center",
            )
            texts.append(txt)

    # Axes
    ax.set_xlabel(
        tf("plot.pca_axis_variance", n=1, pct=var_explained[0]*100),
        fontsize=cfg["font_axes"],
    )
    ax.set_ylabel(
        tf("plot.pca_axis_variance", n=2, pct=var_explained[1]*100),
        fontsize=cfg["font_axes"],
    )
    transform_label = pca_df.attrs.get("transform", "")
    pca_title = cfg["title"]
    if transform_label:
        pca_title += f"  [{transform_label}]"
    ax.set_title(pca_title, fontsize=cfg["font_title"], fontweight="bold", pad=15)

    # Legend
    _apply_legend(ax, loc=legend_loc, fontsize=cfg["font_legend"])

    # Wide margins for labels
    plt.tight_layout(pad=1.5)

    # Slightly expand limits so labels are not clipped
    x_margin = (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.08
    y_margin = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.08
    ax.set_xlim(ax.get_xlim()[0] - x_margin, ax.get_xlim()[1] + x_margin)
    ax.set_ylim(ax.get_ylim()[0] - y_margin, ax.get_ylim()[1] + y_margin)

    return fig


# ═══════════════════════════════════════════════════════════════════════
# MA PLOT
# ═══════════════════════════════════════════════════════════════════════

def prepare_ma_data(
    results_df: pd.DataFrame,
    alpha: float = DESEQ2_DEFAULTS["alpha"],
    log2fc_threshold: float = DESEQ2_DEFAULTS["log2fc_threshold"],
) -> pd.DataFrame:
    """
    Prepare data for MA plot: A (average) vs M (log2FC).

    Parameters
    ----------
    results_df : pd.DataFrame
        DESeq2 results with baseMean, log2FoldChange, padj.
    alpha : float
        Significance threshold.
    log2fc_threshold : float
        |log2FC| threshold.

    Returns
    -------
    pd.DataFrame
        DataFrame with A, M, category (same 4 categories as volcano).
    """
    cfg = MA_PLOT_CONFIG

    ma_df = results_df.dropna(subset=["padj", "baseMean"]).copy()

    # A = log10(baseMean), M = log2FoldChange
    ma_df["A"] = np.log10(ma_df["baseMean"] + 1)
    ma_df["M"] = ma_df["log2FoldChange"]

    # Classification same as volcano
    ma_df["significant"] = (
        (ma_df["padj"] < alpha)
        & (ma_df["log2FoldChange"].abs() > log2fc_threshold)
    )

    ma_df["category"] = "NS"
    up_mask = ma_df["significant"] & (ma_df["log2FoldChange"] > 0)
    down_mask = ma_df["significant"] & (ma_df["log2FoldChange"] < 0)
    ma_df.loc[up_mask, "category"] = "Up"
    ma_df.loc[down_mask, "category"] = "Down"

    return ma_df


def create_ma_plot(
    ma_df: pd.DataFrame,
    test_level: str = "treated",
    reference_level: str = "control",
    legend_loc: str = "upper right",
) -> matplotlib.figure.Figure:
    """
    Generate a professional MA plot.

    Parameters
    ----------
    ma_df : pd.DataFrame
        DataFrame prepared with prepare_ma_data().
    test_level, reference_level : str
        For the Y-axis label.

    Returns
    -------
    matplotlib.figure.Figure
    """
    cfg = MA_PLOT_CONFIG

    with plt.style.context("seaborn-v0_8-whitegrid"):
        fig, ax = plt.subplots(figsize=cfg["figsize"])

    ax.set_facecolor("#FAFAFA")
    fig.patch.set_facecolor("white")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#CCCCCC")
    ax.spines["bottom"].set_color("#CCCCCC")

    ax.yaxis.grid(True, linestyle=":", linewidth=0.5, color="#E0E0E0")
    ax.xaxis.grid(False)
    ax.set_axisbelow(True)

    # Styles per category (same as volcano)
    cat_style = {
        "NS": {"color": cfg["color_ns"], "alpha": cfg["alpha_ns"],
               "size": cfg["size_ns"], "marker": "o", "label": t("plot.ma_not_significant")},
        "Down": {"color": cfg["color_down"], "alpha": cfg["alpha_down"],
                 "size": cfg["size_down"], "marker": "v", "label": t("plot.ma_down_regulated")},
        "Up": {"color": cfg["color_up"], "alpha": cfg["alpha_up"],
               "size": cfg["size_up"], "marker": "^", "label": t("plot.ma_up_regulated")},
    }

    counts_by_cat = {}
    for cat in CATEGORY_ORDER:
        mask = ma_df["category"] == cat
        n = mask.sum()
        counts_by_cat[cat] = n
        if n == 0:
            continue
        style = cat_style[cat]
        ax.scatter(
            ma_df.loc[mask, "A"],
            ma_df.loc[mask, "M"],
            c=style["color"],
            alpha=style["alpha"],
            s=style["size"],
            marker=style["marker"],
            label=f'{style["label"]} ({n:,})',
            edgecolors="white",
            linewidths=0.3,
            zorder=1 if cat == "NS" else 2,
        )

    # Horizontal line at M=0
    ax.axhline(
        0,
        color=cfg["hline_color"],
        linestyle=cfg["hline_style"],
        linewidth=cfg["hline_width"],
        zorder=0,
    )

    ax.set_xlabel(cfg["xlabel"], fontsize=cfg["font_axes"], color="#333333")
    ax.set_ylabel(
        f'{cfg["ylabel"]} ({test_level} vs {reference_level})',
        fontsize=cfg["font_axes"],
        color="#333333",
    )
    ax.set_title(cfg["title"], fontsize=cfg["font_title"], fontweight="bold", pad=15)

    _apply_legend(ax, loc=legend_loc, fontsize=cfg["font_legend"])

    # Summary
    n_total = len(ma_df)
    n_up = counts_by_cat.get("Up", 0)
    n_down = counts_by_cat.get("Down", 0)
    summary_text = tf("plot.ma_summary", n_total=n_total, n_sig=n_up + n_down, n_up=n_up, n_down=n_down)
    ax.annotate(
        summary_text,
        xy=(0.02, 0.97),
        xycoords="axes fraction",
        verticalalignment="top",
        fontsize=cfg["font_annotation"],
        fontstyle="italic",
        color="#555555",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                  edgecolor="#E0E0E0", alpha=0.9),
    )

    plt.tight_layout(pad=1.5)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# HEATMAP
# ═══════════════════════════════════════════════════════════════════════

def prepare_heatmap_data(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    results_df: pd.DataFrame,
    condition_col: str = "condition",
    n_top: int | None = None,
) -> pd.DataFrame:
    """
    Prepare data for a heatmap: top N genes by padj, z-score normalised.

    **Memory considerations:**
    Only the top ``n_top`` genes are kept (default 30).  Z-score
    normalisation is computed with vectorised numpy operations instead
    of ``DataFrame.apply()`` (which creates Python objects per row).

    Parameters
    ----------
    counts_df : pd.DataFrame
        Count matrix (genes × samples).  Only the top genes are read.
    metadata_df : pd.DataFrame
        Metadata with condition column.
    results_df : pd.DataFrame
        DESeq2 results — used to select the top genes.
    condition_col : str
        Condition column name.
    n_top : int or None
        Number of top genes.  Defaults to config value (30).

    Returns
    -------
    pd.DataFrame
        Z-score matrix (genes × samples) for the top genes.
        attrs["conditions"]: list of condition per sample.
        attrs["condition_col"]: name of condition column.
    """
    cfg = HEATMAP_CONFIG
    if n_top is None:
        n_top = cfg.get("n_top_genes", 30)

    # Cap to avoid insane memory use
    n_top = min(n_top, 50)

    # Select top genes by padj
    top_genes = (
        results_df
        .dropna(subset=["padj"])
        .nsmallest(n_top, "padj")
        .index
        .tolist()
    )

    # Only available genes (might be filtered out)
    available = [g for g in top_genes if g in counts_df.index]
    if not available:
        # Fallback: take top by variance
        gene_var = counts_df.var(axis=1)
        available = gene_var.nlargest(min(n_top, len(gene_var))).index.tolist()

    heatmap_counts = counts_df.loc[available]

    # log2(counts + 1) — use float32 to save memory
    log_counts = np.log2(heatmap_counts.values.astype("float32") + 1)

    # Z-score per gene (row) — vectorised
    row_mean = log_counts.mean(axis=1, keepdims=True)
    row_std = log_counts.std(axis=1, keepdims=True)
    row_std[row_std == 0] = 1.0  # avoid division by zero
    zscores = (log_counts - row_mean) / row_std

    zscore_df = pd.DataFrame(
        zscores,
        index=available,
        columns=heatmap_counts.columns,
    )

    # Store conditions as attribute
    conditions = metadata_df.loc[zscore_df.columns, condition_col].tolist()
    zscore_df.attrs["conditions"] = conditions
    zscore_df.attrs["condition_col"] = condition_col

    return zscore_df


def create_heatmap(
    heatmap_df: pd.DataFrame,
) -> matplotlib.figure.Figure:
    """
    Generate a heatmap with hierarchical clustering using seaborn.clustermap.

    Parameters
    ----------
    heatmap_df : pd.DataFrame
        Z-score matrix (genes x samples).
        Must have attrs["conditions"].

    Returns
    -------
    matplotlib.figure.Figure
    """
    import seaborn as sns

    cfg = HEATMAP_CONFIG
    conditions = heatmap_df.attrs.get("conditions", [])

    n_genes = len(heatmap_df.index)
    n_samples = len(heatmap_df.columns)

    # ── Dynamic figure size ─────────────────────────────────────────
    width = max(10, min(18, 4 + n_samples * 0.5))
    height = max(8, min(16, 3 + n_genes * 0.4))

    # ── Colors per condition ────────────────────────────────────────
    col_colors = None
    lut = {}
    if cfg.get("show_condition_bar", True) and conditions:
        unique_conds = list(dict.fromkeys(conditions))
        palette = PCA_PLOT_CONFIG["color_palette"][:len(unique_conds)]
        lut = {c: palette[i] for i, c in enumerate(unique_conds)}
        col_colors = pd.Series(conditions, index=heatmap_df.columns).map(lut)

    # ── Adaptive font size ─────────────────────────────────────────
    gene_fontsize = 10 if n_genes <= 20 else (8 if n_genes <= 35 else 6)
    sample_fontsize = 10 if n_samples <= 15 else (8 if n_samples <= 25 else 6)

    # ── Adaptive clustering ─────────────────────────────────────────
    # Column clustering computes an O(n²) distance matrix.  For 800
    # samples that's ~2.4 GB of float64.  Disable it automatically
    # above a threshold (default 200 samples).
    max_cluster_cols = MEMORY_CONFIG.get("heatmap_cluster_cols_max_samples", 200)
    do_cluster_cols = cfg.get("cluster_cols", True) and n_samples <= max_cluster_cols

    # Also suppress x-axis tick labels for very large sample counts
    show_xticks = n_samples <= 80

    # ── Create clustermap ─────────────────────────────────────────────
    vmax = np.abs(heatmap_df.values).max()

    g = sns.clustermap(
        heatmap_df,
        cmap=cfg.get("cmap", "RdBu_r"),
        vmin=-vmax,
        vmax=vmax,
        center=0,
        row_cluster=cfg.get("cluster_rows", True),
        col_cluster=do_cluster_cols,
        col_colors=col_colors,
        figsize=(width, height),
        dendrogram_ratio=(0.15, 0.15),
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        tree_kws={"linewidths": 1.0, "colors": "#666666"},
        linewidths=0.5 if n_samples <= 100 else 0,
        linecolor="white",
        xticklabels=show_xticks,
        yticklabels=True,
    )

    # ── Gene labels (Y-axis, right side) ─────────────────────────────
    g.ax_heatmap.yaxis.tick_right()
    g.ax_heatmap.yaxis.set_label_position("right")
    g.ax_heatmap.set_yticklabels(
        g.ax_heatmap.get_yticklabels(),
        fontsize=gene_fontsize,
        rotation=0,
    )

    # ── Sample labels (X-axis, bottom) ──────────────────────────────
    if show_xticks:
        rotation = 45 if n_samples <= 15 else 90
        ha = "right" if rotation == 45 else "center"
        g.ax_heatmap.set_xticklabels(
            g.ax_heatmap.get_xticklabels(),
            fontsize=sample_fontsize,
            rotation=rotation,
            ha=ha,
        )

    # ── Title ────────────────────────────────────────────────────────
    g.fig.suptitle(
        cfg.get("title", "Top Differentially Expressed Genes"),
        fontsize=cfg.get("font_title", 14),
        fontweight="bold",
        y=1.02,
    )

    # ── Colorbar ─────────────────────────────────────────────────────
    g.cax.set_ylabel(t("plot.heatmap_zscore_label"), fontsize=10, rotation=90, labelpad=8)
    g.cax.yaxis.set_label_position("left")
    g.cax.tick_params(labelsize=8)

    # ── Condition legend ────────────────────────────────────────────
    if lut:
        from matplotlib.patches import Patch
        legend_handles = [Patch(facecolor=c, edgecolor="gray", label=l)
                         for l, c in lut.items()]
        g.ax_heatmap.legend(
            handles=legend_handles,
            title=t("plot.heatmap_condition_title"),
            loc="upper left",
            bbox_to_anchor=(1.12, 1.0),
            fontsize=9,
            title_fontsize=10,
            frameon=True,
            framealpha=0.95,
            edgecolor="#CCCCCC",
        )

    # ── Adjust margins ─────────────────────────────────────────────
    g.fig.tight_layout(rect=[0.05, 0.08, 0.92, 0.95])

    return g.fig
