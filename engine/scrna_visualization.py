"""
scrna_visualization.py -- Visualization functions for scRNA-seq analysis.

Generates matplotlib figures for each step of the scRNA-seq pipeline:
    - QC violin plots
    - QC scatter plots
    - HVG plot
    - PCA variance ratio (elbow plot)
    - PCA embedding
    - UMAP embedding (colored by cluster, gene, QC metric)
    - Marker gene dotplot
    - Marker gene heatmap
    - Marker gene ranking bar plot

All functions return a matplotlib Figure object for Streamlit display.
"""

import io
from typing import Optional

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from engine.config import SCRNA_CONFIG, THEME, set_theme, THEME_PRESETS
from engine.i18n_labels import get_label, get_label_formatted

# ── Pure i18n layer (no Streamlit dependency) ───────────────────────
_current_lang = "en"


def set_language(lang: str):
    """Set the active language for plot labels."""
    global _current_lang
    _current_lang = lang


def t(key: str) -> str:
    """Return translated label for *key* in the active language."""
    return get_label(key, _current_lang)


def tf(key: str, **kwargs) -> str:
    """Return translated formatted label for *key* in the active language."""
    return get_label_formatted(key, _current_lang, **kwargs)

matplotlib.use("Agg")

_PLOT_CFG = SCRNA_CONFIG["plot"]


# ══════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════

def _cluster_colors(n: int) -> list:
    """Return n colors from the configured palette, cycling if needed."""
    palette = _PLOT_CFG["color_palette"]
    return [palette[i % len(palette)] for i in range(n)]


def _fig_to_bytes(fig: plt.Figure, fmt: str = "png") -> bytes:
    """Serialize a figure to bytes for download."""
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=_PLOT_CFG["dpi"], bbox_inches="tight",
                facecolor=fig.get_facecolor())
    buf.seek(0)
    return buf.getvalue()


def _apply_theme(fig, ax):
    """Apply the cohesive dark theme to a figure and axes."""
    plt.rcParams["font.family"] = THEME["font_family"]
    fig.patch.set_facecolor(THEME["bg"])
    ax.set_facecolor(THEME["surface"])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(THEME["border"])
    ax.spines["bottom"].set_color(THEME["border"])
    ax.yaxis.grid(True, linestyle=THEME["grid_style"],
                  linewidth=THEME["grid_width"], color=THEME["grid_color"])
    ax.xaxis.grid(False)
    ax.set_axisbelow(True)
    ax.tick_params(
        axis="both", which="major",
        labelsize=THEME["font_ticks"], colors=THEME["text_muted"],
        direction="out", length=4, width=0.8,
    )


def _themed_fig(figsize):
    """Create a figure + axes with the dark theme pre-applied."""
    fig, ax = plt.subplots(figsize=figsize)
    _apply_theme(fig, ax)
    return fig, ax


def _apply_legend(ax, loc: str = "best", fontsize: int = 9, **extra_kw) -> None:
    """Create and style a legend with consistent aesthetics.

    Uses bbox_to_anchor to place the legend outside the axes (right side)
    so it never overlaps data points.  Pass loc='none' to suppress the legend.
    """
    if loc == "none":
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()
        return

    legend_kw = dict(
        loc=loc,
        fontsize=fontsize,
        frameon=True,
        framealpha=THEME["legend_alpha"],
        facecolor=THEME["legend_bg"],
        edgecolor=THEME["legend_edge"],
        labelcolor=THEME["text"],
        bbox_to_anchor=(1.02, 1),
        borderaxespad=0,
    )
    legend_kw.update(extra_kw)
    legend = ax.legend(**legend_kw)
    legend.set_zorder(10)


# ══════════════════════════════════════════════════════════════════════
# QC Visualizations
# ══════════════════════════════════════════════════════════════════════

def plot_qc_violin(
    adata: ad.AnnData,
    keys: list[str] | None = None,
) -> plt.Figure:
    """Violin plots of QC metrics (gene counts, total counts, % MT, etc.).

    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix with QC metrics in ``adata.obs``.
    keys : list[str] | None
        QC columns to plot.  Defaults to all available among
        ``n_genes_by_counts``, ``total_counts``, ``pct_counts_mt``,
        ``pct_counts_ribo``, ``pct_counts_hb``.

    Returns
    -------
    plt.Figure
    """
    if keys is None:
        keys = []
        for k in ["n_genes_by_counts", "total_counts", "pct_counts_mt",
                   "pct_counts_ribo", "pct_counts_hb"]:
            if k in adata.obs.columns:
                keys.append(k)

    n = len(keys)
    fig, axes = plt.subplots(1, n, figsize=(_PLOT_CFG["figsize_qc"][0], _PLOT_CFG["figsize_qc"][1]))
    if n == 1:
        axes = [axes]

    fig.patch.set_facecolor(THEME["bg"])
    plt.rcParams["font.family"] = THEME["font_family"]

    for ax, key in zip(axes, keys):
        ax.set_facecolor(THEME["surface"])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color(THEME["border"])
        ax.spines["bottom"].set_color(THEME["border"])
        ax.tick_params(colors=THEME["text_muted"], labelsize=THEME["font_ticks"])

        data = adata.obs[key].dropna().values
        parts = ax.violinplot(data, showmeans=True, showmedians=True)
        for pc in parts["bodies"]:
            pc.set_facecolor(THEME["blue"])
            pc.set_alpha(0.7)
        # Style the lines (means, medians, bars)
        for partname in ("cmeans", "cmedians", "cbars", "cmins", "cmaxes"):
            if partname in parts:
                parts[partname].set_color(THEME["text_muted"])
        ax.set_title(key, fontsize=10, color=THEME["text"])
        ax.set_ylabel("")
        ax.set_xticks([])

    fig.suptitle(t("plot.qc_title"), fontsize=13, fontweight="bold", color=THEME["text"])
    fig.tight_layout()
    return fig


def plot_qc_scatter(
    adata: ad.AnnData,
    x: str = "total_counts",
    y: str = "n_genes_by_counts",
    color: str = "pct_counts_mt",
) -> plt.Figure:
    """Scatter plot of two QC metrics, colored by a third.

    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix with QC metrics in ``adata.obs``.
    x, y : str
        Column names in ``adata.obs`` for the X and Y axes.
    color : str
        Column in ``adata.obs`` used for the colour gradient.
        Falls back to a uniform colour if the column is absent.

    Returns
    -------
    plt.Figure
    """
    fig, ax = _themed_fig((8, 6))

    if color in adata.obs.columns:
        sc = ax.scatter(
            adata.obs[x], adata.obs[y],
            c=adata.obs[color],
            cmap="YlOrRd",
            s=_PLOT_CFG["point_size"],
            alpha=_PLOT_CFG["point_alpha"],
            edgecolors="none",
        )
        cbar = plt.colorbar(sc, ax=ax, label=color, shrink=0.8)
        cbar.ax.yaxis.set_tick_params(color=THEME["text_muted"])
        cbar.ax.yaxis.label.set_color(THEME["text_muted"])
        plt.setp(cbar.ax.yaxis.get_ticklabels(), color=THEME["text_muted"])
    else:
        ax.scatter(
            adata.obs[x], adata.obs[y],
            s=_PLOT_CFG["point_size"],
            alpha=_PLOT_CFG["point_alpha"],
            color=THEME["blue"],
            edgecolors="none",
        )

    ax.set_xlabel(x, fontsize=_PLOT_CFG["font_axes"], color=THEME["text_muted"])
    ax.set_ylabel(y, fontsize=_PLOT_CFG["font_axes"], color=THEME["text_muted"])
    ax.set_title(t("plot.qc_scatter_title"), fontsize=13, fontweight="bold", color=THEME["text"])
    fig.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════
# HVG Plot
# ══════════════════════════════════════════════════════════════════════

def plot_hvg(adata: ad.AnnData) -> plt.Figure:
    """Plot highly variable genes: mean vs dispersion.

    Uses ``adata.raw.var`` when available (after subsetting to HVGs in
    ``scale_data``), so the plot correctly shows HVG vs non-HVG genes.
    """
    fig, ax = _themed_fig((8, 5))

    # Prefer the full gene annotation stored in adata.raw (before HVG subset)
    if adata.raw is not None and "highly_variable" in adata.raw.var.columns:
        var_df = adata.raw.var
    else:
        var_df = adata.var

    hvg_mask = var_df["highly_variable"].values.astype(bool)

    if "dispersions_norm" in var_df.columns and "means" in var_df.columns:
        # cell_ranger / seurat flavors → scatter of mean vs dispersion
        ax.scatter(
            var_df["means"][~hvg_mask],
            var_df["dispersions_norm"][~hvg_mask],
            s=4, alpha=0.3, color=THEME["color_ns"], label=t("plot.hvg_other_label"),
        )
        ax.scatter(
            var_df["means"][hvg_mask],
            var_df["dispersions_norm"][hvg_mask],
            s=6, alpha=0.7, color=THEME["red"], label=t("plot.hvg_label"),
        )
        ax.set_xlabel(t("plot.hvg_xlabel"), fontsize=_PLOT_CFG["font_axes"], color=THEME["text_muted"])
        ax.set_ylabel(t("plot.hvg_ylabel_disp"), fontsize=_PLOT_CFG["font_axes"], color=THEME["text_muted"])

    elif "variances_norm" in var_df.columns and "means" in var_df.columns:
        # seurat_v3 flavor → scatter of mean vs normalized variance
        ax.scatter(
            var_df["means"][~hvg_mask],
            var_df["variances_norm"][~hvg_mask],
            s=4, alpha=0.3, color=THEME["color_ns"], label=t("plot.hvg_other_label"),
        )
        ax.scatter(
            var_df["means"][hvg_mask],
            var_df["variances_norm"][hvg_mask],
            s=6, alpha=0.7, color=THEME["red"], label=t("plot.hvg_label"),
        )
        ax.set_xlabel(t("plot.hvg_xlabel"), fontsize=_PLOT_CFG["font_axes"], color=THEME["text_muted"])
        ax.set_ylabel(t("plot.hvg_ylabel_disp"), fontsize=_PLOT_CFG["font_axes"], color=THEME["text_muted"])

    else:
        # Fallback: bar chart showing HVG count vs other
        n_hvg = int(hvg_mask.sum())
        n_other = int((~hvg_mask).sum())
        ax.bar(
            [f"{t('plot.hvg_label')}\n({n_hvg:,})",
             f"{t('plot.hvg_other_label')}\n({n_other:,})"],
            [n_hvg, n_other],
            color=[THEME["red"], THEME["color_ns"]],
        )
        ax.set_ylabel(t("plot.hvg_ylabel_count"), fontsize=_PLOT_CFG["font_axes"], color=THEME["text_muted"])

    ax.set_title(t("plot.hvg_title"), fontsize=13, fontweight="bold", color=THEME["text"])
    # Only add legend if artists have labels (scatter path has them, bar path doesn't)
    if ax.get_legend_handles_labels()[1]:
        _apply_legend(ax, fontsize=9, bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════
# PCA
# ══════════════════════════════════════════════════════════════════════

def plot_pca_variance(adata: ad.AnnData, n_pcs: int = 30) -> plt.Figure:
    """Elbow plot showing the variance explained by each principal component.

    Parameters
    ----------
    adata : ad.AnnData
        Must have ``adata.uns["pca"]["variance_ratio"]`` computed.
    n_pcs : int
        Maximum number of PCs to display.

    Returns
    -------
    plt.Figure
    """
    fig, ax = _themed_fig((8, 4))
    variance_ratio = adata.uns["pca"]["variance_ratio"]
    n = min(n_pcs, len(variance_ratio))
    ax.bar(range(1, n + 1), variance_ratio[:n], color=THEME["blue"], alpha=0.8)
    ax.set_xlabel(t("plot.pca_elbow_xlabel"), fontsize=_PLOT_CFG["font_axes"], color=THEME["text_muted"])
    ax.set_ylabel(t("plot.pca_elbow_ylabel"), fontsize=_PLOT_CFG["font_axes"], color=THEME["text_muted"])
    ax.set_title(t("plot.pca_elbow_title"), fontsize=13, fontweight="bold", color=THEME["text"])
    fig.tight_layout()
    return fig


def plot_pca_embedding(
    adata: ad.AnnData,
    color: str = "leiden",
    legend_loc: str = "best",
) -> plt.Figure:
    """PCA scatter plot of cells in PC1/PC2 space, colored by a metadata column.

    Categorical columns produce a discrete legend; numeric columns use
    a continuous colour gradient.

    Parameters
    ----------
    adata : ad.AnnData
        Must have ``adata.obsm["X_pca"]`` and ``adata.uns["pca"]``.
    color : str
        Column in ``adata.obs`` to colour cells by.
    legend_loc : str
        Matplotlib legend location (e.g. ``"best"``, ``"none"``).

    Returns
    -------
    plt.Figure
    """
    fig, ax = _themed_fig(_PLOT_CFG["figsize_pca"])

    pca = adata.obsm["X_pca"]

    _has_legend = False

    if color in adata.obs.columns:
        categories = adata.obs[color]
        if hasattr(categories, "cat"):
            cats = categories.cat.categories
            colors = _cluster_colors(len(cats))
            for i, cat in enumerate(cats):
                mask = categories == cat
                ax.scatter(
                    pca[mask, 0], pca[mask, 1],
                    s=_PLOT_CFG["point_size"], alpha=_PLOT_CFG["point_alpha"],
                    color=colors[i], label=str(cat), edgecolors="none",
                )
            _has_legend = legend_loc != "none"
            if _has_legend:
                _apply_legend(
                    ax, loc=legend_loc,
                    fontsize=_PLOT_CFG["font_legend"],
                    title=color, title_fontsize=9,
                    markerscale=2,
                    ncol=max(1, len(cats) // 12),
                )
        else:
            sc = ax.scatter(
                pca[:, 0], pca[:, 1],
                c=categories.values, cmap=_PLOT_CFG["cmap_expression"],
                s=_PLOT_CFG["point_size"], alpha=_PLOT_CFG["point_alpha"],
                edgecolors="none",
            )
            cbar = plt.colorbar(sc, ax=ax, label=color, shrink=0.8)
            cbar.ax.yaxis.set_tick_params(color=THEME["text_muted"])
            cbar.ax.yaxis.label.set_color(THEME["text_muted"])
            plt.setp(cbar.ax.yaxis.get_ticklabels(), color=THEME["text_muted"])
    else:
        ax.scatter(
            pca[:, 0], pca[:, 1],
            s=_PLOT_CFG["point_size"], alpha=_PLOT_CFG["point_alpha"],
            color=THEME["blue"], edgecolors="none",
        )

    var_ratio = adata.uns["pca"]["variance_ratio"]
    ax.set_xlabel(f"PC1 ({var_ratio[0]*100:.1f}%)", fontsize=_PLOT_CFG["font_axes"], color=THEME["text_muted"])
    ax.set_ylabel(f"PC2 ({var_ratio[1]*100:.1f}%)", fontsize=_PLOT_CFG["font_axes"], color=THEME["text_muted"])
    ax.set_title(t("plot.pca_embed_title"), fontsize=13, fontweight="bold", color=THEME["text"])
    fig.tight_layout()
    if _has_legend:
        fig.subplots_adjust(right=0.78)
    return fig


# ══════════════════════════════════════════════════════════════════════
# UMAP
# ══════════════════════════════════════════════════════════════════════

def plot_umap(
    adata: ad.AnnData,
    color: str = "leiden",
    title: Optional[str] = None,
    legend_loc: str = "best",
) -> plt.Figure:
    """UMAP scatter plot colored by clusters, genes, or QC metrics.

    Colour source is resolved in priority order:

    1. **Gene name** — prefers ``adata.raw`` (log-normalised expression
       of all genes) over ``adata.X`` (scaled HVG-only matrix).
    2. **Categorical obs column** — discrete legend per category.
    3. **Numeric obs column** — continuous colour gradient.
    4. **Fallback** — uniform colour if *color* is not found.

    Parameters
    ----------
    adata : ad.AnnData
        Must have ``adata.obsm["X_umap"]`` computed.
    color : str
        Gene name, ``adata.obs`` column, or QC metric to colour by.
    title : str | None
        Custom plot title; defaults to an i18n-derived UMAP title.
    legend_loc : str
        Matplotlib legend location (e.g. ``"best"``, ``"none"``).

    Returns
    -------
    plt.Figure
    """
    fig, ax = _themed_fig(_PLOT_CFG["figsize_umap"])

    umap_coords = adata.obsm["X_umap"]
    _has_legend = False

    # Determine if color refers to a gene — prefer adata.raw (log-normalized
    # expression of ALL genes) over adata.X (scaled HVG-only matrix).
    _is_gene = False
    _gene_source = None
    if adata.raw is not None and color in adata.raw.var_names:
        _is_gene = True
        _gene_source = adata.raw
    elif color in adata.var_names:
        _is_gene = True
        _gene_source = adata

    if _is_gene:
        gene_idx = list(_gene_source.var_names).index(color)
        _X = _gene_source.X
        if hasattr(_X, "toarray"):
            vals = _X[:, gene_idx].toarray().flatten()
        else:
            vals = np.asarray(_X[:, gene_idx]).flatten()
        sc = ax.scatter(
            umap_coords[:, 0], umap_coords[:, 1],
            c=vals, cmap=_PLOT_CFG["cmap_expression"],
            s=_PLOT_CFG["point_size"], alpha=_PLOT_CFG["point_alpha"],
            edgecolors="none",
        )
        cbar = plt.colorbar(sc, ax=ax, label=color, shrink=0.8)
        cbar.ax.yaxis.set_tick_params(color=THEME["text_muted"])
        cbar.ax.yaxis.label.set_color(THEME["text_muted"])
        plt.setp(cbar.ax.yaxis.get_ticklabels(), color=THEME["text_muted"])
    elif color in adata.obs.columns:
        categories = adata.obs[color]
        if hasattr(categories, "cat"):
            cats = categories.cat.categories
            colors = _cluster_colors(len(cats))
            for i, cat in enumerate(cats):
                mask = categories == cat
                ax.scatter(
                    umap_coords[mask, 0], umap_coords[mask, 1],
                    s=_PLOT_CFG["point_size"], alpha=_PLOT_CFG["point_alpha"],
                    color=colors[i], label=str(cat), edgecolors="none",
                )
            _has_legend = legend_loc != "none"
            if _has_legend:
                _apply_legend(
                    ax, loc=legend_loc,
                    fontsize=_PLOT_CFG["font_legend"],
                    title=color, title_fontsize=9,
                    markerscale=2,
                    ncol=max(1, len(cats) // 10),
                )
        else:
            sc = ax.scatter(
                umap_coords[:, 0], umap_coords[:, 1],
                c=categories.values, cmap=_PLOT_CFG["cmap_expression"],
                s=_PLOT_CFG["point_size"], alpha=_PLOT_CFG["point_alpha"],
                edgecolors="none",
            )
            cbar = plt.colorbar(sc, ax=ax, label=color, shrink=0.8)
            cbar.ax.yaxis.set_tick_params(color=THEME["text_muted"])
            cbar.ax.yaxis.label.set_color(THEME["text_muted"])
            plt.setp(cbar.ax.yaxis.get_ticklabels(), color=THEME["text_muted"])
    else:
        ax.scatter(
            umap_coords[:, 0], umap_coords[:, 1],
            s=_PLOT_CFG["point_size"], alpha=_PLOT_CFG["point_alpha"],
            color=THEME["blue"], edgecolors="none",
        )

    ax.set_xlabel(t("plot.umap_xlabel"), fontsize=_PLOT_CFG["font_axes"], color=THEME["text_muted"])
    ax.set_ylabel(t("plot.umap_ylabel"), fontsize=_PLOT_CFG["font_axes"], color=THEME["text_muted"])
    ax.set_title(title or tf("plot.umap_title_prefix", color=color),
                 fontsize=13, fontweight="bold", color=THEME["text"])
    fig.tight_layout()
    if _has_legend:
        fig.subplots_adjust(right=0.78)
    return fig


# ══════════════════════════════════════════════════════════════════════
# Marker Genes
# ══════════════════════════════════════════════════════════════════════

def plot_marker_genes_ranking(
    adata: ad.AnnData,
    n_genes: int = 5,
    n_cols: int = 4,
    groupby: str = "leiden",
) -> plt.Figure:
    """Horizontal bar plot of top marker genes per cluster, ranked by score.

    Parameters
    ----------
    adata : ad.AnnData
        Must have marker gene results in ``adata.uns["rank_genes_groups"]``.
    n_genes : int
        Number of top genes to display per group.
    n_cols : int
        Number of columns in the subplot grid.
    groupby : str
        Observation column used for grouping (e.g. ``"leiden"``).

    Returns
    -------
    plt.Figure
    """
    groups = list(adata.obs[groupby].cat.categories)
    n_groups = len(groups)
    n_rows = (n_groups + n_cols - 1) // n_cols
    colors = _cluster_colors(n_groups)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(3.5 * n_cols, 3 * n_rows),
    )
    axes = np.atleast_2d(axes)
    fig.patch.set_facecolor(THEME["bg"])
    plt.rcParams["font.family"] = THEME["font_family"]

    for idx, group in enumerate(groups):
        row, col = divmod(idx, n_cols)
        ax = axes[row, col]
        ax.set_facecolor(THEME["surface"])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color(THEME["border"])
        ax.spines["bottom"].set_color(THEME["border"])
        ax.tick_params(colors=THEME["text_muted"], labelsize=THEME["font_ticks"])

        try:
            from engine.scrna_pipeline import get_marker_genes_df
            df = get_marker_genes_df(adata, group=group, n_genes=n_genes)
            ax.barh(
                range(len(df)),
                df["scores"].values,
                color=colors[idx], alpha=0.85,
            )
            ax.set_yticks(range(len(df)))
            ax.set_yticklabels(df["names"].values, fontsize=7, color=THEME["text_muted"])
            ax.invert_yaxis()
        except Exception:
            ax.text(0.5, 0.5, t("plot.marker_no_data"), ha="center", va="center",
                    transform=ax.transAxes, color=THEME["text_muted"])

        ax.set_title(str(group), fontsize=10, fontweight="bold", color=THEME["text"])
        ax.set_xlabel(t("plot.marker_score_label"), fontsize=8, color=THEME["text_muted"])

    # Hide unused axes
    for idx in range(n_groups, n_rows * n_cols):
        row, col = divmod(idx, n_cols)
        axes[row, col].axis("off")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.suptitle(tf("plot.marker_ranking_title", groupby=groupby.replace('_', ' ').title()),
                 fontsize=14, fontweight="bold", color=THEME["text"], y=0.98)
    return fig


def plot_dotplot(
    adata: ad.AnnData,
    n_genes: int = 5,
    groupby: str = "leiden",
    marker_genes: Optional[dict] = None,
) -> plt.Figure:
    """Dotplot showing mean expression & fraction of cells expressing
    top marker genes per cluster.

    Uses ``adata.raw`` (log-normalized expression of all genes) when
    available, which gives more meaningful values than the scaled matrix.

    Parameters
    ----------
    adata : AnnData
        Annotated data with ``rank_genes_groups`` in ``.uns``.
    n_genes : int
        Number of top marker genes to show per cluster.
    groupby : str
        Column in ``adata.obs`` to group by (e.g. ``"leiden"``).
    marker_genes : dict, optional
        If provided, a dict mapping group labels → gene lists.
        Otherwise, genes are extracted from ``rank_genes_groups``.
    """
    import scanpy as sc  # local import — only used here

    # ── Resolve which AnnData view to use for expression values ───
    # adata.raw has log-normalized counts for ALL genes (pre-HVG subset);
    # adata.X is the scaled HVG-only matrix — not ideal for a dotplot.
    use_raw = adata.raw is not None
    _ad = adata.raw if use_raw else adata
    available_genes = set(_ad.var_names)

    # ── Build gene list per cluster ───────────────────────────────
    if marker_genes is None and "rank_genes_groups" not in adata.uns:
        fig, ax = plt.subplots(figsize=(6, 3))
        fig.patch.set_facecolor(THEME["bg"])
        ax.set_facecolor(THEME["surface"])
        ax.text(0.5, 0.5, t("plot.marker_no_genes"), ha="center", va="center",
                color=THEME["text_muted"])
        ax.axis("off")
        return fig

    if marker_genes is None:
        _groupby_rgg = adata.uns["rank_genes_groups"]["params"]["groupby"]
        groups = list(adata.obs[_groupby_rgg].cat.categories)
        marker_genes = {}
        for g in groups:
            df = sc.get.rank_genes_groups_df(adata, group=str(g))
            top = [gn for gn in df["names"].tolist() if gn in available_genes][:n_genes]
            marker_genes[str(g)] = top
    else:
        groups = list(marker_genes.keys())

    # Collect unique genes preserving cluster order
    genes = []
    group_boundaries = []  # for vertical separator lines
    for gene_list in marker_genes.values():
        _start = len(genes)
        for g in gene_list:
            if g not in genes and g in available_genes:
                genes.append(g)
        group_boundaries.append(len(genes))

    if not genes:
        fig, ax = plt.subplots(figsize=(6, 3))
        fig.patch.set_facecolor(THEME["bg"])
        ax.set_facecolor(THEME["surface"])
        ax.text(0.5, 0.5, t("plot.marker_no_genes"), ha="center", va="center",
                color=THEME["text_muted"])
        ax.axis("off")
        return fig

    n_g = len(genes)
    n_grp = len(groups)

    # ── Compute mean expression and fraction expressing ───────────
    mean_expr = np.zeros((n_grp, n_g))
    frac_expr = np.zeros((n_grp, n_g))

    for i, group in enumerate(groups):
        mask = (adata.obs[groupby] == group).values
        if use_raw:
            subset = adata.raw[mask, genes]
        else:
            subset = adata[mask, genes]
        X = subset.X.toarray() if hasattr(subset.X, "toarray") else np.asarray(subset.X)
        mean_expr[i] = X.mean(axis=0)
        frac_expr[i] = (X > 0).mean(axis=0)

    # ── Figure layout ─────────────────────────────────────────────
    _w = max(10, 0.5 * n_g + 3)
    _h = max(4, 0.45 * n_grp + 2)
    fig, ax = _themed_fig((_w, _h))

    _DOT_SCALE = 300
    _cmap = _PLOT_CFG["cmap_expression"]
    _norm = plt.Normalize(vmin=mean_expr.min(), vmax=mean_expr.max())

    for i in range(n_grp):
        for j in range(n_g):
            ax.scatter(
                j, i,
                s=frac_expr[i, j] * _DOT_SCALE,
                c=[mean_expr[i, j]],
                cmap=_cmap,
                norm=_norm,
                edgecolors=THEME["text_subtle"],
                linewidths=0.5,
                zorder=3,
            )

    # Cluster separator lines (vertical dashed)
    for bnd in group_boundaries[:-1]:
        ax.axvline(bnd - 0.5, color=THEME["border"], linewidth=0.8, linestyle="--", zorder=1)

    ax.set_xticks(range(n_g))
    ax.set_xticklabels(genes, rotation=90, fontsize=8, fontstyle="italic",
                       color=THEME["text_muted"])
    ax.set_yticks(range(n_grp))
    ax.set_yticklabels(groups, fontsize=9, color=THEME["text_muted"])
    ax.set_xlabel(t("plot.dotplot_xlabel"), fontsize=_PLOT_CFG["font_axes"], color=THEME["text_muted"])
    ax.set_ylabel(t("plot.dotplot_ylabel"), fontsize=_PLOT_CFG["font_axes"], color=THEME["text_muted"])
    ax.set_title(t("plot.dotplot_title"), fontsize=13, fontweight="bold", color=THEME["text"])
    ax.set_xlim(-0.5, n_g - 0.5)
    ax.set_ylim(-0.5, n_grp - 0.5)
    ax.invert_yaxis()

    # ── Colorbar (mean expression) ────────────────────────────────
    sm = plt.cm.ScalarMappable(cmap=_cmap, norm=_norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, label=t("plot.dotplot_cbar_label"),
                        shrink=0.6, pad=0.02)
    cbar.ax.tick_params(labelsize=8, colors=THEME["text_muted"])
    cbar.ax.yaxis.label.set_color(THEME["text_muted"])

    # ── Size legend (fraction of cells) ───────────────────────────
    _legend_fracs = [0.25, 0.50, 0.75, 1.00]
    _legend_handles = []
    for frac in _legend_fracs:
        _legend_handles.append(
            ax.scatter(
                [], [], s=frac * _DOT_SCALE,
                c=THEME["text_subtle"], edgecolors=THEME["text_subtle"], linewidths=0.5,
                label=f"{int(frac * 100)}%",
            )
        )
    leg = ax.legend(
        handles=_legend_handles,
        title=t("plot.dotplot_size_legend"),
        loc="upper left",
        bbox_to_anchor=(1.15, 1.0),
        frameon=True,
        fontsize=8,
        title_fontsize=9,
        labelspacing=1.2,
        handletextpad=1.0,
        facecolor=THEME["legend_bg"],
        edgecolor=THEME["legend_edge"],
        labelcolor=THEME["text"],
    )
    leg.get_title().set_color(THEME["text"])

    fig.tight_layout()
    return fig


def plot_marker_heatmap(
    adata: ad.AnnData,
    n_genes: int = 5,
    groupby: str = "leiden",
) -> plt.Figure:
    """Heatmap of top marker genes per cluster.

    Cells are sorted by cluster and the expression of each top marker
    gene is displayed as a row.  A colour bar above the heatmap
    indicates cluster identity.

    Parameters
    ----------
    adata : ad.AnnData
        Must have marker gene results in ``adata.uns["rank_genes_groups"]``.
    n_genes : int
        Number of top genes per cluster.
    groupby : str
        Observation column used for grouping.

    Returns
    -------
    plt.Figure
    """
    groups = list(adata.obs[groupby].cat.categories)
    from engine.scrna_pipeline import get_marker_genes_df

    # ── Resolve which AnnData view to use for expression values ───
    # adata.raw has log-normalized counts for ALL genes (pre-HVG subset);
    # adata.X is the scaled HVG-only matrix — not ideal for a heatmap.
    use_raw = adata.raw is not None
    _ad = adata.raw if use_raw else adata
    available_genes = set(_ad.var_names)

    # Collect top genes per cluster
    top_genes = []
    for g in groups:
        try:
            df = get_marker_genes_df(adata, group=g, n_genes=n_genes)
            for gene in df["names"].values:
                if gene in available_genes and gene not in top_genes:
                    top_genes.append(gene)
        except Exception:
            pass

    if not top_genes:
        fig, ax = plt.subplots(figsize=(6, 3))
        fig.patch.set_facecolor(THEME["bg"])
        ax.set_facecolor(THEME["surface"])
        ax.text(0.5, 0.5, t("plot.heatmap_no_genes"), ha="center", va="center",
                color=THEME["text_muted"])
        ax.axis("off")
        return fig

    # Sort cells by cluster for clean visualization
    order = adata.obs.sort_values(groupby).index
    subset = _ad[order, top_genes]

    if hasattr(subset.X, "toarray"):
        X = subset.X.toarray()
    else:
        X = np.asarray(subset.X)

    n_cells = X.shape[0]
    n_genes_plot = len(top_genes)

    # Use GridSpec: color-bar | heatmap | colorbar-axis, with sharex
    # so the colorbar stealing space from ax_heat also affects ax_bar.
    fig = plt.figure(figsize=_PLOT_CFG["figsize_heatmap"])
    fig.patch.set_facecolor(THEME["bg"])
    plt.rcParams["font.family"] = THEME["font_family"]

    gs = fig.add_gridspec(
        2, 2,
        height_ratios=[0.03, 1],
        width_ratios=[1, 0.03],
        hspace=0.05,
        wspace=0.03,
    )
    ax_bar = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[1, 0], sharex=ax_bar)
    ax_cbar = fig.add_subplot(gs[1, 1])

    for _ax in [ax_bar, ax_heat, ax_cbar]:
        _ax.set_facecolor(THEME["bg"])

    # ── Cluster color bar (RGB image — no colormap normalisation) ───
    cluster_labels = adata.obs.loc[order, groupby]
    colors = _cluster_colors(len(groups))
    color_map = {g: colors[i] for i, g in enumerate(groups)}
    bar_rgb = np.array(
        [matplotlib.colors.to_rgb(color_map[c]) for c in cluster_labels]
    )
    ax_bar.imshow(
        bar_rgb[np.newaxis, :, :],
        aspect="auto",
        interpolation="nearest",
        extent=(-0.5, n_cells - 0.5, 0, 1),
    )
    ax_bar.set_xticks([])
    ax_bar.set_yticks([])
    ax_bar.set_title(t("plot.heatmap_cluster_title"), fontsize=9, color=THEME["text"])

    # ── Heatmap ─────────────────────────────────────────────────────
    im = ax_heat.imshow(
        X.T,
        aspect="auto",
        cmap=_PLOT_CFG["cmap_heatmap"],
        interpolation="nearest",
        extent=(-0.5, n_cells - 0.5, n_genes_plot - 0.5, -0.5),
    )
    ax_heat.set_yticks(range(n_genes_plot))
    ax_heat.set_yticklabels(top_genes, fontsize=7, color=THEME["text_muted"])
    ax_heat.set_xlabel(t("plot.heatmap_xlabel"), fontsize=10, color=THEME["text_muted"])
    ax_heat.set_xticks([])
    ax_heat.tick_params(colors=THEME["text_muted"])

    # Colorbar in its own dedicated axis → no space stolen from heatmap
    cb = fig.colorbar(im, cax=ax_cbar, label=t("plot.heatmap_cbar_label"))
    cb.ax.yaxis.label.set_color(THEME["text_muted"])
    cb.ax.tick_params(colors=THEME["text_muted"])

    # Hide top-right cell (no colorbar needed for the cluster bar)
    ax_empty = fig.add_subplot(gs[0, 1])
    ax_empty.axis("off")

    fig.suptitle(
        t("plot.heatmap_title"),
        fontsize=13, fontweight="bold", color=THEME["text"], y=1.01,
    )
    # GridSpec with sharex is not compatible with tight_layout; use
    # subplots_adjust instead to avoid the matplotlib warning.
    fig.subplots_adjust(left=0.12, right=0.88, top=0.92, bottom=0.08)
    return fig
