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

from config import SCRNA_CONFIG
from i18n import t, tf

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
    fig.savefig(buf, format=fmt, dpi=_PLOT_CFG["dpi"], bbox_inches="tight")
    buf.seek(0)
    return buf.getvalue()


def _apply_legend(ax, loc: str = "best", fontsize: int = 9, **extra_kw) -> None:
    """Create and style a legend with consistent scientific aesthetics.

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
        framealpha=0.9,
        facecolor="white",
        edgecolor="#E0E0E0",
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
    """Violin plots of QC metrics."""
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

    for ax, key in zip(axes, keys):
        data = adata.obs[key].dropna().values
        parts = ax.violinplot(data, showmeans=True, showmedians=True)
        for pc in parts["bodies"]:
            pc.set_facecolor("#1E88E5")
            pc.set_alpha(0.6)
        ax.set_title(key, fontsize=10)
        ax.set_ylabel("")
        ax.set_xticks([])

    fig.suptitle(t("plot.qc_title"), fontsize=13, fontweight="bold")
    fig.tight_layout()
    return fig


def plot_qc_scatter(
    adata: ad.AnnData,
    x: str = "total_counts",
    y: str = "n_genes_by_counts",
    color: str = "pct_counts_mt",
) -> plt.Figure:
    """Scatter plot of two QC metrics, colored by a third."""
    fig, ax = plt.subplots(figsize=(8, 6))

    if color in adata.obs.columns:
        sc = ax.scatter(
            adata.obs[x], adata.obs[y],
            c=adata.obs[color],
            cmap="YlOrRd",
            s=_PLOT_CFG["point_size"],
            alpha=_PLOT_CFG["point_alpha"],
            edgecolors="none",
        )
        plt.colorbar(sc, ax=ax, label=color, shrink=0.8)
    else:
        ax.scatter(
            adata.obs[x], adata.obs[y],
            s=_PLOT_CFG["point_size"],
            alpha=_PLOT_CFG["point_alpha"],
            color="#1E88E5",
            edgecolors="none",
        )

    ax.set_xlabel(x, fontsize=11)
    ax.set_ylabel(y, fontsize=11)
    ax.set_title(t("plot.qc_scatter_title"), fontsize=13, fontweight="bold")
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
    fig, ax = plt.subplots(figsize=(8, 5))

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
            s=4, alpha=0.3, color="#D5D5D5", label=t("plot.hvg_other_label"),
        )
        ax.scatter(
            var_df["means"][hvg_mask],
            var_df["dispersions_norm"][hvg_mask],
            s=6, alpha=0.6, color="#E53935", label=t("plot.hvg_label"),
        )
        ax.set_xlabel(t("plot.hvg_xlabel"), fontsize=11)
        ax.set_ylabel(t("plot.hvg_ylabel_disp"), fontsize=11)

    elif "variances_norm" in var_df.columns and "means" in var_df.columns:
        # seurat_v3 flavor → scatter of mean vs normalized variance
        ax.scatter(
            var_df["means"][~hvg_mask],
            var_df["variances_norm"][~hvg_mask],
            s=4, alpha=0.3, color="#D5D5D5", label=t("plot.hvg_other_label"),
        )
        ax.scatter(
            var_df["means"][hvg_mask],
            var_df["variances_norm"][hvg_mask],
            s=6, alpha=0.6, color="#E53935", label=t("plot.hvg_label"),
        )
        ax.set_xlabel(t("plot.hvg_xlabel"), fontsize=11)
        ax.set_ylabel(t("plot.hvg_ylabel_disp"), fontsize=11)

    else:
        # Fallback: bar chart showing HVG count vs other
        n_hvg = int(hvg_mask.sum())
        n_other = int((~hvg_mask).sum())
        ax.bar(
            [f"{t('plot.hvg_label')}\n({n_hvg:,})",
             f"{t('plot.hvg_other_label')}\n({n_other:,})"],
            [n_hvg, n_other],
            color=["#E53935", "#D5D5D5"],
        )
        ax.set_ylabel(t("plot.hvg_ylabel_count"), fontsize=11)

    ax.set_title(t("plot.hvg_title"), fontsize=13, fontweight="bold")
    # Only add legend if artists have labels (scatter path has them, bar path doesn't)
    if ax.get_legend_handles_labels()[1]:
        ax.legend(fontsize=9)
    fig.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════
# PCA
# ══════════════════════════════════════════════════════════════════════

def plot_pca_variance(adata: ad.AnnData, n_pcs: int = 30) -> plt.Figure:
    """Elbow plot showing variance explained per PC."""
    fig, ax = plt.subplots(figsize=(8, 4))
    variance_ratio = adata.uns["pca"]["variance_ratio"]
    n = min(n_pcs, len(variance_ratio))
    ax.bar(range(1, n + 1), variance_ratio[:n], color="#1E88E5", alpha=0.7)
    ax.set_xlabel(t("plot.pca_elbow_xlabel"), fontsize=11)
    ax.set_ylabel(t("plot.pca_elbow_ylabel"), fontsize=11)
    ax.set_title(t("plot.pca_elbow_title"), fontsize=13, fontweight="bold")
    fig.tight_layout()
    return fig


def plot_pca_embedding(
    adata: ad.AnnData,
    color: str = "leiden",
    legend_loc: str = "best",
) -> plt.Figure:
    """PCA scatter plot colored by a metadata column."""
    fig, ax = plt.subplots(figsize=_PLOT_CFG["figsize_pca"])

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
                c=categories.values, cmap="viridis",
                s=_PLOT_CFG["point_size"], alpha=_PLOT_CFG["point_alpha"],
                edgecolors="none",
            )
            plt.colorbar(sc, ax=ax, label=color, shrink=0.8)
    else:
        ax.scatter(
            pca[:, 0], pca[:, 1],
            s=_PLOT_CFG["point_size"], alpha=_PLOT_CFG["point_alpha"],
            color="#1E88E5", edgecolors="none",
        )

    var_ratio = adata.uns["pca"]["variance_ratio"]
    ax.set_xlabel(f"PC1 ({var_ratio[0]*100:.1f}%)", fontsize=11)
    ax.set_ylabel(f"PC2 ({var_ratio[1]*100:.1f}%)", fontsize=11)
    ax.set_title(t("plot.pca_embed_title"), fontsize=13, fontweight="bold")
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
    """UMAP scatter plot colored by clusters, genes, or QC metrics."""
    fig, ax = plt.subplots(figsize=_PLOT_CFG["figsize_umap"])

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
        plt.colorbar(sc, ax=ax, label=color, shrink=0.8)
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
                c=categories.values, cmap="viridis",
                s=_PLOT_CFG["point_size"], alpha=_PLOT_CFG["point_alpha"],
                edgecolors="none",
            )
            plt.colorbar(sc, ax=ax, label=color, shrink=0.8)
    else:
        ax.scatter(
            umap_coords[:, 0], umap_coords[:, 1],
            s=_PLOT_CFG["point_size"], alpha=_PLOT_CFG["point_alpha"],
            color="#1E88E5", edgecolors="none",
        )

    ax.set_xlabel(t("plot.umap_xlabel"), fontsize=11)
    ax.set_ylabel(t("plot.umap_ylabel"), fontsize=11)
    ax.set_title(title or tf("plot.umap_title_prefix", color=color), fontsize=13, fontweight="bold")
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
    """Bar plot of top marker genes per group (ranked by score)."""
    groups = list(adata.obs[groupby].cat.categories)
    n_groups = len(groups)
    n_rows = (n_groups + n_cols - 1) // n_cols
    colors = _cluster_colors(n_groups)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(3.5 * n_cols, 3 * n_rows),
    )
    axes = np.atleast_2d(axes)

    for idx, group in enumerate(groups):
        row, col = divmod(idx, n_cols)
        ax = axes[row, col]

        try:
            from scrna_pipeline import get_marker_genes_df
            df = get_marker_genes_df(adata, group=group, n_genes=n_genes)
            ax.barh(
                range(len(df)),
                df["scores"].values,
                color=colors[idx], alpha=0.8,
            )
            ax.set_yticks(range(len(df)))
            ax.set_yticklabels(df["names"].values, fontsize=7)
            ax.invert_yaxis()
        except Exception:
            ax.text(0.5, 0.5, t("plot.marker_no_data"), ha="center", va="center", transform=ax.transAxes)

        ax.set_title(str(group), fontsize=10, fontweight="bold")
        ax.set_xlabel(t("plot.marker_score_label"), fontsize=8)

    # Hide unused axes
    for idx in range(n_groups, n_rows * n_cols):
        row, col = divmod(idx, n_cols)
        axes[row, col].axis("off")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.suptitle(tf("plot.marker_ranking_title", groupby=groupby.replace('_', ' ').title()), fontsize=14, fontweight="bold", y=0.98)
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
        ax.text(0.5, 0.5, t("plot.marker_no_genes"), ha="center", va="center")
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
        ax.text(0.5, 0.5, t("plot.marker_no_genes"), ha="center", va="center")
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
    fig, ax = plt.subplots(figsize=(_w, _h))

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
                edgecolors="grey",
                linewidths=0.5,
                zorder=3,
            )

    # Cluster separator lines (vertical dashed)
    for bnd in group_boundaries[:-1]:
        ax.axvline(bnd - 0.5, color="#CCCCCC", linewidth=0.8, linestyle="--", zorder=1)

    ax.set_xticks(range(n_g))
    ax.set_xticklabels(genes, rotation=90, fontsize=8, fontstyle="italic")
    ax.set_yticks(range(n_grp))
    ax.set_yticklabels(groups, fontsize=9)
    ax.set_xlabel(t("plot.dotplot_xlabel"), fontsize=11)
    ax.set_ylabel(t("plot.dotplot_ylabel"), fontsize=11)
    ax.set_title(t("plot.dotplot_title"), fontsize=13, fontweight="bold")
    ax.set_xlim(-0.5, n_g - 0.5)
    ax.set_ylim(-0.5, n_grp - 0.5)
    ax.invert_yaxis()

    # ── Colorbar (mean expression) ────────────────────────────────
    sm = plt.cm.ScalarMappable(cmap=_cmap, norm=_norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, label=t("plot.dotplot_cbar_label"),
                        shrink=0.6, pad=0.02)
    cbar.ax.tick_params(labelsize=8)

    # ── Size legend (fraction of cells) ───────────────────────────
    _legend_fracs = [0.25, 0.50, 0.75, 1.00]
    _legend_handles = []
    for frac in _legend_fracs:
        _legend_handles.append(
            ax.scatter(
                [], [], s=frac * _DOT_SCALE,
                c="grey", edgecolors="grey", linewidths=0.5,
                label=f"{int(frac * 100)}%",
            )
        )
    ax.legend(
        handles=_legend_handles,
        title=t("plot.dotplot_size_legend"),
        loc="upper left",
        bbox_to_anchor=(1.15, 1.0),
        frameon=True,
        fontsize=8,
        title_fontsize=9,
        labelspacing=1.2,
        handletextpad=1.0,
    )

    fig.tight_layout()
    return fig


def plot_marker_heatmap(
    adata: ad.AnnData,
    n_genes: int = 5,
    groupby: str = "leiden",
) -> plt.Figure:
    """Heatmap of top marker genes per cluster."""
    groups = list(adata.obs[groupby].cat.categories)
    from scrna_pipeline import get_marker_genes_df

    # Collect top genes per cluster
    top_genes = []
    for g in groups:
        try:
            df = get_marker_genes_df(adata, group=g, n_genes=n_genes)
            for gene in df["names"].values:
                if gene in adata.var_names and gene not in top_genes:
                    top_genes.append(gene)
        except Exception:
            pass

    if not top_genes:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, t("plot.heatmap_no_genes"), ha="center", va="center")
        ax.axis("off")
        return fig

    # Sort cells by cluster for clean visualization
    order = adata.obs.sort_values(groupby).index
    subset = adata[order, top_genes]

    if hasattr(subset.X, "toarray"):
        X = subset.X.toarray()
    else:
        X = np.asarray(subset.X)

    n_cells = X.shape[0]
    n_genes_plot = len(top_genes)

    # Use GridSpec: color-bar | heatmap | colorbar-axis, with sharex
    # so the colorbar stealing space from ax_heat also affects ax_bar.
    fig = plt.figure(figsize=_PLOT_CFG["figsize_heatmap"])
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
    ax_bar.set_title(t("plot.heatmap_cluster_title"), fontsize=9)

    # ── Heatmap ─────────────────────────────────────────────────────
    im = ax_heat.imshow(
        X.T,
        aspect="auto",
        cmap=_PLOT_CFG["cmap_heatmap"],
        interpolation="nearest",
        extent=(-0.5, n_cells - 0.5, n_genes_plot - 0.5, -0.5),
    )
    ax_heat.set_yticks(range(n_genes_plot))
    ax_heat.set_yticklabels(top_genes, fontsize=7)
    ax_heat.set_xlabel(t("plot.heatmap_xlabel"), fontsize=10)
    ax_heat.set_xticks([])

    # Colorbar in its own dedicated axis → no space stolen from heatmap
    fig.colorbar(im, cax=ax_cbar, label=t("plot.heatmap_cbar_label"))

    # Hide top-right cell (no colorbar needed for the cluster bar)
    ax_empty = fig.add_subplot(gs[0, 1])
    ax_empty.axis("off")

    fig.suptitle(
        t("plot.heatmap_title"),
        fontsize=13, fontweight="bold", y=1.01,
    )
    # GridSpec with sharex is not compatible with tight_layout; use
    # subplots_adjust instead to avoid the matplotlib warning.
    fig.subplots_adjust(left=0.12, right=0.88, top=0.92, bottom=0.08)
    return fig
