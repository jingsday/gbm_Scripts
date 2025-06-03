import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def dim_heatmap(
    adata,
    dims=[0],
    nfeatures=30,
    cells=None,
    reduction='X_pca',
    disp_min=-2.5,
    disp_max=None,
    balanced=True,
    ncol=None,
    fast=True,
    combine=True,
    layer=None,
    figsize=(4, 4)
):
    """
    Plot heatmaps of top features correlated with PCA (or other) dimensions.
    Features (genes) are shown in rows, cells in columns — like Seurat's DimHeatmap.

    Parameters:
    - adata: AnnData object
    - dims: list of ints, PCA dims to plot (0-based)
    - nfeatures: number of top features per dimension (positive + negative)
    - cells: number of cells to use (or list of indices). Default selects 500 per dim.
    - reduction: str, name of reduction in adata.obsm (default 'X_pca')
    - disp_min, disp_max: float, clamp values for expression
    - balanced: bool, select balanced top positive and negative features/cells
    - ncol: int, number of columns for subplot layout
    - fast: bool, not used here but kept for compatibility
    - combine: bool, if True return combined figure
    - layer: str or None, layer of adata to use for expression data
    - figsize: tuple, size per subplot

    Returns:
    - matplotlib figure if combine=True, else shows plots inline
    """

    if ncol is None:
        ncol = 3 if len(dims) > 2 else len(dims)

    if reduction not in adata.obsm:
        raise ValueError(f"Reduction '{reduction}' not found in adata.obsm. Available: {list(adata.obsm.keys())}")
    embedding = adata.obsm[reduction]

    if cells is None:
        cells = []
        for dim in dims:
            pc_scores = embedding[:, dim]
            if balanced:
                n_half = 250  # 250 + 250 = 500 total cells
                pos_cells = np.argsort(pc_scores)[-n_half:]
                neg_cells = np.argsort(pc_scores)[:n_half]
                dim_cells = np.concatenate([neg_cells, pos_cells])
            else:
                dim_cells = np.argsort(np.abs(pc_scores))[-500:]
            cells.append(dim_cells)
    elif isinstance(cells, (list, np.ndarray)) and not isinstance(cells[0], (list, np.ndarray)):
        cells = [cells] * len(dims)

    if 'PCs' not in adata.varm:
        raise ValueError("Feature loadings not found in adata.varm['PCs']")
    loadings = adata.varm['PCs']

    nrow = (len(dims) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(figsize[0] * ncol, figsize[1] * nrow), squeeze=False)

    def clamp(data, vmin, vmax):
        return np.clip(data, vmin, vmax)

    for i, dim in enumerate(dims):
        feature_scores = loadings[:, dim]
        if balanced:
            n_half = nfeatures // 2
            pos_feats = np.argsort(feature_scores)[-n_half:]
            neg_feats = np.argsort(feature_scores)[:n_half]
            feat_idx = np.concatenate([neg_feats, pos_feats])
        else:
            feat_idx = np.argsort(np.abs(feature_scores))[-nfeatures:]

        features = adata.var_names[feat_idx]
        dim_cells = cells[i]

        if layer is not None:
            expr = adata.layers[layer][dim_cells, :][:, feat_idx]
        else:
            expr = adata.X[dim_cells, :][:, feat_idx]
            if hasattr(expr, "toarray"):
                expr = expr.toarray()

        disp_max_val = disp_max if disp_max is not None else (2.5 if layer is None else 6)
        expr = clamp(expr, disp_min, disp_max_val)

        # Transpose: genes in rows, cells in columns
        expr = expr.T

        ax = axes.flatten()[i]
        sns.heatmap(
            expr,
            ax=ax,
            yticklabels=features,
            xticklabels=False,
            cmap='coolwarm',
            center=0,
            cbar=i == 0
        )
        ax.set_title(f"PC{dim+1}")#        ax.set_title(f"{reduction.upper()} Dim {dim+1}")
        ax.set_xlabel("Cells")
        ax.set_ylabel("Features")

    for j in range(len(dims), nrow * ncol):
        fig.delaxes(axes.flatten()[j])

    plt.tight_layout()

    if combine:
        return fig
    else:
        plt.show()
