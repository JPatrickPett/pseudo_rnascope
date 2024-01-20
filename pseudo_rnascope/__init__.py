"""Approximate missing features from higher dimensionality data neighbours"""
__version__ = "0.0.5"

import numpy as np
import scipy as sp


def get_array(adata, gene_symbol):
    exp_mat = adata.X[:, adata.var_names == gene_symbol]
    if sp.sparse.issparse(exp_mat):
        exp_mat = exp_mat.todense()
    return np.array(exp_mat).flatten()


def scale(x, max_val=255, vmin=None, vmax=None):
    if vmin:
        x[x < vmin] = 0
    if vmax:
        x[x > vmax] = vmax
    return max_val * (x - np.min(x)) / (np.max(x) - np.min(x))


def mix_colors(c1, c2):
    return [round(2 * np.mean([a, b])) for a, b in zip(c1, c2)]


def rgb2hex(rgb_tuple):
    r, g, b = rgb_tuple
    return "#{:02x}{:02x}{:02x}".format(r, g, b)


def hex2rgb(hexcode):
    return tuple(map(ord, hexcode[1:].decode("hex")))


def graph_smooth(expression_matrix, neighbor_matrix):
    """smoothing over the knn graph"""
    return ((neighbor_matrix @ expression_matrix) + expression_matrix) / (
        neighbor_matrix > 0
    ).sum(axis=1)


def add_pseudo_rna_scope(
    adata,
    channel1="gene1",
    channel2="gene2",
    channel1_vmin=None,
    channel1_vmax=None,
    channel1_color=(1, 0, 0),  # RGB
    channel2_vmin=None,
    channel2_vmax=None,
    channel2_color=(0, 1, 0),  # RGB
    auto_range_quantiles=(0.2, 0.9),
    knn_smooth=False,
    na_thr=0,
):
    """
    Add information to an existing `anndata` object, so that two genes can be plotted with two colors. Combined gene expression is represented through additive color mixing.

    Adds a column `pseudo_RNAscope` to the `anndata.obs` dataframe, which can be selected in `scanpy` plotting functions and saves the mixed colors in a corresponding `anndata.uns["pseudo_RNAscope_colors"]` entry.
    Adds a column `pseudo_RNAscope_alpha` to the `anndata.obs` dataframe, which can be passed as alpha values to plotting functions, so that low expression spots are transparent.

    Parameters
    ----------
    adata :
        `anndata` object with basic preprocessing
    channel<x> :
        gene name in `adata.var_names` to plot in this channel
    channel<x>_vmin :
        minimum gene expression values to plot, lower values will be dropped
    channel<x>_vmax :
        maximum gene expression values to plot, larger values will have the same color
    channel<x>_color :
        color to use for this channel as a triple (RGB format, each value in the range [0,1])
    auto_range_quantiles :
        tuple with lower and upper quantiles for automatic selection of `channel<x>_vmin` and `channel<x>_vmax`, respecively. Setting the bounds explicitly overrides this option.
    knn_smooth :
        whether to run nearest-neighbor smoothing on the gene expression values to mitigate sparsity (default: False)

    Returns
    -------
    A dictionary with selected dynamic ranges for plotting (e.g. for adding colorbars).
    """

    for cnl in [channel1, channel2]:
        if cnl not in adata.var_names:
            raise ValueError(f"gene '{cnl}' not in 'adata.var_names'")

    if knn_smooth:
        channel1_vec = np.array(
            graph_smooth(adata[:, channel1].X, adata.obsp["connectivities"])
        ).flatten()
        channel2_vec = np.array(
            graph_smooth(adata[:, channel2].X, adata.obsp["connectivities"])
        ).flatten()
    else:
        channel1_vec = get_array(adata, channel1)
        channel2_vec = get_array(adata, channel2)

    if not channel1_vmin and auto_range_quantiles:
        channel1_vmin = np.quantile(
            channel1_vec[channel1_vec > 0], q=[auto_range_quantiles[0]]
        )
    if not channel1_vmax and auto_range_quantiles:
        channel1_vmax = np.quantile(
            channel1_vec[channel1_vec > 0], q=[auto_range_quantiles[1]]
        )
    if not channel2_vmin and auto_range_quantiles:
        channel2_vmin = np.quantile(
            channel2_vec[channel2_vec > 0], q=[auto_range_quantiles[0]]
        )
    if not channel2_vmax and auto_range_quantiles:
        channel2_vmax = np.quantile(
            channel2_vec[channel2_vec > 0], q=[auto_range_quantiles[1]]
        )

    ### ASSIGN COLORS

    rgb_values = [
        mix_colors(np.array(channel1_color) * x, np.array(channel2_color) * y)
        for x, y in zip(
            scale(channel1_vec, vmin=channel1_vmin, vmax=channel1_vmax),
            scale(channel2_vec, vmin=channel2_vmin, vmax=channel2_vmax),
        )
    ]

    ### STORE ANNDATA COLUMN AND VALUES

    adata.obs["pseudo_RNAscope"] = [
        x if max(y) >= na_thr else np.nan for x, y in zip(adata.obs_names, rgb_values)
    ]
    adata.uns["pseudo_RNAscope_colors"] = [
        rgb2hex(y) for x, y in zip(adata.obs_names, rgb_values) if max(y) >= na_thr
    ]
    adata.obs["pseudo_RNAscope"] = adata.obs["pseudo_RNAscope"].astype("category")
    adata.obs["pseudo_RNAscope_alpha"] = [
        y if max(z) >= na_thr else 0
        for y, z in zip(
            scale(np.array([np.mean(x) for x in rgb_values]), max_val=1), rgb_values
        )
    ]

    return {
        "channel1_vmin": channel1_vmin,
        "channel1_vmax": channel1_vmax,
        "channel2_vmin": channel2_vmin,
        "channel2_vmax": channel2_vmax,
    }
