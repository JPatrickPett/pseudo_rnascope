"""Approximate missing features from higher dimensionality data neighbours"""
__version__ = "0.0.5"

import re
import numpy as np
import scipy as sp


def get_array(adata, gene_symbol):
    exp_mat = adata.X[:, adata.var_names == gene_symbol]
    if sp.sparse.issparse(exp_mat):
        exp_mat = exp_mat.todense()
    return np.array(exp_mat).flatten()


def scale(x, max_val=255, vmin=None, vmax=None):
    x = x.copy().astype(float)
    if vmin:
        x[x < vmin] = np.nan
    if vmax:
        x[x > vmax] = vmax
    return max_val * (x - np.nanmin(x)) / (np.nanmax(x) - np.nanmin(x))


def mix_colors(colors, gamma=4.5):
    assert gamma > 0
    return np.power(sum([c**gamma for c in colors]) / len(colors), 1 / gamma)


def rgb2hex(rgb_tuple):
    r, g, b = np.round(rgb_tuple).astype(int)
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
    channel1_vmin=None,
    channel1_vmax=None,
    channel1_color=(1, 0, 0),  # RGB
    auto_range_quantiles=(0.2, 0.9),
    knn_smooth=False,
    gamma=4.5,
    **kwargs,
):
    """
    Add information to an existing `anndata` object, so that multiple genes can be plotted with different colors. Combined gene expression is represented through additive color mixing. An arbitrary number of channels can be added as **kwargs ("channel<x>" and "channel<x>_color" mandatory, "channel<x>_<vmin|vmax>" optional).

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
    gamma:
        gamma value for gamma correction of mixed colors; larger values make mixed colors appear brighter (default: 4.5)

    Returns
    -------
    A dictionary with selected dynamic ranges for plotting (e.g. for adding colorbars).
    """

    # channel1 info
    channels = [channel1]
    channel_params = {
        channel1: {
            "color": channel1_color,
            "vmin": channel1_vmin,
            "vmax": channel1_vmax,
        }
    }

    # add channel info from kwargs
    channel_kwargs = {k: v for k, v in kwargs.items() if re.match("^channel[0-9]+$", k)}
    channels.extend(channel_kwargs.values())
    channel_params.update(
        {
            {
                "color": kwargs[f"{name}_color"],  # color required
                "vmin": kwargs.get(f"{name}_vmin"),
                "vmax": kwargs.get(f"{name}_vmax"),
            }
            for name, cnl in channel_kwargs.items()
        }
    )

    # checks
    for cnl in channels:
        if cnl not in adata.var_names:
            raise ValueError(f"gene '{cnl}' not in 'adata.var_names'")

    # get expression vectors
    exp_vectors = []
    if knn_smooth:
        for cnl in channels:
            exp_vectors[cnl] = np.array(
                graph_smooth(adata[:, cnl].X, adata.obsp["connectivities"])
            ).flatten()
    else:
        for cnl in channels:
            exp_vectors[cnl] = get_array(adata, cnl)

    for cnl in channels:
        exp_vec = exp_vectors[cnl]
        if not channel_params[cnl]["vmin"] and auto_range_quantiles:
            channel_params[cnl]["vmin"] = np.quantile(
                exp_vec[exp_vec > 0], q=[auto_range_quantiles[0]]
            )
        if not channel_params[cnl]["vmax"] and auto_range_quantiles:
            channel_params[cnl]["vmax"] = np.quantile(
                exp_vec[exp_vec > 0], q=[auto_range_quantiles[1]]
            )

    # assign colors
    rgb_values = [
        mix_colors(
            [
                np.array(channel_params[cnl]["color"]) * (x if not np.isnan(x) else 0)
                for cnl, x in zip(channels, exp_vals_scaled)
            ],
            gamma=gamma,
        )
        for exp_vals_scaled in zip(
            *[
                scale(
                    exp_vectors[cnl],
                    vmin=channel_params[cnl]["vmin"],
                    vmax=channel_params[cnl]["vmax"],
                )
                for cnl in channels
            ]
        )
    ]

    # store information in anndata
    adata.obs["pseudo_RNAscope"] = adata.obs_names.astype("category")
    adata.uns["pseudo_RNAscope_colors"] = [rgb2hex(x) for x in rgb_values]
    adata.obs["pseudo_RNAscope_alpha"] = [
        y for y in scale(np.array([np.mean(x) for x in rgb_values]), max_val=1)
    ]

    return channel_params
