"""Approximate missing features from higher dimensionality data neighbours"""
__version__ = "0.0.5"

import re
import numpy as np
import scipy as sp
import matplotlib as mpl


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
    return np.power(sum([c**gamma for c in colors]) / len(colors), 1 / gamma).round().astype(int)


def rgb2hex(rgb_tuple):
    r, g, b = np.round(rgb_tuple).astype(int)
    return "#{:02x}{:02x}{:02x}".format(r, g, b)


def hex2rgb(hexcode):
    return tuple(map(ord, hexcode[1:].decode("hex")))


def to_decimal(n_list, b):
    """convert each number in n_list of base b to decimal (non-dec numbers represented as lists)"""
    out = []
    for n in n_list:
        result = 0
        for i, num in enumerate(n):
            result += num * b ** i
        out.append(result)
    return out


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
        **{
            cnl: {
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
    exp_vectors = {}
    if knn_smooth:
        for cnl in channels:
            smooth_mat = graph_smooth(adata[:, cnl].X, adata.obsp["connectivities"])
            if sp.sparse.issparse(smooth_mat):
                smooth_mat = smooth_mat.todense()
            exp_vectors[cnl] = np.array(smooth_mat).flatten()
    else:
        for cnl in channels:
            exp_vectors[cnl] = get_array(adata, cnl)

    for cnl in channels:
        exp_vec = exp_vectors[cnl]
        if channel_params[cnl]["vmin"] is None and auto_range_quantiles:
            channel_params[cnl]["vmin"] = np.quantile(
                exp_vec[exp_vec > 0], q=[auto_range_quantiles[0]]
            )
        if channel_params[cnl]["vmax"] is None and auto_range_quantiles:
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
    adata.obs["pseudo_RNAscope_alt"] = [x / 16777215 for x in to_decimal(rgb_values, 256)]  # encode rgb-triplets as decimal float number between 0 and 1
    adata.obs["pseudo_RNAscope_alpha"] = scale(
        np.array([np.max(x) for x in rgb_values]), max_val=1
    )
    
    adata.uns["pseudo_RNAscope_colors"] = [rgb2hex(x) for x in rgb_values]
    
    adata.uns["pseudo_RNAscope"] = {
        "auto_range_quantiles": auto_range_quantiles,
        "knn_smooth": knn_smooth,
        "gamma": gamma,
        "channels": channels,
        "channel_params": channel_params,
        "rgb_values": rgb_values,
        "cmap": RGBcmap('decode_RGB'),
    }

    return channel_params


class RGBcmap(mpl.colors.Colormap):
    """
    matplotlib colormap to decode rgb-triplets saved as decimal [0,1] floating point numbers
    """
    def __init__(self, *args, **kwargs):
        super(RGBcmap, self).__init__(*args, **kwargs)

    def __call__(self, X, alpha=None, bytes=False):
        # TODO: optimize this
        xa = np.array(X, copy=True)
            
        def dec_to_base(n_arr, b):
            out = []
            if n_arr.size==1:
                n_arr = [n_arr]
            for n in n_arr:
                if n == 0:
                    out.append([0])
                    continue
                digits = []
                while n:
                    digits.append(int(n % b))
                    n //= b
                out.append(digits[::-1])
            return out
            
        out = dec_to_base(xa * 16777215, 256)
        for color in out:
            color.reverse()
            color += [0]*(4-len(color))
        rgba = np.array(out, dtype=float)
        rgba /= 255

        if alpha is not None:
            alpha = np.clip(alpha, 0, 1)
            if alpha.shape not in [(), xa.shape]:
                raise ValueError(
                    f"alpha is array-like but its shape {alpha.shape} does "
                    f"not match that of X {xa.shape}")
            rgba[..., -1] = alpha

        if bytes:
            rgba *= 255
        
        if not np.iterable(X):
            rgba = tuple(rgba)
        return rgba





