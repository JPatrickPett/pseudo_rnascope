"""
simple package for plotting multiple genes with different colors and additive colormixing in scanpy or squidpy.
"""

__version__ = "0.0.8"

import re
import numpy as np
import scipy as sp
import matplotlib as mpl


def get_array(adata, annot):
    """
    Extract expression values or annotation for a specified gene or annotation column in an AnnData object.

    Parameters
    ----------
    adata : anndata.AnnData
        An AnnData object containing single-cell or single-nucleus gene expression data.
    annot : str
        The symbol of the gene (in `adata.var_names`) for which expression values are to be extracted or a column in
        `adata.obs` for which values are to be extracted.

    Returns
    -------
    np.ndarray
        A 1-dimensional NumPy array containing the values for `annot` across all cells or nuclei.

    Raises
    ------
    ValueError
        If the provided gene symbol or annotation column is not present in the AnnData object.
    """
    if annot in adata.var_names:
        val_mat = adata.X[:, adata.var_names == annot]
    elif annot in adata.obs.columns:
        val_mat = adata.obs[annot].values
    else:
        raise ValueError(
            f"Gene symbol '{annot}' not found in the provided AnnData object."
        )

    if sp.sparse.issparse(val_mat):
        val_mat = val_mat.todense()

    return np.array(val_mat).flatten()


def scale(x, max_val=255, vmin=None, vmax=None):
    """
    Scale the values of an array to a specified range.

    Parameters
    ----------
    x : np.ndarray
        The input array to be scaled.
    max_val : float, optional
        The maximum value of the scaled range. Default is 255.
    vmin : float or None, optional
        If specified, values below this threshold will be set to NaN before scaling. Default is None.
    vmax : float or None, optional
        If specified, values above this threshold will be set to `vmax` before scaling. Default is None.

    Returns
    -------
    np.ndarray
        The scaled array.

    Notes
    -----
    The function scales the input array `x` to the range [0, `max_val`] based on the specified minimum (`vmin`)
    and maximum (`vmax`) values. If `vmin` or `vmax` is not provided, the minimum and maximum values of the
    input array are used for scaling.

    NaN values are ignored during scaling.

    Examples
    --------
    >>> scaled_array = scale(input_array, max_val=1, vmin=0.5, vmax=2.5)

    """
    x = x.copy().astype(float)

    if vmin:
        x[x < vmin] = np.nan
    if vmax:
        x[x > vmax] = vmax

    return max_val * (x - np.nanmin(x)) / (np.nanmax(x) - np.nanmin(x))


def mix_colors(colors, gamma=10):
    """
    Mix multiple colors by applying gamma correction.

    Parameters
    ----------
    colors : list of tuples
        A list of tuples representing RGB colors. Each tuple should contain three values for the red, green, and blue components.
    gamma : float, optional
        Gamma correction parameter. Must be greater than 0. Default is 10.

    Returns
    -------
    tuple
        A tuple representing the mixed color after applying gamma correction.

    Raises
    ------
    AssertionError
        If the gamma value is not greater than 0.

    Examples
    --------
    >>> mix_colors([(255, 0, 0), (0, 255, 0), (0, 0, 255)])
    (127, 127, 127)

    >>> mix_colors([(255, 0, 0), (0, 255, 0), (0, 0, 255)], gamma=2.2)
    (169, 169, 169)

    """
    assert gamma > 0
    return (
        np.power(sum([c**gamma for c in colors]) / len(colors), 1 / gamma)
        .round()
        .astype(int)
    )


def rgb2hex(rgb_tuple):
    """
    Convert an RGB color represented as a tuple to its hexadecimal representation.

    Parameters
    ----------
    rgb_tuple : tuple
        A tuple representing the RGB color with three values for red, green, and blue components.

    Returns
    -------
    str
        A string representing the hexadecimal color code.

    Examples
    --------
    >>> rgb2hex((255, 0, 128))
    '#ff0080'

    """
    r, g, b = np.round(rgb_tuple).astype(int)
    return "#{:02x}{:02x}{:02x}".format(r, g, b)


def hex2rgb(hexcode):
    """
    Convert a hexadecimal color code to an RGB color represented as a tuple.

    Parameters
    ----------
    hexcode : str
        A string representing the hexadecimal color code.

    Returns
    -------
    tuple
        A tuple containing three values for the red, green, and blue components of the color.

    Examples
    --------
    >>> hex2rgb('#ff0080')
    (255, 0, 128)

    """
    return tuple(map(ord, hexcode[1:].decode("hex")))


def convert_base(number, from_base, to_base):
    """
    Convert a number from one base to another.

    Args:
    - number (List[int]): The list representing the number to be converted.
    - from_base (int): The base of the input number.
    - to_base (int): The desired base for the output.

    Returns:
    - List[int]: The converted number in the specified base.
    """
    # Convert the input number to base 10
    decimal_number = 0
    for digit in number:
        if digit < 0 or digit >= from_base:
            raise ValueError("Invalid digit in the input number")

        decimal_number = decimal_number * from_base + digit

    # Convert the base 10 number to the desired base
    converted_number = []
    while decimal_number > 0:
        remainder = decimal_number % to_base
        converted_number.insert(0, remainder)
        decimal_number //= to_base

    return converted_number if converted_number else [0]


def encode_rgb(rgb_list):
    """
    Encode a list of RGB color values into a list of normalized decimal numbers.

    Parameters
    ----------
    rgb_list : list of tuples
        A list of tuples representing RGB color values. Each tuple should contain three integers for the red, green, and blue components.

    Returns
    -------
    list
        A list of normalized decimal numbers representing the encoded RGB colors.

    Notes
    -----
    This function encodes RGB color values into normalized decimal numbers. Each RGB color value is converted to a decimal
    number using a base conversion approach. The resulting decimal numbers are then normalized to the range [0, 1].

    Examples
    --------
    >>> encode_rgb([(255, 0, 0), (0, 255, 0), (0, 0, 255)])
    [1.5199185323666652e-05, 0.003890991442858663, 0.9960938093718177]

    """
    decimal_numbers = [convert_base(reversed(x), 256, 10) for x in rgb_list]
    decimal_numbers_simple = [int("".join(str(i) for i in x)) for x in decimal_numbers]
    return [x / 16777215 for x in decimal_numbers_simple]


def graph_smooth(expression_matrix, neighbor_matrix):
    """
    Smooth gene expression values over a k-nearest neighbors (knn) graph.

    Parameters
    ----------
    expression_matrix : np.ndarray
        The gene expression matrix where rows represent cells and columns represent genes.
    neighbor_matrix : np.ndarray
        The k-nearest neighbors graph matrix indicating the connectivity between cells.

    Returns
    -------
    np.ndarray
        A smoothed gene expression matrix based on the k-nearest neighbors graph.

    Notes
    -----
    This function performs a smoothing operation over the gene expression matrix using information from the
    k-nearest neighbors graph. It calculates the weighted sum of expression values for each cell and its neighbors,
    then normalizes by the number of neighbors each cell has. The result is a smoothed gene expression matrix.

    Examples
    --------
    >>> smoothed_expression = graph_smooth(expression_matrix, neighbor_matrix)

    """
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
    adata.obs["pseudo_RNAscope_alt"] = encode_rgb(
        rgb_values
    )  # encode rgb-triplets as decimal float number between 0 and 1
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
        "cmap": RGBcmap("decode_RGB"),
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
            if n_arr.size == 1:
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
            color += [0] * (4 - len(color))
        rgba = np.array(out, dtype=float)
        rgba /= 255

        if alpha is not None:
            alpha = np.clip(alpha, 0, 1)
            if alpha.shape not in [(), xa.shape]:
                raise ValueError(
                    f"alpha is array-like but its shape {alpha.shape} does "
                    f"not match that of X {xa.shape}"
                )
            rgba[..., -1] = alpha

        if bytes:
            rgba *= 255

        if not np.iterable(X):
            rgba = tuple(rgba)
        return rgba
