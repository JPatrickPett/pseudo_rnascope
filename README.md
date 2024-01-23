# pseudo-RNAscope

`pseudo_rnascope` is a simple package to plot gene expression of multiple genes as different colors in the same plot using any scanpy plotting function.
The original use-case is to mimic RNAscope with visium data.

![pseudo_rnascope_limb](notebooks/pseudo_rnascope_limb.png?raw=true "Title")
*visualisation of a developing limb with data from: B. Zhang et al. (2023) ‘A human embryonic limb cell atlas resolved in space and time.’ Nature 2023. DOI: 10.1038/s41586-023-00000-0.* 

## TLDR

Run like this...

To plot the expression overlap of two genes, e.g. *FGF2* and *FGFR2*, in spatial locations of visium data

```python
channel1 = "FGF2"
channel2 = "FGFR2"
```

run `add_pseudo_rna_scope` on an existing `anndata` object:

```python
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
from pseudo_rnascope import add_pseudo_rna_scope

ranges = add_pseudo_rna_scope(
    adata,

    channel1 = channel1,
    channel1_color = (1,0,0),  # red (RGB)

    channel2 = channel2,
    channel2_color = (0,1,0),  # green (RGB)

    auto_range_quantiles = (0.2, 0.9),
    knn_smooth = True,
    gamma = 10,
)
```

This will use red for channel1 (*FGF2*) and green for channel2 (*FGFR2*), hence overlaps will appear yellow.

- Dynamic ranges can be set for both channels explicitly, or automatically via upper and lower gene expression quantiles with `auto_range_quantiles`.
- Neareast-neighbor smoothing of expression values can be added to mitigate data sparsity with `knn_smooth=True`.
- Larger `gamma` values make mixed colors appear brighter, while `gamma=1` corresponds to linear color mixing.

Scanpy plotting functions can then be used by specifying the added `adata.obs` column `"pseudo_RNAscope"`.

```python
sc.set_figure_params(figsize=[16,16], dpi=75)

sc.pl.spatial(
    adata, 
    img_key="hires", 
    color='pseudo_RNAscope', 
    size=1.5,
    legend_loc=None,  # turn off legend (necessary)
    alpha = adata.obs["pseudo_RNAscope_alpha"],
    show=False,
)

plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(ranges[channel2]["vmin"], ranges[channel2]["vmax"]), cmap='Greens'),
             orientation='vertical', label=channel2, extend='both', shrink=0.5, pad=-0.04)

plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(ranges[channel1]["vmin"], ranges[channel1]["vmax"]), cmap='Reds'),
             orientation='vertical', label=channel1, extend='both', shrink=0.5, pad=0.03)
```

Corresponding colors for all spots are saved in `adata.uns["pseudo_RNAscope_colors"]` and another `anndata.obs` column called `"pseudo_RNAscope_alpha"` allows to set alpha values, making spots with low color intensity transparent.
Computed colors are also encoded as decimal numbers between 0 and 1, and stored in `adata.obs["pseudo_RNAscope_alt"]`. This can be used as an **alternative** option for plotting.

<details>
<summary><b>alternative</b></summary>

Use values stored in `adata.obs["pseudo_RNAscope_alt"]` and pass `adata.uns["pseudo_RNAscope"]["cmap"]` as a colormap to decode them:

```python
sc.set_figure_params(figsize=[16,16], dpi=75)

sc.pl.spatial(
    adata, 
    img_key="hires", 
    color='pseudo_RNAscope_alt', 
    size=1.5,
    cmap = adata.uns['pseudo_RNAscope']['cmap'],
    vmin=0, # use all values (necessary)
    vmax=1, # use all values (necessary)
    colorbar_loc=None,  # turn off native colorbar, we need separate ones per channel
    alpha = adata.obs.sort_values("pseudo_RNAscope_alt")["pseudo_RNAscope_alpha"],  # since values will be sorted before plotting
    show=False,
)

plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(ranges[channel2]["vmin"], ranges[channel2]["vmax"]), cmap='Greens'),
             orientation='vertical', label=channel2, extend='both', shrink=0.5, pad=-0.04)

plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(ranges[channel1]["vmin"], ranges[channel1]["vmax"]), cmap='Reds'),
             orientation='vertical', label=channel1, extend='both', shrink=0.5, pad=0.03)
```

</details>

## System requirements

<details>
<summary><b>show requirements</b></summary>

### Hardware requirements

`pseudo_rnascope` can run on a standard computer with enough RAM to hold the used datasets in memory.

### Software requirements

**OS requirements**

The package has been tested on:

- macOS Monterey (12.6.7)
- Linux: Ubuntu 18.04.6 bionic

**Python requirements**

A python version `>=3.0` is required. 
Various python libraries are used, listed in `pyproject.toml`. 
Currently only the python scientific stack (`numpy`, `scipy`) is required as a dependency, however, it only works together with an `anndata` object from the `scanpy` package as input and for downstream plotting.
`pseudo_rnascope` and all dependencies can be installed via `pip` (see below).

</details>

## Installation

### Install with pip

**from github**

```bash
pip install git+https://github.com/JPatrickPett/pseudo_rnascope.git
```

**from PyPI**

*will be added*

*(installation time: <1 min)*

## Usage and Documentation

Please refer to the [demo notebook](notebooks/demo.ipynb). There's a function docstring in the source code, which will be rendered on ReadTheDocs once the package goes live.

(*demo running time: around 2 min*)

## Citation

`pseudo_rnascope` is part of the forthcoming manuscript "A multiomic atlas of human early skeletal development" by To, Fei, Pett et al. Stay tuned for details!

