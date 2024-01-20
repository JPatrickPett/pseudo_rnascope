# pseudo-RNAscope

`pseudo_rnascope` is ...

## TLDR

Run like this...

```python
channel1 = "FGF2"
channel2 = "FGFR2"

import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
from pseudo_rnascope import add_pseudo_rna_scope

ranges = add_pseudo_rna_scope(
    adata,
    channel1 = channel1,
    channel2 = channel2,
    auto_range_quantiles = (0.2, 0.9),
    knn_smooth = True,
)

sc.set_figure_params(figsize=[16,16], dpi=75)

sc.pl.spatial(
    adata, 
    img_key="hires", 
    color='pseudo_RNAscope', 
    size=1.5,
    legend_loc=None,
    alpha = adata.obs["pseudo_RNAscope_alpha"],
    show=False,
)

plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(ranges["channel2_vmin"], ranges["channel2_vmax"]), cmap='Reds'),
             orientation='vertical', label=channel2, extend='both', shrink=0.5, pad=-0.04)

plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(ranges["channel1_vmin"], ranges["channel1_vmax"]), cmap='Greens'),
             orientation='vertical', label=channel1, extend='both', shrink=0.5, pad=0.03)
```

## System requirements

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
Currently only `numpy` from the python scientific stack is required as a dependency, however, it only works together with an `anndata` object from the `scanpy` package as input.
`pseudo_rnascope` and all dependencies can be installed via `pip` (see below).

## Installation

*Optional: create and activate a new conda environment (with python<3.12):*
```bash
mamba create -n pseudo_rnascope "python>3.9"
mamba activate pseudo_rnascope
```

### Install with pip

**from PyPI**

*will be added*

**from github**

```bash
pip install git+https://github.com/JPatrickPett/pseudo_rnascope.git
```

*(installation time: around 1 min)*

## Usage and Documentation

Please refer to the [demo notebook](notebooks/demo.ipynb). There's a function docstring in the source code, which will be rendered on ReadTheDocs once the package goes live.

(*demo running time: around 2 min*)

## Citation

`pseudo_rnascope` is part of the forthcoming manuscript "A multiomic atlas of human early skeletal development" by To, Fei, Pett et al. Stay tuned for details!

