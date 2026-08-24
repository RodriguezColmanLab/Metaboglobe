# Metaboglobe

**Metaboglobe** is a Python package to visualize the KEGG metabolic pathway maps with your data. It can be used to visualize metabolic fluxes, enzyme transcript levels, metabolites, or any other values that can be mapped to the KEGG pathway maps. Built on Matplotlib, it is designed to be flexible and easy to use, and it can be used in a Jupyter notebook or in a Python script.

![Metaboglobe example](docs/images/example.png)

## Why Metaboglobe?
Strangely enough, there are only a handful of software packages to visualize metabolic pathways ([Metabolic Atlas](https://metabolicatlas.org/), [Pathview](https://bioconductor.org/packages/release/bioc/html/pathview.html). They may have limitations in the output formats (only low-res PNGs for example), or require uploading your data to a web server. Metaboglobe is designed to work offline, give you the flexibility to customize the plots, and the output can be in any format that Matplotlib supports (PNG, PDF, SVG, etc.).

In addition, Metaboglobe also has built-in support for visualizing Compass results on KEGG pathway maps. [Compass](https://compass-wagnerlab.readthedocs.io/en/latest/) is a software package by Wagner et al. to predict metabolic fluxes from single-cell RNA-seq data. With Metaboglobe, you can easily visualize these predicted fluxes on the KEGG pathway maps, which can help you interpret the results and generate hypotheses about metabolic changes in your single-cell transcriptomic data.

## Installation
Make sure you have Python installed. It might be useful to install Metaboglobe into a separate environment (Venv/Conda/etc.). Then, install with PIP:

    pip install git+https://github.com/RodriguezColmanLab/Metaboglobe.git@main

## Usage

Documentation is available at [rodriguezcolmanlab.github.io/Metaboglobe](https://rodriguezcolmanlab.github.io/Metaboglobe).
