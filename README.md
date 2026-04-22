# Metaboglobe

**Metaboglobe** is a Python package to visualize the KEGG metabolic pathway maps with your data. It can be used to visualize metabolic fluxes, enzyme transcript levels, metabolites, or any other values that can be mapped to the KEGG pathway maps. Built on Matplotlib, it is designed to be flexible and easy to use, and it can be used in a Jupyter notebook or in a Python script.

## Why Metaboglobe?
Strangely enough, there are only a handful of software packages to visualize metabolic pathways [Metabolic Atlas, Pathview]. They may have limitations in the output formats (only low-res PNGs for example), or require uploading your data to a web server. Metaboglobe is designed to work offline, give you the flexibility to customize the plots, and the output can be in any format that Matplotlib supports (PNG, PDF, SVG, etc.).

In addition, Metaboglobe also has built-in support for visualizing Compass results on KEGG pathway maps. Compass is a software package to predict metabolic fluxes from single-cell RNA-seq data [Wagner et al.]. With Metaboglobe, you can easily visualize these predicted fluxes on the KEGG pathway maps, which can help you interpret the results and generate hypotheses about metabolic changes in your single-cell transcriptomic data.

## How to plot a KEGG pathway map with MetaboGlobe
First, download a KEGG pathway map in KGML format from the KEGG website. For example, the Glycolysis pathway map for humans can be downloaded from https://rest.kegg.jp/get/hsa00010/kgml . Then, you can load it using the `load_kegg_map` function:

```python
import metaboglobe.kegg_pathway
# load in the KEGG map for glycolysis
kegg_map = metaboglobe.kegg_pathway.load_kegg_map("path/to/glycolysis_map.xml")
``` 

Plotting the map itself can be done as follows:

```python
import metaboglobe.kegg_plotting
from matplotlib import pyplot as plt

kegg_map = ...  # See above

fig, ax = plt.subplots(figsize=(10, 10))
metaboglobe.kegg_plotting.plot_double_arrows(ax, kegg_map)
plt.show()
```


## How to use with a data table
Let's say you have a pandas DataFrame like this:

```python
import pandas

data_frame = pandas.DataFrame({
    "reactions": ["1.00 * Citrate [c] --> 1.00 * Isocitrate [c] ACO1; IREB2", 
                "1.00 * Isocitrate [c] --> 1.00 * alpha-Ketoglutarate [c] ACO1; IREB2",
                "1.00 * alpha-Ketoglutarate [c] --> 1.00 * Succinyl-CoA [c] OGDH; DLST"],
    "flux_values": [0.5, 1.0, 0.2]
})
```

Then you can add the flux values to the KEGG map using the `insert_values_in_map` function:

```python
import metaboglobe.tabular_data

kegg_map = ...  # See above
data_frame = ...  # See above

metaboglobe.tabular_data.insert_values_in_map(kegg_map, data_frame, reaction_col="reactions", value_col="flux_values")
```

And then plot like above.

Normally your data would not be hardcoded, but come from a data file such as a CSV file. If you're unfamiliar with Pandas, it's easy to load from a CSV file. For example, if you have a CSV file with the same columns as above, you can load it like this:

```python
import pandas
data_frame = pandas.read_csv("path/to/your/data.csv")
```

## How to use with Compass data
[Compass](https://github.com/wagnerlab-berkeley/Compass) is a tool to predict metabolic fluxes from single-cell RNA-seq data. This package provides a convenient way to visualize Compass results on the KEGG metabolic pathway maps.

First, load in the scRNAseq data into an AnnData object, and make sure you have ran Compass on the package. Then, load in your results:

```python
import metaboglobe.compass_data

adata = ...  # load in your AnnData object
metaboglobe.compass_data.add_compass_output(adata, "path/to/compass output/CENTRAL_CARBON_ENERGY", obsm_key="CENTRAL_CARBON_ENERGY")
```

This folder is expected to have a `reactions.tsv` file and a `model.json.gz` file in it. By default, the reaction fluxes will be stored under `adata.obsm[folder_name]`, with `folder_name = CENTRAL_CARBON_ENERGY` in this example. You can specify a different key using the `obsm_key` argument.

Second, load in the KEGG pathway map you want to plot (see above), and apply the coloring to the map:

```python
import metaboglobe.compass_data

kegg_map = ...  # load in the KEGG map for glycolysis
adata = ...  # load in your AnnData object with Compass results

model = metaboglobe.compass_data.load_compass_model("path/to/compass output/CENTRAL_CARBON_ENERGY")
comparison = metaboglobe.compass_data.setup_comparison_to_single_cells(adata, model, groupby="condition",
                                                                       obsm_key="CENTRAL_CARBON_ENERGY",
                                                                       min_percentile=20, max_percentile=80)
comparison.insert_values_in_map(kegg_map, group="condition_a")

# Plot the kegg_map as shown above. It should now be colored according to the Compass fluxes for "condition_a".
```

