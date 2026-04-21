
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

