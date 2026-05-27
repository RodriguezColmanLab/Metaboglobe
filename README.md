# Metaboglobe

**Metaboglobe** is a Python package to visualize the KEGG metabolic pathway maps with your data. It can be used to visualize metabolic fluxes, enzyme transcript levels, metabolites, or any other values that can be mapped to the KEGG pathway maps. Built on Matplotlib, it is designed to be flexible and easy to use, and it can be used in a Jupyter notebook or in a Python script.

## Why Metaboglobe?
Strangely enough, there are only a handful of software packages to visualize metabolic pathways [Metabolic Atlas, Pathview]. They may have limitations in the output formats (only low-res PNGs for example), or require uploading your data to a web server. Metaboglobe is designed to work offline, give you the flexibility to customize the plots, and the output can be in any format that Matplotlib supports (PNG, PDF, SVG, etc.).

In addition, Metaboglobe also has built-in support for visualizing Compass results on KEGG pathway maps. Compass is a software package to predict metabolic fluxes from single-cell RNA-seq data [Wagner et al.]. With Metaboglobe, you can easily visualize these predicted fluxes on the KEGG pathway maps, which can help you interpret the results and generate hypotheses about metabolic changes in your single-cell transcriptomic data.

## Installation
Make sure you have Python installed. It might be useful to install Metaboglobe into a separate environment (Venv/Conda/etc.). Then, install with PIP:

    pip install git+https://github.com/RodriguezColmanLab/Metaboglobe.git@main

## How to plot a KEGG pathway map with MetaboGlobe
First, download a KEGG pathway map in KGML format from the KEGG website. For example, the Glycolysis pathway map for humans can be downloaded from https://rest.kegg.jp/get/hsa00010/kgml . Then, you can load it using the `load_kegg_map` function:

```python
import metaboglobe.kegg_pathway
# load in the KEGG map for glycolysis
kegg_map = metaboglobe.kegg_pathway.load_kegg_map("path/to/glycolysis_map.xml")
``` 

Plotting the map itself can be done as follows:

```python
import metaboglobe.plotting
from matplotlib import pyplot as plt

kegg_map = ...  # See above

fig, ax = plt.subplots(figsize=(10, 10))
metaboglobe.plotting.plot_kegg(ax, kegg_map)
plt.show()
```

You can then color individual reactions or metabolites like this:

```python
from metaboglobe.kegg_species import KeggReactionId, KeggCompoundId
kegg_map = ...  # See above

# Sets the score of KEGG reaction R00703 (lactate -> pyruvate) in the forward direction to 0.5
kegg_map.set_reaction_score(KeggReactionId.create_from_id("R00703").forwards(), 0.5)

# Sets the score of KEGG compound C00022 (pyruvate) to 0.8
kegg_map.set_compound_score(KeggCompoundId.create_from_id("C00022"), 0.8)
```

Since it's not very convenient to set the scores one by one using the KEGG IDs, there are also helper functions to set scores in bulk using reaction equations or compound names. See the next section for how to use these functions.

## How to use with a data table
Let's say you have a pandas DataFrame like this:

```python
import pandas

data_frame = pandas.DataFrame({
    "reactions": ["1.00 * Citrate [c] --> 1.00 * Isocitrate [c]", 
                "1.00 * Isocitrate [c] --> 1.00 * alpha-Ketoglutarate [c]",
                "1.00 * alpha-Ketoglutarate [c] --> 1.00 * Succinyl-CoA [c]"],
    "flux_values": [0.5, 1.0, 0.2]
})
```

Then you can add the flux values to the KEGG map using the `insert_values_in_map` function:

```python
import metaboglobe.tabular_data

kegg_map = ...  # See above
data_frame = ...  # See above

metaboglobe.tabular_data.insert_values_in_map(kegg_map, data_frame, reaction_identifier_col="reactions",
                                              value_col="flux_values")
```

And then plot like above. Instead of reaction equations, you can also use the KEGG reaction IDs (e.g. rn:R00001) in the `reactions` column, ecNums (e.g. ec:1.1.1.1) or IDs from the standard-GEM initiative (e.g. MAR00001). If you use reaction equations, you can also specify enzymes to narrow down possible matches. See the documentation of the function for more details.

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
metaboglobe.compass_data.add_compass_output(adata, "path/to/compass output/")
```

This folder is expected to have a `reactions.tsv` file and a `model.json.gz` file in it. Alternatively, it can contain subfolders with those files, which you will get if you run Module-Compass. By default, the reaction fluxes will be stored under `adata.obsm["compass"]`. You can specify a different key using the `obsm_key` argument.

Second, load in the KEGG pathway map you want to plot (see above), and apply the coloring to the map:

```python
import metaboglobe.compass_data

kegg_map = ...  # load in the KEGG map for glycolysis
adata = ...  # load in your AnnData object with Compass results

model = metaboglobe.compass_data.load_compass_model("path/to/compass output")
comparison = metaboglobe.compass_data.setup_comparison_to_single_cells(adata, model, groupby="condition",
                                                                       min_percentile=20, max_percentile=80)
comparison.insert_values_in_map(kegg_map, group="condition_a")

# Plot the kegg_map as shown above. It should now be colored according to the Compass fluxes for "condition_a".
```

Reactions are first matched to the pathway using the IDs stored in the model, either EC numbers (for the RECON model) or GEM reaction IDs (for the Human1 or similar model), which are translated to KEGG Reaction IDs. However, in some cases this matching fails. For example, the Human1 model models the reaction [R00771](https://www.kegg.jp/entry/R00771), but the glycolysis pathway of KEGG shows the reaction [R13199](https://www.kegg.jp/entry/R13199). These reactions are actually the same, and differ only on the level at which stereoisomers are specified. In these cases matching is done by comparing the product names, reactant names and enzyme gene names. If all of these have at least one match, then the reaction is considered a match.

## How the plotting works
The KEGG map has information on which metabolites, pathway references and enzymes are present on which locations, and how they are connected. What it lacks, is the path of the arrows, or the location of the metabolite labels. To still be able to make a plot, we have created several layout algorithms.

For most reactions, we can draw a U-shaped path, where the metabolites are on both ends of the U, and the enzyme label is at the bottom of the U. The lines of the U are axis-aligned (so along the X or Y axis, no diagonals), and have a rounded corner. The algorithm works by first drawing a bounding box around the two metabolites and the enzyme. If the enzyme is on one of the four sides of the bounding box, then we can draw a U-shaped path. However, if that's not the case, we just draw a straight (non-axis-aligned) line from one metabolite to the enzyme, and then from the enzyme to the other metabolite, with a rounded corner at the metabolite. This is not ideal, but it still allows us to visualize the pathway.

For the metabolite labels, we set up a collision map for the entire plot, blocking the areas where lines or labels are already drawn. Then, we sort all metabolite labels by the text width (largest first). Then, in this order we try to place it in one of the 9 possible locations on or around the metabolite (top, top-right, right, bottom-right, bottom, bottom-left, left, top-left, center). Each location gets a score based on the location itself (right is preferred over bottom-right) and how much it overlaps with blocked areas in the collision map. The location with the best score is chosen, and the collision map is updated to block the area of the label.

## Metabolites
[gem_reactions.csv](src/metaboglobe/data/gem_reactions.tsv) was downloaded from https://github.com/SysBioChalmers/Human-GEM/blob/635f533152dc5f7290ce04d12700eaa882273c3e/model/reactions.tsv .  (Robinson JL, et al. An atlas of human metabolism. Sci. Signal. 13, eaaz1482 (2020). [doi:10.1126/scisignal.aaz1482](https://doi.org/10.1126/scisignal.aaz1482) )

