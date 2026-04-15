from abc import ABC, abstractmethod

import numpy
import os
import gzip
import json

import pandas
import scanpy.get
from anndata import AnnData
from pandas import DataFrame
from typing import NamedTuple

from metaboglobe._util import optimize_for_matching
from metaboglobe.kegg_pathway import KeggMap


_DEFAULT_OBSM_KEY = "compass"


class CompassReaction(NamedTuple):
    reactant_names: list[str]
    product_names: list[str]


class CompassModel:
    """A class for mapping the reaction and metabolite IDs in the Compass output to the names."""
    _species_names_by_id: dict[str, str]
    _reactions_by_id: dict[str, CompassReaction]

    def __init__(self):
        self._species_names_by_id = dict()
        self._reactions_by_id = dict()

    def add_species(self, id: str, name: str):
        self._species_names_by_id[id] = name

    def species(self, id: str) -> str:
        """Gets the species with the given ID. For example, for "MAM01840e" it could return "fructose". Returns
        the ID itself if not found."""
        return optimize_for_matching(self._species_names_by_id.get(id, id))

    def add_reaction(self, id: str, reactant_names: list[str], product_names: list[str]):
        """Adds the reaction with the given ID, reactants and products. The ID will be something like "MAR04363_pos".
        The reactants and products are expected to be list of compound names, like ["fructose", "ATP"] or
         ["fructose-6-phosphate", "ADP"]."""
        self._reactions_by_id[id] = CompassReaction(reactant_names=reactant_names, product_names=product_names)

    def reaction(self, id: str) -> CompassReaction | None:
        """Gets the reaction with the given ID. For example, for "MAR04363_pos". Returns None if not found."""
        return self._reactions_by_id.get(id, None)



class CompassComparison(ABC):
    """For plotting a comparison on a KEGG map."""

    @abstractmethod
    def apply_color(self, kegg_map: KeggMap, group: str):
        """Applies the comparison to the given KEGG map, by coloring the reactions in the map according to the values
        of the given group."""
        return NotImplemented



class _ComparisonToSingleCellValues(CompassComparison):
    _model: CompassModel

    _groupby: str
    _min_percentile: float
    _max_percentile: float
    _obsm_key: str

    _min_values_by_reaction: numpy.ndarray
    _max_values_by_reaction: numpy.ndarray
    _reactions_aggregated: AnnData

    def __init__(self, *, adata: AnnData, model: CompassModel, groupby: str, obsm_key: str, min_percentile: float, max_percentile: float):
        self._model = model

        self._groupby = groupby
        self._min_percentile = min_percentile
        self._max_percentile = max_percentile
        self._obsm_key = obsm_key

        # Calculate scaling values per column
        reaction_scores = adata.obsm[obsm_key]
        if min_percentile >= max_percentile:
            raise ValueError(f"min_percentile '{min_percentile}' must be less than max_percentile '{max_percentile}'")
        self._min_values_by_reaction = numpy.percentile(reaction_scores, min_percentile, axis=0)
        self._max_values_by_reaction = numpy.percentile(reaction_scores, max_percentile, axis=0)

        # Calculate reaction unscaled means
        self._reactions_aggregated = scanpy.get.aggregate(adata, by=self._groupby, func="mean", obsm=self._obsm_key)

    def apply_color(self, kegg_map: KeggMap, group: str):
        reaction_ids = self._reactions_aggregated.var_names

        group_index = self._reactions_aggregated.obs_names.get_loc(group)
        reaction_values = self._reactions_aggregated.layers["mean"][group_index]

        for reaction_id, reaction_value, min_value, max_value in zip(reaction_ids, reaction_values,
                                                                     self._min_values_by_reaction, self._max_values_by_reaction):
            reaction = self._model.reaction(reaction_id)

            match = kegg_map.match_reaction(reaction.reactant_names, reaction.product_names)
            if match is None:
                continue

            scaled_reaction = (reaction_value - min_value) / (max_value - min_value)
            if scaled_reaction < 0:
                scaled_reaction = 0
            elif scaled_reaction > 1:
                scaled_reaction = 1
            kegg_map.set_reaction_score(match, scaled_reaction)


def setup_comparison_to_single_cells(adata: AnnData, model: CompassModel, *, groupby: str, obsm_key: str = _DEFAULT_OBSM_KEY, min_percentile: float = 30, max_percentile: float = 70) -> CompassComparison:
    """Sets up a comparison to single cell values. For every reaction, we calculate the values among all single cells in
    adata, and use the given percentiles for scaling. So the value of min_percentile will become 0, and max_percentile
    1."""
    return _ComparisonToSingleCellValues(adata=adata, model=model, groupby=groupby, obsm_key=obsm_key,
                                         min_percentile=min_percentile, max_percentile=max_percentile)


def load_compass_model(folder: str) -> CompassModel:
    """Reads the model JSON file in the given folder, and stores all the species and reactions in a CompassModel object."""
    file_path = os.path.join(folder, "model.json.gz")
    with gzip.open(file_path, "rt") as handle:
        model_json = handle.read()
    json_object = json.loads(model_json)

    compass_model = CompassModel()
    for species in json_object["species"].values():
        name = species["name"]
        if not name:
            continue
        compass_model.add_species(species["id"], name)

    for reaction in json_object["reactions"].values():
        reactants = reaction["reactants"]
        products = reaction["products"]
        id = reaction["id"]

        reactant_names = [compass_model.species(r) for r in reactants]
        product_names = [compass_model.species(p) for p in products]

        compass_model.add_reaction(id, reactant_names=reactant_names, product_names=product_names)

    return compass_model


def add_compass_output(adata: AnnData, compass_folder: str, *, obsm_key: str = _DEFAULT_OBSM_KEY,
                       microclustering_mapping: DataFrame | None = None, microclustering_column: str = "microclustering"):
    """Reads the Compass output for a given system and adds it to the AnnData object under adata.obsm, by default under
    the "compass" key.

    In case you ran Compass on microclusters, you'll also need to pass a dataframe mapping the microcluster names to
    the cell names (from `adata.obs_names`). The indices of the dataframe are assumed to be the cell names, and the
    values in the column named `microclustering_mapping` ("microclustering" by default) are expected to be the
    names of the microclusters."""

    reaction_scores = pandas.read_csv(os.path.join(compass_folder, "reactions.tsv"), delimiter="\t", index_col=0)

    if microclustering_mapping is not None:
        # Undo the microclustering
        if not microclustering_column in microclustering_mapping.columns:
            raise ValueError(f"Column '{microclustering_column}' not found in microclustering_mapping.")
        reaction_scores = reaction_scores[list(microclustering_mapping[microclustering_column])]
        reaction_scores.columns = microclustering_mapping.index

    adata.obsm[obsm_key] = reaction_scores.T


def color_model(adata: AnnData, kegg_map: KeggMap, *, groupby: str, group: str, reference: str = "rest", obsm_key: str = "compass",
                min_percentile: float = 10, max_percentile: float = 90):
    if not obsm_key in adata.obsm:
        raise ValueError(f"Key '{obsm_key}' not found in adata.obsm")

    reaction_scores = adata.obsm[obsm_key]




