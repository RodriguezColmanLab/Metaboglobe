import gzip
import json
import os
from abc import ABC, abstractmethod
from collections import defaultdict
from importlib.resources import files
from typing import NamedTuple

import numpy
import pandas
import scanpy.get
from anndata import AnnData
from pandas import DataFrame

from metaboglobe._translation_to_kegg_id import map_ecrel_to_kegg_reactions, map_gem_id_to_kegg_reactions
from metaboglobe.kegg_pathway import KeggMap, KeggReactionIdWithReversion
from metaboglobe.kegg_species import KeggReactionId


_DEFAULT_OBSM_KEY = "compass"


class CompassReaction(NamedTuple):
    reactant_names: list[str]
    product_names: list[str]


class CompassModel:
    """A class for mapping the reaction and metabolite IDs in the Compass output to the names."""
    _reactions_by_id: dict[str, list[KeggReactionIdWithReversion]]

    def __init__(self):
        self._reactions_by_id = dict()

    def add_reaction(self, original_id: str, kegg_ids: list[KeggReactionId], reversed: bool):
        """Adds the reaction with the given ID, reactants and products. The original ID will be something like "MAR04363_pos"."""
        self._reactions_by_id[original_id] = [KeggReactionIdWithReversion(kegg_id, reversed) for kegg_id in kegg_ids]

    def reactions(self, id: str) -> list[KeggReactionIdWithReversion]:
        """Gets the reactions with the given ID. For example, for "MAR04363_pos". Returns None if not found."""
        return self._reactions_by_id.get(id, [])

    def update(self, other: "CompassModel"):
        """Adds all the reactions and species from the other model to this model. If there are reactions or species
        with the same ID, they will be overwritten."""
        self._reactions_by_id.update(other._reactions_by_id)


class CompassComparison(ABC):
    """For plotting a comparison on a KEGG map."""

    @abstractmethod
    def insert_values_in_map(self, kegg_map: KeggMap, group: str):
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

    def insert_values_in_map(self, kegg_map: KeggMap, group: str):
        """Applies the values into the given KeggMap. They will be on a scale of 0 to 1."""
        reaction_ids = self._reactions_aggregated.var_names

        group_index = self._reactions_aggregated.obs_names.get_loc(group)
        reaction_values = self._reactions_aggregated.layers["mean"][group_index]

        for reaction_id, reaction_value, min_value, max_value in zip(reaction_ids, reaction_values,
                                                                     self._min_values_by_reaction, self._max_values_by_reaction):
            reactions = self._model.reactions(reaction_id)
            if len(reactions) == 0:
                continue

            if min_value == max_value:
                scaled_reaction = 0.5  # All have the same value, so we just set them to the middle of the scale
            else:
                scaled_reaction = (reaction_value - min_value) / (max_value - min_value)
                if scaled_reaction < 0:
                    scaled_reaction = 0
                elif scaled_reaction > 1:
                    scaled_reaction = 1
            for reaction in reactions:
                kegg_map.set_reaction_score(reaction, scaled_reaction)


def setup_comparison_to_single_cells(adata: AnnData, model: CompassModel, *, groupby: str, obsm_key: str = _DEFAULT_OBSM_KEY, min_percentile: float = 30, max_percentile: float = 70) -> CompassComparison:
    """Sets up a comparison to single cell values. For every reaction, we calculate the values among all single cells in
    adata, and use the given percentiles for scaling. So the value of min_percentile will become 0, and max_percentile
    1."""
    return _ComparisonToSingleCellValues(adata=adata, model=model, groupby=groupby, obsm_key=obsm_key,
                                         min_percentile=min_percentile, max_percentile=max_percentile)


def load_compass_model(folder: str) -> CompassModel:
    """Reads the model JSON file in the given folder, and stores all the species and reactions in a CompassModel object.

    The folder can either be a Compass model folder (containing a model.json.gz file and a reactions.tsv file), or
    a folder containing subfolders, each containing a Compass model folder. The latter is useful if you ran Compass
    on several subsystems using Module-Compass.
    """

    # Case where we have a single model
    if os.path.exists(os.path.join(folder, "model.json.gz")):
        return _load_single_compass_model(folder)

    # Case where each subfolder is a model
    model = None
    for subfolder_name in os.listdir(folder):
        subfolder_path = os.path.join(folder, subfolder_name)
        if os.path.exists(os.path.join(subfolder_path, "model.json.gz")):
            subfolder_model = _load_single_compass_model(subfolder_path)
            if model is None:
                model = subfolder_model
            else:
                model.update(subfolder_model)
    if model is None:
        raise ValueError(f"No Compass model found in folder '{folder}' or its subfolders.")
    return model


def _load_single_compass_model(folder: str) -> CompassModel:
    file_path = os.path.join(folder, "model.json.gz")
    with gzip.open(file_path, "rt") as handle:
        model_json = handle.read()
    json_object = json.loads(model_json)

    compass_model = CompassModel()
    for reaction in json_object["reactions"].values():
        reaction_id = reaction["id"]

        # Check if positive or negative based on suffix, and remove said suffix
        reversed = False
        if reaction_id.endswith("_pos"):
            reaction_id = reaction_id[:-4]
        elif reaction_id.endswith("_neg"):
            reaction_id = reaction_id[:-4]
            reversed = True

        # Translate the ID to the KEGG reaction id
        if "meta" in reaction and "ECNumber" in reaction["meta"]:
            # RECON model
            ec_number = reaction["meta"]["ECNumber"]
            for kegg_id in map_ecrel_to_kegg_reactions(ec_number):
                compass_model.add_reaction(reaction["id"], kegg_id, reversed)
        else:
            # Human1 / Mouse1 model
            kegg_id = map_gem_id_to_kegg_reactions(reaction_id)
            if kegg_id is None:
                continue
            compass_model.add_reaction(reaction["id"], kegg_id, reversed)

    return compass_model


def add_compass_output(adata: AnnData, compass_folder: str, *, obsm_key: str = _DEFAULT_OBSM_KEY,
                       microclustering_mapping: DataFrame | None = None, microclustering_column: str = "microclustering"):
    """Reads the Compass output for a given system and adds it to the AnnData object under adata.obsm, by default under
    the "compass" key.

    In case you ran Compass on microclusters, you'll also need to pass a dataframe mapping the microcluster names to
    the cell names (from `adata.obs_names`). The indices of the dataframe are assumed to be the cell names, and the
    values in the column named `microclustering_mapping` ("microclustering" by default) are expected to be the
    names of the microclusters."""

    if os.path.exists(os.path.join(compass_folder, "reactions.tsv")):
        reaction_penalties = pandas.read_csv(os.path.join(compass_folder, "reactions.tsv"), delimiter="\t", index_col=0)
    else:
        # Search in subfolders
        reaction_penalties = None
        for subfolder_name in os.listdir(compass_folder):
            subfolder_path = os.path.join(compass_folder, subfolder_name)
            if os.path.exists(os.path.join(subfolder_path, "reactions.tsv")):
                # Found a subfolder
                reaction_penalties_of_subfolder = pandas.read_csv(os.path.join(subfolder_path, "reactions.tsv"), delimiter="\t", index_col=0)
                if reaction_penalties is None:
                    reaction_penalties = reaction_penalties_of_subfolder
                else:
                    # Simply concatenate without verifying integrity. We resolve the duplicates later
                    reaction_penalties = pandas.concat([reaction_penalties, reaction_penalties_of_subfolder], axis=0, verify_integrity=False)
        if reaction_penalties is None:
            raise ValueError(f"No reactions.tsv file found in folder '{compass_folder}' or its subfolders.")

        # Remove duplicates (reactions can be part of multiple subsystems)
        reaction_penalties = reaction_penalties[~reaction_penalties.index.duplicated(keep="first")]

    if microclustering_mapping is not None:
        # Undo the microclustering
        if not microclustering_column in microclustering_mapping.columns:
            raise ValueError(f"Column '{microclustering_column}' not found in microclustering_mapping.")
        reaction_penalties = reaction_penalties[list(microclustering_mapping[microclustering_column])]
        reaction_penalties.columns = microclustering_mapping.index

    # To convert to scores, take the minus log of the penalties (as done on
    # https://github.com/wagnerlab-berkeley/Compass/blob/6084e47607bd2258534267a31b6f951fa25b6700/docs/notebooks/compass_analysis.py#L43 )
    adata.obsm[obsm_key] = -numpy.log(reaction_penalties.T + 1)

