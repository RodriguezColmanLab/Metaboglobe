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

from metaboglobe import kegg_pathway
from metaboglobe._translation_to_kegg_id import map_ecrel_to_kegg_reactions, map_gem_id_to_kegg_reactions
from metaboglobe._util import optimize_for_matching
from metaboglobe.kegg_pathway import KeggMap
from metaboglobe.kegg_species import KeggReactionId, KeggReactionIdWithReversion

_DEFAULT_OBSM_KEY = "compass"


class CompassReaction(NamedTuple):
    """A reaction as specified by the Compass model, still needs to be matched to a reaction in the KEGG pathway."""

    kegg_ids: list[KeggReactionId]  # Translated using lookup tables, may be empty

    # If that fails, try to match using reactant and product names
    reactant_names: list[str]
    product_names: list[str]
    enzyme_names: list[str]


class CompassModel:
    """A class for mapping the reaction and metabolite IDs in the Compass output to the names."""

    _species_names_by_id: dict[str, str]
    _reactions_by_id: dict[str, CompassReaction]


    def __init__(self):
        self._reactions_by_id = dict()
        self._species_names_by_id = dict()

    def add_species(self, id: str, name: str):
        self._species_names_by_id[id] = name

    def species(self, id: str) -> str:
        """Gets the species with the given ID. For example, for "MAM01840e" it could return "fructose". Returns
        the ID itself if not found."""
        return self._species_names_by_id.get(id, id)

    def add_reaction(self, original_id: str, compass_reaction: CompassReaction):
        """Adds the reaction with the given ID, reactants and products. The original ID will be something like "MAR04363_pos"."""
        self._reactions_by_id[original_id] = compass_reaction

    def reaction(self, id: str) -> CompassReaction | None:
        """Gets the reactions with the given ID. For example, for "MAR04363_pos". Returns None if not found."""
        return self._reactions_by_id.get(id)

    def update(self, other: "CompassModel"):
        """Adds all the reactions and species from the other model to this model. If there are reactions or species
        with the same ID, they will be overwritten."""
        self._species_names_by_id.update(other._species_names_by_id)
        self._reactions_by_id.update(other._reactions_by_id)


class CompassComparison(ABC):
    """For plotting a comparison on a KEGG map."""

    @abstractmethod
    def insert_values_in_map(self, kegg_map: KeggMap, group: str):
        """Applies the comparison to the given KEGG map, by coloring the reactions in the map according to the values
        of the given group."""
        return NotImplemented


def _match_reactions_to_map(compass_reaction: CompassReaction, kegg_map: KeggMap
                            ) -> list[KeggReactionIdWithReversion]:
    """As described in the README, we use two strategies for matching a Compass reaction to the pathway:

    - first we try using the translated KEGG Reaction ID. However, there are some edge cases where the ID mapping
      fails, and the same reaction in the pathway is not found.
    - second, we try matching by substrate, product and enzyme name. For each of these three, at least one needs to
      match.
    """
    matches = list()
    for kegg_id in compass_reaction.kegg_ids:
        map_reactions = kegg_map.reaction_arrows_by_id(kegg_id)
        if len(map_reactions) > 0:
            # Try matching using KEGG id, we have results
            # Just need to figure out whether our reaction is reversed compared to the reaction on the map, or not
            substrate_names = [optimize_for_matching(name) for name in compass_reaction.reactant_names]
            product_names = [optimize_for_matching(name) for name in compass_reaction.product_names]
            for map_reaction in map_reactions:
                if kegg_pathway.check_for_match(substrate_names, map_reaction.product.compound_id) \
                    or kegg_pathway.check_for_match(product_names, map_reaction.substrate.compound_id):
                    # Appears that reaction is reversed
                    matches.append(KeggReactionIdWithReversion(kegg_id, reversed=True))
                else:
                    matches.append(KeggReactionIdWithReversion(kegg_id, reversed=False))
    if len(matches) > 0:
        return matches  # Succesfull using KEGG id match

    # Try matching by name as fallback
    return list(kegg_map.match_reactions(substrate_names=compass_reaction.reactant_names,
                                         product_names=compass_reaction.product_names,
                                         enzyme_names=compass_reaction.enzyme_names))


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
            compass_reaction = self._model.reaction(reaction_id)
            if compass_reaction is None:
                continue
            reactions = _match_reactions_to_map(compass_reaction, kegg_map)
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

    # Read in the compound names
    for species in json_object["species"].values():
        name = species["name"]
        if not name:
            continue
        compass_model.add_species(species["id"], name)

    # Read in the reactions
    for reaction in json_object["reactions"].values():
        reaction_id = reaction["id"]

        kegg_ids = []

        # Translate the ID to the KEGG reaction id
        if "meta" in reaction and "ECNumber" in reaction["meta"]:
            # RECON model
            ec_number = reaction["meta"]["ECNumber"]
            kegg_ids += map_ecrel_to_kegg_reactions(ec_number)
        else:
            # Human1 / Mouse1 model
            reaction_id_without_direction = reaction["id"]
            if reaction_id_without_direction.endswith("_pos") or reaction_id_without_direction.endswith("_neg"):
                reaction_id_without_direction = reaction_id_without_direction[:-4]
            kegg_ids += map_gem_id_to_kegg_reactions(reaction_id_without_direction)

        reactants = reaction["reactants"].keys()
        products = reaction["products"].keys()
        gene_names = list()
        for gene in reaction["genes"].values():
            gene_names.append(gene["name"])
        reactant_names = [compass_model.species(r) for r in reactants]
        product_names = [compass_model.species(p) for p in products]

        compass_model.add_reaction(reaction_id, CompassReaction(kegg_ids=kegg_ids, reactant_names=reactant_names, product_names=product_names, enzyme_names=gene_names))

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

