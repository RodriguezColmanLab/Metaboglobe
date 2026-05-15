from collections import defaultdict
from enum import Enum, auto
from importlib.resources import files
from typing import NamedTuple, Collection
from xml.etree import ElementTree

import numpy

from metaboglobe._util import optimize_for_matching, get_names_without_stereoisomers, optimize_for_display
from metaboglobe.kegg_species import KeggReactionId, KeggCompoundId, KeggPathwayId
from metaboglobe.math.box_2d import Box2
from metaboglobe.math.vector_2d import Vector2


def _read_accession_number_to_name_mapping() -> dict[KeggCompoundId, list[str]]:
    """Reads the mapping from KEGG accession numbers to names from the "Kegg Metabolism" folder. In the returned mapping,
    the main name is always the first one in the list, and the synonyms are the following ones."""
    file_path = files("metaboglobe.data") / "kegg_compounds.tsv"

    mapping = dict()
    with file_path.open("r") as handle:
        line = handle.readline()
        while line != "":
            if not line.startswith(";"):
                accession_number, names = line.split("\t")
                mapping[KeggCompoundId.create_from_id(accession_number)] = _parse_name_list(names)

            line = handle.readline()
    return mapping


def _parse_name_list(names_unparsed: str) -> list[str]:
    """Puts undesirable names last, so that the most desireable names are more likely to be used for display."""
    desireable_names = list()
    undesireable_names = list()
    for name in names_unparsed.split(";"):
        name = name.strip()

        # Also add the names without stereoisomer specifier
        # (but sometimes those variants are included by KEGG though, so we check for duplicates below)
        undesireable_names += list(get_names_without_stereoisomers(name))

        if "[" in name:
            undesireable_names.append(name)  # Too complex of a name, decrease priority
        else:
            desireable_names.append(name)

    # Return both lists, preserving order, removing any duplicates
    items_seen = set()
    return_list = list()
    for item in desireable_names + undesireable_names:
        if item in items_seen:
            continue
        return_list.append(item)
        items_seen.add(item)
    return return_list


def _parse_kegg_reaction_ids(string: str) -> list[KeggReactionId]:
    """Parses the KEGG reaction IDs from the string, as stored in the xml file. Raises ValueError if parsing fails."""
    parts = string.split(" ")
    return [KeggReactionId(part) for part in parts]


_ACCESSION_NUMBER_TO_NAMES = _read_accession_number_to_name_mapping()


class ReactionType(Enum):
    REVERSIBLE = auto()  # One enzyme, goes both ways
    IRREVERSIBLE = auto()  # Only one way is possible


class KeggCompoundInMap(NamedTuple):
    compound_id: KeggCompoundId
    entry_in_map_id: int
    name: str
    x: float
    y: float


class KeggPathwayReferenceInMap(NamedTuple):
    pathway_id: KeggPathwayId
    entry_in_map_id: int
    name: str
    x: float
    y: float
    width: float
    height: float


class KeggReactionEnzymeInMap(NamedTuple):
    """Represents an enzyme in the KEGG pathway map. Does not state what the substrates and products are."""
    reaction_id: KeggReactionId
    entry_in_map_id: int
    name: str
    x: float
    y: float
    width: float
    height: float


class KeggOtherPathwayRelationInMap(NamedTuple):
    """Represents a relation in the KEGG pathway map."""
    compound: KeggCompoundInMap
    pathway: KeggPathwayReferenceInMap

    @property
    def compound_id(self) -> KeggCompoundId:
        return self.compound.compound_id

    @property
    def pathway_id(self) -> KeggPathwayId:
        return self.pathway.pathway_id


class KeggReactionArrow(NamedTuple):
    """Represents a chemical reaction in the KEGG pathway map. Should contain all the information to draw the arrow."""
    reaction_in_map: KeggReactionEnzymeInMap
    substrate: KeggCompoundInMap
    product: KeggCompoundInMap
    reaction_type: ReactionType

    @property
    def substrate_id(self) -> KeggCompoundId:
        return self.substrate.compound_id

    @property
    def product_id(self) -> KeggCompoundId:
        return self.product.compound_id

    @property
    def reaction_id(self) -> KeggReactionId:
        return self.reaction_in_map.reaction_id


class KeggReactionIdWithReversion(NamedTuple):
    """Represents a reaction in the KEGG pathway map, as well as a flag to indicate whether we are looking at the
    reverse of the reaction."""
    reaction_id: KeggReactionId
    reversed: bool


def get_display_name(kegg_accession_id: KeggCompoundId) -> str:
    """Gets the display name (including LaTeX codes) for the given Kegg Accession ID (like "C00092")."""
    display_names = _ACCESSION_NUMBER_TO_NAMES.get(kegg_accession_id)
    display_name = display_names[0] if display_names is not None else kegg_accession_id.compound_id
    return optimize_for_display(display_name)


def _check_for_match(matching_names: list[str], kegg_identifier: KeggCompoundId) -> bool:
    """Checks if any of the matching_names (like "acetic-acid") matches any of the names behind the KEGG identifier (like "C00003").
    """
    for kegg_name in _ACCESSION_NUMBER_TO_NAMES.get(kegg_identifier, []):
        if optimize_for_matching(kegg_name) in matching_names:
            return True
    return False


class KeggMap:
    """A map of a single KEGG pathway.

    The data strutures are a bit convoluted; this is because metabolites, enzymes, pathways, etc. can appear multiple
    times in the map. So we can't simply index everything by the KEGG id, but have to keep track if local ids within
    the pathway.
    """

    _title: str

    _reaction_enzymes: list[KeggReactionEnzymeInMap]
    _compounds: list[KeggCompoundInMap]
    _pathway_references: list[KeggPathwayReferenceInMap]

    _relations_to_other_pathways: list[KeggOtherPathwayRelationInMap]
    _reaction_arrows: list[KeggReactionArrow]
    _reactions_forward_values: dict[KeggReactionId, float]
    _reactions_backward_values: dict[KeggReactionId, float]
    _reaction_by_compound_names: dict[str, list[KeggReactionArrow]]

    def __init__(self):
        self._title = "Unnamed"

        self._reaction_enzymes = list()
        self._compounds = list()
        self._pathway_references = list()
        self._relations_to_other_pathways = list()

        self._reaction_arrows = list()
        self._reactions_forward_values = dict()
        self._reactions_backward_values = dict()
        self._reaction_by_compound_names = defaultdict(list)

    def box(self) -> Box2:
        """Calculates the axis limits for plotting this pathway."""
        xmax = 0
        ymax = 0
        for entry in self._compounds:
            xmax = max(entry.x, xmax)
            ymax = max(entry.y, ymax)
        for entry in self._pathway_references:
            xmax = max(entry.x + entry.width, xmax)
            ymax = max(entry.y + entry.height, ymax)

        if xmax == 0 or ymax == 0:
            return Box2(Vector2(0, 0), Vector2(1, 1)) # Invalid limits

        # Give some extra space
        xsize = xmax
        ysize = ymax
        xmax += xsize / 10
        ymax += ysize / 10

        return Box2(Vector2(0, 0), Vector2(xmax, ymax))

    def figsize(self) -> tuple[float, float]:
        """Calculates how big the figure size would be (in inches), assuming we only plot this pathway."""
        box = self.box()

        return box.max.x / 100, box.max.y / 100

    @property
    def title(self) -> str:
        """Gets the title of this pathway."""
        return self._title

    @title.setter
    def title(self, title: str):
        if title is None:
            raise ValueError("Title cannot be None")
        self._title = title

    def add_compound(self, compound_id: KeggCompoundId, entry_in_map_id: int, x: float, y: float) -> KeggCompoundInMap:
        """Adds a compound to the map at the given location."""
        compound = KeggCompoundInMap(compound_id, entry_in_map_id, get_display_name(compound_id), x, y)
        self._compounds.append(compound)
        return compound

    def add_reaction_enzyme(self, reaction_id: KeggReactionId, entry_in_map_id: int, name: str, x: float, y: float, width: float, height: float):
        """Adds a reaction enzyme to the map at the given location."""
        enzyme = KeggReactionEnzymeInMap(reaction_id, entry_in_map_id, name, x, y, width, height)
        self._reaction_enzymes.append(enzyme)
        return enzyme

    def add_pathway_reference(self, pathway_id: KeggPathwayId, entry_in_map_id: int, name: str, x: float, y: float, width: float, height: float):
        """Adds a pathway reference to the map at the given location."""
        pathway = KeggPathwayReferenceInMap(pathway_id, entry_in_map_id, name, x, y, width, height)
        self._pathway_references.append(pathway)
        return pathway

    def add_reaction_arrow(self, enzyme: KeggReactionEnzymeInMap, substrate: KeggCompoundInMap, product: KeggCompoundInMap, reaction_type: ReactionType):
        """Add the reaction arrow, connecting a substrate, enzyme and product."""

        reaction = KeggReactionArrow(enzyme, substrate, product, reaction_type)
        self._reaction_arrows.append(reaction)

        # Populate compound name mappings
        for kegg_compound_id in [substrate.compound_id, product.compound_id]:
            names = _ACCESSION_NUMBER_TO_NAMES.get(kegg_compound_id, [])
            for name in names:
                self._reaction_by_compound_names[optimize_for_matching(name)].append(reaction)

    def add_pathway_relation(self, compound: KeggCompoundInMap, pathway: KeggPathwayReferenceInMap):
        """Relations connect entries, identified using the entry IDs."""
        self._relations_to_other_pathways.append(KeggOtherPathwayRelationInMap(compound, pathway))

    def match_reaction(self, substrate_names: list[str], product_names: list[str]) -> KeggReactionIdWithReversion | None:
        """Given a list of substrate names and product names, searches for a reaction in this pathway that matches
        one of the substrate names and one of the product names. The given names are expected to be names like
        "D-Fructose-6P". The method uses all the known names of Kegg to match. If a name is provided without
        specifying the stereoisomer (like "Fructose-6P", without the "D-") it can match to any stereoisomer in the
        pathway. Returns None if there was no match."""
        substrate_names = [optimize_for_matching(name) for name in substrate_names]
        product_names = [optimize_for_matching(name) for name in product_names]

        for substrate_name in substrate_names:
            for reaction in self._reaction_by_compound_names.get(substrate_name, []):
                # Found a relation with at least one match somewhere

                relation_substrate_kegg_accession = reaction.substrate_id
                relation_product_kegg_accession = reaction.product_id

                if (_check_for_match([substrate_name], relation_substrate_kegg_accession)
                        and _check_for_match(product_names, relation_product_kegg_accession)):
                    # Found a match!
                    return KeggReactionIdWithReversion(reaction.reaction_in_map.reaction_id, reversed=False)

                if (reaction.reaction_type != ReactionType.IRREVERSIBLE and
                        (_check_for_match([substrate_name], relation_product_kegg_accession)
                        and _check_for_match(product_names, relation_substrate_kegg_accession))):
                    # Found a reverse match!
                    return KeggReactionIdWithReversion(reaction.reaction_in_map.reaction_id, reversed=True)
        return None

    def _reaction_id_to_str(self, reaction_id: KeggReactionId) -> str:
        for reaction_arrow in self._reaction_arrows:
            if reaction_arrow.reaction_id == reaction_id:
                return (reaction_id.reaction_id + ": "
                        + (_ACCESSION_NUMBER_TO_NAMES[reaction_arrow.substrate_id][0] + " -> "
                        + _ACCESSION_NUMBER_TO_NAMES[reaction_arrow.product_id][0]))
        else:
            return reaction_id.reaction_id + ": (reaction not in map)"

    @property
    def compounds(self) -> Collection[KeggCompoundInMap]:
        """Gets all compounds in the map."""
        return self._compounds

    @property
    def other_pathways(self) -> Collection[KeggPathwayReferenceInMap]:
        """Gets all pathway references in the map."""
        return self._pathway_references

    def other_pathway(self, pathway_id: KeggPathwayId) -> KeggPathwayReferenceInMap:
        """Gets the pathway reference with the given id from the map. Raises ValueError if the pathway is not in the map."""
        return self._pathway_references[pathway_id]

    @property
    def reaction_enzymes(self) -> Collection[KeggReactionEnzymeInMap]:
        """Gets all compounds in the map."""
        return self._reaction_enzymes

    @property
    def reaction_arrows(self) -> Collection[KeggReactionArrow]:
        """Gets all available reaction."""
        return self._reaction_arrows

    @property
    def other_pathway_relations(self) -> Collection[KeggOtherPathwayRelationInMap]:
        """Gets all available relations to other pathways."""
        return self._relations_to_other_pathways

    def set_reaction_score(self, reaction: KeggReactionIdWithReversion, score: float):
        """Adds a reaction score to the map with the given id, which can be used for coloring. Raises ValueError
        if a score was already set for the reaction in the given direction, or if the score is NaN."""
        if numpy.isnan(score):
            raise ValueError(f"Score is NaN for {self._reaction_id_to_str(reaction.reaction_id)}")

        if reaction.reversed:
            if reaction in self._reactions_backward_values:
                raise ValueError(f"Duplicate backward score for {self._reaction_id_to_str(reaction.reaction_id)}")
            self._reactions_backward_values[reaction.reaction_id] = score
        else:
            if reaction in self._reactions_forward_values:
                raise ValueError(f"Duplicate forward score for {self._reaction_id_to_str(reaction.reaction_id)}")
            self._reactions_forward_values[reaction.reaction_id] = score

    def forward_value(self, reaction_id: KeggReactionId) -> float:
        """Gets the forward value for the given reaction, or NaN if not set."""
        return self._reactions_forward_values.get(reaction_id, numpy.nan)

    def backward_value(self, reaction_id: KeggReactionId) -> float:
        """Gets the backward value for the given reaction, or NaN if not set."""
        return self._reactions_backward_values.get(reaction_id, numpy.nan)

    def has_relations_or_reactions(self, entry: KeggCompoundId | KeggReactionId | KeggPathwayId) -> bool:
        """Searches for any reactions for the given entry."""
        if isinstance(entry, KeggCompoundId):
            for reaction in self.reaction_arrows:
                if reaction.substrate_id == entry or reaction.product_id == entry:
                    return True
            for pathway_link in self.other_pathway_relations:
                if pathway_link.compound_id == entry:
                    return True
            return False
        elif isinstance(entry, KeggReactionId):
            for reaction in self.reaction_arrows:
                if reaction.reaction_id == entry:
                    return True
            return False
        elif isinstance(entry, KeggPathwayId):
            for relation in self._relations_to_other_pathways:
                if relation.pathway_id == entry:
                    return True
            return False
        raise ValueError("Expected compound, reaction or pathway id")


def load_kegg_map(path: str) -> KeggMap:
    """Loads a KEGG map. These can be downloaded from webpages like https://rest.kegg.jp/get/hsa00010/kgml ."""
    with open(path) as handle:
        kgml = handle.read()
    root = ElementTree.fromstring(kgml)

    # The KEGG map, plus our own indices using local pathway ids
    kegg_map = KeggMap()
    compound_entries = dict()
    gene_entries = dict()
    map_entries = dict()

    # Read entries
    for entry in root.findall("entry"):
        entry_id = int(entry.attrib["id"])
        graphics = entry.find("graphics")
        if graphics is None or "name" not in graphics.attrib:
            continue
        display_name = graphics.attrib["name"]
        entry_type = entry.attrib["type"]

        if display_name.startswith("TITLE:"):
            # Special-case for how KEGG stores the figure title (it's stored as a map)
            display_name = display_name[len("TITLE:"):]
            kegg_map.title = display_name
            continue

        x = float(graphics.attrib["x"])
        y = float(graphics.attrib["y"])
        width = float(graphics.attrib["width"])
        height = float(graphics.attrib["height"])

        if entry_type == "compound":
            compound_id = KeggCompoundId(entry.attrib["name"])
            compound_in_map = kegg_map.add_compound(compound_id, entry_id, x, y)
            compound_entries[entry_id] = compound_in_map
        elif entry_type == "gene":
            if "reaction" not in entry.attrib:
                continue
            reaction_id = KeggReactionId(entry.attrib["reaction"])
            reaction_enzyme = kegg_map.add_reaction_enzyme(reaction_id, entry_id, display_name, x, y, width, height)
            gene_entries[entry_id] = reaction_enzyme
        elif entry_type == "map":
            pathway_id = KeggPathwayId(entry.attrib["name"])
            pathway = kegg_map.add_pathway_reference(pathway_id, entry_id, display_name, x, y, width, height)
            map_entries[entry_id] = pathway

    # Read reactions
    for reaction in root.findall("reaction"):
        kegg_reaction_enzyme = gene_entries[int(reaction.attrib["id"])]
        substrates = [compound_entries[int(s.attrib['id'])] for s in reaction.findall("substrate")]
        products = [compound_entries[int(p.attrib['id'])] for p in reaction.findall("product")]
        reversible = reaction.attrib.get("type", "reversible") == "reversible"
        reaction_type = ReactionType.REVERSIBLE if reversible else ReactionType.IRREVERSIBLE

        for substrate in substrates:
            for product in products:
                kegg_map.add_reaction_arrow(kegg_reaction_enzyme, substrate, product, reaction_type)

    # Read relations
    for relation in root.findall("relation"):
        if relation.attrib["type"] != "maplink":
            continue

        entry1 = int(relation.attrib["entry1"])
        entry2 = int(relation.attrib["entry2"])
        subtype = relation.find("subtype")
        if subtype is None or subtype.attrib["name"] != "compound":
            raise ValueError("Unexpected subtype,", subtype)
        compound_id = int(subtype.attrib["value"])

        connected_pathway = None
        connected_compounds = list()
        if entry1 in map_entries:
            connected_pathway = map_entries[entry1]
            if entry2 in compound_entries:
                connected_compounds.append(compound_entries[entry2])
        elif entry2 in map_entries:
            connected_pathway = map_entries[entry2]
            if entry1 in compound_entries:
                connected_compounds.append(compound_entries[entry1])
        else:
            continue  # No pathway found
        if compound_id in compound_entries:
            connected_compounds.append(compound_entries[compound_id])

        for compound in connected_compounds:
            kegg_map.add_pathway_relation(compound, connected_pathway)

    return kegg_map
