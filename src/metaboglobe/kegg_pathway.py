from collections import defaultdict
from importlib.resources import files
from xml.etree import ElementTree

import matplotlib.patches
import numpy
import pandas
from enum import Enum, auto
from matplotlib.axes import Axes
from typing import NamedTuple, Collection, Any

from metaboglobe._util import optimize_for_display, wrap_text, optimize_for_matching, MPLColor, \
    get_name_without_enantiomers


def _read_accession_number_to_name_mapping() -> dict[str, list[str]]:
    """Reads the mapping from KEGG accession numbers to names from the "Kegg Metabolism" folder. In the returned mapping,
    the main name is always the first one in the list, and the synonyms are the following ones."""
    file_path = files("metaboglobe.data") / "kegg_compounds.tsv"

    mapping = dict()
    with file_path.open("r") as handle:
        line = handle.readline()
        while line != "":
            if not line.startswith(";"):
                accession_number, names = line.split("\t")
                mapping[accession_number] = _parse_name_list(names)

            line = handle.readline()
    return mapping


def _parse_name_list(names_unparsed: str) -> list[str]:
    """Puts undesirable names last, so that the most desireable names are more likely to be used for display."""
    desireable_names = list()
    undesireable_names = list()
    for name in names_unparsed.split(";"):
        name = name.strip()

        # Also add the names without enantiomer specifier
        # (but sometimes those variants are included by KEGG though, so we check for duplicates below)
        name_without_enantiomers = get_name_without_enantiomers(name)
        undesireable_names.append(name_without_enantiomers)

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





_ACCESSION_NUMBER_TO_NAMES = _read_accession_number_to_name_mapping()


class EntryType(Enum):
    MAP = auto()
    GENE = auto()
    ORTHOLOG = auto()
    COMPOUND = auto()
    TITLE = auto()

    def draw_entry(self, ax: Axes, entry: "KeggEntry"):
        if self == EntryType.TITLE:
            ax.set_title(entry.name)
        elif self == EntryType.MAP:
            display_name = wrap_text(entry.name, int(entry.width // 6))  # Wrap text to fit within the entry box, assuming an average character width of 6 pixels
            rect = matplotlib.patches.Rectangle((entry.x - entry.width / 2, entry.y - entry.height / 2),
                                                entry.width, entry.height, fill=True, color="#c1bcb8")
            ax.add_patch(rect)
            ax.text(entry.x, entry.y, display_name, ha="center", va="center", fontsize=6, zorder=10)
        elif self == EntryType.COMPOUND:
            display_names = _ACCESSION_NUMBER_TO_NAMES.get(entry.name)
            display_name = display_names[0] if display_names is not None else entry.name
            display_name = optimize_for_display(display_name)
            ax.text(entry.x, entry.y, display_name, ha="center", va="center", fontsize=6, zorder=10)
        elif self == EntryType.ORTHOLOG:
            return  # Not interested in orthologs, they just clutter the figure
        elif self == EntryType.GENE:  # Genes
            if "," in entry.name:
                display_name = entry.name.split(",")[0] + ", etc."
            else:
                display_name = entry.name
            ax.text(entry.x, entry.y, display_name, ha="center", va="center", fontsize=6, zorder=10)
        else:
            raise ValueError(f"Unknown entry type: {self}")


class RelationType(Enum):
    REACTION_REVERSIBLE = auto()  # One enzyme, goes both ways
    REACTION_TWO_IRREVERSIBLE = auto()  # Both ways are possible, but using different enzymes
    REACTION_IRREVERSIBLE = auto()  # Only one way is possible

    ECREL = auto()
    MAPLINK = auto()

    def is_reaction(self) -> bool:
        """Checks if this relation represents a reaction between compounds."""
        return (self == RelationType.REACTION_IRREVERSIBLE or self == RelationType.REACTION_REVERSIBLE
                or self == RelationType.REACTION_TWO_IRREVERSIBLE)


class KeggEntry(NamedTuple):
    entry_id: int
    name: str
    x: float
    y: float
    width: float
    height: float
    entry_type: EntryType


class KeggRelation(NamedTuple):

    """Represents a relation in the KEGG pathway map."""
    substrate_id: int
    product_id: int
    relation_type: RelationType

    @property
    def substrate_product_tuple(self) -> tuple[int, int]:
        """Gets a tuple of the substrate and product IDs."""
        return self.substrate_id, self.product_id


class KeggReactionWithReversion(NamedTuple):
    """Represents a reaction in the KEGG pathway map, as well as a flag to indicate whether we are looking at the
    reverse of the reaction."""
    relation: KeggRelation
    reversed: bool


def _check_for_match(matching_names: list[str], kegg_identifier: str) -> bool:
    """Checks if any of the matching_names (like "acetic-acid") matches any of the names behind the KEGG identifier (like "C00003").
    """
    for kegg_name in _ACCESSION_NUMBER_TO_NAMES.get(kegg_identifier, []):
        if optimize_for_matching(kegg_name) in matching_names:
            return True
    return False


class KeggMap:
    _entries_by_id: dict[int, KeggEntry]

    _relations: set[KeggRelation]
    _relations_forward_values: dict[tuple[int, int], float]
    _relations_backward_values: dict[tuple[int, int], float]

    _reaction_by_compound_names: dict[str, list[KeggRelation]]

    def __init__(self):
        self._entries_by_id = dict()

        self._relations = set()
        self._relations_forward_values = dict()
        self._relations_backward_values = dict()
        self._reaction_by_compound_names = defaultdict(list)

    def add_entry(self, entry: KeggEntry):
        """Adds an entry to the plot. Entries are compounds, orthologs, etc."""
        self._entries_by_id[entry.entry_id] = entry

    def _search_relation(self, substrate_id: int, product_id: int) -> KeggRelation | None:
        for relation in self._relations:
            if relation.substrate_id == substrate_id and relation.product_id == product_id:
                return relation
        return None

    def _remove_relation(self, relation: KeggRelation) -> None:
        """Removes a relation, both from the list and from the search indices. Raises KeyError if the relation didn't
        exist."""

        self._relations.remove(relation)
        if relation.relation_type.is_reaction():

            # Remove all substrates and products
            for entry_id in relation.substrate_product_tuple:
                entry = self._entries_by_id[entry_id]
                if entry.entry_type == EntryType.COMPOUND:

                    # Find all names of this entry that we know
                    names = _ACCESSION_NUMBER_TO_NAMES.get(entry.name, [])
                    for name in names:

                        reactions_with_name = self._reaction_by_compound_names[optimize_for_matching(name)]
                        reactions_with_name.remove(relation)

    def add_relation(self, substrate_id: int, product_id: int, relation_type: RelationType):
        """Relations connect entries, identified using the entry IDs."""

        relation = KeggRelation(substrate_id, product_id, relation_type)

        if relation_type.is_reaction():
            # Sometimes reversible reactions are annotated as two irreversible reactions in both directions.
            # In that case, we only want to add one reversible reaction, and ignore the other one.
            reversed_entry = self._search_relation(product_id, substrate_id)
            if reversed_entry is not None:
                # Found the reaction in the other direction, update it to be two-way, but still irreversible (because different enzymes)
                if reversed_entry.relation_type == RelationType.REACTION_IRREVERSIBLE:
                    # We need to remove it. The code below will add a new, two-irreversible reaction
                    self._remove_relation(reversed_entry)
                    relation = KeggRelation(product_id, substrate_id, RelationType.REACTION_TWO_IRREVERSIBLE)
                else:
                    return  # Can leave as-is, nothing to do

        self._relations.add(relation)

        # Populate compound name mappings
        if relation_type.is_reaction():
            for entry_id in [substrate_id, product_id]:
                entry = self._entries_by_id[entry_id]
                if entry.entry_type == EntryType.COMPOUND:
                    names = _ACCESSION_NUMBER_TO_NAMES.get(entry.name, [])
                    for name in names:
                        self._reaction_by_compound_names[optimize_for_matching(name)].append(relation)

    def match_reaction(self, substrate_names: list[str], product_names: list[str]) -> KeggReactionWithReversion | None:
        """Given a list of substrate names and product names, searches for a reaction in this pathway that matches
        one of the substrate names and one of the product names. The given names are expected to be names like
        "D-Fructose-6P". The method uses all the known names of Kegg to match. If a name is provided without
        specifying the stereoisomer (like "Fructose-6P", without the "D-") it can match to any stereoisomer in the
        pathway. Returns None if there was no match."""
        substrate_names = [optimize_for_matching(name) for name in substrate_names]
        product_names = [optimize_for_matching(name) for name in product_names]

        for substrate_name in substrate_names:
            for relation in self._reaction_by_compound_names.get(substrate_name, []):
                # Found a relation with at least one match somewhere

                relation_substrate_kegg_accession = self._entries_by_id[relation.substrate_id].name
                relation_product_kegg_accession = self._entries_by_id[relation.product_id].name

                if (_check_for_match([substrate_name], relation_substrate_kegg_accession)
                        and _check_for_match(product_names, relation_product_kegg_accession)):
                    # Found a match!
                    return KeggReactionWithReversion(relation, reversed=False)

                if (relation.relation_type != RelationType.REACTION_IRREVERSIBLE and
                        (_check_for_match([substrate_name], relation_product_kegg_accession)
                        and _check_for_match(product_names, relation_substrate_kegg_accession))):
                    # Found a reverse match!
                    return KeggReactionWithReversion(relation, reversed=True)
        return None

    def _relation_to_str(self, relation: KeggRelation) -> str:
        return (_ACCESSION_NUMBER_TO_NAMES[self._entries_by_id[relation.substrate_id].name][0] + " -> "
                + _ACCESSION_NUMBER_TO_NAMES[self._entries_by_id[relation.product_id].name][0])

    def entry(self, entry_id: int) -> KeggEntry | None:
        """Gets an entry from the map with the given id (as added using add_entry). Returns None if not found."""
        return self._entries_by_id.get(entry_id, None)

    @property
    def entries(self) -> Collection[KeggEntry]:
        """Gets all available entries."""
        return self._entries_by_id.values()

    @property
    def relations(self) -> Collection[KeggRelation]:
        """Gets all available relations."""
        return self._relations

    def set_reaction_score(self, reaction: KeggReactionWithReversion, score: float):
        """Adds a reaction score to the map with the given id, which can be used for coloring. Raises ValueError
        if a score was already set for the reaction in the given direction, or if the score is NaN."""
        if numpy.isnan(score):
            raise ValueError(f"Score is NaN for {self._relation_to_str(reaction.relation)}")

        if reaction.reversed:
            if reaction in self._relations_backward_values:
                raise ValueError(f"Duplicate backward score for {self._relation_to_str(reaction.relation)}")
            self._relations_backward_values[reaction.relation.substrate_product_tuple] = score
        else:
            if reaction in self._relations_forward_values:
                raise ValueError(f"Duplicate forward score for {self._relation_to_str(reaction.relation)}")
            self._relations_forward_values[reaction.relation.substrate_product_tuple] = score

    def forward_value(self, reaction: KeggRelation) -> float:
        """Gets the forward value for the given reaction, or NaN if not set."""
        return self._relations_forward_values.get(reaction.substrate_product_tuple, numpy.nan)

    def backward_value(self, reaction: KeggRelation) -> float:
        """Gets the backward value for the given reaction, or NaN if not set."""
        return self._relations_backward_values.get(reaction.substrate_product_tuple, numpy.nan)

    def has_relations(self, entry: KeggEntry) -> bool:
        """Searches for any relations for the given entry."""
        entry_id = entry.entry_id
        for relation in self.relations:
            if relation.substrate_id == entry_id or relation.product_id == entry_id:
                return True
        return False


def load_kegg_map(path: str) -> KeggMap:
    """Loads a KEGG map. These can be downloaded from webpages like https://rest.kegg.jp/get/hsa00010/kgml ."""
    with open(path) as handle:
        kgml = handle.read()
    root = ElementTree.fromstring(kgml)
    kegg_map = KeggMap()

    # Read entries
    for entry in root.findall("entry"):
        entry_id = int(entry.attrib["id"])
        graphics = entry.find("graphics")
        if graphics is None:
            continue
        name = graphics.attrib["name"]
        entry_type = EntryType[entry.attrib["type"].upper()]

        if name.startswith("TITLE:"):
            # Special-case for how KEGG stores the figure title (it's stored as a map)
            name = name[len("TITLE:"):]
            entry_type = EntryType.TITLE

        x = float(graphics.attrib["x"])
        y = float(graphics.attrib["y"])
        width = float(graphics.attrib["width"])
        height = float(graphics.attrib["height"])

        kegg_map.add_entry(entry_id, KeggEntry(entry_id, name, x, y, width, height, entry_type))

    # Read reactions
    for reaction in root.findall("reaction"):
        substrates = [int(s.attrib['id']) for s in reaction.findall("substrate")]
        products = [int(p.attrib['id']) for p in reaction.findall("product")]
        reversible = reaction.attrib.get("type", "reversible") == "reversible"
        for substrate in substrates:
            for product in products:
                kegg_map.add_relation(substrate, product,
                                      RelationType.REACTION_REVERSIBLE if reversible else RelationType.REACTION_IRREVERSIBLE)

    # Read other relations
    for relation in root.findall("relation"):
        entry1 = int(relation.attrib["entry1"])
        entry2 = int(relation.attrib["entry2"])
        relation_type = RelationType[relation.attrib["type"].upper()]
        kegg_map.add_relation(entry1, entry2, relation_type)

    return kegg_map
