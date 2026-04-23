from collections import defaultdict
from enum import Enum, auto
from importlib.resources import files
from typing import NamedTuple, Collection, Iterable
from xml.etree import ElementTree

import numpy

from metaboglobe._util import optimize_for_matching, get_names_without_stereoisomers, optimize_for_display


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


_ACCESSION_NUMBER_TO_NAMES = _read_accession_number_to_name_mapping()


class EntryType(Enum):
    MAP = auto()
    GENE = auto()
    ORTHOLOG = auto()
    COMPOUND = auto()
    TITLE = auto()


class ReactionType(Enum):
    REVERSIBLE = auto()  # One enzyme, goes both ways
    IRREVERSIBLE = auto()  # Only one way is possible


class RelationType(Enum):
    ECREL = auto()
    MAPLINK = auto()


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
    from_id: int
    to_id: int
    compound_id: int
    relation_type: RelationType


class KeggReaction(NamedTuple):
    """Represents a chemical reaction in the KEGG pathway map."""
    substrate_id: int
    product_id: int
    gene_id: int
    reaction_type: ReactionType

    @property
    def substrate_product_tuple(self) -> tuple[int, int]:
        """Gets a tuple of the substrate and product IDs."""
        return self.substrate_id, self.product_id


class KeggReactionWithReversion(NamedTuple):
    """Represents a reaction in the KEGG pathway map, as well as a flag to indicate whether we are looking at the
    reverse of the reaction."""
    reaction: KeggReaction
    reversed: bool


def get_display_name(kegg_accession_id: str) -> str:
    """Gets the display name (including LaTeX codes) for the given Kegg Accession ID (like "C00092")."""
    display_names = _ACCESSION_NUMBER_TO_NAMES.get(kegg_accession_id)
    display_name = display_names[0] if display_names is not None else kegg_accession_id
    return optimize_for_display(display_name)


def _check_for_match(matching_names: list[str], kegg_identifier: str) -> bool:
    """Checks if any of the matching_names (like "acetic-acid") matches any of the names behind the KEGG identifier (like "C00003").
    """
    for kegg_name in _ACCESSION_NUMBER_TO_NAMES.get(kegg_identifier, []):
        if optimize_for_matching(kegg_name) in matching_names:
            return True
    return False


class KeggMap:
    _entries_by_id: dict[int, KeggEntry]

    _relations_by_compound_id: dict[int, set[KeggRelation]]

    _reactions: set[KeggReaction]
    _reactions_forward_values: dict[tuple[int, int], float]
    _reactions_backward_values: dict[tuple[int, int], float]
    _reaction_by_compound_names: dict[str, list[KeggReaction]]

    def __init__(self):
        self._entries_by_id = dict()
        self._relations_by_compound_id = defaultdict(set)

        self._reactions = set()
        self._reactions_forward_values = dict()
        self._reactions_backward_values = dict()
        self._reaction_by_compound_names = defaultdict(list)

    def add_entry(self, entry: KeggEntry):
        """Adds an entry to the plot. Entries are compounds, orthologs, etc."""
        self._entries_by_id[entry.entry_id] = entry

    def _search_reaction(self, substrate_id: int, product_id: int) -> KeggReaction | None:
        for reaction in self._reactions:
            if reaction.substrate_id == substrate_id and reaction.product_id == product_id:
                return reaction
        return None

    def add_reaction(self, substrate_id: int, product_id: int, enzyme_id: int, reaction_type: ReactionType):
        """Relations connect compounds, identified using the entry IDs. The enzyme_id must point at a GENE entry."""
        if not self._entries_by_id[substrate_id].entry_type == EntryType.COMPOUND:
            raise ValueError("Substrate ID must be of type COMPOUND")
        if not self._entries_by_id[product_id].entry_type == EntryType.COMPOUND:
            raise ValueError("Product ID must be of type COMPOUND")
        if not self._entries_by_id[enzyme_id].entry_type == EntryType.GENE:
            raise ValueError("Enzyme ID must be of type GENE")

        reaction = KeggReaction(substrate_id, product_id, enzyme_id, reaction_type)
        self._reactions.add(reaction)

        # Populate compound name mappings
        for entry_id in [substrate_id, product_id]:
            entry = self._entries_by_id[entry_id]
            if entry.entry_type == EntryType.COMPOUND:
                names = _ACCESSION_NUMBER_TO_NAMES.get(entry.name, [])
                for name in names:
                    self._reaction_by_compound_names[optimize_for_matching(name)].append(reaction)

    def add_relation(self, from_id: int, to_id: int, compound_id: int, relation_type: RelationType):
        """Relations connect entries, identified using the entry IDs."""
        self._relations_by_compound_id[compound_id].add(KeggRelation(from_id, to_id, compound_id, relation_type))

    def match_reaction(self, substrate_names: list[str], product_names: list[str]) -> KeggReactionWithReversion | None:
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

                relation_substrate_kegg_accession = self._entries_by_id[reaction.substrate_id].name
                relation_product_kegg_accession = self._entries_by_id[reaction.product_id].name

                if (_check_for_match([substrate_name], relation_substrate_kegg_accession)
                        and _check_for_match(product_names, relation_product_kegg_accession)):
                    # Found a match!
                    return KeggReactionWithReversion(reaction, reversed=False)

                if (reaction.reaction_type != ReactionType.IRREVERSIBLE and
                        (_check_for_match([substrate_name], relation_product_kegg_accession)
                        and _check_for_match(product_names, relation_substrate_kegg_accession))):
                    # Found a reverse match!
                    return KeggReactionWithReversion(reaction, reversed=True)
        return None

    def _reaction_to_str(self, reaction: KeggReaction) -> str:
        return (_ACCESSION_NUMBER_TO_NAMES[self._entries_by_id[reaction.substrate_id].name][0] + " -> "
                + _ACCESSION_NUMBER_TO_NAMES[self._entries_by_id[reaction.product_id].name][0])

    def entry(self, entry_id: int) -> KeggEntry:
        """Gets an entry from the map with the given id (as added using add_entry). Raises KeyError if not found."""
        return self._entries_by_id[entry_id]

    @property
    def entries(self) -> Collection[KeggEntry]:
        """Gets all available entries."""
        return self._entries_by_id.values()

    @property
    def reactions(self) -> Collection[KeggReaction]:
        """Gets all available reaction."""
        return self._reactions

    @property
    def relations(self) -> Iterable[KeggRelation]:
        """Gets all available relations."""
        for relations in self._relations_by_compound_id.values():
            yield from relations

    def set_reaction_score(self, reaction: KeggReactionWithReversion, score: float):
        """Adds a reaction score to the map with the given id, which can be used for coloring. Raises ValueError
        if a score was already set for the reaction in the given direction, or if the score is NaN."""
        if numpy.isnan(score):
            raise ValueError(f"Score is NaN for {self._reaction_to_str(reaction.reaction)}")

        if reaction.reversed:
            if reaction in self._reactions_backward_values:
                raise ValueError(f"Duplicate backward score for {self._reaction_to_str(reaction.reaction)}")
            self._reactions_backward_values[reaction.reaction.substrate_product_tuple] = score
        else:
            if reaction in self._reactions_forward_values:
                raise ValueError(f"Duplicate forward score for {self._reaction_to_str(reaction.reaction)}")
            self._reactions_forward_values[reaction.reaction.substrate_product_tuple] = score

    def forward_value(self, reaction: KeggReaction) -> float:
        """Gets the forward value for the given reaction, or NaN if not set."""
        return self._reactions_forward_values.get(reaction.substrate_product_tuple, numpy.nan)

    def backward_value(self, reaction: KeggReaction) -> float:
        """Gets the backward value for the given reaction, or NaN if not set."""
        return self._reactions_backward_values.get(reaction.substrate_product_tuple, numpy.nan)

    def has_relations_or_reactions(self, entry: KeggEntry) -> bool:
        """Searches for any reactions for the given entry."""
        entry_id = entry.entry_id
        if entry.entry_type == EntryType.COMPOUND:
            # Search in reactions and ECRels
            for reaction in self._reactions:
                if reaction.substrate_id == entry_id or reaction.product_id == entry_id:
                    return True
            return len(self._relations_by_compound_id[entry.entry_id]) > 0
        elif entry.entry_type == EntryType.GENE or entry.entry_type == EntryType.MAP:
            # Search in relations
            for relation in self.relations:
                if relation.from_id == entry_id or relation.to_id == entry_id:
                    return True
            return False
        elif entry.entry_type == EntryType.TITLE:
            return True  # Title by definition has a relation to the plot, even if not explicitly specified
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

        kegg_map.add_entry(KeggEntry(entry_id, name, x, y, width, height, entry_type))

    # Read reactions
    for reaction in root.findall("reaction"):
        gene_id = int(reaction.attrib["id"])
        substrates = [int(s.attrib['id']) for s in reaction.findall("substrate")]
        products = [int(p.attrib['id']) for p in reaction.findall("product")]
        reversible = reaction.attrib.get("type", "reversible") == "reversible"
        for substrate in substrates:
            for product in products:
                kegg_map.add_reaction(substrate, product, gene_id,
                                      ReactionType.REVERSIBLE if reversible else ReactionType.IRREVERSIBLE)

    # Read relations
    for relation in root.findall("relation"):
        entry1 = int(relation.attrib["entry1"])
        entry2 = int(relation.attrib["entry2"])
        subtype = relation.find("subtype")
        if subtype is None or subtype.attrib["name"] != "compound":
            raise ValueError("Unexpected subtype,", subtype)
        compound_id = int(subtype.attrib["value"])
        relation_type = RelationType[relation.attrib["type"].upper()]
        kegg_map.add_relation(entry1, entry2, compound_id, relation_type)

    return kegg_map
