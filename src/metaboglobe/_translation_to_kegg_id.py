from collections import defaultdict
from importlib.resources import files

from metaboglobe.kegg_species import KeggReactionId


def _read_gem_to_kegg_reaction_ids() -> dict[str, list[KeggReactionId]]:
    """Reads the mapping from KEGG accession numbers to names from the "Kegg Metabolism" folder. In the returned mapping,
    the main name is always the first one in the list, and the synonyms are the following ones."""
    file_path = files("metaboglobe.data") / "gem_reactions.tsv"

    mapping = defaultdict(list)
    with file_path.open("r") as handle:
        line = handle.readline()
        while line.startswith(";"):
            line = handle.readline()

        # Now we have the first line, read in the header
        headers = line.split("\t")
        gem_reaction_column = headers.index("rxns")
        kegg_reaction_column = headers.index("rxnKEGGID")
        line = handle.readline()

        while line != "":
            if not line.startswith(";"):
                line_split = line.split("\t")
                gem_reaction_id = line_split[gem_reaction_column]
                kegg_reaction_id = line_split[kegg_reaction_column]
                if kegg_reaction_id != "":
                    mapping[gem_reaction_id].append(KeggReactionId.create_from_id(kegg_reaction_id))

            line = handle.readline()
    return mapping


def _read_ecrel_to_kegg_reaction_ids() -> dict[str, list[KeggReactionId]]:
    file_path = files("metaboglobe.data") / "kegg_enzymes.tsv"
    mapping = defaultdict(list)
    with file_path.open("r") as handle:
        line = handle.readline()

        while line != "":
            if not line.startswith(";"):
                ecrel, kegg_id = line.split("\t")
                ecrel = ecrel.strip()
                kegg_id = kegg_id.strip()

                if not ecrel.startswith("ec:"):
                    raise ValueError(f"Expected EC number to start with 'ec:', but got '{ecrel}'")
                if not kegg_id.startswith("rn:"):
                    raise ValueError(f"Expected KEGG reaction ID to start with 'rn:', but got '{kegg_id}'")
                mapping[ecrel[3:]].append(KeggReactionId(kegg_id))
            line = handle.readline()
    return mapping


_GEM_TO_KEGG_REACTION_MAPPING = None
def map_gem_id_to_kegg_reactions(gem_id: str) -> list[KeggReactionId]:
    """Maps an id in the form of MAR04650 to a KEGG reaction ID, for example "rn:R04650". Returns None if no mapping is
    found."""
    global _GEM_TO_KEGG_REACTION_MAPPING
    if _GEM_TO_KEGG_REACTION_MAPPING is None:
        _GEM_TO_KEGG_REACTION_MAPPING = _read_gem_to_kegg_reaction_ids()
    return _GEM_TO_KEGG_REACTION_MAPPING[gem_id]


_ECREL_TO_KEGG_REACTION_MAPPING = None
def map_ecrel_to_kegg_reactions(ecrel: str) -> list[KeggReactionId]:
    """Maps an EC rel in the form of "2.6.1.33" (so without prefix) to a KEGG reaction ID, for example "rn:R04650".
    Returns None if no mapping is found."""
    global _ECREL_TO_KEGG_REACTION_MAPPING
    if _ECREL_TO_KEGG_REACTION_MAPPING is None:
        _ECREL_TO_KEGG_REACTION_MAPPING = _read_ecrel_to_kegg_reaction_ids()
    return _ECREL_TO_KEGG_REACTION_MAPPING[ecrel]
