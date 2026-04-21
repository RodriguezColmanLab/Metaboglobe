import pandas

from metaboglobe.kegg_pathway import KeggMap


def parse_reaction(reaction: str) -> tuple[list[str], list[str]]:
    """Parses reactions in the format:

    1.00 * ATP [c] + 1.00 * Coenzyme A [c] + 1.00 * acetate [c] --> 1.00 * AMP [c] + 1.00 * Diphosphate [c] + 1.00 * Acetyl-CoA [c] AACS; ACSS2

    The plusses are essential to separate the compound names, as well as the arrow "-->". Numbers at the start of a
     compound are ignored. Any part after the square brackets are ignored as well, these tend to be enzymes.
    """
    if "-->" not in reaction:
        raise ValueError(f"Invalid reaction, missing -->: '{reaction}'")

    substrates_str, products_str = reaction.split("-->")
    substrates = _parse_compounds(reaction, substrates_str)
    products = _parse_compounds(reaction, products_str)
    return substrates, products


def _parse_compounds(full_reaction_equation: str, compounds_str: str):
    compounds = list()
    for compound_str in compounds_str.split("+"):
        if "*" not in compound_str:
            raise ValueError(f"Invalid substrate, missing '*' for '{compound_str}' in reaction '{full_reaction_equation}'")
        molecule_count, molecular_formula = compound_str.split("*")

        # The count we just validate, we don't use it
        try:
            float(molecule_count.strip())
        except ValueError:
            raise ValueError(f"Invalid molecule count for '{compound_str}' in reaction '{full_reaction_equation}'")

        if "[" in molecular_formula:
            molecular_formula = molecular_formula[:molecular_formula.index("[")].strip()

        compounds.append(molecular_formula)
    return compounds


def insert_values_in_map(kegg_map: KeggMap, data_frame: pandas.DataFrame, *, reaction_col: str, value_col: str,
                         reversed_col: str | None = None) -> None:

    reactions = data_frame[reaction_col]
    values = data_frame[value_col]
    if reversed_col is None:
        reversions = pandas.Series(False, index=reactions.index)
    else:
        reversions = data_frame[reversed_col]

    for reaction, reversed, value in zip(reactions, reversions, values):
        substrate_names, product_names = parse_reaction(reaction)
        if reversed:
            product_names, substrate_names = substrate_names, product_names
        matched_reaction = kegg_map.match_reaction(substrate_names, product_names)
        if matched_reaction is not None:
            kegg_map.set_reaction_score(matched_reaction, value)
