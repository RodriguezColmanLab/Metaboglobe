import re

import pandas

from metaboglobe import _translation_to_kegg_id
from metaboglobe.kegg_pathway import KeggMap, KeggReactionIdWithReversion
from metaboglobe.kegg_species import KeggReactionId

_ECREL_DETECTION = re.compile(r"(ec:)?([0-9]+\.)+[0-9]+")


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


def insert_values_in_map(kegg_map: KeggMap, data_frame: pandas.DataFrame, *, reaction_identifier_column: str, value_col: str,
                         reversed_col: str | None = None) -> None:
    """Reads values from the given data frame into the given metabolic map. The dataframe must have a column with
    reactions, and a column with values. Optionally, it can also have a column that indicates whether the
    reaction must be reversed or not. (If reversed is true, the corresponding value is registered for the inverse
    reaction instead.)

    The reactions can be:

    - a reaction equation the form

      1.00 * ATP [c] + 1.00 * Coenzyme A [c] + 1.00 * acetate [c] --> 1.00 * AMP [c] + 1.00 * Diphosphate [c] + 1.00 * Acetyl-CoA [c]

      Everything after the first square bracket [ is ignored, so the above equation is equivalent to:

      1 * ATP + 1 * Coenzyme A + 1 * acetate --> 1 * AMP + 1 * Diphosphate + 1 * Acetyl-CoA

      The "-->", "+" and "*" are mandatory though.
      Warning: enzymes are ignored here in the matching, which will be an issue if multiple reactions share the same
      substrate and product in the same map. In that case, the value will be registered for the all reactions that
      matches the equation.
    - an ECNumber identifier: "ec:6.2.1.50". The "ec:" prefix is optional.
    - a KEGG Reaction ID: "rn:R11608". The "rn:" prefix is optional.
    - a GEM reaction ID: "MAR08691". The "MAR" prefix is *mandatory*.
    """

    reactions = data_frame[reaction_identifier_column]
    values = data_frame[value_col]
    if reversed_col is None:
        reversions = pandas.Series(False, index=reactions.index)  # Just assume no reactions are reversed
    else:
        reversions = data_frame[reversed_col]

    for reaction, reversed, value in zip(reactions, reversions, values):
        if not isinstance(reaction, str):
            continue  # Likely NA or something like that

        # Reaction equation
        if "-->" in reaction:
            substrate_names, product_names = parse_reaction(reaction)
            if reversed:
                product_names, substrate_names = substrate_names, product_names
            matched_reactions = list(kegg_map.match_reactions(substrate_names, product_names))
            if len(matched_reactions) == 0 and "L-glutamate(1-)" in reaction and "4-Aminobutanoate" in reaction:
                print(f"No match for {reaction}")
                matched_reactions = list(kegg_map.match_reactions(substrate_names, product_names))
            for matched_reaction in matched_reactions:
                kegg_map.set_reaction_score(matched_reaction, value)

        # ECNums
        elif _ECREL_DETECTION.fullmatch(reaction):
            if reaction.startswith("ec:"):
                reaction = reaction[3:]
            for matched_reaction in _translation_to_kegg_id.map_ecrel_to_kegg_reactions(reaction):
                kegg_map.set_reaction_score(KeggReactionIdWithReversion(matched_reaction, reversed), value)

        # KEGG reaction IDs
        elif reaction.startswith("R"):
            matched_reaction = KeggReactionId.create_from_id(reaction)
            kegg_map.set_reaction_score(KeggReactionIdWithReversion(matched_reaction, reversed), value)
        elif reaction.startswith("rn:"):
            matched_reaction = KeggReactionId(reaction)
            kegg_map.set_reaction_score(KeggReactionIdWithReversion(matched_reaction, reversed), value)

        # GEM IDs
        elif reaction.startswith("MAR"):
            for matched_reaction in _translation_to_kegg_id.map_gem_id_to_kegg_reactions(reaction):
                kegg_map.set_reaction_score(KeggReactionIdWithReversion(matched_reaction, reversed), value)
