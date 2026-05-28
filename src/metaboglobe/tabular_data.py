import re
import numpy
import pandas

from metaboglobe import _translation_to_kegg_id
from metaboglobe.kegg_pathway import KeggMap
from metaboglobe.kegg_species import KeggReactionId, KeggReactionIdWithReversion

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


def insert_values_in_map(kegg_map: KeggMap, data_frame: pandas.DataFrame, *,
                         reaction_identifier_col: str | None = None,
                         value_col: str,
                         reversed_col: str | None = None,
                         enzyme_col: str | None = None,
                         enzyme_col_sep: str = ";") -> None:
    """Reads values from the given data frame into the given metabolic map. The dataframe must have a column with
    reactions and/or enzymes, and a column with values. Optionally, it can also have a column that indicates whether the
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
      matches the equation. You can avoid this by using the enzyme_col parameter.
    - an ECNumber identifier: "ec:6.2.1.50". The "ec:" prefix is optional.
    - a KEGG Reaction ID: "rn:R11608". The "rn:" prefix is optional.
    - a GEM reaction ID: "MAR08691". The "MAR" prefix is *mandatory*.

    The enzymes can be:
    - a string, separated by enzyme_col_sep (";" by default), representing gene names.
    - a list of strings, representing gene names.

    The enzyme_col parameter is only used if reaction_identifier_col is None (then we match solely by enzyme),
    or if it points to a reaction equation (then we further narrow down reaction equations by enzyme).
    However, if reaction_identifier_col  points to a reaction id, the enzyme_col is ignored, as the reaction is already
    fully defined.

    When matching solely by enzyme, multiple enzymes can point to the same reaction. In that case, the absolute highest
    value is used. So if one enzyme has a value of -1.5, and the other of +1, but they belong to the same reaction, only
    the -1.5 is plotted.
    """
    kegg_map.clear_reaction_scores()

    reactions_all = data_frame[reaction_identifier_col] if reaction_identifier_col is not None else None
    values_all = data_frame[value_col]
    if reversed_col is None:
        reversions_all = pandas.Series(False, index=values_all.index)  # Just assume no reactions are reversed
    else:
        reversions_all = data_frame[reversed_col]
    enzymes_all = data_frame[enzyme_col] if enzyme_col is not None else None

    # Match by enzymes
    if reactions_all is None:
        if enzymes_all is None:
            raise ValueError("Reactions DataFrame must have reaction_identifier_col or enzyme_col")
        for i in range(len(enzymes_all)):
            reaction_value = values_all[i]
            reaction_reversed = reversions_all[i]
            enzyme_names = None
            enzymes_str = enzymes_all[i] if enzymes_all is not None else None
            if isinstance(enzymes_str, str):
                enzyme_names = enzymes_str.split(enzyme_col_sep)
            elif isinstance(enzymes_str, list):
                enzyme_names = enzymes_str
            for matched_reaction_id in kegg_map.match_reactions_for_enzymes(enzyme_names):
                matched_reaction = matched_reaction_id.backwards() if reversed_col else matched_reaction_id.forwards()
                existing_value = kegg_map.get_reaction_score(matched_reaction)
                kegg_map.set_reaction_score(matched_reaction, _nan_abs_max(existing_value, reaction_value))
        return

    # Match by reaction identifiers
    for i in range(len(reactions_all)):
        reaction = reactions_all[i]
        if not isinstance(reaction, str):
            continue  # Likely NA or something like that

        reaction_value = values_all[i]
        reaction_reversed = reversions_all[i]

        # Reaction equation
        if "-->" in reaction:
            enzyme_names = None
            enzymes_str = enzymes_all[i] if enzymes_all is not None else None
            if isinstance(enzymes_str, str):
                enzyme_names = enzymes_str.split(enzyme_col_sep)
            elif isinstance(enzymes_str, list):
                enzyme_names = enzymes_str

            substrate_names, product_names = parse_reaction(reaction)
            if reaction_reversed:
                product_names, substrate_names = substrate_names, product_names
            matched_reactions = list(kegg_map.match_reactions(substrate_names=substrate_names, product_names=product_names, enzyme_names=enzyme_names))
            for matched_reaction in matched_reactions:
                kegg_map.set_reaction_score(matched_reaction, reaction_value)

        # ECNums
        elif _ECREL_DETECTION.fullmatch(reaction):
            if reaction.startswith("ec:"):
                reaction = reaction[3:]
            for matched_reaction in _translation_to_kegg_id.map_ecrel_to_kegg_reactions(reaction):
                kegg_map.set_reaction_score(KeggReactionIdWithReversion(matched_reaction, reaction_reversed), reaction_value)

        # KEGG reaction IDs
        elif reaction.startswith("R"):
            matched_reaction = KeggReactionId.create_from_id(reaction)
            kegg_map.set_reaction_score(KeggReactionIdWithReversion(matched_reaction, reaction_reversed), reaction_value)
        elif reaction.startswith("rn:"):
            matched_reaction = KeggReactionId(reaction)
            kegg_map.set_reaction_score(KeggReactionIdWithReversion(matched_reaction, reaction_reversed), reaction_value)

        # GEM IDs
        elif reaction.startswith("MAR"):
            for matched_reaction in _translation_to_kegg_id.map_gem_id_to_kegg_reactions(reaction):
                kegg_map.set_reaction_score(KeggReactionIdWithReversion(matched_reaction, reaction_reversed), reaction_value)


def _nan_abs_max(number1: float, number2: float) -> float:
    if numpy.isnan(number1):
        return number2
    if numpy.isnan(number2):
        return number1
    if abs(number1) > abs(number2):
        return number1
    else:
        return number2