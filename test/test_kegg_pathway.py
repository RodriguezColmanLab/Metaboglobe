from unittest import TestCase

from metaboglobe.kegg_pathway import KeggMap, ReactionType, KeggReactionArrow, KeggReactionIdWithReversion
from metaboglobe.kegg_species import KeggCompoundId, KeggReactionId


class TestKeggMap(TestCase):

    def test_correct_reversability(self):
        kegg_map = KeggMap()

        # Define the four elements
        atp = kegg_map.add_compound(KeggCompoundId.create_from_id("C00002"), 2, 0, 0)
        o2 = kegg_map.add_compound(KeggCompoundId.create_from_id("C00007"), 3, 0, 0)
        catalyst = kegg_map.add_reaction_enzyme(KeggReactionId.create_from_id("R00001"), 5, "ENZYME", 0, 0, 1,1)

        # Define a reaction between atp and o2
        kegg_map.add_reaction_arrow(catalyst, atp, o2, ReactionType.IRREVERSIBLE)

        # Check if we can find back the reaction
        found_reaction = list(kegg_map.match_reactions(["atp"], ["o2"]))
        self.assertEqual(KeggReactionIdWithReversion(KeggReactionId.create_from_id("R00001"), False), found_reaction[0])

        # Also check the other way around (should not be a match)
        found_reaction = list(kegg_map.match_reactions(["o2"], ["atp"]))
        self.assertEqual([], found_reaction)