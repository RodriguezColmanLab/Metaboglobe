from unittest import TestCase

import numpy

from metaboglobe.kegg_pathway import KeggMap, ReactionType, KeggReactionArrow
from metaboglobe.kegg_species import KeggCompoundId, KeggReactionId, KeggReactionIdWithReversion


class TestKeggMap(TestCase):

    def test_reaction_matching(self):
        kegg_map = KeggMap()

        # Define the four elements
        atp = kegg_map.add_compound(KeggCompoundId.create_from_id("C00002"), 2, 0, 0)
        o2 = kegg_map.add_compound(KeggCompoundId.create_from_id("C00007"), 3, 0, 0)
        catalyst = kegg_map.add_reaction_enzyme(KeggReactionId.create_from_id("R00001"), 5, "ENZYME", 0, 0, 1,1)

        # Define a reaction between atp and o2
        kegg_map.add_reaction_arrow(catalyst, atp, o2, ReactionType.IRREVERSIBLE)

        # Check if we can find back the reaction
        found_reaction = list(kegg_map.match_reactions(substrate_names=["atp"], product_names=["o2"]))
        self.assertEqual(KeggReactionIdWithReversion(KeggReactionId.create_from_id("R00001"), False), found_reaction[0])

        # Also check the other way around (should not be a match)
        found_reaction = list(kegg_map.match_reactions(substrate_names=["o2"], product_names=["atp"]))
        self.assertEqual([], found_reaction)

    def test_set_reaction_score(self):
        kegg_map = KeggMap()
        reaction_id = KeggReactionId.create_from_id("R00001")

        # Define the four elements
        atp = kegg_map.add_compound(KeggCompoundId.create_from_id("C00002"), 2, 0, 0)
        o2 = kegg_map.add_compound(KeggCompoundId.create_from_id("C00007"), 3, 0, 0)
        catalyst = kegg_map.add_reaction_enzyme(reaction_id, 5, "ENZYME", 0, 0, 1, 1)

        # Define a reaction between atp and o2
        kegg_map.add_reaction_arrow(catalyst, atp, o2, ReactionType.IRREVERSIBLE)

        # Test setting the forwards score
        self.assertTrue(numpy.isnan(kegg_map.get_reaction_score(reaction_id.forwards())))
        kegg_map.set_reaction_score(reaction_id.forwards(), 1.5)
        self.assertEqual(1.5, kegg_map.get_reaction_score(reaction_id.forwards()))
        self.assertTrue(numpy.isnan(kegg_map.get_reaction_score(reaction_id.backwards())))

    def test_set_compound_score(self):
        kegg_map = KeggMap()
        atp_id = KeggCompoundId.create_from_id("C00002")
        kegg_map.add_compound(atp_id, 2, 0, 0)

        # Test setting a score
        self.assertTrue(numpy.isnan(kegg_map.get_compound_score(atp_id)))
        kegg_map.set_compound_score(atp_id, 1.5)
        self.assertEqual(1.5, kegg_map.get_compound_score(atp_id))
