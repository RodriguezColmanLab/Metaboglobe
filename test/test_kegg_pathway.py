from unittest import TestCase

from metaboglobe.kegg_pathway import KeggMap, KeggEntry, EntryType, RelationType, KeggReactionWithReversion, KeggRelation


def _kegg_compound(name: str, id: int) -> KeggEntry:
    return KeggEntry(name=name, entry_id=id, x=0, y=0, width=1, height=1, entry_type=EntryType.COMPOUND)


class TestKeggMap(TestCase):

    def test_correct_reversability(self):
        """Sometimes KEGG has two irreversible reactions in opposite directions instead of one reversible reaction.
        This test checks that these reactions are corrected to be reversible."""
        kegg_map = KeggMap()

        # Define the four elements
        water = _kegg_compound("C00001", 1)
        atp = _kegg_compound("C00002", 2)
        o2 = _kegg_compound("C00007", 3)
        metal = _kegg_compound("C00050", 4)
        kegg_map.add_entry(water)
        kegg_map.add_entry(atp)
        kegg_map.add_entry(o2)
        kegg_map.add_entry(metal)

        # Define a reaction between atp and o2
        kegg_map.add_relation(2, 3, RelationType.REACTION_IRREVERSIBLE)
        kegg_map.add_relation(3, 2, RelationType.REACTION_IRREVERSIBLE)

        # Check that the reaction that we find is now reversible
        found_reaction = kegg_map.match_reaction(["atp"], ["o2"])
        self.assertEqual(found_reaction.relation, KeggRelation(2, 3, RelationType.REACTION_TWO_IRREVERSIBLE))

        # Also check the other way around
        found_reaction = kegg_map.match_reaction(["o2"], ["atp"])
        self.assertEqual(found_reaction.relation, KeggRelation(2, 3, RelationType.REACTION_TWO_IRREVERSIBLE))