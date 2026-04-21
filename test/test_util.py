from unittest import TestCase

from metaboglobe._util import get_names_without_stereoisomers


class Test(TestCase):
    def test_get_name_without_enantiomers(self):
        self.assertEqual(set(get_names_without_stereoisomers("glucose")), set())

        self.assertEqual(set(get_names_without_stereoisomers("D-fructose")), {"fructose"})

        self.assertEqual(set(get_names_without_stereoisomers("alpha-D-glucose-6-phosphate")),
                         {"alpha-glucose-6-phosphate", "D-glucose-6-phosphate", "glucose-6-phosphate"})

        self.assertEqual(set(get_names_without_stereoisomers("2-(alpha-Hydroxyethyl)thiamine-diphosphate")),
                         {"2-(Hydroxyethyl)thiamine-diphosphate"})
