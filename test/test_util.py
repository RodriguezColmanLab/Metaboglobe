from unittest import TestCase

from metaboglobe._util import get_name_without_enantiomers


class Test(TestCase):
    def test_get_name_without_enantiomers(self):
        self.assertEqual(get_name_without_enantiomers("alpha-D-glucose-6-phosphate"), "glucose-6-phosphate")
        self.assertEqual(get_name_without_enantiomers("glucose"), "glucose")
        self.assertEqual(get_name_without_enantiomers("D-fructose"), "fructose")
        self.assertEqual(get_name_without_enantiomers("2-(alpha-Hydroxyethyl)thiamine-diphosphate"), "2-(Hydroxyethyl)thiamine-diphosphate")
