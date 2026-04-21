from unittest import TestCase

from metaboglobe.tabular_data import parse_reaction


class TestTabularData(TestCase):

    def test_parse_reaction(self):
        self.assertEqual((["ATP", "Coenzyme A", "acetate"], ["AMP", "Diphosphate", "Acetyl-CoA"]),
            parse_reaction("1.00 * ATP [c] + 1.00 * Coenzyme A [c] + 1.00 * acetate [c] --> 1.00 * AMP [c] + 1.00 * Diphosphate [c] + 1.00 * Acetyl-CoA [c] AACS; ACSS2"))

    def test_missing_arrow(self):
        self.assertRaises(ValueError, parse_reaction, "1.00 * ATP")

    def test_missing_counts(self):
        self.assertRaises(ValueError, parse_reaction, "ATP [c] + Coenzyme A [c] + acetate [c] --> AMP [c] + Diphosphate [c] + Acetyl-CoA [c] AACS; ACSS2")

    def test_invalid_counts(self):
        self.assertRaises(ValueError, parse_reaction, "one * ATP [c] + 1.00 * Coenzyme A [c] + 1.00 * acetate [c] --> 1.00 * AMP [c] + 1.00 * Diphosphate [c] + 1.00 * Acetyl-CoA [c] AACS; ACSS2")