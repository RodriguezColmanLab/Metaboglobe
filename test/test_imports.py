from unittest import TestCase


class TestImports(TestCase):

    def test_imports(self):
        """Tests the lazy loading of imports. Can we just use `import metaboglobe` instead of
        `import metaboglobe.tabular_data`?"""
        import metaboglobe
        self.assertEqual((["ATP", "Coenzyme A", "acetate"], ["AMP", "Diphosphate", "Acetyl-CoA"]),
            metaboglobe.tabular_data.parse_reaction("1.00 * ATP [c] + 1.00 * Coenzyme A [c] + 1.00 * acetate [c] --> 1.00 * AMP [c] + 1.00 * Diphosphate [c] + 1.00 * Acetyl-CoA [c] AACS; ACSS2"))
