import math
import unittest

from MSAI.python.easms_statistics import (
    benjamini_hochberg,
    enantiomer_fraction,
    enantiomer_fraction_shift,
    enantiomer_log2_ratio_shift,
    enrichment_fold,
    input_eluate_recovery_fold,
    target_control_enrichment_fold,
)


class EffectSizeSemanticsTests(unittest.TestCase):
    def test_old_enrichment_name_is_compatible_alias(self):
        self.assertIs(enrichment_fold, target_control_enrichment_fold)
        self.assertEqual(target_control_enrichment_fold(10, [1, 2, 3]), 5.0)

    def test_target_control_and_input_eluate_folds_have_explicit_direction(self):
        self.assertEqual(target_control_enrichment_fold(9, [1, 3, 5]), 3.0)
        self.assertEqual(input_eluate_recovery_fold(4, 10), 2.5)
        self.assertEqual(input_eluate_recovery_fold(0, 1, pseudocount=1), 2.0)

    def test_fraction_and_fraction_shift(self):
        self.assertEqual(enantiomer_fraction(3, 1), 0.75)
        positive_shift = enantiomer_fraction_shift(1, 1, 3, 1)
        negative_shift = enantiomer_fraction_shift(3, 1, 1, 1)
        assert positive_shift is not None
        assert negative_shift is not None
        self.assertAlmostEqual(positive_shift, 0.25)
        self.assertAlmostEqual(negative_shift, -0.25)

    def test_zero_denominators_are_undefined_not_infinite(self):
        self.assertIsNone(target_control_enrichment_fold(1, [0, 0]))
        self.assertIsNone(target_control_enrichment_fold(1, []))
        self.assertIsNone(input_eluate_recovery_fold(0, 1))
        self.assertIsNone(enantiomer_fraction(0, 0))
        self.assertIsNone(enantiomer_fraction_shift(0, 0, 1, 1))

    def test_invalid_intensities_and_pseudocounts_are_rejected(self):
        invalid_values = (-1, math.nan, math.inf, -math.inf, "not-a-number")
        for invalid in invalid_values:
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    target_control_enrichment_fold(invalid, [1])  # pyright: ignore[reportArgumentType]
                with self.assertRaises(ValueError):
                    target_control_enrichment_fold(1, [invalid])
                with self.assertRaises(ValueError):
                    enantiomer_fraction(invalid, 1)  # pyright: ignore[reportArgumentType]
        for invalid_pseudocount in (-1, math.nan, math.inf):
            with self.subTest(pseudocount=invalid_pseudocount), self.assertRaises(ValueError):
                input_eluate_recovery_fold(1, 2, invalid_pseudocount)

    def test_log_ratio_requires_positive_pseudocount_and_valid_areas(self):
        self.assertAlmostEqual(enantiomer_log2_ratio_shift(1, 1, 3, 1, 1), 1.0)
        with self.assertRaises(ValueError):
            enantiomer_log2_ratio_shift(1, 1, 3, 1, 0)
        with self.assertRaises(ValueError):
            enantiomer_log2_ratio_shift(1, -1, 3, 1, 1)

    def test_fraction_is_stable_for_large_finite_areas(self):
        self.assertEqual(enantiomer_fraction(1e308, 1e308), 0.5)


class BenjaminiHochbergEdgeTests(unittest.TestCase):
    def test_ties_receive_equal_q_values(self):
        q_values = benjamini_hochberg([0.01, 0.01, 0.5])
        self.assertEqual(q_values, [0.015, 0.015, 0.5])

    def test_all_missing_and_empty_inputs(self):
        self.assertEqual(benjamini_hochberg([None, None]), [None, None])
        self.assertEqual(benjamini_hochberg([]), [])

    def test_order_and_missing_entries_are_preserved(self):
        self.assertEqual(
            benjamini_hochberg([0.04, None, 0.01, 0.03]),
            [0.04, None, 0.03, 0.04],
        )

    def test_invalid_p_values_are_rejected(self):
        for invalid in (-0.01, 1.01, math.nan, math.inf, -math.inf, "bad"):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                benjamini_hochberg([0.1, invalid])


if __name__ == "__main__":
    unittest.main()
