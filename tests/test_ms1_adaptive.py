import math
import unittest
from unittest.mock import patch

from MSAI.python.ms1_peak_picker import PeakPickingConfig, pick_chiral_peaks


class AdaptiveMs1PeakPickingTests(unittest.TestCase):
    def test_clean_low_peak_does_not_need_an_absolute_height_cutoff(self):
        rts = [float(index) for index in range(41)]
        intensities = [
            8.0 + 42.0 * math.exp(-((index - 20) / 3.0) ** 2)
            for index in range(41)
        ]
        result = pick_chiral_peaks(
            rts,
            intensities,
            PeakPickingConfig(
                min_height=None,
                min_prominence_ratio=0.8,
                min_support_scans=3,
                min_distance_scans=5,
            ),
        )

        self.assertEqual(result.chromatographic_status, "single_peak")
        self.assertEqual(result.peaks[0].scan_index, 20)
        self.assertLess(result.peaks[0].intensity, 100.0)
        self.assertGreater(result.peaks[0].prominence_ratio or 0.0, 0.8)
        self.assertIn("absolute_height_disabled", result.review_reasons)

    def test_clear_valley_preserves_close_peaks_and_bounds_raw_refinement(self):
        rts = [float(index) for index in range(25)]
        intensities = [0.0] * 25
        intensities[8:17] = [1.0, 6.0, 10.0, 6.0, 2.0, 7.0, 9.0, 5.0, 1.0]
        result = pick_chiral_peaks(
            rts,
            intensities,
            PeakPickingConfig(
                sg_window=3,
                sg_polyorder=2,
                min_height=None,
                min_prominence_ratio=0.7,
                min_support_scans=2,
                min_distance_scans=10,
                close_peak_max_valley_ratio=0.3,
                refine_radius_scans=5,
            ),
        )

        self.assertEqual(result.chromatographic_status, "double_peak")
        self.assertEqual([peak.scan_index for peak in result.peaks], [10, 14])
        self.assertIn("close_peak_pair", result.review_reasons)

    def test_single_scan_spike_fails_continuous_support(self):
        rts = [float(index) for index in range(21)]
        intensities = [0.0] * 21
        intensities[10] = 20.0
        result = pick_chiral_peaks(
            rts,
            intensities,
            PeakPickingConfig(
                sg_window=3,
                sg_polyorder=2,
                min_height=None,
                min_prominence_ratio=0.5,
                min_support_scans=2,
            ),
        )

        self.assertEqual(result.peaks, ())
        self.assertEqual(result.chromatographic_status, "ambiguous")

    def test_default_configuration_stays_on_legacy_detector(self):
        rts = [float(index) for index in range(100)]
        intensities = [
            220_000 * math.exp(-((index - 25) / 4) ** 2)
            + 180_000 * math.exp(-((index - 70) / 5) ** 2)
            for index in range(100)
        ]
        with patch(
            "MSAI.python.ms1_peak_picker._find_adaptive_peaks",
            side_effect=AssertionError("legacy config entered adaptive detector"),
        ):
            result = pick_chiral_peaks(rts, intensities, PeakPickingConfig())

        self.assertEqual(result.chromatographic_status, "double_peak")
        self.assertEqual([peak.scan_index for peak in result.peaks], [25, 70])
        self.assertTrue(all(peak.prominence is None for peak in result.peaks))

    def test_second_based_distance_uses_irregular_retention_times(self):
        rts = [float(index) for index in range(11)] + [30.0, 31.0, 32.0, 33.0, 34.0]
        intensities = [0.0] * len(rts)
        intensities[8:15] = [1.0, 5.0, 10.0, 2.0, 5.0, 9.0, 1.0]
        result = pick_chiral_peaks(
            rts,
            intensities,
            PeakPickingConfig(
                sg_window=3,
                sg_polyorder=2,
                min_height=None,
                min_prominence_ratio=0.5,
                min_support_scans=1,
                min_distance_scans=20,
                min_distance_sec=10.0,
                refine_radius_scans=0,
            ),
        )

        self.assertEqual([peak.scan_index for peak in result.peaks], [10, 13])


if __name__ == "__main__":
    unittest.main()
