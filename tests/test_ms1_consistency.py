import unittest
from unittest.mock import patch

from MSAI.python.ms1_peak_picker import (
    chromatographic_resolution_fwhm,
    PeakPickResult,
    PeakPickingConfig,
    PickedPeak,
    annotate_peaklist,
    resolution_passes_threshold,
)
from MSAI.python.tune_peak_picker import _call


class Ms1OutputConsistencyTests(unittest.TestCase):
    def test_fwhm_resolution_uses_half_height_conversion(self):
        self.assertAlmostEqual(
            chromatographic_resolution_fwhm(0.0, 10.0, 2.0, 3.0),
            1.17741 * 2,
        )

    def test_annotate_peaklist_outputs_snr_for_each_automatic_peak(self):
        rows = [{"MZ": "100.0"}]
        result = PeakPickResult(
            signal_status="detected",
            chromatographic_status="double_peak",
            peaks=(
                PickedPeak(10, 60.0, 200_000.0, 4.0, 500_000.0, 12.5),
                PickedPeak(20, 120.0, 150_000.0, 5.0, 400_000.0, 8.25),
            ),
            resolution=2.0,
            valley_ratio=0.2,
            max_intensity=200_000.0,
        )

        with (
            patch("MSAI.python.ms1_peak_picker.read_table", return_value=rows),
            patch(
                "MSAI.python.ms1_peak_picker.extract_target_eics_with_config",
                return_value=([60.0], {100.0: [200_000.0]}),
            ),
            patch("MSAI.python.ms1_peak_picker.pick_chiral_peaks", return_value=result),
            patch("MSAI.python.ms1_peak_picker.write_table") as write_table,
        ):
            annotated = annotate_peaklist("peaks.csv", "raw.mzML", "output.csv")

        self.assertEqual(annotated[0]["auto_peak1_snr"], "12.5")
        self.assertEqual(annotated[0]["auto_peak2_snr"], "8.25")
        self.assertEqual(write_table.call_args.args[1], annotated)


class Ms1ResolutionThresholdConsistencyTests(unittest.TestCase):
    @staticmethod
    def _result_without_width_resolution():
        return PeakPickResult(
            signal_status="detected",
            chromatographic_status="double_peak",
            peaks=(
                PickedPeak(10, 60.0, 200.0, None, None, 10.0),
                PickedPeak(20, 120.0, 150.0, None, None, 8.0),
            ),
            resolution=None,
            valley_ratio=0.2,
            max_intensity=200.0,
        )

    def test_zero_and_none_both_disable_resolution_requirement(self):
        result = self._result_without_width_resolution()
        for minimum in (0.0, None):
            with self.subTest(min_resolution=minimum):
                config = PeakPickingConfig(min_resolution=minimum)
                self.assertEqual(config.min_resolution, minimum)
                self.assertTrue(resolution_passes_threshold(result.resolution, minimum))
                self.assertTrue(_call(result, 0.7, minimum, 0.1))

    def test_positive_resolution_requirement_still_rejects_missing_resolution(self):
        result = self._result_without_width_resolution()
        self.assertFalse(resolution_passes_threshold(result.resolution, 0.5))
        self.assertFalse(_call(result, 0.7, 0.5, 0.1))


if __name__ == "__main__":
    unittest.main()
