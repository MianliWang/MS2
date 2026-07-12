import tempfile
import unittest
from pathlib import Path

from MSAI.python.ms1_peak_picker import PeakPickingConfig, PickedPeak, enumerate_peak_candidates
from MSAI.python.ms1_review.classification import (
    background_diagnostics,
    manual_chromatographic_status,
    prediction_diagnostic,
    source_machine_label_interpretation,
)
from PIL import Image

from MSAI.python.ms1_review.exporter import _write_png
from MSAI.python.ms1_review.svg import render_eic_svg


class ReviewClassificationTests(unittest.TestCase):
    def test_candidate_enumerator_keeps_more_than_final_top_two(self):
        rts = [float(index) for index in range(80)]
        signal = [0.0] * 80
        for center in (15, 40, 65):
            signal[center - 1 : center + 2] = [250_000, 500_000, 250_000]
        candidates = enumerate_peak_candidates(
            rts,
            signal,
            PeakPickingConfig(
                sg_window=3,
                min_height=200_000,
                min_prominence_ratio=0.03,
                min_support_scans=2,
                min_distance_sec=3,
            ),
        )
        self.assertEqual(len(candidates), 3)

    def test_manual_status_and_localization_are_separate(self):
        row = {"Peak1": "1", "Peak2": "2"}
        self.assertEqual(manual_chromatographic_status(row), "double_peak")
        self.assertEqual(prediction_diagnostic(row, "double_peak", [1.01, 2.02]), "exact_agreement")
        self.assertEqual(prediction_diagnostic(row, "double_peak", [1.0, 3.0]), "rt_mislocalized")
        self.assertEqual(prediction_diagnostic(row, "single_peak", [1.0]), "class_mismatch")

    def test_low_clean_signal_is_review_only_candidate(self):
        rts = [float(index) for index in range(100)]
        signal = [0.0] * 100
        signal[48:53] = [10, 60_000, 120_000, 60_000, 10]
        result = background_diagnostics(rts, signal, [50 / 60])
        self.assertTrue(result["low_clean_candidate"])

    def test_source_machine_colour_labels_are_not_peak_truth(self):
        self.assertEqual(source_machine_label_interpretation("GREEN"), "expected_double_peak")
        self.assertEqual(source_machine_label_interpretation("yellow"), "expected_single_peak")
        self.assertEqual(source_machine_label_interpretation("RED"), "unusable_no_peak_or_multiple_peak")
        self.assertEqual(source_machine_label_interpretation("CHECK"), "manual_review_needed")


class ReviewSvgTests(unittest.TestCase):
    def test_svg_contains_window_profiles_and_review_markers(self):
        rts = [float(index) for index in range(20)]
        raw = [0.0] * 20
        raw[10] = 1000.0
        peak = PickedPeak(10, 10.0, 1000.0, 2.0, 1000.0, None)
        svg = render_eic_svg(
            compound_id="test<&>", target_mz=300.123456, ppm=3,
            rts_sec=rts, raw=raw, baseline_smooth=raw, experimental_smooth=raw,
            manual_rts_min=[10 / 60], baseline_peaks=(peak,), experimental_peaks=(peak,),
            review_candidates=(),
            metadata={"reference_status":"single_peak","baseline_status":"single_peak","experimental_status":"single_peak","baseline_diagnostic":"exact_agreement","experimental_diagnostic":"exact_agreement","parameter_stability":"stable_class_and_rt","max_intensity":1000,"background_p99":0,"peak_to_background_p99":1000,"max_half_height_support_scans":3,"ms2_status":"not_evaluable","min_height":200000,"baseline_label":"baseline: SG7","experimental_label":"experimental: SG3","baseline_metrics":"P1","experimental_metrics":"P1"},
        )
        self.assertIn("targeted MS1 EIC", svg)
        self.assertIn("baseline: SG7", svg)
        self.assertIn("supplied RT", svg)
        self.assertIn("test&lt;&amp;&gt;", svg)

    def test_png_export_is_a_real_raster_image(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "eic.png"
            _write_png(output, Image.new("RGB", (32, 24), "white"))
            self.assertTrue(output.exists())
            self.assertEqual(output.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")


if __name__ == "__main__":
    unittest.main()
