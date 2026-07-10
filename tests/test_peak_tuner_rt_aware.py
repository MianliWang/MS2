import unittest
from unittest.mock import patch

from MSAI.python.ms1_peak_picker import (
    PeakPickResult,
    PeakPickingConfig,
    PickedPeak,
)
from MSAI.python.tune_peak_picker import (
    _evaluate_rt_aware,
    tune_peak_picker,
)


def _peak(rt_min, intensity=100.0):
    return PickedPeak(1, rt_min * 60, intensity, 5.0, 1000.0, 10.0)


def _result(status, *peaks, valley=0.2, resolution=1.0):
    return PeakPickResult(
        "detected" if peaks else "no_signal",
        status,
        tuple(peaks),
        resolution if len(peaks) == 2 else None,
        valley if len(peaks) == 2 else None,
        max((peak.intensity for peak in peaks), default=0.0),
    )


class RtAwareMetricTests(unittest.TestCase):
    def test_three_classes_require_reviewed_rt_localization(self):
        rows = [
            {"id": "none", "p1": "", "p2": ""},
            {"id": "single", "p1": "2.0", "p2": ""},
            {"id": "close", "p1": "3.0", "p2": "3.1"},
            {"id": "wrong", "p1": "4.0", "p2": "4.4"},
        ]
        picked = {
            "none": _result("not_evaluable"),
            "single": _result("single_peak", _peak(2.05)),
            "close": _result("double_peak", _peak(3.0), _peak(3.1)),
            # Right class, wrong chromatographic location: this must not be a TP.
            "wrong": _result("double_peak", _peak(8.0), _peak(8.4)),
        }

        metrics = _evaluate_rt_aware(
            rows,
            picked,
            "id",
            "p1",
            "p2",
            0.95,
            0.0,
            0.02,
        )

        self.assertAlmostEqual(metrics.accuracy, 0.75)
        self.assertAlmostEqual(metrics.macro_recall, (1.0 + 1.0 + 0.5) / 3)
        self.assertEqual(metrics.mislocalized, 1)
        self.assertEqual(metrics.close_double_rows, 1)
        self.assertEqual(metrics.close_double_recall, 1.0)


class AdaptiveTuningIsolationTests(unittest.TestCase):
    def test_explicit_no_peak_is_used_and_test_is_picked_only_after_selection(self):
        rows = [
            {"id": "cal_none", "mz": "100", "p1": "", "p2": "", "review": "RED"},
            {"id": "cal_single", "mz": "101", "p1": "1", "p2": "", "review": "YELLOW"},
            {"id": "cal_double", "mz": "102", "p1": "1", "p2": "1.1", "review": "GREEN"},
            {"id": "test_single", "mz": "103", "p1": "1", "p2": "", "review": "YELLOW"},
        ]
        traces = {float(value): [float(value)] for value in range(100, 104)}
        calls = []

        def fake_pick(_rts, trace, config):
            marker = int(trace[0])
            calls.append((marker, config))
            if marker == 100:
                return _result("not_evaluable")
            if marker in (101, 103):
                return _result("single_peak", _peak(1.0))
            return _result("double_peak", _peak(1.0), _peak(1.1), valley=0.9)

        core = PeakPickingConfig(
            min_height=None,
            min_prominence_ratio=0.25,
            min_support_scans=2,
            min_distance_sec=6.0,
            close_peak_max_valley_ratio=0.9,
            rank_by_prominence=True,
        )
        with (
            patch("MSAI.python.tune_peak_picker.read_table", return_value=rows),
            patch(
                "MSAI.python.tune_peak_picker.extract_target_eics_with_config",
                return_value=([60.0], traces),
            ),
            patch("MSAI.python.tune_peak_picker._split", side_effect=lambda value: "test" if value.startswith("test") else "calibration"),
            patch("MSAI.python.tune_peak_picker._adaptive_core_grid", return_value=[core]),
            patch("MSAI.python.tune_peak_picker._adaptive_quality_grid", return_value=[(0.95, 0.0, 0.02)]),
            patch("MSAI.python.tune_peak_picker.pick_chiral_peaks", side_effect=fake_pick),
        ):
            report = tune_peak_picker(
                "labels.csv",
                "raw.mzML",
                mz_column="mz",
                peak1_column="p1",
                peak2_column="p2",
                id_column="id",
                adaptive=True,
                review_column="review",
            )

        self.assertEqual([marker for marker, _config in calls], [100, 101, 102, 103])
        self.assertEqual(report["calibration_metrics"]["no_peak_rows"], 1)
        self.assertEqual(report["test_metrics"]["single_peak_recall"], 1.0)
        # Only the holdout call receives the fully selected quality configuration.
        self.assertIsNone(calls[0][1].max_valley_ratio)
        self.assertEqual(calls[-1][1].max_valley_ratio, 0.95)


if __name__ == "__main__":
    unittest.main()
