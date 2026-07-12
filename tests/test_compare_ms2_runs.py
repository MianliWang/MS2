import csv
import json
import tempfile
import unittest
from pathlib import Path

from MSAI.python.compare_ms2_runs import compare_ms2_runs


class CompareMs2RunsTests(unittest.TestCase):
    def test_status_change_queue_keeps_candidate_and_baseline_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fields = [
                "Compound_ID",
                "SGC ID for Component",
                "SGC ID for Pool",
                "MZ",
                "Peak1",
                "Peak2",
                "peak_a_MS2",
                "peak_b_MS2",
                "ms2_diagnostic_status",
                "ms2_issue_codes",
                "ms2_cosine",
                "ms2_entropy_similarity",
                "ms2_matched_peaks",
                "peak_a_explained_intensity",
                "peak_b_explained_intensity",
            ]
            baseline = root / "baseline.csv"
            candidate = root / "candidate.csv"
            common = {
                "Compound_ID": "A",
                "SGC ID for Component": "A",
                "SGC ID for Pool": "P",
                "MZ": "300",
                "Peak1": "1",
                "Peak2": "2",
                "ms2_issue_codes": "",
                "ms2_entropy_similarity": "0.9",
                "ms2_matched_peaks": "6",
                "peak_a_explained_intensity": "0.9",
                "peak_b_explained_intensity": "0.9",
            }
            self._write(
                baseline,
                fields,
                {
                    **common,
                    "peak_a_MS2": "100,10",
                    "peak_b_MS2": "100,10",
                    "ms2_diagnostic_status": "supported_same_compound",
                    "ms2_cosine": "0.95",
                },
            )
            self._write(
                candidate,
                fields,
                {
                    **common,
                    "peak_a_MS2": "100,10;120,20",
                    "peak_b_MS2": "",
                    "ms2_diagnostic_status": "not_evaluable",
                    "ms2_cosine": "",
                },
            )
            output = root / "queue.csv"
            summary = compare_ms2_runs(
                baseline, candidate, output, status_changes_only=True
            )
            self.assertEqual(summary["status_changed_rows"], 1)
            self.assertEqual(summary["queued_rows"], 1)
            rows = list(csv.DictReader(output.open(encoding="utf-8-sig")))
            self.assertEqual(rows[0]["comparison_baseline_status"], "supported_same_compound")
            self.assertEqual(rows[0]["comparison_candidate_status"], "not_evaluable")
            self.assertEqual(rows[0]["comparison_peak_a_spectrum_changed"], "true")
            self.assertEqual(rows[0]["comparison_review_tier"], "manual_parameter_flip")
            payload = json.loads(
                output.with_suffix(".csv.summary.json").read_text(encoding="utf-8")
            )
            self.assertEqual(payload["matched_rows"], 1)

    @staticmethod
    def _write(path, fields, row):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerow(row)


if __name__ == "__main__":
    unittest.main()
