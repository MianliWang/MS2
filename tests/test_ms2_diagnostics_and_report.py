import csv
import tempfile
import unittest
from pathlib import Path

from MSAI.python.chiral_similarity import ChiralPairThresholds, SpectrumSimilarity
from MSAI.python.diagnostic_standards import (
    classify_ms1_reference,
    classify_ms2_diagnostic,
)
from MSAI.python.ms2_review_report import annotate_rows, build_report


class DiagnosticSeparationTests(unittest.TestCase):
    def test_ms1_reference_does_not_claim_ms2_identity(self):
        result = classify_ms1_reference("1.2", "1.8")
        self.assertEqual(result.status, "supplied_double_peak")
        self.assertEqual(result.issue_codes, ())

    def test_supported_and_boundary_ms2_status_are_explicit(self):
        similarity = SpectrumSimilarity(0.75, 6, 7, 7, 0.8, 0.9, 0.72)
        result = classify_ms2_diagnostic(
            similarity,
            "candidate_enantiomer_pair",
            "screen_pass",
            thresholds=ChiralPairThresholds(),
        )
        self.assertEqual(result.status, "supported_same_compound")
        self.assertIn("NEAR_DECISION_BOUNDARY", result.issue_codes)
        self.assertEqual(result.review_priority, "P2_boundary_or_interference")

    def test_low_match_precedes_high_cosine(self):
        similarity = SpectrumSimilarity(0.99, 5, 5, 5, 1.0, 1.0, 0.98)
        result = classify_ms2_diagnostic(
            similarity,
            "insufficient_ms2_evidence",
            "matched_peaks<6",
        )
        self.assertEqual(result.status, "insufficient_evidence")
        self.assertIn("TOO_FEW_MATCHED_FRAGMENTS", result.issue_codes)


class ReportTests(unittest.TestCase):
    def _source_row(self):
        return {
            "Compound_ID": "example",
            "MZ": "300.1",
            "Peak1": "1.0",
            "Peak2": "2.0",
            "peak_a_MS2": "100,1000;120,900;140,800;160,700;180,600;200,500",
            "peak_b_MS2": "100.002,900;120.002,800;140.002,700;160.002,600;180.002,500;200.002,400",
            "peak_a_quality_flags": "ok",
            "peak_b_quality_flags": "ok",
            "dia_window_match_count": "1",
            "ms2_cosine": "0.98",
            "ms2_entropy_similarity": "0.97",
            "ms2_matched_peaks": "6",
            "peak_a_fragment_count": "6",
            "peak_b_fragment_count": "6",
            "peak_a_explained_intensity": "1",
            "peak_b_explained_intensity": "1",
            "enantiomer_pair_status": "candidate_enantiomer_pair",
            "enantiomer_pair_reason": "screen_pass",
        }

    def test_annotation_has_independent_status_columns(self):
        row = annotate_rows([self._source_row()])[0]
        self.assertEqual(row["ms1_reference_status"], "supplied_double_peak")
        self.assertEqual(row["ms2_diagnostic_status"], "supported_same_compound")
        self.assertEqual(row["ms2_review_priority"], "P2_boundary_or_interference")

    def test_shared_rt_and_dia_targets_are_flagged_for_review(self):
        first = self._source_row()
        second = dict(first)
        first["Compound_ID"], second["Compound_ID"] = "first", "second"
        first["dia_window_center"] = second["dia_window_center"] = "305"
        annotated = annotate_rows([first, second])
        self.assertTrue(all("SHARED_RT_DIA_TARGETS" in row["ms2_issue_codes"] for row in annotated))

    def test_static_report_and_manual_queue_are_written(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.csv"
            row = self._source_row()
            with source.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=row.keys())
                writer.writeheader()
                writer.writerow(row)
            template = Path("MSAI/templates/ms2_review_report_shell.html")
            outputs = build_report(source, root / "report", template_path=template)
            report = outputs["report"].read_text(encoding="utf-8")
            self.assertIn('data-contract-section="technical-summary"', report)
            self.assertIn("Peak 1 and Peak 2 mirror spectrum", report)
            self.assertNotIn("DATA_ANALYTICS_HTML_REPORT_RUNTIME", report)
            with outputs["queue"].open(encoding="utf-8-sig", newline="") as handle:
                queue = list(csv.DictReader(handle))
            self.assertEqual(len(queue), 1)
            self.assertIn("manual_ms2_label", queue[0])


if __name__ == "__main__":
    unittest.main()
