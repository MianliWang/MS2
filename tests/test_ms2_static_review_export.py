import csv
import json
import os
import tempfile
import unittest
from pathlib import Path

from MSAI.python.chiral_similarity import align_fragment_spectra
from MSAI.python.ms2_review.acquisition import (
    load_method_profile,
    raw_acquisition_context,
    reconcile_acquisition,
)
from MSAI.python.ms2_review.classification import review_views, spectrum_availability
from MSAI.python.ms2_review.exporter import export_ms2_review
from MSAI.python.ms2_review.model import prepare_mirror_spectrum
from MSAI.python.ms2_review.png import render_ms2_review_png
from MSAI.python.ms2_review.svg import render_ms2_review_svg


def _supported_row():
    return {
        "Compound_ID": "example<&>",
        "SGC ID for Pool": "POOL-A01",
        "Pooled Well": "A1",
        "MZ": "400.1",
        "IG": "GREEN",
        "Peak1": "1.0",
        "Peak2": "2.0",
        "peak_a_MS2": "100,1000;120,900;140,800;160,700;180,600;200,500;250,5",
        "peak_b_MS2": "100.002,900;120.002,800;140.002,700;160.002,600;180.002,500;200.002,400",
        "peak_a_quality_flags": "ok",
        "peak_b_quality_flags": "ok",
        "peak_a_rt_window_start": "52",
        "peak_a_rt_window_end": "68",
        "peak_b_rt_window_start": "112",
        "peak_b_rt_window_end": "128",
        "peak_a_apex_rt": "60",
        "peak_b_apex_rt": "120",
        "peak_a_apex_intensity": "100000",
        "peak_b_apex_intensity": "90000",
        "peak_a_area": "500000",
        "peak_b_area": "450000",
        "peak_a_ms2_count": "4",
        "peak_b_ms2_count": "4",
        "peak_a_consensus_scan_count": "1",
        "peak_b_consensus_scan_count": "1",
        "peak_a_candidate_fragment_count": "100",
        "peak_b_candidate_fragment_count": "100",
        "peak_a_fragment_count": "6",
        "peak_b_fragment_count": "6",
        "dia_window_center": "403",
        "dia_window_lower": "395.5",
        "dia_window_upper": "410.5",
        "dia_window_match_count": "1",
        "rt_separation_sec": "60",
        "ms2_cosine": "0.98",
        "ms2_entropy_similarity": "0.97",
        "ms2_matched_peaks": "6",
        "peak_a_explained_intensity": "0.99",
        "peak_b_explained_intensity": "1",
        "enantiomer_pair_status": "candidate_enantiomer_pair",
        "enantiomer_pair_reason": "screen_pass",
    }


class SpectrumPreparationTests(unittest.TestCase):
    def test_alignment_matches_similarity_filtering(self):
        peak_a, peak_b, matches = align_fragment_spectra(
            "100,1000;150,5", "100.002,900;150.002,4", min_relative_intensity=0.01
        )
        self.assertEqual(peak_a, [(100.0, 1000.0)])
        self.assertEqual(peak_b, [(100.002, 900.0)])
        self.assertEqual(matches, [(0, 0)])

    def test_availability_distinguishes_one_sided_spectra(self):
        self.assertEqual(spectrum_availability({"peak_a_MS2": "100,1"}), "peak1_only")
        self.assertEqual(spectrum_availability({"peak_b_MS2": "100,1"}), "peak2_only")
        self.assertEqual(spectrum_availability({}), "no_spectra")

    def test_views_keep_ms2_status_and_source_independent(self):
        row = _supported_row()
        row.update(
            {
                "ms2_diagnostic_status": "conflicting_spectra",
                "ms2_review_priority": "P1_conflict",
                "ms2_issue_codes": "LOW_COSINE_SIMILARITY",
                "ms1_reference_status": "supplied_double_peak",
            }
        )
        views = review_views(
            row, source_machine_id="IG", source_machine_label="GREEN", source_pool_id="POOL-A01"
        )
        self.assertIn("diagnostic_status/conflicting_spectra", views)
        self.assertIn("source_machine/IG/GREEN", views)
        self.assertIn("cross_stage_followup/ms1_supplied_double_ms2_conflict", views)


class StaticRenderTests(unittest.TestCase):
    def _metadata(self):
        profile = load_method_profile(
            "MSAI/config/acquisition_profiles/adductmlib_chemrxiv_20260129_v1.json"
        )
        return {
            "target_uid": "target",
            "config_hash": "abcdef",
            "standard_id": "MSAI-MS2-DIAGNOSTIC-v1",
            "source_context": "IG label: GREEN (expected_double_peak)",
            "source_pool_id": "POOL-A01",
            "source_pooled_well": "A1",
            "fragment_mz_tol": 0.01,
            "fragment_mz_tol_unit": "Da",
            "min_relative_intensity": 0.01,
            "min_cosine": 0.7,
            "min_matched_peaks": 6,
            "min_explained_intensity": 0.5,
            "dia_windows": [{"lower": 395.5, "upper": 410.5}],
            "acquisition_context": {
                "scope": "entire_raw_file",
                "observed_values": {
                    "activationMethod": [{"value": "HCD", "count": 10}],
                    "collisionEnergy": [{"value": "35.0", "count": 10}],
                    "windowWideness": [{"value": "15.0", "count": 10}],
                },
            },
            "method_profile": profile,
            "acquisition_reconciliation": {"issue_codes": ["RAW_VS_REFERENCE_DIA_WIDTH_DIFFER"]},
        }

    def test_svg_escapes_title_and_contains_acquisition_strip(self):
        row = _supported_row()
        row.update(
            {
                "ms2_diagnostic_status": "supported_same_compound",
                "ms2_review_priority": "P4_supported_audit",
                "ms2_issue_codes": "",
            }
        )
        mirror = prepare_mirror_spectrum(row)
        svg = render_ms2_review_svg(row=row, mirror=mirror, metadata=self._metadata())
        self.assertIn("example&lt;&amp;&gt;", svg)
        self.assertIn("Acquired DIA coverage", svg)
        self.assertIn("matched: thick + circle/triangle", svg)
        self.assertIn("raw observed (entire_raw_file)", svg)
        self.assertIn("not run-confirmed", svg)

    def test_png_is_valid_and_empty_spectrum_is_honest(self):
        row = _supported_row()
        row.update(
            {
                "peak_a_MS2": "",
                "peak_b_MS2": "",
                "ms2_diagnostic_status": "not_evaluable",
                "ms2_review_priority": "P5_not_evaluable",
                "ms2_issue_codes": "NO_ACQUIRED_DIA_WINDOW",
            }
        )
        mirror = prepare_mirror_spectrum(row)
        image = render_ms2_review_png(row=row, mirror=mirror, metadata=self._metadata(), scale=0.5)
        self.assertEqual(image.size, (620, 515))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "review.png"
            image.save(path)
            self.assertEqual(path.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")


class StaticExportTests(unittest.TestCase):
    def test_export_writes_all_rows_hardlinks_and_preserves_annotations(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "result.csv"
            first = _supported_row()
            second = dict(first)
            second.update(
                {
                    "Compound_ID": "no-window",
                    "MZ": "300",
                    "peak_a_MS2": "",
                    "peak_b_MS2": "",
                    "peak_a_quality_flags": "no_matching_DIA_window",
                    "peak_b_quality_flags": "no_matching_DIA_window",
                    "ms2_cosine": "",
                    "ms2_entropy_similarity": "",
                    "ms2_matched_peaks": "",
                    "peak_a_explained_intensity": "",
                    "peak_b_explained_intensity": "",
                    "enantiomer_pair_status": "not_evaluable",
                    "enantiomer_pair_reason": "precursor_outside_acquired_DIA_windows",
                }
            )
            with source.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=first.keys())
                writer.writeheader()
                writer.writerows((first, second))
            sidecar = {
                "parameters": {
                    "fragment_mz_tol": 0.01,
                    "fragment_mz_tol_unit": "Da",
                    "min_relative_intensity": 0.01,
                    "min_cosine": 0.7,
                    "min_matched_peaks": 6,
                    "min_explained_intensity": 0.5,
                },
                "dia_windows": [{"lower": 395.5, "upper": 410.5, "center": 403, "scan_count": 10}],
                "inputs": {},
            }
            Path(str(source) + ".metadata.json").write_text(json.dumps(sidecar), encoding="utf-8")
            output = root / "review"
            summary = export_ms2_review(
                source,
                output,
                standard_path="MSAI/standards/ms2_diagnostic_standard_v1.json",
                method_profile_path="MSAI/config/acquisition_profiles/adductmlib_chemrxiv_20260129_v1.json",
                png_scale=0.25,
            )
            self.assertEqual(summary["images"], 2)
            self.assertEqual(len(list((output / "assets").rglob("*.svg"))), 2)
            self.assertEqual(len(list((output / "assets").rglob("*.png"))), 2)
            no_window_views = list(
                (
                    output
                    / "views"
                    / "acquisition_coverage"
                    / "no_dia_window"
                    / "below_acquired_range"
                ).glob("*.svg")
            )
            self.assertEqual(len(no_window_views), 1)
            manifest = list(
                csv.DictReader(
                    (output / "metadata" / "target_manifest.csv").open(encoding="utf-8-sig")
                )
            )
            asset = output / manifest[1]["svg_asset_path"]
            self.assertTrue(os.path.samefile(asset, no_window_views[0]))
            run_manifest = json.loads(
                (output / "metadata" / "run_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                run_manifest["method_profile"]["profile_id"],
                "ADDUCTMLIB-CHEMRXIV-20260129-v1",
            )
            self.assertIn(
                "REFERENCE_PROFILE_NOT_CONFIRMED_FOR_RUN",
                run_manifest["acquisition_reconciliation"]["issue_codes"],
            )
            annotation = output / "annotations" / "_template" / "review_labels.csv"
            annotation.write_text("preserve-me", encoding="utf-8")
            export_ms2_review(
                source,
                output,
                standard_path="MSAI/standards/ms2_diagnostic_standard_v1.json",
                image_formats=("svg",),
                png_scale=0.25,
            )
            self.assertEqual(annotation.read_text(encoding="utf-8"), "preserve-me")


class AcquisitionProvenanceTests(unittest.TestCase):
    def test_complete_raw_attribute_aggregation_and_reference_reconciliation(self):
        with tempfile.TemporaryDirectory() as directory:
            raw = Path(directory) / "run_2.mzXML"
            raw.write_text(
                '<mzXML><scan msLevel="1"/><scan msLevel="2">'
                '<precursorMz activationMethod="HCD" collisionEnergy="35.0" '
                'windowWideness="15.0">403</precursorMz></scan>'
                '<scan msLevel="2"><precursorMz activationMethod="HCD" '
                'collisionEnergy="35.0" windowWideness="15.0">417</precursorMz>'
                "</scan></mzXML>",
                encoding="utf-8",
            )
            observed = raw_acquisition_context(raw)
            profile = load_method_profile(
                "MSAI/config/acquisition_profiles/adductmlib_chemrxiv_20260129_v1.json"
            )
            reconciliation = reconcile_acquisition(
                observed,
                profile,
                analysis_parameters={
                    "mz_tol": 10,
                    "mz_tol_unit": "ppm",
                    "rt_half_window_sec": 8,
                    "min_fragment_correlation": 0.9,
                },
                dia_windows=[{"lower": 395.5, "upper": 410.5}],
            )
            self.assertEqual(observed["scan_counts"], {"ms1": 1, "ms2": 2})
            self.assertEqual(
                observed["observed_values"]["collisionEnergy"],
                [{"value": "35.0", "count": 2}],
            )
            self.assertIn(
                "RAW_VS_REFERENCE_DIA_WIDTH_DIFFER",
                reconciliation["issue_codes"],
            )
            self.assertIn(
                "RAW_VS_REFERENCE_COLLISION_ENERGY_DIFFER",
                reconciliation["issue_codes"],
            )
            self.assertIn(
                "REFERENCE_THREE_DIA_FILES_CURRENT_INPUT_ONE",
                reconciliation["issue_codes"],
            )
            self.assertIn("CHROMATOGRAPHIC_CORRELATION_0_9", reconciliation["matches"])


if __name__ == "__main__":
    unittest.main()
