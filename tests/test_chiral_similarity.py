import base64
import csv
import math
import struct
import tempfile
import unittest
import xml.etree.ElementTree as ET
import zlib
from pathlib import Path
from unittest.mock import patch

from MSAI.python.chiral_similarity import (
    ChiralPairThresholds,
    SpectrumSimilarity,
    classify_chiral_pair,
    compare_fragment_spectra,
    parse_fragment_string,
)
from MSAI.python.calibrate_thresholds import calibrate_thresholds
from MSAI.python.get_chiral_frag import analyze_chiral_peak_pairs
from MSAI.python.easms_statistics import (
    benjamini_hochberg,
    enantiomer_log2_ratio_shift,
    enrichment_fold,
)
from MSAI.python.ms1_peak_picker import PeakPickingConfig, pick_chiral_peaks, savgol_smooth
from MSAI.python.ms2_core import (
    Spectrum,
    _decode_mzxml_peaks,
    _parse_mzxml_duration,
    build_ms2_index,
    matching_dia_windows,
    read_table,
)


class SimilarityTests(unittest.TestCase):
    def test_identical_spectra_have_unit_cosine(self):
        spectrum = "100,100;120,400;140,900;160,1600;180,2500;200,3600"
        result = compare_fragment_spectra(spectrum, spectrum)
        self.assertAlmostEqual(result.cosine, 1.0)
        self.assertEqual(result.matched_peaks, 6)
        self.assertAlmostEqual(result.peak_a_explained_intensity, 1.0)
        self.assertAlmostEqual(result.peak_b_explained_intensity, 1.0)
        self.assertAlmostEqual(result.entropy_similarity, 1.0)
        self.assertEqual(classify_chiral_pair(result)[0], "candidate_enantiomer_pair")

    def test_peak_matching_is_exclusive(self):
        result = compare_fragment_spectra(
            [(100.000, 100), (100.006, 90)],
            [(100.003, 100)],
            fragment_mz_tol=0.01,
            min_relative_intensity=0,
        )
        self.assertEqual(result.matched_peaks, 1)

    def test_disjoint_spectra_have_zero_entropy_similarity(self):
        result = compare_fragment_spectra(
            [(100, 100), (120, 50)],
            [(200, 100), (220, 50)],
            min_relative_intensity=0,
        )
        self.assertAlmostEqual(result.entropy_similarity, 0.0)

    def test_empty_and_low_match_spectra_are_not_candidates(self):
        empty = compare_fragment_spectra("", "100,100")
        self.assertIsNone(empty.cosine)
        self.assertEqual(classify_chiral_pair(empty)[0], "not_evaluable")
        low_match = SpectrumSimilarity(0.99, 2, 2, 2, 1.0, 1.0)
        self.assertEqual(classify_chiral_pair(low_match)[0], "insufficient_ms2_evidence")

    def test_fragment_parser_drops_bad_values(self):
        self.assertEqual(
            parse_fragment_string("100,10;bad;120,-1;140,nan;160,20"),
            [(100.0, 10.0), (160.0, 20.0)],
        )


class MzXmlFallbackTests(unittest.TestCase):
    def test_decodes_zlib_network_order_64_bit_pairs(self):
        payload = struct.pack(">4d", 100.25, 1234.5, 200.5, 6789.0)
        element = ET.fromstring(
            '<peaks compressionType="zlib" precision="64" byteOrder="network">'
            + base64.b64encode(zlib.compress(payload)).decode("ascii")
            + "</peaks>"
        )
        mz, intensity = _decode_mzxml_peaks(element)
        self.assertEqual(mz, [100.25, 200.5])
        self.assertEqual(intensity, [1234.5, 6789.0])

    def test_parses_retention_time_duration(self):
        self.assertEqual(_parse_mzxml_duration("PT1M2.5S"), 62.5)
        self.assertEqual(_parse_mzxml_duration("PT0.37511S"), 0.37511)


class PairWorkflowTests(unittest.TestCase):
    @staticmethod
    def _spectra():
        fragment_mz = [100, 120, 140, 160, 180, 200]
        fragment_intensity = [3000, 4000, 5000, 6000, 7000, 8000]
        spectra = []
        for center, mass_shift, scale in ((60, 0.0, 1.0), (120, 0.002, 1.5)):
            for rt, apex_scale in ((center - 5, 0.0), (center, 1.0), (center + 5, 0.0)):
                spectra.append(
                    Spectrum(
                        precursor_mz=300,
                        rt=rt,
                        mz=[295] + [mz + mass_shift for mz in fragment_mz],
                        intensity=[10000 * apex_scale]
                        + [value * scale * apex_scale for value in fragment_intensity],
                    )
                )
        return spectra

    def test_per_row_rt_workflow_writes_similarity_and_status(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            peaklist = root / "peaks.csv"
            output = root / "result.csv"
            with peaklist.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["Compound_ID", "MZ", "Peak1", "Peak2"])
                writer.writeheader()
                writer.writerow({"Compound_ID": "test", "MZ": "295", "Peak1": "1", "Peak2": "2"})

            with patch("MSAI.python.get_chiral_frag.load_ms2_spectra", return_value=self._spectra()):
                rows = analyze_chiral_peak_pairs(
                    peaklist,
                    root / "raw.mzXML",
                    output,
                    rt_half_window_sec=8,
                    thresholds=ChiralPairThresholds(),
                    write_metadata=False,
                )

            self.assertTrue(output.exists())
            self.assertEqual(rows[0]["enantiomer_pair_status"], "candidate_enantiomer_pair")
            self.assertEqual(rows[0]["ms2_matched_peaks"], "6")
            self.assertAlmostEqual(float(rows[0]["ms2_cosine"]), 1.0)
            self.assertEqual(rows[0]["peak_a_quality_flags"], "ok")
            self.assertEqual(rows[0]["peak_b_quality_flags"], "ok")
            self.assertEqual(rows[0]["compound_identity_status"], "ms2_supported_same_compound")
            self.assertEqual(rows[0]["chiral_doublet_status"], "candidate_chiral_doublet")
            self.assertTrue(rows[0]["enantioselective_enrichment_status"].startswith("not_assessed"))


class IndexAndPeakPickingTests(unittest.TestCase):
    def test_real_isolation_bounds_and_overlap_are_used(self):
        index = build_ms2_index(
            [
                Spectrum(100, 1, [50], [1], 90, 110),
                Spectrum(120, 1, [50], [1], 109, 131),
            ]
        )
        self.assertEqual(matching_dia_windows(109.5, index, 15), [100, 120])
        self.assertEqual(matching_dia_windows(130, index, 15), [120])

    def test_savgol_preserves_constant_and_picks_two_peaks(self):
        self.assertTrue(all(abs(value - 5) < 1e-9 for value in savgol_smooth([5.0] * 21, 7, 2)))
        rts = [float(index) for index in range(100)]
        intensity = [
            220_000 * math.exp(-((index - 25) / 4) ** 2)
            + 180_000 * math.exp(-((index - 70) / 5) ** 2)
            for index in range(100)
        ]
        result = pick_chiral_peaks(
            rts,
            intensity,
            PeakPickingConfig(min_distance_scans=20),
        )
        self.assertEqual(result.chromatographic_status, "double_peak")
        self.assertEqual([peak.scan_index for peak in result.peaks], [25, 70])
        self.assertGreater(result.resolution, 1)
        self.assertTrue(all((peak.area_fwhm or 0) > 0 for peak in result.peaks))

    def test_blank_csv_rows_are_removed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "rows.csv"
            path.write_text("id,MZ\nA,100\n,\nB,200\n", encoding="utf-8")
            self.assertEqual(len(read_table(path)), 2)


class CalibrationTests(unittest.TestCase):
    def test_calibration_finds_a_perfect_separating_rule(self):
        rows = []
        for label, cosine, entropy, matched, explained in (
            ("positive", 0.95, 0.95, 10, 0.9),
            ("positive", 0.90, 0.90, 8, 0.8),
            ("negative", 0.40, 0.35, 2, 0.2),
            ("negative", 0.55, 0.50, 3, 0.3),
        ):
            rows.append(
                {
                    "truth": label,
                    "ms2_cosine": str(cosine),
                    "ms2_entropy_similarity": str(entropy),
                    "ms2_matched_peaks": str(matched),
                    "peak_a_explained_intensity": str(explained),
                    "peak_b_explained_intensity": str(explained),
                }
            )
        result = calibrate_thresholds(rows, "truth")
        self.assertEqual(result.balanced_accuracy, 1.0)
        self.assertEqual(result.precision, 1.0)


class EasmsStatisticsTests(unittest.TestCase):
    def test_enrichment_and_ratio_shift(self):
        self.assertEqual(enrichment_fold(10, [1, 2, 3]), 5.0)
        self.assertAlmostEqual(enantiomer_log2_ratio_shift(1, 1, 3, 1, 1), 1.0)

    def test_bh_q_values_are_monotone_and_keep_missing(self):
        q_values = benjamini_hochberg([0.01, 0.04, 0.03, None])
        self.assertEqual(q_values, [0.03, 0.04, 0.04, None])


if __name__ == "__main__":
    unittest.main()
