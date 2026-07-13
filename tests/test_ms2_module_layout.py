import unittest
from importlib import import_module

import MSAI.python.chiral_similarity as legacy_similarity
import MSAI.python.diagnostic_standards as legacy_diagnostics
from MSAI.python import ms2_core
from MSAI.python.ms2 import (
    build_ms2_index,
    diagnostics,
    extract_window_result,
    pipeline,
    similarity,
)
from MSAI.python.ms2.extraction import extract_window_result as focused_extract
from MSAI.python.ms2.indexing import build_ms2_index as focused_build_index
from MSAI.python.ms2.models import Eic, Ms2Index, Spectrum


class Ms2ModuleLayoutTests(unittest.TestCase):
    def test_legacy_facade_reexports_identical_objects(self):
        self.assertIs(ms2_core.Spectrum, Spectrum)
        self.assertIs(ms2_core.Eic, Eic)
        self.assertIs(ms2_core.Ms2Index, Ms2Index)
        self.assertIs(ms2_core.build_ms2_index, focused_build_index)
        self.assertIs(ms2_core.extract_window_result, focused_extract)

    def test_package_root_exports_focused_implementations(self):
        self.assertIs(build_ms2_index, focused_build_index)
        self.assertIs(extract_window_result, focused_extract)

    def test_legacy_workflow_facades_reexport_canonical_modules(self):
        legacy_pipeline = import_module("MSAI.python.get_chiral_frag")
        self.assertIs(
            legacy_similarity.compare_fragment_spectra, similarity.compare_fragment_spectra
        )
        self.assertIs(
            legacy_diagnostics.classify_ms2_diagnostic, diagnostics.classify_ms2_diagnostic
        )
        self.assertIs(legacy_pipeline.analyze_chiral_peak_pairs, pipeline.analyze_chiral_peak_pairs)


if __name__ == "__main__":
    unittest.main()
