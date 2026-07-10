import unittest

from MSAI.python import ms2_core
from MSAI.python.ms2 import build_ms2_index, extract_window_result
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


if __name__ == "__main__":
    unittest.main()
