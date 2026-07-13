"""Backward-compatible facade for the modular :mod:`ms2` package.

New code may import focused modules such as ``ms2.raw_io`` or
``ms2.extraction``.  Existing project scripts and third-party callers can keep
using ``ms2_core`` without changing object identity or function signatures.
"""

from __future__ import annotations

try:
    from . import ms2 as _impl
except ImportError:  # direct ``python MSAI/python/script.py`` execution
    import ms2 as _impl  # type: ignore[import-not-found]

Spectrum = _impl.Spectrum
DiaData = _impl.DiaData
Ms2Index = _impl.Ms2Index
Eic = _impl.Eic

build_ms2_index = _impl.build_ms2_index
load_ms2_spectra = _impl.load_ms2_spectra
preclist = _impl.preclist
normalize_rt_window = _impl.normalize_rt_window
mz_bounds = _impl.mz_bounds
summarize_xic_peak = _impl.summarize_xic_peak
format_fragment_string = _impl.format_fragment_string
extract_window_result = _impl.extract_window_result
workspace_paths = _impl.workspace_paths
peak_files = _impl.peak_files
raw_files = _impl.raw_files
matching_raw_file = _impl.matching_raw_file
read_table = _impl.read_table
write_table = _impl.write_table
column = _impl.column
number = _impl.number
matching_dia_windows = _impl.matching_dia_windows
first_dia_window = _impl.first_dia_window
csv_scalar = _impl.csv_scalar
empty_window_result = _impl.empty_window_result
merge_window_result = _impl.merge_window_result
pearson_correlation = _impl.pearson_correlation
raw_acquisition_context = _impl.raw_acquisition_context
sample_sd = _impl.sample_sd

# Private aliases are retained because MS1 and legacy tests historically used
# the XML codecs and low-level helpers from this module.
_ensure_sorted_peaks = _impl._ensure_sorted_peaks
_load_mzxml_stdlib = _impl._load_mzxml_stdlib
_load_mzml_stdlib = _impl._load_mzml_stdlib
_mzml_cv_value = _impl._mzml_cv_value
_mzml_scan_time = _impl._mzml_scan_time
_decode_mzml_arrays = _impl._decode_mzml_arrays
_local_name = _impl._local_name
_parse_mzxml_duration = _impl._parse_mzxml_duration
_decode_mzxml_peaks = _impl._decode_mzxml_peaks
_ms2copy = _impl._ms2copy
_raw_eic = _impl._raw_eic
_merge_candidate_peaks = _impl._merge_candidate_peaks
_clamped_bounds = _impl._clamped_bounds
_string_row = _impl._string_row
_row_has_data = _impl._row_has_data
_empty_xic_summary = _impl._empty_xic_summary
_split_flags = _impl._split_flags
_sd = _impl._sd
_cor = _impl._cor

__all__ = [
    "DiaData",
    "Eic",
    "Ms2Index",
    "Spectrum",
    "build_ms2_index",
    "column",
    "csv_scalar",
    "empty_window_result",
    "extract_window_result",
    "first_dia_window",
    "format_fragment_string",
    "load_ms2_spectra",
    "matching_dia_windows",
    "matching_raw_file",
    "merge_window_result",
    "mz_bounds",
    "normalize_rt_window",
    "number",
    "peak_files",
    "pearson_correlation",
    "preclist",
    "raw_acquisition_context",
    "raw_files",
    "read_table",
    "sample_sd",
    "summarize_xic_peak",
    "workspace_paths",
    "write_table",
]
