"""Composable MS2 primitives.

The package is split by responsibility so raw decoding, indexing, EIC work,
fragment extraction, and table I/O can be tested independently.  Public names
remain available through :mod:`MSAI.python.ms2_core` for legacy callers.
"""

from .chromatograms import (
    _empty_xic_summary,
    _ms2copy,
    _raw_eic,
    _split_flags,
    empty_window_result,
    format_fragment_string,
    merge_window_result,
    summarize_xic_peak,
)
from .extraction import _clamped_bounds, _merge_candidate_peaks, extract_window_result
from .indexing import (
    _ensure_sorted_peaks,
    build_ms2_index,
    first_dia_window,
    matching_dia_windows,
    preclist,
)
from .models import DiaData, Eic, Ms2Index, Spectrum
from .raw_io import (
    _decode_mzml_arrays,
    _load_mzml_stdlib,
    _load_mzxml_stdlib,
    _mzml_cv_value,
    _mzml_scan_time,
    load_ms2_spectra,
)
from .tables import _row_has_data, _string_row, column, read_table, write_table
from .utils import (
    _cor,
    _sd,
    csv_scalar,
    mz_bounds,
    normalize_rt_window,
    number,
    pearson_correlation,
    sample_sd,
)
from .workspace import matching_raw_file, peak_files, raw_files, workspace_paths
from .xml_codec import _decode_mzxml_peaks, _local_name, _parse_mzxml_duration

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
    "raw_files",
    "read_table",
    "sample_sd",
    "summarize_xic_peak",
    "workspace_paths",
    "write_table",
]
