"""手性双峰 fragment extraction。

这个文件只放 GetChiralFrag.R 对应的新逻辑：
- 同一个 precursor m/z 在 peak_a / peak_b 两个 RT window 内分别提取；
- 每个窗口都有独立的 XIC summary 和 MS2 fragments；
- 输出 results/{sample}_chiral_MS2.csv。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

try:
    from .chiral_similarity import (
        ChiralPairThresholds,
        SpectrumSimilarity,
        classify_chiral_pair,
        compare_fragment_spectra,
    )
    from .ms2 import (
        Eic,
        Ms2Index,
        Spectrum,
        build_ms2_index,
        column,
        csv_scalar,
        empty_window_result,
        extract_window_result,
        first_dia_window,
        format_fragment_string,
        load_ms2_spectra,
        matching_dia_windows,
        matching_raw_file,
        mz_bounds,
        normalize_rt_window,
        number,
        peak_files,
        preclist,
        raw_files,
        read_table,
        summarize_xic_peak,
        workspace_paths,
        write_table,
    )
except ImportError:
    from chiral_similarity import (  # type: ignore
        ChiralPairThresholds,
        SpectrumSimilarity,
        classify_chiral_pair,
        compare_fragment_spectra,
    )
    from ms2 import (  # type: ignore
        Eic,
        Ms2Index,
        Spectrum,
        build_ms2_index,
        column,
        csv_scalar,
        empty_window_result,
        extract_window_result,
        first_dia_window,
        format_fragment_string,
        load_ms2_spectra,
        matching_dia_windows,
        matching_raw_file,
        mz_bounds,
        normalize_rt_window,
        number,
        peak_files,
        preclist,
        raw_files,
        read_table,
        summarize_xic_peak,
        workspace_paths,
        write_table,
    )


def GetChiralFrag(
    mz_tol: float,
    DIAisowin: float,
    peak_a_rt_window,
    peak_b_rt_window,
    rt_unit: str = "min",
    mz_tol_unit: str = "legacy_fraction",
    base_path=None,
):
    """兼容 R 函数名的入口；实际实现见 get_chiral_frag()."""

    return get_chiral_frag(
        mz_tol,
        DIAisowin,
        peak_a_rt_window,
        peak_b_rt_window,
        rt_unit,
        mz_tol_unit,
        base_path,
    )


def get_chiral_frag(
    mz_tol: float,
    dia_iso_win: float,
    peak_a_rt_window,
    peak_b_rt_window,
    rt_unit: str = "min",
    mz_tol_unit: str = "legacy_fraction",
    base_path=None,
):
    """手性双峰 fragment extraction。"""

    base = Path(base_path or ".").resolve()
    paths = workspace_paths(base)

    # pyopenms 里的 RT 是秒；用户常用分钟，所以入口处统一转换。
    rt_windows = {
        "peak_a": normalize_rt_window(peak_a_rt_window, rt_unit),
        "peak_b": normalize_rt_window(peak_b_rt_window, rt_unit),
    }

    # 用一个 dummy mz 提前校验 mz tolerance 参数是否合法。
    mz_bounds(1.0, mz_tol, mz_tol_unit)
    if rt_windows["peak_a"][1] > rt_windows["peak_b"][0] and rt_windows["peak_b"][1] > rt_windows["peak_a"][0]:
        warnings.warn("Chiral RT windows overlap; no deconvolution is performed.")

    peakfiles = peak_files(paths["peak"])
    if not peakfiles:
        warnings.warn(f"No peaklist files found in {paths['peak']}; returning None.")
        return None

    msfiles = raw_files(paths["data"])
    last_rows = None
    for k, peakfile in enumerate(peakfiles, start=1):
        print(f"getchiralfragment... {k}")
        rows = read_table(peakfile)
        if not rows:
            continue
        mz_col = column(rows, "mz")
        if mz_col is None:
            warnings.warn(f"{peakfile.name} has no mz/MZ column; skipping.")
            continue

        # 先补齐输出列，保证即使某些 feature 无结果，CSV 结构也稳定。
        _initialize_chiral_columns(rows)

        rawfile = matching_raw_file(peakfile, msfiles)
        if rawfile is None:
            warnings.warn(f"No matching raw MS data found for peaklist {peakfile.name}; skipping.")
            continue

        spectra = load_ms2_spectra(rawfile)
        spectra_index = build_ms2_index(spectra)
        precursor = spectra_index.precursors
        for row in rows:
            mz = number(row.get(mz_col))
            if mz is None:
                continue
            diawin = first_dia_window(mz, spectra_index, dia_iso_win)
            if diawin is None:
                # 没有 DIA window 时，两个手性窗口都写同一个质量标记。
                for prefix in rt_windows:
                    _assign_chiral_result(row, prefix, empty_window_result("no_matching_DIA_window"))
                continue
            for prefix, rt_window in rt_windows.items():
                # 核心算法只处理“一个 feature + 一个 RT window”。
                # 手性双峰就是调用两次同一个 helper。
                result = extract_fragments_for_rt_window(
                    spectra=spectra_index,
                    precurmz=mz,
                    mz_tol=mz_tol,
                    mz_tol_unit=mz_tol_unit,
                    diawin=diawin,
                    rt_window_sec=rt_window,
                )
                _assign_chiral_result(row, prefix, result)

        sample = peakfile.stem
        write_table(paths["results"] / f"{sample}_chiral_MS2.csv", rows)
        last_rows = rows
    return last_rows


def analyze_chiral_peak_pairs(
    peaklist_path,
    raw_path,
    output_path=None,
    *,
    mz_column: str = "MZ",
    peak_a_column: str = "Peak1",
    peak_b_column: str = "Peak2",
    rt_unit: str = "min",
    rt_half_window_sec: float = 8.0,
    mz_tol: float = 10.0,
    mz_tol_unit: str = "ppm",
    dia_iso_win: float = 15.0,
    fragment_mz_tol: float = 0.01,
    fragment_mz_tol_unit: str = "Da",
    min_relative_intensity: float = 0.01,
    min_fragment_intensity: float = 2000.0,
    min_fragment_relative_intensity: float = 0.0,
    min_fragment_correlation: float = 0.9,
    consensus_scans: int = 1,
    thresholds: ChiralPairThresholds | None = None,
    write_metadata: bool = True,
    preloaded_index: Ms2Index | None = None,
):
    """Analyze the per-row ``MZ + Peak1 + Peak2`` format used by this project.

    Unlike :func:`get_chiral_frag`, this entry point does not require one pair
    of global RT windows or filename-based raw-data matching.  Each row gets
    two windows centered on its own picked RT values.  Close peaks receive
    narrower windows automatically so the two extraction windows never
    overlap.
    """

    started = time.perf_counter()
    peaklist_path = Path(peaklist_path).resolve()
    raw_path = Path(raw_path).resolve()
    if output_path is None:
        project_dir = (
            peaklist_path.parent.parent
            if peaklist_path.parent.name.lower() == "peaklist"
            else peaklist_path.parent
        )
        output_path = project_dir / "results" / f"{peaklist_path.stem}_chiral_MS2_similarity.csv"
    output_path = Path(output_path).resolve()
    if not _is_finite_positive(rt_half_window_sec):
        raise ValueError("rt_half_window_sec must be finite and positive.")
    if not _is_finite_positive(dia_iso_win):
        raise ValueError("dia_iso_win must be finite and positive.")

    # Validate tolerance/unit arguments before the expensive raw-data load.
    normalize_rt_window((0, 1), rt_unit)
    mz_bounds(1.0, mz_tol, mz_tol_unit)
    compare_fragment_spectra(
        "",
        "",
        fragment_mz_tol=fragment_mz_tol,
        fragment_mz_tol_unit=fragment_mz_tol_unit,
        min_relative_intensity=min_relative_intensity,
    )
    extract_window_result(
        [],
        1.0,
        mz_tol,
        mz_tol_unit,
        0.0,
        (0.0, 1.0),
        min_fragment_intensity=min_fragment_intensity,
        min_fragment_relative_intensity=min_fragment_relative_intensity,
        min_fragment_correlation=min_fragment_correlation,
        consensus_scans=consensus_scans,
    )
    thresholds = thresholds or ChiralPairThresholds()
    rows = read_table(peaklist_path)
    if not rows:
        raise ValueError(f"Peaklist is empty: {peaklist_path}")
    mz_col = _resolve_column(rows, mz_column)
    peak_a_col = _resolve_column(rows, peak_a_column)
    peak_b_col = _resolve_column(rows, peak_b_column)
    missing = [
        requested
        for requested, resolved in (
            (mz_column, mz_col),
            (peak_a_column, peak_a_col),
            (peak_b_column, peak_b_col),
        )
        if resolved is None
    ]
    if missing:
        raise ValueError(f"Peaklist is missing required column(s): {', '.join(missing)}")

    _initialize_chiral_columns(rows)
    _initialize_similarity_columns(rows)
    if preloaded_index is None:
        spectra = load_ms2_spectra(raw_path)
        if not spectra:
            raise ValueError(f"No MS2 spectra found in raw file: {raw_path}")
        spectra_index = build_ms2_index(spectra)
    else:
        spectra_index = preloaded_index
    precursor = spectra_index.precursors

    for row in rows:
        mz = number(row.get(mz_col))
        peak_a_rt = number(row.get(peak_a_col))
        peak_b_rt = number(row.get(peak_b_col))
        if peak_a_rt is None or peak_b_rt is None:
            _set_status_layers(
                row,
                "not_evaluated",
                "two_peak_retention_times_required",
                chromatographic_status=(
                    "input_single_peak" if peak_a_rt is not None or peak_b_rt is not None else "input_no_peak"
                ),
            )
            continue
        if mz is None:
            _set_status_layers(row, "not_evaluable", "invalid_precursor_mz")
            continue

        peak_a_rt_sec = _rt_to_seconds(peak_a_rt, rt_unit)
        peak_b_rt_sec = _rt_to_seconds(peak_b_rt, rt_unit)
        separation = abs(peak_b_rt_sec - peak_a_rt_sec)
        row["rt_separation_sec"] = csv_scalar(separation)
        if separation == 0:
            _set_status_layers(row, "not_evaluable", "peak_retention_times_are_identical")
            continue

        # Leave a 20% gap between windows.  This avoids cross-contaminating the
        # spectra of close chromatographic peaks while retaining the requested
        # half-width for well-separated peaks.
        half_width = min(float(rt_half_window_sec), separation * 0.4)
        rt_windows = {
            "peak_a": (peak_a_rt_sec - half_width, peak_a_rt_sec + half_width),
            "peak_b": (peak_b_rt_sec - half_width, peak_b_rt_sec + half_width),
        }
        for prefix, rt_window in rt_windows.items():
            row[f"{prefix}_rt_window_start"] = csv_scalar(rt_window[0])
            row[f"{prefix}_rt_window_end"] = csv_scalar(rt_window[1])

        dia_matches = matching_dia_windows(mz, spectra_index, dia_iso_win)
        diawin = dia_matches[0] if dia_matches else None
        if diawin is None:
            for prefix in rt_windows:
                _assign_chiral_result(row, prefix, empty_window_result("no_matching_DIA_window"))
            _set_status_layers(row, "not_evaluable", "precursor_outside_acquired_DIA_windows")
            continue

        dia_data = spectra_index.get(diawin)
        row["dia_window_center"] = csv_scalar(diawin)
        row["dia_window_lower"] = csv_scalar(dia_data.isolation_lower_mz)
        row["dia_window_upper"] = csv_scalar(dia_data.isolation_upper_mz)
        row["dia_window_match_count"] = csv_scalar(len(dia_matches))

        extracted: dict[str, dict] = {}
        for prefix, rt_window in rt_windows.items():
            result = extract_fragments_for_rt_window(
                spectra=spectra_index,
                precurmz=mz,
                mz_tol=mz_tol,
                mz_tol_unit=mz_tol_unit,
                diawin=diawin,
                rt_window_sec=rt_window,
                min_fragment_intensity=min_fragment_intensity,
                min_fragment_relative_intensity=min_fragment_relative_intensity,
                min_fragment_correlation=min_fragment_correlation,
                consensus_scans=consensus_scans,
            )
            extracted[prefix] = result
            _assign_chiral_result(row, prefix, result)

        similarity = compare_fragment_spectra(
            extracted["peak_a"]["fragment_peaks"],
            extracted["peak_b"]["fragment_peaks"],
            fragment_mz_tol=fragment_mz_tol,
            fragment_mz_tol_unit=fragment_mz_tol_unit,
            min_relative_intensity=min_relative_intensity,
        )
        status, reason = classify_chiral_pair(similarity, thresholds)
        _assign_similarity_result(row, similarity, status, reason)

    write_table(output_path, rows)
    if write_metadata:
        _write_run_metadata(
            output_path,
            peaklist_path,
            raw_path,
            rows,
            spectra_index,
            elapsed_seconds=time.perf_counter() - started,
            parameters={
                "mz_column": mz_column,
                "peak_a_column": peak_a_column,
                "peak_b_column": peak_b_column,
                "rt_unit": rt_unit,
                "rt_half_window_sec": rt_half_window_sec,
                "mz_tol": mz_tol,
                "mz_tol_unit": mz_tol_unit,
                "dia_iso_win_fallback": dia_iso_win,
                "fragment_mz_tol": fragment_mz_tol,
                "fragment_mz_tol_unit": fragment_mz_tol_unit,
                "min_relative_intensity": min_relative_intensity,
                "min_fragment_intensity": min_fragment_intensity,
                "min_fragment_relative_intensity": min_fragment_relative_intensity,
                "min_fragment_correlation": min_fragment_correlation,
                "consensus_scans": consensus_scans,
                "min_cosine": thresholds.min_cosine,
                "min_matched_peaks": thresholds.min_matched_peaks,
                "min_explained_intensity": thresholds.min_explained_intensity,
                "min_entropy_similarity": thresholds.min_entropy_similarity,
            },
        )
    return rows


def extract_fragments_for_rt_window(
    spectra: list[Spectrum] | Ms2Index,
    precurmz: float,
    mz_tol: float,
    mz_tol_unit: str,
    diawin: float,
    rt_window_sec: tuple[float, float],
    min_fragment_intensity: float = 2000.0,
    min_fragment_relative_intensity: float = 0.0,
    min_fragment_correlation: float = 0.9,
    consensus_scans: int = 1,
):
    """公开一点的单窗口 helper，便于测试或更细粒度调用。"""

    return extract_window_result(
        spectra,
        precurmz,
        mz_tol,
        mz_tol_unit,
        diawin,
        rt_window_sec,
        min_fragment_intensity=min_fragment_intensity,
        min_fragment_relative_intensity=min_fragment_relative_intensity,
        min_fragment_correlation=min_fragment_correlation,
        consensus_scans=consensus_scans,
    )


def _initialize_chiral_columns(rows: list[dict]):
    """给每一行预先加 peak_a_* / peak_b_* 输出列。"""

    for row in rows:
        for prefix in ("peak_a", "peak_b"):
            row[f"{prefix}_MS2"] = ""
            row[f"{prefix}_apex_rt"] = ""
            row[f"{prefix}_apex_intensity"] = ""
            row[f"{prefix}_area"] = ""
            row[f"{prefix}_ms2_count"] = ""
            row[f"{prefix}_consensus_scan_count"] = ""
            row[f"{prefix}_candidate_fragment_count"] = ""
            row[f"{prefix}_fragment_count"] = ""
            row[f"{prefix}_quality_flags"] = ""


def _initialize_similarity_columns(rows: list[dict]):
    """Add stable scalar columns for windows, similarity, and classification."""

    columns = (
        "peak_a_rt_window_start",
        "peak_a_rt_window_end",
        "peak_b_rt_window_start",
        "peak_b_rt_window_end",
        "rt_separation_sec",
        "dia_window_center",
        "dia_window_lower",
        "dia_window_upper",
        "dia_window_match_count",
        "ms2_cosine",
        "ms2_entropy_similarity",
        "ms2_matched_peaks",
        "peak_a_fragment_count",
        "peak_b_fragment_count",
        "peak_a_explained_intensity",
        "peak_b_explained_intensity",
        "enantiomer_pair_status",
        "enantiomer_pair_reason",
        "chromatographic_status",
        "compound_identity_status",
        "chiral_doublet_status",
        "enantioselective_enrichment_status",
    )
    for row in rows:
        for name in columns:
            row[name] = ""


def _assign_chiral_result(row: dict, prefix: str, result: dict):
    """把一个 window result dict 写回 CSV row。"""

    row[f"{prefix}_MS2"] = result["MS2"]
    row[f"{prefix}_apex_rt"] = csv_scalar(result["apex_rt"])
    row[f"{prefix}_apex_intensity"] = csv_scalar(result["apex_intensity"])
    row[f"{prefix}_area"] = csv_scalar(result["area"])
    row[f"{prefix}_ms2_count"] = csv_scalar(result["ms2_count"])
    row[f"{prefix}_consensus_scan_count"] = csv_scalar(result.get("consensus_scan_count"))
    row[f"{prefix}_candidate_fragment_count"] = csv_scalar(result.get("candidate_fragment_count"))
    row[f"{prefix}_fragment_count"] = csv_scalar(result.get("fragment_count"))
    row[f"{prefix}_quality_flags"] = result["quality_flags"]


def _assign_similarity_result(
    row: dict,
    similarity: SpectrumSimilarity,
    status: str,
    reason: str,
):
    row["ms2_cosine"] = csv_scalar(similarity.cosine)
    row["ms2_entropy_similarity"] = csv_scalar(similarity.entropy_similarity)
    row["ms2_matched_peaks"] = csv_scalar(similarity.matched_peaks)
    row["peak_a_fragment_count"] = csv_scalar(similarity.peak_a_count)
    row["peak_b_fragment_count"] = csv_scalar(similarity.peak_b_count)
    row["peak_a_explained_intensity"] = csv_scalar(similarity.peak_a_explained_intensity)
    row["peak_b_explained_intensity"] = csv_scalar(similarity.peak_b_explained_intensity)
    _set_status_layers(row, status, reason)


def _set_status_layers(
    row: dict,
    status: str,
    reason: str,
    *,
    chromatographic_status: str = "input_double_peak",
):
    """Keep chromatographic, identity, and enrichment conclusions separate."""

    row["enantiomer_pair_status"] = status
    row["enantiomer_pair_reason"] = reason
    row["chromatographic_status"] = chromatographic_status
    if status == "candidate_enantiomer_pair":
        row["compound_identity_status"] = "ms2_supported_same_compound"
        row["chiral_doublet_status"] = "candidate_chiral_doublet"
    elif status == "ms2_not_similar":
        row["compound_identity_status"] = "ms2_conflict"
        row["chiral_doublet_status"] = "not_supported"
    elif status == "insufficient_ms2_evidence":
        row["compound_identity_status"] = "insufficient_ms2_evidence"
        row["chiral_doublet_status"] = "not_evaluable"
    else:
        row["compound_identity_status"] = "not_evaluable"
        row["chiral_doublet_status"] = "not_evaluable"
    row["enantioselective_enrichment_status"] = "not_assessed_requires_input_eluate_replicates"


def _resolve_column(rows: list[dict], requested: str) -> str | None:
    """Prefer an exact header (important when both ``mz`` and ``MZ`` exist)."""

    if requested in rows[0]:
        return requested
    return column(rows, requested)


def _rt_to_seconds(value: float, unit: str) -> float:
    if unit == "min":
        return value * 60
    if unit == "sec":
        return value
    raise ValueError("rt_unit must be either 'min' or 'sec'.")


def _is_finite_positive(value) -> bool:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return False
    return value > 0 and value < float("inf")


def _write_run_metadata(
    output_path: Path,
    peaklist_path: Path,
    raw_path: Path,
    rows: list[dict],
    spectra_index: Ms2Index,
    *,
    elapsed_seconds: float,
    parameters: dict,
):
    counts: dict[str, int] = {}
    for row in rows:
        status = row.get("enantiomer_pair_status", "")
        counts[status] = counts.get(status, 0) + 1
    metadata = {
        "schema_version": "msai-chiral-ms2-v2",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "elapsed_seconds": elapsed_seconds,
        "rows_processed": len(rows),
        "status_counts": counts,
        "parameters": parameters,
        "inputs": {
            "peaklist": _file_fingerprint(peaklist_path),
            "raw": _file_fingerprint(raw_path),
        },
        "dia_windows": [
            {
                "center": center,
                "lower": spectra_index.windows[center].isolation_lower_mz,
                "upper": spectra_index.windows[center].isolation_upper_mz,
                "scan_count": len(spectra_index.windows[center].spectra),
            }
            for center in spectra_index.precursors
        ],
    }
    metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")


def _file_fingerprint(path: Path) -> dict:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    stat = path.stat()
    return {
        "path": str(path),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": digest.hexdigest(),
    }


def _demo():
    """最小自检：覆盖单位转换、mz tolerance、XIC 摘要和 fragment 筛选。"""

    assert normalize_rt_window((1, 2), "min") == (60, 120)
    assert normalize_rt_window((10, 20), "sec") == (10, 20)
    assert mz_bounds(100, 10, "ppm") == (99.999, 100.001)
    assert mz_bounds(100, 0.001, "Da") == (99.999, 100.001)
    assert mz_bounds(100, 0.00001, "legacy_fraction") == (99.999, 100.001)

    eic = Eic(rt=[0, 1, 2], scan=[0, 1, 2], intensity=[0, 10, 0])
    summary = summarize_xic_peak(eic)
    assert summary["apex_rt"] == 1
    assert summary["apex_intensity"] == 10
    assert summary["area"] == 10
    assert summary["quality_flags"] == "ok"
    assert format_fragment_string([], []) == ""
    assert format_fragment_string([50, 75], [1000, 2000]) == "50,1000;75,2000"

    spectra = [
        Spectrum(precursor_mz=300, rt=0, mz=[100, 290, 295], intensity=[0, 0, 0]),
        Spectrum(precursor_mz=300, rt=1, mz=[100, 290, 295], intensity=[5000, 10, 100]),
        Spectrum(precursor_mz=300, rt=2, mz=[100, 290, 295], intensity=[0, 0, 0]),
    ]
    result = extract_fragments_for_rt_window(spectra, 295, 0.01, "Da", 300, (0, 2))
    assert result["MS2"] == "100,5000"


def _main(argv=None):
    parser = argparse.ArgumentParser(description="Extract and compare MS2 spectra for per-row chiral LC peak pairs.")
    parser.add_argument("--peaklist", required=True, help="CSV/XLSX containing MZ, Peak1, and Peak2 columns.")
    parser.add_argument("--raw", required=True, help="Matching mzXML/mzML raw MS file.")
    parser.add_argument("--output", help="Output CSV; defaults to the project results directory.")
    parser.add_argument("--mz-column", default="MZ")
    parser.add_argument("--peak-a-column", default="Peak1")
    parser.add_argument("--peak-b-column", default="Peak2")
    parser.add_argument("--rt-unit", choices=("min", "sec"), default="min")
    parser.add_argument("--rt-half-window-sec", type=float, default=8.0)
    parser.add_argument("--mz-tol", type=float, default=10.0)
    parser.add_argument("--mz-tol-unit", choices=("ppm", "Da", "legacy_fraction"), default="ppm")
    parser.add_argument("--dia-window", type=float, default=15.0)
    parser.add_argument("--fragment-mz-tol", type=float, default=0.01)
    parser.add_argument("--fragment-mz-tol-unit", choices=("Da", "ppm"), default="Da")
    parser.add_argument("--min-relative-intensity", type=float, default=0.01)
    parser.add_argument("--min-fragment-intensity", type=float, default=2000.0)
    parser.add_argument("--min-fragment-relative-intensity", type=float, default=0.0)
    parser.add_argument("--min-fragment-correlation", type=float, default=0.9)
    parser.add_argument("--consensus-scans", type=int, default=1)
    parser.add_argument("--min-cosine", type=float, default=0.7)
    parser.add_argument(
        "--min-entropy-similarity",
        type=float,
        help="Optional calibrated entropy-similarity guardrail; disabled by default.",
    )
    parser.add_argument("--min-matched-peaks", type=int, default=6)
    parser.add_argument("--min-explained-intensity", type=float, default=0.5)
    parser.add_argument("--no-metadata", action="store_true")
    args = parser.parse_args(argv)
    thresholds = ChiralPairThresholds(
        min_cosine=args.min_cosine,
        min_matched_peaks=args.min_matched_peaks,
        min_explained_intensity=args.min_explained_intensity,
        min_entropy_similarity=args.min_entropy_similarity,
    )
    rows = analyze_chiral_peak_pairs(
        args.peaklist,
        args.raw,
        args.output,
        mz_column=args.mz_column,
        peak_a_column=args.peak_a_column,
        peak_b_column=args.peak_b_column,
        rt_unit=args.rt_unit,
        rt_half_window_sec=args.rt_half_window_sec,
        mz_tol=args.mz_tol,
        mz_tol_unit=args.mz_tol_unit,
        dia_iso_win=args.dia_window,
        fragment_mz_tol=args.fragment_mz_tol,
        fragment_mz_tol_unit=args.fragment_mz_tol_unit,
        min_relative_intensity=args.min_relative_intensity,
        min_fragment_intensity=args.min_fragment_intensity,
        min_fragment_relative_intensity=args.min_fragment_relative_intensity,
        min_fragment_correlation=args.min_fragment_correlation,
        consensus_scans=args.consensus_scans,
        thresholds=thresholds,
        write_metadata=not args.no_metadata,
    )
    counts: dict[str, int] = {}
    for row in rows:
        status = row["enantiomer_pair_status"]
        counts[status] = counts.get(status, 0) + 1
    print("; ".join(f"{key}={counts[key]}" for key in sorted(counts)))


if __name__ == "__main__":
    if len(sys.argv) == 1:
        _demo()
        print("ok")
    else:
        _main()
