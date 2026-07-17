"""手性双峰 fragment extraction。

这是第 4--7 步的高层编排器；底层实现分别位于本包的 indexing、
extraction、similarity 和 diagnostics 模块：
- 同一个 precursor m/z 在 peak_a / peak_b 两个 RT window 内分别提取；
- 每个窗口都有独立的 XIC summary 和 MS2 fragments；
- 输出用户指定的分析表和同名 metadata JSON。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
import warnings
from datetime import UTC, datetime
from pathlib import Path

from . import (
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
    raw_acquisition_context,
    raw_files,
    read_table,
    summarize_xic_peak,
    workspace_paths,
    write_table,
)
from .diagnostics import classify_ms1_reference, classify_ms2_diagnostic
from .similarity import (
    ChiralPairThresholds,
    SpectrumSimilarity,
    classify_chiral_pair,
    compare_fragment_spectra,
)


DEFAULT_RT_HALF_WINDOW_SEC = 10.0


def GetChiralFrag(
    mz_tol: float,
    DIAisowin: float,
    peak_a_rt_window,
    peak_b_rt_window,
    rt_unit: str = "min",
    mz_tol_unit: str = "legacy_fraction",
    base_path=None,
):
    """旧式兼容入口；参数原样交给 :func:`get_chiral_frag`。

    新代码不应依赖这个大写名称，而应调用逐行RT的
    :func:`analyze_chiral_peak_pairs`。
    """

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
    """按旧目录布局和两个全局RT窗口批量提取双侧fragment。

    它会扫描``peaklist/``与``data/``并按文件名匹配raw文件，所有目标共享
    ``peak_a_rt_window``和``peak_b_rt_window``。该行为适合兼容历史脚本，
    不适合当前每行具有独立``Peak1/Peak2``的正式E-ASMS流程。
    """

    base = Path(base_path or ".").resolve()
    paths = workspace_paths(base)

    # pyopenms 里的 RT 是秒；用户常用分钟，所以入口处统一转换。
    rt_windows = {
        "peak_a": normalize_rt_window(peak_a_rt_window, rt_unit),
        "peak_b": normalize_rt_window(peak_b_rt_window, rt_unit),
    }

    # 用一个 dummy mz 提前校验 mz tolerance 参数是否合法。
    mz_bounds(1.0, mz_tol, mz_tol_unit)
    if (
        rt_windows["peak_a"][1] > rt_windows["peak_b"][0]
        and rt_windows["peak_b"][1] > rt_windows["peak_a"][0]
    ):
        warnings.warn("Chiral RT windows overlap; no deconvolution is performed.", stacklevel=2)

    peakfiles = peak_files(paths["peak"])
    if not peakfiles:
        warnings.warn(f"No peaklist files found in {paths['peak']}; returning None.", stacklevel=2)
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
            warnings.warn(f"{peakfile.name} has no mz/MZ column; skipping.", stacklevel=2)
            continue

        # 先补齐输出列，保证即使某些 feature 无结果，CSV 结构也稳定。
        _initialize_chiral_columns(rows)

        rawfile = matching_raw_file(peakfile, msfiles)
        if rawfile is None:
            warnings.warn(
                f"No matching raw MS data found for peaklist {peakfile.name}; skipping.",
                stacklevel=2,
            )
            continue

        spectra = load_ms2_spectra(rawfile)
        spectra_index = build_ms2_index(spectra)
        for row in rows:
            mz = number(row.get(mz_col))
            if mz is None:
                continue
            diawin = first_dia_window(mz, spectra_index, dia_iso_win)
            if diawin is None:
                # 没有 DIA window 时，两个手性窗口都写同一个质量标记。
                for prefix in rt_windows:
                    _assign_chiral_result(
                        row, prefix, empty_window_result("no_matching_DIA_window")
                    )
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
    rt_half_window_sec: float = DEFAULT_RT_HALF_WINDOW_SEC,
    mz_tol: float = 10.0,
    mz_tol_unit: str = "ppm",
    precursor_eic_mz_tol: float | None = None,
    precursor_eic_mz_tol_unit: str | None = None,
    fragment_eic_mz_tol: float | None = None,
    fragment_eic_mz_tol_unit: str | None = None,
    dia_iso_win: float = 15.0,
    fragment_mz_tol: float = 0.01,
    fragment_mz_tol_unit: str = "Da",
    min_relative_intensity: float = 0.01,
    min_fragment_intensity: float = 2000.0,
    min_fragment_relative_intensity: float = 0.0,
    min_fragment_correlation: float = 0.9,
    fragment_correlation_mode: str = "full_window",
    correlation_min_relative_intensity: float = 0.05,
    min_correlation_scans: int = 5,
    max_fragment_apex_offset_scans: int = 1,
    min_consecutive_fragment_scans: int = 3,
    consensus_scans: int = 1,
    thresholds: ChiralPairThresholds | None = None,
    write_metadata: bool = True,
    preloaded_index: Ms2Index | None = None,
    method_profile_path: str | Path | None = None,
):
    """分析项目使用的逐行``MZ + Peak1 + Peak2``目标表。

    这是步骤4至7的正式编排边界：读取并索引一个raw文件；逐行检查输入RT和
    真实DIA覆盖；在Peak1/Peak2附近独立重建共洗脱多碎片谱；计算两张谱的
    cosine、entropy、匹配fragment数和双侧解释强度；最后写出分层诊断CSV及
    同名``.metadata.json``。

    ``Peak1/Peak2``默认以分钟输入，底层统一转成秒。两峰很近时，半窗口会
    自动缩小为间距的40%，从而在两个窗口之间保留20%空隙，减少相互污染。
    ``preloaded_index``允许参数扫描复用raw索引以显著减少运行时间。

    返回的是包含新增结果列的行列表。MS2支持相同化合物不等于确认对映体，
    因此色谱、化合物身份和富集结论分别保存在不同状态字段。
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

    resolved_precursor_eic_mz_tol = mz_tol if precursor_eic_mz_tol is None else precursor_eic_mz_tol
    resolved_precursor_eic_mz_tol_unit = (
        mz_tol_unit if precursor_eic_mz_tol_unit is None else precursor_eic_mz_tol_unit
    )
    resolved_fragment_eic_mz_tol = mz_tol if fragment_eic_mz_tol is None else fragment_eic_mz_tol
    resolved_fragment_eic_mz_tol_unit = (
        mz_tol_unit if fragment_eic_mz_tol_unit is None else fragment_eic_mz_tol_unit
    )

    # Validate tolerance/unit arguments before the expensive raw-data load.
    normalize_rt_window((0, 1), rt_unit)
    mz_bounds(1.0, mz_tol, mz_tol_unit)
    mz_bounds(1.0, resolved_precursor_eic_mz_tol, resolved_precursor_eic_mz_tol_unit)
    mz_bounds(1.0, resolved_fragment_eic_mz_tol, resolved_fragment_eic_mz_tol_unit)
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
        precursor_mz_tol=resolved_precursor_eic_mz_tol,
        precursor_mz_tol_unit=resolved_precursor_eic_mz_tol_unit,
        fragment_eic_mz_tol=resolved_fragment_eic_mz_tol,
        fragment_eic_mz_tol_unit=resolved_fragment_eic_mz_tol_unit,
        fragment_correlation_mode=fragment_correlation_mode,
        correlation_min_relative_intensity=correlation_min_relative_intensity,
        min_correlation_scans=min_correlation_scans,
        max_fragment_apex_offset_scans=max_fragment_apex_offset_scans,
        min_consecutive_fragment_scans=min_consecutive_fragment_scans,
    )
    thresholds = thresholds or ChiralPairThresholds()
    rows = read_table(peaklist_path)
    if not rows:
        raise ValueError(f"Peaklist is empty: {peaklist_path}")
    mz_col = _resolve_column(rows, mz_column)
    peak_a_col = _resolve_column(rows, peak_a_column)
    peak_b_col = _resolve_column(rows, peak_b_column)
    if mz_col is None or peak_a_col is None or peak_b_col is None:
        missing = [
            requested
            for requested, resolved in (
                (mz_column, mz_col),
                (peak_a_column, peak_a_col),
                (peak_b_column, peak_b_col),
            )
            if resolved is None
        ]
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
    for row in rows:
        mz = number(row.get(mz_col))
        peak_a_rt = number(row.get(peak_a_col))
        peak_b_rt = number(row.get(peak_b_col))
        ms1_reference = classify_ms1_reference(peak_a_rt, peak_b_rt)
        row["ms1_reference_status"] = ms1_reference.status
        row["ms1_reference_issue_codes"] = ";".join(ms1_reference.issue_codes)
        if peak_a_rt is None or peak_b_rt is None:
            _set_status_layers(
                row,
                "not_evaluated",
                "two_peak_retention_times_required",
                chromatographic_status=(
                    "input_single_peak"
                    if peak_a_rt is not None or peak_b_rt is not None
                    else "input_no_peak"
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

        # 对接近的色谱峰缩窄窗口并保留20%间隔：每侧半宽=0.4×峰间距，
        # 两个全宽合计占80%。这只是避免窗口直接重叠，并不执行色谱解卷积。
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
                precursor_eic_mz_tol=resolved_precursor_eic_mz_tol,
                precursor_eic_mz_tol_unit=resolved_precursor_eic_mz_tol_unit,
                fragment_eic_mz_tol=resolved_fragment_eic_mz_tol,
                fragment_eic_mz_tol_unit=resolved_fragment_eic_mz_tol_unit,
                fragment_correlation_mode=fragment_correlation_mode,
                correlation_min_relative_intensity=correlation_min_relative_intensity,
                min_correlation_scans=min_correlation_scans,
                max_fragment_apex_offset_scans=max_fragment_apex_offset_scans,
                min_consecutive_fragment_scans=min_consecutive_fragment_scans,
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
                "mz_tol_legacy_scope": "shared precursor and fragment EIC default",
                "precursor_eic_mz_tol": resolved_precursor_eic_mz_tol,
                "precursor_eic_mz_tol_unit": resolved_precursor_eic_mz_tol_unit,
                "fragment_eic_mz_tol": resolved_fragment_eic_mz_tol,
                "fragment_eic_mz_tol_unit": resolved_fragment_eic_mz_tol_unit,
                "dia_iso_win_fallback": dia_iso_win,
                "fragment_mz_tol": fragment_mz_tol,
                "fragment_mz_tol_unit": fragment_mz_tol_unit,
                "min_relative_intensity": min_relative_intensity,
                "min_fragment_intensity": min_fragment_intensity,
                "min_fragment_relative_intensity": min_fragment_relative_intensity,
                "min_fragment_correlation": min_fragment_correlation,
                "fragment_correlation_mode": fragment_correlation_mode,
                "correlation_min_relative_intensity": correlation_min_relative_intensity,
                "min_correlation_scans": min_correlation_scans,
                "max_fragment_apex_offset_scans": max_fragment_apex_offset_scans,
                "min_consecutive_fragment_scans": min_consecutive_fragment_scans,
                "consensus_scans": consensus_scans,
                "min_cosine": thresholds.min_cosine,
                "min_matched_peaks": thresholds.min_matched_peaks,
                "min_explained_intensity": thresholds.min_explained_intensity,
                "min_entropy_similarity": thresholds.min_entropy_similarity,
            },
            method_profile_path=method_profile_path,
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
    *,
    precursor_eic_mz_tol: float | None = None,
    precursor_eic_mz_tol_unit: str | None = None,
    fragment_eic_mz_tol: float | None = None,
    fragment_eic_mz_tol_unit: str | None = None,
    fragment_correlation_mode: str = "full_window",
    correlation_min_relative_intensity: float = 0.05,
    min_correlation_scans: int = 5,
    max_fragment_apex_offset_scans: int = 1,
    min_consecutive_fragment_scans: int = 3,
):
    """公开的单RT窗口提取包装，便于测试和细粒度复用。

    参数语义与 :func:`ms2.extraction.extract_window_result` 相同；该包装保留
    旧调用签名，并把独立前体/fragment EIC容差转交给底层实现。
    """

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
        precursor_mz_tol=precursor_eic_mz_tol,
        precursor_mz_tol_unit=precursor_eic_mz_tol_unit,
        fragment_eic_mz_tol=fragment_eic_mz_tol,
        fragment_eic_mz_tol_unit=fragment_eic_mz_tol_unit,
        fragment_correlation_mode=fragment_correlation_mode,
        correlation_min_relative_intensity=correlation_min_relative_intensity,
        min_correlation_scans=min_correlation_scans,
        max_fragment_apex_offset_scans=max_fragment_apex_offset_scans,
        min_consecutive_fragment_scans=min_consecutive_fragment_scans,
    )


def _initialize_chiral_columns(rows: list[dict]):
    """为每一行预建两侧谱、apex、面积、计数和质量标记列。"""

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
    """预建RT窗口、DIA、相似度及分层诊断列，保证CSV schema稳定。"""

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
        "ms1_reference_status",
        "ms1_reference_issue_codes",
        "ms2_diagnostic_status",
        "ms2_issue_codes",
        "ms2_review_priority",
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
    """把一个RT窗口的提取结果写回``peak_a_*``或``peak_b_*``列。"""

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
    """写入数值相似度指标，再调用统一状态分层逻辑。"""

    row["ms2_cosine"] = csv_scalar(similarity.cosine)
    row["ms2_entropy_similarity"] = csv_scalar(similarity.entropy_similarity)
    row["ms2_matched_peaks"] = csv_scalar(similarity.matched_peaks)
    row["peak_a_fragment_count"] = csv_scalar(similarity.peak_a_count)
    row["peak_b_fragment_count"] = csv_scalar(similarity.peak_b_count)
    row["peak_a_explained_intensity"] = csv_scalar(similarity.peak_a_explained_intensity)
    row["peak_b_explained_intensity"] = csv_scalar(similarity.peak_b_explained_intensity)
    _set_status_layers(row, status, reason, similarity=similarity)


def _set_status_layers(
    row: dict,
    status: str,
    reason: str,
    *,
    chromatographic_status: str = "input_double_peak",
    similarity: SpectrumSimilarity | None = None,
):
    """分别记录色谱、MS2身份支持、候选双峰和富集结论。

    当前分析没有input/eluate重复及统计检验，因此无论MS2是否支持同一化合物，
    ``enantioselective_enrichment_status``都明确写为未评价。
    """

    row["enantiomer_pair_status"] = status
    row["enantiomer_pair_reason"] = reason
    row["chromatographic_status"] = chromatographic_status
    diagnostic = classify_ms2_diagnostic(
        similarity,
        status,
        reason,
        peak_a_quality_flags=row.get("peak_a_quality_flags", ""),
        peak_b_quality_flags=row.get("peak_b_quality_flags", ""),
        dia_window_match_count=row.get("dia_window_match_count"),
    )
    row["ms2_diagnostic_status"] = diagnostic.status
    row["ms2_issue_codes"] = ";".join(diagnostic.issue_codes)
    row["ms2_review_priority"] = diagnostic.review_priority
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
    """优先精确匹配表头；只有失败时才进行不区分大小写匹配。"""

    if requested in rows[0]:
        return requested
    return column(rows, requested)


def _rt_to_seconds(value: float, unit: str) -> float:
    """把单个RT值从分钟或秒统一转换为秒。"""

    if unit == "min":
        return value * 60
    if unit == "sec":
        return value
    raise ValueError("rt_unit must be either 'min' or 'sec'.")


def _is_finite_positive(value) -> bool:
    """判断参数是否可转换为严格正的有限浮点数。"""

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
    method_profile_path: str | Path | None = None,
):
    """写出可复现运行所需的JSON sidecar。

    sidecar记录参数、状态计数、Python版本、耗时、输入文件指纹、raw实测采集
    信息、参考方法profile以及逐个DIA窗口边界。参考方法只作为独立来源层，
    永远不会覆盖raw-observed元数据。
    """

    counts: dict[str, int] = {}
    ms2_diagnostic_counts: dict[str, int] = {}
    ms1_reference_counts: dict[str, int] = {}
    for row in rows:
        status = row.get("enantiomer_pair_status", "")
        counts[status] = counts.get(status, 0) + 1
        ms2_status = row.get("ms2_diagnostic_status", "")
        ms2_diagnostic_counts[ms2_status] = ms2_diagnostic_counts.get(ms2_status, 0) + 1
        ms1_status = row.get("ms1_reference_status", "")
        ms1_reference_counts[ms1_status] = ms1_reference_counts.get(ms1_status, 0) + 1
    method_profile = {}
    resolved_method_profile_path = None
    if method_profile_path:
        resolved_method_profile_path = Path(method_profile_path).resolve()
        with resolved_method_profile_path.open(encoding="utf-8") as handle:
            method_profile = json.load(handle)
        if not isinstance(method_profile, dict) or not method_profile.get("profile_id"):
            raise ValueError(
                f"Method profile must define profile_id: {resolved_method_profile_path}"
            )
    metadata = {
        "schema_version": "msai-chiral-ms2-v3-separated-diagnostics",
        "created_utc": datetime.now(UTC).isoformat(),
        "python": platform.python_version(),
        "elapsed_seconds": elapsed_seconds,
        "rows_processed": len(rows),
        "status_counts": counts,
        "ms2_diagnostic_status_counts": ms2_diagnostic_counts,
        "ms1_reference_status_counts": ms1_reference_counts,
        "parameters": parameters,
        "provenance_layers": {
            "analysis_parameters": parameters,
            "raw_observed": raw_acquisition_context(raw_path),
            "reference_method": method_profile,
            "reference_method_path": (
                str(resolved_method_profile_path) if resolved_method_profile_path else None
            ),
            "reference_method_sha256": (
                _file_fingerprint(resolved_method_profile_path)["sha256"]
                if resolved_method_profile_path
                else None
            ),
            "rule": "reference method is contextual evidence and never overwrites raw-observed acquisition metadata",
        },
        "dia_window_provenance": {
            "raw_bounds_window_count": sum(
                spectra_index.windows[center].isolation_lower_mz is not None
                and spectra_index.windows[center].isolation_upper_mz is not None
                for center in spectra_index.precursors
            ),
            "fallback_window_count": sum(
                spectra_index.windows[center].isolation_lower_mz is None
                or spectra_index.windows[center].isolation_upper_mz is None
                for center in spectra_index.precursors
            ),
            "fallback_width_mz": parameters.get("dia_iso_win_fallback"),
            "rule": "raw isolation bounds are authoritative; fallback is used only when either bound is missing",
        },
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
    """返回文件路径、大小、mtime和分块计算的SHA-256指纹。"""

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


def main(argv=None):
    """解析``analyze``命令参数，运行正式逐行MS2分析并打印状态计数。"""

    parser = argparse.ArgumentParser(
        description="Extract and compare MS2 spectra for per-row chiral LC peak pairs."
    )
    parser.add_argument(
        "--peaklist", required=True, help="CSV/XLSX containing MZ, Peak1, and Peak2 columns."
    )
    parser.add_argument("--raw", required=True, help="Matching mzXML/mzML raw MS file.")
    parser.add_argument("--output", help="Output CSV; defaults to the project results directory.")
    parser.add_argument("--mz-column", default="MZ")
    parser.add_argument("--peak-a-column", default="Peak1")
    parser.add_argument("--peak-b-column", default="Peak2")
    parser.add_argument("--rt-unit", choices=("min", "sec"), default="min")
    parser.add_argument(
        "--rt-half-window-sec", type=float, default=DEFAULT_RT_HALF_WINDOW_SEC
    )
    parser.add_argument("--mz-tol", type=float, default=10.0)
    parser.add_argument("--mz-tol-unit", choices=("ppm", "Da", "legacy_fraction"), default="ppm")
    parser.add_argument(
        "--precursor-eic-mz-tol",
        type=float,
        help="Optional precursor EIC tolerance; defaults to --mz-tol.",
    )
    parser.add_argument(
        "--precursor-eic-mz-tol-unit",
        choices=("ppm", "Da", "legacy_fraction"),
        help="Defaults to --mz-tol-unit.",
    )
    parser.add_argument(
        "--fragment-eic-mz-tol",
        type=float,
        help="Optional fragment coelution-EIC/centroid merge tolerance; defaults to --mz-tol.",
    )
    parser.add_argument(
        "--fragment-eic-mz-tol-unit",
        choices=("ppm", "Da", "legacy_fraction"),
        help="Defaults to --mz-tol-unit.",
    )
    parser.add_argument(
        "--dia-window",
        "--dia-window-fallback",
        dest="dia_window",
        type=float,
        default=15.0,
        help="Fallback DIA width used only when raw isolation bounds are missing.",
    )
    parser.add_argument("--fragment-mz-tol", type=float, default=0.01)
    parser.add_argument("--fragment-mz-tol-unit", choices=("Da", "ppm"), default="Da")
    parser.add_argument("--min-relative-intensity", type=float, default=0.01)
    parser.add_argument("--min-fragment-intensity", type=float, default=2000.0)
    parser.add_argument("--min-fragment-relative-intensity", type=float, default=0.0)
    parser.add_argument("--min-fragment-correlation", type=float, default=0.9)
    parser.add_argument(
        "--fragment-correlation-mode",
        choices=("full_window", "active_support"),
        default="full_window",
        help="Use legacy full-window Pearson or a common-zero-resistant active-support Pearson.",
    )
    parser.add_argument("--correlation-min-relative-intensity", type=float, default=0.05)
    parser.add_argument("--min-correlation-scans", type=int, default=5)
    parser.add_argument("--max-fragment-apex-offset-scans", type=int, default=1)
    parser.add_argument("--min-consecutive-fragment-scans", type=int, default=3)
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
    parser.add_argument(
        "--method-profile",
        help="Optional reference-method JSON stored separately from raw-observed metadata.",
    )
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
        precursor_eic_mz_tol=args.precursor_eic_mz_tol,
        precursor_eic_mz_tol_unit=args.precursor_eic_mz_tol_unit,
        fragment_eic_mz_tol=args.fragment_eic_mz_tol,
        fragment_eic_mz_tol_unit=args.fragment_eic_mz_tol_unit,
        dia_iso_win=args.dia_window,
        fragment_mz_tol=args.fragment_mz_tol,
        fragment_mz_tol_unit=args.fragment_mz_tol_unit,
        min_relative_intensity=args.min_relative_intensity,
        min_fragment_intensity=args.min_fragment_intensity,
        min_fragment_relative_intensity=args.min_fragment_relative_intensity,
        min_fragment_correlation=args.min_fragment_correlation,
        fragment_correlation_mode=args.fragment_correlation_mode,
        correlation_min_relative_intensity=args.correlation_min_relative_intensity,
        min_correlation_scans=args.min_correlation_scans,
        max_fragment_apex_offset_scans=args.max_fragment_apex_offset_scans,
        min_consecutive_fragment_scans=args.min_consecutive_fragment_scans,
        consensus_scans=args.consensus_scans,
        thresholds=thresholds,
        write_metadata=not args.no_metadata,
        method_profile_path=args.method_profile,
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
        main()
