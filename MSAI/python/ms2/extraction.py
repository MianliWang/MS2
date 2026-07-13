"""One-feature, one-RT-window DIA fragment extraction workflow."""

from __future__ import annotations

import math
from bisect import bisect_left, bisect_right

from .chromatograms import (
    _ms2copy,
    _raw_eic,
    _split_flags,
    empty_window_result,
    format_fragment_string,
    merge_window_result,
    summarize_xic_peak,
)
from .models import Ms2Index, Spectrum
from .utils import mz_bounds, pearson_correlation, sample_sd


def extract_window_result(
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
    precursor_mz_tol: float | None = None,
    precursor_mz_tol_unit: str | None = None,
    fragment_eic_mz_tol: float | None = None,
    fragment_eic_mz_tol_unit: str | None = None,
    fragment_correlation_mode: str = "full_window",
    correlation_min_relative_intensity: float = 0.05,
    min_correlation_scans: int = 5,
    max_fragment_apex_offset_scans: int = 1,
    min_consecutive_fragment_scans: int = 3,
) -> dict:
    """在一个目标、一个DIA窗口和一个RT窗口内重建共洗脱碎片谱。

    处理步骤：先提取目标前体EIC并确定apex；再从最强的一个或多个consensus
    scans收集低于前体至少10 Da的候选fragment；随后为每个fragment重新提取
    色谱轨迹，并用强度、连续支持和前体/碎片共洗脱相关性筛选。返回值同时
    包含谱字符串、全精度碎片数组、apex、面积、计数和质量标记。

    ``mz_tol``是向后兼容的共享容差。独立的``precursor_mz_tol``与
    ``fragment_eic_mz_tol``允许前体提取和fragment centroid合并分别校准，
    避免修改一个参数时无意改变另一阶段。该函数只重建证据，不判断两张谱
    是否属于同一化合物。
    """

    precursor_mz_tol = mz_tol if precursor_mz_tol is None else precursor_mz_tol
    precursor_mz_tol_unit = mz_tol_unit if precursor_mz_tol_unit is None else precursor_mz_tol_unit
    fragment_eic_mz_tol = mz_tol if fragment_eic_mz_tol is None else fragment_eic_mz_tol
    fragment_eic_mz_tol_unit = (
        mz_tol_unit if fragment_eic_mz_tol_unit is None else fragment_eic_mz_tol_unit
    )

    if min_fragment_intensity < 0 or not math.isfinite(min_fragment_intensity):
        raise ValueError("min_fragment_intensity must be finite and non-negative.")
    if not 0 <= min_fragment_relative_intensity <= 1:
        raise ValueError("min_fragment_relative_intensity must be between 0 and 1.")
    if not -1 <= min_fragment_correlation <= 1:
        raise ValueError("min_fragment_correlation must be between -1 and 1.")
    if consensus_scans < 1:
        raise ValueError("consensus_scans must be at least 1.")
    if fragment_correlation_mode not in {"full_window", "active_support"}:
        raise ValueError("fragment_correlation_mode must be 'full_window' or 'active_support'.")
    if not 0 <= correlation_min_relative_intensity < 1:
        raise ValueError("correlation_min_relative_intensity must be in [0, 1).")
    if min_correlation_scans < 2:
        raise ValueError("min_correlation_scans must be at least 2.")
    if max_fragment_apex_offset_scans < 0:
        raise ValueError("max_fragment_apex_offset_scans must be non-negative.")
    if min_consecutive_fragment_scans < 1:
        raise ValueError("min_consecutive_fragment_scans must be at least 1.")

    dia = _ms2copy(spectra, diawin)
    if not dia.spectra:
        return empty_window_result("no_ms2_scans")

    mzmin, mzmax = _clamped_bounds(precurmz, precursor_mz_tol, precursor_mz_tol_unit, dia.mzrange)
    if mzmin >= mzmax:
        return empty_window_result("mz_out_of_range")

    rtmin = max(dia.rts[0], rt_window_sec[0])
    rtmax = min(dia.rts[-1], rt_window_sec[1])
    if rtmin >= rtmax:
        return empty_window_result("rt_window_out_of_range")

    scan_count = bisect_right(dia.rts, rtmax) - bisect_left(dia.rts, rtmin)
    if scan_count < 1:
        return empty_window_result("no_ms2_scans", 0)

    # “native”是目标前体在当前DIA窗口内的EIC，后续所有fragment都必须
    # 在相同RT采样点上与它比较，才能判断是否真正共同洗脱。
    native = _raw_eic(dia, mzmin, mzmax, rtmin, rtmax)
    peak_summary = summarize_xic_peak(native)
    flags = _split_flags(peak_summary["quality_flags"])
    if not native.scan or not native.intensity:
        return merge_window_result(peak_summary, "", scan_count, [*flags, "no_fragment_scan"])
    if max(native.intensity) <= 0:
        return merge_window_result(peak_summary, "", scan_count, [*flags, "no_precursor_signal"])

    apex_position = max(range(len(native.intensity)), key=lambda index: native.intensity[index])
    apex_scan = native.scan[apex_position]
    if apex_scan >= len(dia.spectra):
        return merge_window_result(peak_summary, "", scan_count, [*flags, "invalid_apex_scan"])

    # consensus_scans>1时合并若干最强前体扫描的候选fragment，提高低丰度
    # 碎片召回；最终强度仍来自完整RT窗口的fragment trace最大值。
    selected_positions = sorted(
        range(len(native.intensity)),
        key=lambda position: native.intensity[position],
        reverse=True,
    )
    selected_scans = [
        native.scan[position] for position in selected_positions if native.intensity[position] > 0
    ][:consensus_scans]
    candidate_peaks = [
        (mz, intensity)
        for scan_index in selected_scans
        for mz, intensity in zip(
            dia.spectra[scan_index].mz,
            dia.spectra[scan_index].intensity,
            strict=True,
        )
        if mz < precurmz - 10 and intensity > 0
    ]
    candidates = _merge_candidate_peaks(
        candidate_peaks, fragment_eic_mz_tol, fragment_eic_mz_tol_unit
    )
    if not candidates:
        return merge_window_result(
            peak_summary,
            "",
            scan_count,
            [*flags, "no_candidate_fragments"],
            consensus_scan_count=len(selected_scans),
        )

    fragment_mz: list[float] = []
    fragment_intensity: list[float] = []
    native_sd = sample_sd(native.intensity)
    correlated: list[tuple[float, float]] = []
    for mz0 in candidates:
        frag_mzmin, frag_mzmax = _clamped_bounds(
            mz0, fragment_eic_mz_tol, fragment_eic_mz_tol_unit, dia.mzrange
        )
        if frag_mzmin >= frag_mzmax:
            continue
        fragment_trace = _raw_eic(dia, frag_mzmin, frag_mzmax, rtmin, rtmax).intensity
        if len(fragment_trace) != len(native.intensity):
            continue
        fragment_sd = sample_sd(fragment_trace)
        if native_sd == 0 or fragment_sd == 0:
            continue
        fragment_max = max(fragment_trace)
        if fragment_max < min_fragment_intensity:
            continue
        correlation = _trace_correlation(
            native.intensity,
            fragment_trace,
            mode=fragment_correlation_mode,
            min_relative_intensity=correlation_min_relative_intensity,
            min_scans=min_correlation_scans,
            max_apex_offset_scans=max_fragment_apex_offset_scans,
            min_consecutive_fragment_scans=min_consecutive_fragment_scans,
        )
        if correlation is not None and correlation > min_fragment_correlation:
            correlated.append((mz0, fragment_max))

    relative_cutoff = (
        max((intensity for _mz, intensity in correlated), default=0.0)
        * min_fragment_relative_intensity
    )
    for fragment, intensity in correlated:
        if intensity >= relative_cutoff:
            fragment_mz.append(fragment)
            fragment_intensity.append(intensity)

    ms2 = format_fragment_string(fragment_mz, fragment_intensity)
    if not ms2:
        flags.append("no_fragments")
    return merge_window_result(
        peak_summary,
        ms2,
        scan_count,
        flags,
        consensus_scan_count=len(selected_scans),
        candidate_fragment_count=len(candidates),
        fragment_count=len(fragment_mz),
        fragment_peaks=list(zip(fragment_mz, fragment_intensity, strict=True)),
    )


def _merge_candidate_peaks(
    peaks: list[tuple[float, float]],
    mz_tol: float,
    mz_tol_unit: str,
) -> list[float]:
    """合并多个apex邻近扫描中重复出现的centroid fragment。

    每个cluster以强度加权m/z为中心，并使用与fragment EIC相同的容差决定
    下一个centroid是否仍属于该cluster。
    """

    if not peaks:
        return []
    merged: list[float] = []
    cluster: list[tuple[float, float]] = []
    for mz, intensity in sorted(peaks):
        if cluster:
            total = sum(value for _mass, value in cluster)
            center = sum(mass * value for mass, value in cluster) / total
            low, high = mz_bounds(center, mz_tol, mz_tol_unit)
            if not low <= mz <= high:
                merged.append(center)
                cluster = []
        cluster.append((mz, intensity))
    if cluster:
        total = sum(value for _mass, value in cluster)
        merged.append(sum(mass * value for mass, value in cluster) / total)
    return merged


def _clamped_bounds(
    mz: float,
    mz_tol: float,
    mz_tol_unit: str,
    mzrange: tuple[float, float],
) -> tuple[float, float]:
    """计算目标m/z容差边界，并裁剪到DIA数据实际m/z范围。"""

    low, high = mz_bounds(mz, mz_tol, mz_tol_unit)
    return max(mzrange[0], low), min(mzrange[1], high)


def _trace_correlation(
    precursor: list[float],
    fragment: list[float],
    *,
    mode: str,
    min_relative_intensity: float,
    min_scans: int,
    max_apex_offset_scans: int,
    min_consecutive_fragment_scans: int,
) -> float | None:
    """计算前体与fragment色谱轨迹的相关性。

    ``full_window``复现旧版整窗Pearson；``active_support``先要求apex位置接近、
    前体有效支持点足够且fragment连续出现，再只在前体活跃区计算，从而减少
    窗口两端大量共同零值导致的虚高相关。
    """

    if len(precursor) != len(fragment) or not precursor:
        return None
    if mode == "full_window":
        return pearson_correlation(precursor, fragment)

    precursor_max = max(precursor, default=0.0)
    fragment_max = max(fragment, default=0.0)
    if precursor_max <= 0 or fragment_max <= 0:
        return None
    precursor_apex = max(range(len(precursor)), key=precursor.__getitem__)
    fragment_apex = max(range(len(fragment)), key=fragment.__getitem__)
    if abs(precursor_apex - fragment_apex) > max_apex_offset_scans:
        return None

    precursor_cutoff = precursor_max * min_relative_intensity
    fragment_cutoff = fragment_max * min_relative_intensity
    support = [index for index, value in enumerate(precursor) if value >= precursor_cutoff]
    if len(support) < min_scans:
        return None
    if (
        _longest_consecutive_run(index for index in support if fragment[index] >= fragment_cutoff)
        < min_consecutive_fragment_scans
    ):
        return None
    return pearson_correlation(
        [precursor[index] for index in support],
        [fragment[index] for index in support],
    )


def _longest_consecutive_run(indices) -> int:
    """返回升序扫描下标中最长的连续区段长度。"""

    best = current = 0
    previous = None
    for index in indices:
        current = current + 1 if previous is not None and index == previous + 1 else 1
        best = max(best, current)
        previous = index
    return best


__all__ = ["extract_window_result"]
