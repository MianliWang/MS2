"""DIA extracted-ion chromatograms and stable result structures."""

from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from itertools import pairwise

from .indexing import _ensure_sorted_peaks
from .models import DiaData, Eic, Ms2Index, Spectrum


def summarize_xic_peak(eic: Eic) -> dict:
    """汇总一条EIC的apex RT、apex强度和梯形积分面积。

    面积直接对实际RT间隔积分，不假设scan间隔完全相等。有效点不足两个时
    面积为空并标记``insufficient_points``，避免把单点强度误作峰面积。
    """

    if not eic.intensity:
        return _empty_xic_summary("empty_xic")
    apex_index = max(
        range(len(eic.intensity)),
        key=lambda index: (
            eic.intensity[index] if math.isfinite(eic.intensity[index]) else -math.inf
        ),
    )
    valid = [
        index
        for index, (rt, intensity) in enumerate(zip(eic.rt, eic.intensity, strict=True))
        if math.isfinite(rt) and math.isfinite(intensity)
    ]
    if not valid:
        return _empty_xic_summary("missing_rt")
    if len(valid) < 2:
        area = None
        flags = "insufficient_points"
    else:
        pairs = sorted((eic.rt[index], eic.intensity[index]) for index in valid)
        area = sum(
            (right[0] - left[0]) * (left[1] + right[1]) / 2 for left, right in pairwise(pairs)
        )
        flags = "ok"
    return {
        "apex_rt": eic.rt[apex_index],
        "apex_intensity": eic.intensity[apex_index],
        "area": area,
        "quality_flags": flags,
    }


def format_fragment_string(fragment_mz, intensity) -> str:
    """将fragment m/z和强度序列化为项目兼容的单元格文本。

    格式为``mz,intensity;mz,intensity``，使用足够有效数字以避免可视化或
    再比较时因过早四舍五入改变匹配结果。
    """

    pairs = [
        (mz, value)
        for mz, value in zip(fragment_mz, intensity, strict=True)
        if mz is not None and value is not None
    ]
    if not pairs:
        return ""
    return ";".join(f"{mz:.12g},{value:.12g}" for mz, value in pairs)


def _ms2copy(spectra: list[Spectrum] | Ms2Index, diawin: float) -> DiaData:
    """取得一个DIA窗口的扫描并按RT排序；已有索引时为常数时间查找。"""

    if isinstance(spectra, Ms2Index):
        return spectra.get(diawin)
    picked = [spectrum for spectrum in spectra if spectrum.precursor_mz == diawin]
    for spectrum in picked:
        _ensure_sorted_peaks(spectrum)
    picked.sort(key=lambda spectrum: spectrum.rt)
    lows = [spectrum.mz[0] for spectrum in picked if spectrum.mz]
    highs = [spectrum.mz[-1] for spectrum in picked if spectrum.mz]
    if not picked or not lows or not highs:
        return DiaData(picked, (math.inf, -math.inf), [spectrum.rt for spectrum in picked])
    return DiaData(picked, (min(lows), max(highs)), [spectrum.rt for spectrum in picked])


def _raw_eic(
    dia: DiaData,
    mzmin: float,
    mzmax: float,
    rtmin: float,
    rtmax: float,
) -> Eic:
    """在闭合m/z和RT范围内提取求和强度EIC。

    先对RT索引、再对每张谱的m/z数组做二分查找，从而避免逐峰全扫描。
    返回的``scan``是DIA窗口内部下标，供后续定位真实apex scan使用。
    """

    rt: list[float] = []
    scan: list[int] = []
    intensity: list[float] = []
    first = bisect_left(dia.rts, rtmin)
    last = bisect_right(dia.rts, rtmax)
    for index in range(first, last):
        spectrum = dia.spectra[index]
        rt.append(spectrum.rt)
        scan.append(index)
        mz_first = bisect_left(spectrum.mz, mzmin)
        mz_last = bisect_right(spectrum.mz, mzmax)
        intensity.append(sum(spectrum.intensity[mz_first:mz_last]))
    return Eic(rt=rt, scan=scan, intensity=intensity)


def _empty_xic_summary(flag: str) -> dict:
    """生成字段完整的空EIC摘要，并保留无法计算的原因。"""

    return {
        "apex_rt": None,
        "apex_intensity": None,
        "area": None,
        "quality_flags": flag,
    }


def empty_window_result(flag: str, ms2_count=None) -> dict:
    """为不可评价的RT/DIA窗口返回字段稳定的空结果。"""

    return {
        "MS2": "",
        "fragment_peaks": [],
        "apex_rt": None,
        "apex_intensity": None,
        "area": None,
        "ms2_count": ms2_count,
        "consensus_scan_count": None,
        "candidate_fragment_count": None,
        "fragment_count": None,
        "quality_flags": flag,
    }


def merge_window_result(
    peak_summary: dict,
    ms2: str,
    ms2_count: int,
    flags: list[str],
    *,
    consensus_scan_count: int | None = None,
    candidate_fragment_count: int | None = None,
    fragment_count: int | None = None,
    fragment_peaks: list[tuple[float, float]] | None = None,
) -> dict:
    """合并EIC摘要、碎片谱、扫描计数和去重后的质量标记。"""

    clean_flags: list[str] = []
    for flag in flags:
        if flag and flag != "ok" and flag not in clean_flags:
            clean_flags.append(flag)
    if not clean_flags:
        clean_flags = ["ok"]
    return {
        "MS2": ms2,
        "fragment_peaks": fragment_peaks or [],
        "apex_rt": peak_summary["apex_rt"],
        "apex_intensity": peak_summary["apex_intensity"],
        "area": peak_summary["area"],
        "ms2_count": ms2_count,
        "consensus_scan_count": consensus_scan_count,
        "candidate_fragment_count": candidate_fragment_count,
        "fragment_count": fragment_count,
        "quality_flags": ";".join(clean_flags),
    }


def _split_flags(flags: str) -> list[str]:
    """把分号分隔质量标记还原为列表；``ok``等价于没有问题。"""

    if not flags or flags == "ok":
        return []
    return flags.split(";")


__all__ = [
    "empty_window_result",
    "format_fragment_string",
    "merge_window_result",
    "summarize_xic_peak",
]
