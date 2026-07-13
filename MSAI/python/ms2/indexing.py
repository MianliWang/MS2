"""DIA-window grouping and lookup."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Iterable
from itertools import pairwise

from .models import DiaData, Ms2Index, Spectrum


def build_ms2_index(spectra: Iterable[Spectrum]) -> Ms2Index:
    """把MS2扫描按DIA中心分组，并为后续批量目标建立一次性索引。

    每个窗口内部按RT排序，同时汇总实际fragment m/z范围和原始隔离边界。
    对同一raw文件分析数百行目标时，应复用返回的``Ms2Index``，而不是为
    每一行重新分组，这是当前速度优化的关键之一。
    """

    grouped: dict[float, list[Spectrum]] = defaultdict(list)
    order: list[float] = []
    for spectrum in spectra:
        _ensure_sorted_peaks(spectrum)
        if spectrum.precursor_mz not in grouped:
            order.append(spectrum.precursor_mz)
        grouped[spectrum.precursor_mz].append(spectrum)

    windows: dict[float, DiaData] = {}
    for precursor_mz, picked in grouped.items():
        picked.sort(key=lambda spectrum: spectrum.rt)
        lows = [spectrum.mz[0] for spectrum in picked if spectrum.mz]
        highs = [spectrum.mz[-1] for spectrum in picked if spectrum.mz]
        mzrange = (min(lows), max(highs)) if lows and highs else (math.inf, -math.inf)
        isolation_lows = [
            spectrum.isolation_lower_mz
            for spectrum in picked
            if spectrum.isolation_lower_mz is not None
        ]
        isolation_highs = [
            spectrum.isolation_upper_mz
            for spectrum in picked
            if spectrum.isolation_upper_mz is not None
        ]
        windows[precursor_mz] = DiaData(
            picked,
            mzrange,
            [spectrum.rt for spectrum in picked],
            min(isolation_lows) if isolation_lows else None,
            max(isolation_highs) if isolation_highs else None,
        )
    return Ms2Index(windows=windows, precursors=order)


def _ensure_sorted_peaks(spectrum: Spectrum) -> None:
    """保证m/z升序且强度仍与其配对，以满足二分截取的前提。"""

    if len(spectrum.mz) != len(spectrum.intensity):
        raise ValueError("Spectrum m/z and intensity arrays must have equal length.")
    if any(left > right for left, right in pairwise(spectrum.mz)):
        pairs = sorted(zip(spectrum.mz, spectrum.intensity, strict=True))
        spectrum.mz = [mz for mz, _intensity in pairs]
        spectrum.intensity = [intensity for _mz, intensity in pairs]


def preclist(spectra: Iterable[Spectrum]) -> list[float]:
    """返回DIA窗口中心并保持它们在原始文件中的首次出现顺序。"""

    seen: list[float] = []
    for spectrum in spectra:
        if spectrum.precursor_mz not in seen:
            seen.append(spectrum.precursor_mz)
    return seen


def matching_dia_windows(
    mz: float,
    precursor: list[float] | Ms2Index,
    dia_iso_win: float,
) -> list[float]:
    """返回覆盖目标m/z的所有已采集窗口，最近中心排在最前。

    传入``Ms2Index``时优先使用raw中的真实隔离边界；传入旧式中心列表时
    才按``dia_iso_win``名义宽度计算。
    """

    if isinstance(precursor, Ms2Index):
        return precursor.matching_windows(mz, dia_iso_win)
    return sorted(
        (value for value in precursor if abs(mz - value) <= 0.5 * dia_iso_win),
        key=lambda value: (abs(value - mz), value),
    )


def first_dia_window(
    mz: float,
    precursor: list[float] | Ms2Index,
    dia_iso_win: float,
) -> float | None:
    """选择覆盖目标的最近DIA窗口；没有覆盖时返回``None``。"""

    matches = matching_dia_windows(mz, precursor, dia_iso_win)
    return matches[0] if matches else None


__all__ = ["build_ms2_index", "first_dia_window", "matching_dia_windows", "preclist"]
