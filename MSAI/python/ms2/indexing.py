"""DIA-window grouping and lookup."""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Iterable

from .models import DiaData, Ms2Index, Spectrum


def build_ms2_index(spectra: Iterable[Spectrum]) -> Ms2Index:
    """Group spectra once instead of rebuilding a DIA copy for every row."""

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
    """Keep m/z and intensity paired while enforcing the bisect invariant."""

    if len(spectrum.mz) != len(spectrum.intensity):
        raise ValueError("Spectrum m/z and intensity arrays must have equal length.")
    if any(left > right for left, right in zip(spectrum.mz, spectrum.mz[1:])):
        pairs = sorted(zip(spectrum.mz, spectrum.intensity))
        spectrum.mz = [mz for mz, _intensity in pairs]
        spectrum.intensity = [intensity for _mz, intensity in pairs]


def preclist(spectra: Iterable[Spectrum]) -> list[float]:
    """Return DIA precursor/window centers in first-seen order."""

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
    """Return all acquired windows containing m/z, nearest center first."""

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
    """Select the closest acquired DIA window containing the feature."""

    matches = matching_dia_windows(mz, precursor, dia_iso_win)
    return matches[0] if matches else None


__all__ = ["build_ms2_index", "first_dia_window", "matching_dia_windows", "preclist"]
