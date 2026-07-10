"""Small data structures shared by the MS2 pipeline."""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass
class Spectrum:
    """Minimal centroided MS2 scan representation."""

    precursor_mz: float
    rt: float
    mz: list[float]
    intensity: list[float]
    isolation_lower_mz: float | None = None
    isolation_upper_mz: float | None = None


@dataclass
class DiaData:
    """MS2 scans and bounds belonging to one DIA isolation window."""

    spectra: list[Spectrum]
    mzrange: tuple[float, float]
    rts: list[float]
    isolation_lower_mz: float | None = None
    isolation_upper_mz: float | None = None


@dataclass
class Ms2Index:
    """DIA-window and retention-time index built once per raw file."""

    windows: dict[float, DiaData]
    precursors: list[float]

    def get(self, precursor_mz: float) -> DiaData:
        return self.windows.get(precursor_mz, DiaData([], (math.inf, -math.inf), []))

    def matching_windows(self, mz: float, fallback_width: float) -> list[float]:
        matches: list[float] = []
        for center in self.precursors:
            window = self.windows[center]
            lower = window.isolation_lower_mz
            upper = window.isolation_upper_mz
            if lower is None or upper is None:
                lower, upper = center - fallback_width / 2, center + fallback_width / 2
            if lower <= mz <= upper:
                matches.append(center)
        return sorted(matches, key=lambda center: (abs(center - mz), center))


@dataclass
class Eic:
    """Extracted-ion chromatogram values indexed by scan and RT seconds."""

    rt: list[float]
    scan: list[int]
    intensity: list[float]


__all__ = ["DiaData", "Eic", "Ms2Index", "Spectrum"]
