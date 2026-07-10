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
) -> dict:
    """Extract coeluting fragments for one feature and one RT window."""

    if min_fragment_intensity < 0 or not math.isfinite(min_fragment_intensity):
        raise ValueError("min_fragment_intensity must be finite and non-negative.")
    if not 0 <= min_fragment_relative_intensity <= 1:
        raise ValueError("min_fragment_relative_intensity must be between 0 and 1.")
    if not -1 <= min_fragment_correlation <= 1:
        raise ValueError("min_fragment_correlation must be between -1 and 1.")
    if consensus_scans < 1:
        raise ValueError("consensus_scans must be at least 1.")

    dia = _ms2copy(spectra, diawin)
    if not dia.spectra:
        return empty_window_result("no_ms2_scans")

    mzmin, mzmax = _clamped_bounds(precurmz, mz_tol, mz_tol_unit, dia.mzrange)
    if mzmin >= mzmax:
        return empty_window_result("mz_out_of_range")

    rtmin = max(dia.rts[0], rt_window_sec[0])
    rtmax = min(dia.rts[-1], rt_window_sec[1])
    if rtmin >= rtmax:
        return empty_window_result("rt_window_out_of_range")

    scan_count = bisect_right(dia.rts, rtmax) - bisect_left(dia.rts, rtmin)
    if scan_count < 1:
        return empty_window_result("no_ms2_scans", 0)

    native = _raw_eic(dia, mzmin, mzmax, rtmin, rtmax)
    peak_summary = summarize_xic_peak(native)
    flags = _split_flags(peak_summary["quality_flags"])
    if not native.scan or not native.intensity:
        return merge_window_result(peak_summary, "", scan_count, flags + ["no_fragment_scan"])
    if max(native.intensity) <= 0:
        return merge_window_result(peak_summary, "", scan_count, flags + ["no_precursor_signal"])

    apex_position = max(range(len(native.intensity)), key=lambda index: native.intensity[index])
    apex_scan = native.scan[apex_position]
    if apex_scan >= len(dia.spectra):
        return merge_window_result(peak_summary, "", scan_count, flags + ["invalid_apex_scan"])

    selected_positions = sorted(
        range(len(native.intensity)),
        key=lambda position: native.intensity[position],
        reverse=True,
    )
    selected_scans = [
        native.scan[position]
        for position in selected_positions
        if native.intensity[position] > 0
    ][:consensus_scans]
    candidate_peaks = [
        (mz, intensity)
        for scan_index in selected_scans
        for mz, intensity in zip(
            dia.spectra[scan_index].mz,
            dia.spectra[scan_index].intensity,
        )
        if mz < precurmz - 10 and intensity > 0
    ]
    candidates = _merge_candidate_peaks(candidate_peaks, mz_tol, mz_tol_unit)
    if not candidates:
        return merge_window_result(
            peak_summary,
            "",
            scan_count,
            flags + ["no_candidate_fragments"],
            consensus_scan_count=len(selected_scans),
        )

    fragment_mz: list[float] = []
    fragment_intensity: list[float] = []
    native_sd = sample_sd(native.intensity)
    correlated: list[tuple[float, float]] = []
    for mz0 in candidates:
        frag_mzmin, frag_mzmax = _clamped_bounds(mz0, mz_tol, mz_tol_unit, dia.mzrange)
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
        if pearson_correlation(native.intensity, fragment_trace) > min_fragment_correlation:
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
        fragment_peaks=list(zip(fragment_mz, fragment_intensity)),
    )


def _merge_candidate_peaks(
    peaks: list[tuple[float, float]],
    mz_tol: float,
    mz_tol_unit: str,
) -> list[float]:
    """Merge repeated centroid peaks from several apex-adjacent scans."""

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
    low, high = mz_bounds(mz, mz_tol, mz_tol_unit)
    return max(mzrange[0], low), min(mzrange[1], high)


__all__ = ["extract_window_result"]
