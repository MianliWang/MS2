"""Targeted MS1 EIC extraction and automatic chiral double-peak picking.

The defaults mirror the documented ChiralRTlib starting point where the
manuscript is explicit (3 ppm, 1e5 height, 20-scan distance).  Parameters that
were not specified in the manuscript remain configurable and are written to
the output for calibration rather than treated as universal constants.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import xml.etree.ElementTree as ET
from array import array
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

try:
    from .ms2 import (
        _decode_mzml_arrays,
        _decode_mzxml_peaks,
        _local_name,
        _mzml_cv_value,
        _mzml_scan_time,
        _parse_mzxml_duration,
        csv_scalar,
        number,
        read_table,
        write_table,
    )
except ImportError:
    from ms2 import (  # type: ignore
        _decode_mzml_arrays,
        _decode_mzxml_peaks,
        _local_name,
        _mzml_cv_value,
        _mzml_scan_time,
        _parse_mzxml_duration,
        csv_scalar,
        number,
        read_table,
        write_table,
    )


@dataclass(frozen=True)
class PeakPickingConfig:
    eic_ppm: float = 3.0
    sg_window: int = 7
    sg_polyorder: int = 2
    min_height: float | None = 100_000.0
    min_distance_scans: int = 20
    refine_radius_scans: int = 3
    min_scan_points: int = 5
    top_k: int = 2
    min_prominence_ratio: float | None = None
    min_support_scans: int | None = None
    min_distance_sec: float | None = None
    close_peak_max_valley_ratio: float | None = None
    rank_by_prominence: bool = False
    max_valley_ratio: float | None = None
    min_resolution: float | None = 0.0
    min_second_peak_ratio: float = 0.0

    def __post_init__(self):
        if self.eic_ppm <= 0 or not math.isfinite(self.eic_ppm):
            raise ValueError("eic_ppm must be finite and positive.")
        if self.sg_window < 3 or self.sg_window % 2 == 0:
            raise ValueError("sg_window must be an odd integer of at least 3.")
        if self.sg_polyorder < 0 or self.sg_polyorder >= self.sg_window:
            raise ValueError("sg_polyorder must be non-negative and smaller than sg_window.")
        if self.min_height is not None and (
            self.min_height < 0 or not math.isfinite(self.min_height)
        ):
            raise ValueError("min_height must be finite and non-negative when provided.")
        if self.min_distance_scans < 1:
            raise ValueError("min_distance_scans must be at least 1.")
        if self.refine_radius_scans < 0 or self.min_scan_points < 3 or self.top_k < 1:
            raise ValueError("Invalid peak-picking scan-count parameter.")
        if self.max_valley_ratio is not None and not 0 <= self.max_valley_ratio <= 1:
            raise ValueError("max_valley_ratio must be between 0 and 1 when provided.")
        if self.min_prominence_ratio is not None and not 0 <= self.min_prominence_ratio <= 1:
            raise ValueError("min_prominence_ratio must be between 0 and 1 when provided.")
        if self.min_support_scans is not None and self.min_support_scans < 1:
            raise ValueError("min_support_scans must be at least 1 when provided.")
        if self.min_distance_sec is not None and (
            self.min_distance_sec < 0 or not math.isfinite(self.min_distance_sec)
        ):
            raise ValueError("min_distance_sec must be finite and non-negative when provided.")
        if self.close_peak_max_valley_ratio is not None and not (
            0 <= self.close_peak_max_valley_ratio <= 1
        ):
            raise ValueError("close_peak_max_valley_ratio must be between 0 and 1 when provided.")
        if self.min_resolution is not None and (
            self.min_resolution < 0 or not math.isfinite(self.min_resolution)
        ):
            raise ValueError("min_resolution must be finite and non-negative when provided.")
        if not 0 <= self.min_second_peak_ratio <= 1:
            raise ValueError("Invalid peak-quality threshold.")

    @property
    def adaptive_detection(self) -> bool:
        """Whether any optional adaptive candidate-detection feature is enabled."""

        return (
            self.min_height is None
            or self.min_prominence_ratio is not None
            or self.min_support_scans is not None
            or self.min_distance_sec is not None
            or self.close_peak_max_valley_ratio is not None
            or self.rank_by_prominence
        )


@dataclass(frozen=True)
class PickedPeak:
    scan_index: int
    rt_sec: float
    intensity: float
    width_sec: float | None
    area_fwhm: float | None
    snr: float | None
    prominence: float | None = None
    prominence_ratio: float | None = None


@dataclass(frozen=True)
class _PeakCandidate:
    scan_index: int
    prominence: float
    prominence_ratio: float
    support_scans: int


@dataclass(frozen=True)
class PeakPickResult:
    signal_status: str
    chromatographic_status: str
    peaks: tuple[PickedPeak, ...]
    resolution: float | None
    valley_ratio: float | None
    max_intensity: float | None
    peak_separation_sec: float | None = None
    review_reasons: tuple[str, ...] = ()


def resolution_passes_threshold(
    resolution: float | None,
    min_resolution: float | None,
) -> bool:
    """Return whether resolution passes; zero/None explicitly disable the threshold."""

    return (
        min_resolution is None
        or min_resolution <= 0
        or (resolution is not None and resolution >= min_resolution)
    )


def extract_target_eics(raw_path, targets) -> tuple[list[float], dict[float, array]]:
    """Parse a raw file once and extract max-intensity MS1 EICs for all targets."""

    return extract_target_eics_with_config(raw_path, targets, PeakPickingConfig())


def extract_target_eics_with_config(
    raw_path,
    targets,
    config: PeakPickingConfig,
) -> tuple[list[float], dict[float, array]]:
    """Configured form of :func:`extract_target_eics`."""

    unique_targets = sorted({float(target) for target in targets if number(target) is not None})
    traces = {target: array("f") for target in unique_targets}
    rts: list[float] = []
    for rt, mz_values, intensity_values in iter_ms1_spectra(Path(raw_path)):
        rts.append(rt)
        for target in unique_targets:
            delta = target * config.eic_ppm / 1_000_000
            first = bisect_left(mz_values, target - delta)
            last = bisect_right(mz_values, target + delta)
            traces[target].append(max(intensity_values[first:last], default=0.0))
    return rts, traces


def iter_ms1_spectra(path: Path):
    """Yield ``(rt_seconds, sorted_mz, intensity)`` from mzML or mzXML."""

    suffix = path.suffix.lower()
    if suffix == ".mzxml":
        yield from _iter_mzxml_ms1(path)
    elif suffix == ".mzml":
        yield from _iter_mzml_ms1(path)
    else:
        raise ValueError(f"Unsupported raw file type: {path.name}")


def _iter_mzxml_ms1(path: Path):
    for _event, element in ET.iterparse(path, events=("end",)):
        if _local_name(element.tag) != "scan":
            continue
        if element.attrib.get("msLevel") == "1":
            peaks_element = next(
                (child for child in element if _local_name(child.tag) == "peaks"),
                None,
            )
            rt = _parse_mzxml_duration(element.attrib.get("retentionTime", ""))
            if peaks_element is not None and rt is not None:
                mz_values, intensity_values = _decode_mzxml_peaks(peaks_element)
                yield rt, mz_values, intensity_values
        element.clear()


def _iter_mzml_ms1(path: Path):
    for _event, element in ET.iterparse(path, events=("end",)):
        if _local_name(element.tag) != "spectrum":
            continue
        ms_level = _mzml_cv_value(element, "MS:1000511")
        if ms_level == "1":
            rt = _mzml_scan_time(element)
            arrays = _decode_mzml_arrays(element)
            if rt is not None and "mz" in arrays and "intensity" in arrays:
                yield rt, arrays["mz"], arrays["intensity"]
        element.clear()

def pick_chiral_peaks(
    rts: list[float],
    intensities,
    config: PeakPickingConfig | None = None,
) -> PeakPickResult:
    """Smooth one EIC, pick/refine up to two peaks, and compute quality metrics."""

    config = config or PeakPickingConfig()
    raw = [float(value) for value in intensities]
    if len(raw) != len(rts) or len(raw) < config.min_scan_points:
        return PeakPickResult("no_signal", "not_evaluable", (), None, None, max(raw, default=None))
    max_intensity = max(raw, default=0.0)
    if config.min_height is not None and max_intensity < config.min_height:
        return PeakPickResult("no_signal", "not_evaluable", (), None, None, max(raw, default=None))

    window = min(config.sg_window, len(raw) if len(raw) % 2 else len(raw) - 1)
    if window < 3:
        return PeakPickResult("no_signal", "not_evaluable", (), None, None, max(raw, default=None))
    polyorder = min(config.sg_polyorder, window - 1)
    smooth = savgol_smooth(raw, window, polyorder)
    baseline = statistics.median(raw)
    adaptive_candidates: dict[int, _PeakCandidate] = {}
    if config.adaptive_detection:
        candidate_details = _find_adaptive_peaks(smooth, rts, baseline, config)
        candidates = [candidate.scan_index for candidate in candidate_details]
        adaptive_candidates = {
            candidate.scan_index: candidate for candidate in candidate_details
        }
    else:
        # Keep the historical absolute-height/scan-distance implementation
        # untouched for existing configurations and calibrated profiles.
        assert config.min_height is not None
        candidates = _find_peaks(smooth, config.min_height, config.min_distance_scans)
    if not candidates:
        if config.adaptive_detection and max_intensity <= baseline:
            return PeakPickResult("no_signal", "not_evaluable", (), None, None, max_intensity)
        reasons = ["no_qualified_peak_candidate"]
        if config.min_height is None:
            reasons.append("absolute_height_disabled")
        return PeakPickResult(
            "detected",
            "ambiguous",
            (),
            None,
            None,
            max(raw),
            review_reasons=tuple(reasons),
        )

    if config.adaptive_detection:
        refined_pairs = _refine_adaptive_candidates(
            raw,
            smooth,
            candidates,
            config.refine_radius_scans,
        )
        if config.rank_by_prominence:
            refined_pairs.sort(
                key=lambda pair: adaptive_candidates[pair[0]].prominence,
                reverse=True,
            )
        else:
            refined_pairs.sort(key=lambda pair: raw[pair[1]], reverse=True)
        refined_pairs = sorted(refined_pairs[: config.top_k], key=lambda pair: pair[1])
    else:
        refined: list[int] = []
        for candidate in candidates:
            first = max(0, candidate - config.refine_radius_scans)
            last = min(len(raw), candidate + config.refine_radius_scans + 1)
            apex = max(range(first, last), key=raw.__getitem__)
            if apex not in refined:
                refined.append(apex)
        refined = sorted(refined, key=raw.__getitem__, reverse=True)[: config.top_k]
        refined.sort()
        refined_pairs = [(index, index) for index in refined]

    mad = statistics.median(abs(value - baseline) for value in raw)
    noise = 1.4826 * mad
    built_peaks: list[PickedPeak] = []
    for candidate_index, index in refined_pairs:
        bounds = _fwhm_bounds(raw, index, baseline)
        width = None if bounds is None else rts[bounds[1]] - rts[bounds[0]]
        area = None if bounds is None else _baseline_corrected_area(rts, raw, bounds, baseline)
        candidate_detail = adaptive_candidates.get(candidate_index)
        built_peaks.append(
            PickedPeak(
                scan_index=index,
                rt_sec=rts[index],
                intensity=raw[index],
                width_sec=width,
                area_fwhm=area,
                snr=None if noise == 0 else max(0.0, (raw[index] - baseline) / noise),
                prominence=None if candidate_detail is None else candidate_detail.prominence,
                prominence_ratio=(
                    None if candidate_detail is None else candidate_detail.prominence_ratio
                ),
            )
        )
    peaks = tuple(built_peaks)
    if len(peaks) < 2:
        reasons = ("absolute_height_disabled",) if config.min_height is None else ()
        return PeakPickResult(
            "detected",
            "single_peak",
            peaks,
            None,
            None,
            max(raw),
            review_reasons=reasons,
        )

    first_peak, second_peak = peaks[0], peaks[1]
    valley = min(raw[first_peak.scan_index : second_peak.scan_index + 1])
    valley_ratio = valley / min(first_peak.intensity, second_peak.intensity)
    resolution = None
    if first_peak.width_sec and second_peak.width_sec:
        resolution = 2 * (second_peak.rt_sec - first_peak.rt_sec) / (
            first_peak.width_sec + second_peak.width_sec
        )
    second_peak_ratio = min(first_peak.intensity, second_peak.intensity) / max(
        first_peak.intensity,
        second_peak.intensity,
    )
    passes_quality = (
        (config.max_valley_ratio is None or valley_ratio <= config.max_valley_ratio)
        and resolution_passes_threshold(resolution, config.min_resolution)
        and second_peak_ratio >= config.min_second_peak_ratio
    )
    status = (
        "double_peak"
        if first_peak.rt_sec < second_peak.rt_sec and passes_quality
        else "ambiguous"
    )
    separation_sec = second_peak.rt_sec - first_peak.rt_sec
    review_reasons: list[str] = []
    if config.min_height is None:
        review_reasons.append("absolute_height_disabled")
    close_for_review = (
        separation_sec <= 2 * config.min_distance_sec
        if config.min_distance_sec is not None
        else abs(second_peak.scan_index - first_peak.scan_index) < config.min_distance_scans
    )
    if close_for_review:
        review_reasons.append("close_peak_pair")
    if config.max_valley_ratio is not None and valley_ratio > config.max_valley_ratio:
        review_reasons.append("shallow_valley")
    if not resolution_passes_threshold(resolution, config.min_resolution):
        review_reasons.append("resolution_below_threshold")
    if second_peak_ratio < config.min_second_peak_ratio:
        review_reasons.append("weak_second_peak")
    return PeakPickResult(
        "detected",
        status,
        peaks,
        resolution,
        valley_ratio,
        max(raw),
        peak_separation_sec=separation_sec,
        review_reasons=tuple(review_reasons),
    )


def savgol_smooth(values: list[float], window: int, polyorder: int) -> list[float]:
    """Dependency-free Savitzky-Golay smoothing at derivative order zero."""

    weights = _savgol_weights(window, polyorder)
    radius = window // 2
    smoothed = list(values)
    for center in range(radius, len(values) - radius):
        smoothed[center] = sum(
            weight * values[center + offset]
            for weight, offset in zip(weights, range(-radius, radius + 1))
        )
    return smoothed


@lru_cache(maxsize=32)
def _savgol_weights(window: int, polyorder: int) -> tuple[float, ...]:
    radius = window // 2
    design = [[float(offset**power) for power in range(polyorder + 1)] for offset in range(-radius, radius + 1)]
    gram = [
        [sum(row[i] * row[j] for row in design) for j in range(polyorder + 1)]
        for i in range(polyorder + 1)
    ]
    inverse = _invert_matrix(gram)
    return tuple(
        sum(inverse[0][power] * row[power] for power in range(polyorder + 1))
        for row in design
    )


def _invert_matrix(matrix: list[list[float]]) -> list[list[float]]:
    size = len(matrix)
    augmented = [row[:] + [float(i == j) for j in range(size)] for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) < 1e-15:
            raise ValueError("Singular Savitzky-Golay design matrix.")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                left - factor * right
                for left, right in zip(augmented[row], augmented[column])
            ]
    return [row[size:] for row in augmented]


def _find_peaks(values: list[float], min_height: float, min_distance: int) -> list[int]:
    candidates = [
        index
        for index in range(1, len(values) - 1)
        if values[index] >= min_height
        and values[index] > values[index - 1]
        and values[index] >= values[index + 1]
    ]
    accepted: list[int] = []
    for index in sorted(candidates, key=values.__getitem__, reverse=True):
        if all(abs(index - other) >= min_distance for other in accepted):
            accepted.append(index)
    return sorted(accepted)


def _find_adaptive_peaks(
    values: list[float],
    rts: list[float],
    baseline: float,
    config: PeakPickingConfig,
) -> list[_PeakCandidate]:
    """Find locally prominent peaks without requiring a universal intensity floor.

    Prominence is measured above the higher of the local contour base and the
    global median baseline.  Dividing it by the baseline-corrected peak height
    makes a clean low peak comparable to a clean high peak.  Continuous support
    is the number of adjacent scans above half prominence.
    """

    local_maxima = [
        index
        for index in range(1, len(values) - 1)
        if values[index] > values[index - 1]
        and values[index] >= values[index + 1]
        and (config.min_height is None or values[index] >= config.min_height)
    ]
    candidates: list[_PeakCandidate] = []
    for index in local_maxima:
        prominence, ratio = _baseline_corrected_prominence(values, index, baseline)
        support = _half_prominence_support(values, index, prominence)
        if config.min_prominence_ratio is not None and ratio < config.min_prominence_ratio:
            continue
        if config.min_support_scans is not None and support < config.min_support_scans:
            continue
        candidates.append(_PeakCandidate(index, prominence, ratio, support))

    rank = (
        (lambda candidate: candidate.prominence)
        if config.rank_by_prominence
        else (lambda candidate: values[candidate.scan_index])
    )
    accepted: list[_PeakCandidate] = []
    for candidate in sorted(candidates, key=rank, reverse=True):
        conflicts = [
            other
            for other in accepted
            if _peaks_are_too_close(candidate.scan_index, other.scan_index, rts, config)
        ]
        if not conflicts or (
            config.close_peak_max_valley_ratio is not None
            and all(
                _baseline_corrected_valley_ratio(
                    values,
                    candidate.scan_index,
                    other.scan_index,
                    baseline,
                )
                <= config.close_peak_max_valley_ratio
                for other in conflicts
            )
        ):
            accepted.append(candidate)
    return sorted(accepted, key=lambda candidate: candidate.scan_index)


def _baseline_corrected_prominence(
    values: list[float],
    index: int,
    baseline: float,
) -> tuple[float, float]:
    peak_height = values[index]
    left_min = peak_height
    for cursor in range(index - 1, -1, -1):
        left_min = min(left_min, values[cursor])
        if values[cursor] > peak_height:
            break
    right_min = peak_height
    for cursor in range(index + 1, len(values)):
        right_min = min(right_min, values[cursor])
        if values[cursor] > peak_height:
            break
    contour_base = max(baseline, left_min, right_min)
    prominence = max(0.0, peak_height - contour_base)
    corrected_height = peak_height - baseline
    ratio = 0.0 if corrected_height <= 0 else min(1.0, prominence / corrected_height)
    return prominence, ratio


def _half_prominence_support(values: list[float], index: int, prominence: float) -> int:
    if prominence <= 0:
        return 0
    threshold = values[index] - prominence / 2
    left = index
    while left > 0 and values[left - 1] >= threshold:
        left -= 1
    right = index
    while right < len(values) - 1 and values[right + 1] >= threshold:
        right += 1
    return right - left + 1


def _peaks_are_too_close(
    first: int,
    second: int,
    rts: list[float],
    config: PeakPickingConfig,
) -> bool:
    if config.min_distance_sec is not None:
        return abs(rts[first] - rts[second]) < config.min_distance_sec
    return abs(first - second) < config.min_distance_scans


def _baseline_corrected_valley_ratio(
    values: list[float],
    first: int,
    second: int,
    baseline: float,
) -> float:
    left, right = sorted((first, second))
    valley = min(values[left : right + 1])
    lower_peak_height = min(values[first], values[second]) - baseline
    if lower_peak_height <= 0:
        return math.inf
    return max(0.0, valley - baseline) / lower_peak_height


def _refine_adaptive_candidates(
    raw: list[float],
    smooth: list[float],
    candidates: list[int],
    radius: int,
) -> list[tuple[int, int]]:
    """Refine candidates on raw data while keeping adjacent peak basins separate."""

    ordered = sorted(candidates)
    valley_boundaries = [
        min(range(left + 1, right), key=smooth.__getitem__)
        for left, right in zip(ordered, ordered[1:])
    ]
    refined: list[tuple[int, int]] = []
    for position, candidate in enumerate(ordered):
        first = max(0, candidate - radius)
        last = min(len(raw), candidate + radius + 1)
        if position > 0:
            first = max(first, valley_boundaries[position - 1] + 1)
        if position < len(valley_boundaries):
            last = min(last, valley_boundaries[position] + 1)
        if first >= last:
            first, last = candidate, candidate + 1
        refined.append((candidate, max(range(first, last), key=raw.__getitem__)))
    return refined


def _fwhm_bounds(values: list[float], apex: int, baseline: float) -> tuple[int, int] | None:
    half_height = baseline + (values[apex] - baseline) / 2
    left = apex
    while left > 0 and values[left] >= half_height:
        left -= 1
    right = apex
    while right < len(values) - 1 and values[right] >= half_height:
        right += 1
    if left == apex or right == apex or right <= left:
        return None
    return left, right


def _baseline_corrected_area(
    rts: list[float],
    values: list[float],
    bounds: tuple[int, int],
    baseline: float,
) -> float:
    left, right = bounds
    corrected = [max(0.0, value - baseline) for value in values[left : right + 1]]
    return sum(
        (rts[left + index + 1] - rts[left + index])
        * (corrected[index] + corrected[index + 1])
        / 2
        for index in range(len(corrected) - 1)
    )


def annotate_peaklist(
    peaklist_path,
    raw_path,
    output_path,
    *,
    mz_column: str = "MZ",
    config: PeakPickingConfig | None = None,
    review_output_path=None,
):
    """Add automatic MS1 peak and quality columns to a CSV/XLSX peaklist."""

    config = config or PeakPickingConfig()
    rows = read_table(Path(peaklist_path))
    if not rows or mz_column not in rows[0]:
        raise ValueError(f"Peaklist must contain an exact {mz_column!r} column.")
    targets = [value for row in rows if (value := number(row.get(mz_column))) is not None]
    rts, traces = extract_target_eics_with_config(raw_path, targets, config)
    for row in rows:
        mz = number(row.get(mz_column))
        result = (
            pick_chiral_peaks(rts, traces[mz], config)
            if mz is not None and mz in traces
            else PeakPickResult("no_signal", "not_evaluable", (), None, None, None)
        )
        row["auto_signal_status"] = result.signal_status
        row["auto_chromatographic_status"] = result.chromatographic_status
        row["auto_peak1_rt_min"] = csv_scalar(result.peaks[0].rt_sec / 60) if result.peaks else ""
        row["auto_peak2_rt_min"] = csv_scalar(result.peaks[1].rt_sec / 60) if len(result.peaks) > 1 else ""
        row["auto_peak1_intensity"] = csv_scalar(result.peaks[0].intensity) if result.peaks else ""
        row["auto_peak2_intensity"] = csv_scalar(result.peaks[1].intensity) if len(result.peaks) > 1 else ""
        row["auto_peak1_width_sec"] = csv_scalar(result.peaks[0].width_sec) if result.peaks else ""
        row["auto_peak2_width_sec"] = csv_scalar(result.peaks[1].width_sec) if len(result.peaks) > 1 else ""
        row["auto_peak1_area_fwhm"] = csv_scalar(result.peaks[0].area_fwhm) if result.peaks else ""
        row["auto_peak2_area_fwhm"] = csv_scalar(result.peaks[1].area_fwhm) if len(result.peaks) > 1 else ""
        row["auto_peak1_snr"] = csv_scalar(result.peaks[0].snr) if result.peaks else ""
        row["auto_peak2_snr"] = csv_scalar(result.peaks[1].snr) if len(result.peaks) > 1 else ""
        row["auto_peak1_prominence"] = csv_scalar(result.peaks[0].prominence) if result.peaks else ""
        row["auto_peak2_prominence"] = csv_scalar(result.peaks[1].prominence) if len(result.peaks) > 1 else ""
        row["auto_peak1_prominence_ratio"] = (
            csv_scalar(result.peaks[0].prominence_ratio) if result.peaks else ""
        )
        row["auto_peak2_prominence_ratio"] = (
            csv_scalar(result.peaks[1].prominence_ratio) if len(result.peaks) > 1 else ""
        )
        total_area = sum(peak.area_fwhm or 0.0 for peak in result.peaks[:2])
        row["auto_peak1_area_fraction"] = (
            csv_scalar((result.peaks[0].area_fwhm or 0.0) / total_area)
            if result.peaks and total_area > 0
            else ""
        )
        row["auto_peak_resolution"] = csv_scalar(result.resolution)
        row["auto_peak_separation_sec"] = csv_scalar(result.peak_separation_sec)
        row["auto_valley_ratio"] = csv_scalar(result.valley_ratio)
        row["auto_max_intensity"] = csv_scalar(result.max_intensity)
        row["auto_manual_review_required"] = "true" if result.review_reasons else "false"
        row["auto_review_reasons"] = ";".join(result.review_reasons)
    write_table(Path(output_path), rows)
    if review_output_path is not None:
        write_table(
            Path(review_output_path),
            [row for row in rows if row["auto_manual_review_required"] == "true"],
        )
    return rows


def load_peak_config(path) -> PeakPickingConfig:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    parameters = dict(payload.get("peak_picking_config", payload))
    parameters.update(payload.get("quality_thresholds", {}))
    return PeakPickingConfig(**parameters)


def _main(argv=None):
    parser = argparse.ArgumentParser(description="Automatically pick up to two chiral peaks from targeted MS1 EICs.")
    parser.add_argument("--peaklist", required=True)
    parser.add_argument("--raw", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--review-output",
        help="Optional CSV containing only rows flagged for manual review.",
    )
    parser.add_argument("--mz-column", default="MZ")
    parser.add_argument("--config-json")
    parser.add_argument("--eic-ppm", type=float, default=3.0)
    parser.add_argument("--sg-window", type=int, default=7)
    parser.add_argument("--sg-polyorder", type=int, default=2)
    parser.add_argument("--min-height", type=float, default=100000.0)
    parser.add_argument("--min-distance-scans", type=int, default=20)
    parser.add_argument("--refine-radius-scans", type=int, default=3)
    parser.add_argument("--min-prominence-ratio", type=float)
    parser.add_argument("--min-support-scans", type=int)
    parser.add_argument("--min-distance-sec", type=float)
    parser.add_argument("--close-peak-max-valley-ratio", type=float)
    parser.add_argument("--rank-by-prominence", action="store_true")
    parser.add_argument("--max-valley-ratio", type=float)
    parser.add_argument("--min-resolution", type=float, default=0.0)
    parser.add_argument("--min-second-peak-ratio", type=float, default=0.0)
    args = parser.parse_args(argv)
    config = (
        load_peak_config(args.config_json)
        if args.config_json
        else PeakPickingConfig(
            eic_ppm=args.eic_ppm,
            sg_window=args.sg_window,
            sg_polyorder=args.sg_polyorder,
            min_height=args.min_height,
            min_distance_scans=args.min_distance_scans,
            refine_radius_scans=args.refine_radius_scans,
            min_prominence_ratio=args.min_prominence_ratio,
            min_support_scans=args.min_support_scans,
            min_distance_sec=args.min_distance_sec,
            close_peak_max_valley_ratio=args.close_peak_max_valley_ratio,
            rank_by_prominence=args.rank_by_prominence,
            max_valley_ratio=args.max_valley_ratio,
            min_resolution=args.min_resolution,
            min_second_peak_ratio=args.min_second_peak_ratio,
        )
    )
    rows = annotate_peaklist(
        args.peaklist,
        args.raw,
        args.output,
        mz_column=args.mz_column,
        config=config,
        review_output_path=args.review_output,
    )
    counts: dict[str, int] = {}
    for row in rows:
        status = row["auto_chromatographic_status"]
        counts[status] = counts.get(status, 0) + 1
    print("; ".join(f"{key}={counts[key]}" for key in sorted(counts)))


if __name__ == "__main__":
    _main()
