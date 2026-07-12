"""MS2 spectrum comparison and evidence-based chiral-pair classification.

The two peaks of an enantiomer pair should have the same precursor mass and,
under an achiral ion source/collision cell, highly similar product-ion spectra.
That evidence is useful for screening peaks already separated by a chiral LC
method, but it cannot assign absolute R/S configuration by itself.
"""

from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from dataclasses import dataclass


@dataclass(frozen=True)
class SpectrumSimilarity:
    """Summary of a one-to-one fragment-ion comparison."""

    cosine: float | None
    matched_peaks: int
    peak_a_count: int
    peak_b_count: int
    peak_a_explained_intensity: float | None
    peak_b_explained_intensity: float | None
    entropy_similarity: float | None = None


@dataclass(frozen=True)
class ChiralPairThresholds:
    """Conservative default screen, configurable for method validation.

    A cosine threshold of 0.7 and six matched ions follow the commonly used
    GNPS spectral-similarity screen.  Explained-intensity guards prevent a
    handful of low-abundance coincidental matches from passing.
    """

    min_cosine: float = 0.7
    min_matched_peaks: int = 6
    min_explained_intensity: float = 0.5
    min_entropy_similarity: float | None = None

    def __post_init__(self):
        if not 0 <= self.min_cosine <= 1:
            raise ValueError("min_cosine must be between 0 and 1.")
        if self.min_matched_peaks < 1:
            raise ValueError("min_matched_peaks must be at least 1.")
        if not 0 <= self.min_explained_intensity <= 1:
            raise ValueError("min_explained_intensity must be between 0 and 1.")
        if self.min_entropy_similarity is not None and not 0 <= self.min_entropy_similarity <= 1:
            raise ValueError("min_entropy_similarity must be between 0 and 1 when provided.")


def parse_fragment_string(value: str | None) -> list[tuple[float, float]]:
    """Parse the project CSV format ``mz,intensity;mz,intensity``.

    Malformed, non-finite, or non-positive peaks are discarded so that one bad
    cell does not abort a batch analysis.
    """

    peaks: list[tuple[float, float]] = []
    for token in str(value or "").split(";"):
        fields = token.split(",")
        if len(fields) != 2:
            continue
        try:
            mz, intensity = (float(field.strip()) for field in fields)
        except ValueError:
            continue
        if math.isfinite(mz) and math.isfinite(intensity) and mz > 0 and intensity > 0:
            peaks.append((mz, intensity))
    return sorted(peaks)


def compare_fragment_spectra(
    peak_a,
    peak_b,
    fragment_mz_tol: float = 0.01,
    fragment_mz_tol_unit: str = "Da",
    min_relative_intensity: float = 0.01,
    intensity_power: float = 0.5,
) -> SpectrumSimilarity:
    """Compare two fragment spectra with exclusive peak matching.

    Intensities are square-root transformed by default before cosine scoring,
    which reduces domination by a single base peak.  A peak may be matched at
    most once; candidate pairs are selected by mass error, then by intensity
    product.  Unmatched peaks remain in each cosine norm.
    """

    if fragment_mz_tol < 0 or not math.isfinite(fragment_mz_tol):
        raise ValueError("fragment_mz_tol must be finite and non-negative.")
    if fragment_mz_tol_unit not in {"Da", "ppm"}:
        raise ValueError("fragment_mz_tol_unit must be either 'Da' or 'ppm'.")
    if not 0 <= min_relative_intensity <= 1:
        raise ValueError("min_relative_intensity must be between 0 and 1.")
    if intensity_power <= 0 or not math.isfinite(intensity_power):
        raise ValueError("intensity_power must be finite and positive.")

    peaks_a = _prepare_peaks(peak_a, min_relative_intensity)
    peaks_b = _prepare_peaks(peak_b, min_relative_intensity)
    if not peaks_a or not peaks_b:
        return SpectrumSimilarity(None, 0, len(peaks_a), len(peaks_b), None, None, None)

    matches = _exclusive_matches(peaks_a, peaks_b, fragment_mz_tol, fragment_mz_tol_unit)

    weighted_a = [intensity**intensity_power for _mz, intensity in peaks_a]
    weighted_b = [intensity**intensity_power for _mz, intensity in peaks_b]
    norm_a = math.sqrt(sum(value * value for value in weighted_a))
    norm_b = math.sqrt(sum(value * value for value in weighted_b))
    cosine = sum(weighted_a[i] * weighted_b[j] for i, j in matches) / (norm_a * norm_b)

    total_a = sum(intensity for _mz, intensity in peaks_a)
    total_b = sum(intensity for _mz, intensity in peaks_b)
    explained_a = sum(peaks_a[i][1] for i, _j in matches) / total_a
    explained_b = sum(peaks_b[j][1] for _i, j in matches) / total_b
    return SpectrumSimilarity(
        cosine=min(1.0, max(0.0, cosine)),
        matched_peaks=len(matches),
        peak_a_count=len(peaks_a),
        peak_b_count=len(peaks_b),
        peak_a_explained_intensity=explained_a,
        peak_b_explained_intensity=explained_b,
        entropy_similarity=_entropy_similarity(peaks_a, peaks_b, matches),
    )


def align_fragment_spectra(
    peak_a,
    peak_b,
    fragment_mz_tol: float = 0.01,
    fragment_mz_tol_unit: str = "Da",
    min_relative_intensity: float = 0.01,
) -> tuple[list[tuple[float, float]], list[tuple[float, float]], list[tuple[int, int]]]:
    """Return the exact filtered one-to-one alignment used for visualization.

    This keeps mirror-spectrum match coloring consistent with the similarity
    calculation instead of independently matching the unfiltered CSV strings.
    """

    if fragment_mz_tol < 0 or not math.isfinite(fragment_mz_tol):
        raise ValueError("fragment_mz_tol must be finite and non-negative.")
    if fragment_mz_tol_unit not in {"Da", "ppm"}:
        raise ValueError("fragment_mz_tol_unit must be either 'Da' or 'ppm'.")
    if not 0 <= min_relative_intensity <= 1:
        raise ValueError("min_relative_intensity must be between 0 and 1.")
    peaks_a = _prepare_peaks(peak_a, min_relative_intensity)
    peaks_b = _prepare_peaks(peak_b, min_relative_intensity)
    matches = _exclusive_matches(peaks_a, peaks_b, fragment_mz_tol, fragment_mz_tol_unit)
    return peaks_a, peaks_b, matches


def classify_chiral_pair(
    similarity: SpectrumSimilarity,
    thresholds: ChiralPairThresholds | None = None,
) -> tuple[str, str]:
    """Classify a separated peak pair without over-claiming R/S identity."""

    thresholds = thresholds or ChiralPairThresholds()
    if similarity.cosine is None:
        return "not_evaluable", "one_or_both_ms2_spectra_are_empty"
    if similarity.matched_peaks < thresholds.min_matched_peaks:
        return (
            "insufficient_ms2_evidence",
            f"matched_peaks<{thresholds.min_matched_peaks}",
        )
    explained_a = similarity.peak_a_explained_intensity or 0.0
    explained_b = similarity.peak_b_explained_intensity or 0.0
    if similarity.cosine < thresholds.min_cosine:
        return "ms2_not_similar", f"cosine<{thresholds.min_cosine:g}"
    if (
        thresholds.min_entropy_similarity is not None
        and (similarity.entropy_similarity or 0.0) < thresholds.min_entropy_similarity
    ):
        return (
            "ms2_not_similar",
            f"entropy_similarity<{thresholds.min_entropy_similarity:g}",
        )
    if min(explained_a, explained_b) < thresholds.min_explained_intensity:
        return (
            "ms2_not_similar",
            f"explained_intensity<{thresholds.min_explained_intensity:g}",
        )
    return (
        "candidate_enantiomer_pair",
        "separated_same_precursor_peaks_with_similar_ms2;R_S_assignment_requires_authentic_enantiomer_standards",
    )


def _prepare_peaks(peak_data, min_relative_intensity: float) -> list[tuple[float, float]]:
    if isinstance(peak_data, str) or peak_data is None:
        peaks = parse_fragment_string(peak_data)
    else:
        peaks = []
        for raw_mz, raw_intensity in peak_data:
            try:
                mz, intensity = float(raw_mz), float(raw_intensity)
            except (TypeError, ValueError):
                continue
            if math.isfinite(mz) and math.isfinite(intensity) and mz > 0 and intensity > 0:
                peaks.append((mz, intensity))
        peaks.sort()
    if not peaks:
        return []
    cutoff = max(intensity for _mz, intensity in peaks) * min_relative_intensity
    return [(mz, intensity) for mz, intensity in peaks if intensity >= cutoff]


def _pair_tolerance(mz_a: float, mz_b: float, tolerance: float, unit: str) -> float:
    if unit == "Da":
        return tolerance
    return ((mz_a + mz_b) / 2) * tolerance / 1_000_000


def _exclusive_matches(
    peaks_a: list[tuple[float, float]],
    peaks_b: list[tuple[float, float]],
    tolerance: float,
    unit: str,
) -> list[tuple[int, int]]:
    """Generate only nearby candidate pairs, then enforce one-to-one use."""

    masses_b = [mz for mz, _intensity in peaks_b]
    candidates: list[tuple[float, float, int, int]] = []
    for i, (mz_a, intensity_a) in enumerate(peaks_a):
        search_tolerance = tolerance if unit == "Da" else mz_a * tolerance / 1_000_000
        first = bisect_left(masses_b, mz_a - search_tolerance * 1.01)
        last = bisect_right(masses_b, mz_a + search_tolerance * 1.01)
        for j in range(first, last):
            mz_b, intensity_b = peaks_b[j]
            error = abs(mz_a - mz_b)
            if error <= _pair_tolerance(mz_a, mz_b, tolerance, unit):
                candidates.append((error, -(intensity_a * intensity_b), i, j))

    used_a: set[int] = set()
    used_b: set[int] = set()
    matches: list[tuple[int, int]] = []
    for _error, _negative_product, i, j in sorted(candidates):
        if i not in used_a and j not in used_b:
            used_a.add(i)
            used_b.add(j)
            matches.append((i, j))
    return matches


def _entropy_similarity(
    peaks_a: list[tuple[float, float]],
    peaks_b: list[tuple[float, float]],
    matches: list[tuple[int, int]],
) -> float:
    """Weighted spectral-entropy similarity on the exclusive alignment."""

    probability_a = _entropy_weight([intensity for _mz, intensity in peaks_a])
    probability_b = _entropy_weight([intensity for _mz, intensity in peaks_b])
    matched_a = {i for i, _j in matches}
    matched_b = {j for _i, j in matches}
    aligned_a = [probability_a[i] for i, _j in matches]
    aligned_b = [probability_b[j] for _i, j in matches]
    aligned_a.extend(probability_a[i] for i in range(len(peaks_a)) if i not in matched_a)
    aligned_b.extend(0.0 for i in range(len(peaks_a)) if i not in matched_a)
    aligned_a.extend(0.0 for j in range(len(peaks_b)) if j not in matched_b)
    aligned_b.extend(probability_b[j] for j in range(len(peaks_b)) if j not in matched_b)

    entropy_a = _shannon_entropy(aligned_a)
    entropy_b = _shannon_entropy(aligned_b)
    mixture = [(left + right) / 2 for left, right in zip(aligned_a, aligned_b)]
    divergence = 2 * _shannon_entropy(mixture) - entropy_a - entropy_b
    return min(1.0, max(0.0, 1.0 - divergence / math.log(4)))


def _entropy_weight(intensities: list[float]) -> list[float]:
    total = sum(intensities)
    probability = [intensity / total for intensity in intensities]
    entropy = _shannon_entropy(probability)
    if entropy < 3:
        power = 0.25 * (1 + entropy)
        probability = [value**power for value in probability]
        weighted_total = sum(probability)
        probability = [value / weighted_total for value in probability]
    return probability


def _shannon_entropy(probability: list[float]) -> float:
    return -sum(value * math.log(value) for value in probability if value > 0)
