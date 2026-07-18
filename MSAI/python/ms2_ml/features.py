"""Fixed, leakage-resistant spectrum features for the shadow classifier."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from typing import Any

try:
    from ..ms2.similarity import compare_fragment_spectra
except ImportError:  # pragma: no cover - direct PYTHONPATH=MSAI/python execution
    from ms2.similarity import compare_fragment_spectra  # type: ignore[import-not-found]
from .contracts import validate_no_leakage

FEATURE_NAMES = (
    "cosine",
    "entropy_similarity",
    "log1p_matched_peaks",
    "min_explained_intensity",
    "explained_intensity_asymmetry",
    "fragment_count_balance",
)
validate_no_leakage(FEATURE_NAMES)


@dataclass(frozen=True)
class SpectrumFeatures:
    """Six interpretable measurements computed directly from two spectra."""

    cosine: float
    entropy_similarity: float
    log1p_matched_peaks: float
    min_explained_intensity: float
    explained_intensity_asymmetry: float
    fragment_count_balance: float

    def as_dict(self) -> dict[str, float]:
        """Return a JSON-safe feature mapping in the frozen feature order."""

        return asdict(self)

    def as_vector(self) -> tuple[float, ...]:
        """Return values in ``FEATURE_NAMES`` order."""

        return tuple(getattr(self, name) for name in FEATURE_NAMES)


def extract_features(spectrum_a: Any, spectrum_b: Any) -> SpectrumFeatures | None:
    """Compute fixed features from two spectra using frozen preprocessing.

    The function exposes no tolerance/filter arguments.  Alignment is always
    exclusive one-to-one at 0.01 Da after the 1% relative-intensity filter,
    matching the frozen extraction/diagnostic contract.  Empty or otherwise
    non-evaluable bilateral input returns ``None`` rather than a negative row.
    """

    similarity = compare_fragment_spectra(
        spectrum_a,
        spectrum_b,
        fragment_mz_tol=0.01,
        fragment_mz_tol_unit="Da",
        min_relative_intensity=0.01,
        intensity_power=0.5,
    )
    if (
        similarity.cosine is None
        or similarity.entropy_similarity is None
        or similarity.peak_a_explained_intensity is None
        or similarity.peak_b_explained_intensity is None
        or similarity.peak_a_count <= 0
        or similarity.peak_b_count <= 0
    ):
        return None
    explained_a = similarity.peak_a_explained_intensity
    explained_b = similarity.peak_b_explained_intensity
    values = SpectrumFeatures(
        cosine=similarity.cosine,
        entropy_similarity=similarity.entropy_similarity,
        log1p_matched_peaks=math.log1p(similarity.matched_peaks),
        min_explained_intensity=min(explained_a, explained_b),
        explained_intensity_asymmetry=abs(explained_a - explained_b),
        fragment_count_balance=(
            min(similarity.peak_a_count, similarity.peak_b_count)
            / max(similarity.peak_a_count, similarity.peak_b_count)
        ),
    )
    if not all(math.isfinite(value) for value in values.as_vector()):
        return None
    return values


def extract_features_from_row(
    row: Mapping[str, Any],
    spectrum_a_column: str = "peak_a_MS2",
    spectrum_b_column: str = "peak_b_MS2",
) -> SpectrumFeatures | None:
    """Read only the two named spectrum cells and compute the fixed features."""

    return extract_features(row.get(spectrum_a_column), row.get(spectrum_b_column))


def feature_vector(
    features: SpectrumFeatures | Mapping[str, Any] | Sequence[float],
) -> tuple[float, ...]:
    """Normalize an exact six-feature record to a finite ordered tuple."""

    if isinstance(features, SpectrumFeatures):
        vector = features.as_vector()
    elif isinstance(features, Mapping):
        missing = [name for name in FEATURE_NAMES if name not in features]
        if missing:
            raise ValueError(f"Missing fixed features: {', '.join(missing)}.")
        vector = tuple(float(features[name]) for name in FEATURE_NAMES)
    else:
        vector = tuple(float(value) for value in features)
    if len(vector) != len(FEATURE_NAMES):
        raise ValueError(f"Expected {len(FEATURE_NAMES)} features, got {len(vector)}.")
    if not all(math.isfinite(value) for value in vector):
        raise ValueError("All spectrum features must be finite.")
    return vector


__all__ = [
    "FEATURE_NAMES",
    "SpectrumFeatures",
    "extract_features",
    "extract_features_from_row",
    "feature_vector",
]
