"""Interpretable threshold baseline for the MS2 shadow classifier.

Only decision thresholds are tuned here.  Spectrum extraction and acquisition
parameters are deliberately absent from this module and cannot enter the grid.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from .contracts import MIN_SPECIFICITY_TARGET, ML_DIAGNOSTIC_STATUSES
from .features import FEATURE_NAMES
from .metrics import BinaryMetrics, binary_metrics

COSINE_GRID = tuple(round(0.50 + index * 0.025, 3) for index in range(19))
MATCHED_PEAKS_GRID = (2, 3, 4, 5, 6, 8, 10)
EXPLAINED_INTENSITY_GRID = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
ENTROPY_SIMILARITY_GRID: tuple[float | None, ...] = (
    None,
    0.5,
    0.6,
    0.7,
    0.75,
    0.8,
    0.85,
    0.9,
)


@dataclass(frozen=True)
class ThresholdParameters:
    """The four trainable, post-extraction decision thresholds."""

    min_cosine: float
    min_matched_peaks: int
    min_explained_intensity: float
    min_entropy_similarity: float | None

    def to_dict(self) -> dict[str, float | int | None]:
        return {
            "min_cosine": self.min_cosine,
            "min_matched_peaks": self.min_matched_peaks,
            "min_explained_intensity": self.min_explained_intensity,
            "min_entropy_similarity": self.min_entropy_similarity,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> ThresholdParameters:
        expected_keys = {
            "min_cosine",
            "min_matched_peaks",
            "min_explained_intensity",
            "min_entropy_similarity",
        }
        if set(value) != expected_keys:
            raise ValueError("Threshold parameters must use the exact four-field schema.")
        matched_peaks = value["min_matched_peaks"]
        if isinstance(matched_peaks, bool) or not isinstance(matched_peaks, int):
            raise ValueError("min_matched_peaks must be an integer from the fixed grid.")
        entropy = value["min_entropy_similarity"]
        if isinstance(entropy, bool):
            raise ValueError("min_entropy_similarity must be null or a fixed-grid number.")
        parameters = cls(
            min_cosine=_finite_float(value["min_cosine"], "min_cosine"),
            min_matched_peaks=matched_peaks,
            min_explained_intensity=_finite_float(
                value["min_explained_intensity"], "min_explained_intensity"
            ),
            min_entropy_similarity=(
                None if entropy is None else _finite_float(entropy, "min_entropy_similarity")
            ),
        )
        validate_threshold_parameters(parameters)
        return parameters


@dataclass(frozen=True)
class ThresholdSearchResult:
    """Selected threshold rule and its development-set diagnostics."""

    parameters: ThresholdParameters
    metrics: BinaryMetrics
    candidates_evaluated: int
    specificity_target: float
    bucket_probabilities: Mapping[str, float]

    def to_dict(self) -> dict[str, object]:
        return {
            "parameters": self.parameters.to_dict(),
            "metrics": self.metrics.as_dict(),
            "candidates_evaluated": self.candidates_evaluated,
            "specificity_target": self.specificity_target,
            "constraint_met": self.metrics.specificity >= self.specificity_target,
            "bucket_probabilities": dict(self.bucket_probabilities),
        }


def iter_threshold_grid() -> Iterable[ThresholdParameters]:
    """Yield the frozen 6,384-combination decision grid in stable order."""

    for cosine, matched, explained, entropy in itertools.product(
        COSINE_GRID,
        MATCHED_PEAKS_GRID,
        EXPLAINED_INTENSITY_GRID,
        ENTROPY_SIMILARITY_GRID,
    ):
        yield ThresholdParameters(cosine, matched, explained, entropy)


def validate_threshold_parameters(parameters: ThresholdParameters) -> ThresholdParameters:
    """Require every serialized threshold to come from the prespecified grid."""

    if (
        parameters.min_cosine not in COSINE_GRID
        or parameters.min_matched_peaks not in MATCHED_PEAKS_GRID
        or parameters.min_explained_intensity not in EXPLAINED_INTENSITY_GRID
        or parameters.min_entropy_similarity not in ENTROPY_SIMILARITY_GRID
    ):
        raise ValueError("Threshold model parameters must come from the fixed 6,384-point grid.")
    return parameters


def threshold_bucket(
    features: Sequence[float] | Mapping[str, float],
    parameters: ThresholdParameters,
) -> str:
    """Map an evaluable spectrum pair to an interpretable evidence bucket."""

    values = _feature_values(features)
    cosine = values["cosine"]
    entropy = values["entropy_similarity"]
    matched_peaks = round(math.expm1(values["log1p_matched_peaks"]))
    explained = values["min_explained_intensity"]
    if matched_peaks < parameters.min_matched_peaks:
        status = "insufficient_evidence"
    elif (
        cosine < parameters.min_cosine
        or explained < parameters.min_explained_intensity
        or (
            parameters.min_entropy_similarity is not None
            and entropy < parameters.min_entropy_similarity
        )
    ):
        status = "conflicting_spectra"
    else:
        status = "supported_same_compound"
    if status not in ML_DIAGNOSTIC_STATUSES:
        raise AssertionError(f"Unexpected threshold status: {status}")
    return status


def threshold_prediction(
    features: Sequence[float] | Mapping[str, float],
    parameters: ThresholdParameters,
) -> int:
    """Return the binary support prediction used during model selection."""

    return int(threshold_bucket(features, parameters) == "supported_same_compound")


def fit_thresholds(
    examples: Sequence[object],
    *,
    specificity_target: float = 0.95,
) -> ThresholdSearchResult:
    """Maximize recall subject to the predeclared specificity constraint."""

    _validate_specificity_target(specificity_target)
    labels = [_example_label(example) for example in examples]
    if not labels or set(labels) != {0, 1}:
        raise ValueError("Threshold fitting requires both eligible truth classes.")

    best: tuple[tuple[float, ...], ThresholdParameters, BinaryMetrics] | None = None
    count = 0
    for parameters in iter_threshold_grid():
        count += 1
        predictions = [
            threshold_prediction(_example_features(example), parameters) for example in examples
        ]
        metrics = binary_metrics(labels, predictions, threshold=0.5)
        if metrics.specificity < specificity_target:
            continue
        key = _threshold_selection_key(parameters, metrics)
        if best is None or key > best[0]:
            best = key, parameters, metrics
    if best is None:
        raise ValueError(
            "No threshold-grid candidate satisfies the predeclared "
            f"specificity target ({specificity_target:.3f})."
        )

    _key, parameters, metrics = best
    probabilities = estimate_bucket_probabilities(examples, parameters)
    return ThresholdSearchResult(
        parameters=parameters,
        metrics=metrics,
        candidates_evaluated=count,
        specificity_target=specificity_target,
        bucket_probabilities=probabilities,
    )


def evaluate_thresholds(
    examples: Sequence[object],
    parameters: ThresholdParameters,
) -> BinaryMetrics:
    """Evaluate a fixed rule without changing any thresholds."""

    labels = [_example_label(example) for example in examples]
    predictions = [
        threshold_prediction(_example_features(example), parameters) for example in examples
    ]
    return binary_metrics(labels, predictions, threshold=0.5)


def estimate_bucket_probabilities(
    examples: Sequence[object],
    parameters: ThresholdParameters,
) -> dict[str, float]:
    """Estimate transparent Laplace-smoothed probabilities for rule buckets."""

    counts = {
        "supported_same_compound": [0, 0],
        "conflicting_spectra": [0, 0],
        "insufficient_evidence": [0, 0],
    }
    for example in examples:
        bucket = threshold_bucket(_example_features(example), parameters)
        label = _example_label(example)
        counts[bucket][0] += label
        counts[bucket][1] += 1
    return {
        bucket: (positive_count + 1) / (total_count + 2)
        for bucket, (positive_count, total_count) in counts.items()
    }


def select_support_probability_threshold(
    labels: Sequence[int],
    probabilities: Sequence[float],
    *,
    specificity_target: float = 0.95,
) -> tuple[float, BinaryMetrics]:
    """Choose the lowest useful support cutoff that maximizes constrained recall."""

    _validate_probability_inputs(labels, probabilities)
    _validate_specificity_target(specificity_target)
    candidates = sorted({0.0, 1.0, *(float(value) for value in probabilities)})
    eligible: list[tuple[tuple[float, ...], float, BinaryMetrics]] = []
    for threshold in candidates:
        metrics = binary_metrics(labels, probabilities, threshold=threshold)
        if metrics.specificity >= specificity_target:
            eligible.append(
                (
                    (
                        metrics.recall,
                        metrics.precision,
                        metrics.f0_5,
                        metrics.specificity,
                        threshold,
                    ),
                    threshold,
                    metrics,
                )
            )
    if not eligible:
        raise ValueError(
            "No probability cutoff satisfies the predeclared "
            f"specificity target ({specificity_target:.3f})."
        )
    _key, threshold, metrics = max(eligible, key=lambda item: item[0])
    return threshold, metrics


def select_conflict_probability_threshold(
    labels: Sequence[int],
    probabilities: Sequence[float],
    support_threshold: float,
    *,
    positive_safety_target: float = 0.95,
) -> float:
    """Select a low-probability conflict cutoff while protecting positives.

    ``positive_safety_target`` bounds the fraction of true positives that may be
    assigned to the conflict tail.  The middle region is intentionally retained
    as abstention/insufficient evidence.
    """

    _validate_probability_inputs(labels, probabilities)
    _validate_specificity_target(positive_safety_target)
    positive_total = sum(labels)
    negative_total = len(labels) - positive_total
    candidates = sorted(
        {0.0, *(float(value) for value in probabilities if value < support_threshold)}
    )
    best: tuple[tuple[float, float], float] | None = None
    for threshold in candidates:
        positive_conflicts = sum(
            label == 1 and probability <= threshold
            for label, probability in zip(labels, probabilities, strict=True)
        )
        negative_conflicts = sum(
            label == 0 and probability <= threshold
            for label, probability in zip(labels, probabilities, strict=True)
        )
        positive_safety = 1.0 - positive_conflicts / positive_total
        if positive_safety < positive_safety_target:
            continue
        negative_capture = negative_conflicts / negative_total
        candidate = (negative_capture, threshold), threshold
        if best is None or candidate[0] > best[0]:
            best = candidate
    if best is None:
        return 0.0
    return min(best[1], math.nextafter(support_threshold, -math.inf))


def _threshold_selection_key(
    parameters: ThresholdParameters,
    metrics: BinaryMetrics,
) -> tuple[float, ...]:
    entropy = (
        -1.0 if parameters.min_entropy_similarity is None else parameters.min_entropy_similarity
    )
    return (
        metrics.recall,
        metrics.precision,
        metrics.f0_5,
        metrics.specificity,
        parameters.min_cosine,
        float(parameters.min_matched_peaks),
        parameters.min_explained_intensity,
        entropy,
    )


def _feature_values(
    features: Sequence[float] | Mapping[str, float],
) -> dict[str, float]:
    if isinstance(features, Mapping):
        values = {str(name): float(features[name]) for name in FEATURE_NAMES}
    else:
        if len(features) != len(FEATURE_NAMES):
            raise ValueError(f"Expected {len(FEATURE_NAMES)} features, received {len(features)}.")
        values = dict(
            zip(
                (str(name) for name in FEATURE_NAMES),
                (float(value) for value in features),
                strict=True,
            )
        )
    if any(not math.isfinite(value) for value in values.values()):
        raise ValueError("Threshold features must all be finite.")
    return values


def _example_features(example: object) -> Sequence[float] | Mapping[str, float]:
    dynamic: Any = example
    if hasattr(dynamic, "features"):
        return dynamic.features
    if isinstance(example, Mapping):
        value = example.get("features")
        if isinstance(value, Sequence | Mapping):
            return value
    raise TypeError("Examples must expose a fixed 'features' vector.")


def _example_label(example: object) -> int:
    dynamic: Any = example
    if hasattr(dynamic, "label"):
        value = dynamic.label
    elif isinstance(example, Mapping):
        value = example.get("label")
    else:
        raise TypeError("Examples must expose a binary 'label'.")
    label = int(value)  # type: ignore[arg-type]
    if label not in {0, 1}:
        raise ValueError(f"Expected binary label 0/1, received {value!r}.")
    return label


def _validate_probability_inputs(
    labels: Sequence[int],
    probabilities: Sequence[float],
) -> None:
    if len(labels) != len(probabilities) or not labels:
        raise ValueError("Labels and probabilities must have equal non-zero length.")
    if {int(value) for value in labels} != {0, 1}:
        raise ValueError("Probability threshold selection requires both truth classes.")
    if any(not math.isfinite(value) or not 0 <= value <= 1 for value in probabilities):
        raise ValueError("Probabilities must be finite values between zero and one.")


def _validate_specificity_target(value: float) -> None:
    if (
        isinstance(value, bool)
        or not isinstance(value, int | float)
        or not math.isfinite(value)
        or not MIN_SPECIFICITY_TARGET <= value <= 1
    ):
        raise ValueError("Specificity/safety target must be in [0.95, 1].")


def _finite_float(value: Any, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number.")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a finite number.") from error
    if not math.isfinite(number):
        raise ValueError(f"{name} must be a finite number.")
    return number


__all__ = [
    "COSINE_GRID",
    "ENTROPY_SIMILARITY_GRID",
    "EXPLAINED_INTENSITY_GRID",
    "MATCHED_PEAKS_GRID",
    "ThresholdParameters",
    "ThresholdSearchResult",
    "estimate_bucket_probabilities",
    "evaluate_thresholds",
    "fit_thresholds",
    "iter_threshold_grid",
    "select_conflict_probability_threshold",
    "select_support_probability_threshold",
    "threshold_bucket",
    "threshold_prediction",
    "validate_threshold_parameters",
]
