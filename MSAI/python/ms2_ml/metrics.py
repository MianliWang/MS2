"""Deterministic binary-classification and calibration metrics."""

from __future__ import annotations

import math
import random
from collections import defaultdict
from collections.abc import Hashable, Sequence
from dataclasses import asdict, dataclass
from itertools import pairwise

from .contracts import MIN_SPECIFICITY_TARGET


@dataclass(frozen=True)
class BinaryMetrics:
    """Binary performance at one prespecified probability threshold."""

    threshold: float
    evaluable_rows: int
    tp: int
    fp: int
    tn: int
    fn: int
    precision: float
    recall: float
    specificity: float
    pr_auc: float
    average_precision: float
    f0_5: float
    brier_score: float
    expected_calibration_error: float

    def as_dict(self) -> dict[str, float | int]:
        """Return a JSON-safe mapping."""

        return asdict(self)


def confusion_matrix(
    labels: Sequence[int],
    scores: Sequence[float],
    threshold: float = 0.5,
) -> tuple[int, int, int, int]:
    """Return ``(tp, fp, tn, fn)`` for ``score >= threshold``."""

    truth, probabilities = _validated(labels, scores)
    if not math.isfinite(threshold):
        raise ValueError("Classification threshold must be finite.")
    tp = fp = tn = fn = 0
    for label, score in zip(truth, probabilities, strict=True):
        prediction = score >= threshold
        if prediction and label == 1:
            tp += 1
        elif prediction:
            fp += 1
        elif label == 0:
            tn += 1
        else:
            fn += 1
    return tp, fp, tn, fn


def binary_metrics(
    labels: Sequence[int],
    scores: Sequence[float],
    threshold: float = 0.5,
    calibration_bins: int = 10,
) -> BinaryMetrics:
    """Compute confusion, discrimination, and probability-quality metrics."""

    truth, probabilities = _validated(labels, scores)
    tp, fp, tn, fn = confusion_matrix(truth, probabilities, threshold)
    precision = _divide(tp, tp + fp)
    recall = _divide(tp, tp + fn)
    specificity = _divide(tn, tn + fp)
    beta_squared = 0.25
    f0_5 = _divide(
        (1 + beta_squared) * precision * recall,
        beta_squared * precision + recall,
    )
    calibration = calibration_table(truth, probabilities, calibration_bins)
    ece = 0.0
    for row in calibration:
        count = row["count"]
        gap = row["calibration_gap"]
        if not isinstance(count, int):
            raise AssertionError("Calibration bin count must be an integer.")
        if gap is not None:
            ece += count / len(truth) * float(gap)
    return BinaryMetrics(
        threshold=float(threshold),
        evaluable_rows=len(truth),
        tp=tp,
        fp=fp,
        tn=tn,
        fn=fn,
        precision=precision,
        recall=recall,
        specificity=specificity,
        pr_auc=pr_auc(truth, probabilities),
        average_precision=average_precision(truth, probabilities),
        f0_5=f0_5,
        brier_score=sum(
            (score - label) ** 2 for label, score in zip(truth, probabilities, strict=True)
        )
        / len(truth),
        expected_calibration_error=ece,
    )


def average_precision(labels: Sequence[int], scores: Sequence[float]) -> float:
    """Return step-integrated average precision with tied scores grouped."""

    truth, probabilities = _validated(labels, scores)
    positives = sum(truth)
    if positives == 0:
        return 0.0
    points = _precision_recall_points(truth, probabilities)
    return sum(
        max(0.0, recall - previous_recall) * precision
        for (previous_recall, _previous_precision), (recall, precision) in pairwise(points)
    )


def pr_auc(labels: Sequence[int], scores: Sequence[float]) -> float:
    """Return trapezoidal area under the grouped precision-recall curve."""

    truth, probabilities = _validated(labels, scores)
    if sum(truth) == 0:
        return 0.0
    points = _precision_recall_points(truth, probabilities)
    return sum(
        max(0.0, recall - previous_recall) * (precision + previous_precision) / 2
        for (previous_recall, previous_precision), (recall, precision) in pairwise(points)
    )


def calibration_table(
    labels: Sequence[int],
    scores: Sequence[float],
    bins: int = 10,
) -> list[dict[str, float | int | None]]:
    """Return equal-width reliability bins, including explicit empty bins."""

    truth, probabilities = _validated(labels, scores)
    if bins < 2:
        raise ValueError("Calibration requires at least two bins.")
    grouped: list[list[tuple[int, float]]] = [[] for _ in range(bins)]
    for label, score in zip(truth, probabilities, strict=True):
        index = min(int(score * bins), bins - 1)
        grouped[index].append((label, score))
    rows: list[dict[str, float | int | None]] = []
    for index, members in enumerate(grouped):
        lower = index / bins
        upper = (index + 1) / bins
        if members:
            mean_prediction = sum(score for _label, score in members) / len(members)
            observed_rate = sum(label for label, _score in members) / len(members)
            gap: float | None = abs(mean_prediction - observed_rate)
        else:
            mean_prediction = observed_rate = gap = None
        rows.append(
            {
                "bin": index,
                "lower_bound": lower,
                "upper_bound": upper,
                "count": len(members),
                "mean_predicted_probability": mean_prediction,
                "observed_positive_rate": observed_rate,
                "calibration_gap": gap,
            }
        )
    return rows


def select_threshold_for_specificity(
    labels: Sequence[int],
    scores: Sequence[float],
    minimum_specificity: float = 0.95,
) -> tuple[float, BinaryMetrics]:
    """Maximize recall subject to a prespecified specificity constraint."""

    truth, probabilities = _validated(labels, scores)
    if (
        isinstance(minimum_specificity, bool)
        or not isinstance(minimum_specificity, int | float)
        or not math.isfinite(minimum_specificity)
        or not MIN_SPECIFICITY_TARGET <= minimum_specificity <= 1
    ):
        raise ValueError(
            f"minimum_specificity must be in [{MIN_SPECIFICITY_TARGET}, 1]; "
            "the model-family safety floor cannot be relaxed."
        )
    if set(truth) != {0, 1}:
        raise ValueError("Threshold selection requires positive and negative truth.")
    candidates = sorted({*probabilities, 1.0})
    scored = [binary_metrics(truth, probabilities, threshold) for threshold in candidates]
    feasible = [metric for metric in scored if metric.specificity >= minimum_specificity]
    if not feasible:
        raise ValueError(
            f"No probability threshold reaches specificity >= {minimum_specificity:.3f}."
        )
    best = max(
        feasible,
        key=lambda metric: (
            metric.recall,
            metric.precision,
            metric.specificity,
            metric.threshold,
        ),
    )
    return best.threshold, best


def group_bootstrap_ci(
    labels: Sequence[int],
    scores: Sequence[float],
    groups: Sequence[Hashable],
    threshold: float = 0.5,
    iterations: int = 1000,
    seed: int = 20260129,
    confidence_level: float = 0.95,
) -> dict[str, dict[str, float | int]]:
    """Bootstrap whole groups and report deterministic percentile intervals.

    Replicates containing only one truth class are skipped because recall or
    specificity would be structurally undefined.  The successful count is
    reported for every metric so a model card can expose unstable intervals.
    """

    truth, probabilities = _validated(labels, scores)
    if len(groups) != len(truth):
        raise ValueError("groups must have the same length as labels and scores.")
    if iterations < 1:
        raise ValueError("iterations must be at least 1.")
    if not 0 < confidence_level < 1:
        raise ValueError("confidence_level must be between 0 and 1.")
    if set(truth) != {0, 1}:
        raise ValueError("Grouped bootstrap requires both truth classes.")
    by_group: dict[Hashable, list[int]] = defaultdict(list)
    for index, group in enumerate(groups):
        try:
            by_group[group].append(index)
        except TypeError as exc:
            raise ValueError("Bootstrap groups must be hashable.") from exc
    ordered_groups = sorted(by_group, key=lambda value: (type(value).__name__, repr(value)))
    if len(ordered_groups) < 2:
        raise ValueError("Grouped bootstrap requires at least two groups.")
    base = binary_metrics(truth, probabilities, threshold)
    names = (
        "precision",
        "recall",
        "specificity",
        "pr_auc",
        "average_precision",
        "f0_5",
        "brier_score",
        "expected_calibration_error",
    )
    sampled: dict[str, list[float]] = {name: [] for name in names}
    rng = random.Random(seed)
    for _iteration in range(iterations):
        chosen = [rng.choice(ordered_groups) for _ in ordered_groups]
        indices = [index for group in chosen for index in by_group[group]]
        replicate_labels = [truth[index] for index in indices]
        if set(replicate_labels) != {0, 1}:
            continue
        replicate = binary_metrics(
            replicate_labels,
            [probabilities[index] for index in indices],
            threshold,
        )
        for name in names:
            sampled[name].append(float(getattr(replicate, name)))
    successful = len(sampled[names[0]])
    if successful == 0:
        raise ValueError("No two-class grouped bootstrap replicate was available.")
    alpha = (1 - confidence_level) / 2
    return {
        name: {
            "estimate": float(getattr(base, name)),
            "lower": _quantile(sampled[name], alpha),
            "upper": _quantile(sampled[name], 1 - alpha),
            "confidence_level": confidence_level,
            "successful_iterations": successful,
        }
        for name in names
    }


def _precision_recall_points(
    labels: tuple[int, ...],
    scores: tuple[float, ...],
) -> list[tuple[float, float]]:
    ranked = sorted(zip(scores, labels, strict=True), key=lambda item: item[0], reverse=True)
    positives = sum(labels)
    tp = fp = 0
    points = [(0.0, 1.0)]
    index = 0
    while index < len(ranked):
        score = ranked[index][0]
        while index < len(ranked) and ranked[index][0] == score:
            if ranked[index][1] == 1:
                tp += 1
            else:
                fp += 1
            index += 1
        points.append((tp / positives, _divide(tp, tp + fp)))
    return points


def _validated(
    labels: Sequence[int],
    scores: Sequence[float],
) -> tuple[tuple[int, ...], tuple[float, ...]]:
    if any(label not in {0, 1} for label in labels):
        raise ValueError("Binary labels must be exactly 0 or 1.")
    truth = tuple(int(label) for label in labels)
    probabilities = tuple(float(score) for score in scores)
    if not truth or len(truth) != len(probabilities):
        raise ValueError("Non-empty labels and scores must have the same length.")
    if any(not math.isfinite(score) or not 0 <= score <= 1 for score in probabilities):
        raise ValueError("Scores must be finite probabilities between 0 and 1.")
    return truth, probabilities


def _divide(numerator: float, denominator: float) -> float:
    return 0.0 if denominator == 0 else numerator / denominator


def _quantile(values: Sequence[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return float(ordered[lower])
    fraction = position - lower
    return float(ordered[lower] * (1 - fraction) + ordered[upper] * fraction)


__all__ = [
    "BinaryMetrics",
    "average_precision",
    "binary_metrics",
    "calibration_table",
    "confusion_matrix",
    "group_bootstrap_ci",
    "pr_auc",
    "select_threshold_for_specificity",
]
