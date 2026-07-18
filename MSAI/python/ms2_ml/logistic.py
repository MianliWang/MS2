"""Dependency-free deterministic StandardScaler + L2 logistic regression."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from .features import FEATURE_NAMES, feature_vector

C_GRID = (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0)
CLASS_WEIGHT_GRID: tuple[None | str, ...] = (None, "balanced")


class LogisticConvergenceError(RuntimeError):
    """Raised when the fixed deterministic Logistic solver is not releaseable."""


@dataclass(frozen=True)
class StandardScaler:
    """Population-mean/variance scaling fitted on training rows only."""

    means: tuple[float, ...]
    scales: tuple[float, ...]
    feature_names: tuple[str, ...] = FEATURE_NAMES

    def __post_init__(self) -> None:
        _validate_feature_names(self.feature_names)
        if len(self.means) != len(FEATURE_NAMES) or len(self.scales) != len(FEATURE_NAMES):
            raise ValueError("StandardScaler dimensions must match the six fixed features.")
        if not all(math.isfinite(value) for value in self.means):
            raise ValueError("StandardScaler means must be finite.")
        if not all(math.isfinite(value) and value > 0 for value in self.scales):
            raise ValueError("StandardScaler scales must be finite and positive.")

    @classmethod
    def fit(
        cls,
        rows: Sequence[Sequence[float]],
        feature_names: Sequence[str] = FEATURE_NAMES,
    ) -> StandardScaler:
        """Fit deterministic population statistics; constant columns use scale 1."""

        names = _validate_feature_names(feature_names)
        matrix = _matrix(rows)
        count = len(matrix)
        means = tuple(sum(row[column] for row in matrix) / count for column in range(len(names)))
        variances = tuple(
            sum((row[column] - means[column]) ** 2 for row in matrix) / count
            for column in range(len(names))
        )
        scales = tuple(math.sqrt(value) if value > 0 else 1.0 for value in variances)
        return cls(means=means, scales=scales, feature_names=names)

    def transform(self, row: Sequence[float]) -> tuple[float, ...]:
        """Scale one raw six-feature vector."""

        values = feature_vector(row)
        return tuple(
            (value - mean) / scale
            for value, mean, scale in zip(values, self.means, self.scales, strict=True)
        )

    def transform_many(self, rows: Sequence[Sequence[float]]) -> tuple[tuple[float, ...], ...]:
        """Scale multiple vectors without changing their order."""

        return tuple(self.transform(row) for row in rows)

    def to_dict(self) -> dict[str, list[float] | list[str]]:
        """Serialize without executable/pickle content."""

        return {
            "feature_names": list(self.feature_names),
            "means": list(self.means),
            "scales": list(self.scales),
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> StandardScaler:
        """Load and validate a JSON scaler payload."""

        return cls(
            means=tuple(float(value) for value in payload.get("means", [])),
            scales=tuple(float(value) for value in payload.get("scales", [])),
            feature_names=tuple(str(value) for value in payload.get("feature_names", [])),
        )


@dataclass(frozen=True)
class L2LogisticModel:
    """Serializable fitted binary probability model."""

    scaler: StandardScaler
    coefficients: tuple[float, ...]
    intercept: float
    c: float
    class_weight: str | None
    feature_names: tuple[str, ...] = FEATURE_NAMES
    iterations: int = 0
    converged: bool = False
    objective: float = 0.0

    def __post_init__(self) -> None:
        names = _validate_feature_names(self.feature_names)
        if names != self.scaler.feature_names:
            raise ValueError("Model and StandardScaler feature names differ.")
        if len(self.coefficients) != len(FEATURE_NAMES):
            raise ValueError("Logistic coefficients must match the six fixed features.")
        if not all(math.isfinite(value) for value in self.coefficients):
            raise ValueError("Logistic coefficients must be finite.")
        if not math.isfinite(self.intercept):
            raise ValueError("Logistic intercept must be finite.")
        _validate_c(self.c)
        _validate_class_weight(self.class_weight)
        if isinstance(self.iterations, bool) or not isinstance(self.iterations, int):
            raise ValueError("iterations must be an integer.")
        if self.iterations < 0:
            raise ValueError("iterations cannot be negative.")
        if not isinstance(self.converged, bool):
            raise ValueError("converged must be a JSON boolean.")
        if not math.isfinite(self.objective):
            raise ValueError("Logistic objective must be finite.")

    def decision_function(self, row: Sequence[float]) -> float:
        """Return the linear log-odds for one raw feature vector."""

        scaled = self.scaler.transform(row)
        return self.intercept + sum(
            coefficient * value
            for coefficient, value in zip(self.coefficients, scaled, strict=True)
        )

    def predict_proba(self, row: Sequence[float]) -> float:
        """Return same-compound probability for one raw feature vector."""

        return _sigmoid(self.decision_function(row))

    def predict_many(self, rows: Sequence[Sequence[float]]) -> tuple[float, ...]:
        """Return probabilities in input order."""

        return tuple(self.predict_proba(row) for row in rows)

    def to_dict(self) -> dict[str, Any]:
        """Serialize all preprocessing and coefficient state as plain JSON data."""

        return {
            "model_type": "standard_scaler_l2_logistic_regression",
            "feature_names": list(self.feature_names),
            "scaler": self.scaler.to_dict(),
            "coefficients": list(self.coefficients),
            "intercept": self.intercept,
            "C": self.c,
            "class_weight": self.class_weight,
            "optimizer": {
                "algorithm": "deterministic_newton_irls",
                "iterations": self.iterations,
                "converged": self.converged,
                "objective": self.objective,
            },
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> L2LogisticModel:
        """Load a model only from validated non-executable JSON values."""

        if payload.get("model_type") != "standard_scaler_l2_logistic_regression":
            raise ValueError("Unsupported logistic model_type.")
        optimizer = payload.get("optimizer")
        if not isinstance(optimizer, Mapping):
            raise ValueError("Logistic payload must contain optimizer metadata.")
        scaler_payload = payload.get("scaler")
        if not isinstance(scaler_payload, Mapping):
            raise ValueError("Logistic payload must contain scaler metadata.")
        coefficients = payload.get("coefficients")
        if not isinstance(coefficients, Sequence) or isinstance(coefficients, str | bytes):
            raise ValueError("Logistic payload coefficients must be a numeric sequence.")
        iterations = optimizer.get("iterations")
        converged = optimizer.get("converged")
        if isinstance(iterations, bool) or not isinstance(iterations, int):
            raise ValueError("Logistic optimizer iterations must be an integer.")
        if not isinstance(converged, bool):
            raise ValueError("Logistic optimizer converged must be a JSON boolean.")
        return cls(
            scaler=StandardScaler.from_dict(scaler_payload),
            coefficients=tuple(float(value) for value in coefficients),
            intercept=_required_float(payload, "intercept"),
            c=_required_float(payload, "C"),
            class_weight=_optional_string(payload.get("class_weight")),
            feature_names=tuple(str(value) for value in payload.get("feature_names", [])),
            iterations=iterations,
            converged=converged,
            objective=_required_float(optimizer, "objective"),
        )


def fit_l2_logistic(
    rows: Sequence[Sequence[float]],
    labels: Sequence[int],
    *,
    c: float = 1.0,
    class_weight: str | None = None,
    feature_names: Sequence[str] = FEATURE_NAMES,
    max_iter: int = 200,
    tolerance: float = 1e-9,
) -> L2LogisticModel:
    """Fit L2 logistic regression with deterministic Newton/IRLS updates.

    L2 regularization excludes the intercept.  Loss and Hessian are normalized
    by total sample weight so C has stable meaning across grouped folds.
    """

    names = _validate_feature_names(feature_names)
    matrix = _matrix(rows)
    truth = _labels(labels, len(matrix))
    if set(truth) != {0, 1}:
        raise ValueError("Logistic fitting requires positive and negative truth.")
    regularization_c = _validate_c(c)
    weighting = _validate_class_weight(class_weight)
    if isinstance(max_iter, bool) or not isinstance(max_iter, int) or max_iter < 1:
        raise ValueError("max_iter must be at least 1.")
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive.")

    scaler = StandardScaler.fit(matrix, names)
    scaled = scaler.transform_many(matrix)
    weights = _sample_weights(truth, weighting)
    positive_fraction = sum(
        weight * label for weight, label in zip(weights, truth, strict=True)
    ) / sum(weights)
    clipped = min(1 - 1e-12, max(1e-12, positive_fraction))
    beta = [math.log(clipped / (1 - clipped)), *([0.0] * len(names))]
    inverse_c = 1.0 / regularization_c
    current_loss = _objective(beta, scaled, truth, weights, inverse_c)
    converged = False
    iterations_run = 0
    for _iteration in range(1, max_iter + 1):
        iterations_run = _iteration
        gradient, hessian = _gradient_hessian(beta, scaled, truth, weights, inverse_c)
        if max(abs(value) for value in gradient) <= tolerance:
            converged = True
            break
        try:
            newton_step = _solve_linear_system(hessian, gradient)
        except RuntimeError as error:
            raise LogisticConvergenceError(
                "Deterministic logistic Hessian was numerically singular."
            ) from error
        directional_decrease = sum(
            gradient[index] * newton_step[index] for index in range(len(beta))
        )
        step_scale = 1.0
        accepted = False
        candidate = beta
        candidate_loss = current_loss
        for _line_search in range(40):
            candidate = [
                value - step_scale * step for value, step in zip(beta, newton_step, strict=True)
            ]
            candidate_loss = _objective(candidate, scaled, truth, weights, inverse_c)
            if candidate_loss <= current_loss - 1e-4 * step_scale * directional_decrease:
                accepted = True
                break
            step_scale *= 0.5
        if not accepted:
            raise LogisticConvergenceError(
                "Deterministic logistic line search failed to reduce loss."
            )
        max_change = max(
            abs(candidate_value - value)
            for candidate_value, value in zip(candidate, beta, strict=True)
        )
        beta = candidate
        current_loss = candidate_loss
        if max_change <= tolerance * (1 + max(abs(value) for value in beta)):
            converged = True
            break
    if not converged:
        raise LogisticConvergenceError(
            f"Deterministic logistic solver did not converge within {max_iter} iterations."
        )
    return L2LogisticModel(
        scaler=scaler,
        coefficients=tuple(beta[1:]),
        intercept=beta[0],
        c=regularization_c,
        class_weight=weighting,
        feature_names=names,
        iterations=iterations_run,
        converged=converged,
        objective=current_loss,
    )


def _gradient_hessian(
    beta: Sequence[float],
    rows: Sequence[Sequence[float]],
    labels: Sequence[int],
    weights: Sequence[float],
    inverse_c: float,
) -> tuple[list[float], list[list[float]]]:
    dimension = len(beta)
    gradient = [0.0] * dimension
    hessian = [[0.0] * dimension for _ in range(dimension)]
    total_weight = sum(weights)
    for row, label, weight in zip(rows, labels, weights, strict=True):
        augmented = (1.0, *row)
        probability = _sigmoid(sum(a * b for a, b in zip(beta, augmented, strict=True)))
        residual = weight * (probability - label) / total_weight
        curvature = weight * probability * (1 - probability) / total_weight
        for first in range(dimension):
            gradient[first] += residual * augmented[first]
            for second in range(first + 1):
                contribution = curvature * augmented[first] * augmented[second]
                hessian[first][second] += contribution
                if first != second:
                    hessian[second][first] += contribution
    for index in range(1, dimension):
        gradient[index] += inverse_c * beta[index]
        hessian[index][index] += inverse_c
    # Numerical floor affects neither the stated L2 penalty nor practical fits;
    # it only protects the unregularized intercept at extreme probabilities.
    for index in range(dimension):
        hessian[index][index] += 1e-12
    return gradient, hessian


def _objective(
    beta: Sequence[float],
    rows: Sequence[Sequence[float]],
    labels: Sequence[int],
    weights: Sequence[float],
    inverse_c: float,
) -> float:
    total_weight = sum(weights)
    data_loss = 0.0
    for row, label, weight in zip(rows, labels, weights, strict=True):
        augmented = (1.0, *row)
        logit = sum(a * b for a, b in zip(beta, augmented, strict=True))
        cross_entropy = max(logit, 0.0) - label * logit + math.log1p(math.exp(-abs(logit)))
        data_loss += weight * cross_entropy
    penalty = 0.5 * inverse_c * sum(value * value for value in beta[1:])
    return data_loss / total_weight + penalty


def _solve_linear_system(matrix: Sequence[Sequence[float]], values: Sequence[float]) -> list[float]:
    """Solve a small dense system by deterministic partial-pivot elimination."""

    size = len(values)
    augmented = [[*matrix[row], float(values[row])] for row in range(size)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) < 1e-18:
            raise RuntimeError("Logistic Hessian is numerically singular.")
        if pivot != column:
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        pivot_value = augmented[column][column]
        for entry in range(column, size + 1):
            augmented[column][entry] /= pivot_value
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor == 0:
                continue
            for entry in range(column, size + 1):
                augmented[row][entry] -= factor * augmented[column][entry]
    return [augmented[row][-1] for row in range(size)]


def _sample_weights(labels: Sequence[int], class_weight: str | None) -> tuple[float, ...]:
    if class_weight is None:
        return (1.0,) * len(labels)
    negatives = labels.count(0)
    positives = labels.count(1)
    count = len(labels)
    by_class = {0: count / (2 * negatives), 1: count / (2 * positives)}
    return tuple(by_class[label] for label in labels)


def _matrix(rows: Sequence[Sequence[float]]) -> tuple[tuple[float, ...], ...]:
    if not rows:
        raise ValueError("At least one feature row is required.")
    return tuple(feature_vector(row) for row in rows)


def _labels(labels: Sequence[int], expected: int) -> tuple[int, ...]:
    if len(labels) != expected:
        raise ValueError("labels must have the same length as feature rows.")
    if any(label not in {0, 1} for label in labels):
        raise ValueError("Logistic labels must be exactly 0 or 1.")
    return tuple(int(label) for label in labels)


def _validate_feature_names(names: Sequence[str]) -> tuple[str, ...]:
    normalized = tuple(str(name) for name in names)
    if normalized != FEATURE_NAMES:
        raise ValueError(f"Logistic feature order is frozen as {FEATURE_NAMES!r}.")
    return normalized


def _validate_c(value: float) -> float:
    numeric = float(value)
    if not math.isfinite(numeric) or numeric <= 0:
        raise ValueError("C must be finite and positive.")
    if not any(math.isclose(numeric, allowed, rel_tol=0.0, abs_tol=1e-15) for allowed in C_GRID):
        raise ValueError(f"C must be selected from {C_GRID!r}.")
    return numeric


def _validate_class_weight(value: str | None) -> str | None:
    if value not in CLASS_WEIGHT_GRID:
        raise ValueError("class_weight must be None or 'balanced'.")
    return value


def _optional_string(value: Any) -> str | None:
    return None if value is None else str(value)


def _required_float(payload: Mapping[str, Any], key: str) -> float:
    value = payload.get(key)
    if not isinstance(value, int | float) or isinstance(value, bool):
        raise ValueError(f"Logistic payload {key!r} must be numeric.")
    return float(value)


def _sigmoid(value: float) -> float:
    if value >= 0:
        inverse = math.exp(-value)
        return 1.0 / (1.0 + inverse)
    exponential = math.exp(value)
    return exponential / (1.0 + exponential)


__all__ = [
    "CLASS_WEIGHT_GRID",
    "C_GRID",
    "L2LogisticModel",
    "LogisticConvergenceError",
    "StandardScaler",
    "fit_l2_logistic",
]
