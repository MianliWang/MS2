"""Strict compatibility wrapper for the MS2 threshold-grid baseline.

Deployment training belongs to ``MSAI.python.cli.ms2 train-shadow`` because it
requires nested batch/compound isolation and a locked final batch.  This module
retains a small import-compatible calibration helper for diagnostics and tests;
it accepts only the explicit truth vocabulary and always recomputes the frozen
features from bilateral spectra.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

if __package__:
    from .ms2_ml.contracts import truth_to_binary, validate_truth_row
    from .ms2_ml.features import extract_features_from_row
    from .ms2_ml.thresholds import fit_thresholds
else:  # pragma: no cover - compatibility for direct script invocation
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from MSAI.python.ms2_ml.contracts import truth_to_binary, validate_truth_row
    from MSAI.python.ms2_ml.features import extract_features_from_row
    from MSAI.python.ms2_ml.thresholds import fit_thresholds


@dataclass(frozen=True)
class CalibrationResult:
    """One constrained threshold-grid result on an eligible labeled set."""

    min_cosine: float
    min_matched_peaks: int
    min_explained_intensity: float
    min_entropy_similarity: float | None
    balanced_accuracy: float
    precision: float
    recall: float
    specificity: float
    pr_auc: float
    average_precision: float
    f0_5: float
    brier_score: float
    expected_calibration_error: float
    tp: int
    fp: int
    tn: int
    fn: int
    evaluable_rows: int
    candidates_evaluated: int
    specificity_target: float


def calibrate_thresholds(
    rows: Sequence[Mapping[str, Any]],
    label_column: str,
    *,
    spectrum_a_column: str = "peak_a_MS2",
    spectrum_b_column: str = "peak_b_MS2",
    specificity_target: float = 0.95,
    positive_values: Sequence[str] | None = None,
    objective: str | None = None,
) -> CalibrationResult:
    """Run only the fixed decision grid on explicit, spectrum-backed truth.

    ``positive_values`` and legacy objectives are accepted in the signature so
    old callers receive an actionable error instead of silently changing label
    semantics.  This helper is not a substitute for grouped deployment
    training and does not write a model artifact.
    """

    if positive_values is not None:
        raise ValueError(
            "Positive aliases are forbidden; use exact truth label 'positive_same_compound'."
        )
    if objective not in {None, "recall_at_specificity"}:
        raise ValueError(
            "The only permitted objective is recall_at_specificity with an explicit "
            "specificity target."
        )
    if not rows or not any(label_column in row for row in rows):
        raise ValueError(f"Calibration data must contain {label_column!r}.")

    examples: list[dict[str, object]] = []
    for row in rows:
        label = validate_truth_row(
            row,
            label_column,
            spectrum_a_column,
            spectrum_b_column,
        )
        binary_label = truth_to_binary(label)
        if binary_label is None:
            continue
        features = extract_features_from_row(row, spectrum_a_column, spectrum_b_column)
        if features is None:
            # ``validate_truth_row`` already catches missing spectra; this is a
            # stronger guard for malformed/non-evaluable spectrum contents.
            raise ValueError("Eligible calibration truth did not yield finite fixed features.")
        examples.append({"label": binary_label, "features": features.as_vector()})
    result = fit_thresholds(examples, specificity_target=specificity_target)
    metrics = result.metrics
    parameters = result.parameters
    return CalibrationResult(
        min_cosine=parameters.min_cosine,
        min_matched_peaks=parameters.min_matched_peaks,
        min_explained_intensity=parameters.min_explained_intensity,
        min_entropy_similarity=parameters.min_entropy_similarity,
        balanced_accuracy=(metrics.recall + metrics.specificity) / 2,
        precision=metrics.precision,
        recall=metrics.recall,
        specificity=metrics.specificity,
        pr_auc=metrics.pr_auc,
        average_precision=metrics.average_precision,
        f0_5=metrics.f0_5,
        brier_score=metrics.brier_score,
        expected_calibration_error=metrics.expected_calibration_error,
        tp=metrics.tp,
        fp=metrics.fp,
        tn=metrics.tn,
        fn=metrics.fn,
        evaluable_rows=metrics.evaluable_rows,
        candidates_evaluated=result.candidates_evaluated,
        specificity_target=result.specificity_target,
    )


def _main(argv: Sequence[str] | None = None) -> int:
    """Fail closed so the legacy CLI cannot bypass grouped shadow training."""

    parser = argparse.ArgumentParser(
        description=(
            "This standalone calibrator is disabled for model training. "
            "Use the grouped train-shadow command."
        )
    )
    parser.add_argument("--input")
    parser.add_argument("--label-column")
    parser.add_argument("--output")
    parser.parse_args(argv)
    parser.error(
        "Standalone calibration cannot write a model artifact. Use: "
        "python -m MSAI.python.cli.ms2 train-shadow --help"
    )


if __name__ == "__main__":
    _main()


__all__ = ["CalibrationResult", "calibrate_thresholds"]
