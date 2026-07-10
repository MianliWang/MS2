"""Tune MS2 decision thresholds on an independent labeled standard set."""

from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

try:
    from .ms2 import number, read_table
except ImportError:
    from ms2 import number, read_table  # type: ignore


@dataclass(frozen=True)
class CalibrationResult:
    min_cosine: float
    min_matched_peaks: int
    min_explained_intensity: float
    min_entropy_similarity: float | None
    balanced_accuracy: float
    precision: float
    recall: float
    specificity: float
    f0_5: float
    tp: int
    fp: int
    tn: int
    fn: int
    evaluable_rows: int


def calibrate_thresholds(
    rows: list[dict],
    label_column: str,
    *,
    positive_values=("1", "true", "positive", "candidate_enantiomer_pair"),
    objective: str = "balanced_accuracy",
) -> CalibrationResult:
    """Grid-search thresholds without treating missing spectra as negatives."""

    if not rows or label_column not in rows[0]:
        raise ValueError(f"Calibration data must contain {label_column!r}.")
    positives = {str(value).strip().lower() for value in positive_values}
    examples = []
    for row in rows:
        cosine = number(row.get("ms2_cosine"))
        matched = number(row.get("ms2_matched_peaks"))
        explained_a = number(row.get("peak_a_explained_intensity"))
        explained_b = number(row.get("peak_b_explained_intensity"))
        entropy = number(row.get("ms2_entropy_similarity"))
        raw_label = str(row.get(label_column, "")).strip().lower()
        if cosine is None or matched is None or explained_a is None or explained_b is None or not raw_label:
            continue
        examples.append(
            (
                cosine,
                int(matched),
                min(explained_a, explained_b),
                entropy,
                raw_label in positives,
            )
        )
    if not examples or len({label for *_metrics, label in examples}) < 2:
        raise ValueError("Calibration requires evaluable labeled positives and negatives.")

    cosine_grid = [round(0.50 + index * 0.025, 3) for index in range(19)]
    matched_grid = [2, 3, 4, 5, 6, 8, 10]
    explained_grid = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    entropy_grid: list[float | None] = [None, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9]
    results = [
        _score(examples, cosine, matched, explained, entropy)
        for cosine, matched, explained, entropy in itertools.product(
            cosine_grid,
            matched_grid,
            explained_grid,
            entropy_grid,
        )
    ]
    if objective not in {"balanced_accuracy", "precision", "f0_5"}:
        raise ValueError("objective must be balanced_accuracy, precision, or f0_5.")
    return max(
        results,
        key=lambda result: (
            getattr(result, objective),
            result.precision,
            result.recall,
            result.min_cosine,
            result.min_matched_peaks,
        ),
    )


def _score(examples, cosine, matched, explained, entropy) -> CalibrationResult:
    tp = fp = tn = fn = 0
    for observed_cosine, observed_matched, observed_explained, observed_entropy, label in examples:
        predicted = (
            observed_cosine >= cosine
            and observed_matched >= matched
            and observed_explained >= explained
            and (entropy is None or (observed_entropy is not None and observed_entropy >= entropy))
        )
        if predicted and label:
            tp += 1
        elif predicted:
            fp += 1
        elif label:
            fn += 1
        else:
            tn += 1
    recall = _divide(tp, tp + fn)
    specificity = _divide(tn, tn + fp)
    precision = _divide(tp, tp + fp)
    balanced_accuracy = (recall + specificity) / 2
    beta_squared = 0.25
    f0_5 = _divide((1 + beta_squared) * precision * recall, beta_squared * precision + recall)
    return CalibrationResult(
        min_cosine=cosine,
        min_matched_peaks=matched,
        min_explained_intensity=explained,
        min_entropy_similarity=entropy,
        balanced_accuracy=balanced_accuracy,
        precision=precision,
        recall=recall,
        specificity=specificity,
        f0_5=f0_5,
        tp=tp,
        fp=fp,
        tn=tn,
        fn=fn,
        evaluable_rows=len(examples),
    )


def _divide(numerator: float, denominator: float) -> float:
    return 0.0 if denominator == 0 else numerator / denominator


def _main(argv=None):
    parser = argparse.ArgumentParser(description="Tune MSAI MS2 thresholds on held-out labeled standards.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--label-column", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--objective", choices=("balanced_accuracy", "precision", "f0_5"), default="balanced_accuracy")
    args = parser.parse_args(argv)
    result = calibrate_thresholds(
        read_table(Path(args.input)),
        args.label_column,
        objective=args.objective,
    )
    Path(args.output).write_text(json.dumps(asdict(result), indent=2), encoding="utf-8")
    print(json.dumps(asdict(result), ensure_ascii=False))


if __name__ == "__main__":
    _main()
