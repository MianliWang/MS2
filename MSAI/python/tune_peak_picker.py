"""Calibrate automatic MS1 double-peak calls against reviewed Peak1/Peak2 labels."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path

try:
    from .ms1_peak_picker import (
        PeakPickingConfig,
        extract_target_eics_with_config,
        pick_chiral_peaks,
        resolution_passes_threshold,
    )
    from .ms2 import number, read_table
except ImportError:
    from ms1_peak_picker import (  # type: ignore
        PeakPickingConfig,
        extract_target_eics_with_config,
        pick_chiral_peaks,
        resolution_passes_threshold,
    )
    from ms2 import number, read_table  # type: ignore


@dataclass(frozen=True)
class PeakCallMetrics:
    balanced_accuracy: float
    precision: float
    recall: float
    specificity: float
    tp: int
    fp: int
    tn: int
    fn: int


@dataclass(frozen=True)
class RtAwarePeakCallMetrics:
    """Three-class metrics that require predicted peaks to match reviewed RTs."""

    accuracy: float
    macro_recall: float
    no_peak_recall: float
    single_peak_recall: float
    double_peak_recall: float
    close_double_recall: float
    rows: int
    correct: int
    mislocalized: int
    no_peak_rows: int
    single_peak_rows: int
    double_peak_rows: int
    close_double_rows: int


def tune_peak_picker(
    peaklist_path,
    raw_path,
    *,
    mz_column="MZ",
    peak1_column="Peak1",
    peak2_column="Peak2",
    id_column="Compound_ID",
    adaptive=False,
    review_column="IG",
    rt_tolerance_min=0.1,
):
    """Tune either the historical binary caller or an RT-aware adaptive caller."""

    if adaptive:
        return _tune_peak_picker_adaptive(
            peaklist_path,
            raw_path,
            mz_column=mz_column,
            peak1_column=peak1_column,
            peak2_column=peak2_column,
            id_column=id_column,
            review_column=review_column,
            rt_tolerance_min=rt_tolerance_min,
        )
    return _tune_peak_picker_legacy(
        peaklist_path,
        raw_path,
        mz_column=mz_column,
        peak1_column=peak1_column,
        peak2_column=peak2_column,
        id_column=id_column,
    )


def _tune_peak_picker_legacy(
    peaklist_path,
    raw_path,
    *,
    mz_column="MZ",
    peak1_column="Peak1",
    peak2_column="Peak2",
    id_column="Compound_ID",
):
    """Tune on deterministic 80% groups and report untouched 20% performance."""

    rows = read_table(Path(peaklist_path))
    required = {mz_column, peak1_column, peak2_column, id_column}
    if not rows or not required <= set(rows[0]):
        raise ValueError(f"Peaklist must contain: {', '.join(sorted(required))}")
    labeled = [
        row
        for row in rows
        if number(row.get(mz_column)) is not None and number(row.get(peak1_column)) is not None
    ]
    targets = [number(row[mz_column]) for row in labeled]
    base_config = PeakPickingConfig(eic_ppm=3.0)
    rts, traces = extract_target_eics_with_config(raw_path, targets, base_config)

    calibration = [row for row in labeled if _split(row[id_column]) != "test"]
    test = [row for row in labeled if _split(row[id_column]) == "test"]
    core_grid = list(itertools.product((5, 7, 9), (100_000.0, 200_000.0, 500_000.0), (15, 20, 25, 30)))
    quality_grid = list(itertools.product((0.5, 0.7, 0.85, 1.0), (0.0, 0.5, 1.0), (0.05, 0.1, 0.2, 0.3)))

    best = None
    cache: dict[tuple, dict[str, object]] = {}
    for sg_window, min_height, min_distance in core_grid:
        config = PeakPickingConfig(
            eic_ppm=3.0,
            sg_window=sg_window,
            sg_polyorder=2,
            min_height=min_height,
            min_distance_scans=min_distance,
        )
        picked = {
            row[id_column]: pick_chiral_peaks(rts, traces[number(row[mz_column])], config)
            for row in labeled
        }
        cache[(sg_window, min_height, min_distance)] = picked
        for max_valley_ratio, min_resolution, min_second_ratio in quality_grid:
            metrics = _evaluate(
                calibration,
                picked,
                id_column,
                peak2_column,
                max_valley_ratio,
                min_resolution,
                min_second_ratio,
            )
            candidate = (
                metrics.balanced_accuracy,
                metrics.precision,
                metrics.recall,
                -sg_window,
                -min_distance,
            )
            if best is None or candidate > best[0]:
                best = (
                    candidate,
                    config,
                    max_valley_ratio,
                    min_resolution,
                    min_second_ratio,
                    metrics,
                )
    assert best is not None
    _rank, config, max_valley_ratio, min_resolution, min_second_ratio, calibration_metrics = best
    picked = cache[(config.sg_window, config.min_height, config.min_distance_scans)]
    test_metrics = _evaluate(
        test,
        picked,
        id_column,
        peak2_column,
        max_valley_ratio,
        min_resolution,
        min_second_ratio,
    )
    rt_validation = _rt_validation(
        test,
        picked,
        id_column,
        peak1_column,
        peak2_column,
        max_valley_ratio,
        min_resolution,
        min_second_ratio,
    )
    return {
        "method": "deterministic_compound_holdout_80_20",
        "calibration_rows": len(calibration),
        "test_rows": len(test),
        "peak_picking_config": asdict(config),
        "quality_thresholds": {
            "max_valley_ratio": max_valley_ratio,
            "min_resolution": min_resolution,
            "min_second_peak_ratio": min_second_ratio,
        },
        "calibration_metrics": asdict(calibration_metrics),
        "test_metrics": asdict(test_metrics),
        "test_rt_validation": rt_validation,
        "warning": "This calibrates agreement with reviewed peak labels, not enantiomer identity.",
    }


def _tune_peak_picker_adaptive(
    peaklist_path,
    raw_path,
    *,
    mz_column,
    peak1_column,
    peak2_column,
    id_column,
    review_column,
    rt_tolerance_min,
):
    """Tune relative/prominence-based detection without looking at holdout calls.

    The raw file is still parsed once for efficiency, but peak picking during the
    grid search is restricted to calibration compounds.  Holdout traces are not
    picked until a complete detection + quality configuration has been selected.
    """

    if rt_tolerance_min <= 0:
        raise ValueError("rt_tolerance_min must be positive.")
    rows = read_table(Path(peaklist_path))
    required = {mz_column, peak1_column, peak2_column, id_column, review_column}
    if not rows or not required <= set(rows[0]):
        raise ValueError(f"Peaklist must contain: {', '.join(sorted(required))}")

    reviewed = []
    for row in rows:
        if number(row.get(mz_column)) is None or not str(row.get(review_column, "")).strip():
            continue
        if not str(row.get(id_column, "")).strip():
            continue
        peak1 = number(row.get(peak1_column))
        peak2 = number(row.get(peak2_column))
        if peak1 is None and peak2 is not None:
            raise ValueError(
                f"Reviewed row {row[id_column]!r} has Peak2 but no Peak1; labels are invalid."
            )
        reviewed.append(row)
    if not reviewed:
        raise ValueError("No explicitly reviewed rows with valid precursor m/z values were found.")
    identifiers = [str(row[id_column]) for row in reviewed]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("Adaptive tuning requires a unique identifier for each reviewed row.")

    # Split before extracting or picking so test membership is fixed before any
    # parameter-dependent operation.
    calibration = [row for row in reviewed if _split(row[id_column]) != "test"]
    test = [row for row in reviewed if _split(row[id_column]) == "test"]
    if not calibration or not test:
        raise ValueError("Deterministic split must contain both calibration and test rows.")

    targets = [float(number(row[mz_column])) for row in reviewed]
    rts, traces = extract_target_eics_with_config(
        raw_path,
        targets,
        PeakPickingConfig(eic_ppm=3.0),
    )

    best = None
    for config in _adaptive_core_grid():
        # Intentionally do not pick test rows inside this loop.
        calibration_picked = {
            row[id_column]: pick_chiral_peaks(rts, traces[float(number(row[mz_column]))], config)
            for row in calibration
        }
        for max_valley_ratio, min_resolution, min_second_ratio in _adaptive_quality_grid():
            metrics = _evaluate_rt_aware(
                calibration,
                calibration_picked,
                id_column,
                peak1_column,
                peak2_column,
                max_valley_ratio,
                min_resolution,
                min_second_ratio,
                rt_tolerance_min=rt_tolerance_min,
            )
            candidate = (
                metrics.macro_recall,
                metrics.accuracy,
                metrics.close_double_recall,
                metrics.double_peak_recall,
                metrics.single_peak_recall,
                metrics.no_peak_recall,
                -config.sg_window,
                -(config.min_support_scans or 0),
            )
            if best is None or candidate > best[0]:
                best = (
                    candidate,
                    config,
                    max_valley_ratio,
                    min_resolution,
                    min_second_ratio,
                    metrics,
                )
    assert best is not None
    (
        _rank,
        detection_config,
        max_valley_ratio,
        min_resolution,
        min_second_ratio,
        calibration_metrics,
    ) = best
    final_config = replace(
        detection_config,
        max_valley_ratio=max_valley_ratio,
        min_resolution=min_resolution,
        min_second_peak_ratio=min_second_ratio,
    )

    # This is the first parameter-dependent peak call on the test compounds.
    test_picked = {
        row[id_column]: pick_chiral_peaks(
            rts,
            traces[float(number(row[mz_column]))],
            final_config,
        )
        for row in test
    }
    test_metrics = _evaluate_rt_aware(
        test,
        test_picked,
        id_column,
        peak1_column,
        peak2_column,
        max_valley_ratio,
        min_resolution,
        min_second_ratio,
        rt_tolerance_min=rt_tolerance_min,
    )
    return {
        "method": "deterministic_compound_holdout_80_20_rt_aware_adaptive",
        "adaptive": True,
        "review_column": review_column,
        "rt_tolerance_min": rt_tolerance_min,
        "close_double_max_separation_sec": 12.0,
        "calibration_rows": len(calibration),
        "test_rows": len(test),
        "peak_picking_config": asdict(final_config),
        "quality_thresholds": {
            "max_valley_ratio": max_valley_ratio,
            "min_resolution": min_resolution,
            "min_second_peak_ratio": min_second_ratio,
        },
        "calibration_metrics": asdict(calibration_metrics),
        "test_metrics": asdict(test_metrics),
        "warning": (
            "This calibrates no/single/double calls and reviewed RT agreement, "
            "not enantiomer identity."
        ),
    }


def _adaptive_core_grid():
    """A bounded grid spanning review-only low peaks and validated height floors."""

    for sg_window, min_height, prominence, support, distance_sec, close_valley in itertools.product(
        (3, 5, 7),
        (None, 100_000.0, 200_000.0),
        (0.03, 0.1, 0.25),
        (2, 3),
        (6.0, 12.0),
        (0.85, 1.0),
    ):
        yield PeakPickingConfig(
            eic_ppm=3.0,
            sg_window=sg_window,
            sg_polyorder=2,
            min_height=min_height,
            min_distance_scans=20,
            min_prominence_ratio=prominence,
            min_support_scans=support,
            min_distance_sec=distance_sec,
            close_peak_max_valley_ratio=close_valley,
            rank_by_prominence=True,
        )


def _adaptive_quality_grid():
    """Include permissive valley/resolution choices for visibly shallow splits."""

    return itertools.product(
        (0.7, 0.85, 0.95, 1.0),
        (0.0,),
        (0.02, 0.05, 0.1, 0.2),
    )


def _split(identifier: str) -> str:
    digest = hashlib.sha256(str(identifier).encode("utf-8")).digest()
    return "test" if digest[0] % 5 == 0 else "calibration"


def _call(result, max_valley_ratio, min_resolution, min_second_ratio) -> bool:
    if result.chromatographic_status != "double_peak" or len(result.peaks) < 2:
        return False
    ratio = min(result.peaks[0].intensity, result.peaks[1].intensity) / max(
        result.peaks[0].intensity,
        result.peaks[1].intensity,
    )
    return (
        (result.valley_ratio is not None and result.valley_ratio <= max_valley_ratio)
        and resolution_passes_threshold(result.resolution, min_resolution)
        and ratio >= min_second_ratio
    )


def _evaluate(
    rows,
    picked,
    id_column,
    peak2_column,
    max_valley_ratio,
    min_resolution,
    min_second_ratio,
) -> PeakCallMetrics:
    tp = fp = tn = fn = 0
    for row in rows:
        truth = number(row.get(peak2_column)) is not None
        prediction = _call(
            picked[row[id_column]],
            max_valley_ratio,
            min_resolution,
            min_second_ratio,
        )
        if prediction and truth:
            tp += 1
        elif prediction:
            fp += 1
        elif truth:
            fn += 1
        else:
            tn += 1
    recall = _divide(tp, tp + fn)
    specificity = _divide(tn, tn + fp)
    precision = _divide(tp, tp + fp)
    return PeakCallMetrics((recall + specificity) / 2, precision, recall, specificity, tp, fp, tn, fn)


def _evaluate_rt_aware(
    rows,
    picked,
    id_column,
    peak1_column,
    peak2_column,
    max_valley_ratio,
    min_resolution,
    min_second_ratio,
    *,
    rt_tolerance_min=0.1,
    close_double_max_sec=12.0,
) -> RtAwarePeakCallMetrics:
    """Score no/single/double calls, counting wrong-RT calls as incorrect.

    A double is localized only when both predicted RTs can be paired to the two
    reviewed RTs within the tolerance.  This prevents a plausible-looking pair
    elsewhere in the chromatogram from being credited as a true positive.
    """

    support = {"no_peak": 0, "single_peak": 0, "double_peak": 0}
    localized = {"no_peak": 0, "single_peak": 0, "double_peak": 0}
    correct = mislocalized = close_rows = close_correct = 0
    for row in rows:
        peak1 = number(row.get(peak1_column))
        peak2 = number(row.get(peak2_column))
        truth = _reviewed_peak_class(peak1, peak2)
        support[truth] += 1
        result = picked[row[id_column]]
        prediction = _predicted_peak_class(
            result,
            max_valley_ratio,
            min_resolution,
            min_second_ratio,
        )

        is_close = (
            truth == "double_peak"
            and peak1 is not None
            and peak2 is not None
            and abs(peak2 - peak1) * 60 <= close_double_max_sec
        )
        if is_close:
            close_rows += 1

        location_ok = False
        if prediction == truth == "no_peak":
            location_ok = True
        elif prediction == truth == "single_peak" and peak1 is not None:
            location_ok = (
                len(result.peaks) == 1
                and abs(result.peaks[0].rt_sec / 60 - peak1) <= rt_tolerance_min
            )
        elif prediction == truth == "double_peak" and peak1 is not None and peak2 is not None:
            location_ok = _double_rts_match(
                tuple(peak.rt_sec / 60 for peak in result.peaks[:2]),
                (peak1, peak2),
                rt_tolerance_min,
            )

        if location_ok:
            correct += 1
            localized[truth] += 1
            if is_close:
                close_correct += 1
        elif prediction == truth and truth != "no_peak":
            mislocalized += 1

    recalls = {
        label: _divide(localized[label], support[label]) for label in support
    }
    supported_recalls = [recalls[label] for label in support if support[label] > 0]
    macro_recall = _divide(sum(supported_recalls), len(supported_recalls))
    return RtAwarePeakCallMetrics(
        accuracy=_divide(correct, len(rows)),
        macro_recall=macro_recall,
        no_peak_recall=recalls["no_peak"],
        single_peak_recall=recalls["single_peak"],
        double_peak_recall=recalls["double_peak"],
        close_double_recall=_divide(close_correct, close_rows),
        rows=len(rows),
        correct=correct,
        mislocalized=mislocalized,
        no_peak_rows=support["no_peak"],
        single_peak_rows=support["single_peak"],
        double_peak_rows=support["double_peak"],
        close_double_rows=close_rows,
    )


def _reviewed_peak_class(peak1, peak2) -> str:
    if peak1 is None and peak2 is None:
        return "no_peak"
    if peak1 is not None and peak2 is None:
        return "single_peak"
    if peak1 is not None and peak2 is not None:
        return "double_peak"
    raise ValueError("Peak2 cannot be reviewed without Peak1.")


def _predicted_peak_class(
    result,
    max_valley_ratio,
    min_resolution,
    min_second_ratio,
) -> str:
    if _call(result, max_valley_ratio, min_resolution, min_second_ratio):
        return "double_peak"
    if result.chromatographic_status == "single_peak" and len(result.peaks) == 1:
        return "single_peak"
    if not result.peaks:
        return "no_peak"
    return "ambiguous"


def _double_rts_match(predicted, reviewed, tolerance) -> bool:
    if len(predicted) != 2:
        return False
    direct = max(abs(predicted[0] - reviewed[0]), abs(predicted[1] - reviewed[1]))
    swapped = max(abs(predicted[0] - reviewed[1]), abs(predicted[1] - reviewed[0]))
    return min(direct, swapped) <= tolerance


def _rt_validation(
    rows,
    picked,
    id_column,
    peak1_column,
    peak2_column,
    max_valley_ratio,
    min_resolution,
    min_second_ratio,
):
    errors: list[float] = []
    positives = 0
    for row in rows:
        if number(row.get(peak2_column)) is None:
            continue
        positives += 1
        result = picked[row[id_column]]
        if not _call(result, max_valley_ratio, min_resolution, min_second_ratio):
            continue
        errors.append(
            max(
                abs(result.peaks[0].rt_sec / 60 - float(row[peak1_column])),
                abs(result.peaks[1].rt_sec / 60 - float(row[peak2_column])),
            )
        )
    errors.sort()
    return {
        "positive_rows": positives,
        "predicted_double_rows": len(errors),
        "both_rt_within_0_1_min": sum(error <= 0.1 for error in errors),
        "median_worst_rt_error_min": None if not errors else errors[len(errors) // 2],
        "p90_worst_rt_error_min": None if not errors else errors[int(0.9 * (len(errors) - 1))],
    }


def _divide(numerator, denominator):
    return 0.0 if denominator == 0 else numerator / denominator


def _main(argv=None):
    parser = argparse.ArgumentParser(description="Calibrate MSAI automatic peak picking on reviewed labels.")
    parser.add_argument("--peaklist", required=True)
    parser.add_argument("--raw", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--mz-column", default="MZ")
    parser.add_argument("--peak1-column", default="Peak1")
    parser.add_argument("--peak2-column", default="Peak2")
    parser.add_argument("--id-column", default="Compound_ID")
    parser.add_argument(
        "--adaptive",
        action="store_true",
        help="Tune relative-prominence detection with RT-aware no/single/double metrics.",
    )
    parser.add_argument(
        "--review-column",
        default="IG",
        help="Non-empty values mark rows as explicitly reviewed (default: IG).",
    )
    parser.add_argument(
        "--rt-tolerance-min",
        type=float,
        default=0.1,
        help="Maximum reviewed-vs-predicted RT error in minutes (default: 0.1).",
    )
    args = parser.parse_args(argv)
    result = tune_peak_picker(
        args.peaklist,
        args.raw,
        mz_column=args.mz_column,
        peak1_column=args.peak1_column,
        peak2_column=args.peak2_column,
        id_column=args.id_column,
        adaptive=args.adaptive,
        review_column=args.review_column,
        rt_tolerance_min=args.rt_tolerance_min,
    )
    Path(args.output).write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    _main()
