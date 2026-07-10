"""Small, dependency-free statistics helpers for future E-ASMS batches.

The functions here only calculate explicitly named effect sizes and adjust
*supplied* p-values.  They deliberately do not infer enrichment labels or
manufacture p-values: those require a documented replicate/control design and
an appropriate statistical test chosen for that design.
"""

from __future__ import annotations

import math
import statistics


def _finite_number(value: float, name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite number.") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite.")
    return number


def _nonnegative_number(value: float, name: str) -> float:
    number = _finite_number(value, name)
    if number < 0:
        raise ValueError(f"{name} must be nonnegative.")
    return number


def _pseudocount(value: float, *, positive: bool = False) -> float:
    number = _nonnegative_number(value, "pseudocount")
    if positive and number == 0:
        raise ValueError("pseudocount must be positive for log-ratio shifts.")
    return number


def _checked_ratio(numerator: float, denominator: float) -> float | None:
    if denominator == 0:
        return None
    result = numerator / denominator
    if not math.isfinite(result):
        raise ValueError("The requested ratio is outside the finite float range.")
    return result


def target_control_enrichment_fold(
    target_intensity: float,
    control_intensities,
    pseudocount: float = 0.0,
) -> float | None:
    """Return ``(target + p) / (mean(controls) + p)``.

    ``None`` means that no controls were supplied or that the denominator is
    zero.  Invalid intensities are rejected instead of being silently removed.
    """

    target = _nonnegative_number(target_intensity, "target_intensity")
    p = _pseudocount(pseudocount)
    controls = [
        _nonnegative_number(value, f"control_intensities[{index}]")
        for index, value in enumerate(control_intensities)
    ]
    if not controls:
        return None
    numerator = target + p
    denominator = statistics.mean(controls) + p
    if not math.isfinite(numerator) or not math.isfinite(denominator):
        raise ValueError("Intensity plus pseudocount must remain finite.")
    return _checked_ratio(numerator, denominator)


# Backwards-compatible name.  New code should use the explicit name above so
# target/control enrichment cannot be confused with eluate/input recovery.
enrichment_fold = target_control_enrichment_fold


def input_eluate_recovery_fold(
    input_intensity: float,
    eluate_intensity: float,
    pseudocount: float = 0.0,
) -> float | None:
    """Return eluate/input recovery as ``(eluate + p) / (input + p)``."""

    input_value = _nonnegative_number(input_intensity, "input_intensity")
    eluate_value = _nonnegative_number(eluate_intensity, "eluate_intensity")
    p = _pseudocount(pseudocount)
    numerator = eluate_value + p
    denominator = input_value + p
    if not math.isfinite(numerator) or not math.isfinite(denominator):
        raise ValueError("Intensity plus pseudocount must remain finite.")
    return _checked_ratio(numerator, denominator)


def enantiomer_fraction(area_1: float, area_2: float, pseudocount: float = 0.0) -> float | None:
    """Return enantiomer-1's fraction of the two-enantiomer signal."""

    first = _nonnegative_number(area_1, "area_1")
    second = _nonnegative_number(area_2, "area_2")
    p = _pseudocount(pseudocount)
    first += p
    second += p
    if not math.isfinite(first) or not math.isfinite(second):
        raise ValueError("Area plus pseudocount must remain finite.")
    scale = max(first, second)
    if scale == 0:
        return None
    # Scaling avoids overflow in first + second for very large finite areas.
    scaled_first = first / scale
    return scaled_first / (scaled_first + second / scale)


def enantiomer_fraction_shift(
    input_area_1: float,
    input_area_2: float,
    eluate_area_1: float,
    eluate_area_2: float,
    pseudocount: float = 0.0,
) -> float | None:
    """Return eluate fraction(enantiomer 1) minus its input fraction."""

    input_fraction = enantiomer_fraction(input_area_1, input_area_2, pseudocount)
    eluate_fraction = enantiomer_fraction(eluate_area_1, eluate_area_2, pseudocount)
    if input_fraction is None or eluate_fraction is None:
        return None
    return eluate_fraction - input_fraction


def enantiomer_log2_ratio_shift(
    input_area_1: float,
    input_area_2: float,
    eluate_area_1: float,
    eluate_area_2: float,
    pseudocount: float,
) -> float:
    first_input = _nonnegative_number(input_area_1, "input_area_1")
    second_input = _nonnegative_number(input_area_2, "input_area_2")
    first_eluate = _nonnegative_number(eluate_area_1, "eluate_area_1")
    second_eluate = _nonnegative_number(eluate_area_2, "eluate_area_2")
    p = _pseudocount(pseudocount, positive=True)
    adjusted = [first_input + p, second_input + p, first_eluate + p, second_eluate + p]
    if not all(math.isfinite(value) for value in adjusted):
        raise ValueError("Area plus pseudocount must remain finite.")
    # Differences of logarithms avoid overflow/underflow in explicit ratios.
    input_log_ratio = math.log2(adjusted[0]) - math.log2(adjusted[1])
    eluate_log_ratio = math.log2(adjusted[2]) - math.log2(adjusted[3])
    return eluate_log_ratio - input_log_ratio


def replicate_summary(values) -> dict:
    clean = [_finite_number(value, f"values[{index}]") for index, value in enumerate(values) if value is not None]
    if not clean:
        return {"n": 0, "mean": None, "sd": None, "cv": None}
    mean = statistics.mean(clean)
    sd = statistics.stdev(clean) if len(clean) > 1 else None
    cv = None if sd is None or mean == 0 else sd / abs(mean)
    return {"n": len(clean), "mean": mean, "sd": sd, "cv": cv}


def benjamini_hochberg(p_values) -> list[float | None]:
    """Adjust supplied p-values with BH; preserve ``None`` entries and order.

    This is only a multiple-testing correction.  It does not calculate p-values
    from intensities, fold changes, or replicates.
    """

    values = list(p_values)
    valid = []
    for index, raw in enumerate(values):
        if raw is None:
            continue
        value = _finite_number(raw, f"p_values[{index}]")
        if not 0 <= value <= 1:
            raise ValueError("p-values must be between 0 and 1.")
        valid.append((value, index))
    count = len(valid)
    output: list[float | None] = [None] * len(values)
    running = 1.0
    for rank_from_end, (value, index) in enumerate(reversed(sorted(valid)), start=1):
        rank = count - rank_from_end + 1
        running = min(running, value * count / rank)
        output[index] = min(1.0, running)
    return output
