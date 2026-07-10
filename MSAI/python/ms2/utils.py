"""Numeric validation and formatting helpers for MS1/MS2 workflows."""

from __future__ import annotations

import math


def normalize_rt_window(rt_window, rt_unit: str) -> tuple[float, float]:
    """Validate an RT window and convert it to seconds."""

    if rt_unit not in {"min", "sec"}:
        raise ValueError("rt_unit must be either 'min' or 'sec'.")
    if len(rt_window) != 2:
        raise ValueError("rt_window must contain exactly two values.")
    start, end = (float(rt_window[0]), float(rt_window[1]))
    if not (math.isfinite(start) and math.isfinite(end)):
        raise ValueError("rt_window must contain finite values.")
    if start >= end:
        raise ValueError("rt_window start must be smaller than rt_window end.")
    return (start * 60, end * 60) if rt_unit == "min" else (start, end)


def mz_bounds(mz: float, mz_tol: float, mz_tol_unit: str) -> tuple[float, float]:
    """Return the inclusive m/z interval for ppm, Da, or legacy tolerance."""

    mz = float(mz)
    mz_tol = float(mz_tol)
    if not (math.isfinite(mz) and math.isfinite(mz_tol)) or mz_tol < 0:
        raise ValueError("mz and mz_tol must be finite, and mz_tol must be non-negative.")
    if mz_tol_unit == "legacy_fraction":
        delta = mz * mz_tol
    elif mz_tol_unit == "ppm":
        delta = mz * mz_tol / 1_000_000
    elif mz_tol_unit == "Da":
        delta = mz_tol
    else:
        raise ValueError("mz_tol_unit must be one of 'legacy_fraction', 'ppm', or 'Da'.")
    return mz - delta, mz + delta


def number(value) -> float | None:
    """Convert a value to a finite float, returning ``None`` on failure."""

    try:
        number_value = float(value)
    except (TypeError, ValueError):
        return None
    return number_value if math.isfinite(number_value) else None


def csv_scalar(value) -> str:
    """Format a scalar for stable CSV output."""

    if value is None:
        return ""
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    return f"{value:g}" if isinstance(value, float) else str(value)


def sample_sd(values: list[float]) -> float:
    """Sample standard deviation, or zero for fewer than two finite values."""

    clean = [value for value in values if math.isfinite(value)]
    if len(clean) < 2:
        return 0.0
    mean = sum(clean) / len(clean)
    return math.sqrt(sum((value - mean) ** 2 for value in clean) / (len(clean) - 1))


def pearson_correlation(left: list[float], right: list[float]) -> float:
    """Pearson correlation over paired finite observations."""

    pairs = [(x, y) for x, y in zip(left, right) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 2:
        return 0.0
    xs, ys = zip(*pairs)
    xmean = sum(xs) / len(xs)
    ymean = sum(ys) / len(ys)
    xdiff = [x - xmean for x in xs]
    ydiff = [y - ymean for y in ys]
    denominator = math.sqrt(sum(x * x for x in xdiff) * sum(y * y for y in ydiff))
    return 0.0 if denominator == 0 else sum(x * y for x, y in zip(xdiff, ydiff)) / denominator


# Private compatibility aliases retained for the old monolithic module.
_sd = sample_sd
_cor = pearson_correlation


__all__ = [
    "csv_scalar",
    "mz_bounds",
    "normalize_rt_window",
    "number",
    "pearson_correlation",
    "sample_sd",
]
