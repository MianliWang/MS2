"""Numeric validation and formatting helpers for MS1/MS2 workflows."""

from __future__ import annotations

import math


def normalize_rt_window(rt_window, rt_unit: str) -> tuple[float, float]:
    """验证二元素RT窗口并统一转换为秒。

    接受 ``min`` 或 ``sec``，拒绝非有限数、反向窗口和零宽窗口。底层读取器
    与提取算法全部使用秒，因此单位转换只应在流程边界进行一次。
    """

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
    """根据容差返回闭区间 ``[mz-delta, mz+delta]``。

    ``ppm``适合质量精度随m/z缩放的仪器设置；``Da``是固定绝对窗口；
    ``legacy_fraction``仅为旧R/Python接口兼容，例如1e-5等同于10 ppm。
    """

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
    """把表格值安全转换成有限浮点数；空值、NaN和无穷大返回``None``。"""

    try:
        number_value = float(value)
    except (TypeError, ValueError):
        return None
    return number_value if math.isfinite(number_value) else None


def csv_scalar(value) -> str:
    """把标量稳定地写成CSV文本；缺失或非有限值写为空字符串。"""

    if value is None:
        return ""
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    return f"{value:g}" if isinstance(value, float) else str(value)


def sample_sd(values: list[float]) -> float:
    """计算有限值的样本标准差；有效点不足两个时返回0。"""

    clean = [value for value in values if math.isfinite(value)]
    if len(clean) < 2:
        return 0.0
    mean = sum(clean) / len(clean)
    return math.sqrt(sum((value - mean) ** 2 for value in clean) / (len(clean) - 1))


def pearson_correlation(left: list[float], right: list[float]) -> float:
    """对成对有限观测计算Pearson相关；点数不足或方差为零时返回0。

    在MS2提取中它衡量候选fragment trace是否与目标前体EIC共同洗脱，
    不是两张最终碎片谱之间的相似度指标。两个trace必须具有相同采样点数；
    长度不一致时主动报错，避免错位数据被``zip``静默截断后产生虚高相关。
    """

    pairs = [
        (x, y) for x, y in zip(left, right, strict=True) if math.isfinite(x) and math.isfinite(y)
    ]
    if len(pairs) < 2:
        return 0.0
    xs, ys = zip(*pairs, strict=True)
    xmean = sum(xs) / len(xs)
    ymean = sum(ys) / len(ys)
    xdiff = [x - xmean for x in xs]
    ydiff = [y - ymean for y in ys]
    denominator = math.sqrt(sum(x * x for x in xdiff) * sum(y * y for y in ydiff))
    return (
        0.0
        if denominator == 0
        else sum(x * y for x, y in zip(xdiff, ydiff, strict=True)) / denominator
    )


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
