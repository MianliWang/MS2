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
    """验证输入可转换为有限浮点数，否则给出带字段名的错误。"""

    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite number.") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite.")
    return number


def _nonnegative_number(value: float, name: str) -> float:
    """验证强度/面积等输入为有限非负数。"""

    number = _finite_number(value, name)
    if number < 0:
        raise ValueError(f"{name} must be nonnegative.")
    return number


def _pseudocount(value: float, *, positive: bool = False) -> float:
    """验证pseudocount；log-ratio场景可要求严格大于零。"""

    number = _nonnegative_number(value, "pseudocount")
    if positive and number == 0:
        raise ValueError("pseudocount must be positive for log-ratio shifts.")
    return number


def _checked_ratio(numerator: float, denominator: float) -> float | None:
    """安全计算有限比值；分母为零返回``None``。"""

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
    """计算``(target+p)/(controls均值+p)``的target/control富集倍数。

    无control或分母为零返回``None``；非法强度会报错而非静默删除。该量与
    eluate/input回收倍数含义不同。
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
    """计算``(eluate+p)/(input+p)``，表示洗脱物相对输入的回收倍数。"""

    input_value = _nonnegative_number(input_intensity, "input_intensity")
    eluate_value = _nonnegative_number(eluate_intensity, "eluate_intensity")
    p = _pseudocount(pseudocount)
    numerator = eluate_value + p
    denominator = input_value + p
    if not math.isfinite(numerator) or not math.isfinite(denominator):
        raise ValueError("Intensity plus pseudocount must remain finite.")
    return _checked_ratio(numerator, denominator)


def enantiomer_fraction(area_1: float, area_2: float, pseudocount: float = 0.0) -> float | None:
    """返回Peak1面积占两峰总面积的比例；双侧均为零时返回``None``。"""

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
    """返回洗脱物Peak1比例减去输入Peak1比例，描述对映体比例方向变化。"""

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
    """计算``log2(Peak1/Peak2)``从input到eluate的变化。

    必须提供正pseudocount处理零面积；使用对数差而非显式相除以减少数值溢出。
    """

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
    """汇总重复值数量、均值、样本SD和CV；单个重复无法计算SD/CV。"""

    clean = [
        _finite_number(value, f"values[{index}]")
        for index, value in enumerate(values)
        if value is not None
    ]
    if not clean:
        return {"n": 0, "mean": None, "sd": None, "cv": None}
    mean = statistics.mean(clean)
    sd = statistics.stdev(clean) if len(clean) > 1 else None
    cv = None if sd is None or mean == 0 else sd / abs(mean)
    return {"n": len(clean), "mean": mean, "sd": sd, "cv": cv}


def benjamini_hochberg(p_values) -> list[float | None]:
    """对外部提供的p值做Benjamini-Hochberg校正并保持原顺序/空值。

    这里只控制多重检验FDR，不会从强度、倍数或重复数据自动制造p值；原始p值
    必须来自与实验设计匹配的统计检验。
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
    for rank_from_end, (value, index) in enumerate(sorted(valid, reverse=True), start=1):
        rank = count - rank_from_end + 1
        running = min(running, value * count / rank)
        output[index] = min(1.0, running)
    return output
