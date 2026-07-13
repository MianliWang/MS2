"""Classification and review-bucket helpers for MS1 EIC images."""

from __future__ import annotations

import math
import statistics

try:
    from ..source_context import (
        folder_component,
    )
    from ..source_context import (
        source_machine_label_interpretation as source_machine_label_interpretation,
    )
except ImportError:
    from source_context import (  # type: ignore[import-not-found]
        folder_component,
    )
    from source_context import (
        source_machine_label_interpretation as source_machine_label_interpretation,
    )


def reference_chromatographic_status(row: dict) -> str:
    """仅依据人工/上游提供的Peak1、Peak2判断参考单峰、双峰或无峰。"""

    has_first = _number(row.get("Peak1")) is not None
    has_second = _number(row.get("Peak2")) is not None
    if has_first and has_second:
        return "double_peak"
    if has_first or has_second:
        return "single_peak"
    return "no_peak"


def manual_chromatographic_status(row: dict) -> str:
    """``reference_chromatographic_status``的旧名称兼容别名。"""

    return reference_chromatographic_status(row)


def prediction_diagnostic(
    row: dict,
    predicted_status: str,
    predicted_rts_min: list[float],
    *,
    rt_tolerance_min: float = 0.1,
) -> str:
    """同时比较自动类别和RT定位，但不把MS1结果当作MS2身份依据。

    类别一致但预测RT超出容差会返回``rt_mislocalized``，防止算法在色谱图
    其他位置找到两个峰却被误计为正确。
    """

    reference_status = reference_chromatographic_status(row)
    reference_rts = [
        value for name in ("Peak1", "Peak2") if (value := _number(row.get(name))) is not None
    ]
    expected_auto = "not_evaluable" if reference_status == "no_peak" else reference_status
    if predicted_status != expected_auto:
        return "class_mismatch"
    if reference_status == "no_peak":
        return "exact_agreement"
    if len(reference_rts) != len(predicted_rts_min):
        return "class_mismatch"
    if all(
        abs(left - right) <= rt_tolerance_min
        for left, right in zip(reference_rts, predicted_rts_min, strict=True)
    ):
        return "exact_agreement"
    return "rt_mislocalized"


def background_diagnostics(
    rts_sec: list[float],
    intensities: list[float],
    protected_rts_min: list[float],
    *,
    exclusion_half_width_sec: float = 20.0,
    absolute_height_floor: float = 200_000.0,
) -> dict[str, float | bool]:
    """汇总峰外背景干净程度，并标记“低但干净”的人工复核候选。

    在保护RT区域之外计算median、MAD、p95/p99、非零比例和峰/background p99；
    低于绝对高度门槛但背景p99很低、非零点少且半高支持充分时标记
    ``low_clean_candidate``，而不是直接归为无信号。
    """

    raw = [max(0.0, float(value)) for value in intensities]
    maximum = max(raw, default=0.0)
    protected_sec = [value * 60 for value in protected_rts_min]
    background = [
        intensity
        for rt, intensity in zip(rts_sec, raw, strict=True)
        if all(abs(rt - center) > exclusion_half_width_sec for center in protected_sec)
    ]
    if not background:
        background = raw
    median = statistics.median(background) if background else 0.0
    mad = statistics.median(abs(value - median) for value in background) if background else 0.0
    p95 = _quantile(background, 0.95)
    p99 = _quantile(background, 0.99)
    nonzero_fraction = sum(value > 0 for value in background) / max(len(background), 1)
    cleanliness_ratio = maximum / max(p99, 1.0)
    half_height_support = _half_height_support(raw)
    low_clean = (
        0 < maximum < absolute_height_floor
        and p99 <= max(1.0, maximum * 0.05)
        and nonzero_fraction <= 0.15
        and half_height_support >= 3
    )
    return {
        "background_median": median,
        "background_mad": mad,
        "background_p95": p95,
        "background_p99": p99,
        "background_nonzero_fraction": nonzero_fraction,
        "peak_to_background_p99": cleanliness_ratio,
        "max_half_height_support_scans": half_height_support,
        "low_clean_candidate": low_clean,
    }


def review_views(
    row: dict,
    *,
    reference_status: str,
    baseline_status: str,
    baseline_diagnostic: str,
    experimental_diagnostic: str,
    low_clean_candidate: bool,
    parameter_stability: str,
    source_machine_id: str = "",
    source_machine_label: str = "",
    source_machine_interpretation: str = "",
) -> list[str]:
    """按参考类别、自动结果、参数稳定性和背景质量返回MS1审核文件夹。"""

    views = [
        f"reference_ms1/supplied_rt/{reference_status}",
        f"auto_ms1/baseline/{baseline_status}",
        f"reference_comparison/baseline/{baseline_diagnostic}",
    ]
    if parameter_stability == "stable_class_and_rt":
        views.append(f"auto_ms1/consensus_stable/{baseline_status}")
    else:
        views.append(f"review_queues/{parameter_stability}")
    first, second = _number(row.get("Peak1")), _number(row.get("Peak2"))
    if first is not None and second is not None and abs(second - first) * 60 < 12:
        views.append("review_queues/close_or_shoulder_double")
    if baseline_status == "ambiguous":
        views.append("review_queues/baseline_ambiguous")
    if baseline_diagnostic != "exact_agreement":
        views.append("review_queues/baseline_vs_reference")
    if experimental_diagnostic != "exact_agreement":
        views.append("review_queues/experimental_vs_reference")
    if baseline_diagnostic == "exact_agreement" and experimental_diagnostic != "exact_agreement":
        views.append("review_queues/experimental_regression")
    if baseline_diagnostic != "exact_agreement" and experimental_diagnostic == "exact_agreement":
        views.append("review_queues/experimental_improvement")
    if low_clean_candidate:
        views.append("review_queues/low_clean_signal")
    if source_machine_id:
        source_id = folder_component(source_machine_id)
        source_label = folder_component(source_machine_label or "UNLABELLED")
        views.append(f"source_machine/{source_id}/{source_label}")
        if source_machine_interpretation == "unusable_no_peak_or_multiple_peak":
            views.append(f"review_queues/source_machine_{source_id}_red_unusable")
        elif source_machine_interpretation == "manual_review_needed":
            views.append(f"review_queues/source_machine_{source_id}_check")
    return list(dict.fromkeys(views))


def _quantile(values: list[float], probability: float) -> float:
    """使用线性插值计算小型背景数组分位数。"""

    if not values:
        return 0.0
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def _half_height_support(values: list[float]) -> int:
    """返回全局最高峰周围连续半高以上scan数量。"""

    if not values or max(values) <= 0:
        return 0
    apex = max(range(len(values)), key=values.__getitem__)
    threshold = values[apex] / 2
    left = apex
    while left > 0 and values[left - 1] >= threshold:
        left -= 1
    right = apex
    while right < len(values) - 1 and values[right + 1] >= threshold:
        right += 1
    return right - left + 1


def _number(value) -> float | None:
    """安全转换有限浮点数；失败时返回``None``。"""

    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None
