"""Folder-view taxonomy for MS2 evidence images."""

from __future__ import annotations

try:
    from ..source_context import folder_component
except ImportError:
    from source_context import folder_component  # type: ignore[import-not-found]


def spectrum_availability(row: dict) -> str:
    """根据两侧MS2字符串是否非空判断谱图可用性。"""

    has_a = bool(str(row.get("peak_a_MS2", "") or "").strip())
    has_b = bool(str(row.get("peak_b_MS2", "") or "").strip())
    if has_a and has_b:
        return "two_spectra"
    if has_a:
        return "peak1_only"
    if has_b:
        return "peak2_only"
    return "no_spectra"


def review_views(
    row: dict,
    *,
    source_machine_id: str = "",
    source_machine_label: str = "",
    source_pool_id: str = "",
) -> list[str]:
    """返回一个目标应出现的全部独立审核视图路径。

    同一张规范图可以通过硬链接同时进入诊断状态、问题代码、机器/来源标签、
    pool和人工审核队列等文件夹。这里的分类不改变核心诊断，只帮助人工按问题
    类型抽查，例如空谱、DIA缺口、近边界、RT过近或fragment严重流失。
    """

    status = folder_component(row.get("ms2_diagnostic_status") or "unknown")
    priority = folder_component(row.get("ms2_review_priority") or "unassigned")
    availability = spectrum_availability(row)
    views = [
        f"diagnostic_status/{status}",
        f"review_priority/{priority}",
        f"spectrum_availability/{availability}",
    ]
    issues = [code for code in str(row.get("ms2_issue_codes", "") or "").split(";") if code]
    views.extend(f"issue_codes/{folder_component(code)}" for code in issues)
    if source_machine_id:
        views.append(
            f"source_machine/{folder_component(source_machine_id)}/"
            f"{folder_component(source_machine_label or 'UNLABELLED')}"
        )
    if source_pool_id:
        views.append(f"source_pool/{folder_component(source_pool_id)}")
    if priority in {"P1_conflict", "P2_boundary_or_interference", "P3_insufficient"}:
        views.append("review_queues/core_manual_review")
    if priority == "P4_supported_audit":
        views.append("review_queues/supported_random_audit")
    if availability in {"peak1_only", "peak2_only"}:
        views.append("review_queues/one_sided_spectrum")
    if "NO_ACQUIRED_DIA_WINDOW" in issues:
        views.append("review_queues/dia_coverage_gap")
    if "EMPTY_ONE_OR_BOTH_SPECTRA" in issues:
        views.append("review_queues/empty_spectrum")
    if "NEAR_DECISION_BOUNDARY" in issues:
        views.append("review_queues/near_decision_boundary")
    if "OVERLAPPING_DIA_WINDOWS" in issues or "SHARED_RT_DIA_TARGETS" in issues:
        views.append("review_queues/possible_dia_interference")
    if status == "conflicting_spectra":
        views.append("review_queues/spectral_conflict")
    if "TOO_FEW_MATCHED_FRAGMENTS" in issues:
        views.append("review_queues/low_matched_fragments")
    if "LOW_EXPLAINED_INTENSITY" in issues:
        views.append("review_queues/low_explained_intensity")
    if "ASYMMETRIC_UNMATCHED_INTENSITY" in issues:
        views.append("review_queues/asymmetric_unmatched_intensity")
    separation = _number(row.get("rt_separation_sec"))
    if separation is not None and separation <= 20:
        views.append("review_queues/close_rt_pair_le20s")
    cosine, matched = _number(row.get("ms2_cosine")), _number(row.get("ms2_matched_peaks"))
    if cosine is not None and cosine >= 0.7 and matched is not None and matched < 6:
        views.append("review_queues/high_similarity_but_sparse_lt6_matches")
    if _extreme_fragment_attrition(row):
        views.append("review_queues/extreme_fragment_attrition_lt5pct")
    if _apex_near_window_edge(row):
        views.append("review_queues/apex_near_rt_window_edge")
    if (
        row.get("ms1_reference_status") == "supplied_double_peak"
        and row.get("ms2_diagnostic_status") == "conflicting_spectra"
    ):
        views.append("cross_stage_followup/ms1_supplied_double_ms2_conflict")
    comparison_tier = str(row.get("comparison_review_tier") or "").strip()
    comparison_level = str(row.get("comparison_change_level") or "").strip()
    if comparison_tier:
        views.append(f"shadow_comparison/review_tier/{folder_component(comparison_tier)}")
    if comparison_level:
        views.append(f"shadow_comparison/change_level/{folder_component(comparison_level)}")
    ml_status = str(row.get("ml_diagnostic_status") or "").strip()
    if ml_status:
        model_id = str(row.get("ml_model_id") or "UNTRAINED").strip()
        v2_status = _normalise_v2_status(row.get("ms2_diagnostic_status"))
        views.extend(
            (
                f"ml_diagnostic_status/{folder_component(ml_status)}",
                f"ml_model/{folder_component(model_id)}",
                "shadow_comparison/transition/"
                f"{folder_component(v2_status)}_to_{folder_component(ml_status)}",
            )
        )
        if v2_status != ml_status:
            views.append("review_queues/v2_ml_disagreement")
        if ml_status == "insufficient_evidence":
            views.append("review_queues/ml_abstention")
        if ml_status == "not_evaluable":
            views.append("review_queues/ml_not_evaluable")
        abstention = str(row.get("ml_abstention_reason") or "").strip()
        if abstention:
            views.append(f"ml_abstention_reason/{folder_component(abstention)}")
        flags = [flag for flag in str(row.get("ml_review_flags") or "").split(";") if flag]
        views.extend(f"ml_review_flags/{folder_component(flag)}" for flag in flags)
    return list(dict.fromkeys(views))


def _normalise_v2_status(value) -> str:
    """Map the v2 no-attempt state onto the four-state ML comparison vocabulary."""

    status = str(value or "unknown").strip()
    return "not_evaluable" if status == "not_attempted" else status


def _extreme_fragment_attrition(row: dict) -> bool:
    """任一侧保留fragment不足候选数5%时返回真。"""

    ratios = []
    for prefix in ("peak_a", "peak_b"):
        candidate = _number(row.get(prefix + "_candidate_fragment_count"))
        retained = _number(row.get(prefix + "_fragment_count"))
        if candidate is not None and candidate > 0 and retained is not None:
            ratios.append(retained / candidate)
    return bool(ratios) and min(ratios) < 0.05


def _apex_near_window_edge(row: dict, tolerance_sec: float = 1.0) -> bool:
    """检测任一提取apex是否位于RT窗口边缘指定秒数内。"""

    for prefix in ("peak_a", "peak_b"):
        apex = _number(row.get(prefix + "_apex_rt"))
        start = _number(row.get(prefix + "_rt_window_start"))
        end = _number(row.get(prefix + "_rt_window_end"))
        if (
            apex is not None
            and start is not None
            and end is not None
            and min(abs(apex - start), abs(apex - end)) <= tolerance_sec
        ):
            return True
    return False


def _number(value):
    """用于审核分类的宽松浮点转换；失败时返回``None``。"""

    try:
        return float(value)
    except (TypeError, ValueError):
        return None
