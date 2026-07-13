"""Create a reproducible review queue for two MS2 processing runs."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict, deque
from datetime import UTC, datetime
from pathlib import Path

try:
    from .ms2 import read_table, write_table
except ImportError:
    from ms2 import read_table, write_table  # type: ignore[import-not-found]


KEY_COLUMNS = ("Compound_ID", "SGC ID for Component", "SGC ID for Pool", "MZ", "Peak1", "Peak2")
SCORE_COLUMNS = (
    "ms2_cosine",
    "ms2_entropy_similarity",
    "ms2_matched_peaks",
    "peak_a_explained_intensity",
    "peak_b_explained_intensity",
)


def compare_ms2_runs(
    baseline_path,
    candidate_path,
    output_path,
    *,
    status_changes_only: bool = False,
) -> dict:
    """比较两次MS2运行并建立可复现的变化审核队列。

    按化合物、pool、MZ和RT组合匹配行，区分状态变化、谱字符串变化及二者同时
    变化，并保留baseline指标以便判断新参数究竟改善还是仅改变了结论。
    """

    baseline_path = Path(baseline_path).resolve()
    candidate_path = Path(candidate_path).resolve()
    output_path = Path(output_path).resolve()
    baseline_rows = read_table(baseline_path)
    candidate_rows = read_table(candidate_path)
    if not baseline_rows or not candidate_rows:
        raise ValueError("Both result tables must contain rows.")

    baseline_index: dict[tuple[str, ...], deque[dict]] = defaultdict(deque)
    for row in baseline_rows:
        baseline_index[_row_key(row)].append(row)

    queue: list[dict] = []
    matched = 0
    status_changes = 0
    spectrum_changes = 0
    for candidate in candidate_rows:
        key = _row_key(candidate)
        if not baseline_index[key]:
            continue
        baseline = baseline_index[key].popleft()
        matched += 1
        status_changed = baseline.get("ms2_diagnostic_status", "") != candidate.get(
            "ms2_diagnostic_status", ""
        )
        peak_a_changed = baseline.get("peak_a_MS2", "") != candidate.get("peak_a_MS2", "")
        peak_b_changed = baseline.get("peak_b_MS2", "") != candidate.get("peak_b_MS2", "")
        evidence_changed = peak_a_changed or peak_b_changed
        status_changes += int(status_changed)
        spectrum_changes += int(evidence_changed)
        if not status_changed and (status_changes_only or not evidence_changed):
            continue

        review_tier, review_flags = _review_tier(baseline, candidate, status_changed)
        row = dict(candidate)
        row.update(
            {
                "comparison_change_level": (
                    "status_and_spectrum"
                    if status_changed and evidence_changed
                    else "status_only"
                    if status_changed
                    else "spectrum_only"
                ),
                "comparison_baseline_status": baseline.get("ms2_diagnostic_status", ""),
                "comparison_candidate_status": candidate.get("ms2_diagnostic_status", ""),
                "comparison_peak_a_spectrum_changed": str(peak_a_changed).lower(),
                "comparison_peak_b_spectrum_changed": str(peak_b_changed).lower(),
                "comparison_baseline_issue_codes": baseline.get("ms2_issue_codes", ""),
                "comparison_review_tier": review_tier,
                "comparison_review_flags": ";".join(review_flags),
                "comparison_matched_peak_delta": _delta(
                    candidate.get("ms2_matched_peaks"), baseline.get("ms2_matched_peaks")
                ),
                "comparison_cosine_delta": _delta(
                    candidate.get("ms2_cosine"), baseline.get("ms2_cosine")
                ),
            }
        )
        for column in SCORE_COLUMNS:
            row[f"comparison_baseline_{column}"] = baseline.get(column, "")
        queue.append(row)

    write_table(output_path, queue)
    summary = {
        "schema_version": "msai-ms2-run-comparison-v1",
        "created_utc": datetime.now(UTC).isoformat(),
        "baseline": str(baseline_path),
        "candidate": str(candidate_path),
        "matched_rows": matched,
        "baseline_only_rows": sum(len(rows) for rows in baseline_index.values()),
        "candidate_only_rows": len(candidate_rows) - matched,
        "status_changed_rows": status_changes,
        "spectrum_changed_rows": spectrum_changes,
        "queued_rows": len(queue),
        "status_changes_only": status_changes_only,
        "output": str(output_path),
    }
    output_path.with_suffix(output_path.suffix + ".summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return summary


def _row_key(row: dict) -> tuple[str, ...]:
    """生成跨两次运行匹配同一目标的复合键。"""

    return tuple(str(row.get(column, "") or "").strip() for column in KEY_COLUMNS)


def _review_tier(baseline: dict, candidate: dict, status_changed: bool) -> tuple[str, list[str]]:
    """依据状态翻转、RT边界、fragment失衡和证据增益分配shadow审核层。"""

    flags: list[str] = []
    separation = _number(candidate.get("rt_separation_sec"))
    if separation is not None and separation < 30:
        flags.append("CLOSE_RT_PAIR_LT30S")
    for label, prefix, supplied in (
        ("PEAK1", "peak_a", "Peak1"),
        ("PEAK2", "peak_b", "Peak2"),
    ):
        apex = _number(candidate.get(prefix + "_apex_rt"))
        supplied_rt = _number(candidate.get(supplied))
        start = _number(candidate.get(prefix + "_rt_window_start"))
        end = _number(candidate.get(prefix + "_rt_window_end"))
        if apex is not None and supplied_rt is not None and abs(apex - supplied_rt * 60) > 6:
            flags.append(label + "_APEX_OFFSET_GT6S")
        if (
            apex is not None
            and start is not None
            and end is not None
            and min(abs(apex - start), abs(apex - end)) < 2
        ):
            flags.append(label + "_APEX_NEAR_WINDOW_EDGE_LT2S")
    count_a = _number(candidate.get("peak_a_fragment_count"))
    count_b = _number(candidate.get("peak_b_fragment_count"))
    if count_a is not None and count_b is not None:
        smaller, larger = min(count_a, count_b), max(count_a, count_b)
        if smaller == 0 < larger or (smaller > 0 and larger / smaller > 4):
            flags.append("FRAGMENT_COUNT_IMBALANCE_GT4X")

    matched = _number(candidate.get("ms2_matched_peaks"))
    baseline_matched = _number(baseline.get("ms2_matched_peaks")) or 0.0
    cosine = _number(candidate.get("ms2_cosine"))
    explained_a = _number(candidate.get("peak_a_explained_intensity"))
    explained_b = _number(candidate.get("peak_b_explained_intensity"))
    strong = (
        candidate.get("ms2_diagnostic_status") == "supported_same_compound"
        and cosine is not None
        and cosine >= 0.9
        and matched is not None
        and matched >= 8
        and matched - baseline_matched >= 3
        and explained_a is not None
        and explained_a >= 0.8
        and explained_b is not None
        and explained_b >= 0.8
        and not flags
    )
    if strong:
        return "strong_shadow_support_candidate", flags
    if candidate.get("ms2_diagnostic_status") == "supported_same_compound":
        return "manual_supported_boundary", flags
    if candidate.get("ms2_diagnostic_status") == "conflicting_spectra":
        return "manual_conflict", flags
    if status_changed:
        return "manual_parameter_flip", flags
    return "stable_status_spectrum_changed", flags


def _number(value) -> float | None:
    """宽松转换比较指标；无法转换时返回``None``。"""

    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _delta(candidate, baseline) -> str:
    """返回candidate减baseline；任一侧缺失时返回空文本。"""

    candidate_number, baseline_number = _number(candidate), _number(baseline)
    if candidate_number is None or baseline_number is None:
        return ""
    return f"{candidate_number - baseline_number:g}"


def main(argv=None):
    """解析baseline/candidate路径并输出运行比较摘要。"""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--status-changes-only", action="store_true")
    args = parser.parse_args(argv)
    summary = compare_ms2_runs(
        args.baseline,
        args.candidate,
        args.output,
        status_changes_only=args.status_changes_only,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
