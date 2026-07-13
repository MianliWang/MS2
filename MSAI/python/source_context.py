"""Neutral source-machine and pool context shared by MS1 and MS2 exports."""

from __future__ import annotations


def source_machine_label_interpretation(label) -> str:
    """解释已确认的IG颜色词汇，但不把来源标签当作分析真值。

    GREEN/YELLOW/RED/CHECK是IG列中的预期/复核标签，不是机器名称。IA、IB等
    其他实验必须保留原始标签，直到实验室确认它们使用完全相同的词汇体系。
    """

    label = str(label or "").strip().upper()
    return {
        "GREEN": "expected_double_peak",
        "YELLOW": "expected_single_peak",
        "RED": "unusable_no_peak_or_multiple_peak",
        "CHECK": "manual_review_needed",
    }.get(label, "unrecognised_or_unlabelled")


def folder_component(value, *, default: str = "unknown", limit: int = 80) -> str:
    """把任意实验、pool或化合物文本清理为安全且长度受限的目录组件。"""

    cleaned = "".join(
        character if character.isalnum() or character in "-_" else "_"
        for character in str(value or "")
    )
    return cleaned[:limit] or default
