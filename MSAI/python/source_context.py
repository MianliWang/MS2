"""Neutral source-machine and pool context shared by MS1 and MS2 exports."""

from __future__ import annotations


def source_machine_label_interpretation(label) -> str:
    """Interpret the confirmed IG colour vocabulary without claiming truth.

    Other source machines (IA, IB, ...) must retain their raw labels until the
    laboratory confirms that they use the same vocabulary.
    """

    label = str(label or "").strip().upper()
    return {
        "GREEN": "expected_double_peak",
        "YELLOW": "expected_single_peak",
        "RED": "unusable_no_peak_or_multiple_peak",
        "CHECK": "manual_review_needed",
    }.get(label, "unrecognised_or_unlabelled")


def folder_component(value, *, default: str = "unknown", limit: int = 80) -> str:
    cleaned = "".join(
        character if character.isalnum() or character in "-_" else "_"
        for character in str(value or "")
    )
    return cleaned[:limit] or default
