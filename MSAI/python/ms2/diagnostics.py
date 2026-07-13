"""Independent MS1-reference and MS2-diagnostic status vocabularies.

MS1 answers whether usable chromatographic peak coordinates were supplied.
MS2 answers whether the fragment evidence from those coordinates supports the
same molecular identity.  Keeping these vocabularies independent prevents an
MS1 double-peak call from being mistaken for an MS2 identity result.
"""

from __future__ import annotations

from dataclasses import dataclass

from .similarity import ChiralPairThresholds, SpectrumSimilarity


@dataclass(frozen=True)
class DiagnosticResult:
    """独立诊断结果：机器可读状态、问题代码和人工审核优先级。"""

    status: str
    issue_codes: tuple[str, ...]
    review_priority: str


def classify_ms1_reference(peak_a_rt, peak_b_rt) -> DiagnosticResult:
    """只描述输入表是否提供了零个、一个或两个RT坐标。

    它不重新检查EIC形状，也不把``supplied_double_peak``当作同一化合物或
    对映体证据。
    """

    has_a = _has_value(peak_a_rt)
    has_b = _has_value(peak_b_rt)
    if has_a and has_b:
        return DiagnosticResult("supplied_double_peak", (), "not_applicable")
    if has_a or has_b:
        return DiagnosticResult("supplied_single_peak", ("MS1_ONLY_ONE_RT",), "not_applicable")
    return DiagnosticResult("no_supplied_peak", ("MS1_NO_RT",), "not_applicable")


def classify_ms2_diagnostic(
    similarity: SpectrumSimilarity | None,
    legacy_status: str,
    legacy_reason: str,
    *,
    thresholds: ChiralPairThresholds | None = None,
    peak_a_quality_flags: str = "",
    peak_b_quality_flags: str = "",
    dia_window_match_count: str | int | float | None = None,
) -> DiagnosticResult:
    """把提取/相似度证据映射为MS2专属状态和显式问题代码。

    判定优先区分未尝试、无DIA覆盖、无效m/z、RT相同和空谱，再根据相似度
    结果区分支持、冲突或证据不足。附加问题代码用于提示窗口重叠、fragment
    数失衡、边界样本等，不会把MS1来源标签混入MS2诊断。
    """

    thresholds = thresholds or ChiralPairThresholds()
    issues: list[str] = []
    reason = str(legacy_reason or "")

    if legacy_status == "not_evaluated":
        issues.append("MS2_NOT_ATTEMPTED_NO_RT_PAIR")
        status = "not_attempted"
    elif "precursor_outside_acquired_DIA_windows" in reason:
        issues.append("NO_ACQUIRED_DIA_WINDOW")
        status = "not_evaluable"
    elif "invalid_precursor_mz" in reason:
        issues.append("INVALID_PRECURSOR_MZ")
        status = "not_evaluable"
    elif "retention_times_are_identical" in reason:
        issues.append("IDENTICAL_RT_COORDINATES")
        status = "not_evaluable"
    elif legacy_status == "candidate_enantiomer_pair":
        status = "supported_same_compound"
    elif legacy_status == "ms2_not_similar":
        status = "conflicting_spectra"
    elif legacy_status == "insufficient_ms2_evidence":
        status = "insufficient_evidence"
    else:
        status = "not_evaluable"

    quality = f"{peak_a_quality_flags};{peak_b_quality_flags}"
    if "no_matching_DIA_window" in quality and "NO_ACQUIRED_DIA_WINDOW" not in issues:
        issues.append("NO_ACQUIRED_DIA_WINDOW")
    if "no_precursor_signal" in quality:
        issues.append("NO_PRECURSOR_SIGNAL")
    if "no_fragments" in quality or "spectra_are_empty" in reason:
        issues.append("EMPTY_ONE_OR_BOTH_SPECTRA")

    if similarity is not None and similarity.cosine is not None:
        if similarity.matched_peaks < thresholds.min_matched_peaks:
            issues.append("TOO_FEW_MATCHED_FRAGMENTS")
        if similarity.cosine < thresholds.min_cosine:
            issues.append("LOW_COSINE_SIMILARITY")
        if (
            thresholds.min_entropy_similarity is not None
            and (similarity.entropy_similarity or 0.0) < thresholds.min_entropy_similarity
        ):
            issues.append("LOW_ENTROPY_SIMILARITY")
        explained = min(
            similarity.peak_a_explained_intensity or 0.0,
            similarity.peak_b_explained_intensity or 0.0,
        )
        if explained < thresholds.min_explained_intensity:
            issues.append("LOW_EXPLAINED_INTENSITY")
        larger = max(similarity.peak_a_count, similarity.peak_b_count)
        smaller = min(similarity.peak_a_count, similarity.peak_b_count)
        if larger >= 8 and (smaller == 0 or larger / max(smaller, 1) >= 3):
            issues.append("FRAGMENT_COUNT_IMBALANCE")
        explained_a = similarity.peak_a_explained_intensity or 0.0
        explained_b = similarity.peak_b_explained_intensity or 0.0
        if (
            similarity.matched_peaks >= thresholds.min_matched_peaks
            and abs(explained_a - explained_b) >= 0.15
            and min(explained_a, explained_b) < 0.85
        ):
            issues.append("ASYMMETRIC_UNMATCHED_INTENSITY")
        if status == "supported_same_compound" and (
            similarity.cosine < thresholds.min_cosine + 0.1
            or similarity.matched_peaks <= thresholds.min_matched_peaks + 1
            or explained < thresholds.min_explained_intensity + 0.1
        ):
            issues.append("NEAR_DECISION_BOUNDARY")

    if dia_window_match_count is not None:
        try:
            if int(float(dia_window_match_count)) > 1:
                issues.append("OVERLAPPING_DIA_WINDOWS")
        except (OverflowError, ValueError):
            pass

    issues = list(dict.fromkeys(issues))
    priority = _review_priority(status, issues)
    return DiagnosticResult(status, tuple(issues), priority)


def legacy_chiral_status(ms2_status: str) -> tuple[str, str, str]:
    """把新MS2状态转换成旧版identity/doublet字段，仅供兼容输出。"""

    if ms2_status == "supported_same_compound":
        return (
            "ms2_supported_same_compound",
            "candidate_chiral_doublet",
            "candidate_enantiomer_pair",
        )
    if ms2_status == "conflicting_spectra":
        return "ms2_conflict", "not_supported", "ms2_not_similar"
    if ms2_status == "insufficient_evidence":
        return "insufficient_ms2_evidence", "not_evaluable", "insufficient_ms2_evidence"
    return "not_evaluable", "not_evaluable", "not_evaluable"


def _review_priority(status: str, issues: list[str]) -> str:
    """依据冲突、干扰/边界、证据不足和支持审计分配P1至P5。"""

    if status == "conflicting_spectra":
        return "P1_conflict"
    if any(
        code in issues
        for code in (
            "OVERLAPPING_DIA_WINDOWS",
            "NEAR_DECISION_BOUNDARY",
            "ASYMMETRIC_UNMATCHED_INTENSITY",
            "SHARED_RT_DIA_TARGETS",
        )
    ):
        return "P2_boundary_or_interference"
    if status == "insufficient_evidence":
        return "P3_insufficient"
    if status == "supported_same_compound":
        return "P4_supported_audit"
    return "P5_not_evaluable"


def _has_value(value) -> bool:
    """判断表格字段是否不是``None``且包含非空白文本。"""

    return value is not None and str(value).strip() != ""
