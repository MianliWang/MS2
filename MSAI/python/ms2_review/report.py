"""Build an MS2-only mirror-spectrum review report and labeling queue."""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import subprocess
import sys
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

try:
    from ..ms2.diagnostics import classify_ms1_reference, classify_ms2_diagnostic
    from ..ms2.similarity import (
        ChiralPairThresholds,
        SpectrumSimilarity,
        align_fragment_spectra,
        parse_fragment_string,
    )
except ImportError:  # Direct execution through the legacy script wrapper.
    from ms2.diagnostics import (  # type: ignore[import-not-found]
        classify_ms1_reference,
        classify_ms2_diagnostic,
    )
    from ms2.similarity import (  # type: ignore[import-not-found]
        ChiralPairThresholds,
        SpectrumSimilarity,
        align_fragment_spectra,
        parse_fragment_string,
    )


STATUS_ORDER = (
    "supported_same_compound",
    "conflicting_spectra",
    "insufficient_evidence",
    "not_evaluable",
    "not_attempted",
)
STATUS_LABELS = {
    "supported_same_compound": "支持",
    "conflicting_spectra": "冲突",
    "insufficient_evidence": "不足",
    "not_evaluable": "不可评价",
    "not_attempted": "未尝试",
}
ISSUE_LABELS = {
    "MS2_NOT_ATTEMPTED_NO_RT_PAIR": "无双 RT",
    "NO_ACQUIRED_DIA_WINDOW": "无 DIA 窗口",
    "EMPTY_ONE_OR_BOTH_SPECTRA": "单侧/双侧空谱",
    "NO_PRECURSOR_SIGNAL": "无前体信号",
    "TOO_FEW_MATCHED_FRAGMENTS": "匹配碎片少",
    "LOW_COSINE_SIMILARITY": "Cosine 低",
    "LOW_ENTROPY_SIMILARITY": "熵相似度低",
    "LOW_EXPLAINED_INTENSITY": "解释强度低",
    "FRAGMENT_COUNT_IMBALANCE": "碎片数失衡",
    "OVERLAPPING_DIA_WINDOWS": "DIA 窗口重叠",
    "NEAR_DECISION_BOUNDARY": "接近边界",
    "ASYMMETRIC_UNMATCHED_INTENSITY": "未匹配强度不对称",
    "SHARED_RT_DIA_TARGETS": "共用 RT/DIA 的多个目标",
    "INVALID_PRECURSOR_MZ": "前体 m/z 无效",
    "IDENTICAL_RT_COORDINATES": "双 RT 坐标相同",
}
PRIORITY_RANK = {
    "P1_conflict": 1,
    "P2_boundary_or_interference": 2,
    "P3_insufficient": 3,
    "P4_supported_audit": 4,
    "P5_not_evaluable": 5,
}


def read_result_rows(path: Path) -> list[dict[str, str]]:
    """以UTF-8/BOM兼容方式读取核心MS2结果CSV。"""

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def annotate_rows(
    rows: list[dict[str, str]],
    thresholds: ChiralPairThresholds | None = None,
    *,
    preserve_input_diagnostics: bool = False,
) -> list[dict]:
    """用当前阈值重建诊断字段并加入人工标注空列。

    还会识别共享相同Peak1/Peak2和DIA窗口的多个Compound_ID，标记潜在
    ``SHARED_RT_DIA_TARGETS``干扰，便于人工优先复核。
    """

    thresholds = thresholds or ChiralPairThresholds()
    annotated: list[dict] = []
    for index, source in enumerate(rows, start=1):
        # CSV列是字符串，但本地审核行会额外加入整数行号，因此这里是有意的
        # 异构表格边界；其余核心计算函数仍保持具体数值类型。
        row: dict[str, Any] = dict(source)
        if preserve_input_diagnostics:
            # Shadow-mode exports are evidence renderings, not a second diagnostic
            # pass.  Keep every supplied v2 field byte-for-byte and add only local
            # review bookkeeping fields when they are absent.
            row["review_row"] = index
            row.setdefault("manual_ms2_label", "")
            row.setdefault("manual_issue_codes", "")
            row.setdefault("manual_notes", "")
            row.setdefault("manual_reviewer", "")
            row.setdefault("manual_review_date", "")
            annotated.append(row)
            continue
        similarity = _similarity_from_row(row)
        ms1 = classify_ms1_reference(row.get("Peak1"), row.get("Peak2"))
        diagnostic = classify_ms2_diagnostic(
            similarity,
            row.get("enantiomer_pair_status", ""),
            row.get("enantiomer_pair_reason", ""),
            thresholds=thresholds,
            peak_a_quality_flags=row.get("peak_a_quality_flags", ""),
            peak_b_quality_flags=row.get("peak_b_quality_flags", ""),
            dia_window_match_count=row.get("dia_window_match_count"),
        )
        row["review_row"] = index
        row["ms1_reference_status"] = ms1.status
        row["ms1_reference_issue_codes"] = ";".join(ms1.issue_codes)
        row["ms2_diagnostic_status"] = diagnostic.status
        row["ms2_issue_codes"] = ";".join(diagnostic.issue_codes)
        row["ms2_review_priority"] = diagnostic.review_priority
        row["manual_ms2_label"] = ""
        row["manual_issue_codes"] = ""
        row["manual_notes"] = ""
        row["manual_reviewer"] = ""
        row["manual_review_date"] = ""
        annotated.append(row)
    if preserve_input_diagnostics:
        return annotated

    shared_groups: dict[tuple[str, str, str], list[dict]] = {}
    for row in annotated:
        if row.get("peak_a_MS2", "").strip() and row.get("peak_b_MS2", "").strip():
            key = (
                str(row.get("Peak1", "")),
                str(row.get("Peak2", "")),
                str(row.get("dia_window_center", "")),
            )
            shared_groups.setdefault(key, []).append(row)
    for group in shared_groups.values():
        compounds = {str(row.get("Compound_ID", "")) for row in group}
        if len(compounds) < 2:
            continue
        for row in group:
            codes = [code for code in row["ms2_issue_codes"].split(";") if code]
            if "SHARED_RT_DIA_TARGETS" not in codes:
                codes.append("SHARED_RT_DIA_TARGETS")
            row["ms2_issue_codes"] = ";".join(codes)
            if (
                PRIORITY_RANK.get(row["ms2_review_priority"], 99)
                > PRIORITY_RANK["P2_boundary_or_interference"]
            ):
                row["ms2_review_priority"] = "P2_boundary_or_interference"
    return annotated


def build_report(
    input_path: Path,
    output_dir: Path,
    *,
    template_path: Path,
    runtime_helper: Path | None = None,
    max_spectra: int = 0,
) -> dict[str, Path]:
    """生成技术报告、谱图库、人工队列、摘要和可嵌入payload。

    ``max_spectra=0``表示包含所有至少一侧有谱的目标；非零值用于快速预览。
    返回各产物路径，便于上层脚本继续发布或归档。
    """

    rows = annotate_rows(read_result_rows(input_path))
    output_dir.mkdir(parents=True, exist_ok=True)
    review_rows = [
        r for r in rows if r.get("peak_a_MS2", "").strip() or r.get("peak_b_MS2", "").strip()
    ]
    review_rows.sort(
        key=lambda r: (
            PRIORITY_RANK.get(r["ms2_review_priority"], 99),
            str(r.get("Compound_ID", "")),
        )
    )
    if max_spectra > 0:
        review_rows = review_rows[:max_spectra]

    queue_path = output_dir / "ms2_manual_review_queue.csv"
    _write_queue(queue_path, review_rows)
    summary_path = output_dir / "ms2_review_summary.json"
    summary = _summary(rows, review_rows, input_path)
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    content, payload = _report_content(rows, review_rows, input_path)
    shell = template_path.read_text(encoding="utf-8")
    shell = shell.replace("{{TITLE}}", "MS2 人工复核与问题诊断")
    shell = shell.replace(
        "{{SOURCE_AND_DATE}}",
        f"{html.escape(input_path.name)} · {datetime.now(UTC).date().isoformat()}",
    )
    shell = shell.replace("{{REPORT_CONTENT}}", content)
    shell_path = output_dir / "ms2_manual_review_report_shell.html"
    payload_path = output_dir / "ms2_manual_review_report_payload.json"
    report_path = output_dir / "ms2_manual_review_report.html"
    gallery_path = output_dir / "ms2_spectrum_gallery.html"
    shell_path.write_text(shell, encoding="utf-8")
    payload_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    if runtime_helper is not None:
        subprocess.run(  # noqa: S603 -- helper路径由本项目代码解析，不来自shell文本。
            [
                sys.executable,
                str(runtime_helper),
                "--input",
                str(shell_path),
                "--payload",
                str(payload_path),
                "--output",
                str(report_path),
            ],
            check=True,
        )
    else:
        report_path.write_text(
            shell.replace("<!-- DATA_ANALYTICS_HTML_REPORT_RUNTIME -->", ""), encoding="utf-8"
        )
    gallery_content = (
        '<article class="reading"><div class="kicker">MS2 mirror spectra</div>'
        '<header data-contract-section="title"><h1>MS2 镜像谱图人工复核册</h1></header>'
        '<p class="deck">按复核优先级排序；Peak 1 向上、Peak 2 向下。蓝/橙为 0.01 Da 内一对一匹配峰，灰色为未匹配峰。</p></article>'
        '<section class="wide"><div class="spectrum-grid">'
        + "".join(_spectrum_card(row, input_path.name, True) for row in review_rows)
        + "</div></section>"
    )
    gallery = template_path.read_text(encoding="utf-8")
    gallery = gallery.replace("{{TITLE}}", "MS2 镜像谱图人工复核册")
    gallery = gallery.replace(
        "{{SOURCE_AND_DATE}}",
        f"{html.escape(input_path.name)} · {datetime.now(UTC).date().isoformat()}",
    )
    gallery = gallery.replace("{{REPORT_CONTENT}}", gallery_content)
    gallery_path.write_text(
        gallery.replace("<!-- DATA_ANALYTICS_HTML_REPORT_RUNTIME -->", ""), encoding="utf-8"
    )
    for priority in PRIORITY_RANK:
        priority_rows = [row for row in review_rows if row["ms2_review_priority"] == priority]
        if not priority_rows:
            continue
        priority_content = (
            '<article class="reading"><div class="kicker">MS2 mirror spectra</div>'
            f'<header data-contract-section="title"><h1>{html.escape(priority)} 谱图复核</h1></header>'
            f'<p class="deck">共 {len(priority_rows)} 条；本页用于分层人工抽查。</p></article>'
            '<section class="wide"><div class="spectrum-grid">'
            + "".join(_spectrum_card(row, input_path.name, True) for row in priority_rows)
            + "</div></section>"
        )
        priority_gallery = template_path.read_text(encoding="utf-8")
        priority_gallery = priority_gallery.replace("{{TITLE}}", f"{priority} MS2 谱图复核")
        priority_gallery = priority_gallery.replace(
            "{{SOURCE_AND_DATE}}", html.escape(input_path.name)
        )
        priority_gallery = priority_gallery.replace("{{REPORT_CONTENT}}", priority_content)
        (output_dir / f"ms2_spectrum_gallery_{priority}.html").write_text(
            priority_gallery.replace("<!-- DATA_ANALYTICS_HTML_REPORT_RUNTIME -->", ""),
            encoding="utf-8",
        )
    return {
        "report": report_path,
        "gallery": gallery_path,
        "queue": queue_path,
        "summary": summary_path,
        "shell": shell_path,
        "payload": payload_path,
    }


def _similarity_from_row(row: dict[str, str]) -> SpectrumSimilarity | None:
    """从结果标量和谱字符串恢复相似度对象；缺少cosine时返回``None``。"""

    cosine = _number(row.get("ms2_cosine"))
    if cosine is None:
        return None
    return SpectrumSimilarity(
        cosine,
        int(_number(row.get("ms2_matched_peaks")) or 0),
        int(
            _number(row.get("peak_a_fragment_count"))
            or len(parse_fragment_string(row.get("peak_a_MS2")))
        ),
        int(
            _number(row.get("peak_b_fragment_count"))
            or len(parse_fragment_string(row.get("peak_b_MS2")))
        ),
        _number(row.get("peak_a_explained_intensity")),
        _number(row.get("peak_b_explained_intensity")),
        _number(row.get("ms2_entropy_similarity")),
    )


def _number(value) -> float | None:
    """把报告字段转换为有限浮点数；失败时返回``None``。"""

    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _write_queue(path: Path, rows: list[dict]) -> None:
    """写出字段受控、含人工标签空列的审核队列CSV。"""

    core = [
        "review_row",
        "Compound_ID",
        "MZ",
        "Peak1",
        "Peak2",
        "ms1_reference_status",
        "ms2_diagnostic_status",
        "ms2_review_priority",
        "ms2_issue_codes",
        "ms2_cosine",
        "ms2_entropy_similarity",
        "ms2_matched_peaks",
        "peak_a_explained_intensity",
        "peak_b_explained_intensity",
        "peak_a_fragment_count",
        "peak_b_fragment_count",
        "dia_window_match_count",
        "peak_a_MS2",
        "peak_b_MS2",
        "manual_ms2_label",
        "manual_issue_codes",
        "manual_notes",
        "manual_reviewer",
        "manual_review_date",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=core, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _summary(rows: list[dict], review_rows: list[dict], source: Path) -> dict:
    """统计状态、问题代码、单/双侧谱数量和审核队列大小。"""

    statuses = Counter(row["ms2_diagnostic_status"] for row in rows)
    issues = Counter(code for row in rows for code in row["ms2_issue_codes"].split(";") if code)
    both = sum(
        bool(r.get("peak_a_MS2", "").strip()) and bool(r.get("peak_b_MS2", "").strip())
        for r in rows
    )
    one = sum(
        bool(r.get("peak_a_MS2", "").strip()) ^ bool(r.get("peak_b_MS2", "").strip()) for r in rows
    )
    return {
        "created_utc": datetime.now(UTC).isoformat(),
        "source": str(source.resolve()),
        "rows": len(rows),
        "status_counts": dict(statuses),
        "issue_counts": dict(issues),
        "rows_with_two_spectra": both,
        "rows_with_one_spectrum": one,
        "manual_queue_rows": len(review_rows),
    }


def _source(value: str, tooltip_id: str, source_name: str) -> str:
    """为报告值生成带来源提示的HTML片段。"""

    return (
        f'<span class="source-tooltip" tabindex="0" aria-describedby="{tooltip_id}">{html.escape(value)}'
        f'<span class="source-tooltip-content" id="{tooltip_id}" role="tooltip">'
        f"Source: local MSAI analysis artifacts<br>File/Dataset: {html.escape(source_name)}</span></span>"
    )


def _report_content(rows: list[dict], review_rows: list[dict], source: Path) -> tuple[str, dict]:
    """组装报告HTML主体以及可供前端再次绘图的数据payload。"""

    source_name = source.name
    statuses = Counter(row["ms2_diagnostic_status"] for row in rows)
    issues = Counter(code for row in rows for code in row["ms2_issue_codes"].split(";") if code)
    comparable = [r for r in rows if _number(r.get("ms2_cosine")) is not None]
    two_spectra = sum(
        bool(r.get("peak_a_MS2", "").strip()) and bool(r.get("peak_b_MS2", "").strip())
        for r in rows
    )
    one_spectrum = sum(
        bool(r.get("peak_a_MS2", "").strip()) ^ bool(r.get("peak_b_MS2", "").strip()) for r in rows
    )
    supported = statuses["supported_same_compound"]
    conflict = statuses["conflicting_spectra"]
    shared_target_rows = issues["SHARED_RT_DIA_TARGETS"]
    asymmetric_rows = issues["ASYMMETRIC_UNMATCHED_INTENSITY"]

    status_data = [{"status": STATUS_LABELS[s], "count": statuses[s]} for s in STATUS_ORDER]
    issue_data = [
        {"issue": ISSUE_LABELS.get(code, code), "count": count}
        for code, count in issues.most_common(10)
    ]
    scatter_data = [
        {
            "compound": r.get("Compound_ID") or f"row-{r['review_row']}",
            "cosine": _number(r.get("ms2_cosine")),
            "entropy": _number(r.get("ms2_entropy_similarity")),
            "status": STATUS_LABELS[r["ms2_diagnostic_status"]],
        }
        for r in comparable
        if _number(r.get("ms2_entropy_similarity")) is not None
    ]

    metrics = "".join(
        f'<div class="metric"><div class="metric-label">{label}</div><div class="metric-value">{_source(str(value), f"metric-{i}", source_name)}</div><div class="metric-note">{note}</div></div>'
        for i, (label, value, note) in enumerate(
            [
                ("总记录", len(rows), "MS2 状态分母"),
                ("双侧谱图", two_spectra, "可画完整镜像谱"),
                ("支持同一化合物", supported, "不是 R/S 指认"),
                ("谱图冲突", conflict, "最高人工复核优先级"),
            ],
            start=1,
        )
    )

    table_rows = [
        (
            "<tr>"
            f"<td>{html.escape(str(r.get('Compound_ID') or r['review_row']))}</td>"
            f"<td>{html.escape(r['ms2_review_priority'])}</td>"
            f"<td>{html.escape(STATUS_LABELS[r['ms2_diagnostic_status']])}</td>"
            f"<td>{_fmt(r.get('ms2_cosine'))}</td><td>{_fmt(r.get('ms2_entropy_similarity'))}</td>"
            f"<td>{_fmt(r.get('ms2_matched_peaks'), 0)}</td>"
            f'<td class="issue">{html.escape(_issue_text(r["ms2_issue_codes"]))}</td></tr>'
        )
        for r in review_rows[:30]
    ]
    review_table = (
        '<div class="table-scroll"><table><thead><tr><th>Compound</th><th>优先级</th><th>MS2 状态</th>'
        "<th>Cosine</th><th>Entropy</th><th>匹配峰</th><th>问题</th></tr></thead><tbody>"
        + "".join(table_rows)
        + "</tbody></table></div>"
    )

    spectra = "".join(
        _spectrum_card(r, source_name, i <= 8) for i, r in enumerate(review_rows, start=1)
    )
    src = f"Source: local MSAI analysis artifacts<br>File/Dataset: {html.escape(source_name)}"
    content = f"""
    <article class="reading">
      <div class="kicker">MS2-only diagnostic · manual review</div>
      <header data-contract-section="title"><h1>MS2 人工复核与问题诊断</h1></header>
      <p class="deck">先把 MS2 身份证据画出来，再用人工标签建立本方法自己的校准标准；MS1 在本报告中只记录 RT 输入上下文，不参与 MS2 通过/失败。</p>
      <section class="summary" data-contract-section="technical-summary"><div class="summary-label">技术摘要</div><div class="summary-body">
        <p>当前标准是可复现的筛查护栏，不是已验证的 gold standard。双谱图先要求至少 6 个一对一匹配碎片，再要求 cosine ≥ 0.7 且两侧解释强度均 ≥ 0.5；entropy 仅报告。</p>
        <p>本批中有 {_source(str(two_spectra), "summary-two", source_name)} 条完整双谱、{_source(str(one_spectrum), "summary-one", source_name)} 条单侧谱。所有有谱记录已进入可填写人工标签的 CSV 队列。</p>
      </div></section>
      <section class="metrics">{metrics}</section>
    </article>

    <section class="wide" data-contract-section="key-findings">
      <h2>主要发现</h2>
      <p>最大问题不是相似度阈值本身，而是采集覆盖：大量前体没有对应 DIA 窗口；进入双谱比较后，主要失败模式是匹配碎片过少，其次才是有足够匹配峰但整体谱形冲突。</p>
      <p>人工查看镜像谱后还能看到两类自动分数不易表达的问题：{_source(str(shared_target_rows), "finding-shared", source_name)} 条记录与另一目标共用同一 RT 对和 DIA 窗口，并出现相似背景峰型；{_source(str(asymmetric_rows), "finding-asym", source_name)} 条记录两侧未匹配强度明显不对称。它们现在只提高复核优先级，不改变自动状态。</p>
      {_chart_card("status-composition", "MS2 独立诊断状态", "分母为全部非空输入记录", _bar_svg(status_data, "status", "count", "#246bce"), src)}
      {_chart_card("issue-ranking", "问题代码频次", "一个记录可以同时具有多个审计问题", _bar_svg(issue_data, "issue", "count", "#d97706"), src)}
      {_chart_card("similarity-map", "可比较双谱的 cosine–entropy 关系", "虚线为当前 cosine 0.7 护栏；entropy 暂不设通过阈值", _scatter_svg(scatter_data), src)}
    </section>

    <article class="reading">
      <section data-contract-section="scope-data-and-metric-definitions"><h2>范围、数据与指标定义</h2>
        <p><code>ms1_reference_status</code> 只说明 Peak1/Peak2 是否存在。<code>ms2_diagnostic_status</code> 只使用两侧 MS2 谱、匹配碎片、cosine、entropy 与解释强度。<code>chiral_doublet_status</code> 保留为兼容整合字段，但不作为本报告分组依据。</p>
        <p>Cosine 使用平方根强度；碎片在 0.01 Da 内一对一匹配；解释强度是已匹配峰占本侧总强度的比例。两侧碎片数严重失衡、重叠 DIA 窗口和接近边界仅触发人工审计，不暗改自动状态。</p>
      </section>
      <section data-contract-section="methodology"><h2>方法与人工复核队列</h2>
        <p>优先级依次为：谱图冲突、边界/干扰、证据不足、支持结果抽审、不可评价。下表显示前 30 条；完整队列在 <code>ms2_manual_review_queue.csv</code>，建议人工填写 <code>manual_ms2_label</code>、问题代码和备注。</p>
        <div class="card">{review_table}</div>
      </section>
    </article>
    <section class="wide"><h2>逐条 MS2 镜像谱图</h2><p>Peak 1 向上、Peak 2 向下；蓝/橙为匹配碎片，灰色为未匹配碎片。默认展开最高优先级的前八条。</p><div class="spectrum-grid">{spectra}</div></section>
    <article class="reading">
      <section data-contract-section="limitations-uncertainty-and-robustness-checks"><h2>局限、不确定性与稳健性检查</h2>
        <ul><li>当前没有独立的 MS2 阳性/阴性 truth labels，因此无法把 0.7、6 个峰和 50% 解释强度称为最优阈值。</li><li>DIA 共隔离会使两个 RT 窗口都混入背景碎片；高相似度可能来自共同干扰，低相似度也可能来自信号不足或峰窗错位。</li><li>单次候选扫描的谱图稳定性尚未验证；下一轮应比较多扫描共识、不同 RT 半窗和共洗脱阈值。</li><li>本报告不把 Peak1/Peak2 命名为 R/S，也不计算输入/洗脱物富集。</li></ul>
      </section>
      <section data-contract-section="recommended-next-steps"><h2>建议的下一轮迭代</h2>
        <ol><li>先盲审所有 P1–P3 镜像谱，再随机抽审 P4；冻结人工标签后才调参。</li><li>优先处理无 DIA 覆盖与空谱：检查离子模式/加合物、前体列、实际 isolation bounds 和 RT 对齐。</li><li>对可比较双谱做 compound-grouped 校准，联合搜索 RT 半窗、共洗脱相关阈值、fragment tolerance、consensus scans 与身份阈值，并保留独立测试批。</li><li>把重叠 DIA 窗口、共用 RT/DIA 的目标、碎片数失衡、未匹配强度不对称与谱图基峰支配纳入明确的 interference QC，而不是只依赖单一 cosine。</li><li>不要直接把六峰规则降到两峰：先用人工标签比较“绝对匹配数”与“匹配比例 + 诊断离子 + 多扫描重复性”的组合，专门解决稀疏但镜像一致的谱图。</li></ol>
      </section>
      <section data-contract-section="further-questions"><h2>后续需要回答的问题</h2><p>哪些化合物有真实同一化合物重复标准？是否有近等质量异构体/空白作为负对照？原始方法的碰撞能量、循环时间和 DIA 窗口是否在整批一致？这些答案决定最终标准能否从“筛查护栏”升级为可报告准确度的诊断标准。</p></section>
      <section class="caveat"><strong>结论边界。</strong> 本报告中的“支持同一化合物”仅是 MS2 一致性证据；手性、绝对构型和富集仍需色谱方法、标准品与输入/洗脱物重复实验独立证明。</section>
    </article>
    """

    payload = {
        "charts": [
            _bar_payload("status-composition", "MS2 独立诊断状态", status_data, "status", "count"),
            _bar_payload("issue-ranking", "问题代码频次", issue_data, "issue", "count"),
            {
                "id": "similarity-map",
                "height": 340,
                "type": "scatter",
                "dataset": {
                    "id": "similarity-map",
                    "title": "Cosine–entropy",
                    "data": scatter_data,
                    "chart_spec": {
                        "id": "similarity-map",
                        "dataset": "similarity-map",
                        "title": "Cosine–entropy",
                        "type": "scatter",
                        "encodings": {
                            "x": {"field": "cosine", "label": "Cosine", "type": "quantitative"},
                            "y": {
                                "field": "entropy",
                                "label": "Entropy similarity",
                                "type": "quantitative",
                            },
                            "color": {"field": "status", "type": "nominal"},
                            "tooltip": [{"field": "compound", "type": "nominal"}],
                        },
                        "xAxisTitle": "Cosine",
                        "yAxisTitle": "Entropy similarity",
                        "valueFormat": "decimal",
                    },
                },
            },
        ]
    }
    return content, payload


def _chart_card(chart_id: str, title: str, subtitle: str, svg: str, source_text: str) -> str:
    """包装一个带标题、说明和来源的报告图表卡片。"""

    host_class = "similarity-map-chart" if chart_id == "similarity-map" else ""
    return f'''<figure class="card source-figure"><div class="card-head"><h3>{html.escape(title)}</h3><p>{html.escape(subtitle)}</p></div><div class="chart-wrap"><div class="{host_class}" data-recharts-chart="{chart_id}"><div class="chart-fallback" data-recharts-fallback>{svg}</div><div data-recharts-live aria-hidden="true"></div></div></div><button type="button" class="source-tooltip">Source<span class="source-tooltip-content" role="tooltip">{source_text}</span></button></figure>'''


def _bar_payload(chart_id: str, title: str, data: list[dict], category: str, value: str) -> dict:
    """生成条形图的结构化payload。"""

    return {
        "id": chart_id,
        "height": 320,
        "type": "bar",
        "dataset": {
            "id": chart_id,
            "title": title,
            "data": data,
            "chart_spec": {
                "id": chart_id,
                "dataset": chart_id,
                "title": title,
                "type": "bar",
                "encodings": {
                    "x": {"field": category, "type": "nominal"},
                    "y": {"field": value, "label": "记录数", "type": "quantitative"},
                },
                "settings": {"orientation": "horizontal", "groupMode": "grouped"},
                "xAxisTitle": "",
                "yAxisTitle": "记录数",
                "valueFormat": "number",
            },
        },
    }


def _bar_svg(data: list[dict], category: str, value: str, color: str) -> str:
    """生成无需JavaScript的内联SVG条形图。"""

    width, height, left, right, top = 960, 330, 245, 60, 24
    plot = width - left - right
    row_h = min(48, (height - top - 25) / max(len(data), 1))
    maximum = max((float(row[value]) for row in data), default=1) or 1
    marks = [f'<line x1="{left}" y1="{top - 4}" x2="{left}" y2="{height - 24}" stroke="#b8c0cc"/>']
    for i, row in enumerate(data):
        y = top + i * row_h
        bar = float(row[value]) / maximum * plot
        marks.append(
            f'<text x="{left - 12}" y="{y + 17}" text-anchor="end" font-size="13" fill="#465166">{html.escape(str(row[category]))}</text>'
        )
        marks.append(
            f'<rect x="{left}" y="{y + 3}" width="{bar:.1f}" height="22" rx="4" fill="{color}" opacity=".88"/>'
        )
        marks.append(
            f'<text x="{left + bar + 8:.1f}" y="{y + 19}" font-size="13" font-weight="700" fill="#172033">{row[value]}</text>'
        )
    return f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="{html.escape(category)} counts">{"".join(marks)}</svg>'


def _scatter_svg(data: list[dict]) -> str:
    """绘制cosine与匹配fragment数的静态散点图。"""

    width, height, left, right, top, bottom = 960, 340, 75, 35, 25, 52
    pw, ph = width - left - right, height - top - bottom
    marks = []
    for tick in (0, 0.25, 0.5, 0.75, 1):
        x, y = left + tick * pw, top + (1 - tick) * ph
        marks.append(
            f'<line x1="{x}" y1="{top}" x2="{x}" y2="{top + ph}" stroke="#e8ebf0"/><text x="{x}" y="{height - 25}" text-anchor="middle" font-size="12" fill="#657085">{tick:g}</text>'
        )
        marks.append(
            f'<line x1="{left}" y1="{y}" x2="{left + pw}" y2="{y}" stroke="#e8ebf0"/><text x="{left - 10}" y="{y + 4}" text-anchor="end" font-size="12" fill="#657085">{tick:g}</text>'
        )
    threshold_x = left + 0.7 * pw
    marks.append(
        f'<line x1="{threshold_x}" y1="{top}" x2="{threshold_x}" y2="{top + ph}" stroke="#bd3d45" stroke-dasharray="6 5"/><text x="{threshold_x + 6}" y="{top + 14}" font-size="12" fill="#bd3d45">cosine 0.7</text>'
    )
    colors = {"支持同一化合物": "#25845b", "谱图冲突": "#bd3d45", "证据不足": "#d97706"}
    for row in data:
        x, y = left + row["cosine"] * pw, top + (1 - row["entropy"]) * ph
        marks.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" fill="{colors.get(row["status"], "#8790a0")}" opacity=".78"><title>{html.escape(row["compound"])}: cosine={row["cosine"]:.3f}, entropy={row["entropy"]:.3f}</title></circle>'
        )
    marks.append(
        f'<text x="{left + pw / 2}" y="{height - 5}" text-anchor="middle" font-size="13" fill="#465166">Cosine</text><text transform="translate(18 {top + ph / 2}) rotate(-90)" text-anchor="middle" font-size="13" fill="#465166">Entropy similarity</text>'
    )
    return f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Cosine versus entropy similarity">{"".join(marks)}</svg>'


def _spectrum_card(row: dict, source_name: str, opened: bool) -> str:
    """生成一个可折叠目标谱图卡片。"""

    compound = str(row.get("Compound_ID") or f"row-{row['review_row']}")
    issues = _issue_text(row["ms2_issue_codes"])
    details_open = " open" if opened else ""
    anchor = "spectrum-" + "".join(ch if ch.isalnum() or ch in "-_" else "-" for ch in compound)
    return f'''<details id="{anchor}" class="spectrum-card"{details_open}><summary>{html.escape(compound)} · {html.escape(STATUS_LABELS[row["ms2_diagnostic_status"]])} <span class="badge">{html.escape(row["ms2_review_priority"])}</span></summary><div class="spectrum-body"><div class="spectrum-meta"><span>m/z {_fmt(row.get("MZ"))}</span><span>RT {_fmt(row.get("Peak1"))} / {_fmt(row.get("Peak2"))} min</span><span>cosine {_fmt(row.get("ms2_cosine"))}</span><span>entropy {_fmt(row.get("ms2_entropy_similarity"))}</span><span>matched {_fmt(row.get("ms2_matched_peaks"), 0)}</span><span class="issue">{html.escape(issues or "无额外审计标记")}</span></div><div class="spectrum-scroll">{_mirror_spectrum_svg(row.get("peak_a_MS2", ""), row.get("peak_b_MS2", ""))}</div><button type="button" class="source-tooltip">Source<span class="source-tooltip-content" role="tooltip">Source: local MSAI analysis artifacts<br>File/Dataset: {html.escape(source_name)}<br>Row: {row["review_row"]}</span></button></div></details>'''


def _mirror_spectrum_svg(a_text: str, b_text: str, tolerance: float = 0.01) -> str:
    """为紧凑报告生成简化镜像谱；正式审核图使用共享model渲染器。"""

    a, b, pairs = align_fragment_spectra(
        a_text,
        b_text,
        fragment_mz_tol=tolerance,
        fragment_mz_tol_unit="Da",
        min_relative_intensity=0.01,
    )
    all_peaks = a + b
    if not all_peaks:
        return '<svg viewBox="0 0 900 220" role="img" aria-label="empty spectra"><text x="450" y="110" text-anchor="middle" fill="#8790a0">No fragment spectrum</text></svg>'
    ma, mb = {i for i, _ in pairs}, {j for _, j in pairs}
    min_mz, max_mz = min(x for x, _ in all_peaks), max(x for x, _ in all_peaks)
    if max_mz - min_mz < 1:
        max_mz = min_mz + 1
    width, height, left, right, mid = 900, 240, 55, 20, 120
    pw, amp = width - left - right, 88

    def scale_x(mz):
        """把fragment m/z映射到紧凑镜像谱横坐标。"""

        return left + (mz - min_mz) / (max_mz - min_mz) * pw

    max_a = max((v for _, v in a), default=1)
    max_b = max((v for _, v in b), default=1)
    marks = [f'<line x1="{left}" y1="{mid}" x2="{width - right}" y2="{mid}" stroke="#8e98a8"/>']
    for i, (mz, intensity) in enumerate(a):
        x, y = scale_x(mz), mid - intensity / max_a * amp
        marks.append(
            f'<line x1="{x:.2f}" y1="{mid}" x2="{x:.2f}" y2="{y:.2f}" stroke="{"#246bce" if i in ma else "#aeb7c5"}" stroke-width="{2 if i in ma else 1}"/>'
        )
    for j, (mz, intensity) in enumerate(b):
        x, y = scale_x(mz), mid + intensity / max_b * amp
        marks.append(
            f'<line x1="{x:.2f}" y1="{mid}" x2="{x:.2f}" y2="{y:.2f}" stroke="{"#d97706" if j in mb else "#aeb7c5"}" stroke-width="{2 if j in mb else 1}"/>'
        )
    for tick in range(5):
        mz = min_mz + (max_mz - min_mz) * tick / 4
        x = scale_x(mz)
        marks.append(
            f'<text x="{x:.1f}" y="{height - 5}" text-anchor="middle" font-size="10" fill="#697487">{mz:.1f}</text>'
        )
    marks.append(
        '<text x="8" y="32" font-size="11" fill="#246bce">Peak 1</text><text x="8" y="212" font-size="11" fill="#d97706">Peak 2</text>'
    )
    return f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Peak 1 and Peak 2 mirror spectrum">{"".join(marks)}</svg>'


def _issue_text(codes: str) -> str:
    """把机器问题代码转换为简短中文显示标签。"""

    return "；".join(ISSUE_LABELS.get(code, code) for code in str(codes or "").split(";") if code)


def _fmt(value, digits: int = 3) -> str:
    """报告数值格式化；缺失值显示为破折号。"""

    number = _number(value)
    if number is None:
        return "—"
    return f"{number:.{digits}f}" if digits else str(round(number))


def main(argv=None) -> None:
    """解析report命令参数并生成完整技术报告。"""

    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--template", type=Path, default=root / "templates" / "ms2_review_report_shell.html"
    )
    parser.add_argument(
        "--runtime-helper", type=Path, help="Optional embed_html_report_runtime.py path."
    )
    parser.add_argument(
        "--max-spectra",
        type=int,
        default=0,
        help="0 includes every row having at least one spectrum.",
    )
    args = parser.parse_args(argv)
    outputs = build_report(
        args.input,
        args.output_dir,
        template_path=args.template,
        runtime_helper=args.runtime_helper,
        max_spectra=args.max_spectra,
    )
    print("\n".join(f"{name}: {path}" for name, path in outputs.items()))


if __name__ == "__main__":
    main()
