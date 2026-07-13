"""Dependency-free SVG rendering for targeted MS1 chromatograms."""

from __future__ import annotations

import html
import math


def render_eic_svg(
    *,
    compound_id: str,
    target_mz: float,
    ppm: float,
    rts_sec: list[float],
    raw: list[float],
    baseline_smooth: list[float],
    experimental_smooth: list[float],
    manual_rts_min: list[float],
    baseline_peaks,
    experimental_peaks,
    review_candidates,
    metadata: dict,
) -> str:
    """渲染同时包含全程、局部放大和log视图的MS1 EIC审核SVG。

    原始信号、baseline/experimental平滑、人工RT、自动峰和高召回候选使用
    不同线型/符号。局部放大优先以人工RT为锚点，避免错误自动峰把真正的近峰
    区域压缩得不可读。
    """

    width, height = 1240, 930
    full_box = (78, 260, 1120, 240)
    focus_box = (78, 550, 1120, 210)
    log_box = (78, 810, 1120, 80)
    rt_min = [value / 60 for value in rts_sec]
    focus = _focus_bounds(rt_min, manual_rts_min, baseline_peaks, experimental_peaks)
    threshold = metadata.get("min_height")

    body = [
        '<rect width="1240" height="930" fill="#ffffff"/>',
        "<style>text{font-family:Segoe UI,Arial,sans-serif}.mono{font-family:Consolas,monospace}</style>",
        f'<text x="42" y="44" font-size="25" font-weight="700" fill="#172033">{html.escape(compound_id)} · targeted MS1 EIC</text>',
        f'<text x="42" y="72" class="mono" font-size="14" fill="#4f5c70">m/z {target_mz:.6f} · {ppm:g} ppm [{target_mz * (1 - ppm / 1e6):.6f}, {target_mz * (1 + ppm / 1e6):.6f}]</text>',
        f'<text x="42" y="98" font-size="13" fill="#25334a">supplied RT reference: {html.escape(str(metadata.get("reference_status", "")))} · source: {html.escape(str(metadata.get("source_machine_context", "")))} · baseline: {html.escape(str(metadata.get("baseline_status", "")))} ({html.escape(str(metadata.get("baseline_diagnostic", "")))}) · experimental: {html.escape(str(metadata.get("experimental_status", "")))} ({html.escape(str(metadata.get("experimental_diagnostic", "")))}) · {html.escape(str(metadata.get("parameter_stability", "")))} · HR candidates={html.escape(str(metadata.get("candidate_count", 0)))}</text>',
        f'<text x="42" y="121" class="mono" font-size="11" fill="#667085">max={_compact(metadata.get("max_intensity"))} · bg p99={_compact(metadata.get("background_p99"))} · peak/bg={_compact(metadata.get("peak_to_background_p99"))} · half-height support={_compact(metadata.get("max_half_height_support_scans"))} scans · MS2 cross-stage={html.escape(str(metadata.get("ms2_status") or "not linked"))}</text>',
        f'<text x="42" y="144" class="mono" font-size="10.5" fill="#246bce">{html.escape(str(metadata.get("baseline_label", "")))}</text>',
        f'<text x="42" y="165" class="mono" font-size="10.5" fill="#246bce">{html.escape(str(metadata.get("baseline_metrics", "")))}</text>',
        f'<text x="42" y="187" class="mono" font-size="10.5" fill="#d97706">{html.escape(str(metadata.get("experimental_label", "")))}</text>',
        f'<text x="42" y="208" class="mono" font-size="10.5" fill="#d97706">{html.escape(str(metadata.get("experimental_metrics", "")))}</text>',
        _legend(metadata),
        '<line x1="42" y1="238" x2="1198" y2="238" stroke="#dfe4eb"/>',
        _panel(
            rt_min,
            raw,
            baseline_smooth,
            experimental_smooth,
            full_box,
            min(rt_min, default=0.0),
            max(rt_min, default=1.0),
            "Full chromatogram · linear intensity (zero baseline)",
            manual_rts_min,
            baseline_peaks,
            experimental_peaks,
            review_candidates,
            threshold,
            False,
        ),
        _panel(
            rt_min,
            raw,
            baseline_smooth,
            experimental_smooth,
            focus_box,
            focus[0],
            focus[1],
            f"Peak-region detail · {focus[0]:.2f}–{focus[1]:.2f} min",
            manual_rts_min,
            baseline_peaks,
            experimental_peaks,
            review_candidates,
            threshold,
            False,
        ),
        _panel(
            rt_min,
            raw,
            baseline_smooth,
            experimental_smooth,
            log_box,
            min(rt_min, default=0.0),
            max(rt_min, default=1.0),
            "Full chromatogram · log10(1 + intensity)",
            manual_rts_min,
            baseline_peaks,
            experimental_peaks,
            review_candidates,
            None,
            True,
            compact=True,
        ),
    ]
    return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="MS1 EIC review for {html.escape(compound_id)}">{"".join(body)}</svg>'


def _panel(
    x_values,
    raw,
    baseline,
    experimental,
    box,
    x_min,
    x_max,
    title,
    manual_rts,
    baseline_peaks,
    experimental_peaks,
    review_candidates,
    threshold,
    log_scale,
    *,
    compact=False,
):
    """绘制一个指定RT范围、线性或log强度的EIC面板。"""

    left, top, plot_width, plot_height = box
    indices = [index for index, value in enumerate(x_values) if x_min <= value <= x_max]
    if not indices:
        return ""
    transform = (
        (lambda value: math.log10(1 + max(0.0, value)))
        if log_scale
        else (lambda value: max(0.0, value))
    )
    maximum = max((transform(raw[i]) for i in indices), default=1.0) or 1.0
    maximum = max(maximum, max((transform(baseline[i]) for i in indices), default=0.0))
    x_span = max(x_max - x_min, 1e-9)

    def sx(value):
        """把保留时间映射到SVG面板横坐标。"""

        return left + (value - x_min) / x_span * plot_width

    def sy(value):
        """把变换后的信号强度映射到SVG面板纵坐标。"""

        return top + plot_height - transform(value) / maximum * plot_height

    parts = [
        f'<text x="{left}" y="{top - 10}" font-size="13" font-weight="650" fill="#25334a">{html.escape(title)}</text>',
        f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="#fbfcfe" stroke="#dfe4eb"/>',
    ]
    for fraction in (0.25, 0.5, 0.75):
        y = top + plot_height * (1 - fraction)
        parts.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_width}" y2="{y:.1f}" stroke="#edf0f4"/>'
        )
    if threshold is not None and transform(float(threshold)) <= maximum:
        y = sy(float(threshold))
        parts.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_width}" y2="{y:.1f}" stroke="#bd3d45" stroke-dasharray="5 5"/><text x="{left + plot_width - 4}" y="{y - 4:.1f}" text-anchor="end" font-size="10" fill="#bd3d45">height floor {_compact(threshold)}</text>'
        )
    selected_x = [x_values[i] for i in indices]
    parts.append(
        _polyline(selected_x, [raw[i] for i in indices], sx, sy, "#a8b2c1", 1.0, plot_width)
    )
    parts.append(
        _polyline(selected_x, [baseline[i] for i in indices], sx, sy, "#246bce", 2.0, plot_width)
    )
    parts.append(
        _polyline(
            selected_x,
            [experimental[i] for i in indices],
            sx,
            sy,
            "#d97706",
            1.5,
            plot_width,
            dash="5 3",
        )
    )
    for index, rt in enumerate(manual_rts):
        if x_min <= rt <= x_max:
            x = sx(rt)
            label_y = top + 13 + (index % 2) * 12
            parts.append(
                f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top + plot_height}" stroke="#172033" stroke-dasharray="3 4"/><text x="{x + 3:.1f}" y="{label_y}" font-size="10" fill="#172033">supplied {rt:.3f}</text>'
            )
    for peak in baseline_peaks:
        rt = peak.rt_sec / 60
        if x_min <= rt <= x_max:
            parts.append(
                f'<circle cx="{sx(rt):.1f}" cy="{sy(peak.intensity):.1f}" r="4" fill="#246bce" stroke="#fff"/>'
            )
    for peak in experimental_peaks:
        rt = peak.rt_sec / 60
        if x_min <= rt <= x_max:
            x, y = sx(rt), sy(peak.intensity)
            parts.append(
                f'<path d="M {x:.1f} {y - 5:.1f} L {x - 5:.1f} {y + 4:.1f} L {x + 5:.1f} {y + 4:.1f} Z" fill="#d97706" stroke="#fff"/>'
            )
    for candidate in review_candidates:
        rt = candidate.rt_sec / 60
        if x_min <= rt <= x_max:
            x, y = sx(rt), sy(candidate.intensity)
            parts.append(
                f'<rect x="{x - 3.5:.1f}" y="{y - 3.5:.1f}" width="7" height="7" fill="none" stroke="#7c3aed" stroke-width="1.4"><title>high-recall candidate RT {rt:.4f}; prominence ratio {candidate.prominence_ratio}; support {candidate.support_scans}</title></rect>'
            )
    tick_count = 6 if not compact else 4
    for index in range(tick_count):
        value = x_min + x_span * index / max(tick_count - 1, 1)
        x = sx(value)
        parts.append(
            f'<text x="{x:.1f}" y="{top + plot_height + 17}" text-anchor="middle" class="mono" font-size="10" fill="#667085">{value:.1f}</text>'
        )
    parts.append(
        f'<text x="{left - 8}" y="{top + 11}" text-anchor="end" class="mono" font-size="9" fill="#667085">{_compact(maximum if not log_scale else (10**maximum - 1))}</text>'
    )
    return "".join(parts)


def _polyline(x_values, y_values, sx, sy, color, stroke_width, pixel_width, dash=None):
    """下采样后生成一条SVG polyline，避免大raw文件产生超大图。"""

    points = _downsample_envelope(x_values, y_values, max(240, int(pixel_width)))
    coords = " ".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in points)
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<polyline points="{coords}" fill="none" stroke="{color}" stroke-width="{stroke_width}" vector-effect="non-scaling-stroke"{dash_attr}/>'


def _downsample_envelope(x_values, y_values, max_points):
    """按bin保留局部最低和最高点，缩小数据同时尽量保留窄峰。"""

    if len(x_values) <= max_points:
        return list(zip(x_values, y_values, strict=True))
    bins = max(1, max_points // 2)
    result = []
    for bin_index in range(bins):
        first = int(bin_index * len(x_values) / bins)
        last = max(first + 1, int((bin_index + 1) * len(x_values) / bins))
        chunk = list(zip(x_values[first:last], y_values[first:last], strict=True))
        low = min(range(len(chunk)), key=lambda i: chunk[i][1])
        high = max(range(len(chunk)), key=lambda i: chunk[i][1])
        result.extend(chunk[index] for index in sorted({low, high}))
    return result


def _focus_bounds(rt_min, manual_rts, baseline_peaks, experimental_peaks):
    """优先围绕人工RT确定局部放大范围；缺失时才使用自动峰。"""

    # Supplied RTs define the review region when available.  Mislocalized
    # automatic peaks remain visible in the full panel but must not stretch the
    # detailed panel until the real close pair becomes unreadable.
    anchors = list(manual_rts)
    if not anchors:
        anchors.extend(peak.rt_sec / 60 for peak in baseline_peaks)
        anchors.extend(peak.rt_sec / 60 for peak in experimental_peaks)
    if not anchors and rt_min:
        anchors = [rt_min[0], rt_min[-1]]
    low = max(min(rt_min, default=0.0), min(anchors, default=0.0) - 0.5)
    high = min(max(rt_min, default=1.0), max(anchors, default=1.0) + 0.5)
    if high - low < 1.0:
        center = (low + high) / 2
        low, high = max(min(rt_min), center - 0.5), min(max(rt_min), center + 0.5)
    return low, high


def _legend(metadata):
    """生成MS1 SVG中各曲线、符号和人工RT的图例。"""

    baseline = html.escape(str(metadata.get("baseline_label", "baseline")).split(":", 1)[0])
    experimental = html.escape(
        str(metadata.get("experimental_label", "experimental")).split(":", 1)[0]
    )
    return (
        '<g transform="translate(510 226)"><line x1="0" y1="0" x2="28" y2="0" stroke="#a8b2c1"/><text x="34" y="4" font-size="11" fill="#667085">raw</text>'
        f'<line x1="85" y1="0" x2="113" y2="0" stroke="#246bce" stroke-width="2"/><circle cx="99" cy="0" r="4" fill="#246bce"/><text x="119" y="4" font-size="11" fill="#667085">{baseline}</text>'
        f'<line x1="225" y1="0" x2="253" y2="0" stroke="#d97706" stroke-dasharray="5 3"/><path d="M239 -5 L234 4 L244 4 Z" fill="#d97706"/><text x="259" y="4" font-size="11" fill="#667085">{experimental}</text>'
        '<rect x="365" y="-4" width="8" height="8" fill="none" stroke="#7c3aed"/><text x="379" y="4" font-size="11" fill="#667085">high-recall candidate</text>'
        '<line x1="535" y1="-8" x2="535" y2="8" stroke="#172033" stroke-dasharray="3 4"/><text x="543" y="4" font-size="11" fill="#667085">supplied RT</text></g>'
    )


def _compact(value) -> str:
    """把MS1强度压缩为适合图内显示的k/M文本。"""

    try:
        number = float(value)
    except (TypeError, ValueError):
        return "—"
    if not math.isfinite(number):
        return "—"
    if abs(number) >= 1_000_000:
        return f"{number / 1_000_000:.2f}M"
    if abs(number) >= 1_000:
        return f"{number / 1_000:.1f}k"
    if abs(number) >= 10:
        return f"{number:.1f}"
    return f"{number:.3f}"
