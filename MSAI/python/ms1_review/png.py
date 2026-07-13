"""Pillow renderer for portable PNG copies of MS1 EIC review images."""

from __future__ import annotations

import math
from functools import lru_cache
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

from .svg import _compact, _downsample_envelope, _focus_bounds


def render_eic_png(
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
    scale: float = 2.0,
) -> Image.Image:
    """用Pillow绘制与SVG信息一致的MS1审核PNG，无需浏览器或Cairo。"""

    width, height = round(1240 * scale), round(930 * scale)
    image = Image.new("RGB", (width, height), "#ffffff")
    draw = ImageDraw.Draw(image)

    def coordinate(value):
        """把逻辑坐标缩放并取整为PNG像素坐标。"""

        return round(value * scale)

    def font(size, *, mono=False, bold=False):
        """按当前输出倍率取得对应字号和字形的字体。"""

        return _font(round(size * scale), mono=mono, bold=bold)

    def text(x, y, value, size=12, color="#172033", mono=False, bold=False):
        """使用统一坐标、字体和颜色约定绘制一段文字。"""

        return draw.text(
            (coordinate(x), coordinate(y)),
            str(value),
            fill=color,
            font=font(size, mono=mono, bold=bold),
        )

    rt_min = [value / 60 for value in rts_sec]
    focus = _focus_bounds(rt_min, manual_rts_min, baseline_peaks, experimental_peaks)
    threshold = metadata.get("min_height")
    text(42, 25, f"{compound_id} · targeted MS1 EIC", 25, bold=True)
    text(
        42,
        58,
        f"m/z {target_mz:.6f} · {ppm:g} ppm [{target_mz * (1 - ppm / 1e6):.6f}, {target_mz * (1 + ppm / 1e6):.6f}]",
        14,
        "#4f5c70",
        mono=True,
    )
    text(
        42,
        85,
        "supplied RT reference: "
        + str(metadata.get("reference_status", ""))
        + " · source: "
        + str(metadata.get("source_machine_context", ""))
        + " · baseline: "
        + str(metadata.get("baseline_status", ""))
        + " ("
        + str(metadata.get("baseline_diagnostic", ""))
        + ")"
        + " · experimental: "
        + str(metadata.get("experimental_status", ""))
        + " ("
        + str(metadata.get("experimental_diagnostic", ""))
        + ")"
        + " · "
        + str(metadata.get("parameter_stability", ""))
        + " · HR candidates="
        + str(metadata.get("candidate_count", 0)),
        11,
        "#25334a",
    )
    text(
        42,
        110,
        f"max={_compact(metadata.get('max_intensity'))} · bg p99={_compact(metadata.get('background_p99'))} · "
        f"peak/bg={_compact(metadata.get('peak_to_background_p99'))} · half-height support={_compact(metadata.get('max_half_height_support_scans'))} scans · "
        f"MS2 cross-stage={metadata.get('ms2_status') or 'not linked'}",
        10,
        "#667085",
        mono=True,
    )
    text(42, 136, metadata.get("baseline_label", ""), 10, "#246bce", mono=True)
    text(42, 157, metadata.get("baseline_metrics", ""), 10, "#246bce", mono=True)
    text(42, 179, metadata.get("experimental_label", ""), 10, "#d97706", mono=True)
    text(42, 200, metadata.get("experimental_metrics", ""), 10, "#d97706", mono=True)
    _draw_legend(draw, coordinate, font)
    draw.line((coordinate(42), coordinate(238), coordinate(1198), coordinate(238)), fill="#dfe4eb")
    _draw_panel(
        draw,
        coordinate,
        font,
        rt_min,
        raw,
        baseline_smooth,
        experimental_smooth,
        (78, 260, 1120, 240),
        min(rt_min, default=0.0),
        max(rt_min, default=1.0),
        "Full chromatogram · linear intensity (zero baseline)",
        manual_rts_min,
        baseline_peaks,
        experimental_peaks,
        review_candidates,
        threshold,
        False,
    )
    _draw_panel(
        draw,
        coordinate,
        font,
        rt_min,
        raw,
        baseline_smooth,
        experimental_smooth,
        (78, 550, 1120, 210),
        focus[0],
        focus[1],
        f"Peak-region detail · {focus[0]:.2f}–{focus[1]:.2f} min",
        manual_rts_min,
        baseline_peaks,
        experimental_peaks,
        review_candidates,
        threshold,
        False,
    )
    _draw_panel(
        draw,
        coordinate,
        font,
        rt_min,
        raw,
        baseline_smooth,
        experimental_smooth,
        (78, 810, 1120, 80),
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
    )
    return image


def _draw_panel(
    draw,
    coordinate,
    font,
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
    """在PNG画布绘制一个EIC面板、阈值、RT标记和候选峰。"""

    left, top, plot_width, plot_height = box
    indices = [index for index, value in enumerate(x_values) if x_min <= value <= x_max]
    if not indices:
        return
    transform = (
        (lambda value: math.log10(1 + max(0.0, value)))
        if log_scale
        else (lambda value: max(0.0, value))
    )
    maximum = max((transform(raw[index]) for index in indices), default=1.0) or 1.0
    maximum = max(maximum, max((transform(baseline[index]) for index in indices), default=0.0))
    span = max(x_max - x_min, 1e-9)

    def sx(value):
        """把保留时间映射到面板横坐标。"""

        return left + (value - x_min) / span * plot_width

    def sy(value):
        """把变换后的信号强度映射到面板纵坐标。"""

        return top + plot_height - transform(value) / maximum * plot_height

    draw.text(
        (coordinate(left), coordinate(top - 20)), title, fill="#25334a", font=font(13, bold=True)
    )
    draw.rectangle(
        (
            coordinate(left),
            coordinate(top),
            coordinate(left + plot_width),
            coordinate(top + plot_height),
        ),
        fill="#fbfcfe",
        outline="#dfe4eb",
    )
    for fraction in (0.25, 0.5, 0.75):
        y = top + plot_height * (1 - fraction)
        draw.line(
            (coordinate(left), coordinate(y), coordinate(left + plot_width), coordinate(y)),
            fill="#edf0f4",
        )
    if threshold is not None and transform(float(threshold)) <= maximum:
        y = sy(float(threshold))
        _dashed_line(draw, coordinate, (left, y), (left + plot_width, y), "#bd3d45")
        _right_text(
            draw,
            coordinate,
            font,
            left + plot_width - 4,
            y - 15,
            f"height floor {_compact(threshold)}",
            10,
            "#bd3d45",
        )
    selected_x = [x_values[index] for index in indices]
    _draw_trace(
        draw,
        coordinate,
        selected_x,
        [raw[index] for index in indices],
        sx,
        sy,
        "#a8b2c1",
        plot_width,
    )
    _draw_trace(
        draw,
        coordinate,
        selected_x,
        [baseline[index] for index in indices],
        sx,
        sy,
        "#246bce",
        plot_width,
        width=2,
    )
    _draw_trace(
        draw,
        coordinate,
        selected_x,
        [experimental[index] for index in indices],
        sx,
        sy,
        "#d97706",
        plot_width,
        dashed=True,
    )
    for index, rt in enumerate(manual_rts):
        if x_min <= rt <= x_max:
            x = sx(rt)
            _dashed_line(draw, coordinate, (x, top), (x, top + plot_height), "#172033", dash=3)
            draw.text(
                (coordinate(x + 3), coordinate(top + 3 + (index % 2) * 12)),
                f"supplied {rt:.3f}",
                fill="#172033",
                font=font(10),
            )
    for peak in baseline_peaks:
        rt = peak.rt_sec / 60
        if x_min <= rt <= x_max:
            x, y = coordinate(sx(rt)), coordinate(sy(peak.intensity))
            draw.ellipse(
                (x - coordinate(4), y - coordinate(4), x + coordinate(4), y + coordinate(4)),
                fill="#246bce",
                outline="#ffffff",
            )
    for peak in experimental_peaks:
        rt = peak.rt_sec / 60
        if x_min <= rt <= x_max:
            x, y = coordinate(sx(rt)), coordinate(sy(peak.intensity))
            draw.polygon(
                (
                    (x, y - coordinate(5)),
                    (x - coordinate(5), y + coordinate(4)),
                    (x + coordinate(5), y + coordinate(4)),
                ),
                fill="#d97706",
                outline="#ffffff",
            )
    for candidate in review_candidates:
        rt = candidate.rt_sec / 60
        if x_min <= rt <= x_max:
            x, y = coordinate(sx(rt)), coordinate(sy(candidate.intensity))
            draw.rectangle(
                (
                    x - coordinate(3.5),
                    y - coordinate(3.5),
                    x + coordinate(3.5),
                    y + coordinate(3.5),
                ),
                outline="#7c3aed",
                width=max(1, coordinate(1.4)),
            )
    tick_count = 4 if compact else 6
    for index in range(tick_count):
        value = x_min + span * index / max(tick_count - 1, 1)
        _center_text(
            draw,
            coordinate,
            font,
            sx(value),
            top + plot_height + 6,
            f"{value:.1f}",
            10,
            "#667085",
            mono=True,
        )
    _right_text(
        draw,
        coordinate,
        font,
        left - 8,
        top + 2,
        _compact(maximum if not log_scale else (10**maximum - 1)),
        9,
        "#667085",
        mono=True,
    )


def _draw_trace(
    draw, coordinate, x_values, y_values, sx, sy, color, pixel_width, *, width=1, dashed=False
):
    """经包络下采样后绘制实线或虚线色谱轨迹。"""

    points = _downsample_envelope(x_values, y_values, max(240, int(pixel_width)))
    scaled = [(coordinate(sx(x)), coordinate(sy(y))) for x, y in points]
    if len(scaled) < 2:
        return
    if not dashed:
        draw.line(scaled, fill=color, width=max(1, coordinate(width)), joint="curve")
        return
    for start in range(0, len(scaled) - 1, 4):
        draw.line(
            scaled[start : min(start + 3, len(scaled))],
            fill=color,
            width=max(1, coordinate(width)),
            joint="curve",
        )


def _dashed_line(draw, coordinate, start, end, color, *, dash=5):
    """用短线段在Pillow画布模拟虚线。"""

    x1, y1 = start
    x2, y2 = end
    length = math.hypot(x2 - x1, y2 - y1)
    if not length:
        return
    step = dash * 2
    offset = 0.0
    while offset < length:
        next_offset = min(offset + dash, length)
        ratio_a, ratio_b = offset / length, next_offset / length
        draw.line(
            (
                coordinate(x1 + (x2 - x1) * ratio_a),
                coordinate(y1 + (y2 - y1) * ratio_a),
                coordinate(x1 + (x2 - x1) * ratio_b),
                coordinate(y1 + (y2 - y1) * ratio_b),
            ),
            fill=color,
        )
        offset += step


def _draw_legend(draw, coordinate, font):
    """绘制raw、两组参数、候选峰和人工RT图例。"""

    y = coordinate(226)
    draw.line((coordinate(510), y, coordinate(538), y), fill="#a8b2c1")
    draw.text((coordinate(544), coordinate(220)), "raw", fill="#667085", font=font(11))
    draw.line((coordinate(595), y, coordinate(623), y), fill="#246bce", width=coordinate(2))
    draw.ellipse(
        (coordinate(605), y - coordinate(4), coordinate(613), y + coordinate(4)), fill="#246bce"
    )
    draw.text((coordinate(629), coordinate(220)), "baseline", fill="#667085", font=font(11))
    _dashed_line(draw, coordinate, (735, 226), (763, 226), "#d97706")
    draw.polygon(
        (
            (coordinate(749), y - coordinate(5)),
            (coordinate(744), y + coordinate(4)),
            (coordinate(754), y + coordinate(4)),
        ),
        fill="#d97706",
    )
    draw.text((coordinate(769), coordinate(220)), "experimental", fill="#667085", font=font(11))
    draw.rectangle(
        (coordinate(875), y - coordinate(4), coordinate(883), y + coordinate(4)), outline="#7c3aed"
    )
    draw.text(
        (coordinate(889), coordinate(220)), "high-recall candidate", fill="#667085", font=font(11)
    )
    _dashed_line(draw, coordinate, (1045, 218), (1045, 234), "#172033", dash=3)
    draw.text((coordinate(1053), coordinate(220)), "supplied RT", fill="#667085", font=font(11))


def _center_text(draw, coordinate, font, x, y, value, size, color, *, mono=False):
    """按文本包围盒在指定x位置居中绘制文字。"""

    box = draw.textbbox((0, 0), value, font=font(size, mono=mono))
    draw.text(
        (coordinate(x) - (box[2] - box[0]) // 2, coordinate(y)),
        value,
        fill=color,
        font=font(size, mono=mono),
    )


def _right_text(draw, coordinate, font, x, y, value, size, color, *, mono=False):
    """按文本包围盒在指定x位置右对齐绘制文字。"""

    box = draw.textbbox((0, 0), value, font=font(size, mono=mono))
    draw.text(
        (coordinate(x) - (box[2] - box[0]), coordinate(y)),
        value,
        fill=color,
        font=font(size, mono=mono),
    )


@lru_cache(maxsize=16)
def _font(size: int, *, mono=False, bold=False):
    """缓存Windows字体；字体不可用时回退Pillow默认字体。"""

    name = "consola.ttf" if mono else "segoeuib.ttf" if bold else "segoeui.ttf"
    path = Path("C:/Windows/Fonts") / name
    try:
        return ImageFont.truetype(str(path), max(1, size))
    except OSError:
        return ImageFont.load_default()
