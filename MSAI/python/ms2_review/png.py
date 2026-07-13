"""Direct Pillow renderer for portable MS2 review PNG files."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

from .model import MirrorSpectrumData
from .svg import (
    _acquisition_label,
    _compact,
    _comparison_label,
    _correlation_label,
    _fmt,
    _method_profile_label,
    _number,
)


def render_ms2_review_png(
    *, row: dict, mirror: MirrorSpectrumData, metadata: dict, scale: float = 2.0
) -> Image.Image:
    """使用Pillow渲染与SVG信息一致的便携PNG审核图。"""

    width, height = round(1240 * scale), round(1030 * scale)
    image = Image.new("RGB", (width, height), "#ffffff")
    draw = ImageDraw.Draw(image)

    def s(value):
        """把逻辑坐标缩放并取整为PNG像素坐标。"""

        return round(value * scale)

    def font(size, mono=False, bold=False):
        """按当前输出倍率取得对应字号和字形的字体。"""

        return _font(round(size * scale), mono=mono, bold=bold)

    def put(x, y, value, size=12, color="#172033", mono=False, bold=False):
        """使用MS2审核图的统一样式绘制一段文字。"""

        return draw.text((s(x), s(y)), str(value), fill=color, font=font(size, mono, bold))

    compound = str(row.get("Compound_ID") or metadata.get("target_uid") or "unknown")
    status = str(row.get("ms2_diagnostic_status") or "unknown")
    priority = str(row.get("ms2_review_priority") or "unassigned")
    issues = str(row.get("ms2_issue_codes") or "none")
    put(42, 24, f"{compound} · targeted MS2 evidence", 26, bold=True)
    put(
        42,
        60,
        f"precursor m/z {_fmt(row.get('MZ'), 6)} · supplied RT {_fmt(row.get('Peak1'), 3)} / {_fmt(row.get('Peak2'), 3)} min · separation {_fmt(row.get('rt_separation_sec'), 1)} s",
        14,
        "#4f5c70",
        mono=True,
    )
    put(
        42,
        88,
        f"source: {metadata.get('source_context', '')} · pool {metadata.get('source_pool_id') or '—'} · well {metadata.get('source_pooled_well') or '—'}",
        11,
        "#4f5c70",
    )
    put(
        42,
        118,
        f"MS2 status: {status} · {priority} · availability: {mirror.availability}",
        15,
        bold=True,
    )
    put(42, 146, "issues: " + _clip(issues + _comparison_label(row), 155), 10, "#8a4b08", mono=True)
    put(
        42,
        170,
        _clip(
            f"analysis: standard {metadata.get('standard_id', '')} · filter ≥{_fmt(float(metadata.get('min_relative_intensity', 0)) * 100, 1)}% base peak · match ±{_fmt(metadata.get('fragment_mz_tol'), 3)} {metadata.get('fragment_mz_tol_unit', '')} · {_correlation_label(metadata)}",
            185,
        ),
        8,
        "#667085",
        mono=True,
    )
    put(
        42,
        190,
        _clip(_acquisition_label(metadata.get("acquisition_context") or {}), 185),
        8,
        "#4f5c70",
        mono=True,
    )
    put(
        42,
        208,
        _clip(
            _method_profile_label(
                metadata.get("method_profile") or {},
                metadata.get("acquisition_reconciliation") or {},
            ),
            185,
        ),
        8,
        "#8a4b08",
        mono=True,
    )
    _meter(
        draw,
        s,
        font,
        55,
        242,
        210,
        "cosine",
        row.get("ms2_cosine"),
        threshold=metadata.get("min_cosine"),
    )
    _meter(draw, s, font, 285, 242, 210, "entropy (reported)", row.get("ms2_entropy_similarity"))
    _meter(
        draw,
        s,
        font,
        515,
        242,
        210,
        "explained Peak 1",
        row.get("peak_a_explained_intensity"),
        threshold=metadata.get("min_explained_intensity"),
    )
    _meter(
        draw,
        s,
        font,
        745,
        242,
        210,
        "explained Peak 2",
        row.get("peak_b_explained_intensity"),
        threshold=metadata.get("min_explained_intensity"),
    )
    _meter(
        draw,
        s,
        font,
        975,
        242,
        210,
        "matched fragments",
        row.get("ms2_matched_peaks"),
        threshold=metadata.get("min_matched_peaks"),
        count=True,
    )
    _mirror(draw, s, font, mirror, issues)
    _dia_coverage(draw, s, font, row.get("MZ"), metadata.get("dia_windows") or [])
    _peak_card(draw, s, font, 55, 785, "Peak 1", "peak_a", row, "#246bce")
    _peak_card(draw, s, font, 630, 785, "Peak 2", "peak_b", row, "#d97706")
    put(
        55,
        997,
        f"DIA center {_fmt(row.get('dia_window_center'), 3)} · bounds [{_fmt(row.get('dia_window_lower'), 3)}, {_fmt(row.get('dia_window_upper'), 3)}] · matching windows {_fmt(row.get('dia_window_match_count'), 0)} · config {metadata.get('config_hash', '')}",
        10,
        "#667085",
        mono=True,
    )
    return image


def _meter(draw, s, font, x, y, width, label, value, *, threshold=None, count=False):
    """在PNG画布绘制指标条和阈值线。"""

    number, threshold_number = _number(value), _number(threshold)
    maximum = max(12.0, number or 0.0, threshold_number or 0.0) if count else 1.0
    fraction = 0.0 if number is None else min(1.0, max(0.0, number / maximum))
    draw.text((s(x), s(y)), label, fill="#4f5c70", font=font(11))
    draw.rounded_rectangle(
        (s(x), s(y + 21), s(x + width), s(y + 47)), radius=s(4), fill="#eef1f5", outline="#dfe4eb"
    )
    if fraction > 0:
        draw.rounded_rectangle(
            (s(x), s(y + 21), s(x + width * fraction), s(y + 47)), radius=s(4), fill="#246bce"
        )
    if threshold_number is not None:
        threshold_x = x + width * min(1.0, threshold_number / maximum)
        _dashed_line(draw, s, threshold_x, y + 21, threshold_x, y + 47, "#172033")
    display = "—" if number is None else str(round(number)) if count else f"{number:.3f}"
    _center_text(draw, s, font, x + width / 2, y + 29, display, 11, "#172033", bold=True)


def _mirror(draw, s, font, mirror, issues):
    """在PNG画布绘制双侧镜像fragment谱。"""

    left, right, top, bottom, mid = 70, 1190, 330, 690, 510
    width, amplitude = right - left, 140

    def sx(mz):
        """把fragment m/z映射到镜像谱横坐标。"""

        return left + (mz - mirror.min_mz) / (mirror.max_mz - mirror.min_mz) * width

    draw.text(
        (s(70), s(300)),
        "Normalized mirror spectrum · each side scaled to its own base peak",
        fill="#25334a",
        font=font(14, bold=True),
    )
    draw.rectangle((s(left), s(top), s(right), s(bottom)), fill="#fbfcfe", outline="#dfe4eb")
    draw.line((s(left), s(mid), s(right), s(mid)), fill="#8792a3")
    draw.text((s(78), s(338)), "Peak 1 ↑", fill="#246bce", font=font(12))
    draw.text((s(78), s(656)), "Peak 2 ↓", fill="#d97706", font=font(12))
    if not mirror.peak_a and not mirror.peak_b:
        _center_text(draw, s, font, 630, 480, "No retained fragment spectrum", 17, "#8790a0")
        _center_text(draw, s, font, 630, 512, _clip(issues, 100), 10, "#8790a0", mono=True)
    for index, (mz, intensity) in enumerate(mirror.peak_a):
        x, y = sx(mz), mid - intensity / mirror.max_a * amplitude
        matched = index in mirror.matched_a
        color, line_width = ("#246bce", s(2.3)) if matched else ("#aeb7c5", s(1))
        draw.line((s(x), s(mid), s(x), s(y)), fill=color, width=max(1, line_width))
        if matched:
            draw.ellipse(
                (s(x) - s(2.8), s(y) - s(2.8), s(x) + s(2.8), s(y) + s(2.8)), fill="#246bce"
            )
    for index, (mz, intensity) in enumerate(mirror.peak_b):
        x, y = sx(mz), mid + intensity / mirror.max_b * amplitude
        matched = index in mirror.matched_b
        color, line_width = ("#d97706", s(2.3)) if matched else ("#aeb7c5", s(1))
        draw.line((s(x), s(mid), s(x), s(y)), fill=color, width=max(1, line_width))
        if matched:
            draw.polygon(
                ((s(x), s(y + 3.5)), (s(x - 3.5), s(y - 2.5)), (s(x + 3.5), s(y - 2.5))),
                fill="#d97706",
            )
    for tick in range(7):
        mz = mirror.min_mz + (mirror.max_mz - mirror.min_mz) * tick / 6
        x = sx(mz)
        draw.line((s(x), s(bottom), s(x), s(bottom + 5)), fill="#8792a3")
        _center_text(draw, s, font, x, bottom + 7, f"{mz:.1f}", 10, "#667085", mono=True)
    _right_text(
        draw,
        s,
        font,
        right,
        338,
        "matched: thick + circle/triangle · unmatched: thin grey",
        10,
        "#667085",
    )


def _peak_card(draw, s, font, x, y, title, prefix, row, color):
    """绘制一侧提取质量卡片。"""

    draw.rounded_rectangle(
        (s(x), s(y), s(x + 550), s(y + 190)), radius=s(9), fill="#fbfcfe", outline="#dfe4eb"
    )
    draw.rounded_rectangle((s(x), s(y), s(x + 6), s(y + 190)), radius=s(3), fill=color)
    draw.text((s(x + 18), s(y + 12)), title, fill=color, font=font(15, bold=True))
    supplied = row.get("Peak1" if prefix == "peak_a" else "Peak2")
    lines = [
        f"supplied RT {_fmt(supplied, 3)} min · extraction [{_fmt(row.get(prefix + '_rt_window_start'), 1)}, {_fmt(row.get(prefix + '_rt_window_end'), 1)}] s",
        f"apex RT {_fmt(row.get(prefix + '_apex_rt'), 1)} s · apex intensity {_compact(row.get(prefix + '_apex_intensity'))} · area {_compact(row.get(prefix + '_area'))}",
        f"MS2 scans {_fmt(row.get(prefix + '_ms2_count'), 0)} · consensus scans {_fmt(row.get(prefix + '_consensus_scan_count'), 0)}",
        f"candidate fragments {_fmt(row.get(prefix + '_candidate_fragment_count'), 0)} · retained fragments {_fmt(row.get(prefix + '_fragment_count'), 0)}",
        "quality: " + _clip(str(row.get(prefix + "_quality_flags") or "—"), 82),
    ]
    for index, line in enumerate(lines):
        draw.text(
            (s(x + 18), s(y + 48 + index * 25)), line, fill="#4f5c70", font=font(10, mono=True)
        )


def _dia_coverage(draw, s, font, target_mz, windows):
    """绘制目标m/z相对于真实DIA窗口的位置。"""

    if not windows:
        draw.text(
            (s(70), s(732)),
            "Acquired DIA-window metadata unavailable",
            fill="#8790a0",
            font=font(11),
        )
        return
    lowers = [float(window["lower"]) for window in windows]
    uppers = [float(window["upper"]) for window in windows]
    target = _number(target_mz)
    minimum = min(lowers + ([target - 5] if target is not None else []))
    maximum = max(uppers + ([target + 5] if target is not None else []))
    span = max(maximum - minimum, 1)

    def sx(value):
        """把DIA窗口边界映射到coverage面板横坐标。"""

        return 70 + (value - minimum) / span * 1120

    covered = target is not None and any(
        lower <= target <= upper for lower, upper in zip(lowers, uppers, strict=True)
    )
    draw.text((s(70), s(717)), "Acquired DIA coverage", fill="#25334a", font=font(11, bold=True))
    _right_text(
        draw,
        s,
        font,
        1190,
        717,
        f"{min(lowers):.1f}–{max(uppers):.1f} m/z · n={len(windows)}",
        10,
        "#667085",
    )
    draw.rectangle((s(70), s(741), s(1190), s(759)), fill="#f3f5f8", outline="#dfe4eb")
    for index, (lower, upper) in enumerate(zip(lowers, uppers, strict=True)):
        draw.rectangle(
            (s(sx(lower)), s(741), s(sx(upper)), s(759)),
            fill="#dce9fb" if index % 2 == 0 else "#eef3fb",
            outline="#8aaee3",
        )
    if target is not None:
        x = sx(target)
        draw.line((s(x), s(736), s(x), s(764)), fill="#d97706", width=max(1, s(2)))
        label = "inside acquired window" if covered else "outside acquired windows"
        draw.text(
            (s(min(max(x + 5, 75), 1040)), s(765)),
            f"target {target:.4f} · {label}",
            fill="#8a4b08",
            font=font(10, mono=True),
        )


def _dashed_line(draw, s, x1, y1, x2, y2, color):
    """在MS2 PNG画布上绘制虚线。"""

    segments = 8
    for index in range(0, segments, 2):
        a, b = index / segments, (index + 1) / segments
        draw.line(
            (
                s(x1 + (x2 - x1) * a),
                s(y1 + (y2 - y1) * a),
                s(x1 + (x2 - x1) * b),
                s(y1 + (y2 - y1) * b),
            ),
            fill=color,
        )


def _center_text(draw, s, font, x, y, value, size, color, *, mono=False, bold=False):
    """居中绘制MS2审核图文字。"""

    active = font(size, mono, bold)
    box = draw.textbbox((0, 0), value, font=active)
    draw.text((s(x) - (box[2] - box[0]) // 2, s(y)), value, fill=color, font=active)


def _right_text(draw, s, font, x, y, value, size, color):
    """右对齐绘制MS2审核图文字。"""

    active = font(size)
    box = draw.textbbox((0, 0), value, font=active)
    draw.text((s(x) - (box[2] - box[0]), s(y)), value, fill=color, font=active)


def _clip(value: str, limit: int) -> str:
    """将过长图内文本截断并添加省略号。"""

    return value if len(value) <= limit else value[: limit - 1] + "…"


@lru_cache(maxsize=16)
def _font(size: int, mono=False, bold=False):
    """缓存所需Windows字体并在缺失时安全回退。"""

    name = "consola.ttf" if mono else "segoeuib.ttf" if bold else "segoeui.ttf"
    try:
        return ImageFont.truetype(str(Path("C:/Windows/Fonts") / name), max(1, size))
    except OSError:
        return ImageFont.load_default()
