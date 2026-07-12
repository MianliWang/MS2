"""Dependency-free SVG renderer for one MS2 pair-evidence review image."""

from __future__ import annotations

import html
import math

from .model import MirrorSpectrumData


def render_ms2_review_svg(*, row: dict, mirror: MirrorSpectrumData, metadata: dict) -> str:
    width, height = 1240, 1030
    compound = str(row.get("Compound_ID") or metadata.get("target_uid") or "unknown")
    status = str(row.get("ms2_diagnostic_status") or "unknown")
    priority = str(row.get("ms2_review_priority") or "unassigned")
    issues = str(row.get("ms2_issue_codes") or "none")
    comparison = _comparison_label(row)
    mz = _number(row.get("MZ"))
    peak1, peak2 = _number(row.get("Peak1")), _number(row.get("Peak2"))
    body = [
        '<rect width="1240" height="1030" fill="#ffffff"/>',
        '<style>text{font-family:Segoe UI,Arial,sans-serif}.mono{font-family:Consolas,monospace}</style>',
        f'<text x="42" y="45" font-size="26" font-weight="700" fill="#172033">{html.escape(compound)} · targeted MS2 evidence</text>',
        f'<text x="42" y="75" class="mono" font-size="14" fill="#4f5c70">precursor m/z {_fmt(mz, 6)} · supplied RT {_fmt(peak1, 3)} / {_fmt(peak2, 3)} min · separation {_fmt(row.get("rt_separation_sec"), 1)} s</text>',
        f'<text x="42" y="102" font-size="12" fill="#4f5c70">source: {html.escape(str(metadata.get("source_context", "")))} · pool {html.escape(str(metadata.get("source_pool_id") or "—"))} · well {html.escape(str(metadata.get("source_pooled_well") or "—"))}</text>',
        f'<text x="42" y="132" font-size="15" font-weight="650" fill="#172033">MS2 status: {html.escape(status)} · {html.escape(priority)} · availability: {html.escape(mirror.availability)}</text>',
        f'<text x="42" y="158" class="mono" font-size="11" fill="#8a4b08">issues: {html.escape(issues)}{html.escape(comparison)}</text>',
        f'<text x="42" y="181" class="mono" font-size="9.5" fill="#667085">analysis: standard {html.escape(str(metadata.get("standard_id", "")))} · filter ≥{_fmt(float(metadata.get("min_relative_intensity", 0))*100, 1)}% base peak · match ±{_fmt(metadata.get("fragment_mz_tol"), 3)} {html.escape(str(metadata.get("fragment_mz_tol_unit", "")))} · {html.escape(_correlation_label(metadata))}</text>',
        f'<text x="42" y="201" class="mono" font-size="9.5" fill="#4f5c70">{html.escape(_acquisition_label(metadata.get("acquisition_context") or {}))}</text>',
        f'<text x="42" y="219" class="mono" font-size="9.5" fill="#8a4b08">{html.escape(_method_profile_label(metadata.get("method_profile") or {}, metadata.get("acquisition_reconciliation") or {}))}</text>',
        _meter_svg(55, 242, 210, "cosine", row.get("ms2_cosine"), threshold=metadata.get("min_cosine")),
        _meter_svg(285, 242, 210, "entropy (reported)", row.get("ms2_entropy_similarity")),
        _meter_svg(515, 242, 210, "explained Peak 1", row.get("peak_a_explained_intensity"), threshold=metadata.get("min_explained_intensity")),
        _meter_svg(745, 242, 210, "explained Peak 2", row.get("peak_b_explained_intensity"), threshold=metadata.get("min_explained_intensity")),
        _meter_svg(975, 242, 210, "matched fragments", row.get("ms2_matched_peaks"), threshold=metadata.get("min_matched_peaks"), count=True),
        _mirror_svg(mirror, issues),
        _dia_coverage_svg(mz, metadata.get("dia_windows") or []),
        _peak_card_svg(55, 785, "Peak 1", "peak_a", row, "#246bce"),
        _peak_card_svg(630, 785, "Peak 2", "peak_b", row, "#d97706"),
        f'<text x="55" y="1008" class="mono" font-size="10.5" fill="#667085">DIA center {_fmt(row.get("dia_window_center"), 3)} · bounds [{_fmt(row.get("dia_window_lower"), 3)}, {_fmt(row.get("dia_window_upper"), 3)}] · matching windows {_fmt(row.get("dia_window_match_count"), 0)} · config {html.escape(str(metadata.get("config_hash", "")))}</text>',
    ]
    return f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="MS2 evidence review for {html.escape(compound)}">{"".join(body)}</svg>'


def _meter_svg(x, y, width, label, value, *, threshold=None, count=False):
    number = _number(value)
    threshold_number = _number(threshold)
    maximum = max(12.0, number or 0.0, threshold_number or 0.0) if count else 1.0
    fraction = 0.0 if number is None else min(1.0, max(0.0, number / maximum))
    threshold_fraction = None if threshold_number is None else min(1.0, threshold_number / maximum)
    threshold_mark = "" if threshold_fraction is None else (
        f'<line x1="{x+width*threshold_fraction:.1f}" y1="{y+21}" x2="{x+width*threshold_fraction:.1f}" y2="{y+47}" stroke="#172033" stroke-dasharray="3 2"/>'
    )
    display = "—" if number is None else str(int(round(number))) if count else f"{number:.3f}"
    return (
        f'<text x="{x}" y="{y}" font-size="11" fill="#4f5c70">{html.escape(label)}</text>'
        f'<rect x="{x}" y="{y+21}" width="{width}" height="26" rx="4" fill="#eef1f5" stroke="#dfe4eb"/>'
        f'<rect x="{x}" y="{y+21}" width="{width*fraction:.1f}" height="26" rx="4" fill="#246bce"/>{threshold_mark}'
        f'<text x="{x+width/2}" y="{y+39}" text-anchor="middle" font-size="11" font-weight="650" fill="#172033">{display}</text>'
    )


def _mirror_svg(mirror: MirrorSpectrumData, issues: str) -> str:
    left, right, top, bottom, mid = 70, 1190, 330, 690, 510
    plot_width, amplitude = right - left, 140
    scale_x = lambda mz: left + (mz - mirror.min_mz) / (mirror.max_mz - mirror.min_mz) * plot_width
    marks = [
        '<text x="70" y="316" font-size="14" font-weight="650" fill="#25334a">Normalized mirror spectrum · each side scaled to its own base peak</text>',
        f'<rect x="{left}" y="{top}" width="{plot_width}" height="{bottom-top}" fill="#fbfcfe" stroke="#dfe4eb"/>',
        f'<line x1="{left}" y1="{mid}" x2="{right}" y2="{mid}" stroke="#8792a3"/>',
        '<text x="78" y="351" font-size="12" fill="#246bce">Peak 1 ↑</text>',
        '<text x="78" y="675" font-size="12" fill="#d97706">Peak 2 ↓</text>',
    ]
    if not mirror.peak_a and not mirror.peak_b:
        marks.append(f'<text x="630" y="495" text-anchor="middle" font-size="17" fill="#8790a0">No retained fragment spectrum</text><text x="630" y="525" text-anchor="middle" font-size="11" fill="#8790a0">{html.escape(issues)}</text>')
    for index, (mz, intensity) in enumerate(mirror.peak_a):
        x, y = scale_x(mz), mid - intensity / mirror.max_a * amplitude
        matched = index in mirror.matched_a
        color, stroke = ("#246bce", 2.3) if matched else ("#aeb7c5", 1.0)
        marks.append(f'<line x1="{x:.2f}" y1="{mid}" x2="{x:.2f}" y2="{y:.2f}" stroke="{color}" stroke-width="{stroke}"><title>Peak 1 m/z {mz:.5f}; intensity {intensity:.6g}; {"matched" if matched else "unmatched"}</title></line>')
        if matched:
            marks.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.8" fill="#246bce"/>')
    for index, (mz, intensity) in enumerate(mirror.peak_b):
        x, y = scale_x(mz), mid + intensity / mirror.max_b * amplitude
        matched = index in mirror.matched_b
        color, stroke = ("#d97706", 2.3) if matched else ("#aeb7c5", 1.0)
        marks.append(f'<line x1="{x:.2f}" y1="{mid}" x2="{x:.2f}" y2="{y:.2f}" stroke="{color}" stroke-width="{stroke}"><title>Peak 2 m/z {mz:.5f}; intensity {intensity:.6g}; {"matched" if matched else "unmatched"}</title></line>')
        if matched:
            marks.append(f'<path d="M {x:.2f} {y+3.5:.2f} L {x-3.5:.2f} {y-2.5:.2f} L {x+3.5:.2f} {y-2.5:.2f} Z" fill="#d97706"/>')
    for tick in range(7):
        mz = mirror.min_mz + (mirror.max_mz - mirror.min_mz) * tick / 6
        x = scale_x(mz)
        marks.append(f'<line x1="{x:.1f}" y1="{bottom}" x2="{x:.1f}" y2="{bottom+5}" stroke="#8792a3"/><text x="{x:.1f}" y="{bottom+20}" text-anchor="middle" class="mono" font-size="10" fill="#667085">{mz:.1f}</text>')
    marks.append(f'<text x="{right}" y="350" text-anchor="end" font-size="10" fill="#667085">matched: thick + circle/triangle · unmatched: thin grey</text>')
    return "".join(marks)


def _peak_card_svg(x, y, title, prefix, row, color):
    width, height = 550, 190
    lines = [
        f"supplied RT {_fmt(row.get('Peak1' if prefix == 'peak_a' else 'Peak2'), 3)} min · extraction [{_fmt(row.get(prefix+'_rt_window_start'), 1)}, {_fmt(row.get(prefix+'_rt_window_end'), 1)}] s",
        f"apex RT {_fmt(row.get(prefix+'_apex_rt'), 1)} s · apex intensity {_compact(row.get(prefix+'_apex_intensity'))} · area {_compact(row.get(prefix+'_area'))}",
        f"MS2 scans {_fmt(row.get(prefix+'_ms2_count'), 0)} · consensus scans {_fmt(row.get(prefix+'_consensus_scan_count'), 0)}",
        f"candidate fragments {_fmt(row.get(prefix+'_candidate_fragment_count'), 0)} · retained fragments {_fmt(row.get(prefix+'_fragment_count'), 0)}",
        "quality: " + str(row.get(prefix + "_quality_flags") or "—")[:88],
    ]
    text = "".join(
        f'<text x="{x+18}" y="{y+58+index*25}" class="mono" font-size="10.5" fill="#4f5c70">{html.escape(line)}</text>'
        for index, line in enumerate(lines)
    )
    return (
        f'<rect x="{x}" y="{y}" width="{width}" height="{height}" rx="9" fill="#fbfcfe" stroke="#dfe4eb"/>'
        f'<rect x="{x}" y="{y}" width="6" height="{height}" rx="3" fill="{color}"/>'
        f'<text x="{x+18}" y="{y+30}" font-size="15" font-weight="650" fill="{color}">{title}</text>{text}'
    )


def _dia_coverage_svg(target_mz, windows):
    if not windows:
        return '<text x="70" y="748" font-size="11" fill="#8790a0">Acquired DIA-window metadata unavailable</text>'
    lowers = [float(window["lower"]) for window in windows]
    uppers = [float(window["upper"]) for window in windows]
    target = _number(target_mz)
    minimum = min(lowers + ([target - 5] if target is not None else []))
    maximum = max(uppers + ([target + 5] if target is not None else []))
    span = max(maximum - minimum, 1)
    sx = lambda value: 70 + (value - minimum) / span * 1120
    covered = target is not None and any(lower <= target <= upper for lower, upper in zip(lowers, uppers))
    marks = [
        '<text x="70" y="732" font-size="11" font-weight="650" fill="#25334a">Acquired DIA coverage</text>',
        '<rect x="70" y="741" width="1120" height="18" fill="#f3f5f8" stroke="#dfe4eb"/>',
    ]
    for index, (lower, upper) in enumerate(zip(lowers, uppers)):
        marks.append(f'<rect x="{sx(lower):.2f}" y="741" width="{max(1.0, sx(upper)-sx(lower)):.2f}" height="18" fill="{"#dce9fb" if index%2 == 0 else "#eef3fb"}" stroke="#8aaee3" stroke-width="0.6"/>')
    if target is not None:
        x = sx(target)
        label = "inside acquired window" if covered else "outside acquired windows"
        marks.append(f'<line x1="{x:.2f}" y1="736" x2="{x:.2f}" y2="764" stroke="#d97706" stroke-width="2"/><text x="{min(max(x+5, 75), 1040):.2f}" y="778" class="mono" font-size="10" fill="#8a4b08">target {target:.4f} · {label}</text>')
    marks.append(f'<text x="1190" y="732" text-anchor="end" class="mono" font-size="10" fill="#667085">{min(lowers):.1f}–{max(uppers):.1f} m/z · n={len(windows)}</text>')
    return "".join(marks)


def _number(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _fmt(value, digits=3):
    number = _number(value)
    if number is None:
        return "—"
    return str(int(round(number))) if digits == 0 else f"{number:.{digits}f}"


def _compact(value):
    number = _number(value)
    if number is None:
        return "—"
    if abs(number) >= 1_000_000:
        return f"{number/1_000_000:.2f}M"
    if abs(number) >= 1_000:
        return f"{number/1_000:.1f}k"
    return f"{number:.3g}"


def _acquisition_label(context: dict) -> str:
    if not context:
        return "raw observed: acquisition metadata unavailable"
    return (
        f"raw observed ({context.get('scope', 'available portion')}): "
        f"activation={_observed_summary(context, 'activationMethod')} · "
        f"collisionEnergy field={_observed_summary(context, 'collisionEnergy')} · "
        f"windowWideness={_observed_summary(context, 'windowWideness')}"
    )


def _method_profile_label(profile: dict, reconciliation: dict) -> str:
    if not profile:
        return "reference method: not attached"
    dia = ((profile.get("acquisition") or {}).get("dia") or {})
    processing = profile.get("processing_reference") or {}
    hcd = "/".join(_fmt(value, 0) for value in dia.get("stepped_hcd_percent", [])) or "—"
    flags = len(reconciliation.get("issue_codes") or [])
    return (
        f"reference {profile.get('profile_id', 'method')} (not run-confirmed): "
        f"DIA {_fmt(dia.get('nominal_isolation_width_mz'), 1)} m/z; preprocessing {_fmt(dia.get('preprocessing_fixed_isolation_width_da'), 1)} Da · "
        f"stepped HCD {hcd}% · {_fmt(dia.get('files_per_sample'), 0)} files/sample · "
        f"RT ±{_fmt(processing.get('precursor_fragment_rt_alignment_sec'), 0)} s · r>{_fmt(processing.get('min_chromatographic_correlation'), 2)} · provenance flags={flags}"
    )


def _correlation_label(metadata: dict) -> str:
    mode = str(metadata.get("fragment_correlation_mode") or "full_window")
    base = (
        f"RT half-window≤{_fmt(metadata.get('rt_half_window_sec'), 1)}s · "
        f"coelution r>{_fmt(metadata.get('min_fragment_correlation'), 2)} ({mode})"
    )
    if mode != "active_support":
        return base
    return (
        base
        + f"; precursor≥{_fmt(float(metadata.get('correlation_min_relative_intensity', 0.05))*100, 1)}% apex"
        + f", scans≥{_fmt(metadata.get('min_correlation_scans'), 0)}"
        + f", apexΔ≤{_fmt(metadata.get('max_fragment_apex_offset_scans'), 0)}"
        + f", consecutive≥{_fmt(metadata.get('min_consecutive_fragment_scans'), 0)}"
    )


def _observed_summary(context: dict, key: str) -> str:
    values = (context.get("observed_values") or {}).get(key) or []
    if values:
        return "/".join(
            f"{item.get('value', '—')} (n={item.get('count', '—')})" for item in values
        )
    return str(context.get(key, "—"))


def _comparison_label(row: dict) -> str:
    baseline = str(row.get("comparison_baseline_status") or "")
    candidate = str(row.get("comparison_candidate_status") or "")
    if not baseline and not candidate:
        return ""
    changed = "/".join(
        name
        for name, field in (
            ("P1", "comparison_peak_a_spectrum_changed"),
            ("P2", "comparison_peak_b_spectrum_changed"),
        )
        if str(row.get(field) or "").lower() == "true"
    ) or "none"
    tier = str(row.get("comparison_review_tier") or "")
    tier_label = f"; tier={tier}" if tier else ""
    return f" · shadow: {baseline or '—'} → {candidate or '—'}; spectra={changed}{tier_label}"
