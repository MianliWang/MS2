"""Batch exporter for canonical MS1 EIC SVGs and review-folder views."""

from __future__ import annotations

import csv
import hashlib
import html
import json
import os
import shutil
import statistics
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

try:
    from ..ms1_peak_picker import (
        PeakPickingConfig,
        enumerate_peak_candidates,
        extract_target_eics_with_config,
        load_peak_config,
        pick_chiral_peaks,
        savgol_smooth,
    )
    from ..ms2 import number, read_table
except ImportError:
    from ms1_peak_picker import (  # type: ignore
        PeakPickingConfig,
        enumerate_peak_candidates,
        extract_target_eics_with_config,
        load_peak_config,
        pick_chiral_peaks,
        savgol_smooth,
    )
    from ms2 import number, read_table  # type: ignore
from .classification import (
    background_diagnostics,
    prediction_diagnostic,
    reference_chromatographic_status,
    review_views,
    source_machine_label_interpretation,
)
from .svg import render_eic_svg
from .png import render_eic_png


def export_ms1_review(
    peaklist_path,
    raw_path,
    output_dir,
    *,
    baseline_config_path=None,
    experimental_config_path=None,
    candidate_config_path=None,
    ms2_result_path=None,
    mz_column: str = "MZ",
    source_machine_label_column: str | None = "IG",
    rt_tolerance_min: float = 0.1,
    image_formats: tuple[str, ...] = ("svg", "png"),
    png_scale: float = 2.0,
) -> dict:
    peaklist_path, raw_path, output_dir = map(Path, (peaklist_path, raw_path, output_dir))
    image_formats = _normalise_image_formats(image_formats)
    if png_scale <= 0:
        raise ValueError("png_scale must be positive")
    baseline_config = load_peak_config(baseline_config_path) if baseline_config_path else PeakPickingConfig()
    experimental_config = (
        load_peak_config(experimental_config_path) if experimental_config_path else baseline_config
    )
    candidate_config = (
        load_peak_config(candidate_config_path) if candidate_config_path else experimental_config
    )
    rows = read_table(peaklist_path)
    targets = [value for row in rows if (value := number(row.get(mz_column))) is not None]
    rts, traces = extract_target_eics_with_config(raw_path, targets, baseline_config)
    ms2_by_compound = _ms2_lookup(ms2_result_path)
    _prepare_output(output_dir)
    config_hash = _config_hash(baseline_config, experimental_config, candidate_config)
    run_id = _safe_name(raw_path.stem)
    dataset_id = _safe_name(peaklist_path.stem)
    images_dir = output_dir / "assets" / "ms1" / f"cfg-{config_hash}" / run_id
    views_dir = output_dir / "views"
    metadata_dir = output_dir / "metadata"
    asset_dirs = {image_format: images_dir / image_format for image_format in image_formats}
    for asset_dir in asset_dirs.values():
        asset_dir.mkdir(parents=True, exist_ok=True)
    views_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    manifest: list[dict] = []
    view_index: list[dict] = []
    view_files: dict[str, list[dict[str, str]]] = {}
    median_scan_interval = statistics.median(
        right - left for left, right in zip(rts, rts[1:])
    ) if len(rts) > 1 else 0.0
    for row_index, row in enumerate(rows, start=1):
        mz = number(row.get(mz_column))
        if mz is None or mz not in traces:
            continue
        raw = [float(value) for value in traces[mz]]
        baseline_result = pick_chiral_peaks(rts, raw, baseline_config)
        experimental_result = pick_chiral_peaks(rts, raw, experimental_config)
        review_candidates = enumerate_peak_candidates(rts, raw, candidate_config)
        baseline_smooth = _smooth(raw, baseline_config)
        experimental_smooth = _smooth(raw, experimental_config)
        reference_rts = [value for name in ("Peak1", "Peak2") if (value := number(row.get(name))) is not None]
        reference_status = reference_chromatographic_status(row)
        source_machine_id = str(source_machine_label_column or "").strip().upper()
        source_machine_label = str(row.get(source_machine_label_column, "") or "").strip() if source_machine_label_column else ""
        source_machine_interpretation = source_machine_label_interpretation(source_machine_label)
        baseline_rts = [peak.rt_sec / 60 for peak in baseline_result.peaks]
        experimental_rts = [peak.rt_sec / 60 for peak in experimental_result.peaks]
        baseline_diagnostic = prediction_diagnostic(
            row, baseline_result.chromatographic_status, baseline_rts, rt_tolerance_min=rt_tolerance_min
        )
        experimental_diagnostic = prediction_diagnostic(
            row, experimental_result.chromatographic_status, experimental_rts, rt_tolerance_min=rt_tolerance_min
        )
        stability = _parameter_stability(
            baseline_result.chromatographic_status,
            baseline_rts,
            experimental_result.chromatographic_status,
            experimental_rts,
        )
        protected = reference_rts or baseline_rts or experimental_rts
        background = background_diagnostics(
            rts,
            raw,
            protected,
            absolute_height_floor=baseline_config.min_height or 0.0,
        )
        compound_id = str(row.get("Compound_ID") or row.get("SGC ID for Component") or f"row-{row_index}")
        ms2_status = ms2_by_compound.get(compound_id, "")
        target_uid = _target_uid(dataset_id, run_id, compound_id, mz)
        filename_base = (
            f"{_safe_name(compound_id)}__{run_id}__mz{mz:.5f}"
            f"__t-{target_uid}__cfg-{config_hash}"
        )
        metadata = {
            "reference_status": reference_status,
            "source_machine_context": (
                f"{source_machine_id} label: {source_machine_label or '—'}"
                f" ({source_machine_interpretation})"
                if source_machine_id else "no source-machine label"
            ),
            "baseline_status": baseline_result.chromatographic_status,
            "experimental_status": experimental_result.chromatographic_status,
            "baseline_diagnostic": baseline_diagnostic,
            "experimental_diagnostic": experimental_diagnostic,
            "parameter_stability": stability,
            "ms2_status": ms2_status,
            "min_height": baseline_config.min_height,
            "max_intensity": max(raw, default=0.0),
            "baseline_label": _config_label("baseline", baseline_config, median_scan_interval),
            "experimental_label": _config_label("experimental", experimental_config, median_scan_interval),
            "baseline_metrics": _result_metrics(baseline_result),
            "experimental_metrics": _result_metrics(experimental_result),
            "candidate_label": _config_label("high-recall candidates", candidate_config, median_scan_interval),
            "candidate_count": len(review_candidates),
            **background,
        }
        svg = render_eic_svg(
                compound_id=compound_id,
                target_mz=mz,
                ppm=baseline_config.eic_ppm,
                rts_sec=rts,
                raw=raw,
                baseline_smooth=baseline_smooth,
                experimental_smooth=experimental_smooth,
                manual_rts_min=reference_rts,
                baseline_peaks=baseline_result.peaks,
                experimental_peaks=experimental_result.peaks,
                review_candidates=review_candidates,
                metadata=metadata,
        )
        canonical_paths: dict[str, Path] = {}
        if "svg" in image_formats:
            canonical_paths["svg"] = asset_dirs["svg"] / f"{filename_base}.svg"
            canonical_paths["svg"].write_text(svg, encoding="utf-8")
        if "png" in image_formats:
            canonical_paths["png"] = asset_dirs["png"] / f"{filename_base}.png"
            png = render_eic_png(
                compound_id=compound_id,
                target_mz=mz,
                ppm=baseline_config.eic_ppm,
                rts_sec=rts,
                raw=raw,
                baseline_smooth=baseline_smooth,
                experimental_smooth=experimental_smooth,
                manual_rts_min=reference_rts,
                baseline_peaks=baseline_result.peaks,
                experimental_peaks=experimental_result.peaks,
                review_candidates=review_candidates,
                metadata=metadata,
                scale=png_scale,
            )
            _write_png(canonical_paths["png"], png)
        artifact_sha256 = hashlib.sha256(svg.encode("utf-8")).hexdigest()
        png_sha256 = (
            hashlib.sha256(canonical_paths["png"].read_bytes()).hexdigest()
            if "png" in canonical_paths
            else ""
        )
        views = review_views(
            row,
            reference_status=reference_status,
            baseline_status=baseline_result.chromatographic_status,
            baseline_diagnostic=baseline_diagnostic,
            experimental_diagnostic=experimental_diagnostic,
            low_clean_candidate=bool(background["low_clean_candidate"]),
            parameter_stability=stability,
            source_machine_id=source_machine_id,
            source_machine_label=source_machine_label,
            source_machine_interpretation=source_machine_interpretation,
        )
        if baseline_result.resolution is not None and baseline_result.resolution < 1:
            views.append("review_queues/low_resolution_fwhm_lt1")
        if len(baseline_result.peaks) >= 2:
            height_ratio = min(peak.intensity for peak in baseline_result.peaks[:2]) / max(
                peak.intensity for peak in baseline_result.peaks[:2]
            )
            if height_ratio < 0.1:
                views.append("review_queues/weak_second_peak_ratio_lt0.1")
        if "shallow_valley" in baseline_result.review_reasons:
            views.append("review_queues/shallow_valley")
        max_intensity = max(raw, default=0.0)
        if 200_000 <= max_intensity < 500_000:
            views.append("audit/intensity_0.2_to_0.5M")
        elif 500_000 <= max_intensity < 1_000_000:
            views.append("audit/intensity_0.5_to_1M")
        if len(review_candidates) >= 12:
            views.append("review_queues/high_recall_candidate_explosion")
        views = list(dict.fromkeys(views))
        for view in views:
            destination_dir = views_dir / view
            destination_dir.mkdir(parents=True, exist_ok=True)
            view_artifacts = {}
            for image_format, canonical_path in canonical_paths.items():
                filename = f"{filename_base}.{image_format}"
                destination = destination_dir / filename
                _hardlink(canonical_path, destination)
                view_artifacts[image_format] = filename
                view_index.append(
                    {
                        "target_uid": target_uid,
                        "view": view,
                        "format": image_format,
                        "view_path": destination.relative_to(output_dir).as_posix(),
                        "asset_path": canonical_path.relative_to(output_dir).as_posix(),
                    }
                )
            view_files.setdefault(view, []).append(view_artifacts)
        result_fields = {}
        result_fields.update(_result_fields("baseline", baseline_result))
        result_fields.update(_result_fields("experimental", experimental_result))
        manifest.append(
            {
                "target_uid": target_uid,
                "artifact_sha256": artifact_sha256,
                "png_sha256": png_sha256,
                "dataset_id": dataset_id,
                "run_id": run_id,
                "config_hash": config_hash,
                "row_index": row_index,
                "compound_id": compound_id,
                "target_mz": f"{mz:.8g}",
                "eic_ppm": baseline_config.eic_ppm,
                "svg_asset_path": (
                    canonical_paths["svg"].relative_to(output_dir).as_posix()
                    if "svg" in canonical_paths else ""
                ),
                "png_asset_path": (
                    canonical_paths["png"].relative_to(output_dir).as_posix()
                    if "png" in canonical_paths else ""
                ),
                "reference_status": reference_status,
                "supplied_peak1_rt_min": row.get("Peak1", ""),
                "supplied_peak2_rt_min": row.get("Peak2", ""),
                "source_machine_id": source_machine_id,
                "source_machine_label_column": source_machine_label_column or "",
                "source_machine_label": source_machine_label,
                "source_machine_interpretation": source_machine_interpretation,
                "source_pool_id": row.get("SGC ID for Pool", ""),
                "source_pooled_well": row.get("Pooled Well", ""),
                "baseline_status": baseline_result.chromatographic_status,
                "baseline_vs_reference": baseline_diagnostic,
                "baseline_review_reasons": ";".join(baseline_result.review_reasons),
                "experimental_status": experimental_result.chromatographic_status,
                "experimental_vs_reference": experimental_diagnostic,
                "experimental_review_reasons": ";".join(experimental_result.review_reasons),
                "parameter_stability": stability,
                "consensus_auto_status": (
                    baseline_result.chromatographic_status
                    if stability == "stable_class_and_rt" else "manual_review_required"
                ),
                "supplied_close_double": str("review_queues/close_or_shoulder_double" in views).lower(),
                "low_clean_candidate": str(background["low_clean_candidate"]).lower(),
                "background_p99": background["background_p99"],
                "background_nonzero_fraction": background["background_nonzero_fraction"],
                "peak_to_background_p99": background["peak_to_background_p99"],
                "max_half_height_support_scans": background["max_half_height_support_scans"],
                "cross_stage_ms2_diagnostic_status": ms2_status,
                "high_recall_candidate_count": len(review_candidates),
                "high_recall_candidate_rts_min": ";".join(
                    f"{candidate.rt_sec/60:.6g}" for candidate in review_candidates
                ),
                "high_recall_candidate_details": ";".join(
                    f"rt={candidate.rt_sec/60:.6g}|prom={_value(candidate.prominence)}|ratio={_value(candidate.prominence_ratio)}|support={candidate.support_scans or ''}"
                    for candidate in review_candidates
                ),
                "high_recall_pair_covers_supplied_double": str(
                    _candidate_pair_covers_reference(row, review_candidates, rt_tolerance_min)
                ).lower(),
                "views": ";".join(views),
                **result_fields,
            }
        )

    _write_manifest(metadata_dir / "target_manifest.csv", manifest)
    _write_manifest(metadata_dir / "view_index.csv", view_index)
    _write_annotation_template(output_dir / "annotations" / "_template" / "review_labels.csv", manifest)
    _write_galleries(output_dir, view_files)
    summary = _summary(
        manifest,
        view_files,
        peaklist_path,
        raw_path,
        baseline_config,
        experimental_config,
        candidate_config,
    )
    summary["config_hash"] = config_hash
    summary["dataset_id"] = dataset_id
    summary["run_id"] = run_id
    (metadata_dir / "generation_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (metadata_dir / "run_manifest.json").write_text(
        json.dumps(
            {
                "dataset_id": dataset_id,
                "run_id": run_id,
                "config_hash": config_hash,
                "source_peaklist": str(peaklist_path.resolve()),
                "source_raw": str(raw_path.resolve()),
                "baseline_config": baseline_config.__dict__,
                "experimental_config": experimental_config.__dict__,
                "candidate_config": candidate_config.__dict__,
                "image_formats": image_formats,
                "png_scale": png_scale if "png" in image_formats else None,
                "source_machine_label_column": source_machine_label_column,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    (output_dir / "README.md").write_text(_readme(summary), encoding="utf-8")
    return summary


def _smooth(raw: list[float], config: PeakPickingConfig) -> list[float]:
    window = min(config.sg_window, len(raw) if len(raw) % 2 else len(raw) - 1)
    return savgol_smooth(raw, window, min(config.sg_polyorder, window - 1)) if window >= 3 else raw[:]


def _prepare_output(output_dir: Path) -> None:
    resolved = output_dir.resolve()
    if len(resolved.parts) < 3 or resolved.name in {"", ".", ".."}:
        raise ValueError(f"Unsafe output directory: {resolved}")
    output_dir.mkdir(parents=True, exist_ok=True)
    # Generated assets are reproducible.  Human annotations are deliberately
    # outside this cleanup list and survive every rerun.
    for name in ("assets", "images", "views", "galleries", "metadata"):
        target = output_dir / name
        if target.exists():
            shutil.rmtree(target)


def _hardlink(source: Path, destination: Path) -> None:
    if destination.exists():
        destination.unlink()
    try:
        os.link(source, destination)
    except OSError as error:
        raise OSError(
            f"Could not create NTFS hardlink {destination}. "
            "No silent copy fallback is used; open galleries/index.html instead."
        ) from error


def _ms2_lookup(path) -> dict[str, str]:
    if not path or not Path(path).exists():
        return {}
    rows = read_table(Path(path))
    return {
        str(row.get("Compound_ID", "")): str(
            row.get("ms2_diagnostic_status") or row.get("compound_identity_status") or row.get("enantiomer_pair_status") or ""
        )
        for row in rows
        if row.get("Compound_ID")
    }


def _write_manifest(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_galleries(output_dir: Path, views: dict[str, list[dict[str, str]]]) -> None:
    galleries = output_dir / "galleries"
    galleries.mkdir(parents=True, exist_ok=True)
    links = []
    for view, artifacts in sorted(views.items()):
        name = view.replace("/", "__") + ".html"
        view_path = "../views/" + view + "/"
        cards = "".join(_gallery_card(view_path, artifact) for artifact in artifacts)
        document = _gallery_document(view, len(artifacts), cards)
        (galleries / name).write_text(document, encoding="utf-8")
        links.append(f'<li><a href="{html.escape(name)}">{html.escape(view)}</a> ({len(artifacts)})</li>')
    (galleries / "index.html").write_text(
        '<!doctype html><meta charset="utf-8"><title>MS1 EIC review galleries</title><style>body{font:16px Segoe UI,Arial;max-width:900px;margin:40px auto;color:#172033}li{margin:8px}</style><h1>MS1 EIC review galleries</h1><ul>'
        + "".join(links) + "</ul>", encoding="utf-8"
    )


def _gallery_document(title: str, count: int, cards: str) -> str:
    return f'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>{html.escape(title)}</title><style>body{{font:14px Segoe UI,Arial;background:#eef1f5;color:#172033;margin:24px}}h1{{margin-left:1%}}main{{display:grid;grid-template-columns:repeat(auto-fit,minmax(560px,1fr));gap:16px}}figure{{margin:0;background:#fff;border:1px solid #dfe4eb;border-radius:12px;overflow:hidden}}img{{display:block;width:100%;height:auto}}figcaption{{padding:8px 12px;color:#667085}}@media(max-width:650px){{main{{grid-template-columns:1fr}}}}</style></head><body><h1>{html.escape(title)} <small>n={count}</small></h1><main>{cards}</main></body></html>'''


def _gallery_card(view_path: str, artifacts: dict[str, str]) -> str:
    preview = artifacts.get("png") or artifacts.get("svg")
    primary = artifacts.get("svg") or preview
    if not preview or not primary:
        return ""
    links = " · ".join(
        f'<a href="{html.escape(view_path + filename)}">{image_format.upper()}</a>'
        for image_format, filename in sorted(artifacts.items())
    )
    return (
        f'<figure><a href="{html.escape(view_path + primary)}"><img loading="lazy" '
        f'src="{html.escape(view_path + preview)}"></a><figcaption>{html.escape(primary)} '
        f'({links})</figcaption></figure>'
    )


def _summary(manifest, views, peaklist, raw, baseline, experimental, candidate):
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_peaklist": str(peaklist.resolve()),
        "source_raw": str(raw.resolve()),
        "images": len(manifest),
        "reference_status_counts": dict(Counter(row["reference_status"] for row in manifest)),
        "baseline_status_counts": dict(Counter(row["baseline_status"] for row in manifest)),
        "baseline_vs_reference_counts": dict(Counter(row["baseline_vs_reference"] for row in manifest)),
        "experimental_vs_reference_counts": dict(Counter(row["experimental_vs_reference"] for row in manifest)),
        "parameter_stability_counts": dict(Counter(row["parameter_stability"] for row in manifest)),
        "stable_reference_agreement": {
            "stable_rows": sum(row["parameter_stability"] == "stable_class_and_rt" for row in manifest),
            "stable_exact_reference_rows": sum(
                row["parameter_stability"] == "stable_class_and_rt"
                and row["baseline_vs_reference"] == "exact_agreement"
                for row in manifest
            ),
        },
        "low_clean_candidates": sum(row["low_clean_candidate"] == "true" for row in manifest),
        "view_counts": {name: len(files) for name, files in sorted(views.items())},
        "baseline_config": baseline.__dict__,
        "experimental_config": experimental.__dict__,
        "candidate_config": candidate.__dict__,
        "high_recall_candidate_pair_coverage": sum(
            row["high_recall_pair_covers_supplied_double"] == "true" for row in manifest
        ),
    }


def _readme(summary: dict) -> str:
    return f"""# MS1 EIC review export

Generated {summary['images']} canonical chromatogram images in SVG and PNG
formats.  `assets/` is the single canonical set.  `views/` contains NTFS
hard-linked folder views by supplied RT
reference, automatic status, parameter stability, and focused review reason.
Open `galleries/index.html` for browser-based thumbnail review.

`Peak1`/`Peak2` are called a supplied/legacy reference, not final manual truth.
Copy `annotations/_template/review_labels.csv` to
`annotations/<reviewer>/round_1/labels.csv` before entering labels.  The
exporter never deletes `annotations/`.  MS2 status is present only as a
cross-stage manifest field and never changes an MS1 folder classification.

The conservative consensus accepts a row only when baseline and experimental
models have the same class and every picked RT agrees within 0.05 min.  Other
rows are routed to `review_queues/unstable_*`; this is an abstaining review
layer, not a new truth label.
"""


def _normalise_image_formats(image_formats) -> tuple[str, ...]:
    formats = tuple(dict.fromkeys(str(image_format).lower() for image_format in image_formats))
    invalid = set(formats).difference({"svg", "png"})
    if invalid or not formats:
        raise ValueError("image_formats must contain svg and/or png")
    return formats


def _write_png(path: Path, image) -> None:
    """Write a deterministic RGB PNG from the same EIC rendering inputs."""

    # Fast, modest compression keeps a full 525-target review export practical.
    # SVG remains the archival/vector artifact; PNG is a convenience copy for
    # file-browser thumbnails, PowerPoint, and reviewers without SVG support.
    image.save(path, format="PNG", optimize=False, compress_level=3)


def _peak_rt(result, index: int) -> str:
    return f"{result.peaks[index].rt_sec / 60:.6g}" if len(result.peaks) > index else ""


def _safe_name(value: str) -> str:
    cleaned = "".join(character if character.isalnum() or character in "-_" else "_" for character in value)
    return cleaned[:80] or "unknown"


def _config_hash(
    baseline: PeakPickingConfig,
    experimental: PeakPickingConfig,
    candidate: PeakPickingConfig,
) -> str:
    payload = json.dumps(
        {
            "baseline": baseline.__dict__,
            "experimental": experimental.__dict__,
            "candidate": candidate.__dict__,
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:10]


def _target_uid(dataset_id: str, run_id: str, compound_id: str, mz: float) -> str:
    payload = f"{dataset_id}|{run_id}|{compound_id}|{mz:.8f}"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]


def _candidate_pair_covers_reference(row: dict, candidates, tolerance_min: float) -> bool:
    first, second = number(row.get("Peak1")), number(row.get("Peak2"))
    if first is None or second is None:
        return False
    candidate_rts = [candidate.rt_sec / 60 for candidate in candidates]
    return any(
        abs(left - first) <= tolerance_min and abs(right - second) <= tolerance_min
        for left in candidate_rts
        for right in candidate_rts
        if left != right
    )


def _parameter_stability(
    baseline_status: str,
    baseline_rts: list[float],
    experimental_status: str,
    experimental_rts: list[float],
    *,
    rt_tolerance_min: float = 0.05,
) -> str:
    if baseline_status != experimental_status:
        return "unstable_class_across_smoothing"
    if baseline_status == "ambiguous":
        return "ambiguous_in_both_models"
    if len(baseline_rts) != len(experimental_rts) or any(
        abs(left - right) > rt_tolerance_min
        for left, right in zip(baseline_rts, experimental_rts)
    ):
        return "unstable_rt_across_smoothing"
    return "stable_class_and_rt"


def _config_label(name: str, config: PeakPickingConfig, scan_interval_sec: float) -> str:
    smooth_width = config.sg_window * scan_interval_sec
    distance = (
        f"{config.min_distance_sec:g}s"
        if config.min_distance_sec is not None
        else f"{config.min_distance_scans} scans"
    )
    return (
        f"{name}: SG{config.sg_window}/p{config.sg_polyorder} (~{smooth_width:.1f}s), "
        f"minH={_compact(config.min_height)}, distance={distance}, "
        f"valley<={_compact(config.max_valley_ratio)}, second>={_compact(config.min_second_peak_ratio)}"
    )


def _result_metrics(result) -> str:
    peaks = []
    for index, peak in enumerate(result.peaks[:2], start=1):
        peaks.append(
            f"P{index} RT={peak.rt_sec/60:.3f}, I={_compact(peak.intensity)}, "
            f"FWHM={_compact(peak.width_sec)}s, area={_compact(peak.area_fwhm)}, "
            f"SNR={_compact(peak.snr)}, prom={_compact(peak.prominence)}"
        )
    pair = (
        f"sep={_compact(result.peak_separation_sec)}s, Rs(FWHM)={_compact(result.resolution)}, "
        f"valley={_compact(result.valley_ratio)}"
    )
    return ("; ".join(peaks) or "no picked peak") + "; " + pair


def _result_fields(prefix: str, result) -> dict:
    fields = {
        f"{prefix}_peak_resolution_fwhm": _value(result.resolution),
        f"{prefix}_peak_separation_sec": _value(result.peak_separation_sec),
        f"{prefix}_valley_ratio": _value(result.valley_ratio),
    }
    for index in range(2):
        peak = result.peaks[index] if len(result.peaks) > index else None
        number = index + 1
        fields.update(
            {
                f"{prefix}_peak{number}_rt_min": "" if peak is None else _value(peak.rt_sec / 60),
                f"{prefix}_peak{number}_intensity": "" if peak is None else _value(peak.intensity),
                f"{prefix}_peak{number}_fwhm_sec": "" if peak is None else _value(peak.width_sec),
                f"{prefix}_peak{number}_area_fwhm": "" if peak is None else _value(peak.area_fwhm),
                f"{prefix}_peak{number}_snr": "" if peak is None else _value(peak.snr),
                f"{prefix}_peak{number}_prominence": "" if peak is None else _value(peak.prominence),
                f"{prefix}_peak{number}_prominence_ratio": "" if peak is None else _value(peak.prominence_ratio),
            }
        )
    if len(result.peaks) >= 2 and max(result.peaks[0].intensity, result.peaks[1].intensity) > 0:
        fields[f"{prefix}_second_peak_height_ratio"] = _value(
            min(result.peaks[0].intensity, result.peaks[1].intensity)
            / max(result.peaks[0].intensity, result.peaks[1].intensity)
        )
    else:
        fields[f"{prefix}_second_peak_height_ratio"] = ""
    return fields


def _write_annotation_template(path: Path, manifest: list[dict]) -> None:
    if path.exists():
        return
    rows = [
        {
            "target_uid": row["target_uid"],
            "artifact_sha256": row["artifact_sha256"],
            "reviewer_id": "",
            "review_round": "1",
            "manual_peak_class": "",
            "confidence": "",
            "peak1_rt_min": "",
            "peak2_rt_min": "",
            "qc_flags": "",
            "notes": "",
            "reviewed_utc": "",
        }
        for row in manifest
    ]
    _write_manifest(path, rows)


def _value(value) -> str:
    return "" if value is None else f"{float(value):.8g}"


def _compact(value) -> str:
    if value is None:
        return "—"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if abs(number) >= 1_000_000:
        return f"{number/1_000_000:.2f}M"
    if abs(number) >= 1_000:
        return f"{number/1_000:.1f}k"
    return f"{number:.3g}"
