"""Batch exporter for canonical MS2 SVG/PNG evidence images and review views."""

from __future__ import annotations

import csv
import hashlib
import html
import json
import os
import shutil
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

try:
    from ..ms2.similarity import ChiralPairThresholds
    from ..source_context import folder_component, source_machine_label_interpretation
    from .report import annotate_rows, read_result_rows
except ImportError:
    from ms2.similarity import ChiralPairThresholds  # type: ignore[import-not-found]
    from ms2_review.report import annotate_rows, read_result_rows  # type: ignore[import-not-found]
    from source_context import (  # type: ignore[import-not-found]
        folder_component,
        source_machine_label_interpretation,
    )

from .acquisition import (
    load_method_profile,
    raw_acquisition_context,
    reconcile_acquisition,
)
from .classification import review_views, spectrum_availability
from .model import prepare_mirror_spectrum
from .png import render_ms2_review_png
from .svg import render_ms2_review_svg


def export_ms2_review(
    input_path,
    output_dir,
    *,
    standard_path=None,
    sidecar_path=None,
    method_profile_path=None,
    source_machine_label_column: str | None = "IG",
    image_formats: tuple[str, ...] = ("svg", "png"),
    png_scale: float = 2.0,
) -> dict:
    """批量生成可追溯的逐目标MS2人工复核资料。

    函数读取分析CSV及metadata sidecar，沿用当次运行的fragment匹配容差和
    判定阈值，为每行生成唯一SVG/PNG资产；随后通过硬链接将同一资产放入
    状态、问题、机器、pool和审核优先级视图。最终写出target manifest、
    view index、人工标注模板、HTML gallery和运行来源信息。

    ``standard_path``描述诊断规则版本，``method_profile_path``仅提供独立方法
    上下文。二者都不会覆盖raw中观察到的DIA边界。
    """

    input_path, output_dir = Path(input_path), Path(output_dir)
    image_formats = _normalise_formats(image_formats)
    if png_scale <= 0:
        raise ValueError("png_scale must be positive")
    sidecar_path = Path(sidecar_path) if sidecar_path else Path(str(input_path) + ".metadata.json")
    sidecar = _read_json(sidecar_path) if sidecar_path.exists() else {}
    standard = _read_json(Path(standard_path)) if standard_path else {}
    method_profile = load_method_profile(method_profile_path)
    parameters = dict(sidecar.get("parameters") or {})
    thresholds = ChiralPairThresholds(
        min_cosine=float(parameters.get("min_cosine", 0.7)),
        min_matched_peaks=int(parameters.get("min_matched_peaks", 6)),
        min_explained_intensity=float(parameters.get("min_explained_intensity", 0.5)),
        min_entropy_similarity=parameters.get("min_entropy_similarity"),
    )
    fragment_mz_tol = float(parameters.get("fragment_mz_tol", 0.01))
    fragment_mz_tol_unit = str(parameters.get("fragment_mz_tol_unit", "Da"))
    min_relative_intensity = float(parameters.get("min_relative_intensity", 0.01))
    rows = annotate_rows(read_result_rows(input_path), thresholds)
    standard_id = str(standard.get("standard_id") or "MSAI-MS2-DIAGNOSTIC-v1")
    dataset_id = folder_component(input_path.stem)
    raw_path = str(((sidecar.get("inputs") or {}).get("raw") or {}).get("path") or "")
    acquisition_context = raw_acquisition_context(raw_path)
    run_id = folder_component(Path(raw_path).stem if raw_path else dataset_id)
    dia_windows = list(sidecar.get("dia_windows") or [])
    acquisition_reconciliation = reconcile_acquisition(
        acquisition_context,
        method_profile,
        analysis_parameters=parameters,
        dia_windows=dia_windows,
        raw_input_count=1 if raw_path else 0,
    )
    config_hash = _config_hash(standard_id, parameters, thresholds, method_profile)

    _prepare_output(output_dir)
    asset_root = (
        output_dir
        / "assets"
        / "ms2"
        / folder_component(standard_id)
        / f"cfg-{config_hash}"
        / run_id
    )
    asset_dirs = {image_format: asset_root / image_format for image_format in image_formats}
    for directory in asset_dirs.values():
        directory.mkdir(parents=True, exist_ok=True)
    views_root = output_dir / "views"
    metadata_root = output_dir / "metadata"
    views_root.mkdir(parents=True, exist_ok=True)
    metadata_root.mkdir(parents=True, exist_ok=True)

    manifest: list[dict] = []
    view_index: list[dict] = []
    view_files: dict[str, list[dict[str, str]]] = {}
    for row_index, row in enumerate(rows, start=1):
        compound = str(
            row.get("Compound_ID") or row.get("SGC ID for Component") or f"row-{row_index}"
        )
        target_mz = _number(row.get("MZ"))
        source_machine_id = str(source_machine_label_column or "").strip().upper()
        source_machine_label = (
            str(row.get(source_machine_label_column, "") or "").strip()
            if source_machine_label_column
            else ""
        )
        interpretation = (
            source_machine_label_interpretation(source_machine_label)
            if source_machine_id == "IG"
            else "unconfirmed_for_source_machine"
        )
        source_pool_id = str(row.get("SGC ID for Pool") or "")
        source_pooled_well = str(row.get("Pooled Well") or "")
        target_uid = _target_uid(
            dataset_id, compound, source_pool_id, target_mz, row.get("Peak1"), row.get("Peak2")
        )
        filename_base = (
            f"{folder_component(compound)}__{run_id}__mz{_filename_mz(target_mz)}"
            f"__t-{target_uid}__cfg-{config_hash}"
        )
        mirror = prepare_mirror_spectrum(
            row,
            fragment_mz_tol=fragment_mz_tol,
            fragment_mz_tol_unit=fragment_mz_tol_unit,
            min_relative_intensity=min_relative_intensity,
        )
        source_context = (
            f"{source_machine_id} label: {source_machine_label or '—'} ({interpretation})"
            if source_machine_id
            else "no source-machine label"
        )
        render_metadata = {
            "target_uid": target_uid,
            "config_hash": config_hash,
            "standard_id": standard_id,
            "source_context": source_context,
            "source_pool_id": source_pool_id,
            "source_pooled_well": source_pooled_well,
            "fragment_mz_tol": fragment_mz_tol,
            "fragment_mz_tol_unit": fragment_mz_tol_unit,
            "min_relative_intensity": min_relative_intensity,
            "min_cosine": thresholds.min_cosine,
            "min_matched_peaks": thresholds.min_matched_peaks,
            "min_explained_intensity": thresholds.min_explained_intensity,
            "rt_half_window_sec": parameters.get("rt_half_window_sec", 10.0),
            "min_fragment_correlation": parameters.get("min_fragment_correlation", 0.9),
            "fragment_correlation_mode": parameters.get("fragment_correlation_mode", "full_window"),
            "correlation_min_relative_intensity": parameters.get(
                "correlation_min_relative_intensity", 0.05
            ),
            "min_correlation_scans": parameters.get("min_correlation_scans", 5),
            "max_fragment_apex_offset_scans": parameters.get("max_fragment_apex_offset_scans", 1),
            "min_consecutive_fragment_scans": parameters.get("min_consecutive_fragment_scans", 3),
            "dia_windows": dia_windows,
            "acquisition_context": acquisition_context,
            "method_profile": method_profile,
            "acquisition_reconciliation": acquisition_reconciliation,
        }
        svg = render_ms2_review_svg(row=row, mirror=mirror, metadata=render_metadata)
        canonical: dict[str, Path] = {}
        if "svg" in image_formats:
            canonical["svg"] = asset_dirs["svg"] / f"{filename_base}.svg"
            canonical["svg"].write_text(svg, encoding="utf-8")
        if "png" in image_formats:
            canonical["png"] = asset_dirs["png"] / f"{filename_base}.png"
            image = render_ms2_review_png(
                row=row, mirror=mirror, metadata=render_metadata, scale=png_scale
            )
            image.save(canonical["png"], format="PNG", optimize=False, compress_level=3)
        svg_sha256 = hashlib.sha256(svg.encode("utf-8")).hexdigest()
        png_sha256 = _sha256_file(canonical["png"]) if "png" in canonical else ""

        views = review_views(
            row,
            source_machine_id=source_machine_id,
            source_machine_label=source_machine_label,
            source_pool_id=source_pool_id,
        )
        coverage_view = _coverage_view(row, dia_windows)
        if coverage_view:
            views.append(coverage_view)
        views = list(dict.fromkeys(views))
        for view in views:
            destination_dir = views_root / view
            destination_dir.mkdir(parents=True, exist_ok=True)
            artifacts = {}
            for image_format, source in canonical.items():
                destination = destination_dir / f"{filename_base}.{image_format}"
                _hardlink(source, destination)
                artifacts[image_format] = destination.name
                view_index.append(
                    {
                        "target_uid": target_uid,
                        "view": view,
                        "format": image_format,
                        "view_path": destination.relative_to(output_dir).as_posix(),
                        "asset_path": source.relative_to(output_dir).as_posix(),
                    }
                )
            view_files.setdefault(view, []).append(artifacts)
        manifest.append(
            {
                "target_uid": target_uid,
                "artifact_sha256": svg_sha256,
                "png_sha256": png_sha256,
                "dataset_id": dataset_id,
                "run_id": run_id,
                "standard_id": standard_id,
                "config_hash": config_hash,
                "source_row": row_index,
                "compound_id": compound,
                "target_mz": "" if target_mz is None else f"{target_mz:.8g}",
                "supplied_peak1_rt_min": row.get("Peak1", ""),
                "supplied_peak2_rt_min": row.get("Peak2", ""),
                "source_machine_id": source_machine_id,
                "source_machine_label": source_machine_label,
                "source_machine_interpretation": interpretation,
                "source_pool_id": source_pool_id,
                "source_pooled_well": source_pooled_well,
                "ms1_reference_status": row.get("ms1_reference_status", ""),
                "ms2_diagnostic_status": row.get("ms2_diagnostic_status", ""),
                "ms2_review_priority": row.get("ms2_review_priority", ""),
                "ms2_issue_codes": row.get("ms2_issue_codes", ""),
                "spectrum_availability": spectrum_availability(row),
                "ms2_cosine": row.get("ms2_cosine", ""),
                "ms2_entropy_similarity": row.get("ms2_entropy_similarity", ""),
                "ms2_matched_peaks": row.get("ms2_matched_peaks", ""),
                "peak_a_explained_intensity": row.get("peak_a_explained_intensity", ""),
                "peak_b_explained_intensity": row.get("peak_b_explained_intensity", ""),
                "peak_a_fragment_count": row.get("peak_a_fragment_count", ""),
                "peak_b_fragment_count": row.get("peak_b_fragment_count", ""),
                "peak_a_candidate_fragment_count": row.get("peak_a_candidate_fragment_count", ""),
                "peak_b_candidate_fragment_count": row.get("peak_b_candidate_fragment_count", ""),
                "peak_a_quality_flags": row.get("peak_a_quality_flags", ""),
                "peak_b_quality_flags": row.get("peak_b_quality_flags", ""),
                "dia_window_lower": row.get("dia_window_lower", ""),
                "dia_window_upper": row.get("dia_window_upper", ""),
                "dia_window_match_count": row.get("dia_window_match_count", ""),
                "method_profile_id": method_profile.get("profile_id", ""),
                "acquisition_reconciliation_status": acquisition_reconciliation.get("status", ""),
                "acquisition_issue_codes": ";".join(
                    acquisition_reconciliation.get("issue_codes", [])
                ),
                "rt_half_window_sec": parameters.get("rt_half_window_sec", ""),
                "fragment_correlation_mode": parameters.get(
                    "fragment_correlation_mode", "full_window"
                ),
                "svg_asset_path": canonical["svg"].relative_to(output_dir).as_posix()
                if "svg" in canonical
                else "",
                "png_asset_path": canonical["png"].relative_to(output_dir).as_posix()
                if "png" in canonical
                else "",
                "views": ";".join(views),
            }
        )

    _write_csv(metadata_root / "target_manifest.csv", manifest)
    _write_csv(metadata_root / "view_index.csv", view_index)
    _write_annotation_template(
        output_dir / "annotations" / "_template" / "review_labels.csv", manifest
    )
    _write_galleries(output_dir, view_files)
    summary = _summary(manifest, view_files, input_path, sidecar_path if sidecar else None)
    summary.update(
        {
            "standard_id": standard_id,
            "config_hash": config_hash,
            "dataset_id": dataset_id,
            "run_id": run_id,
            "method_profile_id": method_profile.get("profile_id"),
            "acquisition_reconciliation_status": acquisition_reconciliation.get("status"),
            "acquisition_issue_codes": acquisition_reconciliation.get("issue_codes", []),
        }
    )
    (metadata_root / "generation_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    run_manifest = {
        "created_utc": datetime.now(UTC).isoformat(),
        "dataset_id": dataset_id,
        "run_id": run_id,
        "standard_id": standard_id,
        "config_hash": config_hash,
        "source_result": str(input_path.resolve()),
        "source_result_sha256": _sha256_file(input_path),
        "source_sidecar": str(sidecar_path.resolve()) if sidecar else None,
        "source_sidecar_sha256": _sha256_file(sidecar_path) if sidecar else None,
        "sidecar_provenance_warning": None
        if sidecar
        else "sidecar missing; explicit defaults were used",
        "source_inputs": sidecar.get("inputs", {}),
        "parameters": parameters,
        "thresholds": thresholds.__dict__,
        "dia_windows": dia_windows,
        "acquisition_context": acquisition_context,
        "method_profile": method_profile,
        "method_profile_path": str(Path(method_profile_path).resolve())
        if method_profile_path
        else None,
        "method_profile_sha256": _sha256_file(Path(method_profile_path))
        if method_profile_path
        else None,
        "acquisition_reconciliation": acquisition_reconciliation,
        "image_formats": image_formats,
        "png_scale": png_scale if "png" in image_formats else None,
        "source_machine_label_column": source_machine_label_column,
    }
    (metadata_root / "run_manifest.json").write_text(
        json.dumps(run_manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (output_dir / "README.md").write_text(
        _readme(summary, bool(sidecar), bool(method_profile)), encoding="utf-8"
    )
    return summary


def _prepare_output(output_dir: Path) -> None:
    """清理本导出器拥有的旧子目录并重新建立空输出根。"""

    resolved = output_dir.resolve()
    if len(resolved.parts) < 3 or resolved.name in {"", ".", ".."}:
        raise ValueError(f"Unsafe output directory: {resolved}")
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in ("assets", "views", "galleries", "metadata"):
        target = output_dir / name
        if target.exists():
            shutil.rmtree(target)


def _hardlink(source: Path, destination: Path) -> None:
    """优先创建硬链接复用规范资产；文件系统不支持时回退到复制。"""

    if destination.exists():
        destination.unlink()
    try:
        os.link(source, destination)
    except OSError as error:
        raise OSError(
            f"Could not create hardlink {destination}; no silent copy fallback is used."
        ) from error


def _write_csv(path: Path, rows: list[dict]) -> None:
    """以所有记录键的稳定并集写出UTF-8 CSV。"""

    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_annotation_template(path: Path, manifest: list[dict]) -> None:
    """创建仅含目标标识和空人工标签字段的审核模板。"""

    if path.exists():
        return
    rows = [
        {
            "target_uid": row["target_uid"],
            "artifact_sha256": row["artifact_sha256"],
            "reviewer_id": "",
            "review_round": "1",
            "manual_ms2_identity_label": "",
            "confidence": "",
            "peak1_spectrum_quality": "",
            "peak2_spectrum_quality": "",
            "manual_issue_codes": "",
            "diagnostic_fragment_notes": "",
            "notes": "",
            "reviewed_utc": "",
        }
        for row in manifest
    ]
    _write_csv(path, rows)


def _write_galleries(output_dir: Path, views: dict[str, list[dict[str, str]]]) -> None:
    """为每个分类视图生成可浏览SVG/PNG缩略图的HTML页面。"""

    galleries = output_dir / "galleries"
    galleries.mkdir(parents=True, exist_ok=True)
    links = []
    for view, artifacts in sorted(views.items()):
        name = view.replace("/", "__") + ".html"
        prefix = "../views/" + view + "/"
        cards = "".join(_gallery_card(prefix, artifact) for artifact in artifacts)
        (galleries / name).write_text(
            _gallery_document(view, len(artifacts), cards), encoding="utf-8"
        )
        links.append(
            f'<li><a href="{html.escape(name)}">{html.escape(view)}</a> ({len(artifacts)})</li>'
        )
    (galleries / "index.html").write_text(
        '<!doctype html><meta charset="utf-8"><title>MS2 evidence review galleries</title>'
        "<style>body{font:16px Segoe UI,Arial;max-width:960px;margin:40px auto;color:#172033}li{margin:8px}</style>"
        "<h1>MS2 evidence review galleries</h1><p>Folders are non-exclusive evidence views; MS1 peak shape never determines MS2 status.</p><ul>"
        + "".join(links)
        + "</ul>",
        encoding="utf-8",
    )


def _gallery_card(prefix: str, artifacts: dict[str, str]) -> str:
    """生成一个优先显示PNG、并链接其他格式的gallery卡片。"""

    preview, primary = (
        artifacts.get("png") or artifacts.get("svg"),
        artifacts.get("svg") or artifacts.get("png"),
    )
    if not preview or not primary:
        return ""
    links = " · ".join(
        f'<a href="{html.escape(prefix + filename)}">{image_format.upper()}</a>'
        for image_format, filename in sorted(artifacts.items())
    )
    return f'<figure><a href="{html.escape(prefix + primary)}"><img loading="lazy" src="{html.escape(prefix + preview)}"></a><figcaption>{html.escape(primary)} ({links})</figcaption></figure>'


def _gallery_document(title: str, count: int, cards: str) -> str:
    """把卡片包装成独立、无需服务器的HTML文档。"""

    return f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>{html.escape(title)}</title><style>body{{font:14px Segoe UI,Arial;background:#eef1f5;color:#172033;margin:24px}}main{{display:grid;grid-template-columns:repeat(auto-fit,minmax(560px,1fr));gap:16px}}figure{{margin:0;background:#fff;border:1px solid #dfe4eb;border-radius:12px;overflow:hidden}}img{{display:block;width:100%;height:auto}}figcaption{{padding:8px 12px;color:#667085;overflow-wrap:anywhere}}@media(max-width:650px){{main{{grid-template-columns:1fr}}}}</style></head><body><h1>{html.escape(title)} <small>n={count}</small></h1><main>{cards}</main></body></html>"""


def _coverage_view(row: dict, dia_windows: list[dict]) -> str:
    """根据目标m/z是否位于sidecar真实DIA边界内返回coverage视图。"""

    issues = set(str(row.get("ms2_issue_codes", "") or "").split(";"))
    if "NO_ACQUIRED_DIA_WINDOW" not in issues or not dia_windows:
        return (
            "acquisition_coverage/within_acquired_window"
            if _number(row.get("dia_window_match_count"))
            else ""
        )
    target = _number(row.get("MZ"))
    lower = min(float(window["lower"]) for window in dia_windows)
    upper = max(float(window["upper"]) for window in dia_windows)
    if target is not None and target < lower:
        return "acquisition_coverage/no_dia_window/below_acquired_range"
    if target is not None and target > upper:
        return "acquisition_coverage/no_dia_window/above_acquired_range"
    return "acquisition_coverage/no_dia_window/internal_gap"


def _summary(manifest, views, input_path, sidecar_path):
    """汇总资产数量、状态分布、谱可用性和来源路径。"""

    return {
        "created_utc": datetime.now(UTC).isoformat(),
        "source": str(input_path.resolve()),
        "sidecar": str(sidecar_path.resolve()) if sidecar_path else None,
        "images": len(manifest),
        "status_counts": dict(Counter(row["ms2_diagnostic_status"] for row in manifest)),
        "priority_counts": dict(Counter(row["ms2_review_priority"] for row in manifest)),
        "spectrum_availability_counts": dict(
            Counter(row["spectrum_availability"] for row in manifest)
        ),
        "issue_counts": dict(
            Counter(code for row in manifest for code in row["ms2_issue_codes"].split(";") if code)
        ),
        "view_counts": {view: len(files) for view, files in sorted(views.items())},
    }


def _readme(summary: dict, has_sidecar: bool, has_method_profile: bool) -> str:
    """生成导出目录内面向人工审核者的说明文件。"""

    provenance = (
        "The matching metadata sidecar supplied extraction parameters and acquired DIA windows."
        if has_sidecar
        else "No metadata sidecar was found; inspect the provenance warning before scientific use."
    )
    method = (
        "A reference-method profile is shown separately from raw-observed acquisition fields; discrepancies are intentional review warnings."
        if has_method_profile
        else "No external reference-method profile was attached to this export."
    )
    return f"""# MS2 static evidence review export

Generated {summary["images"]} canonical SVG images and matching PNG copies,
including honest diagnostic placeholders for records without spectra. `views/`
contains non-exclusive NTFS hardlink views by MS2 diagnostic status, review
priority, spectrum availability, issue code, acquisition coverage, source
machine, and pool. Open `galleries/index.html` for thumbnail review.

{provenance}

{method}

Copy `annotations/_template/review_labels.csv` to
`annotations/<reviewer>/round_1/labels.csv`; reruns never delete `annotations/`.
Allowed identity labels should remain same-compound support/conflict,
insufficient, not evaluable, or uncertain. This workflow does not assign R/S.
"""


def _normalise_formats(image_formats) -> tuple[str, ...]:
    """去重并验证请求的图像格式，仅接受SVG和PNG。"""

    formats = tuple(dict.fromkeys(str(value).lower() for value in image_formats))
    if not formats or set(formats).difference({"svg", "png"}):
        raise ValueError("image_formats must contain svg and/or png")
    return formats


def _read_json(path: Path) -> dict:
    """读取JSON对象；顶层不是mapping时显式报错。"""

    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _config_hash(
    standard_id: str,
    parameters: dict,
    thresholds: ChiralPairThresholds,
    method_profile: dict | None = None,
) -> str:
    """由诊断标准、分析参数、阈值和方法profile生成短配置指纹。"""

    payload = json.dumps(
        {
            "standard_id": standard_id,
            "parameters": parameters,
            "thresholds": thresholds.__dict__,
            "method_profile": method_profile or {},
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode()).hexdigest()[:10]


def _target_uid(dataset, compound, pool, mz, peak1, peak2) -> str:
    """由数据集和目标关键字段生成稳定、去重的目标ID。"""

    payload = f"{dataset}|{compound}|{pool}|{mz}|{peak1}|{peak2}"
    return hashlib.sha256(payload.encode()).hexdigest()[:12]


def _filename_mz(value) -> str:
    """把目标m/z格式化为适合文件名且可读的文本。"""

    return "unknown" if value is None else f"{value:.5f}"


def _number(value):
    """把审核字段宽松转换为浮点数；失败时返回``None``。"""

    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _sha256_file(path: Path) -> str:
    """分块计算文件SHA-256，用于manifest来源校验。"""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
