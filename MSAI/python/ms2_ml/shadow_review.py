"""Build an auditable v2-versus-shadow-ML disagreement package."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from ..ms2_review.exporter import export_ms2_review
from .contracts import (
    ML_DIAGNOSTIC_STATUSES,
    ML_OUTPUT_FIELDS,
    SHADOW_COMPARISON_FIELDS,
    V2_DIAGNOSTIC_STATUSES,
)
from .review_html import render_shadow_disagreement_html

V2_STANDARD_ID = "MSAI-MS2-DIAGNOSTIC-v2"


def build_shadow_review(
    shadow_csv,
    output_dir,
    *,
    model_artifact_path=None,
    sidecar_path=None,
    standard_path=None,
    method_profile_path=None,
    image_formats: tuple[str, ...] = ("svg", "png"),
    png_scale: float = 2.0,
) -> dict[str, Any]:
    """Export every v2/ML status disagreement as CSV, HTML, PNG, and SVG.

    The source row remains a full v2 record with its five append-only ML fields.
    ``not_attempted`` is mapped to ``not_evaluable`` only for the comparison;
    its original CSV value is preserved in all evidence artifacts.
    """

    source_path = Path(shadow_csv).resolve()
    destination = Path(output_dir).resolve()
    fields, rows = _read_csv(source_path)
    _validate_shadow_fields(fields)

    resolved_standard = (
        Path(standard_path).resolve()
        if standard_path
        else Path(__file__).resolve().parents[2] / "standards" / "ms2_diagnostic_standard_v2.json"
    )
    standard = _read_json(resolved_standard)
    if standard.get("standard_id") != V2_STANDARD_ID:
        raise ValueError(f"Shadow disagreement review requires {V2_STANDARD_ID}.")

    model_path = Path(model_artifact_path).resolve() if model_artifact_path else None
    row_model_ids = sorted({str(row.get("ml_model_id") or "").strip() for row in rows})
    if "" in row_model_ids:
        raise ValueError("Shadow CSV rows require a non-empty ml_model_id.")
    if len(row_model_ids) > 1:
        raise ValueError("Shadow CSV may contain predictions from exactly one model ID.")
    model_id = "UNTRAINED" if not row_model_ids else row_model_ids[0]
    if model_path:
        from .shadow import load_model_artifact, predict_shadow

        validated_artifact = load_model_artifact(model_path)
        artifact_model_id = str(validated_artifact["trainable_decision"]["model_id"])
        if row_model_ids and row_model_ids != [artifact_model_id]:
            raise ValueError("Shadow CSV model IDs do not match the attached model artifact.")
        model_id = artifact_model_id
        _validate_row_predictions(rows, validated_artifact, predict_shadow)
    elif model_id == "UNTRAINED":
        from .shadow import predict_shadow

        _validate_row_predictions(rows, None, predict_shadow)

    comparison_rows: list[dict[str, Any]] = []
    transition_counts: Counter[str] = Counter()
    for source_row in rows:
        v2_status = str(source_row.get("ms2_diagnostic_status") or "").strip()
        ml_status = str(source_row.get("ml_diagnostic_status") or "").strip()
        if v2_status not in V2_DIAGNOSTIC_STATUSES:
            raise ValueError(f"Invalid v2 diagnostic status: {v2_status!r}.")
        if ml_status not in ML_DIAGNOSTIC_STATUSES:
            raise ValueError(f"Invalid ML diagnostic status: {ml_status!r}.")
        normalized = _normalise_v2_status(v2_status)
        transition_counts[f"{normalized}->{ml_status}"] += 1
        if normalized == ml_status:
            continue
        row = dict(source_row)
        row.update(
            {
                "comparison_change_level": "v2_ml_status",
                "comparison_baseline_status": v2_status,
                "comparison_v2_normalized_status": normalized,
                "comparison_candidate_status": ml_status,
                "comparison_review_tier": "v2_ml_disagreement",
                "comparison_review_flags": source_row.get("ml_review_flags", ""),
            }
        )
        comparison_rows.append(row)

    destination.mkdir(parents=True, exist_ok=True)
    disagreement_csv = destination / "v2_ml_disagreements.csv"
    comparison_fields = [*fields, *SHADOW_COMPARISON_FIELDS]
    _write_csv(disagreement_csv, comparison_fields, comparison_rows)

    resolved_sidecar = (
        Path(sidecar_path).resolve() if sidecar_path else Path(str(source_path) + ".metadata.json")
    )
    render_summary = export_ms2_review(
        disagreement_csv,
        destination,
        standard_path=resolved_standard,
        sidecar_path=resolved_sidecar if resolved_sidecar.exists() else None,
        method_profile_path=method_profile_path,
        image_formats=image_formats,
        png_scale=png_scale,
        preserve_input_diagnostics=True,
        model_artifact_path=model_artifact_path,
    )
    manifest_path = destination / "metadata" / "target_manifest.csv"
    manifest = _read_manifest(manifest_path)
    if len(manifest) != len(comparison_rows):
        raise AssertionError("Review manifest does not contain one row per disagreement.")

    disagreement_html = destination / "v2_ml_disagreements.html"
    disagreement_html.write_text(
        render_shadow_disagreement_html(comparison_rows, manifest), encoding="utf-8"
    )

    summary = {
        "schema_version": "msai-ms2-shadow-disagreement-review-v1",
        "created_utc": datetime.now(UTC).isoformat(),
        "source_shadow_csv": str(source_path),
        "source_shadow_sha256": _sha256_file(source_path),
        "standard": {
            "standard_id": V2_STANDARD_ID,
            "path": str(resolved_standard),
            "sha256": _sha256_file(resolved_standard),
        },
        "model": {
            "model_id": model_id,
            "path": str(model_path) if model_path else None,
            "sha256": _sha256_file(model_path) if model_path else None,
            "artifact_attached": model_path is not None,
            "provenance_warning": (
                None
                if model_path or model_id == "UNTRAINED"
                else "trained model artifact was not attached to this review export"
            ),
        },
        "rows": len(rows),
        "agreements": len(rows) - len(comparison_rows),
        "disagreements": len(comparison_rows),
        "transition_counts": dict(sorted(transition_counts.items())),
        "outputs": {
            "csv": str(disagreement_csv),
            "html": str(disagreement_html),
            "target_manifest": str(manifest_path) if manifest_path.exists() else None,
            "png_images": len(list((destination / "assets").rglob("*.png"))),
            "svg_images": len(list((destination / "assets").rglob("*.svg"))),
        },
        "render_summary": render_summary,
    }
    summary_path = destination / "shadow_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary


def _validate_shadow_fields(fields: list[str]) -> None:
    if len(fields) != len(set(fields)):
        raise ValueError("Shadow CSV header fields must be unique.")
    missing = [
        field for field in ("ms2_diagnostic_status", *ML_OUTPUT_FIELDS) if field not in fields
    ]
    if missing:
        raise ValueError("Shadow CSV is missing required fields: " + ", ".join(missing))
    if tuple(fields[-len(ML_OUTPUT_FIELDS) :]) != ML_OUTPUT_FIELDS:
        raise ValueError("Shadow CSV must append the exact five ML output fields at the end.")


def _validate_row_predictions(rows, artifact, predictor) -> None:
    for index, row in enumerate(rows):
        expected = predictor(row, artifact)
        observed = {field: str(row.get(field) or "") for field in ML_OUTPUT_FIELDS}
        if observed != expected:
            raise ValueError(
                f"Shadow CSV row {index} outputs do not match its declared model prediction."
            )


def _normalise_v2_status(status: str) -> str:
    return "not_evaluable" if status == "not_attempted" else status


def _read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("Shadow CSV must contain a header.")
        return list(reader.fieldnames), list(reader)


def _read_manifest(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, fields, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def _read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {path}.")
    return payload


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


__all__ = ["build_shadow_review"]
