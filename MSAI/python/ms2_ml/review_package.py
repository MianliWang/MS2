"""Build a compact, reproducible ZIP for manual MS2 shadow review."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import tempfile
import zipfile
import zlib
from html.parser import HTMLParser
from pathlib import Path, PurePosixPath
from typing import Any
from xml.etree import ElementTree

from ..ms2_review.acquisition import reconcile_acquisition
from ..ms2_review.model import prepare_mirror_spectrum
from ..ms2_review.render_context import build_review_render_metadata, review_config_hash
from ..ms2_review.svg import render_ms2_review_svg
from ..source_context import folder_component, source_machine_label_interpretation
from .contracts import (
    ML_DIAGNOSTIC_STATUSES,
    ML_OUTPUT_FIELDS,
    SHADOW_COMPARISON_FIELDS,
    V2_DIAGNOSTIC_STATUSES,
    validate_frozen_source,
)
from .provenance import canonical_json_sha256, canonical_table_sha256
from .review_html import render_shadow_disagreement_html

PACKAGE_SCHEMA_VERSION = "msai-ms2-shadow-manual-review-package-v1"
PACKAGE_ROOT = "MS2_manual_review"
REQUIRED_REVIEW_FILES = {
    "html": "v2_ml_disagreements.html",
    "queue": "v2_ml_disagreements.csv",
    "labels": "annotations/_template/review_labels.csv",
    "manifest": "metadata/target_manifest.csv",
    "run_manifest": "metadata/run_manifest.json",
    "summary": "shadow_summary.json",
}
_IMAGE_SUFFIXES = frozenset({".png", ".svg"})
_V2_STANDARD_ID = "MSAI-MS2-DIAGNOSTIC-v2"
_REVIEW_LABEL_FIELDS = (
    "target_uid",
    "artifact_sha256",
    "reviewer_id",
    "review_round",
    "manual_ms2_identity_label",
    "confidence",
    "peak1_spectrum_quality",
    "peak2_spectrum_quality",
    "manual_issue_codes",
    "diagnostic_fragment_notes",
    "notes",
    "reviewed_utc",
)


class _ImageReferenceParser(HTMLParser):
    """Collect local PNG/SVG references without interpreting arbitrary HTML."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.references: set[str] = set()

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        del tag
        for name, value in attrs:
            if name not in {"href", "src"} or not value:
                continue
            reference = value.split("#", 1)[0].split("?", 1)[0]
            if PurePosixPath(reference).suffix.lower() in _IMAGE_SUFFIXES:
                self.references.add(reference)


def package_shadow_review(
    shadow_csv: str | Path,
    review_dir: str | Path,
    output_path: str | Path | None = None,
    *,
    model_artifact_path: str | Path | None = None,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Package one final shadow CSV and its canonical manual-review images.

    The archive deliberately excludes duplicated ``views/`` hardlinks,
    galleries, and implementation artifacts.  Every image referenced by the
    review HTML is required to appear exactly once in the canonical manifest
    and is verified against its recorded SHA-256 before the ZIP is written.
    """

    result_path = Path(shadow_csv).resolve()
    review_root = Path(review_dir).resolve()
    destination = (
        Path(output_path).resolve()
        if output_path is not None
        else (Path.cwd() / "share" / f"{result_path.stem}_manual_review.zip").resolve()
    )
    if destination.suffix.lower() != ".zip":
        raise ValueError("Manual-review package output must use a .zip suffix.")
    if not result_path.is_file():
        raise FileNotFoundError(f"Final shadow CSV does not exist: {result_path}")
    if not review_root.is_dir():
        raise FileNotFoundError(f"Shadow review directory does not exist: {review_root}")
    if destination.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists; pass overwrite=True to replace it: {destination}"
        )

    required = {
        name: _required_file(review_root, relative)
        for name, relative in REQUIRED_REVIEW_FILES.items()
    }
    sidecar_path = Path(str(result_path) + ".metadata.json")
    if not sidecar_path.is_file():
        raise FileNotFoundError(f"Final result metadata sidecar does not exist: {sidecar_path}")

    result_bytes = result_path.read_bytes()
    sidecar_bytes = sidecar_path.read_bytes()
    metadata = _read_json_bytes(sidecar_bytes, sidecar_path)
    review_summary = _read_json_file(required["summary"])
    run_manifest = _read_json_file(required["run_manifest"])
    result_fields, result_rows = _read_csv_bytes(result_bytes, result_path)
    _validate_result_fields(result_fields)
    _validate_result_sidecar(metadata, result_bytes, result_rows)
    standard_sha256 = _validate_source_provenance(metadata, result_fields, result_rows)
    _validate_review_summary(
        review_summary,
        result_bytes,
        result_rows,
        standard_sha256,
    )
    render_context = _authoritative_render_context(metadata)

    manifest_fields, manifest_rows = _read_csv_bytes(
        required["manifest"].read_bytes(), required["manifest"]
    )
    queue_bytes = required["queue"].read_bytes()
    queue_fields, queue_rows = _read_csv_bytes(queue_bytes, required["queue"])
    labels_bytes = required["labels"].read_bytes()
    label_fields, label_rows = _read_csv_bytes(labels_bytes, required["labels"])
    _validate_review_tables(
        result_fields,
        result_rows,
        manifest_fields,
        manifest_rows,
        queue_fields,
        queue_rows,
        label_fields,
        label_rows,
    )
    _validate_run_manifest(
        run_manifest,
        queue_bytes,
        sidecar_bytes,
        metadata,
        review_summary,
        render_context,
    )
    asset_entries = _validated_assets(
        review_root,
        manifest_rows,
        queue_rows,
        render_context,
        png_scale=_positive_number(run_manifest.get("png_scale"), "png_scale"),
    )
    html_bytes = required["html"].read_bytes()
    expected_html = render_shadow_disagreement_html(queue_rows, manifest_rows).encode("utf-8")
    if html_bytes != expected_html:
        raise ValueError("Review HTML is not the exact deterministic rendering of the queue.")
    html_references = _html_image_references(html_bytes)
    if html_references != set(asset_entries):
        missing = sorted(set(asset_entries) - html_references)
        unexpected = sorted(html_references - set(asset_entries))
        raise ValueError(
            "Review HTML/image manifest mismatch: "
            f"missing_references={missing}, unexpected_references={unexpected}."
        )

    _validate_cross_file_counts(review_summary, result_rows, manifest_rows, asset_entries)
    model, decision = _model_summary(
        metadata,
        review_summary,
        result_rows,
        _source_decision_parameters(metadata),
        model_artifact_path,
    )
    entries: dict[str, bytes] = {
        f"final_result/{result_path.name}": result_bytes,
        f"final_result/{sidecar_path.name}": sidecar_bytes,
        "START_HERE.html": html_bytes,
        "review_queue.csv": queue_bytes,
        "review_labels.csv": _render_review_labels(manifest_rows, label_rows),
        **asset_entries,
    }
    entries["README_FIRST.txt"] = _render_readme(
        result_path.name,
        len(result_rows),
        len(manifest_rows),
        model,
        decision,
    )
    entries["package_summary.json"] = _render_summary(
        result_path,
        result_bytes,
        len(result_rows),
        len(manifest_rows),
        asset_entries,
        model,
        decision,
    )
    entries["manifest_sha256.csv"] = _render_hash_manifest(entries)
    _validate_entry_names(entries)

    destination.parent.mkdir(parents=True, exist_ok=True)
    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
    )
    os.close(file_descriptor)
    temporary = Path(temporary_name)
    try:
        _write_archive(temporary, entries)
        _validate_archive(temporary, entries, html_references)
        if destination.exists() and not overwrite:
            raise FileExistsError(
                f"Output appeared while packaging; refusing to replace it: {destination}"
            )
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)

    return {
        "schema_version": PACKAGE_SCHEMA_VERSION,
        "output": str(destination),
        "sha256": _sha256_bytes(destination.read_bytes()),
        "size_bytes": destination.stat().st_size,
        "result_rows": len(result_rows),
        "review_rows": len(manifest_rows),
        "png_images": sum(path.endswith(".png") for path in asset_entries),
        "svg_images": sum(path.endswith(".svg") for path in asset_entries),
        "model_id": model["model_id"],
        "trained": model["trained"],
    }


def _required_file(root: Path, relative: str) -> Path:
    path = root / relative
    if not path.is_file():
        raise FileNotFoundError(f"Required shadow-review file is missing: {path}")
    return path


def _read_json_file(path: Path) -> dict[str, Any]:
    return _read_json_bytes(path.read_bytes(), path)


def _read_json_bytes(data: bytes, label: Path) -> dict[str, Any]:
    value = json.loads(data.decode("utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {label}.")
    return value


def _read_csv_bytes(data: bytes, label: Path) -> tuple[list[str], list[dict[str, str]]]:
    with io.StringIO(data.decode("utf-8-sig"), newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"CSV must contain a header: {label}")
        return list(reader.fieldnames), list(reader)


def _validate_result_fields(fields: list[str]) -> None:
    if len(fields) != len(set(fields)):
        raise ValueError("Final shadow CSV contains duplicate column names.")
    required = {"ms2_diagnostic_status", *ML_OUTPUT_FIELDS}
    missing = sorted(required - set(fields))
    if missing:
        raise ValueError("Final shadow CSV is missing required fields: " + ", ".join(missing))
    if tuple(fields[-len(ML_OUTPUT_FIELDS) :]) != ML_OUTPUT_FIELDS:
        raise ValueError("Final shadow CSV must append the exact five ML fields at its tail.")


def _validate_result_sidecar(
    metadata: dict[str, Any], result_bytes: bytes, rows: list[dict[str, str]]
) -> None:
    inference = metadata.get("shadow_inference")
    if not isinstance(inference, dict):
        raise ValueError("Final result metadata is missing shadow_inference provenance.")
    expected_counts: dict[str, int] = {}
    for row in rows:
        status = str(row.get("ml_diagnostic_status") or "")
        expected_counts[status] = expected_counts.get(status, 0) + 1
    if inference.get("schema_version") != "msai-ms2-shadow-output-v1":
        raise ValueError("Unsupported final shadow-result metadata schema.")
    if inference.get("output_result_sha256") != _sha256_bytes(result_bytes):
        raise ValueError("Final result metadata does not fingerprint the supplied CSV.")
    if inference.get("rows") != len(rows):
        raise ValueError("Final result metadata row count does not match the supplied CSV.")
    if inference.get("status_counts") != expected_counts:
        raise ValueError("Final result metadata status counts do not match the supplied CSV.")
    if inference.get("append_only_fields") != list(ML_OUTPUT_FIELDS):
        raise ValueError("Final result metadata does not record the exact five appended ML fields.")


def _validate_source_provenance(
    metadata: dict[str, Any],
    result_fields: list[str],
    result_rows: list[dict[str, str]],
) -> str:
    inference = metadata["shadow_inference"]
    source_fields = result_fields[: -len(ML_OUTPUT_FIELDS)]
    source_rows = [{field: row.get(field, "") for field in source_fields} for row in result_rows]
    expected_result_content_hash = inference.get("source_result_content_sha256")
    if expected_result_content_hash is not None:
        if expected_result_content_hash != canonical_table_sha256(source_fields, source_rows):
            raise ValueError("Final result does not match the recorded source-result content hash.")
    else:
        reconstructed_result = _render_csv_bytes(source_fields, source_rows)
        if inference.get("source_result_sha256") != _sha256_bytes(reconstructed_result):
            raise ValueError(
                "Legacy final result cannot be reconstructed to the recorded source-result SHA-256."
            )

    source_sidecar = dict(metadata)
    source_sidecar.pop("shadow_inference", None)
    expected_sidecar_content_hash = inference.get("source_sidecar_content_sha256")
    if expected_sidecar_content_hash is not None:
        if expected_sidecar_content_hash != canonical_json_sha256(source_sidecar):
            raise ValueError(
                "Final metadata does not match the recorded source-sidecar content hash."
            )
    else:
        reconstructed_sidecar = json.dumps(
            source_sidecar,
            indent=2,
            ensure_ascii=False,
        ).encode("utf-8")
        if inference.get("source_sidecar_sha256") != _sha256_bytes(reconstructed_sidecar):
            raise ValueError(
                "Legacy final metadata cannot be reconstructed to the recorded "
                "source-sidecar SHA-256."
            )

    standard_value = inference.get("diagnostic_standard")
    expected_standard_hash = inference.get("diagnostic_standard_sha256")
    if not isinstance(standard_value, str) or not standard_value:
        raise ValueError("Final result metadata is missing the diagnostic-standard path.")
    if not isinstance(expected_standard_hash, str) or len(expected_standard_hash) != 64:
        raise ValueError("Final result metadata has an invalid diagnostic-standard SHA-256.")
    standard_path = Path(standard_value).resolve()
    if not standard_path.is_file():
        raise FileNotFoundError(f"Recorded diagnostic standard does not exist: {standard_path}")
    standard_bytes = standard_path.read_bytes()
    if _sha256_bytes(standard_bytes) != expected_standard_hash:
        raise ValueError("Recorded diagnostic standard does not match its SHA-256.")
    standard = _read_json_bytes(standard_bytes, standard_path)
    validate_frozen_source(source_sidecar, standard)
    return expected_standard_hash


def _validate_review_summary(
    summary: dict[str, Any],
    result_bytes: bytes,
    result_rows: list[dict[str, str]],
    standard_sha256: str,
) -> None:
    if summary.get("schema_version") != "msai-ms2-shadow-disagreement-review-v1":
        raise ValueError("Unsupported shadow-review summary schema.")
    if summary.get("source_shadow_sha256") != _sha256_bytes(result_bytes):
        raise ValueError("Shadow-review summary does not fingerprint the supplied final CSV.")
    if summary.get("rows") != len(result_rows):
        raise ValueError("Shadow-review summary row count does not match the final CSV.")
    standard = summary.get("standard")
    if (
        not isinstance(standard, dict)
        or standard.get("standard_id") != _V2_STANDARD_ID
        or standard.get("sha256") != standard_sha256
    ):
        raise ValueError("Shadow-review summary does not match the frozen diagnostic standard.")

    expected_transitions: dict[str, int] = {}
    for row in result_rows:
        v2 = str(row.get("ms2_diagnostic_status") or "").strip()
        ml = str(row.get("ml_diagnostic_status") or "").strip()
        if v2 not in V2_DIAGNOSTIC_STATUSES:
            raise ValueError(f"Invalid v2 diagnostic status in final result: {v2!r}.")
        if ml not in ML_DIAGNOSTIC_STATUSES:
            raise ValueError(f"Invalid ML diagnostic status in final result: {ml!r}.")
        transition = f"{_normalise_v2_status(v2)}->{ml}"
        expected_transitions[transition] = expected_transitions.get(transition, 0) + 1
    if summary.get("transition_counts") != dict(sorted(expected_transitions.items())):
        raise ValueError("Shadow-review transition counts do not match the final result.")


def _validate_run_manifest(
    manifest: dict[str, Any],
    queue_bytes: bytes,
    sidecar_bytes: bytes,
    metadata: dict[str, Any],
    review_summary: dict[str, Any],
    render_context: dict[str, Any],
) -> None:
    if manifest.get("source_result_sha256") != _sha256_bytes(queue_bytes):
        raise ValueError("Review run manifest does not fingerprint the disagreement queue.")
    if manifest.get("source_sidecar_sha256") != _sha256_bytes(sidecar_bytes):
        raise ValueError("Review run manifest does not fingerprint the final result sidecar.")
    for key in ("dataset_id", "run_id", "standard_id", "config_hash"):
        if manifest.get(key) != render_context[key]:
            raise ValueError(
                f"Review run manifest field {key!r} does not match sidecar-derived provenance."
            )
    if manifest.get("preserve_input_diagnostics") is not True:
        raise ValueError("Shadow review must preserve input diagnostic fields.")
    if manifest.get("source_machine_label_column") != render_context["source_machine_label_column"]:
        raise ValueError(
            "Review run manifest source-machine column does not match the review contract."
        )
    image_formats = manifest.get("image_formats")
    if not isinstance(image_formats, list) or set(image_formats) != {"svg", "png"}:
        raise ValueError("Manual-review packaging requires SVG and PNG render formats.")
    _positive_number(manifest.get("png_scale"), "png_scale")
    for key in (
        "parameters",
        "thresholds",
        "dia_windows",
        "acquisition_context",
        "method_profile",
        "acquisition_reconciliation",
    ):
        if manifest.get(key) != render_context[key]:
            raise ValueError(
                f"Review run manifest field {key!r} does not match final-sidecar provenance."
            )

    inference = metadata["shadow_inference"]
    model = inference.get("model")
    run_model = manifest.get("model_artifact")
    summary_model = review_summary.get("model")
    if not isinstance(model, dict) or not isinstance(run_model, dict):
        raise ValueError("Review run manifest is missing model provenance.")
    expected_model = {
        "model_id": model.get("model_id"),
        "sha256": model.get("sha256"),
    }
    if {key: run_model.get(key) for key in expected_model} != expected_model:
        raise ValueError("Review run manifest model provenance does not match the result.")
    if (
        not isinstance(summary_model, dict)
        or {key: summary_model.get(key) for key in expected_model} != expected_model
    ):
        raise ValueError("Shadow-review summary model provenance does not match the result.")


def _authoritative_render_context(metadata: dict[str, Any]) -> dict[str, Any]:
    """Derive every scientific render input from the final result sidecar."""

    parameters = _required_mapping(metadata, "parameters")
    thresholds = _source_decision_parameters(metadata)
    dia_windows = metadata.get("dia_windows", [])
    if not isinstance(dia_windows, list) or not all(
        isinstance(window, dict) for window in dia_windows
    ):
        raise ValueError("Final sidecar dia_windows must be a list of mappings.")

    layers_value = metadata.get("provenance_layers")
    if layers_value is None:
        layers: dict[str, Any] = {}
    elif isinstance(layers_value, dict):
        layers = layers_value
    else:
        raise ValueError("Final sidecar provenance_layers must be a mapping.")

    inputs_value = metadata.get("inputs")
    if inputs_value is None:
        inputs: dict[str, Any] = {}
    elif isinstance(inputs_value, dict):
        inputs = inputs_value
    else:
        raise ValueError("Final sidecar inputs must be a mapping.")
    raw_entry = inputs.get("raw") or {}
    if not isinstance(raw_entry, dict):
        raise ValueError("Final sidecar inputs.raw must be a mapping.")
    raw_path = str(raw_entry.get("path") or "")

    raw_observed = layers.get("raw_observed")
    if raw_observed is None:
        if raw_path:
            raise ValueError(
                "Final sidecar with a raw input must include provenance_layers.raw_observed."
            )
        acquisition_context: dict[str, Any] = {}
    elif isinstance(raw_observed, dict):
        acquisition_context = dict(raw_observed)
    else:
        raise ValueError("Final sidecar raw_observed provenance must be a mapping.")

    reference_method = layers.get("reference_method")
    if reference_method is None:
        method_profile: dict[str, Any] = {}
    elif isinstance(reference_method, dict):
        method_profile = dict(reference_method)
    else:
        raise ValueError("Final sidecar reference_method provenance must be a mapping.")

    inference = _required_mapping(metadata, "shadow_inference")
    model = _required_mapping(inference, "model")
    model_context = {
        "model_id": str(model.get("model_id") or ""),
        "sha256": model.get("sha256") or "",
    }
    dataset_id = folder_component(Path(REQUIRED_REVIEW_FILES["queue"]).stem)
    run_id = folder_component(Path(raw_path).stem if raw_path else dataset_id)
    reconciliation = reconcile_acquisition(
        acquisition_context,
        method_profile,
        analysis_parameters=parameters,
        dia_windows=dia_windows,
        raw_input_count=1 if raw_path else 0,
    )
    config_hash = review_config_hash(
        standard_id=_V2_STANDARD_ID,
        parameters=parameters,
        thresholds=thresholds,
        method_profile=method_profile,
        model_context=model_context,
    )
    return {
        "dataset_id": dataset_id,
        "run_id": run_id,
        "standard_id": _V2_STANDARD_ID,
        "config_hash": config_hash,
        "source_machine_label_column": "IG",
        "parameters": parameters,
        "thresholds": thresholds,
        "dia_windows": dia_windows,
        "acquisition_context": acquisition_context,
        "method_profile": method_profile,
        "acquisition_reconciliation": reconciliation,
        "model_context": model_context,
    }


def _render_csv_bytes(fields: list[str], rows: list[dict[str, str]], *, bom: bool = False) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, extrasaction="raise")
    writer.writeheader()
    writer.writerows(rows)
    prefix = "\ufeff" if bom else ""
    return (prefix + output.getvalue()).encode("utf-8")


def _validated_assets(
    root: Path,
    rows: list[dict[str, str]],
    queue_rows: list[dict[str, str]],
    render_context: dict[str, Any],
    *,
    png_scale: float,
) -> dict[str, bytes]:
    entries: dict[str, bytes] = {}
    for index, (row, queue) in enumerate(zip(rows, queue_rows, strict=True), start=1):
        uid = _validate_manifest_identity(row, queue, render_context, index)
        relatives: dict[str, str] = {}
        png_dimensions: tuple[int, int] | None = None
        for path_field, hash_field, suffix in (
            ("png_asset_path", "png_sha256", ".png"),
            ("svg_asset_path", "artifact_sha256", ".svg"),
        ):
            relative = _safe_relative_path(row.get(path_field, ""), suffix)
            relatives[suffix] = relative
            if relative in entries:
                raise ValueError(f"Duplicate canonical review asset in manifest: {relative}")
            path = (root / relative).resolve()
            try:
                path.relative_to(root)
            except ValueError as error:
                raise ValueError(f"Review asset escapes its root: {relative}") from error
            if not path.is_file():
                raise FileNotFoundError(f"Canonical review asset is missing: {path}")
            data = path.read_bytes()
            dimensions = _validate_image_payload(data, suffix, relative)
            if suffix == ".png":
                png_dimensions = dimensions
            expected_hash = str(row.get(hash_field) or "")
            if len(expected_hash) != 64 or expected_hash != _sha256_bytes(data):
                raise ValueError(f"Review asset SHA-256 mismatch: {relative}")
            entries[relative] = data
        png_stem = PurePosixPath(relatives[".png"]).stem
        svg_stem = PurePosixPath(relatives[".svg"]).stem
        if png_stem != svg_stem:
            raise ValueError(f"Canonical PNG/SVG stems do not match for target {uid}.")
        if f"__t-{uid}__" not in png_stem:
            raise ValueError(f"Canonical asset stem does not contain target_uid {uid}.")
        expected_svg, expected_png = _render_expected_assets(
            queue,
            uid,
            render_context,
            png_scale=png_scale,
        )
        if entries[relatives[".svg"]] != expected_svg:
            raise ValueError(f"Canonical SVG is not the deterministic rendering for target {uid}.")
        expected_dimensions = expected_png.size
        if png_dimensions != expected_dimensions:
            raise ValueError(
                f"Canonical PNG dimensions do not match render scale for target {uid}."
            )
        _validate_png_pixels(entries[relatives[".png"]], expected_png, relatives[".png"])
    return entries


def _render_expected_assets(
    queue: dict[str, str],
    target_uid: str,
    render_context: dict[str, Any],
    *,
    png_scale: float,
) -> tuple[bytes, Any]:
    machine_id = str(render_context["source_machine_label_column"] or "").strip().upper()
    machine_label = str(queue.get(machine_id) or "").strip() if machine_id else ""
    interpretation = (
        source_machine_label_interpretation(machine_label)
        if machine_id == "IG"
        else "unconfirmed_for_source_machine"
    )
    metadata = build_review_render_metadata(
        target_uid=target_uid,
        config_hash=str(render_context["config_hash"]),
        standard_id=str(render_context["standard_id"]),
        source_machine_id=machine_id,
        source_machine_label=machine_label,
        source_machine_interpretation=interpretation,
        source_pool_id=str(queue.get("SGC ID for Pool") or ""),
        source_pooled_well=str(queue.get("Pooled Well") or ""),
        parameters=render_context["parameters"],
        thresholds=render_context["thresholds"],
        dia_windows=render_context["dia_windows"],
        acquisition_context=render_context["acquisition_context"],
        method_profile=render_context["method_profile"],
        acquisition_reconciliation=render_context["acquisition_reconciliation"],
        model_context=render_context["model_context"],
    )
    mirror = prepare_mirror_spectrum(
        queue,
        fragment_mz_tol=float(metadata["fragment_mz_tol"]),
        fragment_mz_tol_unit=str(metadata["fragment_mz_tol_unit"]),
        min_relative_intensity=float(metadata["min_relative_intensity"]),
    )
    svg = render_ms2_review_svg(row=queue, mirror=mirror, metadata=metadata).encode("utf-8")
    from ..ms2_review.png import render_ms2_review_png

    png = render_ms2_review_png(
        row=queue,
        mirror=mirror,
        metadata=metadata,
        scale=png_scale,
    )
    return svg, png


def _validate_png_pixels(data: bytes, expected: Any, relative: str) -> None:
    """Compare decoded pixels with the final-sidecar/queue-derived rendering."""

    from PIL import Image, UnidentifiedImageError  # type: ignore[import-not-found]

    try:
        with Image.open(io.BytesIO(data)) as observed:
            observed.load()
            observed_rgb = observed.convert("RGB")
    except (OSError, UnidentifiedImageError) as error:
        raise ValueError(f"Canonical PNG cannot be decoded: {relative}") from error
    expected_rgb = expected.convert("RGB")
    if observed_rgb.size != expected_rgb.size or observed_rgb.tobytes() != expected_rgb.tobytes():
        raise ValueError(f"Canonical PNG is not the deterministic rendering: {relative}")


def _required_mapping(value: dict[str, Any], key: str) -> dict[str, Any]:
    candidate = value.get(key)
    if not isinstance(candidate, dict):
        raise ValueError(f"Review run manifest field {key!r} must be a mapping.")
    return candidate


def _positive_number(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ValueError(f"Review run manifest {name} must be a positive number.")
    numeric = float(value)
    if not 0 < numeric <= 8:
        raise ValueError(f"Review run manifest {name} must be in (0, 8].")
    return numeric


def _validate_manifest_identity(
    manifest: dict[str, str],
    queue: dict[str, str],
    render_context: dict[str, Any],
    index: int,
) -> str:
    for field in ("dataset_id", "run_id", "standard_id", "config_hash"):
        if str(manifest.get(field) or "") != str(render_context.get(field) or ""):
            raise ValueError(
                f"Target manifest field {field!r} does not match sidecar-derived provenance."
            )
    if manifest.get("standard_id") != _V2_STANDARD_ID:
        raise ValueError("Target manifest has the wrong diagnostic standard ID.")

    compound = str(queue.get("Compound_ID") or queue.get("SGC ID for Component") or f"row-{index}")
    pool = str(queue.get("SGC ID for Pool") or "")
    target_mz = _number(queue.get("MZ"))
    peak1 = queue.get("Peak1")
    peak2 = queue.get("Peak2")
    if str(manifest.get("compound_id") or "") != compound:
        raise ValueError("Target manifest compound_id does not match the review queue.")
    if str(manifest.get("source_pool_id") or "") != pool:
        raise ValueError("Target manifest source_pool_id does not match the review queue.")
    pooled_well = str(queue.get("Pooled Well") or "")
    if str(manifest.get("source_pooled_well") or "") != pooled_well:
        raise ValueError("Target manifest source_pooled_well does not match the review queue.")
    machine_id = str(render_context["source_machine_label_column"] or "").strip().upper()
    machine_label = str(queue.get(machine_id) or "").strip() if machine_id else ""
    interpretation = (
        source_machine_label_interpretation(machine_label)
        if machine_id == "IG"
        else "unconfirmed_for_source_machine"
    )
    expected_machine = {
        "source_machine_id": machine_id,
        "source_machine_label": machine_label,
        "source_machine_interpretation": interpretation,
    }
    for field, expected in expected_machine.items():
        if str(manifest.get(field) or "") != expected:
            raise ValueError(f"Target manifest field {field!r} does not match the review queue.")
    expected_mz = "" if target_mz is None else f"{target_mz:.8g}"
    if str(manifest.get("target_mz") or "") != expected_mz:
        raise ValueError("Target manifest target_mz does not match the review queue.")
    if str(manifest.get("supplied_peak1_rt_min") or "") != str(peak1 or ""):
        raise ValueError("Target manifest Peak1 does not match the review queue.")
    if str(manifest.get("supplied_peak2_rt_min") or "") != str(peak2 or ""):
        raise ValueError("Target manifest Peak2 does not match the review queue.")

    payload = f"{render_context['dataset_id']}|{compound}|{pool}|{target_mz}|{peak1}|{peak2}"
    expected_uid = hashlib.sha256(payload.encode()).hexdigest()[:12]
    uid = str(manifest.get("target_uid") or "")
    if uid != expected_uid:
        raise ValueError("Target manifest target_uid does not match its source identity fields.")
    return uid


def _number(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _validate_image_payload(data: bytes, suffix: str, relative: str) -> tuple[int, int] | None:
    if suffix == ".png":
        return _validate_png_payload(data, relative)
    if b"<!doctype" in data.lower():
        raise ValueError(f"Canonical SVG must not contain a document type: {relative}")
    try:
        root = ElementTree.fromstring(data)
    except ElementTree.ParseError as error:
        raise ValueError(f"Canonical SVG is not valid XML: {relative}") from error
    if root.tag.rsplit("}", 1)[-1] != "svg":
        raise ValueError(f"Canonical SVG root element is not <svg>: {relative}")
    _validate_static_svg(root, relative)
    return None


def _validate_png_payload(data: bytes, relative: str) -> tuple[int, int]:
    if not data.startswith(b"\x89PNG\r\n\x1a\n"):
        raise ValueError(f"Canonical PNG has an invalid signature: {relative}")
    offset = 8
    width = height = bit_depth = color_type = None
    idat_parts: list[bytes] = []
    seen_ihdr = seen_iend = idat_closed = False
    while offset < len(data):
        if offset + 12 > len(data):
            raise ValueError(f"Canonical PNG has a truncated chunk: {relative}")
        length = int.from_bytes(data[offset : offset + 4], "big")
        chunk_type = data[offset + 4 : offset + 8]
        chunk_end = offset + 12 + length
        if chunk_end > len(data):
            raise ValueError(f"Canonical PNG has an invalid chunk length: {relative}")
        payload = data[offset + 8 : offset + 8 + length]
        expected_crc = int.from_bytes(data[offset + 8 + length : chunk_end], "big")
        if zlib.crc32(chunk_type + payload) & 0xFFFFFFFF != expected_crc:
            raise ValueError(f"Canonical PNG chunk CRC mismatch: {relative}")
        if not seen_ihdr and chunk_type != b"IHDR":
            raise ValueError(f"Canonical PNG must begin with IHDR: {relative}")
        if chunk_type == b"IHDR":
            if seen_ihdr or length != 13:
                raise ValueError(f"Canonical PNG has an invalid IHDR chunk: {relative}")
            seen_ihdr = True
            width = int.from_bytes(payload[0:4], "big")
            height = int.from_bytes(payload[4:8], "big")
            bit_depth, color_type, compression, filter_method, interlace = payload[8:13]
            valid_depths = {
                0: {1, 2, 4, 8, 16},
                2: {8, 16},
                3: {1, 2, 4, 8},
                4: {8, 16},
                6: {8, 16},
            }
            if (
                width <= 0
                or height <= 0
                or color_type not in valid_depths
                or bit_depth not in valid_depths[color_type]
                or compression != 0
                or filter_method != 0
                or interlace != 0
            ):
                raise ValueError(f"Canonical PNG has unsupported IHDR values: {relative}")
        elif chunk_type == b"IDAT":
            if not seen_ihdr or idat_closed:
                raise ValueError(f"Canonical PNG has invalid IDAT ordering: {relative}")
            idat_parts.append(payload)
        elif chunk_type == b"IEND":
            if length != 0 or seen_iend:
                raise ValueError(f"Canonical PNG has an invalid IEND chunk: {relative}")
            seen_iend = True
            offset = chunk_end
            break
        else:
            if idat_parts:
                idat_closed = True
            if chunk_type[:1].isupper() and chunk_type not in {b"PLTE"}:
                raise ValueError(f"Canonical PNG has an unknown critical chunk: {relative}")
        offset = chunk_end
    if not seen_ihdr or not seen_iend or offset != len(data) or not idat_parts:
        raise ValueError(f"Canonical PNG is incomplete or has trailing data: {relative}")
    assert width is not None and height is not None
    assert bit_depth is not None and color_type is not None
    channels = {0: 1, 2: 3, 3: 1, 4: 2, 6: 4}[color_type]
    row_bytes = (width * channels * bit_depth + 7) // 8
    expected_size = (row_bytes + 1) * height
    if expected_size > 256 * 1024 * 1024:
        raise ValueError(f"Canonical PNG decoded payload is too large: {relative}")
    decompressor = zlib.decompressobj()
    try:
        decoded = decompressor.decompress(b"".join(idat_parts), expected_size + 1)
        if decompressor.unconsumed_tail:
            raise ValueError(f"Canonical PNG decoded payload exceeds IHDR dimensions: {relative}")
        decoded += decompressor.flush()
    except zlib.error as error:
        raise ValueError(f"Canonical PNG IDAT stream is invalid: {relative}") from error
    if (
        not decompressor.eof
        or decompressor.unused_data
        or len(decoded) != expected_size
        or any(decoded[row * (row_bytes + 1)] > 4 for row in range(height))
    ):
        raise ValueError(f"Canonical PNG scanlines do not match IHDR dimensions: {relative}")
    return width, height


def _validate_static_svg(root: ElementTree.Element, relative: str) -> None:
    forbidden_elements = {"script", "foreignobject", "iframe", "object", "embed"}
    unsafe_tokens = ("javascript:", "data:", "@import", "expression(")
    for element in root.iter():
        tag = element.tag.rsplit("}", 1)[-1].casefold()
        if tag in forbidden_elements:
            raise ValueError(f"Canonical SVG contains active content: {relative}")
        if tag == "style":
            style = str(element.text or "").casefold()
            if any(token in style for token in (*unsafe_tokens, "url(")):
                raise ValueError(f"Canonical SVG style contains external content: {relative}")
        for raw_name, raw_value in element.attrib.items():
            name = raw_name.rsplit("}", 1)[-1].casefold()
            value = str(raw_value).strip().casefold()
            if name.startswith("on") or any(token in value for token in unsafe_tokens):
                raise ValueError(f"Canonical SVG contains an unsafe attribute: {relative}")
            if name in {"href", "src"} and value and not value.startswith("#"):
                raise ValueError(f"Canonical SVG contains an external reference: {relative}")
            if "url(" in value and not (value.startswith("url(#") and value.endswith(")")):
                raise ValueError(f"Canonical SVG contains an external URL: {relative}")


def _safe_relative_path(value: str, suffix: str) -> str:
    candidate = PurePosixPath(str(value).replace("\\", "/"))
    if not value or candidate.is_absolute() or ".." in candidate.parts:
        raise ValueError(f"Unsafe review asset path: {value!r}")
    normalised = candidate.as_posix()
    if not candidate.parts or candidate.parts[0] != "assets":
        raise ValueError(f"Canonical review assets must be under assets/: {normalised}")
    if candidate.suffix.lower() != suffix:
        raise ValueError(f"Review asset must end with {suffix}: {normalised}")
    return normalised


def _html_image_references(data: bytes) -> set[str]:
    parser = _ImageReferenceParser()
    parser.feed(data.decode("utf-8"))
    references = {_safe_html_reference(reference) for reference in parser.references}
    return references


def _safe_html_reference(value: str) -> str:
    candidate = PurePosixPath(value.replace("\\", "/"))
    if candidate.is_absolute() or ".." in candidate.parts:
        raise ValueError(f"Unsafe image reference in review HTML: {value!r}")
    return candidate.as_posix()


def _validate_cross_file_counts(
    summary: dict[str, Any],
    result_rows: list[dict[str, str]],
    manifest_rows: list[dict[str, str]],
    assets: dict[str, bytes],
) -> None:
    disagreements = summary.get("disagreements")
    agreements = summary.get("agreements")
    if (
        isinstance(disagreements, bool)
        or not isinstance(disagreements, int)
        or disagreements < 0
        or disagreements != len(manifest_rows)
    ):
        raise ValueError("Review manifest count does not match summary disagreements.")
    if (
        isinstance(agreements, bool)
        or not isinstance(agreements, int)
        or agreements < 0
        or agreements + disagreements != len(result_rows)
    ):
        raise ValueError("Review summary agreement/disagreement counts are inconsistent.")
    png_count = sum(path.endswith(".png") for path in assets)
    svg_count = sum(path.endswith(".svg") for path in assets)
    if png_count != len(manifest_rows) or svg_count != len(manifest_rows):
        raise ValueError("Manual-review package requires one canonical PNG and SVG per row.")


def _validate_review_tables(
    result_fields: list[str],
    result_rows: list[dict[str, str]],
    manifest_fields: list[str],
    manifest_rows: list[dict[str, str]],
    queue_fields: list[str],
    queue_rows: list[dict[str, str]],
    label_fields: list[str],
    label_rows: list[dict[str, str]],
) -> None:
    if queue_fields != [*result_fields, *SHADOW_COMPARISON_FIELDS]:
        raise ValueError(
            "Review queue header must be the exact final-result header plus comparison fields."
        )
    if len(manifest_fields) != len(set(manifest_fields)):
        raise ValueError("Target manifest contains duplicate column names.")
    required_manifest_fields = {
        "target_uid",
        "artifact_sha256",
        "png_sha256",
        "dataset_id",
        "run_id",
        "standard_id",
        "config_hash",
        "source_row",
        "compound_id",
        "target_mz",
        "supplied_peak1_rt_min",
        "supplied_peak2_rt_min",
        "source_machine_id",
        "source_machine_label",
        "source_machine_interpretation",
        "source_pool_id",
        "source_pooled_well",
        "svg_asset_path",
        "png_asset_path",
        "ms2_diagnostic_status",
        *ML_OUTPUT_FIELDS,
    }
    missing_manifest = sorted(required_manifest_fields - set(manifest_fields))
    if missing_manifest:
        raise ValueError(
            "Target manifest is missing required fields: " + ", ".join(missing_manifest)
        )
    if tuple(label_fields) != _REVIEW_LABEL_FIELDS:
        raise ValueError("Review label template header does not match the fixed schema.")

    for row in result_rows:
        v2_status = str(row.get("ms2_diagnostic_status") or "").strip()
        ml_status = str(row.get("ml_diagnostic_status") or "").strip()
        if v2_status not in V2_DIAGNOSTIC_STATUSES:
            raise ValueError(f"Invalid v2 diagnostic status in final result: {v2_status!r}.")
        if ml_status not in ML_DIAGNOSTIC_STATUSES:
            raise ValueError(f"Invalid ML diagnostic status in final result: {ml_status!r}.")
    expected_queue = [
        row
        for row in result_rows
        if _normalise_v2_status(str(row.get("ms2_diagnostic_status") or ""))
        != str(row.get("ml_diagnostic_status") or "")
    ]
    if len(queue_rows) != len(manifest_rows) or len(queue_rows) != len(expected_queue):
        raise ValueError("Review queue and target manifest row counts do not match.")
    for expected, queue in zip(expected_queue, queue_rows, strict=True):
        if any(
            str(queue.get(field) or "") != str(expected.get(field) or "") for field in result_fields
        ):
            raise ValueError(
                "Review queue is not the exact ordered disagreement subset of the CSV."
            )
        expected_comparison = _comparison_values(expected)
        if any(
            str(queue.get(field) or "") != expected_comparison[field]
            for field in SHADOW_COMPARISON_FIELDS
        ):
            raise ValueError("Review queue comparison fields do not match recomputed values.")
    if not ({"Compound_ID", "SGC ID for Component"} & set(queue_fields)):
        raise ValueError("Review queue is missing a compound identifier column.")
    if len(label_rows) != len(manifest_rows):
        raise ValueError("Review label template and target manifest row counts do not match.")

    manual_fields = {
        "reviewer_id",
        "manual_ms2_identity_label",
        "confidence",
        "peak1_spectrum_quality",
        "peak2_spectrum_quality",
        "manual_issue_codes",
        "diagnostic_fragment_notes",
        "notes",
        "reviewed_utc",
    }

    seen_uids: set[str] = set()
    for index, (manifest, queue, label) in enumerate(
        zip(manifest_rows, queue_rows, label_rows, strict=True), start=1
    ):
        uid = str(manifest.get("target_uid") or "")
        artifact_hash = str(manifest.get("artifact_sha256") or "")
        if not uid or uid in seen_uids:
            raise ValueError("Target manifest contains a missing or duplicate target_uid.")
        seen_uids.add(uid)
        if len(artifact_hash) != 64:
            raise ValueError(f"Target manifest has an invalid artifact SHA-256 for {uid}.")
        if label.get("target_uid") != uid or label.get("artifact_sha256") != artifact_hash:
            raise ValueError("Review label template does not match target_uid/artifact_sha256.")
        populated = sorted(field for field in manual_fields if str(label.get(field) or "").strip())
        if populated:
            raise ValueError(
                "The packaged review label template must not contain prior manual values: "
                + ", ".join(populated)
            )
        try:
            source_row = int(str(manifest.get("source_row") or ""))
        except ValueError as error:
            raise ValueError(f"Target manifest has an invalid source_row for {uid}.") from error
        if source_row != index:
            raise ValueError("Target manifest source_row values must match review queue order.")
        queue_compound = str(
            queue.get("Compound_ID") or queue.get("SGC ID for Component") or ""
        ).strip()
        if queue_compound != str(manifest.get("compound_id") or "").strip():
            raise ValueError("Review queue compound IDs do not match the target manifest.")
        for field in ("ms2_diagnostic_status", *ML_OUTPUT_FIELDS):
            if str(queue.get(field) or "") != str(manifest.get(field) or ""):
                raise ValueError(f"Review queue field {field!r} does not match target manifest.")


def _comparison_values(row: dict[str, str]) -> dict[str, str]:
    v2_status = str(row.get("ms2_diagnostic_status") or "").strip()
    ml_status = str(row.get("ml_diagnostic_status") or "").strip()
    return {
        "comparison_change_level": "v2_ml_status",
        "comparison_baseline_status": v2_status,
        "comparison_v2_normalized_status": _normalise_v2_status(v2_status),
        "comparison_candidate_status": ml_status,
        "comparison_review_tier": "v2_ml_disagreement",
        "comparison_review_flags": str(row.get("ml_review_flags") or ""),
    }


def _normalise_v2_status(status: str) -> str:
    return "not_evaluable" if status == "not_attempted" else status


def _model_summary(
    metadata: dict[str, Any],
    summary: dict[str, Any],
    rows: list[dict[str, str]],
    source_parameters: dict[str, Any],
    model_artifact_path: str | Path | None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    inference = metadata.get("shadow_inference")
    if not isinstance(inference, dict) or not isinstance(inference.get("model"), dict):
        raise ValueError("Final result metadata is missing shadow_inference.model provenance.")
    model = inference["model"]
    trained = model.get("trained")
    model_id = model.get("model_id")
    if type(trained) is not bool or not isinstance(model_id, str) or not model_id:
        raise ValueError("Invalid shadow model provenance in final result metadata.")
    row_model_ids = {str(row.get("ml_model_id") or "") for row in rows}
    if rows and row_model_ids != {model_id}:
        raise ValueError("Final result rows do not match the metadata model ID.")
    _validate_row_model_outputs(rows)
    summary_model = summary.get("model")
    if not isinstance(summary_model, dict) or summary_model.get("model_id") != model_id:
        raise ValueError("Shadow-review summary does not match the final result model ID.")
    if not trained:
        if model_artifact_path is not None:
            raise ValueError("UNTRAINED results must not be packaged with a model artifact.")
        if model_id != "UNTRAINED":
            raise ValueError("Untrained shadow metadata must use model_id=UNTRAINED.")
        if model.get("path") not in {None, ""} or model.get("sha256") not in {None, ""}:
            raise ValueError("Untrained shadow metadata must not claim a model artifact.")
        if summary_model.get("artifact_attached") is not False:
            raise ValueError("Untrained review summary must record artifact_attached=false.")
        _validate_predictions(rows, None)
        return (
            {"model_id": model_id, "trained": False, "selected_model": None},
            {"type": "untrained_v2_default_thresholds", **source_parameters},
        )

    if model_id == "UNTRAINED":
        raise ValueError("Trained shadow metadata cannot use model_id=UNTRAINED.")
    if model_artifact_path is None:
        raise ValueError("Trained shadow packaging requires the validated model artifact.")
    from .shadow import load_model_artifact

    model_path = Path(model_artifact_path).resolve()
    artifact = load_model_artifact(model_path)
    decision_payload = artifact["trainable_decision"]
    artifact_hash = _sha256_bytes(model_path.read_bytes())
    if decision_payload.get("model_id") != model_id:
        raise ValueError("Attached model artifact does not match the result model ID.")
    if model.get("sha256") != artifact_hash:
        raise ValueError("Final result metadata does not match the attached model SHA-256.")
    if summary_model.get("artifact_attached") is not True:
        raise ValueError("Trained review summary must record artifact_attached=true.")
    if summary_model.get("sha256") != artifact_hash:
        raise ValueError("Shadow-review summary does not match the attached model SHA-256.")
    _validate_predictions(rows, artifact)
    selected = str(decision_payload["selected_model"])
    if selected == "threshold":
        effective = {
            "type": "trained_threshold",
            **dict(decision_payload["model"]["parameters"]),
        }
    else:
        fitted = decision_payload["model"]
        boundaries = decision_payload["decision_thresholds"]
        effective = {
            "type": "trained_l2_logistic_regression",
            "C": fitted["C"],
            "class_weight": fitted["class_weight"],
            "conflict_max_probability": boundaries["conflict_max_probability"],
            "support_min_probability": boundaries["support_min_probability"],
        }
    return (
        {
            "model_id": model_id,
            "trained": True,
            "selected_model": selected,
            "artifact_sha256": artifact_hash,
        },
        effective,
    )


def _validate_row_model_outputs(rows: list[dict[str, str]]) -> None:
    for row in rows:
        status = str(row.get("ml_diagnostic_status") or "")
        if status not in ML_DIAGNOSTIC_STATUSES:
            raise ValueError(f"Invalid ML diagnostic status in final result: {status!r}.")
        probability = str(row.get("ml_same_compound_probability") or "").strip()
        if not probability:
            continue
        try:
            numeric = float(probability)
        except ValueError as error:
            raise ValueError("Final result contains a non-numeric ML probability.") from error
        if not 0 <= numeric <= 1:
            raise ValueError("Final result ML probabilities must be in [0, 1].")


def _validate_predictions(rows: list[dict[str, str]], artifact: dict[str, Any] | None) -> None:
    from .shadow import predict_shadow

    for row in rows:
        expected = predict_shadow(row, artifact)
        if any(str(row.get(field) or "") != expected[field] for field in ML_OUTPUT_FIELDS):
            mode = "trained" if artifact is not None else "UNTRAINED"
            raise ValueError(f"Final result does not match recomputed {mode} shadow predictions.")


def _source_decision_parameters(metadata: dict[str, Any]) -> dict[str, Any]:
    parameters = metadata.get("parameters")
    if not isinstance(parameters, dict):
        raise ValueError("Final result metadata is missing analysis parameters.")
    keys = (
        "min_cosine",
        "min_matched_peaks",
        "min_explained_intensity",
        "min_entropy_similarity",
    )
    if any(key not in parameters for key in keys):
        raise ValueError("Final result metadata is missing downstream decision parameters.")
    return {key: parameters[key] for key in keys}


def _render_review_labels(
    manifest_rows: list[dict[str, str]], label_rows: list[dict[str, str]]
) -> bytes:
    fields = [
        "target_uid",
        "artifact_sha256",
        "compound_id",
        "png_asset_path",
        "svg_asset_path",
        "current_ms2_status",
        "ml_diagnostic_status",
        "truth_label",
        "confidence",
        "peak1_spectrum_quality",
        "peak2_spectrum_quality",
        "manual_issue_codes",
        "diagnostic_fragment_notes",
        "notes",
        "reviewer_id",
        "review_round",
        "reviewed_utc",
    ]
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for manifest, label in zip(manifest_rows, label_rows, strict=True):
        writer.writerow(
            {
                "target_uid": manifest["target_uid"],
                "artifact_sha256": manifest["artifact_sha256"],
                "compound_id": manifest.get("compound_id", ""),
                "png_asset_path": manifest.get("png_asset_path", ""),
                "svg_asset_path": manifest.get("svg_asset_path", ""),
                "current_ms2_status": manifest.get("ms2_diagnostic_status", ""),
                "ml_diagnostic_status": manifest.get("ml_diagnostic_status", ""),
                "truth_label": "",
                "confidence": label.get("confidence", ""),
                "peak1_spectrum_quality": label.get("peak1_spectrum_quality", ""),
                "peak2_spectrum_quality": label.get("peak2_spectrum_quality", ""),
                "manual_issue_codes": label.get("manual_issue_codes", ""),
                "diagnostic_fragment_notes": label.get("diagnostic_fragment_notes", ""),
                "notes": label.get("notes", ""),
                "reviewer_id": label.get("reviewer_id", ""),
                "review_round": label.get("review_round", "1"),
                "reviewed_utc": label.get("reviewed_utc", ""),
            }
        )
    return ("\ufeff" + output.getvalue()).encode("utf-8")


def _render_readme(
    result_name: str,
    result_rows: int,
    review_rows: int,
    model: dict[str, Any],
    decision: dict[str, Any],
) -> bytes:
    if model["trained"]:
        training_note = (
            f"结果来自通过artifact完整性校验的已训练模型 {model['model_id']}，"
            f"选择的模型类型为{model['selected_model']}。"
        )
    else:
        training_note = (
            "当前数据没有独立ground truth；本包为UNTRAINED feasibility结果，"
            "下列阈值是当前经验默认值，不是训练得到的最佳参数。"
        )
    decision_lines = _decision_readme_lines(decision)
    text = f"""MS2人工复核包
==============

范围
----
- 最终结果：final_result/{result_name}（{result_rows}行）
- 可视化人工复核队列：{review_rows}条，每条包含PNG和SVG
- {training_note}

当前下游判定配置
----------------
{decision_lines}

使用方法
--------
1. 打开 START_HERE.html 浏览全部人工复核图。
2. 在 review_labels.csv 中填写truth_label、confidence、quality和notes；该表已包含
   target_uid、Compound_ID和对应PNG/SVG路径，不应依靠行顺序猜测映射。
3. review_queue.csv 保存自动状态、谱图指标和比较上下文，便于筛选排序。
4. final_result 目录保留完整最终CSV及其参数/provenance metadata。

人工标签词表
------------
- positive_same_compound
- negative_different_or_interference
- uncertain
- not_evaluable

科学边界
--------
这些图只评价两个保留时间位置的MS2证据是否支持同一化合物，不能直接证明
其为对映体，也不进行R/S指认。论文Method参数保持固定，不由本流程训练。
图片显式展示当前自动status、指标和阈值，因此本包适合audit/triage，不是盲法
ground-truth收集包；用于训练前仍需独立盲审和adjudication。
"""
    return text.encode("utf-8-sig")


def _decision_readme_lines(decision: dict[str, Any]) -> str:
    kind = decision["type"]
    if kind in {"untrained_v2_default_thresholds", "trained_threshold"}:
        entropy = decision["min_entropy_similarity"]
        entropy_text = "关闭" if entropy is None else str(entropy)
        prefix = "经验默认阈值" if kind.startswith("untrained") else "训练选定阈值"
        return "\n".join(
            (
                f"- 类型：{prefix}",
                f"- cosine >= {decision['min_cosine']}",
                f"- matched fragments >= {decision['min_matched_peaks']}",
                f"- 双侧 explained intensity >= {decision['min_explained_intensity']}",
                f"- entropy threshold：{entropy_text}",
            )
        )
    class_weight = decision["class_weight"]
    return "\n".join(
        (
            "- 类型：训练选定的StandardScaler + L2 Logistic Regression",
            f"- C = {decision['C']}",
            f"- class_weight = {class_weight}",
            f"- conflict probability <= {decision['conflict_max_probability']}",
            f"- support probability >= {decision['support_min_probability']}",
        )
    )


def _render_summary(
    result_path: Path,
    result_bytes: bytes,
    result_rows: int,
    review_rows: int,
    assets: dict[str, bytes],
    model: dict[str, Any],
    decision: dict[str, Any],
) -> bytes:
    value = {
        "schema_version": PACKAGE_SCHEMA_VERSION,
        "final_result": {
            "name": result_path.name,
            "rows": result_rows,
            "sha256": _sha256_bytes(result_bytes),
        },
        "review": {
            "rows": review_rows,
            "png_images": sum(path.endswith(".png") for path in assets),
            "svg_images": sum(path.endswith(".svg") for path in assets),
        },
        "model": model,
        "decision_configuration": decision,
    }
    return json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True).encode("utf-8")


def _render_hash_manifest(entries: dict[str, bytes]) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(["relative_path", "size_bytes", "sha256"])
    for path, data in sorted(entries.items()):
        writer.writerow([path, len(data), _sha256_bytes(data)])
    return ("\ufeff" + output.getvalue()).encode("utf-8")


def _validate_entry_names(entries: dict[str, bytes]) -> None:
    lowered: set[str] = set()
    for name in entries:
        candidate = PurePosixPath(name)
        if candidate.is_absolute() or ".." in candidate.parts:
            raise ValueError(f"Unsafe archive path: {name!r}")
        folded = name.casefold()
        if folded in lowered:
            raise ValueError(f"Case-insensitive duplicate archive path: {name}")
        lowered.add(folded)


def _write_archive(path: Path, entries: dict[str, bytes]) -> None:
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for relative, data in sorted(entries.items()):
            info = zipfile.ZipInfo(f"{PACKAGE_ROOT}/{relative}", (1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, data, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)


def _validate_archive(path: Path, entries: dict[str, bytes], html_references: set[str]) -> None:
    expected = {f"{PACKAGE_ROOT}/{relative}" for relative in entries}
    with zipfile.ZipFile(path) as archive:
        if archive.testzip() is not None:
            raise ValueError("ZIP CRC validation failed.")
        names = archive.namelist()
        if set(names) != expected or len(names) != len(expected):
            raise ValueError("ZIP inventory does not match the selected review files.")
        if len({name.casefold() for name in names}) != len(names):
            raise ValueError("ZIP contains case-insensitive duplicate paths.")
        for relative, data in entries.items():
            if archive.read(f"{PACKAGE_ROOT}/{relative}") != data:
                raise ValueError(f"ZIP payload differs from source: {relative}")
        for reference in html_references:
            if f"{PACKAGE_ROOT}/{reference}" not in expected:
                raise ValueError(f"START_HERE.html has a broken image reference: {reference}")


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


__all__ = ["PACKAGE_SCHEMA_VERSION", "package_shadow_review"]
