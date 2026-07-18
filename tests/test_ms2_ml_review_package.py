import base64
import csv
import hashlib
import json
import tempfile
import unittest
import zipfile
from contextlib import redirect_stdout
from io import BytesIO, StringIO
from pathlib import Path
from unittest.mock import patch

from MSAI.python.cli.ms2 import main as ms2_main
from MSAI.python.ms2_ml.contracts import (
    FROZEN_SOURCE_PARAMETERS,
    ML_OUTPUT_FIELDS,
    SHADOW_COMPARISON_FIELDS,
)
from MSAI.python.ms2_ml.provenance import canonical_json_sha256, canonical_table_sha256
from MSAI.python.ms2_ml.review_html import render_shadow_disagreement_html
from MSAI.python.ms2_ml.review_package import PACKAGE_ROOT, package_shadow_review
from MSAI.python.ms2_review.acquisition import reconcile_acquisition
from MSAI.python.ms2_review.model import prepare_mirror_spectrum
from MSAI.python.ms2_review.png import render_ms2_review_png
from MSAI.python.ms2_review.render_context import (
    build_review_render_metadata,
    review_config_hash,
)
from MSAI.python.ms2_review.svg import render_ms2_review_svg
from MSAI.python.source_context import source_machine_label_interpretation

PNG_1X1 = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII="
)
STANDARD_PATH = (
    Path(__file__).resolve().parents[1] / "MSAI" / "standards" / "ms2_diagnostic_standard_v2.json"
)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _csv_bytes(fields: list[str], rows: list[dict[str, str]], *, bom: bool = False) -> bytes:
    output = StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
    return (("\ufeff" if bom else "") + output.getvalue()).encode("utf-8")


def _fixture(root: Path) -> tuple[Path, Path]:
    result = root / "shadow.csv"
    source_fields = [
        "Compound_ID",
        "SGC ID for Pool",
        "MZ",
        "Peak1",
        "Peak2",
        "peak_a_MS2",
        "peak_b_MS2",
        "ms2_diagnostic_status",
    ]
    result_fields = [*source_fields, *ML_OUTPUT_FIELDS]
    source_row = {
        "Compound_ID": "CMP-1",
        "SGC ID for Pool": "POOL-1",
        "MZ": "100.123456",
        "Peak1": "1.25",
        "Peak2": "1.75",
        "peak_a_MS2": "100,1000;120,900",
        "peak_b_MS2": "100.002,950;120.002,850",
        "ms2_diagnostic_status": "supported_same_compound",
    }
    result_row = {
        **source_row,
        "ml_same_compound_probability": "",
        "ml_diagnostic_status": "not_evaluable",
        "ml_model_id": "UNTRAINED",
        "ml_abstention_reason": "NO_INDEPENDENT_TRUTH_MODEL",
        "ml_review_flags": "SHADOW_ONLY;FEASIBILITY_ONLY;NO_INDEPENDENT_TRUTH",
    }
    _write_csv(result, result_fields, [result_row])

    parameters = {
        **FROZEN_SOURCE_PARAMETERS,
        "min_cosine": 0.7,
        "min_matched_peaks": 6,
        "min_explained_intensity": 0.5,
        "min_entropy_similarity": None,
    }
    source_metadata = {"parameters": parameters}
    source_sidecar_bytes = json.dumps(source_metadata, indent=2, ensure_ascii=False).encode("utf-8")
    metadata = {
        **source_metadata,
        "shadow_inference": {
            "schema_version": "msai-ms2-shadow-output-v1",
            "source_result_sha256": _sha256(_csv_bytes(source_fields, [source_row])),
            "source_result_content_sha256": canonical_table_sha256(source_fields, [source_row]),
            "source_sidecar_sha256": _sha256(source_sidecar_bytes),
            "source_sidecar_content_sha256": canonical_json_sha256(source_metadata),
            "diagnostic_standard": str(STANDARD_PATH),
            "diagnostic_standard_sha256": _sha256(STANDARD_PATH.read_bytes()),
            "output_result_sha256": _sha256(result.read_bytes()),
            "append_only_fields": list(ML_OUTPUT_FIELDS),
            "rows": 1,
            "status_counts": {"not_evaluable": 1},
            "model": {
                "model_id": "UNTRAINED",
                "trained": False,
                "path": None,
                "sha256": None,
            },
        },
    }
    sidecar_path = Path(str(result) + ".metadata.json")
    sidecar_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

    review = root / "review"
    dataset_id = "v2_ml_disagreements"
    run_id = dataset_id
    thresholds = {
        key: parameters[key]
        for key in (
            "min_cosine",
            "min_matched_peaks",
            "min_explained_intensity",
            "min_entropy_similarity",
        )
    }
    model_context = {
        "path": None,
        "sha256": None,
        "model_id": "UNTRAINED",
    }
    acquisition_reconciliation = reconcile_acquisition(
        {},
        {},
        analysis_parameters=parameters,
        dia_windows=[],
        raw_input_count=0,
    )
    config_hash = review_config_hash(
        standard_id="MSAI-MS2-DIAGNOSTIC-v2",
        parameters=parameters,
        thresholds=thresholds,
        method_profile={},
        model_context=model_context,
    )
    target_mz = float(source_row["MZ"])
    uid_payload = (
        f"{dataset_id}|CMP-1|POOL-1|{target_mz}|{source_row['Peak1']}|{source_row['Peak2']}"
    )
    target_uid = hashlib.sha256(uid_payload.encode()).hexdigest()[:12]
    stem = f"CMP-1__{run_id}__mz{target_mz:.5f}__t-{target_uid}__cfg-{config_hash}"
    png_relative = f"assets/run/png/{stem}.png"
    svg_relative = f"assets/run/svg/{stem}.svg"
    png_path = review / png_relative
    svg_path = review / svg_relative
    png_path.parent.mkdir(parents=True)
    svg_path.parent.mkdir(parents=True)
    queue_row = {
        **result_row,
        "comparison_change_level": "v2_ml_status",
        "comparison_baseline_status": "supported_same_compound",
        "comparison_v2_normalized_status": "supported_same_compound",
        "comparison_candidate_status": "not_evaluable",
        "comparison_review_tier": "v2_ml_disagreement",
        "comparison_review_flags": result_row["ml_review_flags"],
    }
    render_metadata = build_review_render_metadata(
        target_uid=target_uid,
        config_hash=config_hash,
        standard_id="MSAI-MS2-DIAGNOSTIC-v2",
        source_machine_id="IG",
        source_machine_label="",
        source_machine_interpretation=source_machine_label_interpretation(""),
        source_pool_id=source_row["SGC ID for Pool"],
        source_pooled_well="",
        parameters=parameters,
        thresholds=thresholds,
        dia_windows=[],
        acquisition_context={},
        method_profile={},
        acquisition_reconciliation=acquisition_reconciliation,
        model_context=model_context,
    )
    mirror = prepare_mirror_spectrum(
        queue_row,
        fragment_mz_tol=float(render_metadata["fragment_mz_tol"]),
        fragment_mz_tol_unit=str(render_metadata["fragment_mz_tol_unit"]),
        min_relative_intensity=float(render_metadata["min_relative_intensity"]),
    )
    svg_bytes = render_ms2_review_svg(
        row=queue_row,
        mirror=mirror,
        metadata=render_metadata,
    ).encode("utf-8")
    svg_path.write_bytes(svg_bytes)
    png_image = render_ms2_review_png(
        row=queue_row,
        mirror=mirror,
        metadata=render_metadata,
        scale=0.25,
    )
    png_buffer = BytesIO()
    png_image.save(png_buffer, format="PNG", optimize=False, compress_level=3)
    png_bytes = png_buffer.getvalue()
    png_path.write_bytes(png_bytes)
    queue_path = review / "v2_ml_disagreements.csv"
    _write_csv(
        queue_path,
        [*result_fields, *SHADOW_COMPARISON_FIELDS],
        [queue_row],
    )
    _write_csv(
        review / "annotations/_template/review_labels.csv",
        [
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
        ],
        [
            {
                "target_uid": target_uid,
                "artifact_sha256": _sha256(svg_bytes),
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
        ],
    )
    manifest_fields = [
        "target_uid",
        "source_row",
        "compound_id",
        "artifact_sha256",
        "png_sha256",
        "svg_asset_path",
        "png_asset_path",
        "dataset_id",
        "run_id",
        "standard_id",
        "config_hash",
        "target_mz",
        "supplied_peak1_rt_min",
        "supplied_peak2_rt_min",
        "source_machine_id",
        "source_machine_label",
        "source_machine_interpretation",
        "source_pool_id",
        "source_pooled_well",
        "ms2_diagnostic_status",
        "ml_same_compound_probability",
        "ml_diagnostic_status",
        "ml_model_id",
        "ml_abstention_reason",
        "ml_review_flags",
    ]
    manifest_row = {
        "target_uid": target_uid,
        "source_row": "1",
        "compound_id": "CMP-1",
        "artifact_sha256": _sha256(svg_bytes),
        "png_sha256": _sha256(png_bytes),
        "svg_asset_path": svg_relative,
        "png_asset_path": png_relative,
        "dataset_id": dataset_id,
        "run_id": run_id,
        "standard_id": "MSAI-MS2-DIAGNOSTIC-v2",
        "config_hash": config_hash,
        "target_mz": f"{target_mz:.8g}",
        "supplied_peak1_rt_min": source_row["Peak1"],
        "supplied_peak2_rt_min": source_row["Peak2"],
        "source_machine_id": "IG",
        "source_machine_label": "",
        "source_machine_interpretation": source_machine_label_interpretation(""),
        "source_pool_id": source_row["SGC ID for Pool"],
        "source_pooled_well": "",
        **{
            key: queue_row[key]
            for key in (
                "ms2_diagnostic_status",
                "ml_same_compound_probability",
                "ml_diagnostic_status",
                "ml_model_id",
                "ml_abstention_reason",
                "ml_review_flags",
            )
        },
    }
    _write_csv(
        review / "metadata/target_manifest.csv",
        manifest_fields,
        [manifest_row],
    )
    (review / "v2_ml_disagreements.html").write_text(
        render_shadow_disagreement_html([queue_row], [manifest_row]), encoding="utf-8"
    )
    summary = {
        "schema_version": "msai-ms2-shadow-disagreement-review-v1",
        "source_shadow_sha256": _sha256(result.read_bytes()),
        "standard": {
            "standard_id": "MSAI-MS2-DIAGNOSTIC-v2",
            "path": str(STANDARD_PATH),
            "sha256": _sha256(STANDARD_PATH.read_bytes()),
        },
        "rows": 1,
        "agreements": 0,
        "disagreements": 1,
        "transition_counts": {"supported_same_compound->not_evaluable": 1},
        "model": {
            "model_id": "UNTRAINED",
            "artifact_attached": False,
            "path": None,
            "sha256": None,
        },
    }
    (review / "shadow_summary.json").write_text(json.dumps(summary), encoding="utf-8")
    run_manifest = {
        "dataset_id": dataset_id,
        "run_id": run_id,
        "standard_id": "MSAI-MS2-DIAGNOSTIC-v2",
        "config_hash": config_hash,
        "source_result_sha256": _sha256(queue_path.read_bytes()),
        "source_sidecar_sha256": _sha256(sidecar_path.read_bytes()),
        "parameters": parameters,
        "thresholds": thresholds,
        "dia_windows": [],
        "acquisition_context": {},
        "method_profile": {},
        "acquisition_reconciliation": acquisition_reconciliation,
        "image_formats": ["svg", "png"],
        "png_scale": 0.25,
        "source_machine_label_column": "IG",
        "preserve_input_diagnostics": True,
        "model_artifact": model_context,
    }
    (review / "metadata/run_manifest.json").write_text(json.dumps(run_manifest), encoding="utf-8")
    return result, review


def _read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def _sync_review_provenance(result: Path, review: Path) -> None:
    result_fields, result_rows = _read_csv(result)
    source_fields = result_fields[: -len(ML_OUTPUT_FIELDS)]
    source_rows = [{field: row.get(field, "") for field in source_fields} for row in result_rows]
    sidecar_path = Path(str(result) + ".metadata.json")
    metadata = json.loads(sidecar_path.read_text(encoding="utf-8"))
    inference = metadata["shadow_inference"]
    inference["source_result_sha256"] = _sha256(_csv_bytes(source_fields, source_rows))
    inference["source_result_content_sha256"] = canonical_table_sha256(source_fields, source_rows)
    source_metadata = dict(metadata)
    source_metadata.pop("shadow_inference")
    inference["source_sidecar_sha256"] = _sha256(
        json.dumps(source_metadata, indent=2, ensure_ascii=False).encode("utf-8")
    )
    inference["source_sidecar_content_sha256"] = canonical_json_sha256(source_metadata)
    inference["output_result_sha256"] = _sha256(result.read_bytes())
    inference["rows"] = len(result_rows)
    status_counts: dict[str, int] = {}
    for row in result_rows:
        status = row["ml_diagnostic_status"]
        status_counts[status] = status_counts.get(status, 0) + 1
    inference["status_counts"] = status_counts
    sidecar_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

    queue_rows = []
    for row in result_rows:
        v2 = row["ms2_diagnostic_status"]
        normalized = "not_evaluable" if v2 == "not_attempted" else v2
        ml = row["ml_diagnostic_status"]
        if normalized == ml:
            continue
        queue_rows.append(
            {
                **row,
                "comparison_change_level": "v2_ml_status",
                "comparison_baseline_status": v2,
                "comparison_v2_normalized_status": normalized,
                "comparison_candidate_status": ml,
                "comparison_review_tier": "v2_ml_disagreement",
                "comparison_review_flags": row["ml_review_flags"],
            }
        )
    queue_path = review / "v2_ml_disagreements.csv"
    _write_csv(
        queue_path,
        [*result_fields, *SHADOW_COMPARISON_FIELDS],
        queue_rows,
    )

    manifest_path = review / "metadata/target_manifest.csv"
    manifest_fields, manifest_rows = _read_csv(manifest_path)
    for manifest, queue in zip(manifest_rows, queue_rows, strict=True):
        for field in ("ms2_diagnostic_status", *ML_OUTPUT_FIELDS):
            manifest[field] = queue[field]
    _write_csv(manifest_path, manifest_fields, manifest_rows)
    (review / "v2_ml_disagreements.html").write_text(
        render_shadow_disagreement_html(queue_rows, manifest_rows), encoding="utf-8"
    )

    summary_path = review / "shadow_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["source_shadow_sha256"] = _sha256(result.read_bytes())
    summary["rows"] = len(result_rows)
    summary["disagreements"] = len(queue_rows)
    summary["agreements"] = len(result_rows) - len(queue_rows)
    transitions: dict[str, int] = {}
    for row in result_rows:
        v2 = row["ms2_diagnostic_status"]
        normalized = "not_evaluable" if v2 == "not_attempted" else v2
        key = f"{normalized}->{row['ml_diagnostic_status']}"
        transitions[key] = transitions.get(key, 0) + 1
    summary["transition_counts"] = dict(sorted(transitions.items()))
    model = inference["model"]
    summary["model"] = {
        "model_id": model["model_id"],
        "artifact_attached": bool(model["trained"]),
        "path": model.get("path"),
        "sha256": model.get("sha256"),
    }
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    run_path = review / "metadata/run_manifest.json"
    run_manifest = json.loads(run_path.read_text(encoding="utf-8"))
    run_manifest["source_result_sha256"] = _sha256(queue_path.read_bytes())
    run_manifest["source_sidecar_sha256"] = _sha256(sidecar_path.read_bytes())
    run_manifest["parameters"] = metadata["parameters"]
    run_manifest["thresholds"] = {
        key: metadata["parameters"][key]
        for key in (
            "min_cosine",
            "min_matched_peaks",
            "min_explained_intensity",
            "min_entropy_similarity",
        )
    }
    run_manifest["model_artifact"] = {
        "path": model.get("path"),
        "sha256": model.get("sha256"),
        "model_id": model["model_id"],
    }
    run_manifest["config_hash"] = review_config_hash(
        standard_id=run_manifest["standard_id"],
        parameters=run_manifest["parameters"],
        thresholds=run_manifest["thresholds"],
        method_profile=run_manifest["method_profile"],
        model_context=run_manifest["model_artifact"],
    )
    for manifest, queue in zip(manifest_rows, queue_rows, strict=True):
        manifest["config_hash"] = run_manifest["config_hash"]
        machine_id = str(run_manifest["source_machine_label_column"] or "").upper()
        machine_label = str(queue.get(machine_id) or "").strip()
        interpretation = source_machine_label_interpretation(machine_label)
        manifest["source_machine_id"] = machine_id
        manifest["source_machine_label"] = machine_label
        manifest["source_machine_interpretation"] = interpretation
        manifest["source_pooled_well"] = str(queue.get("Pooled Well") or "")
        render_metadata = build_review_render_metadata(
            target_uid=manifest["target_uid"],
            config_hash=run_manifest["config_hash"],
            standard_id=run_manifest["standard_id"],
            source_machine_id=machine_id,
            source_machine_label=machine_label,
            source_machine_interpretation=interpretation,
            source_pool_id=str(queue.get("SGC ID for Pool") or ""),
            source_pooled_well=str(queue.get("Pooled Well") or ""),
            parameters=run_manifest["parameters"],
            thresholds=run_manifest["thresholds"],
            dia_windows=run_manifest["dia_windows"],
            acquisition_context=run_manifest["acquisition_context"],
            method_profile=run_manifest["method_profile"],
            acquisition_reconciliation=run_manifest["acquisition_reconciliation"],
            model_context=run_manifest["model_artifact"],
        )
        mirror = prepare_mirror_spectrum(
            queue,
            fragment_mz_tol=float(render_metadata["fragment_mz_tol"]),
            fragment_mz_tol_unit=str(render_metadata["fragment_mz_tol_unit"]),
            min_relative_intensity=float(render_metadata["min_relative_intensity"]),
        )
        svg_bytes = render_ms2_review_svg(
            row=queue,
            mirror=mirror,
            metadata=render_metadata,
        ).encode("utf-8")
        (review / manifest["svg_asset_path"]).write_bytes(svg_bytes)
        manifest["artifact_sha256"] = _sha256(svg_bytes)
        png_image = render_ms2_review_png(
            row=queue,
            mirror=mirror,
            metadata=render_metadata,
            scale=run_manifest["png_scale"],
        )
        png_buffer = BytesIO()
        png_image.save(png_buffer, format="PNG", optimize=False, compress_level=3)
        png_bytes = png_buffer.getvalue()
        (review / manifest["png_asset_path"]).write_bytes(png_bytes)
        manifest["png_sha256"] = _sha256(png_bytes)
    _write_csv(manifest_path, manifest_fields, manifest_rows)
    label_path = review / "annotations/_template/review_labels.csv"
    label_fields, label_rows = _read_csv(label_path)
    for manifest, label in zip(manifest_rows, label_rows, strict=True):
        label["artifact_sha256"] = manifest["artifact_sha256"]
    _write_csv(label_path, label_fields, label_rows)
    run_path.write_text(json.dumps(run_manifest), encoding="utf-8")


class ShadowReviewPackageTests(unittest.TestCase):
    def test_package_is_compact_integrity_checked_and_reproducible(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, review = _fixture(root)
            first = root / "first.zip"
            second = root / "second.zip"

            first_summary = package_shadow_review(result, review, first)
            second_summary = package_shadow_review(result, review, second)

            self.assertEqual(first_summary["sha256"], second_summary["sha256"])
            self.assertFalse(first_summary["trained"])
            self.assertEqual(first_summary["review_rows"], 1)
            self.assertEqual(first_summary["png_images"], 1)
            self.assertEqual(first_summary["svg_images"], 1)
            with zipfile.ZipFile(first) as archive:
                self.assertIsNone(archive.testzip())
                names = set(archive.namelist())
                self.assertIn(f"{PACKAGE_ROOT}/START_HERE.html", names)
                self.assertIn(f"{PACKAGE_ROOT}/final_result/shadow.csv", names)
                png_names = sorted(name for name in names if name.endswith(".png"))
                svg_names = sorted(name for name in names if name.endswith(".svg"))
                self.assertEqual(len(png_names), 1)
                self.assertEqual(len(svg_names), 1)
                self.assertIn("__t-", png_names[0])
                self.assertFalse(any("/views/" in name or "/galleries/" in name for name in names))
                readme = archive.read(f"{PACKAGE_ROOT}/README_FIRST.txt").decode("utf-8-sig")
                self.assertIn("UNTRAINED feasibility", readme)
                self.assertIn("cosine >= 0.7", readme)
                package_summary = json.loads(archive.read(f"{PACKAGE_ROOT}/package_summary.json"))
                self.assertEqual(package_summary["final_result"]["rows"], 1)
                self.assertEqual(package_summary["model"]["model_id"], "UNTRAINED")
                self.assertEqual(
                    package_summary["decision_configuration"]["type"],
                    "untrained_v2_default_thresholds",
                )
                with archive.open(f"{PACKAGE_ROOT}/review_labels.csv") as handle:
                    labels = list(
                        csv.DictReader(StringIO(handle.read().decode("utf-8-sig"), newline=""))
                    )
                self.assertEqual(labels[0]["compound_id"], "CMP-1")
                self.assertEqual(labels[0]["truth_label"], "")
                self.assertEqual(f"{PACKAGE_ROOT}/{labels[0]['png_asset_path']}", png_names[0])

    def test_corrupt_asset_is_rejected_before_output_is_written(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, review = _fixture(root)
            _fields, manifest_rows = _read_csv(review / "metadata/target_manifest.csv")
            (review / manifest_rows[0]["png_asset_path"]).write_bytes(PNG_1X1 + b"changed")
            output = root / "invalid.zip"

            with self.assertRaisesRegex(ValueError, "trailing data"):
                package_shadow_review(result, review, output)
            self.assertFalse(output.exists())

    def test_unsafe_asset_path_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, review = _fixture(root)
            manifest = review / "metadata/target_manifest.csv"
            with manifest.open(encoding="utf-8-sig", newline="") as handle:
                rows = list(csv.DictReader(handle))
            rows[0]["png_asset_path"] = "../CMP-1.png"
            _write_csv(manifest, list(rows[0]), rows)

            with self.assertRaisesRegex(ValueError, "Unsafe review asset path"):
                package_shadow_review(result, review, root / "unsafe.zip")

    def test_existing_output_requires_explicit_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, review = _fixture(root)
            output = root / "review.zip"
            initial = package_shadow_review(result, review, output)

            with self.assertRaises(FileExistsError):
                package_shadow_review(result, review, output)
            replaced = package_shadow_review(result, review, output, overwrite=True)
            self.assertEqual(initial["sha256"], replaced["sha256"])

    def test_queue_and_label_template_must_match_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, review = _fixture(root)
            queue = review / "v2_ml_disagreements.csv"
            with queue.open(encoding="utf-8-sig", newline="") as handle:
                fields = list(csv.DictReader(handle).fieldnames or [])
            _write_csv(queue, fields, [])
            with self.assertRaisesRegex(ValueError, "queue and target manifest row counts"):
                package_shadow_review(result, review, root / "queue-mismatch.zip")

            result, review = _fixture(root / "label-case")
            labels_path = review / "annotations/_template/review_labels.csv"
            with labels_path.open(encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                label_fields = list(reader.fieldnames or [])
                labels = list(reader)
            labels[0]["target_uid"] = "wrong-uid"
            _write_csv(labels_path, label_fields, labels)
            with self.assertRaisesRegex(ValueError, "target_uid/artifact_sha256"):
                package_shadow_review(result, review, root / "label-mismatch.zip")

    def test_broken_html_reference_and_symlink_escape_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, review = _fixture(root)
            html = review / "v2_ml_disagreements.html"
            _fields, manifest_rows = _read_csv(review / "metadata/target_manifest.csv")
            png_name = Path(manifest_rows[0]["png_asset_path"]).name
            html.write_text(
                html.read_text(encoding="utf-8").replace(png_name, "missing.png"),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "exact deterministic rendering"):
                package_shadow_review(result, review, root / "broken-html.zip")

            result, review = _fixture(root / "symlink-case")
            outside = root / "outside.png"
            outside.write_bytes(PNG_1X1)
            link = review / "assets/run/png/link.png"
            link.symlink_to(outside)
            manifest_path = review / "metadata/target_manifest.csv"
            with manifest_path.open(encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                manifest_fields = list(reader.fieldnames or [])
                manifest_rows = list(reader)
            manifest_rows[0]["png_asset_path"] = "assets/run/png/link.png"
            _write_csv(manifest_path, manifest_fields, manifest_rows)
            html_path = review / "v2_ml_disagreements.html"
            html_path.write_text(
                html_path.read_text(encoding="utf-8").replace("CMP-1.png", "link.png"),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "escapes its root"):
                package_shadow_review(result, review, root / "symlink.zip")

    def test_untrained_semantic_contract_is_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, review = _fixture(root)
            with result.open(encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                fields = list(reader.fieldnames or [])
                rows = list(reader)
            rows[0]["ml_review_flags"] = "SHADOW_ONLY;FEASIBILITY_ONLY"
            _write_csv(result, fields, rows)
            metadata_path = Path(str(result) + ".metadata.json")
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            metadata["shadow_inference"]["output_result_sha256"] = _sha256(result.read_bytes())
            metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
            for path in (
                review / "v2_ml_disagreements.csv",
                review / "metadata/target_manifest.csv",
            ):
                with path.open(encoding="utf-8-sig", newline="") as handle:
                    reader = csv.DictReader(handle)
                    table_fields = list(reader.fieldnames or [])
                    table_rows = list(reader)
                table_rows[0]["ml_review_flags"] = rows[0]["ml_review_flags"]
                _write_csv(path, table_fields, table_rows)
            summary_path = review / "shadow_summary.json"
            self.assertTrue(summary_path.is_file())
            _sync_review_provenance(result, review)

            with self.assertRaisesRegex(ValueError, "recomputed UNTRAINED shadow predictions"):
                package_shadow_review(result, review, root / "bad-untrained.zip")

    def test_trained_package_uses_attached_model_decision_not_source_defaults(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, review = _fixture(root)
            with result.open(encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                fields = list(reader.fieldnames or [])
                rows = list(reader)
            rows[0].update(
                {
                    "ms2_diagnostic_status": "insufficient_evidence",
                    "ml_same_compound_probability": "0.91",
                    "ml_diagnostic_status": "supported_same_compound",
                    "ml_model_id": "MODEL-1",
                    "ml_abstention_reason": "",
                    "ml_review_flags": "SHADOW_ONLY;INTERPRETABLE_THRESHOLD_MODEL",
                }
            )
            _write_csv(result, fields, rows)
            model_path = root / "model.json"
            model_path.write_text("{}", encoding="utf-8")
            model_hash = _sha256(model_path.read_bytes())
            metadata_path = Path(str(result) + ".metadata.json")
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            metadata["shadow_inference"]["model"] = {
                "model_id": "MODEL-1",
                "trained": True,
                "path": str(model_path),
                "sha256": model_hash,
            }
            metadata["shadow_inference"]["output_result_sha256"] = _sha256(result.read_bytes())
            metadata["shadow_inference"]["status_counts"] = {"supported_same_compound": 1}
            metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
            queue_path = review / "v2_ml_disagreements.csv"
            with queue_path.open(encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                queue_fields = list(reader.fieldnames or [])
                queue_rows = list(reader)
            for field in (
                "ms2_diagnostic_status",
                "ml_same_compound_probability",
                "ml_diagnostic_status",
                "ml_model_id",
                "ml_abstention_reason",
                "ml_review_flags",
            ):
                queue_rows[0][field] = rows[0][field]
            _write_csv(queue_path, queue_fields, queue_rows)
            manifest_path = review / "metadata/target_manifest.csv"
            with manifest_path.open(encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                manifest_fields = list(reader.fieldnames or [])
                manifest_rows = list(reader)
            for field in (
                "ms2_diagnostic_status",
                "ml_same_compound_probability",
                "ml_diagnostic_status",
                "ml_model_id",
                "ml_abstention_reason",
                "ml_review_flags",
            ):
                manifest_rows[0][field] = rows[0][field]
            _write_csv(manifest_path, manifest_fields, manifest_rows)
            summary_path = review / "shadow_summary.json"
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            summary["model"] = {
                "model_id": "MODEL-1",
                "artifact_attached": True,
                "path": str(model_path),
                "sha256": model_hash,
            }
            summary_path.write_text(json.dumps(summary), encoding="utf-8")
            _sync_review_provenance(result, review)
            loaded = {
                "trainable_decision": {
                    "model_id": "MODEL-1",
                    "selected_model": "threshold",
                    "model": {
                        "parameters": {
                            "min_cosine": 0.85,
                            "min_matched_peaks": 8,
                            "min_explained_intensity": 0.6,
                            "min_entropy_similarity": 0.75,
                        }
                    },
                }
            }
            output = root / "trained.zip"
            predicted = {
                field: rows[0][field]
                for field in (
                    "ml_same_compound_probability",
                    "ml_diagnostic_status",
                    "ml_model_id",
                    "ml_abstention_reason",
                    "ml_review_flags",
                )
            }
            with (
                patch(
                    "MSAI.python.ms2_ml.shadow.load_model_artifact",
                    return_value=loaded,
                ),
                patch(
                    "MSAI.python.ms2_ml.shadow.predict_shadow",
                    return_value=predicted,
                ),
            ):
                package_shadow_review(
                    result,
                    review,
                    output,
                    model_artifact_path=model_path,
                )
            with zipfile.ZipFile(output) as archive:
                package_summary = json.loads(archive.read(f"{PACKAGE_ROOT}/package_summary.json"))
            self.assertEqual(package_summary["decision_configuration"]["min_cosine"], 0.85)
            self.assertEqual(package_summary["model"]["selected_model"], "threshold")

    def test_router_exposes_package_command_help(self):
        output = StringIO()
        with redirect_stdout(output), self.assertRaises(SystemExit) as raised:
            ms2_main(["package-shadow-review", "--help"])
        self.assertEqual(raised.exception.code, 0)
        self.assertIn("canonical PNG/SVG", output.getvalue())

    def test_result_header_requires_exact_ml_tail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, review = _fixture(root)
            fields, rows = _read_csv(result)
            fields.remove("ml_model_id")
            fields.insert(2, "ml_model_id")
            _write_csv(result, fields, rows)

            with self.assertRaisesRegex(ValueError, "exact five ML fields at its tail"):
                package_shadow_review(result, review, root / "bad-header.zip")

    def test_invalid_v2_status_and_forged_comparison_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, review = _fixture(root)
            fields, rows = _read_csv(result)
            rows[0]["ms2_diagnostic_status"] = "invented_status"
            _write_csv(result, fields, rows)
            _sync_review_provenance(result, review)

            with self.assertRaisesRegex(ValueError, "Invalid v2 diagnostic status"):
                package_shadow_review(result, review, root / "bad-status.zip")

            result, review = _fixture(root / "comparison-case")
            queue_path = review / "v2_ml_disagreements.csv"
            queue_fields, queue_rows = _read_csv(queue_path)
            queue_rows[0]["comparison_candidate_status"] = "conflicting_spectra"
            _write_csv(queue_path, queue_fields, queue_rows)
            with self.assertRaisesRegex(ValueError, "comparison fields"):
                package_shadow_review(result, review, root / "bad-comparison.zip")

    def test_source_sidecar_parameter_tampering_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, review = _fixture(root)
            sidecar_path = Path(str(result) + ".metadata.json")
            metadata = json.loads(sidecar_path.read_text(encoding="utf-8"))
            metadata["parameters"]["min_cosine"] = 0.99
            sidecar_path.write_text(
                json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
            )

            with self.assertRaisesRegex(ValueError, "source-sidecar content hash"):
                package_shadow_review(result, review, root / "tampered-sidecar.zip")

    def test_content_hashes_allow_noncanonical_original_byte_formatting(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, review = _fixture(root)
            sidecar_path = Path(str(result) + ".metadata.json")
            metadata = json.loads(sidecar_path.read_text(encoding="utf-8"))
            metadata["shadow_inference"]["source_result_sha256"] = "a" * 64
            metadata["shadow_inference"]["source_sidecar_sha256"] = "b" * 64
            sidecar_path.write_text(json.dumps(metadata), encoding="utf-8")
            run_path = review / "metadata/run_manifest.json"
            run_manifest = json.loads(run_path.read_text(encoding="utf-8"))
            run_manifest["source_sidecar_sha256"] = _sha256(sidecar_path.read_bytes())
            run_path.write_text(json.dumps(run_manifest), encoding="utf-8")

            summary = package_shadow_review(result, review, root / "content-hashes.zip")

            self.assertEqual(summary["result_rows"], 1)

    def test_html_body_and_target_uid_tampering_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, review = _fixture(root)
            html_path = review / "v2_ml_disagreements.html"
            html_path.write_text(
                html_path.read_text(encoding="utf-8").replace("CMP-1", "FORGED", 1),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "exact deterministic rendering"):
                package_shadow_review(result, review, root / "tampered-html.zip")

            result, review = _fixture(root / "uid-case")
            manifest_path = review / "metadata/target_manifest.csv"
            manifest_fields, manifest_rows = _read_csv(manifest_path)
            manifest_rows[0]["target_uid"] = "000000000000"
            _write_csv(manifest_path, manifest_fields, manifest_rows)
            label_path = review / "annotations/_template/review_labels.csv"
            label_fields, label_rows = _read_csv(label_path)
            label_rows[0]["target_uid"] = "000000000000"
            _write_csv(label_path, label_fields, label_rows)
            with self.assertRaisesRegex(ValueError, "target_uid does not match"):
                package_shadow_review(result, review, root / "tampered-uid.zip")

    def test_forged_or_active_image_payloads_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, review = _fixture(root / "png-case")
            manifest_path = review / "metadata/target_manifest.csv"
            manifest_fields, manifest_rows = _read_csv(manifest_path)
            forged_png = b"\x89PNG\r\n\x1a\nnot-a-png"
            (review / manifest_rows[0]["png_asset_path"]).write_bytes(forged_png)
            manifest_rows[0]["png_sha256"] = _sha256(forged_png)
            _write_csv(manifest_path, manifest_fields, manifest_rows)
            with self.assertRaisesRegex(ValueError, "truncated chunk"):
                package_shadow_review(result, review, root / "forged-png.zip")

            result, review = _fixture(root / "valid-but-wrong-png-case")
            manifest_path = review / "metadata/target_manifest.csv"
            manifest_fields, manifest_rows = _read_csv(manifest_path)
            png_path = review / manifest_rows[0]["png_asset_path"]
            from PIL import Image  # type: ignore[import-not-found]

            with Image.open(png_path) as source:
                forged_image = source.convert("RGB")
            forged_image.putpixel((0, 0), (255, 0, 0))
            forged_image.save(png_path, format="PNG", optimize=False, compress_level=3)
            manifest_rows[0]["png_sha256"] = _sha256(png_path.read_bytes())
            _write_csv(manifest_path, manifest_fields, manifest_rows)
            with self.assertRaisesRegex(ValueError, "PNG is not the deterministic rendering"):
                package_shadow_review(result, review, root / "wrong-pixels.zip")

            result, review = _fixture(root / "script-case")
            manifest_path = review / "metadata/target_manifest.csv"
            manifest_fields, manifest_rows = _read_csv(manifest_path)
            active_svg = b'<svg xmlns="http://www.w3.org/2000/svg"><script>1</script></svg>'
            (review / manifest_rows[0]["svg_asset_path"]).write_bytes(active_svg)
            manifest_rows[0]["artifact_sha256"] = _sha256(active_svg)
            _write_csv(manifest_path, manifest_fields, manifest_rows)
            label_path = review / "annotations/_template/review_labels.csv"
            label_fields, label_rows = _read_csv(label_path)
            label_rows[0]["artifact_sha256"] = _sha256(active_svg)
            _write_csv(label_path, label_fields, label_rows)
            with self.assertRaisesRegex(ValueError, "active content"):
                package_shadow_review(result, review, root / "active-svg.zip")

            result, review = _fixture(root / "unrelated-case")
            manifest_path = review / "metadata/target_manifest.csv"
            manifest_fields, manifest_rows = _read_csv(manifest_path)
            unrelated_svg = (
                b'<svg xmlns="http://www.w3.org/2000/svg" width="1240" height="1030"></svg>'
            )
            (review / manifest_rows[0]["svg_asset_path"]).write_bytes(unrelated_svg)
            manifest_rows[0]["artifact_sha256"] = _sha256(unrelated_svg)
            _write_csv(manifest_path, manifest_fields, manifest_rows)
            label_path = review / "annotations/_template/review_labels.csv"
            label_fields, label_rows = _read_csv(label_path)
            label_rows[0]["artifact_sha256"] = _sha256(unrelated_svg)
            _write_csv(label_path, label_fields, label_rows)
            with self.assertRaisesRegex(ValueError, "deterministic rendering"):
                package_shadow_review(result, review, root / "unrelated-svg.zip")

    def test_render_context_and_manifest_display_fields_are_sidecar_queue_derived(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, review = _fixture(root / "context-case")
            run_path = review / "metadata/run_manifest.json"
            run_manifest = json.loads(run_path.read_text(encoding="utf-8"))
            run_manifest["method_profile"] = {
                "profile_id": "FORGED-PAPER-METHOD",
                "acquisition": {"dia": {"nominal_isolation_width_mz": 999}},
            }
            run_manifest["config_hash"] = review_config_hash(
                standard_id=run_manifest["standard_id"],
                parameters=run_manifest["parameters"],
                thresholds=run_manifest["thresholds"],
                method_profile=run_manifest["method_profile"],
                model_context=run_manifest["model_artifact"],
            )
            run_path.write_text(json.dumps(run_manifest), encoding="utf-8")
            manifest_path = review / "metadata/target_manifest.csv"
            manifest_fields, manifest_rows = _read_csv(manifest_path)
            manifest_rows[0]["config_hash"] = run_manifest["config_hash"]
            _write_csv(manifest_path, manifest_fields, manifest_rows)
            with self.assertRaisesRegex(ValueError, "sidecar-derived provenance"):
                package_shadow_review(result, review, root / "forged-context.zip")

            result, review = _fixture(root / "manifest-case")
            manifest_path = review / "metadata/target_manifest.csv"
            manifest_fields, manifest_rows = _read_csv(manifest_path)
            manifest_rows[0]["source_machine_label"] = "RED"
            manifest_rows[0]["source_machine_interpretation"] = "unexpected_single_peak"
            _write_csv(manifest_path, manifest_fields, manifest_rows)
            with self.assertRaisesRegex(ValueError, "source_machine_label.*review queue"):
                package_shadow_review(result, review, root / "forged-manifest.zip")

    def test_review_run_manifest_is_bound_to_queue_and_sidecar(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, review = _fixture(root)
            run_path = review / "metadata/run_manifest.json"
            run_manifest = json.loads(run_path.read_text(encoding="utf-8"))
            run_manifest["source_sidecar_sha256"] = "0" * 64
            run_path.write_text(json.dumps(run_manifest), encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "final result sidecar"):
                package_shadow_review(result, review, root / "bad-run-manifest.zip")

    def test_zero_disagreement_package_is_empty_but_valid(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, review = _fixture(root)
            result_fields, result_rows = _read_csv(result)
            result_rows[0]["ms2_diagnostic_status"] = "not_evaluable"
            _write_csv(result, result_fields, result_rows)
            manifest_path = review / "metadata/target_manifest.csv"
            manifest_fields, _manifest_rows = _read_csv(manifest_path)
            _write_csv(manifest_path, manifest_fields, [])
            label_path = review / "annotations/_template/review_labels.csv"
            label_fields, _label_rows = _read_csv(label_path)
            _write_csv(label_path, label_fields, [])
            _sync_review_provenance(result, review)

            output = root / "zero.zip"
            summary = package_shadow_review(result, review, output)

            self.assertEqual(summary["review_rows"], 0)
            self.assertEqual(summary["png_images"], 0)
            self.assertEqual(summary["svg_images"], 0)
            with zipfile.ZipFile(output) as archive:
                self.assertIsNone(archive.testzip())
                self.assertFalse(
                    any(name.endswith((".png", ".svg")) for name in archive.namelist())
                )
                labels = archive.read(f"{PACKAGE_ROOT}/review_labels.csv").decode("utf-8-sig")
                self.assertEqual(len(list(csv.DictReader(StringIO(labels)))), 0)


if __name__ == "__main__":
    unittest.main()
