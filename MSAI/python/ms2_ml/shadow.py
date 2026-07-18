"""Append-only shadow inference for the interpretable MS2 classifier."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import Counter
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from .contracts import (
    MIN_SPECIFICITY_TARGET,
    ML_DIAGNOSTIC_STATUSES,
    ML_OUTPUT_FIELDS,
    MODEL_SCHEMA_VERSION,
    TRAINABLE_DECISION_KEYS,
    canonical_model_id,
    validate_frozen_source,
    validate_trainable_artifact,
)
from .features import FEATURE_NAMES, extract_features_from_row, feature_vector
from .logistic import C_GRID, CLASS_WEIGHT_GRID, L2LogisticModel
from .provenance import canonical_json_sha256, canonical_table_sha256
from .thresholds import (
    COSINE_GRID,
    ENTROPY_SIMILARITY_GRID,
    EXPLAINED_INTENSITY_GRID,
    MATCHED_PEAKS_GRID,
    ThresholdParameters,
    threshold_bucket,
)

_TRAINABLE_SECTION = "trainable_decision"
_SELECTED_MODELS = frozenset({"threshold", "logistic_regression"})
_THRESHOLD_BUCKETS = frozenset(
    {
        "supported_same_compound",
        "conflicting_spectra",
        "insufficient_evidence",
    }
)
_THRESHOLD_GRID_SIZE = (
    len(COSINE_GRID)
    * len(MATCHED_PEAKS_GRID)
    * len(EXPLAINED_INTENSITY_GRID)
    * len(ENTROPY_SIMILARITY_GRID)
)
_FOLD_KEYS = frozenset(
    {
        "name",
        "train_row_ids",
        "validation_row_ids",
        "purged_row_ids",
        "train_batches",
        "validation_batches",
        "train_compounds",
        "validation_compounds",
    }
)
_EVALUATION_METRIC_KEYS = frozenset(
    {
        "threshold",
        "evaluable_rows",
        "tp",
        "fp",
        "tn",
        "fn",
        "precision",
        "recall",
        "specificity",
        "pr_auc",
        "average_precision",
        "f0_5",
        "brier_score",
        "expected_calibration_error",
    }
)


def load_model_artifact(path) -> dict[str, Any]:
    """Load and validate the sole allowed training-artifact section."""

    model_path = Path(path)
    with model_path.open(encoding="utf-8") as handle:
        artifact = json.load(handle)
    return validate_model_artifact(artifact)


def validate_model_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    """Deeply validate an in-memory artifact using the production loader contract."""

    if not isinstance(artifact, Mapping):
        raise ValueError("Model artifact must be a JSON object.")
    payload = validate_trainable_artifact(artifact)
    _validate_model_payload(payload)
    return {_TRAINABLE_SECTION: dict(payload)}


def predict_shadow(
    row: Mapping[str, Any],
    model_artifact: Mapping[str, Any] | None = None,
) -> dict[str, str]:
    """Return exactly five shadow fields without consulting any v2 rule label.

    When independent truth has not produced a model artifact, the explicit
    ``UNTRAINED`` sentinel exercises the output plumbing but emits neither a
    probability nor a same/different-compound conclusion.
    """

    missing = _missing_spectrum_flags(row)
    if model_artifact is None:
        flags = ["SHADOW_ONLY", "FEASIBILITY_ONLY", "NO_INDEPENDENT_TRUTH", *missing]
        return _output(
            probability="",
            status="not_evaluable",
            model_id="UNTRAINED",
            abstention_reason="NO_INDEPENDENT_TRUTH_MODEL",
            flags=flags,
        )

    payload = validate_trainable_artifact(model_artifact)
    _validate_model_payload(payload)
    model_id = str(payload["model_id"])
    if missing:
        reason = "MISSING_BILATERAL_SPECTRA" if len(missing) == 2 else missing[0]
        return _output(
            probability="",
            status="not_evaluable",
            model_id=model_id,
            abstention_reason=reason,
            flags=["SHADOW_ONLY", *missing],
        )

    features = extract_features_from_row(row)
    if features is None:
        return _output(
            probability="",
            status="not_evaluable",
            model_id=model_id,
            abstention_reason="INVALID_OR_NON_EVALUABLE_SPECTRA",
            flags=["SHADOW_ONLY", "INVALID_SPECTRUM_FEATURES"],
        )
    vector = feature_vector(features)
    selected = str(payload["selected_model"])
    model = payload["model"]
    flags = ["SHADOW_ONLY"]
    if _outside_training_range(vector, payload.get("training_feature_ranges")):
        flags.append("FEATURE_OUTSIDE_TRAINING_RANGE")

    if selected == "threshold":
        parameters = ThresholdParameters.from_dict(model["parameters"])
        status = threshold_bucket(vector, parameters)
        bucket_probabilities = model["bucket_probabilities"]
        probability = _probability(bucket_probabilities.get(status), f"bucket {status!r}")
        flags.append("INTERPRETABLE_THRESHOLD_MODEL")
        reason = "TOO_FEW_MATCHED_FRAGMENTS" if status == "insufficient_evidence" else ""
    else:
        probability = _logistic_probability(vector, model)
        decision_thresholds = payload["decision_thresholds"]
        conflict = _probability(
            decision_thresholds.get("conflict_max_probability"),
            "conflict_max_probability",
        )
        support = _probability(
            decision_thresholds.get("support_min_probability"),
            "support_min_probability",
        )
        if conflict >= support:
            raise ValueError("Conflict probability threshold must be below support threshold.")
        if probability <= conflict:
            status, reason = "conflicting_spectra", ""
        elif probability >= support:
            status, reason = "supported_same_compound", ""
        else:
            status, reason = "insufficient_evidence", "PROBABILITY_ABSTENTION_INTERVAL"
        flags.append("LOGISTIC_REGRESSION_MODEL")
        if min(abs(probability - conflict), abs(probability - support)) <= 0.05:
            flags.append("NEAR_ML_DECISION_BOUNDARY")

    return _output(
        probability=f"{probability:.12g}",
        status=status,
        model_id=model_id,
        abstention_reason=reason,
        flags=flags,
    )


def apply_shadow_rows(
    rows: Sequence[Mapping[str, Any]],
    model_artifact: Mapping[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Copy input rows and append exactly the five ML output fields."""

    if model_artifact is not None:
        payload = validate_trainable_artifact(model_artifact)
        _validate_model_payload(payload)
    output: list[dict[str, Any]] = []
    for source in rows:
        existing = [field for field in ML_OUTPUT_FIELDS if field in source]
        if existing:
            raise ValueError(
                "Shadow input already contains ML output fields; refusing to overwrite: "
                + ", ".join(existing)
            )
        row = dict(source)
        prediction = predict_shadow(source, model_artifact)
        if tuple(prediction) != ML_OUTPUT_FIELDS:
            raise AssertionError("Shadow predictor did not return the fixed five-field schema.")
        row.update(prediction)
        output.append(row)
    return output


def apply_shadow_csv(
    input_path,
    output_path,
    *,
    model_artifact_path=None,
    source_sidecar_path=None,
    standard_path=None,
) -> dict[str, Any]:
    """Apply shadow inference after validating the source sidecar and v2 standard."""

    source_path = Path(input_path).resolve()
    destination_path = Path(output_path).resolve()
    if source_path == destination_path:
        raise ValueError("Shadow output must not overwrite its v2 source CSV.")
    fields, rows = _read_csv(source_path)
    existing = [field for field in ML_OUTPUT_FIELDS if field in fields]
    if existing:
        raise ValueError(
            "Shadow input already contains ML output fields; refusing to overwrite: "
            + ", ".join(existing)
        )
    source_sidecar_path = (
        Path(source_sidecar_path).resolve()
        if source_sidecar_path
        else Path(str(source_path) + ".metadata.json")
    )
    if not source_sidecar_path.is_file():
        raise ValueError(f"Source metadata sidecar does not exist: {source_sidecar_path}")
    source_sidecar = _read_json_mapping(source_sidecar_path, "Source metadata sidecar")
    standard_path = (
        Path(standard_path).resolve()
        if standard_path
        else Path(__file__).resolve().parents[2] / "standards" / "ms2_diagnostic_standard_v2.json"
    )
    if not standard_path.is_file():
        raise ValueError(f"Diagnostic standard does not exist: {standard_path}")
    standard = _read_json_mapping(standard_path, "Diagnostic standard")
    validate_frozen_source(source_sidecar, standard)

    artifact = load_model_artifact(model_artifact_path) if model_artifact_path else None
    output_rows = apply_shadow_rows(rows, artifact)
    output_fields = [*fields, *ML_OUTPUT_FIELDS]
    _write_csv(destination_path, output_fields, output_rows)

    model_path = Path(model_artifact_path).resolve() if model_artifact_path else None
    model_id = (
        str(artifact[_TRAINABLE_SECTION]["model_id"]) if artifact is not None else "UNTRAINED"
    )
    summary = {
        "schema_version": "msai-ms2-shadow-output-v1",
        "created_utc": datetime.now(UTC).isoformat(),
        "source_result": str(source_path),
        "source_result_sha256": _sha256_file(source_path),
        "source_result_content_sha256": canonical_table_sha256(fields, rows),
        "source_sidecar": str(source_sidecar_path),
        "source_sidecar_sha256": _sha256_file(source_sidecar_path),
        "source_sidecar_content_sha256": canonical_json_sha256(source_sidecar),
        "diagnostic_standard": str(standard_path),
        "diagnostic_standard_sha256": _sha256_file(standard_path),
        "output_result": str(destination_path),
        "output_result_sha256": _sha256_file(destination_path),
        "append_only_fields": list(ML_OUTPUT_FIELDS),
        "rows": len(output_rows),
        "status_counts": dict(Counter(str(row["ml_diagnostic_status"]) for row in output_rows)),
        "model": {
            "model_id": model_id,
            "path": str(model_path) if model_path else None,
            "sha256": _sha256_file(model_path) if model_path else None,
            "trained": artifact is not None,
        },
    }
    sidecar = Path(str(destination_path) + ".metadata.json")
    combined_metadata = dict(source_sidecar)
    combined_metadata["shadow_inference"] = summary
    sidecar.write_text(
        json.dumps(combined_metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return summary


def _validate_model_payload(payload: Mapping[str, Any]) -> None:
    if set(payload) != TRAINABLE_DECISION_KEYS:
        missing = sorted(TRAINABLE_DECISION_KEYS - set(payload))
        unexpected = sorted(set(payload) - TRAINABLE_DECISION_KEYS)
        details = []
        if missing:
            details.append("missing=" + ",".join(missing))
        if unexpected:
            details.append("unexpected=" + ",".join(unexpected))
        raise ValueError("Model artifact trainable_decision schema mismatch: " + "; ".join(details))
    if payload.get("schema_version") != MODEL_SCHEMA_VERSION:
        raise ValueError(f"Model artifact schema_version must be {MODEL_SCHEMA_VERSION!r}.")
    selected = str(payload.get("selected_model"))
    if selected not in _SELECTED_MODELS:
        raise ValueError(f"Unsupported selected_model: {selected!r}.")
    feature_names = tuple(str(name) for name in payload.get("feature_names") or ())
    if feature_names != FEATURE_NAMES:
        raise ValueError("Model artifact feature_names do not match the frozen feature order.")
    model = payload.get("model")
    if not isinstance(model, Mapping) or model.get("type") != selected:
        raise ValueError("Serialized model type must match selected_model.")
    thresholds = payload.get("decision_thresholds")
    if not isinstance(thresholds, Mapping) or set(thresholds) != {
        "conflict_max_probability",
        "support_min_probability",
    }:
        raise ValueError("decision_thresholds must use the exact two-field schema.")
    if selected == "threshold":
        if set(model) != {"type", "parameters", "bucket_probabilities"}:
            raise ValueError("Threshold model must use the exact serialized schema.")
        if not isinstance(model.get("parameters"), Mapping):
            raise ValueError("Threshold model is missing parameters.")
        ThresholdParameters.from_dict(model["parameters"])
        probabilities = model.get("bucket_probabilities")
        if not isinstance(probabilities, Mapping) or set(probabilities) != _THRESHOLD_BUCKETS:
            raise ValueError("Threshold bucket_probabilities must use the exact three buckets.")
        for status in sorted(_THRESHOLD_BUCKETS):
            _probability(probabilities.get(status), f"bucket {status!r}")
        conflict = _probability(
            thresholds.get("conflict_max_probability"),
            "conflict_max_probability",
        )
        support = _probability(
            thresholds.get("support_min_probability"),
            "support_min_probability",
        )
        if conflict != 0.0 or support != 1.0:
            raise ValueError("Threshold models require decision-threshold sentinels 0.0 and 1.0.")
    else:
        _validate_logistic_model_schema(model)
        conflict = _probability(
            thresholds.get("conflict_max_probability"),
            "conflict_max_probability",
        )
        support = _probability(
            thresholds.get("support_min_probability"),
            "support_min_probability",
        )
        if conflict >= support:
            raise ValueError("Conflict probability threshold must be below support threshold.")
        _logistic_probability((0.0,) * len(FEATURE_NAMES), model)

    _validate_training_feature_ranges(payload.get("training_feature_ranges"))
    _validate_training_evidence(payload, selected)
    observed_model_id = payload.get("model_id")
    expected_model_id = canonical_model_id(payload)
    if observed_model_id != expected_model_id:
        raise ValueError(
            "Model artifact model_id does not match the canonical trainable payload hash."
        )


def _validate_logistic_model_schema(model: Mapping[str, Any]) -> None:
    expected = {
        "type",
        "model_type",
        "feature_names",
        "scaler",
        "coefficients",
        "intercept",
        "C",
        "class_weight",
        "optimizer",
    }
    if set(model) != expected:
        raise ValueError("Logistic model must use the exact serialized schema.")
    scaler = model.get("scaler")
    if not isinstance(scaler, Mapping) or set(scaler) != {
        "feature_names",
        "means",
        "scales",
    }:
        raise ValueError("Logistic scaler must use the exact serialized schema.")
    optimizer = model.get("optimizer")
    if not isinstance(optimizer, Mapping) or set(optimizer) != {
        "algorithm",
        "iterations",
        "converged",
        "objective",
    }:
        raise ValueError("Logistic optimizer must use the exact serialized schema.")
    if optimizer.get("algorithm") != "deterministic_newton_irls":
        raise ValueError("Unsupported logistic optimizer algorithm.")
    if optimizer.get("converged") is not True:
        raise ValueError("Logistic optimizer must record converged=true for deployment.")
    iterations = optimizer.get("iterations")
    if (
        isinstance(iterations, bool)
        or not isinstance(iterations, int)
        or not 1 <= iterations <= 200
    ):
        raise ValueError("Logistic optimizer iterations must be an integer in [1, 200].")
    _finite_number(optimizer.get("objective"), "logistic optimizer objective")


def _validate_training_feature_ranges(value: Any) -> None:
    if not isinstance(value, Mapping) or set(value) != set(FEATURE_NAMES):
        raise ValueError("training_feature_ranges must cover the six fixed features exactly.")
    for name in FEATURE_NAMES:
        bounds = value.get(name)
        if not isinstance(bounds, Mapping) or set(bounds) != {"minimum", "maximum"}:
            raise ValueError(f"Training feature range for {name!r} must contain minimum/maximum.")
        minimum = _finite_number(bounds.get("minimum"), f"{name} minimum")
        maximum = _finite_number(bounds.get("maximum"), f"{name} maximum")
        if minimum > maximum:
            raise ValueError(f"Training feature range for {name!r} is reversed.")


def _validate_training_evidence(payload: Mapping[str, Any], selected: str) -> None:
    mapping_sections = (
        "candidates",
        "selection",
        "development_validation",
        "locked_final_evaluation",
        "split_manifest",
        "annotation_summary",
        "provenance",
        "model_card",
    )
    for name in mapping_sections:
        value = payload.get(name)
        if not isinstance(value, Mapping) or not value:
            raise ValueError(f"Model artifact requires non-empty training evidence {name!r}.")

    selection = payload["selection"]
    development = payload["development_validation"]
    specificity_target = _validate_specificity_evidence(selection, development)
    split_info = _validate_split_evidence(payload["split_manifest"])
    outer_metrics = _validate_development_evidence(
        development,
        split_info["outer_folds"],
        specificity_target,
    )
    _validate_candidate_evidence(payload["candidates"], specificity_target)
    _validate_selected_model_matches_final_tuning(payload, selected)
    _validate_selection_evidence(
        selection,
        selected,
        outer_metrics,
        specificity_target,
    )

    locked = payload["locked_final_evaluation"]
    if set(locked) != {"batch", "used_once_after_model_selection", "evaluation"}:
        raise ValueError("Locked final evidence must use the exact three-field schema.")
    if locked.get("used_once_after_model_selection") is not True:
        raise ValueError("Locked final evidence must confirm one-time post-selection use.")
    if not str(locked.get("batch") or "").strip() or not isinstance(
        locked.get("evaluation"), Mapping
    ):
        raise ValueError("Locked final evidence requires a batch and evaluation mapping.")
    locked_batch = next(iter(split_info["locked_batches"]))
    if locked.get("batch") != locked_batch:
        raise ValueError("Locked final evaluation batch does not match split evidence.")
    _validate_evaluation_record(
        locked["evaluation"],
        "locked final",
        expected_rows=len(split_info["locked_row_ids"]),
        expected_compounds=len(split_info["locked_compounds"]),
    )

    provenance = payload["provenance"]
    if provenance.get("truth_contract") != "independent_external_labels_required":
        raise ValueError("Model provenance must require independent external truth labels.")
    if provenance.get("frozen_source_contract_validated") is not True:
        raise ValueError("Model provenance must confirm frozen source-contract validation.")
    _validate_provenance_evidence(
        provenance,
        eligible_batches=split_info["eligible_batches"],
        locked_batch=locked_batch,
    )

    annotation = payload["annotation_summary"]
    expected_annotation_keys = {
        "input_rows",
        "positive_same_compound",
        "negative_different_or_interference",
        "uncertain",
        "not_evaluable",
        "eligible_rows",
    }
    if set(annotation) != expected_annotation_keys:
        raise ValueError("Annotation evidence must use the exact truth-count schema.")
    for name in expected_annotation_keys:
        count = annotation.get(name)
        if isinstance(count, bool) or not isinstance(count, int) or count < 0:
            raise ValueError(f"Annotation evidence {name!r} must be a non-negative integer.")
    for label in (
        "positive_same_compound",
        "negative_different_or_interference",
    ):
        count = annotation.get(label)
        if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
            raise ValueError(f"Annotation evidence requires a positive {label!r} count.")
    if annotation["eligible_rows"] != (
        annotation["positive_same_compound"] + annotation["negative_different_or_interference"]
    ):
        raise ValueError("eligible_rows must equal the two fit-eligible truth counts.")
    if annotation["input_rows"] != sum(
        annotation[name]
        for name in (
            "positive_same_compound",
            "negative_different_or_interference",
            "uncertain",
            "not_evaluable",
        )
    ):
        raise ValueError("input_rows must equal all four explicit truth-label counts.")
    if annotation["eligible_rows"] != len(split_info["eligible_row_ids"]):
        raise ValueError("Annotation eligible_rows does not match the split row universe.")

    limitations = payload.get("limitations")
    if (
        not isinstance(limitations, Sequence)
        or isinstance(limitations, str | bytes)
        or not limitations
    ):
        raise ValueError("Model artifact requires a non-empty limitations list.")


def _validate_specificity_evidence(
    selection: Mapping[str, Any],
    development: Mapping[str, Any],
) -> float:
    selection_target = _specificity_target(selection.get("specificity_target"), "selection")
    development_target = _specificity_target(
        development.get("specificity_target"),
        "development_validation",
    )
    if selection_target != development_target:
        raise ValueError("Selection and development specificity_target values must match exactly.")
    return selection_target


def _specificity_target(value: Any, label: str) -> float:
    target = _finite_number(value, f"{label} specificity_target")
    if not MIN_SPECIFICITY_TARGET <= target <= 1.0:
        raise ValueError(f"{label} specificity_target must be in [{MIN_SPECIFICITY_TARGET}, 1].")
    return target


def _validate_candidate_evidence(candidates: Mapping[str, Any], target: float) -> None:
    if set(candidates) != {"threshold", "logistic_regression"}:
        raise ValueError("Training candidates must document both permitted model classes.")
    threshold = candidates.get("threshold")
    logistic = candidates.get("logistic_regression")
    if not isinstance(threshold, Mapping) or set(threshold) != {
        "grid_size",
        "final_development_tuning",
    }:
        raise ValueError("Threshold candidate evidence must use the exact training schema.")
    if threshold.get("grid_size") != _THRESHOLD_GRID_SIZE:
        raise ValueError("Threshold candidate evidence must record the fixed 6,384-point grid.")
    _validate_threshold_tuning(threshold.get("final_development_tuning"), target)

    if not isinstance(logistic, Mapping) or set(logistic) != {
        "c_grid",
        "class_weight_grid",
        "final_development_tuning",
    }:
        raise ValueError("Logistic candidate evidence must use the exact training schema.")
    if logistic.get("c_grid") != list(C_GRID) or logistic.get("class_weight_grid") != list(
        CLASS_WEIGHT_GRID
    ):
        raise ValueError("Logistic candidate evidence must record the fixed C/weight grids.")
    _validate_logistic_tuning(logistic.get("final_development_tuning"), target)


def _validate_selected_model_matches_final_tuning(
    payload: Mapping[str, Any], selected: str
) -> None:
    candidates = payload["candidates"]
    model = payload["model"]
    decision_thresholds = payload["decision_thresholds"]
    if selected == "threshold":
        tuning = candidates["threshold"]["final_development_tuning"]
        if dict(model["parameters"]) != dict(tuning["parameters"]):
            raise ValueError("Deployed threshold parameters do not match final tuning evidence.")
        if dict(model["bucket_probabilities"]) != dict(tuning["bucket_probabilities"]):
            raise ValueError(
                "Deployed threshold bucket probabilities do not match final tuning evidence."
            )
        return

    tuning = candidates["logistic_regression"]["final_development_tuning"]
    if _finite_number(model.get("C"), "deployed Logistic C") != _finite_number(
        tuning.get("c"), "final tuning Logistic C"
    ) or model.get("class_weight") != tuning.get("class_weight"):
        raise ValueError("Deployed Logistic configuration does not match final tuning evidence.")
    if _probability(
        decision_thresholds.get("support_min_probability"),
        "deployed support threshold",
    ) != _probability(
        tuning.get("support_threshold"), "final tuning support threshold"
    ) or _probability(
        decision_thresholds.get("conflict_max_probability"),
        "deployed conflict threshold",
    ) != _probability(tuning.get("conflict_threshold"), "final tuning conflict threshold"):
        raise ValueError(
            "Deployed Logistic probability boundaries do not match final tuning evidence."
        )


def _validate_threshold_tuning(value: Any, target: float) -> None:
    expected = {
        "parameters",
        "metrics",
        "candidates_evaluated",
        "specificity_target",
        "constraint_met",
        "bucket_probabilities",
    }
    if not isinstance(value, Mapping) or set(value) != expected:
        raise ValueError("Threshold tuning evidence must use the exact serialized schema.")
    parameters = value.get("parameters")
    if not isinstance(parameters, Mapping):
        raise ValueError("Threshold tuning evidence requires parameters.")
    ThresholdParameters.from_dict(parameters)
    if value.get("candidates_evaluated") != _THRESHOLD_GRID_SIZE:
        raise ValueError("Threshold tuning must evaluate all 6,384 fixed candidates.")
    if _specificity_target(value.get("specificity_target"), "threshold tuning") != target:
        raise ValueError("Threshold tuning specificity_target does not match development.")
    if value.get("constraint_met") is not True:
        raise ValueError("Threshold tuning must satisfy the specificity constraint.")
    metrics = value.get("metrics")
    validated_metrics = _validate_metric_record(metrics, "threshold tuning")
    specificity = float(validated_metrics["specificity"])
    if specificity < target:
        raise ValueError("Threshold tuning metrics violate the specificity target.")
    probabilities = value.get("bucket_probabilities")
    if not isinstance(probabilities, Mapping) or set(probabilities) != _THRESHOLD_BUCKETS:
        raise ValueError("Threshold tuning bucket probabilities use an invalid schema.")
    for status in _THRESHOLD_BUCKETS:
        _probability(probabilities.get(status), f"threshold tuning bucket {status!r}")


def _validate_logistic_tuning(value: Any, target: float) -> None:
    expected = {
        "c",
        "class_weight",
        "support_threshold",
        "conflict_threshold",
        "metrics",
        "candidates_evaluated",
        "candidates_attempted",
        "convergence_failures",
    }
    if not isinstance(value, Mapping) or set(value) != expected:
        raise ValueError("Logistic tuning evidence must use the exact serialized schema.")
    c_value = _finite_number(value.get("c"), "Logistic tuning C")
    if c_value not in C_GRID or value.get("class_weight") not in CLASS_WEIGHT_GRID:
        raise ValueError("Logistic tuning must select from the fixed C/weight grids.")
    attempted = value.get("candidates_attempted")
    evaluated = value.get("candidates_evaluated")
    expected_attempted = len(C_GRID) * len(CLASS_WEIGHT_GRID)
    if attempted != expected_attempted:
        raise ValueError("Logistic tuning must attempt the fixed 14 candidates.")
    if (
        isinstance(evaluated, bool)
        or not isinstance(evaluated, int)
        or not 1 <= evaluated <= expected_attempted
    ):
        raise ValueError("Logistic tuning must evaluate at least one converged candidate.")
    failures = value.get("convergence_failures")
    if (
        not isinstance(failures, Sequence)
        or isinstance(failures, str | bytes)
        or any(not isinstance(item, str) or not item.strip() for item in failures)
        or len(failures) != expected_attempted - evaluated
    ):
        raise ValueError("Logistic tuning convergence-failure evidence is inconsistent.")
    conflict = _probability(value.get("conflict_threshold"), "Logistic conflict threshold")
    support = _probability(value.get("support_threshold"), "Logistic support threshold")
    if conflict >= support:
        raise ValueError("Logistic tuning conflict threshold must be below support threshold.")
    metrics = value.get("metrics")
    validated_metrics = _validate_metric_record(metrics, "Logistic tuning")
    specificity = float(validated_metrics["specificity"])
    if specificity < target:
        raise ValueError("Logistic tuning metrics violate the specificity target.")


def _validate_evaluation_record(
    value: Any,
    label: str,
    *,
    expected_rows: int,
    expected_compounds: int,
) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} evaluation must be a mapping.")
    if not {"metrics", "row_count", "compound_count"} <= set(value):
        raise ValueError(f"{label} evaluation is missing metrics or group counts.")
    row_count = value.get("row_count")
    compound_count = value.get("compound_count")
    if isinstance(row_count, bool) or not isinstance(row_count, int) or row_count != expected_rows:
        raise ValueError(f"{label} evaluation row_count does not match split evidence.")
    if (
        isinstance(compound_count, bool)
        or not isinstance(compound_count, int)
        or compound_count != expected_compounds
        or not 1 <= compound_count <= row_count
    ):
        raise ValueError(f"{label} evaluation compound_count does not match split evidence.")
    metrics = _validate_metric_record(
        value.get("metrics"),
        label,
        expected_rows=row_count,
    )
    return metrics


def _validate_metric_record(
    value: Any,
    label: str,
    *,
    expected_rows: int | None = None,
) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or set(value) != _EVALUATION_METRIC_KEYS:
        raise ValueError(f"{label} metrics must use the exact all-row metric schema.")
    count_names = ("evaluable_rows", "tp", "fp", "tn", "fn")
    counts: dict[str, int] = {}
    for name in count_names:
        count = value.get(name)
        if isinstance(count, bool) or not isinstance(count, int) or count < 0:
            raise ValueError(f"{label} metric {name!r} must be a non-negative integer.")
        counts[name] = count
    evaluable_rows = counts["evaluable_rows"]
    if evaluable_rows < 1:
        raise ValueError(f"{label} metrics require at least one evaluable row.")
    if expected_rows is not None and evaluable_rows != expected_rows:
        raise ValueError(f"{label} confusion counts do not match evaluation row_count.")
    if sum(counts[name] for name in ("tp", "fp", "tn", "fn")) != evaluable_rows:
        raise ValueError(f"{label} confusion counts do not match evaluable_rows.")
    numeric = {
        name: _probability(value.get(name), f"{label} metric {name}")
        for name in _EVALUATION_METRIC_KEYS - set(count_names)
    }
    expected_operating = {
        "precision": _ratio(counts["tp"], counts["tp"] + counts["fp"]),
        "recall": _ratio(counts["tp"], counts["tp"] + counts["fn"]),
        "specificity": _ratio(counts["tn"], counts["tn"] + counts["fp"]),
    }
    expected_operating["f0_5"] = _ratio(
        1.25 * expected_operating["precision"] * expected_operating["recall"],
        0.25 * expected_operating["precision"] + expected_operating["recall"],
    )
    for name, expected in expected_operating.items():
        if not math.isclose(numeric[name], expected, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"{label} metric {name!r} is inconsistent with confusion counts.")
    return value


def _ratio(numerator: float, denominator: float) -> float:
    return 0.0 if denominator == 0 else numerator / denominator


def _validate_selection_evidence(
    selection: Mapping[str, Any],
    selected: str,
    outer_metrics: Mapping[str, Mapping[str, Mapping[str, Any]]],
    target: float,
) -> None:
    expected = {
        "selected_model",
        "threshold_wins_ties",
        "criterion",
        "reliably_better",
        "specificity_target",
        "batch_deltas",
    }
    if set(selection) != expected:
        raise ValueError("Selection evidence must use the exact serialized schema.")
    if selection.get("selected_model") != selected:
        raise ValueError("Selection evidence does not match selected_model.")
    if selection.get("threshold_wins_ties") is not True:
        raise ValueError("Selection evidence must confirm that threshold wins ties.")
    reliable = selection.get("reliably_better")
    if not isinstance(reliable, bool) or reliable != (selected == "logistic_regression"):
        raise ValueError("Selection reliability flag does not match selected_model.")
    if _specificity_target(selection.get("specificity_target"), "selection") != target:
        raise ValueError("Selection specificity_target does not match development.")
    if not isinstance(selection.get("criterion"), str) or not selection["criterion"].strip():
        raise ValueError("Selection evidence requires a non-empty criterion.")
    deltas = selection.get("batch_deltas")
    if (
        not isinstance(deltas, Sequence)
        or isinstance(deltas, str | bytes)
        or len(deltas) != len(outer_metrics)
    ):
        raise ValueError("Selection batch_deltas must cover every outer validation batch.")
    observed_batches: set[str] = set()
    logistic_is_reliable = True
    for index, raw_delta in enumerate(deltas):
        expected_delta = {
            "batch",
            "recall_delta",
            "threshold_recall",
            "logistic_recall",
            "logistic_specificity",
        }
        if not isinstance(raw_delta, Mapping) or set(raw_delta) != expected_delta:
            raise ValueError(f"Selection batch delta {index} has an invalid schema.")
        batch = str(raw_delta.get("batch") or "").strip()
        if not batch or batch in observed_batches:
            raise ValueError("Selection batch_deltas require unique non-empty batches.")
        observed_batches.add(batch)
        threshold_recall = _probability(raw_delta.get("threshold_recall"), "threshold batch recall")
        logistic_recall = _probability(raw_delta.get("logistic_recall"), "Logistic batch recall")
        logistic_specificity = _probability(
            raw_delta.get("logistic_specificity"), "Logistic batch specificity"
        )
        if batch not in outer_metrics:
            raise ValueError("Selection batch_deltas reference an unknown validation batch.")
        expected_threshold = outer_metrics[batch]["threshold"]
        expected_logistic = outer_metrics[batch]["logistic_regression"]
        if (
            not math.isclose(
                threshold_recall,
                float(expected_threshold["recall"]),
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            or not math.isclose(
                logistic_recall,
                float(expected_logistic["recall"]),
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            or not math.isclose(
                logistic_specificity,
                float(expected_logistic["specificity"]),
                rel_tol=0.0,
                abs_tol=1e-12,
            )
        ):
            raise ValueError("Selection batch metrics do not match outer-fold evaluations.")
        recall_delta = _finite_number(raw_delta.get("recall_delta"), "batch recall delta")
        if not math.isclose(
            recall_delta,
            logistic_recall - threshold_recall,
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise ValueError("Selection recall_delta does not match recorded recalls.")
        logistic_is_reliable &= logistic_specificity >= target and recall_delta > 0
    if observed_batches != set(outer_metrics):
        raise ValueError("Selection batch_deltas do not match outer validation batches.")
    if reliable != logistic_is_reliable:
        raise ValueError("Selection reliability does not follow the prespecified batch rule.")


def _validate_development_evidence(
    development: Mapping[str, Any],
    split_outer_folds: Mapping[str, Mapping[str, Any]],
    target: float,
) -> dict[str, dict[str, Mapping[str, Any]]]:
    if set(development) != {"strategy", "specificity_target", "outer_folds", "pooled"}:
        raise ValueError("Development evidence must use the exact serialized schema.")
    if development.get("strategy") != "nested_leave_one_batch_out_with_compound_purge":
        raise ValueError("Development evidence must use nested grouped leave-one-batch-out CV.")
    if _specificity_target(development.get("specificity_target"), "development") != target:
        raise ValueError("Development specificity_target is inconsistent.")
    outer_folds = development.get("outer_folds")
    if (
        not isinstance(outer_folds, Sequence)
        or isinstance(outer_folds, str | bytes)
        or len(outer_folds) != len(split_outer_folds)
    ):
        raise ValueError("Development evidence must match every split outer fold.")
    observed_names: set[str] = set()
    validation_batches: set[str] = set()
    outer_metrics: dict[str, dict[str, Mapping[str, Any]]] = {}
    for index, record in enumerate(outer_folds):
        expected = {
            "fold",
            "validation_batch",
            "validation_compounds",
            "threshold",
            "logistic_regression",
        }
        if not isinstance(record, Mapping) or set(record) != expected:
            raise ValueError(f"Development outer fold {index} has an invalid schema.")
        name = str(record.get("fold") or "").strip()
        if not name or name in observed_names or name not in split_outer_folds:
            raise ValueError("Development outer fold names must uniquely match split evidence.")
        observed_names.add(name)
        split_fold = split_outer_folds[name]
        split_validation_batches = split_fold["validation_batches"]
        batch = str(record.get("validation_batch") or "").strip()
        if {batch} != split_validation_batches:
            raise ValueError("Development validation_batch does not match split evidence.")
        validation_batches.add(batch)
        compounds = _string_set(
            record.get("validation_compounds"),
            f"development.outer[{index}].validation_compounds",
        )
        if compounds != split_fold["validation_compounds"]:
            raise ValueError("Development validation compounds do not match split evidence.")
        batch_metrics: dict[str, Mapping[str, Any]] = {}
        for model_name, validator in (
            ("threshold", _validate_threshold_tuning),
            ("logistic_regression", _validate_logistic_tuning),
        ):
            candidate = record.get(model_name)
            if not isinstance(candidate, Mapping) or set(candidate) != {"tuning", "evaluation"}:
                raise ValueError(f"Development {model_name} evidence has an invalid schema.")
            validator(candidate.get("tuning"), target)
            batch_metrics[model_name] = _validate_evaluation_record(
                candidate.get("evaluation"),
                f"development {model_name} fold {name}",
                expected_rows=len(split_fold["validation_row_ids"]),
                expected_compounds=len(split_fold["validation_compounds"]),
            )
        outer_metrics[batch] = batch_metrics
    if observed_names != set(split_outer_folds):
        raise ValueError("Development evidence does not cover all split outer folds.")
    pooled = development.get("pooled")
    if not isinstance(pooled, Mapping) or set(pooled) != _SELECTED_MODELS:
        raise ValueError("Development pooled evidence must cover both model classes.")
    development_rows = set().union(
        *(fold["validation_row_ids"] for fold in split_outer_folds.values())
    )
    development_compounds = set().union(
        *(fold["validation_compounds"] for fold in split_outer_folds.values())
    )
    for model_name in sorted(_SELECTED_MODELS):
        _validate_evaluation_record(
            pooled[model_name],
            f"development pooled {model_name}",
            expected_rows=len(development_rows),
            expected_compounds=len(development_compounds),
        )
    return outer_metrics


def _validate_provenance_evidence(
    provenance: Mapping[str, Any],
    *,
    eligible_batches: set[str],
    locked_batch: str,
) -> None:
    expected = {
        "training_input",
        "source_sidecar",
        "diagnostic_standard",
        "seed",
        "label_column",
        "compound_group_column",
        "batch_group_column",
        "row_id_column",
        "locked_final_batch",
        "truth_contract",
        "frozen_source_contract_validated",
        "batch_source_contract_validated",
        "batch_sources",
    }
    if set(provenance) != expected:
        raise ValueError("Model provenance must use the exact release-evidence schema.")
    seed = provenance.get("seed")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("Model provenance requires a non-negative integer seed.")
    for name in ("training_input", "source_sidecar", "diagnostic_standard"):
        entry = provenance.get(name)
        if not isinstance(entry, Mapping) or set(entry) != {"path", "sha256"}:
            raise ValueError(f"Model provenance requires {name!r} path/hash evidence.")
        path = entry.get("path")
        digest = entry.get("sha256")
        if not isinstance(path, str) or not path.strip():
            raise ValueError(f"Model provenance {name!r} path must be non-empty.")
        if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
            raise ValueError(f"Model provenance {name!r} sha256 must be lowercase hex.")
    for name in ("label_column", "compound_group_column", "batch_group_column", "row_id_column"):
        if not isinstance(provenance.get(name), str) or not provenance[name].strip():
            raise ValueError(f"Model provenance {name!r} must be a non-empty string.")
    if provenance.get("locked_final_batch") != locked_batch:
        raise ValueError("Model provenance locked_final_batch does not match split evidence.")
    if provenance.get("batch_source_contract_validated") is not True:
        raise ValueError("Model provenance must confirm per-batch source-contract validation.")
    batch_sources = provenance.get("batch_sources")
    if not isinstance(batch_sources, Mapping) or set(batch_sources) != eligible_batches:
        raise ValueError("Model provenance batch_sources must exactly cover eligible batches.")
    for batch, raw_entry in batch_sources.items():
        if not isinstance(raw_entry, Mapping) or set(raw_entry) != {
            "path",
            "sha256",
            "proof_type",
        }:
            raise ValueError(f"Model provenance batch_sources[{batch!r}] has an invalid schema.")
        path = raw_entry.get("path")
        digest = raw_entry.get("sha256")
        proof_type = raw_entry.get("proof_type")
        if not isinstance(path, str) or not path.strip():
            raise ValueError(f"Model provenance batch_sources[{batch!r}] path must be non-empty.")
        if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
            raise ValueError(
                f"Model provenance batch_sources[{batch!r}] sha256 must be lowercase hex."
            )
        if proof_type not in {"source_file", "embedded_canonical_json"}:
            raise ValueError(f"Model provenance batch_sources[{batch!r}] proof_type is invalid.")
        if proof_type == "embedded_canonical_json" and path != f"<embedded:{batch}>":
            raise ValueError(
                f"Embedded batch source path for {batch!r} must use the canonical sentinel."
            )


def _validate_split_evidence(split: Mapping[str, Any]) -> dict[str, Any]:
    expected = {
        "strategy",
        "eligible_batches",
        "development_row_ids",
        "locked_final_row_ids",
        "locked_compound_overlap_purged_row_ids",
        "locked_compound_overlap_purged_batches",
        "development_batches",
        "locked_final_batches",
        "development_compounds",
        "locked_final_compounds",
        "outer_folds",
        "final_development_folds",
    }
    if set(split) != expected:
        raise ValueError("Split evidence must use the exact nested grouped-CV schema.")
    if split.get("strategy") != "locked_batch_then_nested_leave_one_batch_out":
        raise ValueError("Split evidence must use the locked grouped nested-CV strategy.")
    eligible_batches = _string_set(split.get("eligible_batches"), "eligible_batches")
    development_batches = _string_set(split.get("development_batches"), "development_batches")
    locked_batches = _string_set(split.get("locked_final_batches"), "locked_final_batches")
    if len(locked_batches) != 1:
        raise ValueError("Split evidence requires exactly one locked final batch.")
    if len(development_batches) < 3:
        raise ValueError("Nested grouped CV requires at least three development batches.")
    locked_purged_batches = _string_set(
        split.get("locked_compound_overlap_purged_batches"),
        "locked_compound_overlap_purged_batches",
        allow_empty=True,
    )
    if development_batches | locked_batches | locked_purged_batches != eligible_batches:
        raise ValueError(
            "eligible_batches must equal development, locked, and locked-purge batches."
        )
    development_compounds = _string_set(split.get("development_compounds"), "development_compounds")
    locked_compounds = _string_set(split.get("locked_final_compounds"), "locked_final_compounds")
    _require_disjoint(development_batches, locked_batches, "development/locked batches")
    _require_disjoint(
        development_compounds,
        locked_compounds,
        "development/locked compounds",
    )
    development_rows = _string_set(split.get("development_row_ids"), "development_row_ids")
    locked_rows = _string_set(split.get("locked_final_row_ids"), "locked_final_row_ids")
    locked_purged_rows = _string_set(
        split.get("locked_compound_overlap_purged_row_ids"),
        "locked_compound_overlap_purged_row_ids",
        allow_empty=True,
    )
    _require_pairwise_disjoint(
        (development_rows, locked_rows, locked_purged_rows),
        "development/locked/purged row IDs",
    )
    eligible_row_ids = development_rows | locked_rows | locked_purged_rows
    outer_folds = split.get("outer_folds")
    if (
        not isinstance(outer_folds, Sequence)
        or isinstance(outer_folds, str | bytes)
        or len(outer_folds) != len(development_batches)
    ):
        raise ValueError("Split evidence requires one outer fold per development batch.")
    parsed_outer: dict[str, Mapping[str, Any]] = {}
    outer_validation_batches: set[str] = set()
    for index, fold in enumerate(outer_folds):
        parsed = _validate_fold_manifest(
            fold,
            f"outer[{index}]",
            parent_rows=development_rows,
            parent_batches=development_batches,
            parent_compounds=development_compounds,
            with_inner=True,
        )
        name = parsed["name"]
        if name in parsed_outer:
            raise ValueError("Outer fold names must be unique.")
        parsed_outer[name] = parsed
        validation_batch = next(iter(parsed["validation_batches"]))
        if validation_batch in outer_validation_batches:
            raise ValueError("Outer validation batches must be held out exactly once.")
        outer_validation_batches.add(validation_batch)
        inner_folds = fold.get("inner_folds")
        if (
            not isinstance(inner_folds, Sequence)
            or isinstance(inner_folds, str | bytes)
            or len(inner_folds) < 2
        ):
            raise ValueError(f"Outer fold {index} must record at least two inner folds.")
        inner_validation_batches: set[str] = set()
        for inner_index, inner_fold in enumerate(inner_folds):
            inner = _validate_fold_manifest(
                inner_fold,
                f"outer[{index}].inner[{inner_index}]",
                parent_rows=parsed["train_row_ids"],
                parent_batches=parsed["train_batches"],
                parent_compounds=parsed["train_compounds"],
                with_inner=False,
            )
            validation = next(iter(inner["validation_batches"]))
            if validation in inner_validation_batches:
                raise ValueError("Inner validation batches must be held out exactly once.")
            inner_validation_batches.add(validation)
        if inner_validation_batches != parsed["train_batches"]:
            raise ValueError("Inner validation folds must cover the outer training batches.")
    if outer_validation_batches != development_batches:
        raise ValueError("Outer validation folds must cover every development batch.")

    final_folds = split.get("final_development_folds")
    if (
        not isinstance(final_folds, Sequence)
        or isinstance(final_folds, str | bytes)
        or len(final_folds) != len(development_batches)
    ):
        raise ValueError("Final development evidence requires one fold per development batch.")
    final_validation_batches: set[str] = set()
    for index, fold in enumerate(final_folds):
        parsed = _validate_fold_manifest(
            fold,
            f"final_development[{index}]",
            parent_rows=development_rows,
            parent_batches=development_batches,
            parent_compounds=development_compounds,
            with_inner=False,
        )
        validation = next(iter(parsed["validation_batches"]))
        if validation in final_validation_batches:
            raise ValueError("Final development validation batches must be unique.")
        final_validation_batches.add(validation)
    if final_validation_batches != development_batches:
        raise ValueError("Final development folds must cover every development batch.")
    return {
        "eligible_batches": eligible_batches,
        "eligible_row_ids": eligible_row_ids,
        "locked_batches": locked_batches,
        "locked_row_ids": locked_rows,
        "locked_compounds": locked_compounds,
        "outer_folds": parsed_outer,
    }


def _validate_fold_manifest(
    value: Any,
    label: str,
    *,
    parent_rows: set[str],
    parent_batches: set[str],
    parent_compounds: set[str],
    with_inner: bool,
) -> dict[str, Any]:
    expected = set(_FOLD_KEYS)
    if with_inner:
        expected.add("inner_folds")
    if not isinstance(value, Mapping) or set(value) != expected:
        raise ValueError(f"Split evidence {label!r} fold schema is invalid.")
    name = value.get("name")
    if not isinstance(name, str) or not name.strip():
        raise ValueError(f"Split evidence {label!r} requires a non-empty fold name.")
    train_rows = _string_set(value.get("train_row_ids"), f"{label}.train_row_ids")
    validation_rows = _string_set(value.get("validation_row_ids"), f"{label}.validation_row_ids")
    purged_rows = _string_set(
        value.get("purged_row_ids"), f"{label}.purged_row_ids", allow_empty=True
    )
    _require_pairwise_disjoint((train_rows, validation_rows, purged_rows), f"{label} row roles")
    if train_rows | validation_rows | purged_rows != parent_rows:
        raise ValueError(f"Split evidence {label!r} rows do not partition the parent universe.")
    train_batches = _string_set(value.get("train_batches"), f"{label}.train_batches")
    validation_batches = _string_set(value.get("validation_batches"), f"{label}.validation_batches")
    if len(validation_batches) != 1:
        raise ValueError(f"Split evidence {label!r} must hold out exactly one batch.")
    _require_disjoint(train_batches, validation_batches, f"{label} batches")
    if not train_batches | validation_batches <= parent_batches:
        raise ValueError(f"Split evidence {label!r} references a batch outside its parent.")
    train_compounds = _string_set(value.get("train_compounds"), f"{label}.train_compounds")
    validation_compounds = _string_set(
        value.get("validation_compounds"), f"{label}.validation_compounds"
    )
    _require_disjoint(train_compounds, validation_compounds, f"{label} compounds")
    if not train_compounds | validation_compounds <= parent_compounds:
        raise ValueError(f"Split evidence {label!r} references a compound outside its parent.")
    return {
        "name": name,
        "train_row_ids": train_rows,
        "validation_row_ids": validation_rows,
        "purged_row_ids": purged_rows,
        "train_batches": train_batches,
        "validation_batches": validation_batches,
        "train_compounds": train_compounds,
        "validation_compounds": validation_compounds,
    }


def _string_set(value: Any, label: str, *, allow_empty: bool = False) -> set[str]:
    if isinstance(value, Sequence) and not isinstance(value, str | bytes):
        normalized = list(value)
    else:
        normalized = []
    if (
        (not normalized and not allow_empty)
        or any(not isinstance(item, str) or not item.strip() for item in normalized)
        or len(normalized) != len(set(normalized))
    ):
        requirement = "a unique string list" if allow_empty else "a non-empty unique string list"
        raise ValueError(f"Split evidence {label!r} must be {requirement}.")
    return set(normalized)


def _require_pairwise_disjoint(values: Sequence[set[str]], label: str) -> None:
    observed: set[str] = set()
    for value in values:
        overlap = observed & value
        if overlap:
            raise ValueError(f"Split evidence overlaps {label}: {', '.join(sorted(overlap))}.")
        observed.update(value)


def _require_disjoint(first: set[str], second: set[str], label: str) -> None:
    overlap = sorted(first & second)
    if overlap:
        raise ValueError(f"Split evidence overlaps {label}: {', '.join(overlap)}.")


def _logistic_probability(vector: Sequence[float], model: Mapping[str, Any]) -> float:
    try:
        fitted = L2LogisticModel.from_dict(model)
        probability = fitted.predict_proba(vector)
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("Invalid serialized logistic model.") from error
    return _probability(probability, "logistic probability")


def _outside_training_range(vector: Sequence[float], ranges: Any) -> bool:
    """Flag extrapolation without changing the prediction or its status."""

    if not isinstance(ranges, Mapping):
        return False
    for name, value in zip(FEATURE_NAMES, vector, strict=True):
        bounds = ranges.get(name)
        if not isinstance(bounds, Mapping):
            continue
        try:
            minimum = float(bounds["minimum"])
            maximum = float(bounds["maximum"])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(minimum) and math.isfinite(maximum) and not minimum <= value <= maximum:
            return True
    return False


def _missing_spectrum_flags(row: Mapping[str, Any]) -> list[str]:
    flags: list[str] = []
    for field, flag in (
        ("peak_a_MS2", "MISSING_PEAK_A_SPECTRUM"),
        ("peak_b_MS2", "MISSING_PEAK_B_SPECTRUM"),
    ):
        if not str(row.get(field) or "").strip():
            flags.append(flag)
    return flags


def _probability(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"Invalid {label} probability: {value!r}.")
    try:
        probability = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid {label} probability: {value!r}.") from error
    if not math.isfinite(probability) or not 0.0 <= probability <= 1.0:
        raise ValueError(f"Invalid {label} probability: {value!r}.")
    return probability


def _finite_number(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"Invalid finite {label}: {value!r}.")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid finite {label}: {value!r}.") from error
    if not math.isfinite(number):
        raise ValueError(f"Invalid finite {label}: {value!r}.")
    return number


def _output(
    *,
    probability: str,
    status: str,
    model_id: str,
    abstention_reason: str,
    flags: Sequence[str],
) -> dict[str, str]:
    if status not in ML_DIAGNOSTIC_STATUSES:
        raise ValueError(f"Invalid ML diagnostic status: {status!r}.")
    values = (
        probability,
        status,
        model_id,
        abstention_reason,
        ";".join(dict.fromkeys(flag for flag in flags if flag)),
    )
    return dict(zip(ML_OUTPUT_FIELDS, values, strict=True))


def _read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("Shadow source CSV must contain a header.")
        return list(reader.fieldnames), list(reader)


def _read_json_mapping(path: Path, label: str) -> dict[str, Any]:
    try:
        with path.open(encoding="utf-8") as handle:
            value = json.load(handle)
    except json.JSONDecodeError as error:
        raise ValueError(f"{label} must be valid JSON: {path}") from error
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object.")
    return value


def _write_csv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


__all__ = [
    "apply_shadow_csv",
    "apply_shadow_rows",
    "load_model_artifact",
    "predict_shadow",
    "validate_model_artifact",
]
