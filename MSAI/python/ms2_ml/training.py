"""Leakage-resistant training and evaluation for the MS2 shadow classifier.

The pipeline is intentionally narrow: one fixed threshold grid and one
StandardScaler + L2 Logistic Regression candidate.  It refuses row-random
splits, holds one batch untouched for final evaluation, and purges compounds
across every train/validation boundary.
"""

from __future__ import annotations

import hashlib
import json
import math
import random
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .contracts import (
    MIN_SPECIFICITY_TARGET,
    MODEL_SCHEMA_VERSION,
    NEGATIVE_TRUTH_LABEL,
    POSITIVE_TRUTH_LABEL,
    canonical_model_id,
    make_trainable_artifact,
    normalize_truth_label,
    truth_to_binary,
    validate_frozen_source,
    validate_truth_row,
)
from .features import FEATURE_NAMES, extract_features_from_row
from .logistic import (
    C_GRID,
    CLASS_WEIGHT_GRID,
    L2LogisticModel,
    LogisticConvergenceError,
    fit_l2_logistic,
)
from .metrics import binary_metrics, calibration_table, group_bootstrap_ci
from .splits import (
    FoldSplit,
    GroupedExample,
    assert_group_isolation,
    build_leave_one_batch_out_folds,
    build_locked_batch_split,
    validate_binary_fold,
)
from .thresholds import (
    ThresholdSearchResult,
    fit_thresholds,
    select_conflict_probability_threshold,
    select_support_probability_threshold,
    threshold_bucket,
)

DEFAULT_SEED = 20_260_129
DEFAULT_SPECIFICITY_TARGET = MIN_SPECIFICITY_TARGET
DEFAULT_BOOTSTRAP_ITERATIONS = 1_000

_SUPPORT_STATUS = "supported_same_compound"
_CONFLICT_STATUS = "conflicting_spectra"
_ABSTENTION_STATUS = "insufficient_evidence"
_TRAINING_DECISION_STATUSES = (
    _SUPPORT_STATUS,
    _CONFLICT_STATUS,
    _ABSTENTION_STATUS,
)
_PROBABILITY_METRIC_NAMES = (
    "pr_auc",
    "average_precision",
    "brier_score",
    "expected_calibration_error",
)


@dataclass(frozen=True)
class LogisticTuningResult:
    """Nested-CV choice for the only permitted statistical model class."""

    c: float
    class_weight: str | None
    support_threshold: float
    conflict_threshold: float
    metrics: Mapping[str, Any]
    oof_labels: tuple[int, ...]
    oof_probabilities: tuple[float, ...]
    oof_indices: tuple[int, ...]
    candidates_evaluated: int
    candidates_attempted: int
    convergence_failures: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "c": self.c,
            "class_weight": self.class_weight,
            "support_threshold": self.support_threshold,
            "conflict_threshold": self.conflict_threshold,
            "metrics": dict(self.metrics),
            "candidates_evaluated": self.candidates_evaluated,
            "candidates_attempted": self.candidates_attempted,
            "convergence_failures": list(self.convergence_failures),
        }


def train_shadow_model(
    rows: Sequence[Mapping[str, Any]],
    *,
    final_batch: str,
    source_sidecar: Mapping[str, Any],
    diagnostic_standard: Mapping[str, Any],
    label_column: str = "truth_label",
    compound_column: str = "Compound_ID",
    batch_column: str = "batch_id",
    row_id_column: str = "target_uid",
    spectrum_a_column: str = "peak_a_MS2",
    spectrum_b_column: str = "peak_b_MS2",
    input_path: str | Path | None = None,
    source_sidecar_path: str | Path | None = None,
    diagnostic_standard_path: str | Path | None = None,
    provenance_discrepancies: Sequence[str] = (),
    seed: int = DEFAULT_SEED,
    specificity_target: float = DEFAULT_SPECIFICITY_TARGET,
    bootstrap_iterations: int = DEFAULT_BOOTSTRAP_ITERATIONS,
) -> dict[str, dict[str, Any]]:
    """Train both permitted baselines and evaluate once on a locked batch.

    The returned object is a deterministic JSON artifact whose sole top-level
    key is ``trainable_decision``.  No file is written by this function.
    """

    final_batch_name = str(final_batch).strip()
    _validate_training_options(final_batch_name, seed, specificity_target, bootstrap_iterations)
    if (spectrum_a_column, spectrum_b_column) != ("peak_a_MS2", "peak_b_MS2"):
        raise ValueError(
            "Deployment training spectrum columns are frozen as 'peak_a_MS2' and 'peak_b_MS2'."
        )
    # Bind every optional provenance path to the exact in-memory content used
    # below.  Merely hashing a caller-supplied path would let a programmatic
    # caller train on different rows/JSON while the artifact claimed the file.
    training_input_provenance = _provenance_entry(input_path, rows, source_kind="table")
    source_sidecar_provenance = _provenance_entry(
        source_sidecar_path,
        source_sidecar,
        source_kind="json",
    )
    diagnostic_standard_provenance = _provenance_entry(
        diagnostic_standard_path,
        diagnostic_standard,
        source_kind="json",
    )
    examples, annotation_summary = _prepare_examples(
        rows,
        label_column=label_column,
        compound_column=compound_column,
        batch_column=batch_column,
        row_id_column=row_id_column,
        spectrum_a_column=spectrum_a_column,
        spectrum_b_column=spectrum_b_column,
    )
    # Validate the immutable upstream contract for every eligible batch, but
    # store only source hashes in the trainable artifact so Method/extraction
    # values cannot be overwritten.  A single sidecar cannot prove a
    # multi-batch truth set was extracted consistently.
    batch_source_provenance = _validate_batch_source_proof(
        source_sidecar,
        sorted({example.batch_id for example in examples}),
        diagnostic_standard,
        single_source_path=source_sidecar_path,
    )

    locked = build_locked_batch_split(examples, final_batch_name)
    _validate_locked_split(examples, locked.development_indices, locked.final_indices)
    development_indices = tuple(locked.development_indices)
    final_indices = tuple(locked.final_indices)

    outer_folds = tuple(
        build_leave_one_batch_out_folds(
            examples,
            indices=development_indices,
            name_prefix="outer",
        )
    )
    if len(outer_folds) < 2:
        raise ValueError("Nested grouped CV requires at least two evaluable outer folds.")

    outer_records: list[dict[str, Any]] = []
    split_manifest_outer: list[dict[str, Any]] = []
    threshold_outer_statuses: list[str] = []
    threshold_outer_probabilities: list[float] = []
    logistic_outer_scores: list[float] = []
    logistic_outer_statuses: list[str] = []
    outer_labels: list[int] = []
    outer_groups: list[str] = []
    for outer_index, fold in enumerate(outer_folds):
        validate_binary_fold(examples, fold)
        assert_group_isolation(examples, fold.train_indices, fold.validation_indices)
        inner_folds = _nested_folds(
            examples,
            fold.train_indices,
            name_prefix=f"inner_outer_{outer_index}",
        )
        threshold_tuning = _tune_threshold_nested(examples, inner_folds, specificity_target)
        logistic_tuning = _tune_logistic_nested(examples, inner_folds, specificity_target)

        validation_examples = _select(examples, fold.validation_indices)
        validation_labels = [example.label for example in validation_examples]
        validation_groups = [example.compound_id for example in validation_examples]
        threshold_statuses = [
            threshold_bucket(example.features, threshold_tuning.parameters)
            for example in validation_examples
        ]
        threshold_probabilities = [
            threshold_tuning.bucket_probabilities[
                threshold_bucket(example.features, threshold_tuning.parameters)
            ]
            for example in validation_examples
        ]
        logistic_model = _fit_logistic(
            examples,
            fold.train_indices,
            logistic_tuning.c,
            logistic_tuning.class_weight,
        )
        logistic_scores = logistic_model.predict_many(
            [example.features for example in validation_examples]
        )
        logistic_statuses = _logistic_statuses(
            logistic_scores,
            conflict_threshold=logistic_tuning.conflict_threshold,
            support_threshold=logistic_tuning.support_threshold,
        )

        threshold_evaluation = _threshold_evaluation_record(
            validation_labels,
            threshold_statuses,
            threshold_probabilities,
            validation_groups,
            seed=seed + outer_index * 101 + 1,
            bootstrap_iterations=bootstrap_iterations,
        )
        logistic_evaluation = _evaluation_record(
            validation_labels,
            logistic_scores,
            validation_groups,
            statuses=logistic_statuses,
            decision_policy=(
                "probability <= conflict threshold is conflicting; probability >= "
                "support threshold is supported; the middle interval abstains"
            ),
            seed=seed + outer_index * 101 + 2,
            bootstrap_iterations=bootstrap_iterations,
        )
        validation_batch = _single_validation_batch(examples, fold)
        outer_records.append(
            {
                "fold": fold.name,
                "validation_batch": validation_batch,
                "validation_compounds": sorted(set(validation_groups)),
                "threshold": {
                    "tuning": threshold_tuning.to_dict(),
                    "evaluation": threshold_evaluation,
                },
                "logistic_regression": {
                    "tuning": logistic_tuning.to_dict(),
                    "evaluation": logistic_evaluation,
                },
            }
        )
        split_manifest_outer.append(
            {
                **_fold_manifest(examples, fold),
                "inner_folds": [_fold_manifest(examples, inner) for inner in inner_folds],
            }
        )
        outer_labels.extend(validation_labels)
        threshold_outer_statuses.extend(threshold_statuses)
        threshold_outer_probabilities.extend(threshold_probabilities)
        logistic_outer_scores.extend(logistic_scores)
        logistic_outer_statuses.extend(logistic_statuses)
        outer_groups.extend(validation_groups)

    reliability = _logistic_reliability(outer_records, specificity_target)
    selected_model = "logistic_regression" if reliability["reliably_better"] else "threshold"

    final_tuning_folds = _nested_folds(
        examples,
        development_indices,
        name_prefix="final_development",
    )
    final_threshold_tuning = _tune_threshold_nested(
        examples,
        final_tuning_folds,
        specificity_target,
    )
    final_logistic_tuning = _tune_logistic_nested(
        examples,
        final_tuning_folds,
        specificity_target,
    )
    final_examples = _select(examples, final_indices)
    final_labels = [example.label for example in final_examples]
    final_groups = [example.compound_id for example in final_examples]
    selected_statuses: list[str]

    if selected_model == "logistic_regression":
        selected_fitted_model = _fit_logistic(
            examples,
            development_indices,
            final_logistic_tuning.c,
            final_logistic_tuning.class_weight,
        )
        selected_scores = selected_fitted_model.predict_many(
            [example.features for example in final_examples]
        )
        selected_statuses = _logistic_statuses(
            selected_scores,
            conflict_threshold=final_logistic_tuning.conflict_threshold,
            support_threshold=final_logistic_tuning.support_threshold,
        )
        decision_thresholds = {
            "conflict_max_probability": final_logistic_tuning.conflict_threshold,
            "support_min_probability": final_logistic_tuning.support_threshold,
        }
        model_payload = {
            "type": "logistic_regression",
            **selected_fitted_model.to_dict(),
        }
    else:
        selected_statuses = [
            threshold_bucket(example.features, final_threshold_tuning.parameters)
            for example in final_examples
        ]
        selected_scores = [
            final_threshold_tuning.bucket_probabilities[
                threshold_bucket(example.features, final_threshold_tuning.parameters)
            ]
            for example in final_examples
        ]
        model_payload = {
            "type": "threshold",
            "parameters": final_threshold_tuning.parameters.to_dict(),
            "bucket_probabilities": dict(final_threshold_tuning.bucket_probabilities),
        }
        # Status comes directly from an interpretable rule bucket.  These
        # sentinels document that no probability cutoffs control that status.
        decision_thresholds = {
            "conflict_max_probability": 0.0,
            "support_min_probability": 1.0,
        }

    if selected_model == "threshold":
        locked_evaluation = _threshold_evaluation_record(
            final_labels,
            selected_statuses,
            selected_scores,
            final_groups,
            seed=seed + 900_001,
            bootstrap_iterations=bootstrap_iterations,
        )
    else:
        locked_evaluation = _evaluation_record(
            final_labels,
            selected_scores,
            final_groups,
            statuses=selected_statuses,
            decision_policy=(
                "probability <= conflict threshold is conflicting; probability >= "
                "support threshold is supported; the middle interval abstains"
            ),
            seed=seed + 900_001,
            bootstrap_iterations=bootstrap_iterations,
        )
    selected_model_training_ranges = _feature_ranges(_select(examples, development_indices))

    provenance = {
        "training_input": training_input_provenance,
        "source_sidecar": source_sidecar_provenance,
        "diagnostic_standard": diagnostic_standard_provenance,
        "seed": seed,
        "label_column": label_column,
        "compound_group_column": compound_column,
        "batch_group_column": batch_column,
        "row_id_column": row_id_column,
        "locked_final_batch": final_batch_name,
        "truth_contract": "independent_external_labels_required",
        "frozen_source_contract_validated": True,
        "batch_source_contract_validated": True,
        "batch_sources": batch_source_provenance,
    }
    decision_without_id: dict[str, Any] = {
        "schema_version": MODEL_SCHEMA_VERSION,
        "selected_model": selected_model,
        "feature_names": list(FEATURE_NAMES),
        "model": model_payload,
        "decision_thresholds": decision_thresholds,
        "training_feature_ranges": selected_model_training_ranges,
        "candidates": {
            "threshold": {
                "grid_size": final_threshold_tuning.candidates_evaluated,
                "final_development_tuning": final_threshold_tuning.to_dict(),
            },
            "logistic_regression": {
                "c_grid": list(C_GRID),
                "class_weight_grid": list(CLASS_WEIGHT_GRID),
                "final_development_tuning": final_logistic_tuning.to_dict(),
            },
        },
        "selection": reliability,
        "development_validation": {
            "strategy": "nested_leave_one_batch_out_with_compound_purge",
            "specificity_target": specificity_target,
            "outer_folds": outer_records,
            "pooled": {
                "threshold": _threshold_evaluation_record(
                    outer_labels,
                    threshold_outer_statuses,
                    threshold_outer_probabilities,
                    outer_groups,
                    seed=seed + 800_001,
                    bootstrap_iterations=bootstrap_iterations,
                ),
                "logistic_regression": _evaluation_record(
                    outer_labels,
                    logistic_outer_scores,
                    outer_groups,
                    statuses=logistic_outer_statuses,
                    decision_policy=(
                        "fold-specific nested-CV conflict/support thresholds; middle "
                        "probabilities abstain"
                    ),
                    seed=seed + 800_002,
                    bootstrap_iterations=bootstrap_iterations,
                ),
            },
        },
        "locked_final_evaluation": {
            "batch": final_batch_name,
            "used_once_after_model_selection": True,
            "evaluation": locked_evaluation,
        },
        "split_manifest": {
            "strategy": "locked_batch_then_nested_leave_one_batch_out",
            "eligible_batches": _group_values(examples, range(len(examples)), "batch_id"),
            "development_row_ids": _row_ids(examples, development_indices),
            "locked_final_row_ids": _row_ids(examples, final_indices),
            "locked_compound_overlap_purged_row_ids": _row_ids(
                examples,
                locked.purged_indices,
            ),
            "locked_compound_overlap_purged_batches": _group_values(
                examples,
                locked.purged_indices,
                "batch_id",
            ),
            "development_batches": _group_values(examples, development_indices, "batch_id"),
            "locked_final_batches": _group_values(examples, final_indices, "batch_id"),
            "development_compounds": _group_values(
                examples,
                development_indices,
                "compound_id",
            ),
            "locked_final_compounds": _group_values(examples, final_indices, "compound_id"),
            "outer_folds": split_manifest_outer,
            "final_development_folds": [
                _fold_manifest(examples, fold) for fold in final_tuning_folds
            ],
        },
        "annotation_summary": annotation_summary,
        "provenance": provenance,
        "model_card": build_model_card(
            selected_model=selected_model,
            annotation_summary=annotation_summary,
            development_validation=outer_records,
            locked_final_evaluation=locked_evaluation,
            provenance=provenance,
            provenance_discrepancies=provenance_discrepancies,
        ),
        "limitations": [
            "Shadow mode only; the existing v2 diagnostic result is preserved.",
            "No identity claim is valid beyond the independently labeled batches represented here.",
            "The locked final batch is evaluation-only and was not used for model selection.",
            "Upstream Method and spectrum-extraction settings are frozen outside this artifact.",
        ],
    }
    model_id = canonical_model_id(decision_without_id)
    artifact = make_trainable_artifact(
        {
            "model_id": model_id,
            **decision_without_id,
        }
    )
    # Apply the same deep production validator used by inference before the
    # locked-final artifact can be returned or written.
    from .shadow import validate_model_artifact

    validate_model_artifact(artifact)
    return artifact


def build_model_card(
    *,
    selected_model: str,
    annotation_summary: Mapping[str, int],
    development_validation: Sequence[Mapping[str, Any]],
    locked_final_evaluation: Mapping[str, Any],
    provenance: Mapping[str, Any],
    provenance_discrepancies: Sequence[str] = (),
) -> dict[str, Any]:
    """Build the model-card section embedded in every training artifact."""

    return {
        "intended_use": (
            "Shadow-mode prioritization of whether two independently extracted MS2 "
            "spectra support the same-compound hypothesis."
        ),
        "claim_boundary": (
            "This classifier does not assign absolute stereochemistry, replace the v2 "
            "diagnostic standard, or establish compound identity without orthogonal evidence."
        ),
        "selected_model": selected_model,
        "training_data": {
            "annotation_counts": dict(annotation_summary),
            "truth_source": provenance["training_input"],
            "grouping": {
                "compound_column": provenance["compound_group_column"],
                "batch_column": provenance["batch_group_column"],
                "locked_final_batch": provenance["locked_final_batch"],
            },
        },
        "metrics": {
            "nested_grouped_cv_batches": [
                {
                    "batch": record["validation_batch"],
                    "threshold": {
                        "all_row_support_call_metrics": record["threshold"]["evaluation"][
                            "all_row_support_call_metrics"
                        ],
                        "selective_non_abstained_metrics": record["threshold"]["evaluation"][
                            "selective_non_abstained_metrics"
                        ],
                        "emitted_probability_quality": record["threshold"]["evaluation"][
                            "probability_metrics"
                        ],
                        "emitted_probability_calibration": record["threshold"]["evaluation"][
                            "calibration"
                        ],
                    },
                    "logistic_regression": {
                        "all_row_support_call_metrics": record["logistic_regression"]["evaluation"][
                            "all_row_support_call_metrics"
                        ],
                        "selective_non_abstained_metrics": record["logistic_regression"][
                            "evaluation"
                        ]["selective_non_abstained_metrics"],
                        "probability_quality_all_rows": record["logistic_regression"]["evaluation"][
                            "probability_metrics"
                        ],
                        "calibration": record["logistic_regression"]["evaluation"]["calibration"],
                    },
                }
                for record in development_validation
            ],
            "locked_final": {
                "all_row_support_call_metrics": locked_final_evaluation[
                    "all_row_support_call_metrics"
                ],
                "selective_non_abstained_metrics": locked_final_evaluation[
                    "selective_non_abstained_metrics"
                ],
                "probability_quality_all_rows": locked_final_evaluation["probability_metrics"],
                "emitted_probability_calibration": locked_final_evaluation["calibration"],
            },
        },
        "abstention": {
            "status": _ABSTENTION_STATUS,
            "threshold_policy": (
                "threshold supported/conflicting buckets are decisions; "
                "insufficient_evidence abstains"
            ),
            "logistic_regression_policy": (
                "probability >= support boundary is supported; probability <= conflict "
                "boundary is conflicting; the middle interval abstains"
            ),
            "reporting_scope": (
                "coverage and selective metrics accompany all-row support-call metrics; "
                "PR-AUC, average precision, Brier score, ECE, and calibration use all "
                "evaluable probabilities including abstentions"
            ),
            "nested_grouped_cv_batches": [
                {
                    "batch": record["validation_batch"],
                    "threshold": _model_card_abstention(record["threshold"]["evaluation"]),
                    "logistic_regression": _model_card_abstention(
                        record["logistic_regression"]["evaluation"]
                    ),
                }
                for record in development_validation
            ],
            "locked_final": _model_card_abstention(locked_final_evaluation),
        },
        "limitations": [
            "Performance applies only to the independently labeled batches and chemistry represented.",
            "Middle-range Logistic Regression probabilities abstain as insufficient evidence.",
            "Missing bilateral spectra are not evaluable and are never negative examples.",
            "Threshold baseline is retained unless Logistic Regression improves every outer batch.",
        ],
        "provenance_discrepancies": [str(value) for value in provenance_discrepancies],
        "current_data_validation": (
            "not_validated_by_this_model_card; current unlabeled feasibility data cannot "
            "supply independent performance estimates"
        ),
    }


def train_and_write_shadow_model(
    rows: Sequence[Mapping[str, Any]],
    output_path: str | Path,
    **kwargs: Any,
) -> dict[str, dict[str, Any]]:
    """Train fully, then atomically write one validated deterministic artifact."""

    artifact = train_shadow_model(rows, **kwargs)
    from .shadow import validate_model_artifact

    validate_model_artifact(artifact)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    temporary.write_text(
        json.dumps(artifact, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(output)
    return artifact


def _prepare_examples(
    rows: Sequence[Mapping[str, Any]],
    *,
    label_column: str,
    compound_column: str,
    batch_column: str,
    row_id_column: str,
    spectrum_a_column: str,
    spectrum_b_column: str,
) -> tuple[list[GroupedExample], dict[str, int]]:
    if not rows or not any(label_column in row for row in rows):
        raise ValueError(
            "No independent truth annotations are available; feasibility data cannot train a model."
        )
    raw_labels = [str(row.get(label_column, "")).strip() for row in rows]
    if not any(raw_labels):
        raise ValueError(
            "No independent truth annotations are available; feasibility data cannot train a model."
        )

    counts = {
        "input_rows": len(rows),
        "positive_same_compound": 0,
        "negative_different_or_interference": 0,
        "uncertain": 0,
        "not_evaluable": 0,
        "eligible_rows": 0,
    }
    examples: list[GroupedExample] = []
    seen_row_ids: set[str] = set()
    for index, row in enumerate(rows):
        label_name = validate_truth_row(
            row,
            label_column,
            spectrum_a_column,
            spectrum_b_column,
        )
        counts[normalize_truth_label(label_name)] += 1
        binary_label = truth_to_binary(label_name)
        if binary_label is None:
            continue
        compound_id = str(row.get(compound_column, "")).strip()
        batch_id = str(row.get(batch_column, "")).strip()
        if not compound_id or not batch_id:
            raise ValueError(
                f"Eligible row {index} requires non-empty {compound_column!r} and {batch_column!r}."
            )
        row_id = str(row.get(row_id_column, "")).strip() or f"row-{index:06d}"
        if row_id in seen_row_ids:
            raise ValueError(f"Duplicate training row ID: {row_id!r}.")
        seen_row_ids.add(row_id)
        extracted = extract_features_from_row(row, spectrum_a_column, spectrum_b_column)
        if extracted is None:
            raise ValueError(
                f"Eligible truth row {row_id!r} has bilateral cells but no evaluable fixed features."
            )
        examples.append(
            GroupedExample(
                row_id=row_id,
                compound_id=compound_id,
                batch_id=batch_id,
                label=binary_label,
                features=extracted.as_vector(),
            )
        )
        counts["eligible_rows"] += 1
    if not examples or {example.label for example in examples} != {0, 1}:
        raise ValueError(
            "Training requires evaluable independent positive_same_compound and "
            "negative_different_or_interference labels."
        )
    return examples, counts


def _validate_batch_source_proof(
    source_sidecar: Mapping[str, Any],
    eligible_batches: Sequence[str],
    diagnostic_standard: Mapping[str, Any],
    *,
    single_source_path: str | Path | None,
) -> dict[str, dict[str, str]]:
    batches = sorted({str(batch) for batch in eligible_batches})
    if not batches:
        raise ValueError("At least one eligible batch is required for source validation.")
    manifest = source_sidecar.get("batch_sources")
    if not isinstance(manifest, Mapping):
        if len(batches) > 1:
            raise ValueError(
                "Multi-batch training requires source_sidecar.batch_sources with one "
                "embedded frozen-source sidecar for every eligible batch."
            )
        entry: dict[str, Any] = {"sidecar": source_sidecar}
        if single_source_path is not None:
            entry["path"] = str(single_source_path)
        return {
            batches[0]: _validate_batch_source_entry(
                batches[0],
                entry,
                diagnostic_standard,
                relative_to=None,
            )
        }

    normalized_manifest = {str(key): value for key, value in manifest.items()}
    if len(normalized_manifest) != len(manifest):
        raise ValueError("batch_sources contains duplicate batch IDs after string normalization.")
    manifest_batches = set(normalized_manifest)
    missing = sorted(set(batches).difference(manifest_batches))
    extra = sorted(manifest_batches.difference(batches))
    if missing or extra:
        details: list[str] = []
        if missing:
            details.append("missing eligible batches: " + ", ".join(missing))
        if extra:
            details.append("unexpected batches: " + ", ".join(extra))
        raise ValueError("batch_sources must exactly match eligible batches; " + "; ".join(details))
    return {
        batch: _validate_batch_source_entry(
            batch,
            normalized_manifest[batch],
            diagnostic_standard,
            relative_to=(
                Path(single_source_path).resolve().parent
                if single_source_path is not None
                else Path.cwd()
            ),
        )
        for batch in batches
    }


def _validate_batch_source_entry(
    batch: str,
    raw_entry: Any,
    diagnostic_standard: Mapping[str, Any],
    *,
    relative_to: Path | None,
) -> dict[str, str]:
    if not isinstance(raw_entry, Mapping):
        raise ValueError(f"batch_sources[{batch!r}] must be a mapping.")
    path_value = str(raw_entry.get("path") or "").strip()
    source_path: Path | None = None
    if path_value:
        candidate = Path(path_value)
        if not candidate.is_absolute() and relative_to is not None:
            candidate = relative_to / candidate
        source_path = candidate.resolve()
        if not source_path.is_file():
            raise ValueError(
                f"Batch source path for {batch!r} does not exist or is not a file: {source_path}"
            )

    embedded = raw_entry.get("sidecar")
    if embedded is None:
        direct_embedded = {
            key: value for key, value in raw_entry.items() if str(key) not in {"path", "sha256"}
        }
        embedded = direct_embedded or None
    if embedded is not None and not isinstance(embedded, Mapping):
        raise ValueError(f"batch_sources[{batch!r}].sidecar must be an embedded sidecar mapping.")
    if embedded is None and source_path is None:
        raise ValueError(f"batch_sources[{batch!r}] requires an embedded sidecar or a source path.")
    embedded_digest = ""
    if embedded is not None:
        validate_frozen_source(embedded, diagnostic_standard)
        embedded_digest = _sha256_json(embedded)

    declared_digest_value = raw_entry.get("sha256")
    declared_digest = ""
    if declared_digest_value is not None:
        declared_digest = str(declared_digest_value).strip().lower()
        if len(declared_digest) != 64 or any(
            character not in "0123456789abcdef" for character in declared_digest
        ):
            raise ValueError(f"batch_sources[{batch!r}].sha256 must be a 64-character hex digest.")

    if source_path is not None:
        with source_path.open(encoding="utf-8") as handle:
            file_sidecar = json.load(handle)
        if not isinstance(file_sidecar, Mapping):
            raise ValueError(f"Batch source file for {batch!r} must contain a JSON object.")
        validate_frozen_source(file_sidecar, diagnostic_standard)
        if embedded is not None and _sha256_json(file_sidecar) != embedded_digest:
            raise ValueError(
                f"Embedded sidecar for batch {batch!r} does not match its source file."
            )
        observed_digest = _sha256_file(source_path)
        proof_type = "source_file"
        resolved_path = str(source_path)
    else:
        if embedded is None:
            raise AssertionError("Embedded source unexpectedly missing.")
        observed_digest = embedded_digest
        proof_type = "embedded_canonical_json"
        resolved_path = path_value or f"<embedded:{batch}>"

    if declared_digest and declared_digest != observed_digest:
        raise ValueError(
            f"Source SHA-256 mismatch for batch {batch!r}: expected "
            f"{declared_digest}, observed {observed_digest}."
        )
    return {
        "path": resolved_path,
        "sha256": observed_digest,
        "proof_type": proof_type,
    }


def _nested_folds(
    examples: Sequence[GroupedExample],
    indices: Sequence[int],
    *,
    name_prefix: str,
) -> tuple[FoldSplit, ...]:
    folds = tuple(
        build_leave_one_batch_out_folds(
            examples,
            indices=tuple(indices),
            name_prefix=name_prefix,
        )
    )
    if len(folds) < 2:
        raise ValueError(
            f"{name_prefix} requires at least two batch-held-out folds after compound purging."
        )
    for fold in folds:
        validate_binary_fold(examples, fold)
        assert_group_isolation(examples, fold.train_indices, fold.validation_indices)
    return folds


def _tune_threshold_nested(
    examples: Sequence[GroupedExample],
    folds: Sequence[FoldSplit],
    specificity_target: float,
) -> ThresholdSearchResult:
    # A threshold rule has no fitted coefficients.  Pooling only inner
    # validation partitions evaluates every grid point out-of-batch while the
    # fold manifest proves batch and compound isolation.
    validation_indices = tuple(index for fold in folds for index in fold.validation_indices)
    if len(validation_indices) != len(set(validation_indices)):
        raise AssertionError("Nested validation rows must occur in exactly one fold.")
    return fit_thresholds(
        _select(examples, validation_indices),
        specificity_target=specificity_target,
    )


def _tune_logistic_nested(
    examples: Sequence[GroupedExample],
    folds: Sequence[FoldSplit],
    specificity_target: float,
) -> LogisticTuningResult:
    best: tuple[tuple[float, ...], LogisticTuningResult] | None = None
    candidate_count = 0
    converged_count = 0
    convergence_failures: list[str] = []
    for class_weight in CLASS_WEIGHT_GRID:
        for c in C_GRID:
            candidate_count += 1
            indexed_predictions: list[tuple[int, int, float]] = []
            candidate_failed = False
            for fold in folds:
                try:
                    model = _fit_logistic(examples, fold.train_indices, c, class_weight)
                except LogisticConvergenceError as error:
                    weight_name = "none" if class_weight is None else class_weight
                    convergence_failures.append(
                        f"C={c:g};class_weight={weight_name};fold={fold.name};error={error}"
                    )
                    candidate_failed = True
                    break
                for index, probability in zip(
                    fold.validation_indices,
                    model.predict_many(
                        [examples[row_index].features for row_index in fold.validation_indices]
                    ),
                    strict=True,
                ):
                    indexed_predictions.append((index, examples[index].label, probability))
            if candidate_failed:
                # A numerical candidate that does not converge in every inner
                # fold is not eligible for scoring or selection.
                continue
            converged_count += 1
            indexed_predictions.sort(key=lambda value: value[0])
            indices = tuple(index for index, _label, _probability in indexed_predictions)
            if len(indices) != len(set(indices)):
                raise AssertionError("Nested logistic OOF rows must be unique.")
            labels = tuple(label for _index, label, _probability in indexed_predictions)
            probabilities = tuple(
                probability for _index, _label, probability in indexed_predictions
            )
            try:
                support_threshold, metrics = select_support_probability_threshold(
                    labels,
                    probabilities,
                    specificity_target=specificity_target,
                )
            except ValueError:
                continue
            conflict_threshold = select_conflict_probability_threshold(
                labels,
                probabilities,
                support_threshold,
                positive_safety_target=specificity_target,
            )
            result = LogisticTuningResult(
                c=float(c),
                class_weight=class_weight,
                support_threshold=support_threshold,
                conflict_threshold=conflict_threshold,
                metrics=metrics.as_dict(),
                oof_labels=labels,
                oof_probabilities=probabilities,
                oof_indices=indices,
                candidates_evaluated=0,
                candidates_attempted=0,
                convergence_failures=(),
            )
            # Prefer constrained recall, then probability quality, then the
            # simpler unweighted/smaller-C solution for exact ties.
            key = (
                metrics.recall,
                metrics.precision,
                metrics.f0_5,
                metrics.average_precision,
                -metrics.brier_score,
                float(class_weight is None),
                -float(c),
            )
            if best is None or key > best[0]:
                best = key, result
    if best is None:
        raise ValueError(
            "No StandardScaler + L2 Logistic Regression configuration satisfies "
            f"specificity >= {specificity_target:.3f} in nested grouped CV."
        )
    selected = best[1]
    return LogisticTuningResult(
        c=selected.c,
        class_weight=selected.class_weight,
        support_threshold=selected.support_threshold,
        conflict_threshold=selected.conflict_threshold,
        metrics=selected.metrics,
        oof_labels=selected.oof_labels,
        oof_probabilities=selected.oof_probabilities,
        oof_indices=selected.oof_indices,
        candidates_evaluated=converged_count,
        candidates_attempted=candidate_count,
        convergence_failures=tuple(convergence_failures),
    )


def _fit_logistic(
    examples: Sequence[GroupedExample],
    indices: Sequence[int],
    c: float,
    class_weight: str | None,
) -> L2LogisticModel:
    selected = _select(examples, indices)
    model = fit_l2_logistic(
        [example.features for example in selected],
        [example.label for example in selected],
        c=c,
        class_weight=class_weight,
        feature_names=FEATURE_NAMES,
    )
    if not model.converged:
        raise LogisticConvergenceError(
            "StandardScaler + L2 Logistic Regression did not converge; "
            "the candidate is ineligible for training or deployment."
        )
    return model


def _logistic_reliability(
    outer_records: Sequence[Mapping[str, Any]],
    specificity_target: float,
) -> dict[str, Any]:
    batch_deltas: list[dict[str, Any]] = []
    for record in outer_records:
        threshold_metrics = record["threshold"]["evaluation"]["metrics"]
        logistic_metrics = record["logistic_regression"]["evaluation"]["metrics"]
        batch_deltas.append(
            {
                "batch": record["validation_batch"],
                "recall_delta": logistic_metrics["recall"] - threshold_metrics["recall"],
                "threshold_recall": threshold_metrics["recall"],
                "logistic_recall": logistic_metrics["recall"],
                "logistic_specificity": logistic_metrics["specificity"],
            }
        )
    reliable = (
        len(batch_deltas) >= 2
        and all(item["logistic_specificity"] >= specificity_target for item in batch_deltas)
        and all(item["recall_delta"] > 0 for item in batch_deltas)
    )
    return {
        "selected_model": "logistic_regression" if reliable else "threshold",
        "threshold_wins_ties": True,
        "criterion": (
            "logistic specificity must meet the target and logistic recall must be "
            "strictly greater than threshold recall in every outer validation batch"
        ),
        "reliably_better": reliable,
        "specificity_target": specificity_target,
        "batch_deltas": batch_deltas,
    }


def _evaluation_record(
    labels: Sequence[int],
    scores: Sequence[float],
    groups: Sequence[str],
    *,
    statuses: Sequence[str],
    decision_policy: str,
    seed: int,
    bootstrap_iterations: int,
) -> dict[str, Any]:
    if not (len(labels) == len(scores) == len(groups) == len(statuses)) or not labels:
        raise ValueError(
            "Evaluation labels, probabilities, groups, and statuses must have equal "
            "non-zero length."
        )
    normalized_statuses = tuple(str(status) for status in statuses)
    unexpected = sorted(set(normalized_statuses).difference(_TRAINING_DECISION_STATUSES))
    if unexpected:
        raise ValueError("Unexpected training evaluation statuses: " + ", ".join(unexpected))
    support_calls = tuple(float(status == _SUPPORT_STATUS) for status in normalized_statuses)
    support_metrics = binary_metrics(labels, support_calls, threshold=0.5)
    probability_metrics = binary_metrics(labels, scores, threshold=0.5)
    combined_metrics = support_metrics.as_dict()
    probability_metric_values = probability_metrics.as_dict()
    for name in _PROBABILITY_METRIC_NAMES:
        combined_metrics[name] = probability_metric_values[name]
    abstention = _abstention_report(labels, normalized_statuses)
    abstention_bootstrap = _abstention_group_bootstrap_ci(
        labels,
        normalized_statuses,
        groups,
        iterations=bootstrap_iterations,
        seed=seed + 2,
    )
    return {
        "decision_policy": decision_policy,
        "metrics": combined_metrics,
        "metrics_note": (
            "Confusion, precision, recall, specificity, and F0.5 are all-row "
            "support-call metrics. PR-AUC, average precision, Brier score, and ECE "
            "use all emitted probabilities, including abstentions."
        ),
        "all_row_support_call_metrics": _operating_metrics(support_metrics.as_dict()),
        "all_row_support_call_metrics_note": (
            "All eligible rows remain in the denominator; supported_same_compound is "
            "the positive call, while conflicting_spectra and insufficient_evidence "
            "are non-support calls. Consult selective metrics and coverage before "
            "interpreting recall."
        ),
        "probability_metrics": probability_metric_values,
        "probability_metrics_note": (
            "PR-AUC, average precision, Brier score, ECE, and calibration are computed "
            "from every evaluable probability, including rows whose status abstains."
        ),
        **abstention,
        "calibration": calibration_table(labels, scores),
        "bootstrap_95_ci": group_bootstrap_ci(
            labels,
            support_calls,
            groups,
            threshold=0.5,
            iterations=bootstrap_iterations,
            seed=seed,
        ),
        "bootstrap_note": (
            "This group bootstrap follows the all-row binary support call; use only "
            "its operating-point precision/recall/specificity/F0.5 intervals."
        ),
        "probability_bootstrap_95_ci": group_bootstrap_ci(
            labels,
            scores,
            groups,
            threshold=0.5,
            iterations=bootstrap_iterations,
            seed=seed + 1,
        ),
        "probability_bootstrap_seed": seed + 1,
        "probability_bootstrap_note": (
            "This all-row probability bootstrap is authoritative for PR-AUC, average "
            "precision, Brier score, and ECE; its 0.5 confusion cutoff is descriptive only."
        ),
        "abstention_bootstrap_95_ci": abstention_bootstrap,
        "selective_non_abstained_bootstrap_95_ci": abstention_bootstrap[
            "selective_non_abstained_metrics"
        ],
        "abstention_bootstrap_seed": seed + 2,
        "abstention_bootstrap_note": (
            "Compound_ID groups are resampled as whole units. Each interval records "
            "its actual successful replicate count; class-specific replicates lacking "
            "that truth class are skipped only for that class."
        ),
        "bootstrap_seed": seed,
        "bootstrap_iterations": bootstrap_iterations,
        "row_count": len(labels),
        "compound_count": len(set(groups)),
    }


def _threshold_evaluation_record(
    labels: Sequence[int],
    statuses: Sequence[str],
    probabilities: Sequence[float],
    groups: Sequence[str],
    *,
    seed: int,
    bootstrap_iterations: int,
) -> dict[str, Any]:
    """Evaluate rule statuses while calibrating all emitted probabilities."""

    return _evaluation_record(
        labels,
        probabilities,
        groups,
        statuses=statuses,
        decision_policy=(
            "supported_same_compound and conflicting_spectra are decisions; "
            "insufficient_evidence abstains"
        ),
        seed=seed,
        bootstrap_iterations=bootstrap_iterations,
    )


def _logistic_statuses(
    probabilities: Sequence[float],
    *,
    conflict_threshold: float,
    support_threshold: float,
) -> list[str]:
    if not 0 <= conflict_threshold < support_threshold <= 1:
        raise ValueError("Logistic conflict threshold must be below support threshold.")
    return [
        _CONFLICT_STATUS
        if probability <= conflict_threshold
        else _SUPPORT_STATUS
        if probability >= support_threshold
        else _ABSTENTION_STATUS
        for probability in probabilities
    ]


def _abstention_report(
    labels: Sequence[int],
    statuses: Sequence[str],
) -> dict[str, Any]:
    if len(labels) != len(statuses) or not labels:
        raise ValueError("Abstention reporting requires equal non-empty labels and statuses.")
    status_truth_counts = {
        status: {
            POSITIVE_TRUTH_LABEL: 0,
            NEGATIVE_TRUTH_LABEL: 0,
        }
        for status in _TRAINING_DECISION_STATUSES
    }
    for label, status in zip(labels, statuses, strict=True):
        truth_name = POSITIVE_TRUTH_LABEL if label == 1 else NEGATIVE_TRUTH_LABEL
        status_truth_counts[status][truth_name] += 1

    abstained = sum(status == _ABSTENTION_STATUS for status in statuses)
    row_count = len(labels)
    class_specific: dict[str, dict[str, float | int]] = {}
    for label_value, truth_name in (
        (1, POSITIVE_TRUTH_LABEL),
        (0, NEGATIVE_TRUTH_LABEL),
    ):
        class_indices = [index for index, label in enumerate(labels) if label == label_value]
        class_abstained = sum(statuses[index] == _ABSTENTION_STATUS for index in class_indices)
        class_rows = len(class_indices)
        class_specific[truth_name] = {
            "row_count": class_rows,
            "decided_count": class_rows - class_abstained,
            "abstained_count": class_abstained,
            "decision_coverage": _safe_divide(class_rows - class_abstained, class_rows),
            "abstention_rate": _safe_divide(class_abstained, class_rows),
        }

    selective_indices = [
        index for index, status in enumerate(statuses) if status != _ABSTENTION_STATUS
    ]
    selective_labels = [labels[index] for index in selective_indices]
    selective_predictions = [int(statuses[index] == _SUPPORT_STATUS) for index in selective_indices]
    return {
        "status_truth_counts": status_truth_counts,
        "status_rates": {
            status: sum(observed == status for observed in statuses) / row_count
            for status in _TRAINING_DECISION_STATUSES
        },
        "decision_coverage": (row_count - abstained) / row_count,
        "abstention_rate": abstained / row_count,
        "class_specific_coverage": class_specific,
        "selective_non_abstained_metrics": _selective_operating_metrics(
            selective_labels,
            selective_predictions,
        ),
        "selective_non_abstained_metrics_note": (
            "Confusion, precision, recall, and specificity exclude only "
            "insufficient_evidence abstentions; coverage reports the excluded fraction."
        ),
    }


def _selective_operating_metrics(
    labels: Sequence[int],
    predictions: Sequence[int],
) -> dict[str, float | int]:
    tp = fp = tn = fn = 0
    for label, prediction in zip(labels, predictions, strict=True):
        if prediction == 1 and label == 1:
            tp += 1
        elif prediction == 1:
            fp += 1
        elif label == 0:
            tn += 1
        else:
            fn += 1
    precision = _safe_divide(tp, tp + fp)
    recall = _safe_divide(tp, tp + fn)
    return {
        "evaluable_rows": len(labels),
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "specificity": _safe_divide(tn, tn + fp),
        "f0_5": _safe_divide(1.25 * precision * recall, 0.25 * precision + recall),
    }


def _operating_metrics(metrics: Mapping[str, Any]) -> dict[str, Any]:
    return {
        name: metrics[name]
        for name in (
            "evaluable_rows",
            "tp",
            "fp",
            "tn",
            "fn",
            "precision",
            "recall",
            "specificity",
            "f0_5",
        )
    }


def _model_card_abstention(evaluation: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "status_truth_counts": evaluation["status_truth_counts"],
        "status_rates": evaluation["status_rates"],
        "decision_coverage": evaluation["decision_coverage"],
        "abstention_rate": evaluation["abstention_rate"],
        "class_specific_coverage": evaluation["class_specific_coverage"],
        "selective_non_abstained_metrics": evaluation["selective_non_abstained_metrics"],
        "selective_non_abstained_bootstrap_95_ci": evaluation[
            "selective_non_abstained_bootstrap_95_ci"
        ],
        "compound_group_bootstrap_95_ci": evaluation["abstention_bootstrap_95_ci"],
        "bootstrap_seed": evaluation["abstention_bootstrap_seed"],
        "bootstrap_iterations": evaluation["bootstrap_iterations"],
    }


def _abstention_group_bootstrap_ci(
    labels: Sequence[int],
    statuses: Sequence[str],
    groups: Sequence[str],
    *,
    iterations: int,
    seed: int,
) -> dict[str, Any]:
    if len(groups) != len(labels):
        raise ValueError("Abstention bootstrap groups must match evaluation rows.")
    by_group: dict[str, list[int]] = {}
    for index, group in enumerate(groups):
        by_group.setdefault(str(group), []).append(index)
    ordered_groups = sorted(by_group)
    if len(ordered_groups) < 2:
        raise ValueError("Abstention group bootstrap requires at least two compounds.")

    base = _abstention_report(labels, statuses)
    overall_samples = {"decision_coverage": [], "abstention_rate": []}
    status_samples = {status: [] for status in _TRAINING_DECISION_STATUSES}
    class_samples = {
        truth_name: {"decision_coverage": [], "abstention_rate": []}
        for truth_name in (POSITIVE_TRUTH_LABEL, NEGATIVE_TRUTH_LABEL)
    }
    selective_samples = {name: [] for name in ("precision", "recall", "specificity", "f0_5")}
    rng = random.Random(seed)
    for _iteration in range(iterations):
        chosen = [rng.choice(ordered_groups) for _group in ordered_groups]
        indices = [index for group in chosen for index in by_group[group]]
        replicate_labels = [labels[index] for index in indices]
        replicate_statuses = [statuses[index] for index in indices]
        report = _abstention_report(replicate_labels, replicate_statuses)
        for name in overall_samples:
            overall_samples[name].append(float(report[name]))
        for status in _TRAINING_DECISION_STATUSES:
            status_samples[status].append(float(report["status_rates"][status]))
        for label_value, truth_name in (
            (1, POSITIVE_TRUTH_LABEL),
            (0, NEGATIVE_TRUTH_LABEL),
        ):
            if label_value not in replicate_labels:
                continue
            for name in class_samples[truth_name]:
                class_samples[truth_name][name].append(
                    float(report["class_specific_coverage"][truth_name][name])
                )
        selective = report["selective_non_abstained_metrics"]
        if selective["tp"] + selective["fn"] > 0 and selective["tn"] + selective["fp"] > 0:
            for name in selective_samples:
                selective_samples[name].append(float(selective[name]))

    return {
        **{
            name: _bootstrap_interval(float(base[name]), values, iterations)
            for name, values in overall_samples.items()
        },
        "status_rates": {
            status: _bootstrap_interval(
                float(base["status_rates"][status]),
                values,
                iterations,
            )
            for status, values in status_samples.items()
        },
        "class_specific_coverage": {
            truth_name: {
                name: _bootstrap_interval(
                    float(base["class_specific_coverage"][truth_name][name]),
                    values,
                    iterations,
                )
                for name, values in metrics.items()
            }
            for truth_name, metrics in class_samples.items()
        },
        "selective_non_abstained_metrics": {
            name: _bootstrap_interval(
                float(base["selective_non_abstained_metrics"][name]),
                values,
                iterations,
            )
            for name, values in selective_samples.items()
        },
    }


def _bootstrap_interval(
    estimate: float,
    values: Sequence[float],
    requested_iterations: int,
) -> dict[str, float | int | None]:
    if not values:
        lower = upper = None
    else:
        lower = _quantile(values, 0.025)
        upper = _quantile(values, 0.975)
    return {
        "estimate": estimate,
        "lower": lower,
        "upper": upper,
        "confidence_level": 0.95,
        "successful_iterations": len(values),
        "requested_iterations": requested_iterations,
    }


def _quantile(values: Sequence[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower_index = math.floor(position)
    upper_index = math.ceil(position)
    if lower_index == upper_index:
        return float(ordered[lower_index])
    fraction = position - lower_index
    return float(ordered[lower_index] * (1 - fraction) + ordered[upper_index] * fraction)


def _safe_divide(numerator: float | int, denominator: float | int) -> float:
    return 0.0 if denominator == 0 else numerator / denominator


def _fold_manifest(
    examples: Sequence[GroupedExample],
    fold: FoldSplit,
) -> dict[str, Any]:
    return {
        "name": fold.name,
        "train_row_ids": _row_ids(examples, fold.train_indices),
        "validation_row_ids": _row_ids(examples, fold.validation_indices),
        "purged_row_ids": _row_ids(examples, fold.purged_indices),
        "train_batches": _group_values(examples, fold.train_indices, "batch_id"),
        "validation_batches": _group_values(
            examples,
            fold.validation_indices,
            "batch_id",
        ),
        "train_compounds": _group_values(examples, fold.train_indices, "compound_id"),
        "validation_compounds": _group_values(
            examples,
            fold.validation_indices,
            "compound_id",
        ),
    }


def _validate_locked_split(
    examples: Sequence[GroupedExample],
    development_indices: Sequence[int],
    final_indices: Sequence[int],
) -> None:
    if not development_indices or not final_indices:
        raise ValueError("Locked split requires non-empty development and final sets.")
    development = _select(examples, development_indices)
    final = _select(examples, final_indices)
    if {example.label for example in development} != {0, 1}:
        raise ValueError("Development data must contain both eligible truth classes.")
    if {example.label for example in final} != {0, 1}:
        raise ValueError("Locked final batch must contain both eligible truth classes.")
    development_batches = {example.batch_id for example in development}
    final_batches = {example.batch_id for example in final}
    development_compounds = {example.compound_id for example in development}
    final_compounds = {example.compound_id for example in final}
    if development_batches & final_batches or development_compounds & final_compounds:
        raise AssertionError("Locked split must be disjoint by both batch and compound.")


def _single_validation_batch(
    examples: Sequence[GroupedExample],
    fold: FoldSplit,
) -> str:
    batches = {examples[index].batch_id for index in fold.validation_indices}
    if len(batches) != 1:
        raise AssertionError("Leave-one-batch-out fold must have one validation batch.")
    return next(iter(batches))


def _feature_ranges(examples: Sequence[GroupedExample]) -> dict[str, dict[str, float]]:
    return {
        name: {
            "minimum": min(example.features[index] for example in examples),
            "maximum": max(example.features[index] for example in examples),
        }
        for index, name in enumerate(FEATURE_NAMES)
    }


def _select(
    examples: Sequence[GroupedExample],
    indices: Sequence[int],
) -> list[GroupedExample]:
    return [examples[index] for index in indices]


def _row_ids(
    examples: Sequence[GroupedExample],
    indices: Sequence[int],
) -> list[str]:
    return [examples[index].row_id for index in indices]


def _group_values(
    examples: Sequence[GroupedExample],
    indices: Sequence[int],
    attribute: str,
) -> list[str]:
    return sorted({str(getattr(examples[index], attribute)) for index in indices})


def _provenance_entry(
    path: str | Path | None,
    content: Any,
    *,
    source_kind: str,
) -> dict[str, str]:
    if path is None:
        return {"path": "<in-memory>", "sha256": _sha256_json(content)}

    source = Path(path).resolve()
    if not source.is_file():
        raise ValueError(f"Provenance source does not exist or is not a file: {source}")
    if source_kind == "table":
        from ..ms2.tables import read_table

        observed: Any = read_table(source)
        expected = [dict(row) for row in content]
    elif source_kind == "json":
        with source.open(encoding="utf-8") as handle:
            observed = json.load(handle)
        expected = content
    else:  # pragma: no cover - all callers use one of the two fixed kinds
        raise AssertionError(f"Unsupported provenance source kind: {source_kind}")
    if observed != expected:
        raise ValueError(
            f"Provenance source content does not match the in-memory {source_kind} data: {source}"
        )
    return {"path": str(source), "sha256": _sha256_file(source)}


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_json(value: Any) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _validate_training_options(
    final_batch: str,
    seed: int,
    specificity_target: float,
    bootstrap_iterations: int,
) -> None:
    if not str(final_batch).strip():
        raise ValueError("A locked final batch must be declared before training.")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("seed must be a non-negative integer.")
    if (
        isinstance(specificity_target, bool)
        or not isinstance(specificity_target, int | float)
        or not math.isfinite(specificity_target)
        or not MIN_SPECIFICITY_TARGET <= specificity_target <= 1
    ):
        raise ValueError(
            "specificity_target must be in [0.95, 1]; training may not relax the "
            "prespecified minimum-specificity safety target."
        )
    if (
        isinstance(bootstrap_iterations, bool)
        or not isinstance(bootstrap_iterations, int)
        or bootstrap_iterations < 1
    ):
        raise ValueError("bootstrap_iterations must be a positive integer.")


__all__ = [
    "DEFAULT_BOOTSTRAP_ITERATIONS",
    "DEFAULT_SEED",
    "DEFAULT_SPECIFICITY_TARGET",
    "MODEL_SCHEMA_VERSION",
    "LogisticTuningResult",
    "build_model_card",
    "train_and_write_shadow_model",
    "train_shadow_model",
]
