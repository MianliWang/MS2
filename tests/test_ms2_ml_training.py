import csv
import hashlib
import json
import math
import tempfile
import unittest
from pathlib import Path
from typing import Any
from unittest.mock import patch

from MSAI.python.ms2_ml.contracts import FROZEN_SOURCE_PARAMETERS
from MSAI.python.ms2_ml.features import FEATURE_NAMES
from MSAI.python.ms2_ml.logistic import (
    C_GRID,
    CLASS_WEIGHT_GRID,
    L2LogisticModel,
    LogisticConvergenceError,
    fit_l2_logistic,
)
from MSAI.python.ms2_ml.metrics import (
    binary_metrics,
    calibration_table,
    group_bootstrap_ci,
    select_threshold_for_specificity,
)
from MSAI.python.ms2_ml.shadow import predict_shadow
from MSAI.python.ms2_ml.splits import GroupedExample
from MSAI.python.ms2_ml.thresholds import fit_thresholds, iter_threshold_grid
from MSAI.python.ms2_ml.training import (
    _evaluation_record,
    train_and_write_shadow_model,
    train_shadow_model,
)

POSITIVE_A = "100,1000;120,900;140,800;160,700;180,600;200,500"
POSITIVE_B = "100.002,950;120.002,850;140.002,750;160.002,650;180.002,550;200.002,450"
NEGATIVE_A = POSITIVE_A
NEGATIVE_B = "101,1000;121,900;141,800;161,700;181,600;201,500"


def _standard():
    return json.loads(
        Path("MSAI/standards/ms2_diagnostic_standard_v2.json").read_text(encoding="utf-8")
    )


def _sidecar():
    return {"parameters": dict(FROZEN_SOURCE_PARAMETERS)}


def _batch_source_manifest(
    batches=("BATCH-A", "BATCH-B", "BATCH-C", "BATCH-FINAL"),
) -> dict[str, Any]:
    return {"batch_sources": {batch: {"sidecar": _sidecar()} for batch in batches}}


def _independent_rows():
    rows = []
    for batch in ("BATCH-A", "BATCH-B", "BATCH-C", "BATCH-FINAL"):
        rows.extend(
            (
                {
                    "target_uid": f"{batch}-positive",
                    "Compound_ID": f"{batch}-P",
                    "batch_id": batch,
                    "truth_label": "positive_same_compound",
                    "peak_a_MS2": POSITIVE_A,
                    "peak_b_MS2": POSITIVE_B,
                },
                {
                    "target_uid": f"{batch}-negative",
                    "Compound_ID": f"{batch}-N",
                    "batch_id": batch,
                    "truth_label": "negative_different_or_interference",
                    "peak_a_MS2": NEGATIVE_A,
                    "peak_b_MS2": NEGATIVE_B,
                },
            )
        )
    return rows


class ThresholdAndLogisticTests(unittest.TestCase):
    def test_threshold_grid_is_the_prespecified_6384_combinations(self):
        grid = tuple(iter_threshold_grid())
        self.assertEqual(len(grid), 19 * 7 * 6 * 8)
        self.assertEqual(grid[0].min_cosine, 0.5)
        self.assertEqual(grid[-1].min_cosine, 0.95)
        self.assertEqual(grid[0].min_matched_peaks, 2)
        self.assertEqual(grid[-1].min_matched_peaks, 10)
        self.assertIsNone(grid[0].min_entropy_similarity)
        self.assertEqual(grid[-1].min_entropy_similarity, 0.9)

    def test_threshold_search_maximizes_recall_under_specificity_constraint(self):
        examples = (
            GroupedExample("p1", "P1", "A", 1, (0.95, 0.95, math.log1p(8), 0.9, 0.0, 1.0)),
            GroupedExample("p2", "P2", "B", 1, (0.85, 0.85, math.log1p(6), 0.8, 0.0, 1.0)),
            GroupedExample("n1", "N1", "A", 0, (0.40, 0.40, math.log1p(1), 0.2, 0.5, 0.5)),
            GroupedExample("n2", "N2", "B", 0, (0.45, 0.45, math.log1p(2), 0.2, 0.5, 0.5)),
        )
        result = fit_thresholds(examples, specificity_target=0.95)
        self.assertEqual(result.candidates_evaluated, 6384)
        self.assertEqual(result.metrics.recall, 1.0)
        self.assertEqual(result.metrics.specificity, 1.0)
        self.assertEqual(
            set(result.bucket_probabilities),
            {"supported_same_compound", "conflicting_spectra", "insufficient_evidence"},
        )
        with self.assertRaisesRegex(
            ValueError, r"Specificity/safety target must be in \[0.95, 1\]"
        ):
            fit_thresholds(examples, specificity_target=0.949)

    def test_standard_scaler_l2_logistic_is_serializable_and_deterministic(self):
        rows = (
            (0.95, 0.95, math.log1p(8), 0.9, 0.0, 1.0),
            (0.85, 0.85, math.log1p(6), 0.8, 0.05, 0.9),
            (0.20, 0.20, math.log1p(1), 0.1, 0.8, 0.2),
            (0.30, 0.30, math.log1p(2), 0.2, 0.7, 0.3),
        )
        labels = (1, 1, 0, 0)
        first = fit_l2_logistic(rows, labels, c=1.0, class_weight="balanced")
        second = fit_l2_logistic(rows, labels, c=1.0, class_weight="balanced")
        self.assertEqual(first, second)
        self.assertEqual(tuple(first.to_dict()["feature_names"]), FEATURE_NAMES)
        restored = L2LogisticModel.from_dict(first.to_dict())
        self.assertEqual(restored, first)
        probabilities = first.predict_many(rows)
        self.assertGreater(min(probabilities[:2]), max(probabilities[2:]))
        self.assertEqual(C_GRID, (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0))
        self.assertEqual(CLASS_WEIGHT_GRID, (None, "balanced"))
        with self.assertRaises(ValueError):
            fit_l2_logistic(rows, labels, c=0.05)
        with self.assertRaisesRegex(LogisticConvergenceError, "did not converge"):
            fit_l2_logistic(
                rows,
                labels,
                c=1.0,
                class_weight="balanced",
                max_iter=1,
                tolerance=1e-15,
            )


class MetricTests(unittest.TestCase):
    def test_known_probability_vector_reports_all_required_metrics(self):
        labels = (1, 1, 0, 0)
        scores = (0.9, 0.8, 0.2, 0.1)
        result = binary_metrics(labels, scores)
        self.assertEqual((result.tp, result.fp, result.tn, result.fn), (2, 0, 2, 0))
        self.assertEqual(result.precision, 1.0)
        self.assertEqual(result.recall, 1.0)
        self.assertEqual(result.specificity, 1.0)
        self.assertEqual(result.average_precision, 1.0)
        self.assertEqual(result.pr_auc, 1.0)
        self.assertEqual(result.f0_5, 1.0)
        self.assertAlmostEqual(result.brier_score, 0.025)
        self.assertTrue(calibration_table(labels, scores))

    def test_public_specificity_selector_cannot_relax_the_model_family_floor(self):
        with self.assertRaisesRegex(ValueError, r"minimum_specificity must be in \[0.95, 1\]"):
            select_threshold_for_specificity(
                (1, 1, 0, 0),
                (0.9, 0.8, 0.2, 0.1),
                minimum_specificity=0.949,
            )

    def test_group_bootstrap_ci_is_seed_deterministic(self):
        labels = (1, 1, 0, 0)
        scores = (0.9, 0.8, 0.2, 0.1)
        groups = ("P1", "P2", "N1", "N2")
        first = group_bootstrap_ci(labels, scores, groups, iterations=40, seed=7)
        second = group_bootstrap_ci(labels, scores, groups, iterations=40, seed=7)
        self.assertEqual(first, second)
        self.assertIsInstance(first["recall"]["successful_iterations"], int)
        self.assertLessEqual(first["recall"]["lower"], first["recall"]["upper"])

    def test_abstention_report_separates_all_row_and_selective_metrics(self):
        result = _evaluation_record(
            labels=(1, 1, 0, 0),
            scores=(0.9, 0.6, 0.1, 0.4),
            groups=("P1", "P2", "N1", "N2"),
            statuses=(
                "supported_same_compound",
                "insufficient_evidence",
                "conflicting_spectra",
                "insufficient_evidence",
            ),
            decision_policy="test policy",
            seed=17,
            bootstrap_iterations=40,
        )
        self.assertEqual(result["decision_coverage"], 0.5)
        self.assertEqual(result["abstention_rate"], 0.5)
        self.assertEqual(
            result["status_truth_counts"]["insufficient_evidence"],
            {
                "positive_same_compound": 1,
                "negative_different_or_interference": 1,
            },
        )
        self.assertEqual(
            result["class_specific_coverage"]["positive_same_compound"]["decision_coverage"],
            0.5,
        )
        self.assertEqual(
            result["all_row_support_call_metrics"]["recall"],
            0.5,
        )
        selective = result["selective_non_abstained_metrics"]
        self.assertEqual((selective["tp"], selective["tn"]), (1, 1))
        self.assertEqual(selective["recall"], 1.0)
        self.assertEqual(selective["specificity"], 1.0)
        self.assertAlmostEqual(result["probability_metrics"]["brier_score"], 0.085)
        bootstrap = result["abstention_bootstrap_95_ci"]
        self.assertEqual(bootstrap["decision_coverage"]["successful_iterations"], 40)
        self.assertEqual(bootstrap["decision_coverage"]["requested_iterations"], 40)
        self.assertGreater(
            bootstrap["class_specific_coverage"]["positive_same_compound"]["decision_coverage"][
                "successful_iterations"
            ],
            0,
        )
        selective_bootstrap = result["selective_non_abstained_bootstrap_95_ci"]
        self.assertEqual(selective_bootstrap["recall"]["estimate"], selective["recall"])
        self.assertEqual(
            selective_bootstrap,
            bootstrap["selective_non_abstained_metrics"],
        )
        self.assertGreater(selective_bootstrap["specificity"]["successful_iterations"], 0)


class EndToEndTrainingTests(unittest.TestCase):
    def test_nested_grouped_training_is_deterministic_and_threshold_wins_ties(self):
        rows = _independent_rows()
        first = train_shadow_model(
            rows,
            final_batch="BATCH-FINAL",
            source_sidecar=_batch_source_manifest(),
            diagnostic_standard=_standard(),
            bootstrap_iterations=20,
        )
        second = train_shadow_model(
            rows,
            final_batch="BATCH-FINAL",
            source_sidecar=_batch_source_manifest(),
            diagnostic_standard=_standard(),
            bootstrap_iterations=20,
        )
        self.assertEqual(first, second)
        self.assertEqual(set(first), {"trainable_decision"})
        decision = first["trainable_decision"]
        self.assertEqual(decision["selected_model"], "threshold")
        self.assertTrue(decision["selection"]["threshold_wins_ties"])
        self.assertFalse(decision["selection"]["reliably_better"])
        self.assertEqual(decision["candidates"]["threshold"]["grid_size"], 6384)
        self.assertEqual(decision["candidates"]["logistic_regression"]["c_grid"], list(C_GRID))
        self.assertEqual(decision["split_manifest"]["locked_final_batches"], ["BATCH-FINAL"])
        self.assertFalse(
            set(decision["split_manifest"]["development_compounds"])
            & set(decision["split_manifest"]["locked_final_compounds"])
        )
        for fold in decision["split_manifest"]["outer_folds"]:
            self.assertFalse(set(fold["train_batches"]) & set(fold["validation_batches"]))
            self.assertFalse(set(fold["train_compounds"]) & set(fold["validation_compounds"]))
        metrics = decision["locked_final_evaluation"]["evaluation"]["metrics"]
        for name in (
            "precision",
            "recall",
            "specificity",
            "pr_auc",
            "average_precision",
            "f0_5",
            "brier_score",
            "expected_calibration_error",
        ):
            self.assertIn(name, metrics)
        self.assertIn("model_card", decision)
        self.assertIn("limitations", decision["model_card"])
        self.assertIn("abstention", decision["model_card"])
        evaluation_records = [
            candidate["evaluation"]
            for fold in decision["development_validation"]["outer_folds"]
            for candidate in (fold["threshold"], fold["logistic_regression"])
        ]
        evaluation_records.extend(decision["development_validation"]["pooled"].values())
        evaluation_records.append(decision["locked_final_evaluation"]["evaluation"])
        for evaluation in evaluation_records:
            self.assertIn("status_truth_counts", evaluation)
            self.assertIn("decision_coverage", evaluation)
            self.assertIn("abstention_rate", evaluation)
            self.assertIn("class_specific_coverage", evaluation)
            self.assertIn("all_row_support_call_metrics", evaluation)
            self.assertIn("selective_non_abstained_metrics", evaluation)
            self.assertIn("selective_non_abstained_bootstrap_95_ci", evaluation)
            self.assertIn("abstention_bootstrap_95_ci", evaluation)
            self.assertAlmostEqual(
                evaluation["decision_coverage"] + evaluation["abstention_rate"],
                1.0,
            )
        self.assertNotIn("created_utc", json.dumps(first))

    def test_selected_logistic_branch_produces_a_loadable_shadow_artifact(self):
        def force_consistent_logistic_selection(outer_records, specificity_target):
            deltas = []
            for record in outer_records:
                threshold = record["threshold"]["evaluation"]["metrics"]
                logistic = record["logistic_regression"]["evaluation"]["metrics"]
                positive_count = threshold["tp"] + threshold["fn"]
                negative_count = logistic["tn"] + logistic["fp"]
                threshold.update(
                    {
                        "tp": 0,
                        "fn": positive_count,
                        "precision": 0.0,
                        "recall": 0.0,
                        "f0_5": 0.0,
                    }
                )
                logistic.update(
                    {
                        "tp": positive_count,
                        "fn": 0,
                        "tn": negative_count,
                        "fp": 0,
                        "precision": 1.0,
                        "recall": 1.0,
                        "specificity": 1.0,
                        "f0_5": 1.0,
                    }
                )
                deltas.append(
                    {
                        "batch": record["validation_batch"],
                        "recall_delta": 1.0,
                        "threshold_recall": 0.0,
                        "logistic_recall": 1.0,
                        "logistic_specificity": 1.0,
                    }
                )
            return {
                "selected_model": "logistic_regression",
                "threshold_wins_ties": True,
                "criterion": "test-only forced branch with internally bound outer evidence",
                "reliably_better": True,
                "specificity_target": specificity_target,
                "batch_deltas": deltas,
            }

        with patch(
            "MSAI.python.ms2_ml.training._logistic_reliability",
            side_effect=force_consistent_logistic_selection,
        ):
            artifact = train_shadow_model(
                _independent_rows(),
                final_batch="BATCH-FINAL",
                source_sidecar=_batch_source_manifest(),
                diagnostic_standard=_standard(),
                bootstrap_iterations=5,
            )
        self.assertEqual(
            artifact["trainable_decision"]["selected_model"],
            "logistic_regression",
        )
        prediction = predict_shadow(_independent_rows()[0], artifact)
        self.assertEqual(
            prediction["ml_model_id"],
            artifact["trainable_decision"]["model_id"],
        )
        self.assertIn(
            prediction["ml_diagnostic_status"],
            {
                "supported_same_compound",
                "conflicting_spectra",
                "insufficient_evidence",
            },
        )

    def test_locked_final_batch_name_is_normalized_once_for_deployment(self):
        artifact = train_shadow_model(
            _independent_rows(),
            final_batch="  BATCH-FINAL  ",
            source_sidecar=_batch_source_manifest(),
            diagnostic_standard=_standard(),
            bootstrap_iterations=5,
        )
        decision = artifact["trainable_decision"]
        self.assertEqual(decision["provenance"]["locked_final_batch"], "BATCH-FINAL")
        self.assertEqual(decision["locked_final_evaluation"]["batch"], "BATCH-FINAL")
        prediction = predict_shadow(_independent_rows()[0], artifact)
        self.assertEqual(prediction["ml_model_id"], decision["model_id"])

    def test_manifest_records_real_paths_hashes_seed_and_groups(self):
        rows = _independent_rows()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "truth.csv"
            with input_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            sidecar_path = root / "truth.csv.metadata.json"
            source_manifest = _batch_source_manifest()
            sidecar_path.write_text(json.dumps(source_manifest), encoding="utf-8")
            standard_path = Path("MSAI/standards/ms2_diagnostic_standard_v2.json")
            artifact = train_shadow_model(
                rows,
                final_batch="BATCH-FINAL",
                source_sidecar=source_manifest,
                diagnostic_standard=_standard(),
                input_path=input_path,
                source_sidecar_path=sidecar_path,
                diagnostic_standard_path=standard_path,
                seed=1234,
                bootstrap_iterations=10,
            )
            provenance = artifact["trainable_decision"]["provenance"]
            self.assertEqual(provenance["seed"], 1234)
            self.assertEqual(provenance["training_input"]["path"], str(input_path.resolve()))
            for source in ("training_input", "source_sidecar", "diagnostic_standard"):
                self.assertEqual(len(provenance[source]["sha256"]), 64)
            self.assertEqual(
                set(provenance["batch_sources"]),
                {
                    "BATCH-A",
                    "BATCH-B",
                    "BATCH-C",
                    "BATCH-FINAL",
                },
            )
            for source in provenance["batch_sources"].values():
                self.assertTrue(source["path"].startswith("<embedded:"))
                self.assertEqual(len(source["sha256"]), 64)
                self.assertEqual(source["proof_type"], "embedded_canonical_json")

    def test_provenance_paths_must_match_the_exact_training_inputs(self):
        rows = _independent_rows()
        source_manifest = _batch_source_manifest()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "truth.csv"
            with input_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            sidecar_path = root / "source.json"
            sidecar_path.write_text(json.dumps(source_manifest), encoding="utf-8")
            standard_path = root / "standard.json"
            standard_path.write_text(json.dumps(_standard()), encoding="utf-8")

            changed_rows = [dict(row) for row in rows]
            changed_rows[0]["truth_label"] = "negative_different_or_interference"
            with self.assertRaisesRegex(ValueError, "does not match.*table data"):
                train_shadow_model(
                    changed_rows,
                    final_batch="BATCH-FINAL",
                    source_sidecar=source_manifest,
                    diagnostic_standard=_standard(),
                    input_path=input_path,
                    source_sidecar_path=sidecar_path,
                    diagnostic_standard_path=standard_path,
                    bootstrap_iterations=5,
                )

            changed_manifest = _batch_source_manifest()
            changed_manifest["unexpected"] = True
            with self.assertRaisesRegex(ValueError, "does not match.*json data"):
                train_shadow_model(
                    rows,
                    final_batch="BATCH-FINAL",
                    source_sidecar=changed_manifest,
                    diagnostic_standard=_standard(),
                    input_path=input_path,
                    source_sidecar_path=sidecar_path,
                    diagnostic_standard_path=standard_path,
                    bootstrap_iterations=5,
                )

            missing_path = root / "missing.json"
            with self.assertRaisesRegex(ValueError, "does not exist or is not a file"):
                train_shadow_model(
                    rows,
                    final_batch="BATCH-FINAL",
                    source_sidecar=source_manifest,
                    diagnostic_standard=_standard(),
                    input_path=input_path,
                    source_sidecar_path=sidecar_path,
                    diagnostic_standard_path=missing_path,
                    bootstrap_iterations=5,
                )

    def test_multi_batch_training_requires_complete_per_batch_source_proof(self):
        rows = _independent_rows()
        with self.assertRaisesRegex(ValueError, "Multi-batch training requires"):
            train_shadow_model(
                rows,
                final_batch="BATCH-FINAL",
                source_sidecar=_sidecar(),
                diagnostic_standard=_standard(),
                bootstrap_iterations=5,
            )
        incomplete = _batch_source_manifest(("BATCH-A", "BATCH-B", "BATCH-C"))
        with self.assertRaisesRegex(ValueError, "missing eligible batches: BATCH-FINAL"):
            train_shadow_model(
                rows,
                final_batch="BATCH-FINAL",
                source_sidecar=incomplete,
                diagnostic_standard=_standard(),
                bootstrap_iterations=5,
            )
        extra = _batch_source_manifest((*incomplete["batch_sources"], "BATCH-FINAL", "EXTRA"))
        with self.assertRaisesRegex(ValueError, "unexpected batches: EXTRA"):
            train_shadow_model(
                rows,
                final_batch="BATCH-FINAL",
                source_sidecar=extra,
                diagnostic_standard=_standard(),
                bootstrap_iterations=5,
            )

    def test_batch_source_paths_are_relative_to_manifest_and_hash_checked(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_directory = root / "batch-sidecars"
            source_directory.mkdir()
            entries = {}
            for batch in ("BATCH-A", "BATCH-B", "BATCH-C", "BATCH-FINAL"):
                source_path = source_directory / f"{batch}.json"
                source_path.write_text(json.dumps(_sidecar()), encoding="utf-8")
                entries[batch] = {
                    "path": str(source_path.relative_to(root)),
                    "sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
                }
            manifest = {"batch_sources": entries}
            manifest_path = root / "batch-source-manifest.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            artifact = train_shadow_model(
                _independent_rows(),
                final_batch="BATCH-FINAL",
                source_sidecar=manifest,
                source_sidecar_path=manifest_path,
                diagnostic_standard=_standard(),
                bootstrap_iterations=5,
            )
            recorded = artifact["trainable_decision"]["provenance"]["batch_sources"]
            for batch, source in recorded.items():
                self.assertEqual(
                    source["path"], str((source_directory / f"{batch}.json").resolve())
                )
                self.assertEqual(source["sha256"], entries[batch]["sha256"])
                self.assertEqual(source["proof_type"], "source_file")

            missing = _batch_source_manifest()
            missing["batch_sources"]["BATCH-A"] = {"path": "does-not-exist.json"}
            manifest_path.write_text(json.dumps(missing), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "does not exist or is not a file"):
                train_shadow_model(
                    _independent_rows(),
                    final_batch="BATCH-FINAL",
                    source_sidecar=missing,
                    source_sidecar_path=manifest_path,
                    diagnostic_standard=_standard(),
                    bootstrap_iterations=5,
                )

    def test_specificity_target_cannot_relax_below_95_percent(self):
        with self.assertRaisesRegex(ValueError, r"specificity_target must be in \[0.95, 1\]"):
            train_shadow_model(
                _independent_rows(),
                final_batch="BATCH-FINAL",
                source_sidecar=_batch_source_manifest(),
                diagnostic_standard=_standard(),
                specificity_target=0.949,
                bootstrap_iterations=5,
            )

    def test_deployment_training_spectrum_columns_are_frozen(self):
        with self.assertRaisesRegex(ValueError, "spectrum columns are frozen"):
            train_shadow_model(
                _independent_rows(),
                final_batch="BATCH-FINAL",
                source_sidecar=_batch_source_manifest(),
                diagnostic_standard=_standard(),
                spectrum_a_column="alternate_a",
                bootstrap_iterations=5,
            )

    def test_unlabelled_feasibility_data_fails_before_artifact_write(self):
        rows = [
            {
                "target_uid": "unlabelled",
                "Compound_ID": "X",
                "batch_id": "ONLY-BATCH",
                "truth_label": "",
                "peak_a_MS2": POSITIVE_A,
                "peak_b_MS2": POSITIVE_B,
            }
        ]
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "must-not-exist.json"
            with self.assertRaisesRegex(ValueError, "No independent truth"):
                train_and_write_shadow_model(
                    rows,
                    output,
                    final_batch="ONLY-BATCH",
                    source_sidecar=_sidecar(),
                    diagnostic_standard=_standard(),
                    bootstrap_iterations=5,
                )
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
