import csv
import json
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path

from MSAI.python.ms2_ml.contracts import (
    FROZEN_SOURCE_PARAMETERS,
    ML_OUTPUT_FIELDS,
    canonical_model_id,
)
from MSAI.python.ms2_ml.features import FEATURE_NAMES
from MSAI.python.ms2_ml.review_package import package_shadow_review
from MSAI.python.ms2_ml.shadow import (
    apply_shadow_csv,
    apply_shadow_rows,
    predict_shadow,
)
from MSAI.python.ms2_ml.shadow_review import build_shadow_review

SUPPORTED_A = "100,1000;120,900;140,800;160,700;180,600;200,500"
SUPPORTED_B = "100.002,950;120.002,850;140.002,750;160.002,650;180.002,550;200.002,450"
CONFLICT_A = "100,100;120,100;140,100;300,10000"
CONFLICT_B = "100.002,100;120.002,100;140.002,100;400,10000"
SPARSE_A = "100,1000;120,900;140,800"
SPARSE_B = "100.002,1000;221,900;241,800"


def _threshold_artifact():
    return _complete_artifact(
        {
            "schema_version": "ms2_shadow_trainable_decision_v1",
            "selected_model": "threshold",
            "feature_names": list(FEATURE_NAMES),
            "model": {
                "type": "threshold",
                "parameters": {
                    "min_cosine": 0.7,
                    "min_matched_peaks": 3,
                    "min_explained_intensity": 0.5,
                    "min_entropy_similarity": None,
                },
                "bucket_probabilities": {
                    "supported_same_compound": 0.9,
                    "conflicting_spectra": 0.1,
                    "insufficient_evidence": 0.4,
                },
            },
            "decision_thresholds": {
                "conflict_max_probability": 0.0,
                "support_min_probability": 1.0,
            },
            "training_feature_ranges": {},
        }
    )


def _logistic_artifact():
    return _complete_artifact(
        {
            "schema_version": "ms2_shadow_trainable_decision_v1",
            "selected_model": "logistic_regression",
            "feature_names": list(FEATURE_NAMES),
            "model": {
                "type": "logistic_regression",
                "model_type": "standard_scaler_l2_logistic_regression",
                "feature_names": list(FEATURE_NAMES),
                "scaler": {
                    "feature_names": list(FEATURE_NAMES),
                    "means": [0.0] * len(FEATURE_NAMES),
                    "scales": [1.0] * len(FEATURE_NAMES),
                },
                "coefficients": [10.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                "intercept": -5.0,
                "C": 1.0,
                "class_weight": None,
                "optimizer": {
                    "algorithm": "deterministic_newton_irls",
                    "iterations": 1,
                    "converged": True,
                    "objective": 0.0,
                },
            },
            "decision_thresholds": {
                "conflict_max_probability": 0.2,
                "support_min_probability": 0.8,
            },
            "training_feature_ranges": {},
        }
    )


def _complete_artifact(payload):
    selected = payload["selected_model"]
    batches = ("A", "B", "C")
    batch_rows = {batch: [f"{batch}-P", f"{batch}-N"] for batch in (*batches, "FINAL")}
    batch_compounds = {batch: [f"{batch}-P", f"{batch}-N"] for batch in (*batches, "FINAL")}

    def fold_manifest(train_batches, validation_batch, name):
        return {
            "name": name,
            "train_row_ids": [row for batch in train_batches for row in batch_rows[batch]],
            "validation_row_ids": list(batch_rows[validation_batch]),
            "purged_row_ids": [],
            "train_batches": list(train_batches),
            "validation_batches": [validation_batch],
            "train_compounds": [
                compound for batch in train_batches for compound in batch_compounds[batch]
            ],
            "validation_compounds": list(batch_compounds[validation_batch]),
        }

    def threshold_tuning():
        return {
            "parameters": {
                "min_cosine": 0.7,
                "min_matched_peaks": 3,
                "min_explained_intensity": 0.5,
                "min_entropy_similarity": None,
            },
            "metrics": evaluation(
                positives_supported=1,
                positive_count=1,
                negative_count=1,
            )["metrics"],
            "candidates_evaluated": 6384,
            "specificity_target": 0.95,
            "constraint_met": True,
            "bucket_probabilities": {
                "supported_same_compound": 0.9,
                "conflicting_spectra": 0.1,
                "insufficient_evidence": 0.4,
            },
        }

    def logistic_tuning():
        return {
            "c": 1.0,
            "class_weight": None,
            "support_threshold": 0.8,
            "conflict_threshold": 0.2,
            "metrics": evaluation(
                positives_supported=1,
                positive_count=1,
                negative_count=1,
            )["metrics"],
            "candidates_evaluated": 14,
            "candidates_attempted": 14,
            "convergence_failures": [],
        }

    def evaluation(*, positives_supported, positive_count, negative_count):
        tp = positives_supported
        fn = positive_count - positives_supported
        tn = negative_count
        fp = 0
        precision = 1.0 if tp else 0.0
        recall = tp / positive_count
        specificity = 1.0
        f0_5 = 0.0 if not tp else 1.25 * precision * recall / (0.25 * precision + recall)
        row_count = positive_count + negative_count
        return {
            "metrics": {
                "threshold": 0.5,
                "evaluable_rows": row_count,
                "tp": tp,
                "fp": fp,
                "tn": tn,
                "fn": fn,
                "precision": precision,
                "recall": recall,
                "specificity": specificity,
                "pr_auc": 1.0,
                "average_precision": 1.0,
                "f0_5": f0_5,
                "brier_score": 0.0,
                "expected_calibration_error": 0.0,
            },
            "row_count": row_count,
            "compound_count": row_count,
        }

    outer_manifests = []
    outer_evaluations = []
    for validation_batch in batches:
        train_batches = tuple(batch for batch in batches if batch != validation_batch)
        outer = fold_manifest(train_batches, validation_batch, f"outer:{validation_batch}")
        outer["inner_folds"] = [
            fold_manifest(
                tuple(batch for batch in train_batches if batch != inner_validation),
                inner_validation,
                f"inner_outer_{validation_batch}:{inner_validation}",
            )
            for inner_validation in train_batches
        ]
        outer_manifests.append(outer)
        outer_evaluations.append(
            {
                "fold": outer["name"],
                "validation_batch": validation_batch,
                "validation_compounds": list(batch_compounds[validation_batch]),
                "threshold": {
                    "tuning": threshold_tuning(),
                    "evaluation": evaluation(
                        positives_supported=0 if selected == "logistic_regression" else 1,
                        positive_count=1,
                        negative_count=1,
                    ),
                },
                "logistic_regression": {
                    "tuning": logistic_tuning(),
                    "evaluation": evaluation(
                        positives_supported=1,
                        positive_count=1,
                        negative_count=1,
                    ),
                },
            }
        )

    reliable = selected == "logistic_regression"
    batch_deltas = [
        {
            "batch": batch,
            "recall_delta": 1.0 if reliable else 0.0,
            "threshold_recall": 0.0 if reliable else 1.0,
            "logistic_recall": 1.0,
            "logistic_specificity": 1.0,
        }
        for batch in batches
    ]
    payload["training_feature_ranges"] = {
        name: {"minimum": -10.0, "maximum": 10.0} for name in FEATURE_NAMES
    }
    payload.update(
        {
            "candidates": {
                "threshold": {
                    "grid_size": 6384,
                    "final_development_tuning": threshold_tuning(),
                },
                "logistic_regression": {
                    "c_grid": [0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0],
                    "class_weight_grid": [None, "balanced"],
                    "final_development_tuning": logistic_tuning(),
                },
            },
            "selection": {
                "selected_model": selected,
                "threshold_wins_ties": True,
                "criterion": "synthetic fixture follows the fixed every-batch rule",
                "reliably_better": reliable,
                "specificity_target": 0.95,
                "batch_deltas": batch_deltas,
            },
            "development_validation": {
                "strategy": "nested_leave_one_batch_out_with_compound_purge",
                "specificity_target": 0.95,
                "outer_folds": outer_evaluations,
                "pooled": {
                    "threshold": evaluation(
                        positives_supported=0 if reliable else 3,
                        positive_count=3,
                        negative_count=3,
                    ),
                    "logistic_regression": evaluation(
                        positives_supported=3,
                        positive_count=3,
                        negative_count=3,
                    ),
                },
            },
            "locked_final_evaluation": {
                "batch": "FINAL",
                "used_once_after_model_selection": True,
                "evaluation": evaluation(
                    positives_supported=1,
                    positive_count=1,
                    negative_count=1,
                ),
            },
            "split_manifest": {
                "strategy": "locked_batch_then_nested_leave_one_batch_out",
                "eligible_batches": ["A", "B", "C", "FINAL"],
                "development_row_ids": [row for batch in batches for row in batch_rows[batch]],
                "locked_final_row_ids": list(batch_rows["FINAL"]),
                "locked_compound_overlap_purged_row_ids": [],
                "locked_compound_overlap_purged_batches": [],
                "development_batches": list(batches),
                "locked_final_batches": ["FINAL"],
                "development_compounds": [
                    compound for batch in batches for compound in batch_compounds[batch]
                ],
                "locked_final_compounds": ["FINAL-P", "FINAL-N"],
                "outer_folds": outer_manifests,
                "final_development_folds": [
                    fold_manifest(
                        tuple(batch for batch in batches if batch != validation_batch),
                        validation_batch,
                        f"final_development:{validation_batch}",
                    )
                    for validation_batch in batches
                ],
            },
            "annotation_summary": {
                "input_rows": 8,
                "eligible_rows": 8,
                "positive_same_compound": 4,
                "negative_different_or_interference": 4,
                "uncertain": 0,
                "not_evaluable": 0,
            },
            "provenance": {
                "truth_contract": "independent_external_labels_required",
                "frozen_source_contract_validated": True,
                "batch_source_contract_validated": True,
                "seed": 7,
                "label_column": "truth_label",
                "compound_group_column": "Compound_ID",
                "batch_group_column": "batch_id",
                "row_id_column": "target_uid",
                "locked_final_batch": "FINAL",
                "training_input": {"path": "synthetic.csv", "sha256": "0" * 64},
                "source_sidecar": {"path": "synthetic.json", "sha256": "1" * 64},
                "diagnostic_standard": {"path": "v2.json", "sha256": "2" * 64},
                "batch_sources": {
                    batch: {
                        "path": f"<embedded:{batch}>",
                        "sha256": str(index) * 64,
                        "proof_type": "embedded_canonical_json",
                    }
                    for index, batch in enumerate(("A", "B", "C", "FINAL"), start=3)
                },
            },
            "model_card": {"intended_use": "structurally complete inference fixture"},
            "limitations": ["Synthetic fixture; not an approved release artifact."],
        }
    )
    payload["model_id"] = canonical_model_id(payload)
    return {"trainable_decision": payload}


def _review_row(compound="example"):
    return {
        "Compound_ID": compound,
        "SGC ID for Pool": "POOL-A01",
        "Pooled Well": "A1",
        "MZ": "400.1",
        "IG": "GREEN",
        "Peak1": "1.0",
        "Peak2": "2.0",
        "peak_a_MS2": SUPPORTED_A,
        "peak_b_MS2": SUPPORTED_B,
        "peak_a_quality_flags": "ok",
        "peak_b_quality_flags": "ok",
        "peak_a_rt_window_start": "50",
        "peak_a_rt_window_end": "70",
        "peak_b_rt_window_start": "110",
        "peak_b_rt_window_end": "130",
        "peak_a_apex_rt": "60",
        "peak_b_apex_rt": "120",
        "peak_a_apex_intensity": "100000",
        "peak_b_apex_intensity": "90000",
        "peak_a_area": "500000",
        "peak_b_area": "450000",
        "peak_a_ms2_count": "4",
        "peak_b_ms2_count": "4",
        "peak_a_consensus_scan_count": "1",
        "peak_b_consensus_scan_count": "1",
        "peak_a_candidate_fragment_count": "6",
        "peak_b_candidate_fragment_count": "6",
        "peak_a_fragment_count": "6",
        "peak_b_fragment_count": "6",
        "dia_window_center": "403",
        "dia_window_lower": "395.5",
        "dia_window_upper": "410.5",
        "dia_window_match_count": "1",
        "rt_separation_sec": "60",
        "ms2_cosine": "0.99",
        "ms2_entropy_similarity": "0.99",
        "ms2_matched_peaks": "6",
        "peak_a_explained_intensity": "1",
        "peak_b_explained_intensity": "1",
        "ms1_reference_status": "supplied_double_peak",
        "ms1_reference_issue_codes": "",
        "ms2_diagnostic_status": "supported_same_compound",
        "ms2_issue_codes": "NEAR_DECISION_BOUNDARY",
        "ms2_review_priority": "P2_boundary_or_interference",
        "enantiomer_pair_status": "candidate_enantiomer_pair",
        "enantiomer_pair_reason": "screen_pass",
    }


def _sidecar():
    parameters: dict[str, object] = dict(FROZEN_SOURCE_PARAMETERS)
    parameters.update(
        {
            "min_cosine": 0.7,
            "min_matched_peaks": 6,
            "min_explained_intensity": 0.5,
            "min_entropy_similarity": None,
        }
    )
    return {
        "parameters": parameters,
        "inputs": {},
        "dia_windows": [{"center": 403, "lower": 395.5, "upper": 410.5, "scan_count": 10}],
    }


def _write_csv(path, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _rehash(artifact):
    payload = artifact["trainable_decision"]
    payload["model_id"] = canonical_model_id(payload)
    return artifact


class ShadowPredictionTests(unittest.TestCase):
    def test_untrained_feasibility_never_emits_a_probability_or_identity_call(self):
        prediction = predict_shadow({"peak_a_MS2": SUPPORTED_A, "peak_b_MS2": SUPPORTED_B})
        self.assertEqual(tuple(prediction), ML_OUTPUT_FIELDS)
        self.assertEqual(prediction["ml_same_compound_probability"], "")
        self.assertEqual(prediction["ml_diagnostic_status"], "not_evaluable")
        self.assertEqual(prediction["ml_model_id"], "UNTRAINED")
        self.assertEqual(prediction["ml_abstention_reason"], "NO_INDEPENDENT_TRUTH_MODEL")
        self.assertIn("FEASIBILITY_ONLY", prediction["ml_review_flags"])

    def test_threshold_model_emits_supported_conflicting_insufficient_and_missing(self):
        cases = (
            (SUPPORTED_A, SUPPORTED_B, "supported_same_compound", "0.9"),
            (CONFLICT_A, CONFLICT_B, "conflicting_spectra", "0.1"),
            (SPARSE_A, SPARSE_B, "insufficient_evidence", "0.4"),
            (SUPPORTED_A, "", "not_evaluable", ""),
        )
        for spectrum_a, spectrum_b, status, probability in cases:
            with self.subTest(status=status):
                result = predict_shadow(
                    {"peak_a_MS2": spectrum_a, "peak_b_MS2": spectrum_b},
                    _threshold_artifact(),
                )
                self.assertEqual(result["ml_diagnostic_status"], status)
                self.assertEqual(result["ml_same_compound_probability"], probability)

    def test_logistic_high_low_and_middle_probability_abstention(self):
        middle_a = "100,100;200,100"
        middle_b = "100.002,100;300,100"
        cases = (
            (SUPPORTED_A, SUPPORTED_B, "supported_same_compound"),
            (CONFLICT_A, "101,100;121,100;141,100;401,10000", "conflicting_spectra"),
            (middle_a, middle_b, "insufficient_evidence"),
        )
        for spectrum_a, spectrum_b, status in cases:
            with self.subTest(status=status):
                result = predict_shadow(
                    {"peak_a_MS2": spectrum_a, "peak_b_MS2": spectrum_b},
                    _logistic_artifact(),
                )
                self.assertEqual(result["ml_diagnostic_status"], status)
                if status == "insufficient_evidence":
                    self.assertEqual(
                        result["ml_abstention_reason"], "PROBABILITY_ABSTENTION_INTERVAL"
                    )

    def test_trained_inference_rejects_incomplete_forged_or_off_grid_artifacts(self):
        row = {"peak_a_MS2": SUPPORTED_A, "peak_b_MS2": SUPPORTED_B}
        minimal = {
            "trainable_decision": {
                "schema_version": "ms2_shadow_trainable_decision_v1",
                "model_id": "forged",
                "selected_model": "threshold",
                "feature_names": list(FEATURE_NAMES),
                "model": _threshold_artifact()["trainable_decision"]["model"],
                "decision_thresholds": {
                    "conflict_max_probability": 0.0,
                    "support_min_probability": 1.0,
                },
            }
        }
        with self.assertRaisesRegex(ValueError, "schema mismatch"):
            predict_shadow(row, minimal)

        forged_id = _threshold_artifact()
        forged_id["trainable_decision"]["model_id"] = "ms2-shadow-forged"
        with self.assertRaisesRegex(ValueError, "canonical.*hash"):
            predict_shadow(row, forged_id)

        wrong_schema = _threshold_artifact()
        wrong_schema["trainable_decision"]["schema_version"] = "v2"
        with self.assertRaisesRegex(ValueError, "schema_version"):
            predict_shadow(row, wrong_schema)

        off_grid = _threshold_artifact()
        off_grid["trainable_decision"]["model"]["parameters"]["min_cosine"] = 0.71
        with self.assertRaisesRegex(ValueError, "fixed 6,384-point grid"):
            predict_shadow(row, off_grid)

        mismatched_threshold = _threshold_artifact()
        mismatched_threshold["trainable_decision"]["model"]["parameters"]["min_cosine"] = 0.95
        with self.assertRaisesRegex(ValueError, "parameters do not match final tuning"):
            predict_shadow(row, _rehash(mismatched_threshold))

        mismatched_buckets = _threshold_artifact()
        mismatched_buckets["trainable_decision"]["model"]["bucket_probabilities"][
            "supported_same_compound"
        ] = 0.8
        with self.assertRaisesRegex(ValueError, "bucket probabilities do not match"):
            predict_shadow(row, _rehash(mismatched_buckets))

        mismatched_logistic_c = _logistic_artifact()
        mismatched_logistic_c["trainable_decision"]["model"]["C"] = 10.0
        with self.assertRaisesRegex(ValueError, "Logistic configuration does not match"):
            predict_shadow(row, _rehash(mismatched_logistic_c))

        mismatched_logistic_boundary = _logistic_artifact()
        mismatched_logistic_boundary["trainable_decision"]["decision_thresholds"][
            "support_min_probability"
        ] = 0.85
        with self.assertRaisesRegex(ValueError, "probability boundaries do not match"):
            predict_shadow(row, _rehash(mismatched_logistic_boundary))

        forged_selection = _threshold_artifact()
        forged_delta = forged_selection["trainable_decision"]["selection"]["batch_deltas"][0]
        forged_delta.update({"threshold_recall": 0.5, "logistic_recall": 0.5, "recall_delta": 0.0})
        with self.assertRaisesRegex(ValueError, "do not match outer-fold evaluations"):
            predict_shadow(row, _rehash(forged_selection))

        forged_locked_count = _threshold_artifact()
        forged_locked_count["trainable_decision"]["locked_final_evaluation"]["evaluation"][
            "row_count"
        ] = 3
        with self.assertRaisesRegex(ValueError, "row_count does not match split evidence"):
            predict_shadow(row, _rehash(forged_locked_count))

        missing_locked_metric = _threshold_artifact()
        missing_locked_metric["trainable_decision"]["locked_final_evaluation"]["evaluation"][
            "metrics"
        ].pop("recall")
        with self.assertRaisesRegex(ValueError, "exact all-row metric schema"):
            predict_shadow(row, _rehash(missing_locked_metric))

        missing_tuning_metric = _threshold_artifact()
        missing_tuning_metric["trainable_decision"]["candidates"]["threshold"][
            "final_development_tuning"
        ]["metrics"].pop("recall")
        with self.assertRaisesRegex(ValueError, "exact all-row metric schema"):
            predict_shadow(row, _rehash(missing_tuning_metric))

        missing_batch_proof = _threshold_artifact()
        missing_batch_proof["trainable_decision"]["provenance"].pop("batch_sources")
        with self.assertRaisesRegex(ValueError, "release-evidence schema|batch_sources"):
            predict_shadow(row, _rehash(missing_batch_proof))

        relaxed_specificity = _threshold_artifact()
        relaxed_specificity["trainable_decision"]["selection"]["specificity_target"] = 0.949
        relaxed_specificity["trainable_decision"]["development_validation"][
            "specificity_target"
        ] = 0.949
        with self.assertRaisesRegex(ValueError, r"specificity_target must be in \[0.95, 1\]"):
            predict_shadow(row, _rehash(relaxed_specificity))

        forged_inner = _threshold_artifact()
        forged_inner["trainable_decision"]["split_manifest"]["outer_folds"][0]["inner_folds"] = [
            {"name": "forged"}
        ]
        with self.assertRaisesRegex(ValueError, "at least two inner folds|fold schema is invalid"):
            predict_shadow(row, _rehash(forged_inner))

        leaking_inner = _threshold_artifact()
        leaking_fold = leaking_inner["trainable_decision"]["split_manifest"]["outer_folds"][0][
            "inner_folds"
        ][0]
        leaking_fold["validation_row_ids"] = ["A-P", "A-N"]
        leaking_fold["validation_batches"] = ["A"]
        leaking_fold["validation_compounds"] = ["A-P", "A-N"]
        with self.assertRaisesRegex(ValueError, "partition the parent universe|outside its parent"):
            predict_shadow(row, _rehash(leaking_inner))

        for bad_convergence in (False, "false"):
            forged_logistic = _logistic_artifact()
            forged_logistic["trainable_decision"]["model"]["optimizer"]["converged"] = (
                bad_convergence
            )
            with (
                self.subTest(converged=bad_convergence),
                self.assertRaisesRegex(ValueError, "converged=true"),
            ):
                predict_shadow(row, _rehash(forged_logistic))
        for bad_iterations in (0, 201):
            forged_logistic = _logistic_artifact()
            forged_logistic["trainable_decision"]["model"]["optimizer"]["iterations"] = (
                bad_iterations
            )
            with (
                self.subTest(iterations=bad_iterations),
                self.assertRaisesRegex(ValueError, r"iterations must be an integer in \[1, 200\]"),
            ):
                predict_shadow(row, _rehash(forged_logistic))

    def test_apply_rows_preserves_source_and_appends_exactly_five_fields(self):
        source = _review_row()
        original = deepcopy(source)
        combined = apply_shadow_rows([source], _threshold_artifact())[0]
        self.assertEqual(source, original)
        self.assertEqual(
            {key: combined[key] for key in original},
            original,
        )
        self.assertEqual(tuple(combined)[-5:], ML_OUTPUT_FIELDS)
        with self.assertRaisesRegex(ValueError, "refusing to overwrite"):
            apply_shadow_rows([combined], _threshold_artifact())


class ShadowCsvAndReviewTests(unittest.TestCase):
    def test_csv_sidecar_and_every_disagreement_png_svg_are_written(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "v2.csv"
            disagreement = _review_row("disagreement")
            agreement = _review_row("agreement")
            # The same spectra produce ML support. Only the first row disagrees.
            disagreement["ms2_diagnostic_status"] = "conflicting_spectra"
            disagreement["ms2_issue_codes"] = "LOW_COSINE_SIMILARITY"
            disagreement["ms2_review_priority"] = "P1_conflict"
            _write_csv(source, [disagreement, agreement])
            source_sidecar = Path(str(source) + ".metadata.json")
            sidecar_payload = _sidecar()
            method_profile_path = Path(
                "MSAI/config/acquisition_profiles/adductmlib_chemrxiv_20260129_v1.json"
            )
            sidecar_payload["provenance_layers"] = {
                "reference_method": json.loads(method_profile_path.read_text(encoding="utf-8"))
            }
            source_sidecar.write_text(json.dumps(sidecar_payload), encoding="utf-8")
            model = root / "model.json"
            model.write_text(json.dumps(_threshold_artifact()), encoding="utf-8")
            shadow = root / "shadow.csv"
            summary = apply_shadow_csv(
                source,
                shadow,
                model_artifact_path=model,
                source_sidecar_path=source_sidecar,
            )
            self.assertEqual(summary["rows"], 2)
            with shadow.open(encoding="utf-8-sig", newline="") as handle:
                shadow_rows = list(csv.DictReader(handle))
            self.assertEqual(list(shadow_rows[0])[-5:], list(ML_OUTPUT_FIELDS))
            self.assertEqual(shadow_rows[0]["ms2_issue_codes"], "LOW_COSINE_SIMILARITY")
            combined_sidecar = json.loads(
                Path(str(shadow) + ".metadata.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                combined_sidecar["shadow_inference"]["model"]["model_id"],
                _threshold_artifact()["trainable_decision"]["model_id"],
            )
            self.assertEqual(
                len(combined_sidecar["shadow_inference"]["source_result_content_sha256"]),
                64,
            )
            self.assertEqual(
                len(combined_sidecar["shadow_inference"]["source_sidecar_content_sha256"]),
                64,
            )

            output = root / "review"
            review = build_shadow_review(
                shadow,
                output,
                model_artifact_path=model,
                standard_path="MSAI/standards/ms2_diagnostic_standard_v2.json",
                method_profile_path=method_profile_path,
                image_formats=("svg", "png"),
                png_scale=0.25,
            )
            self.assertEqual(review["disagreements"], 1)
            self.assertEqual(review["outputs"]["svg_images"], 1)
            self.assertEqual(review["outputs"]["png_images"], 1)
            self.assertEqual(len(list((output / "assets").rglob("*.svg"))), 1)
            self.assertEqual(len(list((output / "assets").rglob("*.png"))), 1)
            with (output / "v2_ml_disagreements.csv").open(
                encoding="utf-8-sig", newline=""
            ) as handle:
                disagreements = list(csv.DictReader(handle))
            self.assertEqual(disagreements[0]["ms2_diagnostic_status"], "conflicting_spectra")
            self.assertEqual(disagreements[0]["ms2_issue_codes"], "LOW_COSINE_SIMILARITY")
            self.assertEqual(disagreements[0]["ms2_review_priority"], "P1_conflict")
            with (output / "metadata" / "target_manifest.csv").open(
                encoding="utf-8-sig", newline=""
            ) as handle:
                manifest = list(csv.DictReader(handle))
            for field in ML_OUTPUT_FIELDS:
                self.assertIn(field, manifest[0])
            self.assertTrue((output / manifest[0]["svg_asset_path"]).is_file())
            self.assertTrue((output / manifest[0]["png_asset_path"]).is_file())
            html = (output / "v2_ml_disagreements.html").read_text(encoding="utf-8")
            self.assertIn("conflicting_spectra", html)
            self.assertIn("supported_same_compound", html)
            self.assertIn(_threshold_artifact()["trainable_decision"]["model_id"], html)
            packaged = package_shadow_review(
                shadow,
                output,
                root / "review.zip",
                model_artifact_path=model,
            )
            self.assertEqual(packaged["review_rows"], 1)

    def test_shadow_review_rejects_invalid_model_before_writing_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "v2.csv"
            row = _review_row("disagreement")
            row["ms2_diagnostic_status"] = "conflicting_spectra"
            _write_csv(source, [row])
            sidecar = Path(str(source) + ".metadata.json")
            sidecar.write_text(json.dumps(_sidecar()), encoding="utf-8")
            valid_model = root / "valid-model.json"
            valid_model.write_text(json.dumps(_threshold_artifact()), encoding="utf-8")
            shadow = root / "shadow.csv"
            apply_shadow_csv(
                source,
                shadow,
                model_artifact_path=valid_model,
                source_sidecar_path=sidecar,
            )
            invalid = _threshold_artifact()
            invalid["trainable_decision"]["candidates"].pop("logistic_regression")
            invalid_model = root / "invalid-model.json"
            invalid_model.write_text(json.dumps(invalid), encoding="utf-8")
            output = root / "must-not-exist"
            with self.assertRaisesRegex(ValueError, "both permitted model classes"):
                build_shadow_review(
                    shadow,
                    output,
                    model_artifact_path=invalid_model,
                    standard_path="MSAI/standards/ms2_diagnostic_standard_v2.json",
                )
            self.assertFalse(output.exists())

    def test_shadow_review_binds_every_row_to_one_declared_model_prediction(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            artifact = _threshold_artifact()
            model = root / "model.json"
            model.write_text(json.dumps(artifact), encoding="utf-8")

            mixed_rows = apply_shadow_rows(
                [_review_row("first"), _review_row("second")],
                artifact,
            )
            mixed_rows[1]["ml_model_id"] = "ms2-shadow-different-model"
            mixed = root / "mixed.csv"
            _write_csv(mixed, mixed_rows)
            mixed_output = root / "mixed-review"
            with self.assertRaisesRegex(ValueError, "exactly one model ID"):
                build_shadow_review(mixed, mixed_output)
            self.assertFalse(mixed_output.exists())

            forged_rows = apply_shadow_rows([_review_row("forged")], artifact)
            forged_rows[0]["ml_same_compound_probability"] = "0.1"
            forged = root / "forged.csv"
            _write_csv(forged, forged_rows)
            forged_output = root / "forged-review"
            with self.assertRaisesRegex(ValueError, "do not match.*model prediction"):
                build_shadow_review(
                    forged,
                    forged_output,
                    model_artifact_path=model,
                )
            self.assertFalse(forged_output.exists())

    def test_apply_csv_validates_sidecar_and_v2_standard_inside_the_api(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "v2.csv"
            _write_csv(source, [_review_row()])

            missing_sidecar_output = root / "missing-sidecar.csv"
            with self.assertRaisesRegex(ValueError, "sidecar does not exist"):
                apply_shadow_csv(source, missing_sidecar_output)
            self.assertFalse(missing_sidecar_output.exists())

            sidecar = root / "v2.csv.metadata.json"
            invalid_sidecar = _sidecar()
            invalid_sidecar["parameters"]["min_correlation_scans"] = 4
            sidecar.write_text(json.dumps(invalid_sidecar), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "min_correlation_scans"):
                apply_shadow_csv(source, root / "invalid-sidecar.csv", source_sidecar_path=sidecar)

            sidecar.write_text(json.dumps(_sidecar()), encoding="utf-8")
            standard = json.loads(
                Path("MSAI/standards/ms2_diagnostic_standard_v2.json").read_text(encoding="utf-8")
            )
            standard["preprocessing"]["fragment_matching"] = "many_to_one"
            standard_path = root / "invalid-standard.json"
            standard_path.write_text(json.dumps(standard), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "fragment_matching"):
                apply_shadow_csv(
                    source,
                    root / "invalid-standard.csv",
                    source_sidecar_path=sidecar,
                    standard_path=standard_path,
                )

    def test_zero_disagreement_still_writes_header_html_and_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "v2.csv"
            _write_csv(source, [_review_row("agreement")])
            source_sidecar = Path(str(source) + ".metadata.json")
            source_sidecar.write_text(json.dumps(_sidecar()), encoding="utf-8")
            model = root / "model.json"
            model.write_text(json.dumps(_threshold_artifact()), encoding="utf-8")
            shadow = root / "shadow.csv"
            apply_shadow_csv(
                source,
                shadow,
                model_artifact_path=model,
                source_sidecar_path=source_sidecar,
            )
            output = root / "review"
            summary = build_shadow_review(
                shadow,
                output,
                model_artifact_path=model,
                image_formats=("svg",),
                png_scale=0.25,
            )
            self.assertEqual(summary["disagreements"], 0)
            csv_path = output / "v2_ml_disagreements.csv"
            self.assertTrue(csv_path.is_file())
            with csv_path.open(encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                self.assertIn("comparison_candidate_status", reader.fieldnames or [])
                self.assertEqual(list(reader), [])
            html = (output / "v2_ml_disagreements.html").read_text(encoding="utf-8")
            self.assertIn("No v2/ML diagnostic disagreements", html)
            self.assertTrue((output / "shadow_summary.json").is_file())
            self.assertFalse(list((output / "assets").rglob("*.svg")))

    def test_untrained_zero_disagreement_and_zero_row_outputs_package_cleanly(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cases = {
                "agreement": [_review_row("agreement")],
                "empty": [],
            }
            cases["agreement"][0]["ms2_diagnostic_status"] = "not_evaluable"
            for name, rows in cases.items():
                with self.subTest(name=name):
                    case_root = root / name
                    case_root.mkdir()
                    source = case_root / "v2.csv"
                    fields = list(_review_row())
                    with source.open("w", encoding="utf-8-sig", newline="") as handle:
                        writer = csv.DictWriter(handle, fieldnames=fields)
                        writer.writeheader()
                        writer.writerows(rows)
                    source_sidecar = Path(str(source) + ".metadata.json")
                    source_sidecar.write_text(json.dumps(_sidecar()), encoding="utf-8")
                    shadow = case_root / "shadow.csv"
                    apply_shadow_csv(
                        source,
                        shadow,
                        source_sidecar_path=source_sidecar,
                    )
                    output = case_root / "review"
                    summary = build_shadow_review(shadow, output)
                    self.assertEqual(summary["disagreements"], 0)
                    run_manifest = json.loads(
                        (output / "metadata/run_manifest.json").read_text(encoding="utf-8")
                    )
                    self.assertEqual(run_manifest["model_artifact"]["model_id"], "UNTRAINED")
                    package = package_shadow_review(
                        shadow,
                        output,
                        case_root / "review.zip",
                    )
                    self.assertEqual(package["review_rows"], 0)


if __name__ == "__main__":
    unittest.main()
