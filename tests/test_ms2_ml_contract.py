import io
import json
import math
import unittest
from contextlib import redirect_stdout
from copy import deepcopy
from pathlib import Path

from MSAI.python.ms2_ml.cli import train_main
from MSAI.python.ms2_ml.contracts import (
    ALLOWED_SOURCE_PARAMETER_KEYS,
    ELIGIBLE_TRUTH_LABELS,
    FROZEN_CLOSE_PEAK_RULE,
    FROZEN_SOURCE_PARAMETERS,
    FROZEN_STANDARD_PREPROCESSING,
    LEAKAGE_DENYLIST,
    ML_DIAGNOSTIC_STATUSES,
    ML_OUTPUT_FIELDS,
    SOURCE_DECISION_GUARDRAIL_KEYS,
    TRUTH_LABELS,
    filter_labeled_rows,
    make_trainable_artifact,
    normalize_truth_label,
    validate_frozen_source,
    validate_no_leakage,
    validate_trainable_artifact,
    validate_truth_row,
)
from MSAI.python.ms2_ml.features import (
    FEATURE_NAMES,
    extract_features_from_row,
)
from MSAI.python.ms2_ml.splits import (
    GroupedExample,
    assert_group_isolation,
    build_leave_one_batch_out_folds,
    build_locked_batch_split,
)

SPECTRUM_A = "100,1000;120,800;140,600;160,400"
SPECTRUM_B = "100.002,900;120.002,700;140.002,500;160.002,300"


def _sidecar(**overrides):
    parameters = dict(FROZEN_SOURCE_PARAMETERS)
    parameters.update(overrides)
    return {"parameters": parameters}


def _example(row_id, compound, batch, label):
    return GroupedExample(
        row_id=row_id,
        compound_id=compound,
        batch_id=batch,
        label=label,
        features=(0.9, 0.9, math.log1p(6), 0.8, 0.0, 1.0),
    )


class TruthAndFeatureContractTests(unittest.TestCase):
    def test_truth_vocabulary_is_exact_and_exclusions_are_explicit(self):
        self.assertEqual(
            TRUTH_LABELS,
            (
                "positive_same_compound",
                "negative_different_or_interference",
                "uncertain",
                "not_evaluable",
            ),
        )
        self.assertEqual(
            ELIGIBLE_TRUTH_LABELS,
            ("positive_same_compound", "negative_different_or_interference"),
        )
        self.assertEqual(normalize_truth_label("positive_same_compound"), TRUTH_LABELS[0])
        for alias in ("positive", "negative", "1", "candidate_enantiomer_pair", ""):
            with self.subTest(alias=alias), self.assertRaises(ValueError):
                normalize_truth_label(alias)

    def test_missing_bilateral_spectrum_is_never_an_eligible_negative(self):
        row = {
            "truth": "negative_different_or_interference",
            "peak_a_MS2": SPECTRUM_A,
            "peak_b_MS2": "",
        }
        with self.assertRaisesRegex(ValueError, "requires bilateral spectra"):
            validate_truth_row(row, "truth")
        row["truth"] = "not_evaluable"
        self.assertEqual(validate_truth_row(row, "truth"), "not_evaluable")

    def test_uncertain_and_not_evaluable_are_excluded_from_fitting(self):
        base = {"peak_a_MS2": SPECTRUM_A, "peak_b_MS2": SPECTRUM_B}
        rows = [
            {**base, "truth": "positive_same_compound"},
            {**base, "truth": "negative_different_or_interference"},
            {**base, "truth": "uncertain"},
            {**base, "truth": "not_evaluable"},
        ]
        self.assertEqual(len(filter_labeled_rows(rows, "truth")), 2)

    def test_features_are_recomputed_from_spectra_and_ignore_rule_outputs(self):
        first = {
            "peak_a_MS2": SPECTRUM_A,
            "peak_b_MS2": SPECTRUM_B,
            "ms2_diagnostic_status": "supported_same_compound",
            "enantiomer_pair_status": "candidate_enantiomer_pair",
            "ms2_issue_codes": "",
            "IG": "GREEN",
        }
        second = dict(first)
        second.update(
            {
                "ms2_diagnostic_status": "conflicting_spectra",
                "enantiomer_pair_status": "ms2_not_similar",
                "ms2_issue_codes": "LOW_COSINE_SIMILARITY",
                "IG": "RED",
            }
        )
        features_a = extract_features_from_row(first)
        features_b = extract_features_from_row(second)
        self.assertIsNotNone(features_a)
        self.assertEqual(features_a, features_b)
        self.assertEqual(tuple(features_a.as_dict()), FEATURE_NAMES)  # type: ignore[union-attr]
        self.assertIsNone(extract_features_from_row({"peak_a_MS2": SPECTRUM_A}))
        for leaked in (
            "ms2_diagnostic_status",
            "enantiomer_pair_status",
            "IG",
            "truth_label",
            "positive_same_compound",
            "ml_same_compound_probability",
            "ml_diagnostic_status",
        ):
            with self.subTest(leaked=leaked), self.assertRaises(ValueError):
                validate_no_leakage((*FEATURE_NAMES, leaked))


class ParameterAndArtifactLockTests(unittest.TestCase):
    def test_frozen_source_accepts_only_rt10_v2_contract(self):
        standard = json.loads(
            Path("MSAI/standards/ms2_diagnostic_standard_v2.json").read_text(encoding="utf-8")
        )
        self.assertEqual(validate_frozen_source(_sidecar(), standard), FROZEN_SOURCE_PARAMETERS)
        for key, bad_value in (
            ("rt_half_window_sec", 8),
            ("min_fragment_correlation", 0.85),
            ("fragment_mz_tol", 0.02),
            ("min_fragment_intensity", 1000),
            ("dia_iso_win_fallback", 12),
            ("correlation_min_relative_intensity", 0.1),
            ("min_correlation_scans", 3),
            ("max_fragment_apex_offset_scans", 2),
            ("min_consecutive_fragment_scans", 2),
            ("consensus_scans", 3),
        ):
            with self.subTest(key=key), self.assertRaises(ValueError):
                validate_frozen_source(_sidecar(**{key: bad_value}), standard)

    def test_explicit_eic_fields_are_all_or_none_before_legacy_fallback(self):
        standard = json.loads(
            Path("MSAI/standards/ms2_diagnostic_standard_v2.json").read_text(encoding="utf-8")
        )
        legacy_parameters = dict(FROZEN_SOURCE_PARAMETERS)
        for key in (
            "precursor_eic_mz_tol",
            "precursor_eic_mz_tol_unit",
            "fragment_eic_mz_tol",
            "fragment_eic_mz_tol_unit",
        ):
            legacy_parameters.pop(key)
        self.assertEqual(
            validate_frozen_source({"parameters": legacy_parameters}, standard),
            FROZEN_SOURCE_PARAMETERS,
        )
        partial = dict(FROZEN_SOURCE_PARAMETERS)
        partial.pop("fragment_eic_mz_tol_unit")
        with self.assertRaisesRegex(ValueError, "all-or-none"):
            validate_frozen_source({"parameters": partial}, standard)

    def test_source_parameters_are_closed_world_and_duplicate_views_must_match(self):
        standard = json.loads(
            Path("MSAI/standards/ms2_diagnostic_standard_v2.json").read_text(encoding="utf-8")
        )
        unknown = _sidecar(future_extraction_filter=999)
        with self.assertRaisesRegex(ValueError, "unknown analysis/extraction parameters"):
            validate_frozen_source(unknown, standard)

        full_parameters: dict[str, object] = dict(FROZEN_SOURCE_PARAMETERS)
        full_parameters.update(
            {
                "min_cosine": 0.7,
                "min_matched_peaks": 6,
                "min_explained_intensity": 0.5,
                "min_entropy_similarity": None,
            }
        )
        self.assertEqual(
            validate_frozen_source({"parameters": full_parameters}, standard),
            FROZEN_SOURCE_PARAMETERS,
        )
        partial_guardrails = dict(full_parameters)
        partial_guardrails.pop("min_entropy_similarity")
        with self.assertRaisesRegex(ValueError, "guardrail fields must be supplied all-or-none"):
            validate_frozen_source({"parameters": partial_guardrails}, standard)
        invalid_guardrails = dict(full_parameters)
        invalid_guardrails["min_matched_peaks"] = 0
        with self.assertRaisesRegex(ValueError, "positive integer"):
            validate_frozen_source({"parameters": invalid_guardrails}, standard)

        conflict = {
            "parameters": dict(FROZEN_SOURCE_PARAMETERS),
            "provenance_layers": {
                "analysis_parameters": {
                    **FROZEN_SOURCE_PARAMETERS,
                    "rt_half_window_sec": 8.0,
                }
            },
        }
        with self.assertRaisesRegex(ValueError, "parameters conflict"):
            validate_frozen_source(conflict, standard)

    def test_every_v2_preprocessing_field_is_frozen(self):
        standard = json.loads(
            Path("MSAI/standards/ms2_diagnostic_standard_v2.json").read_text(encoding="utf-8")
        )
        self.assertEqual(standard["preprocessing"], FROZEN_STANDARD_PREPROCESSING)
        for key, original in FROZEN_STANDARD_PREPROCESSING.items():
            changed = deepcopy(standard)
            changed["preprocessing"][key] = (
                original + 0.1 if isinstance(original, float) else f"{original}-changed"
            )
            with (
                self.subTest(key=key),
                self.assertRaisesRegex(ValueError, "Frozen v2 preprocessing"),
            ):
                validate_frozen_source(_sidecar(), changed)
        changed = deepcopy(standard)
        changed["preprocessing"]["unexpected"] = 1
        with self.assertRaisesRegex(ValueError, "schema mismatch"):
            validate_frozen_source(_sidecar(), changed)

    def test_training_artifact_has_one_section_and_no_method_keys_recursively(self):
        artifact = make_trainable_artifact(
            {"schema_version": "msai-ms2-shadow-model-v1", "model": {"type": "threshold"}}
        )
        self.assertEqual(set(artifact), {"trainable_decision"})
        self.assertEqual(validate_trainable_artifact(artifact)["model"]["type"], "threshold")
        with self.assertRaisesRegex(ValueError, "top level"):
            validate_trainable_artifact({**artifact, "method": {}})
        with self.assertRaisesRegex(ValueError, "non-trainable"):
            validate_trainable_artifact(
                {"trainable_decision": {"model": {"rt_half_window_sec": 10}}}
            )
        for forbidden in (
            {"method": {}},
            {"extraction": {}},
            {"acquisition": {}},
            {"dia": {}},
            {"model": {"correlation_min_relative_intensity": 0.05}},
            {"model": {"fragment_matching": "exclusive_one_to_one"}},
            {"model": {"fragment_correlation_operator": "pearson"}},
            {"model": {"precursor_eic_tolerance": 10.0}},
            {"model": {"fragment_alignment_tolerance_unit": "Da"}},
            {"model": {"intensity_power": 0.5}},
            {"frozen_input_generation": {}},
            {"parameter_provenance": {}},
            {"guardrails": {}},
        ):
            with (
                self.subTest(forbidden=forbidden),
                self.assertRaisesRegex(ValueError, "non-trainable"),
            ):
                validate_trainable_artifact({"trainable_decision": {"model": {"probe": forbidden}}})
        with self.assertRaisesRegex(ValueError, "unsupported.*root keys"):
            validate_trainable_artifact({"trainable_decision": {"future_knob": 1}})

    def test_machine_readable_shadow_contract_matches_code_and_no_v3_is_added(self):
        contract = json.loads(
            Path("MSAI/standards/ms2_ml_shadow_contract_v1.json").read_text(encoding="utf-8")
        )
        self.assertEqual(tuple(contract["truth_labels"]), TRUTH_LABELS)
        self.assertEqual(tuple(contract["ml_statuses"]), ML_DIAGNOSTIC_STATUSES)
        self.assertEqual(tuple(contract["shadow_output_fields"]), ML_OUTPUT_FIELDS)
        self.assertEqual(tuple(contract["features_in_order"]), FEATURE_NAMES)
        source_schema = contract["frozen_input_generation"]["source_analysis_parameter_schema"]
        self.assertEqual(set(source_schema["allowed_keys"]), ALLOWED_SOURCE_PARAMETER_KEYS)
        self.assertEqual(
            set(source_schema["downstream_decision_guardrail_atomic_group"]),
            SOURCE_DECISION_GUARDRAIL_KEYS,
        )
        extraction = contract["frozen_input_generation"][
            "extraction_never_trained_in_this_model_family"
        ]
        self.assertEqual(extraction["close_peak_window_rule"], FROZEN_CLOSE_PEAK_RULE)
        self.assertEqual(
            {value.casefold() for value in contract["forbidden_feature_sources"]},
            LEAKAGE_DENYLIST,
        )
        self.assertFalse(list(Path("MSAI/standards").glob("ms2_diagnostic_standard_v3*")))

    def test_training_cli_exposes_no_method_or_extraction_tuning_flags(self):
        output = io.StringIO()
        with redirect_stdout(output), self.assertRaises(SystemExit) as raised:
            train_main(["--help"])
        self.assertEqual(raised.exception.code, 0)
        help_text = output.getvalue()
        self.assertIn("--specificity-target", help_text)
        for forbidden in (
            "--rt-half-window",
            "--min-fragment-correlation",
            "--fragment-mz-tol",
            "--min-fragment-intensity",
            "--consensus-scans",
            "--dia-iso-win",
        ):
            self.assertNotIn(forbidden, help_text)


class GroupIsolationContractTests(unittest.TestCase):
    def test_locked_batch_and_cv_purge_compound_overlap(self):
        examples = (
            _example("a-pos", "A", "batch-a", 1),
            _example("a-neg", "B", "batch-a", 0),
            _example("b-pos", "C", "batch-b", 1),
            _example("b-neg", "D", "batch-b", 0),
            _example("c-pos", "E", "batch-c", 1),
            _example("c-neg", "F", "batch-c", 0),
            _example("overlap", "LOCKED", "batch-c", 1),
            _example("final-pos", "LOCKED", "batch-final", 1),
            _example("final-neg", "G", "batch-final", 0),
        )
        split = build_locked_batch_split(examples, "batch-final")
        self.assertEqual(split.purged_indices, (6,))
        assert_group_isolation(examples, split.development_indices, split.final_indices)
        folds = build_leave_one_batch_out_folds(examples, split.development_indices)
        self.assertEqual(len(folds), 3)
        for fold in folds:
            assert_group_isolation(examples, fold.train_indices, fold.validation_indices)

    def test_row_random_or_single_batch_fallback_does_not_exist(self):
        examples = (
            _example("one", "A", "batch-a", 1),
            _example("two", "B", "batch-a", 0),
        )
        with self.assertRaisesRegex(ValueError, "at least two eligible batches"):
            build_leave_one_batch_out_folds(examples)


if __name__ == "__main__":
    unittest.main()
