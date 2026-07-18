"""Safety contracts for the MS2 shadow classifier.

This module deliberately separates immutable extraction/Method parameters from
the trainable decision layer.  Nothing here infers truth from the existing v2
rule output: supervised labels must use the four explicit values below.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Iterable, Mapping
from typing import Any

try:
    from ..ms2.similarity import parse_fragment_string
except ImportError:  # pragma: no cover - direct PYTHONPATH=MSAI/python execution
    from ms2.similarity import parse_fragment_string  # type: ignore[import-not-found]

POSITIVE_TRUTH_LABEL = "positive_same_compound"
NEGATIVE_TRUTH_LABEL = "negative_different_or_interference"
UNCERTAIN_TRUTH_LABEL = "uncertain"
NOT_EVALUABLE_TRUTH_LABEL = "not_evaluable"

TRUTH_LABELS = (
    POSITIVE_TRUTH_LABEL,
    NEGATIVE_TRUTH_LABEL,
    UNCERTAIN_TRUTH_LABEL,
    NOT_EVALUABLE_TRUTH_LABEL,
)
ELIGIBLE_TRUTH_LABELS = (POSITIVE_TRUTH_LABEL, NEGATIVE_TRUTH_LABEL)

ML_DIAGNOSTIC_STATUSES = (
    "supported_same_compound",
    "conflicting_spectra",
    "insufficient_evidence",
    "not_evaluable",
)
V2_DIAGNOSTIC_STATUSES = (
    "supported_same_compound",
    "conflicting_spectra",
    "insufficient_evidence",
    "not_evaluable",
    "not_attempted",
)
ML_OUTPUT_FIELDS = (
    "ml_same_compound_probability",
    "ml_diagnostic_status",
    "ml_model_id",
    "ml_abstention_reason",
    "ml_review_flags",
)
SHADOW_COMPARISON_FIELDS = (
    "comparison_change_level",
    "comparison_baseline_status",
    "comparison_v2_normalized_status",
    "comparison_candidate_status",
    "comparison_review_tier",
    "comparison_review_flags",
)

# Any future row-based feature builder must apply this case-insensitive denylist.
# The current feature builder is safer still: it accepts only two spectra.
LEAKAGE_DENYLIST = frozenset(
    {
        "ms2_diagnostic_status",
        "enantiomer_pair_status",
        "ms2_issue_codes",
        "ig",
        "supported_same_compound",
        "conflicting_spectra",
        "candidate_enantiomer_pair",
        "ms2_not_similar",
        "insufficient_ms2_evidence",
        "compound_identity_status",
        "truth_label",
        "positive_same_compound",
        "negative_different_or_interference",
        "uncertain",
        "not_evaluable",
        "ml_same_compound_probability",
        "ml_diagnostic_status",
        "ml_model_id",
        "ml_abstention_reason",
        "ml_review_flags",
    }
)

FROZEN_CLOSE_PEAK_RULE = "half_width=min(10.0, 0.4*rt_separation_seconds)"
FROZEN_SOURCE_PARAMETERS: dict[str, float | int | str] = {
    "mz_column": "MZ",
    "peak_a_column": "Peak1",
    "peak_b_column": "Peak2",
    "rt_unit": "min",
    "rt_half_window_sec": 10.0,
    "mz_tol": 10.0,
    "mz_tol_unit": "ppm",
    "mz_tol_legacy_scope": "shared precursor and fragment EIC default",
    "precursor_eic_mz_tol": 10.0,
    "precursor_eic_mz_tol_unit": "ppm",
    "fragment_eic_mz_tol": 10.0,
    "fragment_eic_mz_tol_unit": "ppm",
    "dia_iso_win_fallback": 15.0,
    "fragment_mz_tol": 0.01,
    "fragment_mz_tol_unit": "Da",
    "min_relative_intensity": 0.01,
    "min_fragment_intensity": 2000.0,
    "min_fragment_relative_intensity": 0.0,
    "min_fragment_correlation": 0.9,
    "fragment_correlation_mode": "full_window",
    "correlation_min_relative_intensity": 0.05,
    "min_correlation_scans": 5,
    "max_fragment_apex_offset_scans": 1,
    "min_consecutive_fragment_scans": 3,
    "consensus_scans": 1,
}

# These four values are emitted beside the extraction parameters by the v2
# pipeline, but they act only after spectra have been generated.  They are
# therefore permitted in a source sidecar without being frozen: this model
# family may replace them using independent truth.  Every other key in the
# analysis-parameter mapping is rejected so a future extraction knob cannot
# silently enter an already-trained model family.
SOURCE_DECISION_GUARDRAIL_KEYS = frozenset(
    {
        "min_cosine",
        "min_matched_peaks",
        "min_explained_intensity",
        "min_entropy_similarity",
    }
)
ALLOWED_SOURCE_PARAMETER_KEYS = frozenset(FROZEN_SOURCE_PARAMETERS) | (
    SOURCE_DECISION_GUARDRAIL_KEYS
)

FROZEN_STANDARD_PREPROCESSING: dict[str, float | str] = {
    "precursor_eic_tolerance_ppm": 10.0,
    "rt_half_window_seconds": 10.0,
    "close_peak_window_rule": FROZEN_CLOSE_PEAK_RULE,
    "fragment_coelution_min_pearson": 0.9,
    "fragment_coelution_operator": ">",
    "spectrum_min_relative_intensity": 0.01,
    "intensity_transform_for_cosine": "square_root",
    "fragment_alignment_tolerance_da": 0.01,
    "fragment_matching": "exclusive_one_to_one",
}

MODEL_SCHEMA_VERSION = "ms2_shadow_trainable_decision_v1"
MIN_SPECIFICITY_TARGET = 0.95
TRAINABLE_DECISION_KEYS = frozenset(
    {
        "schema_version",
        "model_id",
        "selected_model",
        "feature_names",
        "model",
        "decision_thresholds",
        "training_feature_ranges",
        "candidates",
        "selection",
        "development_validation",
        "locked_final_evaluation",
        "split_manifest",
        "annotation_summary",
        "provenance",
        "model_card",
        "limitations",
    }
)

_EXPLICIT_EIC_KEYS = frozenset(
    {
        "precursor_eic_mz_tol",
        "precursor_eic_mz_tol_unit",
        "fragment_eic_mz_tol",
        "fragment_eic_mz_tol_unit",
    }
)
_SOURCE_DECISION_GUARDRAIL_KEYS = SOURCE_DECISION_GUARDRAIL_KEYS
_FORBIDDEN_SECTION_KEYS = frozenset(
    {
        "method",
        "methods",
        "method_parameters",
        "reference_method",
        "extraction",
        "extraction_parameters",
        "acquisition",
        "acquisition_parameters",
        "dia",
        "dia_parameters",
        "preprocessing",
        "frozen_input_generation",
        "paper_method_never_trained",
        "extraction_never_trained_in_this_model_family",
        "parameter_provenance",
        "guardrails",
    }
)

# Keys that may describe, but must never become part of, a trainable decision
# payload.  Provenance should record hashes/IDs of the immutable source sidecar
# and standard instead of copying these values into the model artifact.
NON_TRAINABLE_PARAMETER_KEYS = frozenset(
    {
        *FROZEN_SOURCE_PARAMETERS,
        *FROZEN_STANDARD_PREPROCESSING,
        "rt_half_window_seconds",
        "precursor_eic_tolerance_ppm",
        "precursor_eic_tolerance",
        "precursor_eic_tolerance_unit",
        "fragment_eic_tolerance_ppm",
        "fragment_eic_tolerance",
        "fragment_eic_tolerance_unit",
        "fragment_alignment_tolerance_da",
        "fragment_alignment_tolerance",
        "fragment_alignment_tolerance_unit",
        "fragment_coelution_min_pearson",
        "fragment_correlation_operator",
        "fragment_coelution_operator",
        "close_peak_window_rule",
        "anti_overlap_rule",
        "intensity_power",
        "consensus_scan_count",
        "mz_tol",
        "mz_tol_unit",
        "dia_iso_win",
        "dia_iso_win_fallback",
        "dia_window",
        "collision_energy",
        "collisionenergy",
        "activation_method",
        "activationmethod",
        "acquisition_parameters",
        "method_parameters",
        "extraction_parameters",
        *_FORBIDDEN_SECTION_KEYS,
    }
)


def normalize_truth_label(value: Any) -> str:
    """Return an exact truth label, rejecting blanks, aliases, and rule labels."""

    label = str(value).strip()
    if label not in TRUTH_LABELS:
        allowed = ", ".join(TRUTH_LABELS)
        raise ValueError(f"Unknown truth label {label!r}; expected one of: {allowed}.")
    return label


def truth_to_binary(label: Any) -> int | None:
    """Map eligible truth to 1/0 and explicitly excluded truth to ``None``."""

    normalized = normalize_truth_label(label)
    if normalized == POSITIVE_TRUTH_LABEL:
        return 1
    if normalized == NEGATIVE_TRUTH_LABEL:
        return 0
    return None


def validate_truth_row(
    row: Mapping[str, Any],
    label_column: str = "truth_label",
    spectrum_a_column: str = "peak_a_MS2",
    spectrum_b_column: str = "peak_b_MS2",
) -> str:
    """Validate one annotation without treating a missing spectrum as negative.

    An eligible positive or negative annotation is a contradiction when either
    independently extracted spectrum is empty.  ``uncertain`` and
    ``not_evaluable`` rows are valid but remain excluded from fitting.
    """

    if label_column not in row:
        raise ValueError(f"Truth data must contain {label_column!r}.")
    label = normalize_truth_label(row[label_column])
    if label in ELIGIBLE_TRUTH_LABELS:
        missing = [
            name
            for name in (spectrum_a_column, spectrum_b_column)
            if not _has_valid_spectrum(row.get(name))
        ]
        if missing:
            columns = ", ".join(missing)
            raise ValueError(
                f"Eligible truth label {label!r} requires bilateral spectra; missing: {columns}."
            )
    return label


def filter_labeled_rows(
    rows: Iterable[Mapping[str, Any]],
    label_column: str = "truth_label",
    spectrum_a_column: str = "peak_a_MS2",
    spectrum_b_column: str = "peak_b_MS2",
) -> list[dict[str, Any]]:
    """Validate all annotations and return copies of only binary-truth rows."""

    eligible: list[dict[str, Any]] = []
    for row in rows:
        label = validate_truth_row(row, label_column, spectrum_a_column, spectrum_b_column)
        if label in ELIGIBLE_TRUTH_LABELS:
            eligible.append(dict(row))
    return eligible


def validate_no_leakage(feature_names: Iterable[str]) -> tuple[str, ...]:
    """Reject any requested feature whose name is a known derived-rule leak."""

    names = tuple(str(name) for name in feature_names)
    leaked = sorted({name for name in names if name.casefold() in LEAKAGE_DENYLIST})
    if leaked:
        raise ValueError(f"Forbidden label/rule leakage features: {', '.join(leaked)}.")
    return names


def validate_frozen_source(
    sidecar: Mapping[str, Any],
    standard: Mapping[str, Any] | None = None,
) -> dict[str, float | int | str]:
    """Validate that spectra came from the frozen 10-second extraction contract.

    ``sidecar`` may be either the run metadata (with ``parameters``), its
    ``provenance_layers.analysis_parameters`` view, or the parameter mapping
    itself.  The optional v2 ``standard`` additionally locks the close-peak
    anti-overlap rule, which run metadata does not currently duplicate.
    """

    parameters = _analysis_parameters(sidecar)
    unexpected_parameters = sorted(set(parameters).difference(ALLOWED_SOURCE_PARAMETER_KEYS))
    if unexpected_parameters:
        raise ValueError(
            "Source sidecar contains unknown analysis/extraction parameters; "
            "a model-family contract update is required before use: "
            + ", ".join(unexpected_parameters)
        )
    explicit_eic_present = _EXPLICIT_EIC_KEYS.intersection(parameters)
    if explicit_eic_present and explicit_eic_present != _EXPLICIT_EIC_KEYS:
        missing = sorted(_EXPLICIT_EIC_KEYS - explicit_eic_present)
        raise ValueError(
            "Explicit precursor/fragment EIC fields must be supplied all-or-none; missing: "
            + ", ".join(missing)
        )

    decision_guardrails_present = _SOURCE_DECISION_GUARDRAIL_KEYS.intersection(parameters)
    if decision_guardrails_present and decision_guardrails_present != (
        _SOURCE_DECISION_GUARDRAIL_KEYS
    ):
        missing = sorted(_SOURCE_DECISION_GUARDRAIL_KEYS - decision_guardrails_present)
        raise ValueError(
            "Downstream diagnostic guardrail fields must be supplied all-or-none; missing: "
            + ", ".join(missing)
        )
    if decision_guardrails_present:
        _validate_source_decision_guardrails(parameters)

    resolved: dict[str, float | int | str] = {}
    for key, expected in FROZEN_SOURCE_PARAMETERS.items():
        value = parameters.get(key)
        # Older sidecars record a single shared 10 ppm EIC tolerance.  It is
        # acceptable only when all four explicit precursor/fragment fields are
        # absent; a partially explicit contract is ambiguous and rejected.
        if value is None and not explicit_eic_present and key in _EXPLICIT_EIC_KEYS:
            value = parameters.get("mz_tol_unit" if key.endswith("_unit") else "mz_tol")
        if value is None:
            raise ValueError(f"Source sidecar is missing frozen parameter {key!r}.")
        if not _same_parameter(value, expected):
            raise ValueError(
                f"Frozen source parameter {key!r} must be {expected!r}, got {value!r}."
            )
        resolved[key] = expected

    if standard is not None:
        if standard.get("standard_id") != "MSAI-MS2-DIAGNOSTIC-v2":
            raise ValueError("Shadow ML requires diagnostic standard MSAI-MS2-DIAGNOSTIC-v2.")
        preprocessing = standard.get("preprocessing")
        if not isinstance(preprocessing, Mapping):
            raise ValueError("Diagnostic standard must contain a preprocessing mapping.")
        if set(preprocessing) != set(FROZEN_STANDARD_PREPROCESSING):
            missing = sorted(set(FROZEN_STANDARD_PREPROCESSING) - set(preprocessing))
            unexpected = sorted(set(preprocessing) - set(FROZEN_STANDARD_PREPROCESSING))
            details = []
            if missing:
                details.append("missing=" + ",".join(missing))
            if unexpected:
                details.append("unexpected=" + ",".join(unexpected))
            raise ValueError("Diagnostic v2 preprocessing schema mismatch: " + "; ".join(details))
        for key, expected in FROZEN_STANDARD_PREPROCESSING.items():
            value = preprocessing.get(key)
            if not _same_parameter(value, expected):
                raise ValueError(
                    f"Frozen v2 preprocessing parameter {key!r} must be "
                    f"{expected!r}, got {value!r}."
                )
    return resolved


def make_trainable_artifact(payload: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    """Wrap and validate a JSON model payload under its sole allowed top-level key."""

    artifact: dict[str, dict[str, Any]] = {"trainable_decision": dict(payload)}
    validate_trainable_artifact(artifact)
    return artifact


def validate_trainable_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    """Enforce the one-section artifact and reject copied immutable parameters."""

    if set(artifact) != {"trainable_decision"}:
        raise ValueError("Training artifact top level must contain only 'trainable_decision'.")
    payload = artifact.get("trainable_decision")
    if not isinstance(payload, Mapping):
        raise ValueError("'trainable_decision' must be a mapping.")
    unexpected = sorted(set(payload) - TRAINABLE_DECISION_KEYS)
    if unexpected:
        raise ValueError(
            "Training artifact contains unsupported trainable_decision root keys: "
            + ", ".join(unexpected)
        )
    forbidden = sorted(_find_forbidden_keys(payload))
    if forbidden:
        raise ValueError(
            "Training artifact contains non-trainable Method/extraction keys: "
            + ", ".join(forbidden)
        )
    return dict(payload)


def canonical_model_id(payload: Mapping[str, Any]) -> str:
    """Return the deterministic ID used by trained shadow artifacts."""

    decision_without_id = {str(key): value for key, value in payload.items() if key != "model_id"}
    try:
        encoded = json.dumps(
            decision_without_id,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise ValueError("Model payload must be canonical JSON data.") from error
    return "ms2-shadow-" + hashlib.sha256(encoded).hexdigest()[:20]


def _analysis_parameters(sidecar: Mapping[str, Any]) -> Mapping[str, Any]:
    parameters: Mapping[str, Any] | None = None
    analysis: Mapping[str, Any] | None = None
    if "parameters" in sidecar:
        raw_parameters = sidecar.get("parameters")
        if not isinstance(raw_parameters, Mapping):
            raise ValueError("Source sidecar 'parameters' must be a mapping.")
        parameters = raw_parameters
    if "provenance_layers" in sidecar:
        provenance = sidecar.get("provenance_layers")
        if not isinstance(provenance, Mapping):
            raise ValueError("Source sidecar 'provenance_layers' must be a mapping.")
        if "analysis_parameters" in provenance:
            raw_analysis = provenance.get("analysis_parameters")
            if not isinstance(raw_analysis, Mapping):
                raise ValueError(
                    "Source sidecar provenance_layers.analysis_parameters must be a mapping."
                )
            analysis = raw_analysis
    if parameters is not None and analysis is not None and dict(parameters) != dict(analysis):
        raise ValueError(
            "Source sidecar parameters conflict with provenance_layers.analysis_parameters."
        )
    if parameters is not None:
        return parameters
    if analysis is not None:
        return analysis
    return sidecar


def _same_parameter(observed: Any, expected: float | int | str) -> bool:
    if isinstance(expected, str):
        return isinstance(observed, str) and observed == expected
    if isinstance(observed, bool):
        return False
    try:
        numeric = float(observed)
    except (TypeError, ValueError):
        return False
    return math.isfinite(numeric) and math.isclose(
        numeric,
        float(expected),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def _validate_source_decision_guardrails(parameters: Mapping[str, Any]) -> None:
    cosine = _bounded_source_number(parameters.get("min_cosine"), "min_cosine")
    explained = _bounded_source_number(
        parameters.get("min_explained_intensity"),
        "min_explained_intensity",
    )
    matched = parameters.get("min_matched_peaks")
    if isinstance(matched, bool) or not isinstance(matched, int) or matched < 1:
        raise ValueError(
            "Source decision guardrail 'min_matched_peaks' must be a positive integer."
        )
    entropy = parameters.get("min_entropy_similarity")
    if entropy is not None:
        _bounded_source_number(entropy, "min_entropy_similarity")
    # Bind local names so type-checkers and future refactors cannot accidentally
    # remove validation of either bounded field.
    if not (0.0 <= cosine <= 1.0 and 0.0 <= explained <= 1.0):  # pragma: no cover
        raise AssertionError("Validated source decision guardrail left its bounds.")


def _bounded_source_number(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ValueError(f"Source decision guardrail {name!r} must be numeric in [0, 1].")
    numeric = float(value)
    if not math.isfinite(numeric) or not 0 <= numeric <= 1:
        raise ValueError(f"Source decision guardrail {name!r} must be numeric in [0, 1].")
    return numeric


def _has_valid_spectrum(value: Any) -> bool:
    if isinstance(value, str) or value is None:
        return bool(parse_fragment_string(value))
    try:
        return any(
            math.isfinite(float(mz))
            and math.isfinite(float(intensity))
            and float(mz) > 0
            and float(intensity) > 0
            for mz, intensity in value
        )
    except (TypeError, ValueError):
        return False


def _find_forbidden_keys(value: Any, prefix: str = "") -> set[str]:
    found: set[str] = set()
    if isinstance(value, Mapping):
        for raw_key, child in value.items():
            key = str(raw_key)
            path = f"{prefix}.{key}" if prefix else key
            canonical_key = _canonical_key(key)
            if canonical_key in NON_TRAINABLE_PARAMETER_KEYS or _is_forbidden_section_key(
                canonical_key
            ):
                found.add(path)
            found.update(_find_forbidden_keys(child, path))
    elif isinstance(value, list | tuple):
        for index, child in enumerate(value):
            found.update(_find_forbidden_keys(child, f"{prefix}[{index}]"))
    return found


def _canonical_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.casefold()).strip("_")


def _is_forbidden_section_key(key: str) -> bool:
    if key in _FORBIDDEN_SECTION_KEYS:
        return True
    tokens = ("method", "extraction", "acquisition", "dia")
    suffixes = ("parameters", "settings", "config", "profile")
    return any(key == f"{token}_{suffix}" for token in tokens for suffix in suffixes)


__all__ = [
    "ALLOWED_SOURCE_PARAMETER_KEYS",
    "ELIGIBLE_TRUTH_LABELS",
    "FROZEN_CLOSE_PEAK_RULE",
    "FROZEN_SOURCE_PARAMETERS",
    "FROZEN_STANDARD_PREPROCESSING",
    "LEAKAGE_DENYLIST",
    "MIN_SPECIFICITY_TARGET",
    "ML_DIAGNOSTIC_STATUSES",
    "ML_OUTPUT_FIELDS",
    "MODEL_SCHEMA_VERSION",
    "NEGATIVE_TRUTH_LABEL",
    "NON_TRAINABLE_PARAMETER_KEYS",
    "NOT_EVALUABLE_TRUTH_LABEL",
    "POSITIVE_TRUTH_LABEL",
    "SHADOW_COMPARISON_FIELDS",
    "SOURCE_DECISION_GUARDRAIL_KEYS",
    "TRAINABLE_DECISION_KEYS",
    "TRUTH_LABELS",
    "UNCERTAIN_TRUTH_LABEL",
    "V2_DIAGNOSTIC_STATUSES",
    "canonical_model_id",
    "filter_labeled_rows",
    "make_trainable_artifact",
    "normalize_truth_label",
    "truth_to_binary",
    "validate_frozen_source",
    "validate_no_leakage",
    "validate_trainable_artifact",
    "validate_truth_row",
]
