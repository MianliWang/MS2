"""Acquisition provenance helpers for MS2 review artifacts.

The raw file, a manuscript method, and the parameters used by the analysis are
different evidence layers.  This module keeps them separate and reports
disagreements without allowing a reference method to overwrite raw metadata.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

try:
    from ..ms2.raw_metadata import raw_acquisition_context
except ImportError:
    from ms2.raw_metadata import raw_acquisition_context  # type: ignore


def load_method_profile(path: str | Path | None) -> dict:
    """Load a versioned reference-method profile, or return an empty mapping."""

    if not path:
        return {}
    profile_path = Path(path)
    with profile_path.open(encoding="utf-8") as handle:
        profile = json.load(handle)
    if not isinstance(profile, dict) or not str(profile.get("profile_id") or "").strip():
        raise ValueError(f"Method profile must define profile_id: {profile_path}")
    return profile


def reconcile_acquisition(
    observed_raw: dict,
    reference_method: dict,
    *,
    analysis_parameters: dict | None = None,
    dia_windows: list[dict] | None = None,
    raw_input_count: int = 1,
) -> dict:
    """Compare provenance layers without claiming that they describe one run."""

    if not reference_method:
        return {
            "status": "raw_only" if observed_raw else "metadata_unavailable",
            "issue_codes": [],
            "matches": [],
            "comparisons": [],
        }

    parameters = analysis_parameters or {}
    dia = ((reference_method.get("acquisition") or {}).get("dia") or {})
    processing = reference_method.get("processing_reference") or {}
    applicability = reference_method.get("applicability") or {}
    issue_codes: list[str] = []
    matches: list[str] = []
    comparisons: list[dict] = []

    if not applicability.get("exact_current_run_confirmed", False):
        issue_codes.append("REFERENCE_PROFILE_NOT_CONFIRMED_FOR_RUN")

    internal_conflicts = reference_method.get("known_internal_conflicts") or []
    if internal_conflicts:
        issue_codes.append("REFERENCE_INTERNAL_DIA_WIDTH_CONFLICT")

    raw_widths = _numeric_observed_values(observed_raw, "windowWideness")
    method_width = _number(dia.get("nominal_isolation_width_mz"))
    preprocessing_width = _number(dia.get("preprocessing_fixed_isolation_width_da"))
    if raw_widths and any(
        value is not None and not any(_close(raw, value) for raw in raw_widths)
        for value in (method_width, preprocessing_width)
    ):
        issue_codes.append("RAW_VS_REFERENCE_DIA_WIDTH_DIFFER")
    comparisons.append(
        {
            "field": "dia_isolation_width",
            "raw_observed": raw_widths,
            "reference_method_nominal_mz": method_width,
            "reference_preprocessing_fixed_da": preprocessing_width,
            "interpretation": "different method versions or unresolved manuscript wording",
        }
    )

    raw_collision = _numeric_observed_values(observed_raw, "collisionEnergy")
    reference_collision = [
        value for value in (_number(item) for item in dia.get("stepped_hcd_percent", []))
        if value is not None
    ]
    if raw_collision and reference_collision and set(raw_collision) != set(reference_collision):
        issue_codes.append("RAW_VS_REFERENCE_COLLISION_ENERGY_DIFFER")
    comparisons.append(
        {
            "field": "collision_energy",
            "raw_observed": raw_collision,
            "reference_stepped_hcd_percent": reference_collision,
            "interpretation": "raw field and documented stepped method are retained separately",
        }
    )

    raw_activation = _observed_values(observed_raw, "activationMethod")
    reference_activation = str(dia.get("activation_method") or "")
    if reference_activation and raw_activation and all(
        value.upper() == reference_activation.upper() for value in raw_activation
    ):
        matches.append("ACTIVATION_METHOD_HCD")

    expected_files = _integer(dia.get("files_per_sample"))
    if expected_files and raw_input_count < expected_files:
        issue_codes.append("REFERENCE_THREE_DIA_FILES_CURRENT_INPUT_ONE")
    comparisons.append(
        {
            "field": "dia_file_count",
            "current_raw_inputs": raw_input_count,
            "reference_files_per_sample": expected_files,
            "interpretation": "verify whether companion DIA files belong to this sample",
        }
    )

    reference_correlation = _number(processing.get("min_chromatographic_correlation"))
    analysis_correlation = _number(parameters.get("min_fragment_correlation"))
    if (
        reference_correlation is not None
        and analysis_correlation is not None
        and _close(reference_correlation, analysis_correlation)
    ):
        matches.append("CHROMATOGRAPHIC_CORRELATION_0_9")
    comparisons.append(
        {
            "field": "chromatographic_correlation",
            "analysis": analysis_correlation,
            "reference": reference_correlation,
        }
    )

    reference_rt = _number(processing.get("precursor_fragment_rt_alignment_sec"))
    analysis_rt = _number(parameters.get("rt_half_window_sec"))
    if reference_rt is not None and analysis_rt is not None and not _close(reference_rt, analysis_rt):
        issue_codes.append("REFERENCE_RT_ALIGNMENT_DIFFERS_FROM_EXTRACTION_WINDOW")
    comparisons.append(
        {
            "field": "rt_seconds",
            "analysis_extraction_half_window": analysis_rt,
            "reference_alignment_limit": reference_rt,
            "interpretation": "related sensitivity parameter, not necessarily the same operation",
        }
    )

    lab_ppm = _number(((processing.get("mass_tolerance_ppm") or {}).get("value")))
    precursor_value = parameters.get("precursor_eic_mz_tol", parameters.get("mz_tol"))
    precursor_unit = parameters.get("precursor_eic_mz_tol_unit", parameters.get("mz_tol_unit"))
    fragment_value = parameters.get("fragment_eic_mz_tol", parameters.get("mz_tol"))
    fragment_unit = parameters.get("fragment_eic_mz_tol_unit", parameters.get("mz_tol_unit"))
    analysis_ppm = _number(precursor_value) if str(precursor_unit or "") == "ppm" else None
    fragment_analysis_ppm = (
        _number(fragment_value) if str(fragment_unit or "") == "ppm" else None
    )
    if lab_ppm is not None and analysis_ppm is not None and not _close(lab_ppm, analysis_ppm):
        issue_codes.append("LAB_NOTE_3PPM_NOT_APPLIED_TO_CURRENT_MS2_EXTRACTION")
    comparisons.append(
        {
            "field": "mass_tolerance_ppm",
            "analysis_precursor_eic": analysis_ppm,
            "analysis_fragment_eic": fragment_analysis_ppm,
            "author_lab_note": lab_ppm,
            "interpretation": "lab-note scope is unconfirmed; shadow validation required",
        }
    )

    if dia_windows:
        bounds = [
            (_number(window.get("lower")), _number(window.get("upper")))
            for window in dia_windows
        ]
        finite = [(lower, upper) for lower, upper in bounds if lower is not None and upper is not None]
        if finite:
            comparisons.append(
                {
                    "field": "current_acquired_dia_coverage",
                    "lower": min(lower for lower, _upper in finite),
                    "upper": max(upper for _lower, upper in finite),
                    "window_count": len(finite),
                    "source": "analysis sidecar derived from raw isolation bounds",
                }
            )

    return {
        "status": "reference_only_with_discrepancies" if issue_codes else "compatible_reference",
        "profile_id": reference_method.get("profile_id"),
        "issue_codes": list(dict.fromkeys(issue_codes)),
        "matches": list(dict.fromkeys(matches)),
        "comparisons": comparisons,
    }


def _observed_values(context: dict, key: str) -> list[str]:
    values = (context.get("observed_values") or {}).get(key) or []
    return [str(item.get("value")) for item in values if item.get("value") is not None]


def _numeric_observed_values(context: dict, key: str) -> list[float]:
    return [value for value in (_number(item) for item in _observed_values(context, key)) if value is not None]


def _number(value) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _integer(value) -> int | None:
    number = _number(value)
    return None if number is None else int(number)


def _close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=1e-9, abs_tol=1e-9)


__all__ = ["load_method_profile", "raw_acquisition_context", "reconcile_acquisition"]
