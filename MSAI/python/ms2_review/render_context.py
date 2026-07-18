"""Shared deterministic metadata assembly for MS2 review renderers."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from typing import Any


def build_review_render_metadata(
    *,
    target_uid: str,
    config_hash: str,
    standard_id: str,
    source_machine_id: str,
    source_machine_label: str,
    source_machine_interpretation: str,
    source_pool_id: str,
    source_pooled_well: str,
    parameters: Mapping[str, Any],
    thresholds: Mapping[str, Any],
    dia_windows: Sequence[Mapping[str, Any]],
    acquisition_context: Mapping[str, Any],
    method_profile: Mapping[str, Any],
    acquisition_reconciliation: Mapping[str, Any],
    model_context: Mapping[str, Any],
) -> dict[str, Any]:
    """Build the exact context consumed by both SVG and PNG renderers."""

    source_context = (
        f"{source_machine_id} label: {source_machine_label or '—'} "
        f"({source_machine_interpretation})"
        if source_machine_id
        else "no source-machine label"
    )
    return {
        "target_uid": target_uid,
        "config_hash": config_hash,
        "standard_id": standard_id,
        "source_context": source_context,
        "source_pool_id": source_pool_id,
        "source_pooled_well": source_pooled_well,
        "fragment_mz_tol": float(parameters.get("fragment_mz_tol", 0.01)),
        "fragment_mz_tol_unit": str(parameters.get("fragment_mz_tol_unit", "Da")),
        "min_relative_intensity": float(parameters.get("min_relative_intensity", 0.01)),
        "min_cosine": float(thresholds.get("min_cosine", 0.7)),
        "min_matched_peaks": int(thresholds.get("min_matched_peaks", 6)),
        "min_explained_intensity": float(thresholds.get("min_explained_intensity", 0.5)),
        "rt_half_window_sec": parameters.get("rt_half_window_sec", 10.0),
        "min_fragment_correlation": parameters.get("min_fragment_correlation", 0.9),
        "fragment_correlation_mode": parameters.get("fragment_correlation_mode", "full_window"),
        "correlation_min_relative_intensity": parameters.get(
            "correlation_min_relative_intensity", 0.05
        ),
        "min_correlation_scans": parameters.get("min_correlation_scans", 5),
        "max_fragment_apex_offset_scans": parameters.get("max_fragment_apex_offset_scans", 1),
        "min_consecutive_fragment_scans": parameters.get("min_consecutive_fragment_scans", 3),
        "dia_windows": list(dia_windows),
        "acquisition_context": dict(acquisition_context),
        "method_profile": dict(method_profile),
        "acquisition_reconciliation": dict(acquisition_reconciliation),
        "ml_model_id": model_context.get("model_id") or "",
        "ml_model_sha256": model_context.get("sha256") or "",
    }


def review_config_hash(
    *,
    standard_id: str,
    parameters: Mapping[str, Any],
    thresholds: Mapping[str, Any],
    method_profile: Mapping[str, Any],
    model_context: Mapping[str, Any],
) -> str:
    """Return the stable review-render configuration fingerprint.

    A model artifact SHA-256 already commits to its complete trainable decision,
    so the render configuration needs only the declared model ID and digest.
    Keeping this function beside the shared metadata builder lets exporters and
    package validators independently derive the same value.
    """

    payload = json.dumps(
        {
            "standard_id": standard_id,
            "parameters": dict(parameters),
            "thresholds": dict(thresholds),
            "method_profile": dict(method_profile),
            "model": {
                "model_id": model_context.get("model_id") or "",
                "sha256": model_context.get("sha256") or "",
                # Preserve the v1 fingerprint shape. The artifact SHA-256 above
                # commits to the complete trainable decision without trusting a
                # mutable copy from review metadata.
                "trainable_decision": {},
            },
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode()).hexdigest()[:10]


__all__ = ["build_review_render_metadata", "review_config_hash"]
