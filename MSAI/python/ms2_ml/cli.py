"""Command-line entry points for the MS2 ML shadow workflow."""

from __future__ import annotations

import argparse
import json
from collections.abc import Mapping
from pathlib import Path

from ..ms2.tables import read_table
from .contracts import validate_frozen_source


def train_main(argv=None) -> int:
    """Fit the prespecified baselines on independent grouped truth data."""

    from .training import train_and_write_shadow_model

    parser = argparse.ArgumentParser(
        description=(
            "Train threshold and L2 Logistic Regression MS2 decision models. "
            "Method and extraction parameters are fixed and are not CLI options."
        )
    )
    parser.add_argument("--input", required=True, help="Independently labelled CSV/XLSX")
    parser.add_argument(
        "--source-sidecar",
        required=True,
        help="Per-batch frozen-source manifest JSON (batch_sources)",
    )
    parser.add_argument("--output", required=True, help="trainable_decision JSON artifact")
    parser.add_argument("--standard", default=str(_default_standard()))
    parser.add_argument("--label-column", default="truth_label")
    parser.add_argument("--compound-column", default="Compound_ID")
    parser.add_argument("--batch-column", default="batch_id")
    parser.add_argument("--row-id-column", default="target_uid")
    parser.add_argument("--locked-final-batch", required=True)
    parser.add_argument("--seed", type=int, default=20260129)
    parser.add_argument("--specificity-target", type=float, default=0.95)
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    args = parser.parse_args(argv)

    input_path = Path(args.input)
    sidecar_path = Path(args.source_sidecar)
    standard_path = Path(args.standard)
    sidecar = _read_mapping(sidecar_path)
    standard = _read_mapping(standard_path)
    # ``train_shadow_model`` validates every eligible batch entry in the
    # manifest.  Validating this top-level manifest as if it were one run would
    # incorrectly let a single sidecar stand in for a multi-batch truth set.
    artifact = train_and_write_shadow_model(
        read_table(input_path),
        Path(args.output),
        final_batch=args.locked_final_batch,
        source_sidecar=sidecar,
        diagnostic_standard=standard,
        label_column=args.label_column,
        compound_column=args.compound_column,
        batch_column=args.batch_column,
        row_id_column=args.row_id_column,
        input_path=input_path,
        source_sidecar_path=sidecar_path,
        diagnostic_standard_path=standard_path,
        seed=args.seed,
        specificity_target=args.specificity_target,
        bootstrap_iterations=args.bootstrap_iterations,
    )
    decision = artifact["trainable_decision"]
    print(
        json.dumps(
            {
                "model_id": decision["model_id"],
                "selected_model": decision["selected_model"],
                "output": str(Path(args.output).resolve()),
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    return 0


def apply_main(argv=None) -> int:
    """Append the five shadow fields while preserving every v2 field."""

    from .shadow import apply_shadow_csv

    parser = argparse.ArgumentParser(
        description="Apply an MS2 model in shadow mode; omit --model for feasibility-only output."
    )
    parser.add_argument("--input", required=True, help="Existing v2 result CSV")
    parser.add_argument("--output", required=True, help="Combined v2 plus ML shadow CSV")
    parser.add_argument("--model", help="trainable_decision JSON; omitted means UNTRAINED")
    parser.add_argument("--source-sidecar", help="Defaults to <input>.metadata.json")
    parser.add_argument("--standard", default=str(_default_standard()))
    args = parser.parse_args(argv)

    input_path = Path(args.input)
    sidecar_path = (
        Path(args.source_sidecar)
        if args.source_sidecar
        else Path(str(input_path) + ".metadata.json")
    )
    sidecar = _read_mapping(sidecar_path)
    standard = _read_mapping(Path(args.standard))
    validate_frozen_source(sidecar, standard)
    summary = apply_shadow_csv(
        input_path,
        Path(args.output),
        model_artifact_path=args.model,
        source_sidecar_path=sidecar_path,
        standard_path=args.standard,
    )
    if summary is None:
        summary = {
            "output": str(Path(args.output).resolve()),
            "mode": "trained_shadow" if args.model else "feasibility_only_untrained",
        }
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


def review_main(argv=None) -> int:
    """Build the v2-versus-ML disagreement table, HTML, PNG, and SVG."""

    parser = argparse.ArgumentParser(description="Export all v2/ML shadow disagreements.")
    parser.add_argument("--input", required=True, help="Combined v2 plus ML shadow CSV")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--model", help="Optional trainable_decision model artifact")
    parser.add_argument("--source-sidecar", help="Defaults to the source shadow CSV sidecar")
    parser.add_argument("--standard", default=str(_default_standard()))
    parser.add_argument("--method-profile")
    parser.add_argument("--formats", nargs="+", choices=("svg", "png"), default=("svg", "png"))
    parser.add_argument("--png-scale", type=float, default=2.0)
    args = parser.parse_args(argv)

    # Keep ``shadow-review --help`` available in minimal environments.  Pillow
    # is needed only after argument parsing, when PNG rendering is requested.
    from .shadow_review import build_shadow_review

    summary = build_shadow_review(
        args.input,
        args.output_dir,
        model_artifact_path=args.model,
        sidecar_path=args.source_sidecar,
        standard_path=args.standard,
        method_profile_path=args.method_profile,
        image_formats=tuple(args.formats),
        png_scale=args.png_scale,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


def package_review_main(argv=None) -> int:
    """Build a compact, integrity-checked ZIP for human shadow review."""

    parser = argparse.ArgumentParser(
        description=(
            "Package one final shadow CSV plus the canonical PNG/SVG review assets. "
            "Duplicated views, galleries, and implementation files are excluded."
        )
    )
    parser.add_argument("--input", required=True, help="Final combined v2 plus ML shadow CSV")
    parser.add_argument(
        "--review-dir",
        required=True,
        help="Directory previously produced by shadow-review",
    )
    parser.add_argument(
        "--output",
        help="ZIP path; defaults to share/<input-stem>_manual_review.zip",
    )
    parser.add_argument(
        "--model",
        help="Required for trained results; deeply validated but not copied into the ZIP",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace an existing output ZIP after all inputs pass validation",
    )
    args = parser.parse_args(argv)

    from .review_package import package_shadow_review

    summary = package_shadow_review(
        args.input,
        args.review_dir,
        args.output,
        model_artifact_path=args.model,
        overwrite=args.force,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


def _default_standard() -> Path:
    return Path(__file__).resolve().parents[2] / "standards" / "ms2_diagnostic_standard_v2.json"


def _read_mapping(path: Path) -> dict:
    with Path(path).open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, Mapping):
        raise ValueError(f"Expected a JSON object: {path}")
    return dict(value)


__all__ = ["apply_main", "package_review_main", "review_main", "train_main"]
