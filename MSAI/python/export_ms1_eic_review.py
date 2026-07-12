"""CLI facade for the modular MS1 EIC image-review exporter."""

from __future__ import annotations

import argparse
import json

try:
    from .ms1_review import export_ms1_review
except ImportError:
    from ms1_review import export_ms1_review  # type: ignore


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peaklist", required=True)
    parser.add_argument("--raw", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--baseline-config", required=True)
    parser.add_argument("--experimental-config", required=True)
    parser.add_argument("--candidate-config", required=True)
    parser.add_argument("--ms2-results")
    parser.add_argument("--mz-column", default="MZ")
    parser.add_argument(
        "--source-machine-label-column",
        default="IG",
        help="Column containing the source machine's colour/result label (default: IG).",
    )
    parser.add_argument("--rt-tolerance-min", type=float, default=0.1)
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("svg", "png"),
        default=("svg", "png"),
        help="Canonical image formats to export (default: svg png).",
    )
    parser.add_argument(
        "--png-scale",
        type=float,
        default=2.0,
        help="Raster scale relative to the 1240x930 SVG viewBox (default: 2).",
    )
    args = parser.parse_args(argv)
    summary = export_ms1_review(
        args.peaklist,
        args.raw,
        args.output_dir,
        baseline_config_path=args.baseline_config,
        experimental_config_path=args.experimental_config,
        candidate_config_path=args.candidate_config,
        ms2_result_path=args.ms2_results,
        mz_column=args.mz_column,
        source_machine_label_column=args.source_machine_label_column,
        rt_tolerance_min=args.rt_tolerance_min,
        image_formats=tuple(args.formats),
        png_scale=args.png_scale,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
