"""CLI for standalone MS2 SVG/PNG evidence images and review-folder views."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

try:
    from .ms2_review import export_ms2_review
except ImportError:
    from ms2_review import export_ms2_review  # type: ignore[import-not-found]


def main(argv=None):
    """解析独立MS2图像导出参数并打印生成摘要JSON。"""

    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--standard", default=str(root / "standards" / "ms2_diagnostic_standard_v1.json")
    )
    parser.add_argument("--sidecar", help="Defaults to <input>.metadata.json")
    parser.add_argument(
        "--method-profile",
        help="Optional reference-method JSON; displayed separately from raw-observed metadata.",
    )
    parser.add_argument("--source-machine-label-column", default="IG")
    parser.add_argument("--formats", nargs="+", choices=("svg", "png"), default=("svg", "png"))
    parser.add_argument("--png-scale", type=float, default=2.0)
    args = parser.parse_args(argv)
    summary = export_ms2_review(
        args.input,
        args.output_dir,
        standard_path=args.standard,
        sidecar_path=args.sidecar,
        method_profile_path=args.method_profile,
        source_machine_label_column=args.source_machine_label_column,
        image_formats=tuple(args.formats),
        png_scale=args.png_scale,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
