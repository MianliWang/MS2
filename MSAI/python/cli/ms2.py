"""Single command router for the supported MS2 workflow.

Examples::

    python -m MSAI.python.cli.ms2 analyze --help
    python -m MSAI.python.cli.ms2 report --help
    python -m MSAI.python.cli.ms2 export-review --help
    python -m MSAI.python.cli.ms2 train-shadow --help
    python -m MSAI.python.cli.ms2 apply-shadow --help
    python -m MSAI.python.cli.ms2 shadow-review --help
    python -m MSAI.python.cli.ms2 package-shadow-review --help
"""

from __future__ import annotations

import sys


def analyze_main(argv):
    """延迟加载并执行核心Peak1/Peak2 MS2分析命令。"""

    from ..ms2.pipeline import main

    return main(argv)


def report_main(argv):
    """延迟加载并执行紧凑HTML/CSV技术报告命令。"""

    from ..ms2_review.report import main

    return main(argv)


def export_review_main(argv):
    """延迟加载并执行逐目标SVG/PNG人工复核导出命令。"""

    from ..export_ms2_spectrum_review import main

    return main(argv)


def train_shadow_main(argv):
    """延迟加载独立truth、分组交叉验证的shadow训练命令。"""

    from ..ms2_ml.cli import train_main

    return train_main(argv)


def apply_shadow_main(argv):
    """延迟加载只追加五个ML字段的shadow推理命令。"""

    from ..ms2_ml.cli import apply_main

    return apply_main(argv)


def shadow_review_main(argv):
    """延迟加载v2/ML分歧CSV、HTML和谱图导出命令。"""

    from ..ms2_ml.cli import review_main

    return review_main(argv)


def package_shadow_review_main(argv):
    """延迟加载可重复、完整性校验的人工复核ZIP打包命令。"""

    from ..ms2_ml.cli import package_review_main

    return package_review_main(argv)


COMMANDS = {
    "analyze": analyze_main,
    "report": report_main,
    "export-review": export_review_main,
    "train-shadow": train_shadow_main,
    "apply-shadow": apply_shadow_main,
    "shadow-review": shadow_review_main,
    "package-shadow-review": package_shadow_review_main,
}


def main(argv=None):
    """统一路由MS2分析、人工复核和ML shadow子命令。"""

    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in {"-h", "--help"}:
        names = "\n".join(
            (
                "  analyze        extract and compare Peak1/Peak2 DIA-MS2 spectra",
                "  report         build the compact HTML/CSV technical report",
                "  export-review  build per-target SVG/PNG review folders",
                "  train-shadow   fit fixed threshold/logistic baselines on grouped truth",
                "  apply-shadow   append five ML shadow fields without changing v2",
                "  shadow-review  export all v2/ML disagreements as CSV/HTML/PNG/SVG",
                "  package-shadow-review  package final CSV and canonical review images",
            )
        )
        print(f"Usage: python -m MSAI.python.cli.ms2 <command> [options]\n\nCommands:\n{names}")
        return 0
    command = argv.pop(0)
    handler = COMMANDS.get(command)
    if handler is None:
        choices = ", ".join(COMMANDS)
        raise SystemExit(f"Unknown MS2 command {command!r}; choose one of: {choices}")
    return handler(argv)


if __name__ == "__main__":
    main()
