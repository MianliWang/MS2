"""Single command router for the supported MS2 workflow.

Examples::

    python -m MSAI.python.cli.ms2 analyze --help
    python -m MSAI.python.cli.ms2 report --help
    python -m MSAI.python.cli.ms2 export-review --help
"""

from __future__ import annotations

import sys


def analyze_main(argv):
    from ..ms2.pipeline import main

    return main(argv)


def report_main(argv):
    from ..ms2_review.report import main

    return main(argv)


def export_review_main(argv):
    from ..export_ms2_spectrum_review import main

    return main(argv)


COMMANDS = {
    "analyze": analyze_main,
    "report": report_main,
    "export-review": export_review_main,
}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in {"-h", "--help"}:
        names = "\n".join(
            (
                "  analyze        extract and compare Peak1/Peak2 DIA-MS2 spectra",
                "  report         build the compact HTML/CSV technical report",
                "  export-review  build per-target SVG/PNG review folders",
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
