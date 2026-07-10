"""Project path discovery and strict raw-file matching."""

from __future__ import annotations

import re
from pathlib import Path


def workspace_paths(base: Path) -> dict[str, Path]:
    """Return the legacy project's data, peak-list, and result paths."""

    return {
        "data": base / "data",
        "peak": base / "peaklist",
        "results": base / "results",
    }


def peak_files(path: Path) -> list[Path]:
    if not path.exists():
        return []
    return sorted(
        candidate
        for candidate in path.iterdir()
        if candidate.is_file() and candidate.suffix.lower() in {".csv", ".xlsx"}
    )


def raw_files(path: Path) -> list[Path]:
    if not path.exists():
        return []
    return sorted(
        candidate
        for candidate in path.iterdir()
        if candidate.is_file() and candidate.suffix.lower() in {".mzml", ".mzxml"}
    )


def matching_raw_file(peakfile: Path, msfiles: list[Path]) -> Path | None:
    """Match one peak list to exactly one raw file by stem or well token."""

    stem = peakfile.stem
    matches = [path for path in msfiles if stem in path.stem]
    if not matches:
        tokens = re.findall(r"[A-Za-z]\d{2,}", stem)
        matches = [
            path
            for path in msfiles
            if any(token.lower() in path.stem.lower() for token in tokens)
        ]
    if not matches:
        return None
    if len(matches) > 1:
        names = ", ".join(path.name for path in matches)
        raise ValueError(
            f"Multiple raw MS data files matched peaklist {peakfile.name}: {names}. "
            "Use an explicit manifest or split the peaklist by sample/well."
        )
    return matches[0]


__all__ = ["matching_raw_file", "peak_files", "raw_files", "workspace_paths"]
