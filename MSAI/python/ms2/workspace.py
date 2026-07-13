"""Project path discovery and strict raw-file matching."""

from __future__ import annotations

import re
from pathlib import Path


def workspace_paths(base: Path) -> dict[str, Path]:
    """返回旧目录约定中的``data``、``peaklist``和``results``路径。

    新的显式CLI并不依赖该布局；此函数主要服务旧版``GetFrag``接口。
    """

    return {
        "data": base / "data",
        "peak": base / "peaklist",
        "results": base / "results",
    }


def peak_files(path: Path) -> list[Path]:
    """列出目录中可作为peaklist的CSV/XLSX文件，按路径稳定排序。"""

    if not path.exists():
        return []
    return sorted(
        candidate
        for candidate in path.iterdir()
        if candidate.is_file() and candidate.suffix.lower() in {".csv", ".xlsx"}
    )


def raw_files(path: Path) -> list[Path]:
    """列出目录中的mzML/mzXML原始文件，按路径稳定排序。"""

    if not path.exists():
        return []
    return sorted(
        candidate
        for candidate in path.iterdir()
        if candidate.is_file() and candidate.suffix.lower() in {".mzml", ".mzxml"}
    )


def matching_raw_file(peakfile: Path, msfiles: list[Path]) -> Path | None:
    """按文件stem或well token为一个peaklist匹配唯一原始文件。

    无匹配返回``None``；多个候选直接报错，避免把一个pool悄悄路由到错误
    raw文件。正式多文件运行更推荐使用显式manifest。
    """

    stem = peakfile.stem
    matches = [path for path in msfiles if stem in path.stem]
    if not matches:
        tokens = re.findall(r"[A-Za-z]\d{2,}", stem)
        matches = [
            path for path in msfiles if any(token.lower() in path.stem.lower() for token in tokens)
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
