"""CSV/XLSX table I/O used by analysis and calibration commands."""

from __future__ import annotations

import csv
import importlib
from pathlib import Path
from typing import Any, cast


def read_table(path: Path) -> list[dict[str, str]]:
    """读取CSV/XLSX为字符串字典列表，仅删除完全空白记录。

    CSV使用``utf-8-sig``兼容Excel BOM；XLSX按活动工作表读取，并要求可选
    依赖``openpyxl``。这里不擅自删除重复行或解释业务字段。
    """

    path = Path(path)
    if path.suffix.lower() == ".csv":
        with path.open(newline="", encoding="utf-8-sig") as handle:
            return [
                clean for row in csv.DictReader(handle) if _row_has_data(clean := _string_row(row))
            ]
    if path.suffix.lower() == ".xlsx":
        try:
            openpyxl = cast(Any, importlib.import_module("openpyxl"))
        except ImportError as exc:
            raise RuntimeError(
                "XLSX peaklists need openpyxl: python -m pip install openpyxl"
            ) from exc
        workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
        sheet = workbook.active
        if sheet is None:
            return []
        values = sheet.iter_rows(values_only=True)
        try:
            header_row = next(values)
        except StopIteration:
            return []
        headers = ["" if value is None else str(value) for value in header_row]
        return [
            clean
            for row in values
            if _row_has_data(
                clean := dict(
                    zip(
                        headers,
                        ["" if value is None else str(value) for value in row],
                        strict=False,
                    )
                )
            )
        ]
    raise ValueError(f"Unsupported peaklist type: {path.name}")


def _string_row(row: dict[str | None, str | None]) -> dict[str, str]:
    """移除空表头并把``None``单元格规范为空字符串。"""

    return {
        str(key): "" if value is None else str(value)
        for key, value in row.items()
        if key is not None
    }


def _row_has_data(row: dict) -> bool:
    """判断一行是否至少包含一个非空白值。"""

    return any(str(value).strip() for value in row.values())


def write_table(path: Path, rows: list[dict]) -> None:
    """写出CSV，并按首次出现顺序保留所有行中出现过的字段。

    这样即使部分目标没有MS2结果，输出仍保持完整、稳定的列结构。
    """

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def column(rows: list[dict], name: str) -> str | None:
    """返回第一个不区分大小写的列名匹配；不存在时返回``None``。"""

    wanted = name.lower()
    for key in rows[0]:
        if key.lower() == wanted:
            return key
    return None


__all__ = ["column", "read_table", "write_table"]
