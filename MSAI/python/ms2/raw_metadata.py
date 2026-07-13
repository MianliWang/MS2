"""Lightweight scan-wide acquisition metadata aggregation for XML raw files."""

from __future__ import annotations

import mmap
import re
from collections import Counter
from pathlib import Path

_RAW_ATTRIBUTE_PATTERN = re.compile(
    rb'\b(msLevel|collisionEnergy|activationMethod|windowWideness)="([^"]+)"'
)


def raw_acquisition_context(path: str | Path | None) -> dict:
    """汇总整个raw XML中可直接观察到的关键采集属性。

    该函数使用内存映射扫描属性，不解码峰数组，因此适合在报告阶段快速统计
    MS1/MS2 scan数、碰撞能、activation method和window width。返回内容只
    代表原始文件记录；论文方法和实验室笔记必须在其他provenance层保存。
    """

    if not path:
        return {}
    raw_path = Path(path)
    if not raw_path.exists() or raw_path.stat().st_size == 0:
        return {}

    counts: dict[str, Counter[str]] = {
        "msLevel": Counter(),
        "collisionEnergy": Counter(),
        "activationMethod": Counter(),
        "windowWideness": Counter(),
    }
    with (
        raw_path.open("rb") as handle,
        mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ) as mapped,
    ):
        for match in _RAW_ATTRIBUTE_PATTERN.finditer(mapped):
            key = match.group(1).decode("ascii")
            value = match.group(2).decode("utf-8", errors="replace")
            counts[key][value] += 1

    observed_values = {
        key: [
            {"value": value, "count": count}
            for value, count in sorted(counter.items(), key=lambda item: item[0])
        ]
        for key, counter in counts.items()
        if counter
    }
    context: dict = {
        "source": "raw_xml",
        "scope": "entire_raw_file",
        "observed_values": observed_values,
    }
    for key in ("collisionEnergy", "activationMethod", "windowWideness"):
        if len(counts[key]) == 1:
            context[key] = next(iter(counts[key]))
    ms1_count = counts["msLevel"].get("1", 0)
    ms2_count = counts["msLevel"].get("2", 0)
    context["scan_counts"] = {"ms1": ms1_count, "ms2": ms2_count}
    if ms1_count:
        context["ms2_per_ms1"] = ms2_count / ms1_count
    return context


__all__ = ["raw_acquisition_context"]
