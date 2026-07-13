"""Small data structures shared by the MS2 pipeline."""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass
class Spectrum:
    """一张质心化MS2扫描的最小内存表示。

    ``precursor_mz``通常是DIA隔离窗口中心而不一定是某个真实前体峰；
    ``rt``统一使用秒。``mz``与``intensity``必须等长并保持一一对应。
    隔离上下界优先来自原始文件，缺失时才允许上层使用名义窗口宽度回退。
    """

    precursor_mz: float
    rt: float
    mz: list[float]
    intensity: list[float]
    isolation_lower_mz: float | None = None
    isolation_upper_mz: float | None = None


@dataclass
class DiaData:
    """同一DIA隔离窗口内按保留时间排序的一组MS2扫描。

    ``mzrange``描述这些扫描实际出现过的fragment m/z范围；``rts``是为了
    使用二分查找快速截取RT窗口而保存的并行索引。
    """

    spectra: list[Spectrum]
    mzrange: tuple[float, float]
    rts: list[float]
    isolation_lower_mz: float | None = None
    isolation_upper_mz: float | None = None


@dataclass
class Ms2Index:
    """每个原始文件只建立一次的DIA窗口/保留时间索引。

    ``windows``以隔离窗口中心为键，``precursors``保留窗口首次出现顺序。
    参数扫描或批量处理多行目标时复用这个对象，可避免反复解析原始文件。
    """

    windows: dict[float, DiaData]
    precursors: list[float]

    def get(self, precursor_mz: float) -> DiaData:
        """返回指定DIA中心的数据；不存在时返回结构稳定的空窗口。"""

        return self.windows.get(precursor_mz, DiaData([], (math.inf, -math.inf), []))

    def matching_windows(self, mz: float, fallback_width: float) -> list[float]:
        """返回覆盖目标m/z的窗口中心，按离目标的距离由近到远排序。

        原始文件给出隔离上下界时以真实边界为准；只有边界缺失时才用
        ``center ± fallback_width/2``。因此fallback宽度不会覆盖真实采集信息。
        """

        matches: list[float] = []
        for center in self.precursors:
            window = self.windows[center]
            lower = window.isolation_lower_mz
            upper = window.isolation_upper_mz
            if lower is None or upper is None:
                lower, upper = center - fallback_width / 2, center + fallback_width / 2
            if lower <= mz <= upper:
                matches.append(center)
        return sorted(matches, key=lambda center: (abs(center - mz), center))


@dataclass
class Eic:
    """提取离子流图（EIC/XIC）。

    ``rt``单位为秒；``scan``是所属DIA窗口内部的扫描下标；``intensity``
    是对应m/z容差内的求和强度。三个数组应具有相同长度。
    """

    rt: list[float]
    scan: list[int]
    intensity: list[float]


__all__ = ["DiaData", "Eic", "Ms2Index", "Spectrum"]
