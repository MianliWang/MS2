"""Prepared mirror-spectrum data shared by the SVG and PNG renderers."""

from __future__ import annotations

from dataclasses import dataclass

try:
    from ..ms2.similarity import align_fragment_spectra
except ImportError:
    from ms2.similarity import align_fragment_spectra  # type: ignore[import-not-found]


@dataclass(frozen=True)
class MirrorSpectrumData:
    """SVG与PNG共用的镜像谱绘图数据。

    ``matches``及``matched_a/b``来自核心相似度模块的同一套过滤和匹配，
    ``max_a/max_b``用于让两侧分别按各自base peak归一化显示。
    """

    peak_a: tuple[tuple[float, float], ...]
    peak_b: tuple[tuple[float, float], ...]
    matches: tuple[tuple[int, int], ...]
    matched_a: frozenset[int]
    matched_b: frozenset[int]
    min_mz: float
    max_mz: float
    max_a: float
    max_b: float

    @property
    def availability(self) -> str:
        """返回双侧、仅Peak1、仅Peak2或双侧都无谱的可用性分类。"""

        if self.peak_a and self.peak_b:
            return "two_spectra"
        if self.peak_a:
            return "peak1_only"
        if self.peak_b:
            return "peak2_only"
        return "no_spectra"


def prepare_mirror_spectrum(
    row: dict,
    *,
    fragment_mz_tol: float = 0.01,
    fragment_mz_tol_unit: str = "Da",
    min_relative_intensity: float = 0.01,
) -> MirrorSpectrumData:
    """从一行结果准备镜像谱坐标、匹配集合和显示范围。

    这里调用``align_fragment_spectra``，确保图片中的匹配颜色与CSV相似度
    计算使用完全相同的m/z容差和相对强度过滤。
    """

    peak_a, peak_b, matches = align_fragment_spectra(
        row.get("peak_a_MS2"),
        row.get("peak_b_MS2"),
        fragment_mz_tol=fragment_mz_tol,
        fragment_mz_tol_unit=fragment_mz_tol_unit,
        min_relative_intensity=min_relative_intensity,
    )
    all_peaks = peak_a + peak_b
    minimum = min((mz for mz, _intensity in all_peaks), default=0.0)
    maximum = max((mz for mz, _intensity in all_peaks), default=1.0)
    if maximum - minimum < 1:
        maximum = minimum + 1
    return MirrorSpectrumData(
        tuple(peak_a),
        tuple(peak_b),
        tuple(matches),
        frozenset(index for index, _ in matches),
        frozenset(index for _, index in matches),
        minimum,
        maximum,
        max((intensity for _mz, intensity in peak_a), default=1.0),
        max((intensity for _mz, intensity in peak_b), default=1.0),
    )
