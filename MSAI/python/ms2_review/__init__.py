"""Static, per-target MS2 evidence images and human-review folder views."""

from .acquisition import (
    load_method_profile,
    raw_acquisition_context,
    reconcile_acquisition,
)
from .classification import review_views, spectrum_availability
from .model import MirrorSpectrumData, prepare_mirror_spectrum
from .report import annotate_rows, build_report, read_result_rows


def export_ms2_review(*args, **kwargs):
    """仅在实际导出图像时延迟加载Pillow/SVG exporter。

    这样只运行核心MS2分析或读取报告数据时不必安装PNG相关依赖。
    """

    from .exporter import export_ms2_review as _export_ms2_review

    return _export_ms2_review(*args, **kwargs)


__all__ = [
    "MirrorSpectrumData",
    "annotate_rows",
    "build_report",
    "export_ms2_review",
    "load_method_profile",
    "prepare_mirror_spectrum",
    "raw_acquisition_context",
    "read_result_rows",
    "reconcile_acquisition",
    "review_views",
    "spectrum_availability",
]
