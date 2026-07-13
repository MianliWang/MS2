"""兼容层：旧代码如果 import MSAI.python.fragments，仍然可以继续使用。

新代码建议直接看：
- get_frag.py: 旧版单峰 GetFrag。
- get_chiral_frag.py: 新版手性双峰 GetChiralFrag。
- ms2/: 单峰/双峰共享的数据结构、读取、索引、EIC 和 DIA helper。
"""

from __future__ import annotations

try:
    from .chiral_similarity import (
        ChiralPairThresholds,
        classify_chiral_pair,
        compare_fragment_spectra,
    )
    from .get_chiral_frag import (
        GetChiralFrag,
        analyze_chiral_peak_pairs,
        extract_fragments_for_rt_window,
        get_chiral_frag,
    )
    from .get_frag import GetFrag, get_frag
    from .ms2 import (
        Eic,
        Spectrum,
        format_fragment_string,
        mz_bounds,
        normalize_rt_window,
        summarize_xic_peak,
    )
except ImportError:
    from chiral_similarity import (  # type: ignore[import-not-found]
        ChiralPairThresholds,
        classify_chiral_pair,
        compare_fragment_spectra,
    )
    from get_chiral_frag import (  # type: ignore[import-not-found]
        GetChiralFrag,
        analyze_chiral_peak_pairs,
        extract_fragments_for_rt_window,
        get_chiral_frag,
    )
    from get_frag import GetFrag, get_frag  # type: ignore[import-not-found]
    from ms2 import (  # type: ignore[import-not-found]
        Eic,
        Spectrum,
        format_fragment_string,
        mz_bounds,
        normalize_rt_window,
        summarize_xic_peak,
    )

__all__ = [
    "ChiralPairThresholds",
    "Eic",
    "GetChiralFrag",
    "GetFrag",
    "Spectrum",
    "analyze_chiral_peak_pairs",
    "classify_chiral_pair",
    "compare_fragment_spectra",
    "extract_fragments_for_rt_window",
    "format_fragment_string",
    "get_chiral_frag",
    "get_frag",
    "mz_bounds",
    "normalize_rt_window",
    "summarize_xic_peak",
]


def _demo():
    """保留原来的 `python MSAI/python/fragments.py` 自检入口。"""

    assert normalize_rt_window((1, 2), "min") == (60, 120)
    assert normalize_rt_window((10, 20), "sec") == (10, 20)
    assert mz_bounds(100, 10, "ppm") == (99.999, 100.001)
    assert mz_bounds(100, 0.001, "Da") == (99.999, 100.001)
    assert mz_bounds(100, 0.00001, "legacy_fraction") == (99.999, 100.001)

    eic = Eic(rt=[0, 1, 2], scan=[0, 1, 2], intensity=[0, 10, 0])
    summary = summarize_xic_peak(eic)
    assert summary["apex_rt"] == 1
    assert summary["apex_intensity"] == 10
    assert summary["area"] == 10
    assert summary["quality_flags"] == "ok"
    assert format_fragment_string([], []) == ""
    assert format_fragment_string([50, 75], [1000, 2000]) == "50,1000;75,2000"

    spectra = [
        Spectrum(precursor_mz=300, rt=0, mz=[100, 290, 295], intensity=[0, 0, 0]),
        Spectrum(precursor_mz=300, rt=1, mz=[100, 290, 295], intensity=[5000, 10, 100]),
        Spectrum(precursor_mz=300, rt=2, mz=[100, 290, 295], intensity=[0, 0, 0]),
    ]
    result = extract_fragments_for_rt_window(spectra, 295, 0.01, "Da", 300, (0, 2))
    assert result["MS2"] == "100,5000"


if __name__ == "__main__":
    _demo()
    print("ok")
