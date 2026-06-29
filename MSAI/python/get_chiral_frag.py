"""手性双峰 fragment extraction。

这个文件只放 GetChiralFrag.R 对应的新逻辑：
- 同一个 precursor m/z 在 peak_a / peak_b 两个 RT window 内分别提取；
- 每个窗口都有独立的 XIC summary 和 MS2 fragments；
- 输出 results/{sample}_chiral_MS2.csv。
"""

from __future__ import annotations

import warnings
from pathlib import Path

try:
    from .ms2_core import (
        Eic,
        Spectrum,
        column,
        csv_scalar,
        empty_window_result,
        extract_window_result,
        first_dia_window,
        format_fragment_string,
        load_ms2_spectra,
        matching_raw_file,
        mz_bounds,
        normalize_rt_window,
        number,
        peak_files,
        preclist,
        raw_files,
        read_table,
        summarize_xic_peak,
        workspace_paths,
        write_table,
    )
except ImportError:
    from ms2_core import (  # type: ignore
        Eic,
        Spectrum,
        column,
        csv_scalar,
        empty_window_result,
        extract_window_result,
        first_dia_window,
        format_fragment_string,
        load_ms2_spectra,
        matching_raw_file,
        mz_bounds,
        normalize_rt_window,
        number,
        peak_files,
        preclist,
        raw_files,
        read_table,
        summarize_xic_peak,
        workspace_paths,
        write_table,
    )


def GetChiralFrag(
    mz_tol: float,
    DIAisowin: float,
    peak_a_rt_window,
    peak_b_rt_window,
    rt_unit: str = "min",
    mz_tol_unit: str = "legacy_fraction",
    base_path=None,
):
    """兼容 R 函数名的入口；实际实现见 get_chiral_frag()."""

    return get_chiral_frag(
        mz_tol,
        DIAisowin,
        peak_a_rt_window,
        peak_b_rt_window,
        rt_unit,
        mz_tol_unit,
        base_path,
    )


def get_chiral_frag(
    mz_tol: float,
    dia_iso_win: float,
    peak_a_rt_window,
    peak_b_rt_window,
    rt_unit: str = "min",
    mz_tol_unit: str = "legacy_fraction",
    base_path=None,
):
    """手性双峰 fragment extraction。"""

    base = Path(base_path or ".").resolve()
    paths = workspace_paths(base)

    # pyopenms 里的 RT 是秒；用户常用分钟，所以入口处统一转换。
    rt_windows = {
        "peak_a": normalize_rt_window(peak_a_rt_window, rt_unit),
        "peak_b": normalize_rt_window(peak_b_rt_window, rt_unit),
    }

    # 用一个 dummy mz 提前校验 mz tolerance 参数是否合法。
    mz_bounds(1.0, mz_tol, mz_tol_unit)
    if rt_windows["peak_a"][1] > rt_windows["peak_b"][0] and rt_windows["peak_b"][1] > rt_windows["peak_a"][0]:
        warnings.warn("Chiral RT windows overlap; no deconvolution is performed.")

    peakfiles = peak_files(paths["peak"])
    if not peakfiles:
        warnings.warn(f"No peaklist files found in {paths['peak']}; returning None.")
        return None

    msfiles = raw_files(paths["data"])
    last_rows = None
    for k, peakfile in enumerate(peakfiles, start=1):
        print(f"getchiralfragment... {k}")
        rows = read_table(peakfile)
        if not rows:
            continue
        mz_col = column(rows, "mz")
        if mz_col is None:
            warnings.warn(f"{peakfile.name} has no mz/MZ column; skipping.")
            continue

        # 先补齐输出列，保证即使某些 feature 无结果，CSV 结构也稳定。
        _initialize_chiral_columns(rows)

        rawfile = matching_raw_file(peakfile, msfiles)
        if rawfile is None:
            warnings.warn(f"No matching raw MS data found for peaklist {peakfile.name}; skipping.")
            continue

        spectra = load_ms2_spectra(rawfile)
        precursor = preclist(spectra)
        for row in rows:
            mz = number(row.get(mz_col))
            if mz is None:
                continue
            diawin = first_dia_window(mz, precursor, dia_iso_win)
            if diawin is None:
                # 没有 DIA window 时，两个手性窗口都写同一个质量标记。
                for prefix in rt_windows:
                    _assign_chiral_result(row, prefix, empty_window_result("no_matching_DIA_window"))
                continue
            for prefix, rt_window in rt_windows.items():
                # 核心算法只处理“一个 feature + 一个 RT window”。
                # 手性双峰就是调用两次同一个 helper。
                result = extract_fragments_for_rt_window(
                    spectra=spectra,
                    precurmz=mz,
                    mz_tol=mz_tol,
                    mz_tol_unit=mz_tol_unit,
                    diawin=diawin,
                    rt_window_sec=rt_window,
                )
                _assign_chiral_result(row, prefix, result)

        sample = peakfile.stem
        write_table(paths["results"] / f"{sample}_chiral_MS2.csv", rows)
        last_rows = rows
    return last_rows


def extract_fragments_for_rt_window(
    spectra: list[Spectrum],
    precurmz: float,
    mz_tol: float,
    mz_tol_unit: str,
    diawin: float,
    rt_window_sec: tuple[float, float],
):
    """公开一点的单窗口 helper，便于测试或更细粒度调用。"""

    return extract_window_result(spectra, precurmz, mz_tol, mz_tol_unit, diawin, rt_window_sec)


def _initialize_chiral_columns(rows: list[dict]):
    """给每一行预先加 peak_a_* / peak_b_* 输出列。"""

    for row in rows:
        for prefix in ("peak_a", "peak_b"):
            row[f"{prefix}_MS2"] = ""
            row[f"{prefix}_apex_rt"] = ""
            row[f"{prefix}_apex_intensity"] = ""
            row[f"{prefix}_area"] = ""
            row[f"{prefix}_ms2_count"] = ""
            row[f"{prefix}_quality_flags"] = ""


def _assign_chiral_result(row: dict, prefix: str, result: dict):
    """把一个 window result dict 写回 CSV row。"""

    row[f"{prefix}_MS2"] = result["MS2"]
    row[f"{prefix}_apex_rt"] = csv_scalar(result["apex_rt"])
    row[f"{prefix}_apex_intensity"] = csv_scalar(result["apex_intensity"])
    row[f"{prefix}_area"] = csv_scalar(result["area"])
    row[f"{prefix}_ms2_count"] = csv_scalar(result["ms2_count"])
    row[f"{prefix}_quality_flags"] = result["quality_flags"]


def _demo():
    """最小自检：覆盖单位转换、mz tolerance、XIC 摘要和 fragment 筛选。"""

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
