"""旧版单峰 fragment extraction。

这个文件只放 GetFrag.R 对应的逻辑：
- peaklist 里每个 feature 使用 mz + rt；
- rt 按分钟解释，并使用 rt +/- 10 秒；
- 输出 results/{sample}_MS2.csv。
"""

from __future__ import annotations

import warnings
from pathlib import Path

try:
    from .ms2_core import (
        column,
        extract_window_result,
        first_dia_window,
        load_ms2_spectra,
        matching_raw_file,
        number,
        peak_files,
        preclist,
        raw_files,
        read_table,
        workspace_paths,
        write_table,
    )
except ImportError:
    from ms2_core import (  # type: ignore
        column,
        extract_window_result,
        first_dia_window,
        load_ms2_spectra,
        matching_raw_file,
        number,
        peak_files,
        preclist,
        raw_files,
        read_table,
        workspace_paths,
        write_table,
    )


def GetFrag(mz_tol: float, DIAisowin: float, RTwin=None, base_path=None):
    """兼容 R 函数名的入口；实际实现见 get_frag()."""

    return get_frag(mz_tol, DIAisowin, RTwin, base_path)


def get_frag(mz_tol: float, dia_iso_win: float, rt_win=None, base_path=None):
    """旧版单 RT 中心点 fragment extraction。

    rt_win 参数只保留接口兼容；旧 R 核心逻辑里也没有真正使用它。
    """

    base = Path(base_path or ".").resolve()
    paths = workspace_paths(base)

    peakfiles = peak_files(paths["peak"])
    if not peakfiles:
        warnings.warn(f"No peaklist files found in {paths['peak']}; returning None.")
        return None

    msfiles = raw_files(paths["data"])
    last_rows = None
    for k, peakfile in enumerate(peakfiles, start=1):
        print(f"getfragment... {k}")
        rows = read_table(peakfile)
        if not rows:
            continue

        # 兼容 R 版小写 mz，也兼容当前 CSV 里的大写 MZ。
        mz_col = column(rows, "mz")
        rt_col = column(rows, "rt")
        if mz_col is None:
            warnings.warn(f"{peakfile.name} has no mz/MZ column; skipping.")
            continue
        if rt_col is None:
            warnings.warn(f"{peakfile.name} has no rt column; skipping.")
            continue

        for row in rows:
            row["MS2"] = "0"

        rawfile = matching_raw_file(peakfile, msfiles)
        if rawfile is None:
            warnings.warn(f"No matching raw MS data found for peaklist {peakfile.name}; skipping.")
            continue

        spectra = load_ms2_spectra(rawfile)
        precursor = preclist(spectra)
        for row in rows:
            mz = number(row.get(mz_col))
            rt = number(row.get(rt_col))
            if mz is None or rt is None:
                continue

            diawin = first_dia_window(mz, precursor, dia_iso_win)
            if diawin is None:
                continue

            # legacy 行为：peaklist rt 单位是分钟，转换成秒后取 +/- 10 秒。
            fragments = _extract_fragments(
                spectra=spectra,
                precurmz=mz,
                mz_tol=mz_tol,
                mz_tol_unit="legacy_fraction",
                diawin=diawin,
                rt_window_sec=(rt * 60 - 10, rt * 60 + 10),
                legacy_last_scan_none=True,
            )
            if fragments is None:
                continue
            row["MS2"] = fragments

        sample = peakfile.stem
        write_table(paths["results"] / f"{sample}_MS2.csv", rows)
        last_rows = rows
    return last_rows


def _extract_fragments(
    spectra,
    precurmz: float,
    mz_tol: float,
    mz_tol_unit: str,
    diawin: float,
    rt_window_sec: tuple[float, float],
    legacy_last_scan_none: bool,
):
    """旧 GetFrag 使用的薄包装。

    单峰版只需要 MS2 fragment 字符串，不需要双峰版的 summary dict。
    legacy_last_scan_none 用来保留 R 版 apex scan 落在最后一个 scan 时返回 NULL 的行为。
    """

    result = extract_window_result(spectra, precurmz, mz_tol, mz_tol_unit, diawin, rt_window_sec)
    if legacy_last_scan_none and "invalid_apex_scan" in result["quality_flags"]:
        return None
    return result["MS2"]
