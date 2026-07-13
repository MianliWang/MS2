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
    from .ms2 import (
        build_ms2_index,
        column,
        extract_window_result,
        first_dia_window,
        load_ms2_spectra,
        matching_raw_file,
        number,
        peak_files,
        raw_files,
        read_table,
        workspace_paths,
        write_table,
    )
except ImportError:
    from ms2 import (  # type: ignore[import-not-found]
        build_ms2_index,
        column,
        extract_window_result,
        first_dia_window,
        load_ms2_spectra,
        matching_raw_file,
        number,
        peak_files,
        raw_files,
        read_table,
        workspace_paths,
        write_table,
    )


def GetFrag(mz_tol: float, DIAisowin: float, RTwin=None, base_path=None):
    """旧式大写名称兼容入口；实际逻辑见 :func:`get_frag`。"""

    return get_frag(mz_tol, DIAisowin, RTwin, base_path)


def get_frag(mz_tol: float, dia_iso_win: float, rt_win=None, base_path=None):
    """按旧目录布局执行单RT中心点fragment extraction。

    每行使用``mz + rt``，RT按分钟解释并固定提取``±10秒``。``rt_win``仅为
    历史签名兼容，实际上不参与计算。当前Peak1/Peak2正式流程不应调用此函数。
    """

    base = Path(base_path or ".").resolve()
    paths = workspace_paths(base)

    peakfiles = peak_files(paths["peak"])
    if not peakfiles:
        warnings.warn(f"No peaklist files found in {paths['peak']}; returning None.", stacklevel=2)
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
            warnings.warn(f"{peakfile.name} has no mz/MZ column; skipping.", stacklevel=2)
            continue
        if rt_col is None:
            warnings.warn(f"{peakfile.name} has no rt column; skipping.", stacklevel=2)
            continue

        for row in rows:
            row["MS2"] = "0"

        rawfile = matching_raw_file(peakfile, msfiles)
        if rawfile is None:
            warnings.warn(
                f"No matching raw MS data found for peaklist {peakfile.name}; skipping.",
                stacklevel=2,
            )
            continue

        spectra = load_ms2_spectra(rawfile)
        spectra_index = build_ms2_index(spectra)
        for row in rows:
            mz = number(row.get(mz_col))
            rt = number(row.get(rt_col))
            if mz is None or rt is None:
                continue

            diawin = first_dia_window(mz, spectra_index, dia_iso_win)
            if diawin is None:
                continue

            # legacy 行为：peaklist rt 单位是分钟，转换成秒后取 +/- 10 秒。
            fragments = _extract_fragments(
                spectra=spectra_index,
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
    """把新单窗口结果压缩为旧GetFrag只需要的fragment字符串。

    ``legacy_last_scan_none``保留历史边界行为；它不应传播到新分析接口。
    """

    result = extract_window_result(spectra, precurmz, mz_tol, mz_tol_unit, diawin, rt_window_sec)
    if legacy_last_scan_none and "invalid_apex_scan" in result["quality_flags"]:
        return None
    return result["MS2"]
