"""MSAI fragment extraction 的 MS2 核心工具。

这个文件只放单峰和双峰都会用到的部分：
- 简化后的 MS2 数据结构；
- mzML/mzXML 读取；
- DIA window 和 EIC/XIC 提取；
- m/z tolerance、相关性、CSV 读写等小工具。
"""

from __future__ import annotations

import csv
import importlib
import math
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, cast


@dataclass
class Spectrum:
    """一个 MS2 scan 的最小表示。

    R 代码直接操作 xcmsRaw 的 S4 slots；Python 版不复制那个复杂对象，
    只保留本算法真正需要的字段：
    - precursor_mz: 当前 MS2 scan 对应的 DIA precursor/window m/z。
    - rt: scan retention time，单位是秒。
    - mz/intensity: 这个 scan 里的峰表。
    """

    precursor_mz: float
    rt: float
    mz: list[float]
    intensity: list[float]


@dataclass
class DiaData:
    """某一个 DIA window 下的 MS2 scans。

    这对应 R 代码里的 ms2copy() 返回值：把同一个 precursor window 的
    MS2 scan 抽出来，后续在 fragment m/z 维度上做 EIC。
    """

    spectra: list[Spectrum]
    mzrange: tuple[float, float]


@dataclass
class Eic:
    """Extracted Ion Chromatogram / XIC。

    每个点对应一个 scan：
    - rt: scan 时间，秒。
    - scan: 在当前 DiaData.spectra 里的下标。
    - intensity: 指定 mz window 内所有 peak intensity 的和。
    """

    rt: list[float]
    scan: list[int]
    intensity: list[float]


def load_ms2_spectra(path: Path) -> list[Spectrum]:
    """读取 mzML/mzXML 并只保留 MS2 spectra。

    pyopenms 是可选依赖，所以这里动态导入：
    - 没有安装 pyopenms 时，普通 helper 自检仍能跑；
    - 只有真正读 raw data 时才报清晰错误。
    """

    try:
        oms = cast(Any, importlib.import_module("pyopenms"))
    except ImportError as exc:
        raise RuntimeError("Raw mzML/mzXML loading needs pyopenms: python -m pip install pyopenms") from exc

    exp = oms.MSExperiment()
    suffix = path.suffix.lower()
    if suffix == ".mzml":
        oms.MzMLFile().load(str(path), exp)
    elif suffix == ".mzxml":
        oms.MzXMLFile().load(str(path), exp)
    else:
        raise ValueError(f"Unsupported raw file type: {path.name}")

    spectra: list[Spectrum] = []
    for spec in exp:
        # 本算法只使用 MS2；MS1 spectra 对 fragment extraction 没有参与。
        if spec.getMSLevel() != 2:
            continue
        precursors = spec.getPrecursors()
        if not precursors:
            continue
        mz, intensity = spec.get_peaks()
        spectra.append(
            Spectrum(
                precursor_mz=float(precursors[0].getMZ()),
                rt=float(spec.getRT()),
                mz=[float(x) for x in mz],
                intensity=[float(x) for x in intensity],
            )
        )
    return spectra


def preclist(spectra: Iterable[Spectrum]) -> list[float]:
    """提取出现过的 DIA precursor/window m/z，保留原始出现顺序。"""

    seen: list[float] = []
    for spec in spectra:
        if spec.precursor_mz not in seen:
            seen.append(spec.precursor_mz)
    return seen


def normalize_rt_window(rt_window, rt_unit: str) -> tuple[float, float]:
    """校验 RT window，并统一转换成秒。"""

    if rt_unit not in {"min", "sec"}:
        raise ValueError("rt_unit must be either 'min' or 'sec'.")
    if len(rt_window) != 2:
        raise ValueError("rt_window must contain exactly two values.")
    start, end = (float(rt_window[0]), float(rt_window[1]))
    if not (math.isfinite(start) and math.isfinite(end)):
        raise ValueError("rt_window must contain finite values.")
    if start >= end:
        raise ValueError("rt_window start must be smaller than rt_window end.")
    if rt_unit == "min":
        return (start * 60, end * 60)
    return (start, end)


def mz_bounds(mz: float, mz_tol: float, mz_tol_unit: str) -> tuple[float, float]:
    """根据 tolerance 单位计算 mz 下界和上界。

    legacy_fraction 保留旧 R 代码行为：mz +/- mz * mz_tol。
    ppm 和 Da 是 GetChiralFrag.R 里新增的更明确单位。
    """

    mz = float(mz)
    mz_tol = float(mz_tol)
    if not (math.isfinite(mz) and math.isfinite(mz_tol)) or mz_tol < 0:
        raise ValueError("mz and mz_tol must be finite, and mz_tol must be non-negative.")
    if mz_tol_unit == "legacy_fraction":
        delta = mz * mz_tol
    elif mz_tol_unit == "ppm":
        delta = mz * mz_tol / 1_000_000
    elif mz_tol_unit == "Da":
        delta = mz_tol
    else:
        raise ValueError("mz_tol_unit must be one of 'legacy_fraction', 'ppm', or 'Da'.")
    return (mz - delta, mz + delta)


def summarize_xic_peak(eic: Eic):
    """总结一条 EIC/XIC 的峰顶和面积。"""

    if not eic.intensity:
        return _empty_xic_summary("empty_xic")
    apex_index = max(
        range(len(eic.intensity)),
        key=lambda i: eic.intensity[i] if math.isfinite(eic.intensity[i]) else -math.inf,
    )
    valid = [i for i, (rt, inten) in enumerate(zip(eic.rt, eic.intensity)) if math.isfinite(rt) and math.isfinite(inten)]
    if not valid:
        return _empty_xic_summary("missing_rt")
    if len(valid) < 2:
        area = None
        flags = "insufficient_points"
    else:
        pairs = sorted((eic.rt[i], eic.intensity[i]) for i in valid)
        area = sum((b[0] - a[0]) * (a[1] + b[1]) / 2 for a, b in zip(pairs, pairs[1:]))
        flags = "ok"
    return {
        "apex_rt": eic.rt[apex_index],
        "apex_intensity": eic.intensity[apex_index],
        "area": area,
        "quality_flags": flags,
    }


def format_fragment_string(fragment_mz, intensity) -> str:
    """把 fragment mz 和强度拼成旧 CSV cell 格式。"""

    pairs = [(mz, inten) for mz, inten in zip(fragment_mz, intensity) if mz is not None and inten is not None]
    if not pairs:
        return ""
    return ";".join(f"{mz:g},{inten:g}" for mz, inten in pairs)


def extract_window_result(
    spectra: list[Spectrum],
    precurmz: float,
    mz_tol: float,
    mz_tol_unit: str,
    diawin: float,
    rt_window_sec: tuple[float, float],
):
    """核心提取逻辑：一个 feature + 一个 RT window -> 一个结果 dict。

    单峰和双峰入口都调用这里，区别只在 RT window 来自哪里：
    - 单峰：peaklist rt +/- 10 秒；
    - 双峰：用户传入 peak_a / peak_b 两个 RT window。
    """

    # 只保留目标 DIA window 的 MS2 scans，相当于 R 里的 ms2copy()。
    dia = _ms2copy(spectra, diawin)
    if not dia.spectra:
        return empty_window_result("no_ms2_scans")

    # precursor EIC 的 m/z window，同时裁剪到该 DIA window 实际 mz 范围内。
    mzmin, mzmax = _clamped_bounds(precurmz, mz_tol, mz_tol_unit, dia.mzrange)
    if mzmin >= mzmax:
        return empty_window_result("mz_out_of_range")

    # 用户给的 RT window 也裁剪到 raw data 实际覆盖的 RT 范围内。
    rtmin = max(min(s.rt for s in dia.spectra), rt_window_sec[0])
    rtmax = min(max(s.rt for s in dia.spectra), rt_window_sec[1])
    if rtmin >= rtmax:
        return empty_window_result("rt_window_out_of_range")

    # ms2_count 用于告诉用户这个窗口内实际有多少 MS2 scans。
    scan_count = sum(rtmin <= s.rt <= rtmax for s in dia.spectra)
    if scan_count < 1:
        return empty_window_result("no_ms2_scans", 0)

    # native EIC 是 precursor m/z 在当前 RT window 内的 chromatogram。
    native = _raw_eic(dia, mzmin, mzmax, rtmin, rtmax)
    peak_summary = summarize_xic_peak(native)
    flags = _split_flags(peak_summary["quality_flags"])
    if not native.scan or not native.intensity:
        return merge_window_result(peak_summary, "", scan_count, flags + ["no_fragment_scan"])

    # 找 native EIC 的峰顶 scan，然后从这个 scan 的 MS2 peak list 里找候选 fragment。
    apex_pos = max(range(len(native.intensity)), key=lambda i: native.intensity[i])
    apex_scan = native.scan[apex_pos]
    if apex_scan >= len(dia.spectra) - 1:
        return merge_window_result(peak_summary, "", scan_count, flags + ["invalid_apex_scan"])

    apex_spec = dia.spectra[apex_scan]
    # legacy 规则：fragment m/z 必须比 precursor m/z 至少小 10 Da。
    candidates = [mz for mz in apex_spec.mz if mz < precurmz - 10]
    if not candidates:
        return merge_window_result(peak_summary, "", scan_count, flags + ["no_candidate_fragments"])

    fragment_mz: list[float] = []
    fragment_intensity: list[float] = []
    native_sd = _sd(native.intensity)
    for mz0 in candidates:
        # 对每个 candidate 用同样的 mz tolerance 提取 fragment EIC。
        frag_mzmin, frag_mzmax = _clamped_bounds(mz0, mz_tol, mz_tol_unit, dia.mzrange)
        if frag_mzmin >= frag_mzmax:
            continue
        frag = _raw_eic(dia, frag_mzmin, frag_mzmax, rtmin, rtmax).intensity
        if len(frag) != len(native.intensity):
            continue

        # flat trace 的相关性没有意义；R 版用 sd == 0 跳过。
        frag_sd = _sd(frag)
        if native_sd == 0 or frag_sd == 0:
            continue
        frag_max = max(frag)

        # legacy intensity threshold。
        if frag_max < 2000:
            continue

        # 共洗脱判定：fragment EIC 和 native EIC 的 Pearson correlation > 0.9。
        if _cor(native.intensity, frag) > 0.9:
            fragment_mz.append(mz0)
            fragment_intensity.append(frag_max)

    ms2 = format_fragment_string(fragment_mz, fragment_intensity)
    if not ms2:
        flags.append("no_fragments")
    return merge_window_result(peak_summary, ms2, scan_count, flags)


def _ms2copy(spectra: list[Spectrum], diawin: float) -> DiaData:
    """复制/筛选某个 DIA precursor window 的 MS2 scans。"""

    picked = [spec for spec in spectra if spec.precursor_mz == diawin]
    all_mz = [mz for spec in picked for mz in spec.mz]
    if not picked or not all_mz:
        return DiaData(picked, (math.inf, -math.inf))
    return DiaData(picked, (min(all_mz), max(all_mz)))


def _raw_eic(dia: DiaData, mzmin: float, mzmax: float, rtmin: float, rtmax: float) -> Eic:
    """在给定 mz/RT window 内提取 EIC。"""

    rt: list[float] = []
    scan: list[int] = []
    intensity: list[float] = []
    for i, spec in enumerate(dia.spectra):
        if not (rtmin <= spec.rt <= rtmax):
            continue
        rt.append(spec.rt)
        scan.append(i)
        intensity.append(sum(y for x, y in zip(spec.mz, spec.intensity) if mzmin <= x <= mzmax))
    return Eic(rt=rt, scan=scan, intensity=intensity)


def _clamped_bounds(mz: float, mz_tol: float, mz_tol_unit: str, mzrange: tuple[float, float]) -> tuple[float, float]:
    """计算 mz tolerance window，并裁剪到实际 mzrange 内。"""

    low, high = mz_bounds(mz, mz_tol, mz_tol_unit)
    return max(mzrange[0], low), min(mzrange[1], high)


def workspace_paths(base: Path):
    """生成旧项目约定的三个工作目录路径。"""

    return {
        "data": base / "data",
        "peak": base / "peaklist",
        "results": base / "results",
    }


def peak_files(path: Path) -> list[Path]:
    """列出支持的 peaklist 文件。"""

    if not path.exists():
        return []
    return sorted(p for p in path.iterdir() if p.is_file() and p.suffix.lower() in {".csv", ".xlsx"})


def raw_files(path: Path) -> list[Path]:
    """列出支持的 raw MS 文件。"""

    if not path.exists():
        return []
    return sorted(p for p in path.iterdir() if p.is_file() and p.suffix.lower() in {".mzml", ".mzxml"})


def matching_raw_file(peakfile: Path, msfiles: list[Path]) -> Path | None:
    """用 peaklist 文件名匹配 raw data。

    第一优先级：peaklist stem 是 raw file stem 的子串。
    第二优先级：从文件名提取 A01/A11 这类 plate/well token，再匹配 raw file。
    """

    stem = peakfile.stem
    matches = [p for p in msfiles if stem in p.stem]
    if not matches:
        tokens = re.findall(r"[A-Za-z]\d{2,}", stem)
        matches = [p for p in msfiles if any(token.lower() in p.stem.lower() for token in tokens)]
    if not matches:
        return None
    if len(matches) > 1:
        warnings.warn(f"Multiple raw MS data files matched peaklist {peakfile.name}; using first match: {matches[0].name}")
    return matches[0]


def read_table(path: Path) -> list[dict[str, str]]:
    """读取 CSV/XLSX peaklist，并统一成 list[dict[str, str]]。"""

    if path.suffix.lower() == ".csv":
        with path.open(newline="", encoding="utf-8-sig") as handle:
            return [_string_row(row) for row in csv.DictReader(handle)]
    if path.suffix.lower() == ".xlsx":
        try:
            openpyxl = cast(Any, importlib.import_module("openpyxl"))
        except ImportError as exc:
            raise RuntimeError("XLSX peaklists need openpyxl: python -m pip install openpyxl") from exc
        workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
        sheet = workbook.active
        if sheet is None:
            return []
        rows = sheet.iter_rows(values_only=True)
        headers = ["" if value is None else str(value) for value in next(rows)]
        return [dict(zip(headers, ["" if value is None else str(value) for value in row])) for row in rows]
    raise ValueError(f"Unsupported peaklist type: {path.name}")


def _string_row(row: dict[str | None, str | None]) -> dict[str, str]:
    """把 csv.DictReader 的一行清洗成普通字符串 dict。"""

    return {str(key): "" if value is None else str(value) for key, value in row.items() if key is not None}


def write_table(path: Path, rows: list[dict]):
    """写 CSV，并保留 rows 中出现过的所有列。"""

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
    """大小写不敏感地查找列名，例如 mz 可以匹配 MZ。"""

    wanted = name.lower()
    for key in rows[0].keys():
        if key.lower() == wanted:
            return key
    return None


def number(value) -> float | None:
    """把单元格值转成有限 float；失败时返回 None。"""

    try:
        number_value = float(value)
    except (TypeError, ValueError):
        return None
    return number_value if math.isfinite(number_value) else None


def first_dia_window(mz: float, precursor: list[float], dia_iso_win: float) -> float | None:
    """找到第一个距离 feature mz 不超过半个 DIA window 宽度的 precursor。"""

    for value in precursor:
        if abs(mz - value) <= 0.5 * dia_iso_win:
            return value
    return None


def csv_scalar(value) -> str:
    """把 None/NaN/float 转成适合 CSV cell 的字符串。"""

    if value is None:
        return ""
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    return f"{value:g}" if isinstance(value, float) else str(value)


def _empty_xic_summary(flag: str):
    """EIC 层失败或为空时的统一摘要结构。"""

    return {
        "apex_rt": None,
        "apex_intensity": None,
        "area": None,
        "quality_flags": flag,
    }


def empty_window_result(flag: str, ms2_count=None):
    """整个 feature/window 无法提取时的统一结果结构。"""

    return {
        "MS2": "",
        "apex_rt": None,
        "apex_intensity": None,
        "area": None,
        "ms2_count": ms2_count,
        "quality_flags": flag,
    }


def merge_window_result(peak_summary: dict, ms2: str, ms2_count: int, flags: list[str]):
    """合并 XIC summary、fragment 字符串和质量标记。"""

    clean_flags = []
    for flag in flags:
        if flag and flag != "ok" and flag not in clean_flags:
            clean_flags.append(flag)
    if not clean_flags:
        clean_flags = ["ok"]
    return {
        "MS2": ms2,
        "apex_rt": peak_summary["apex_rt"],
        "apex_intensity": peak_summary["apex_intensity"],
        "area": peak_summary["area"],
        "ms2_count": ms2_count,
        "quality_flags": ";".join(clean_flags),
    }


def _split_flags(flags: str) -> list[str]:
    """把 CSV 里的分号 flag 字符串拆回列表。"""

    if not flags or flags == "ok":
        return []
    return flags.split(";")


def _sd(values: list[float]) -> float:
    """样本标准差；长度不足或全无效时返回 0。"""

    clean = [x for x in values if math.isfinite(x)]
    if len(clean) < 2:
        return 0.0
    mean = sum(clean) / len(clean)
    return math.sqrt(sum((x - mean) ** 2 for x in clean) / (len(clean) - 1))


def _cor(left: list[float], right: list[float]) -> float:
    """Pearson correlation，忽略非有限值。"""

    pairs = [(x, y) for x, y in zip(left, right) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 2:
        return 0.0
    xs, ys = zip(*pairs)
    xmean = sum(xs) / len(xs)
    ymean = sum(ys) / len(ys)
    xdiff = [x - xmean for x in xs]
    ydiff = [y - ymean for y in ys]
    denom = math.sqrt(sum(x * x for x in xdiff) * sum(y * y for y in ydiff))
    return 0.0 if denom == 0 else sum(x * y for x, y in zip(xdiff, ydiff)) / denom
