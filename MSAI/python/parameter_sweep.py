"""One-factor sensitivity analysis with one shared in-memory raw index."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

try:
    from .get_chiral_frag import analyze_chiral_peak_pairs
    from .ms2 import build_ms2_index, load_ms2_spectra, write_table
except ImportError:
    from get_chiral_frag import analyze_chiral_peak_pairs  # type: ignore[import-not-found]
    from ms2 import build_ms2_index, load_ms2_spectra, write_table  # type: ignore[import-not-found]


DEFAULT_VARIANTS = {
    "mz_tol": [2.5, 3.0, 5.0, 10.0],
    "precursor_eic_mz_tol": [3.0, 5.0],
    "fragment_eic_mz_tol": [3.0, 5.0],
    "rt_half_window_sec": [4.0, 8.0, 10.0, 12.0],
    "min_fragment_correlation": [0.7, 0.8, 0.9, 0.95],
    "consensus_scans": [1, 2, 3, 5],
    "fragment_mz_tol": [0.003, 0.005, 0.01],
}


def run_sensitivity_sweep(
    peaklist_path,
    raw_path,
    output_dir,
    *,
    mz_column="MZ",
    peak_a_column="Peak1",
    peak_b_column="Peak2",
):
    """执行单因素MS2敏感性扫描并复用一个内存raw索引。

    输出展示状态和双侧可用谱数量对参数的敏感程度；没有独立人工真值时，它
    不能估计accuracy，也不应以“得到更多supported”作为唯一优化目标。
    """

    output_dir = Path(output_dir).resolve()
    load_started = time.perf_counter()
    index = build_ms2_index(load_ms2_spectra(Path(raw_path)))
    raw_load_seconds = time.perf_counter() - load_started
    baseline = {
        "mz_tol": 10.0,
        "rt_half_window_sec": 10.0,
        "min_fragment_correlation": 0.9,
        "consensus_scans": 1,
        "fragment_mz_tol": 0.01,
    }
    summaries: list[dict] = []
    seen: set[tuple] = set()
    for parameter, values in DEFAULT_VARIANTS.items():
        for value in values:
            settings = dict(baseline)
            settings[parameter] = value
            key = tuple(sorted(settings.items()))
            if key in seen:
                continue
            seen.add(key)
            label = f"{parameter}_{value:g}"
            output = output_dir / f"{label}.csv"
            started = time.perf_counter()
            rows = analyze_chiral_peak_pairs(
                peaklist_path,
                raw_path,
                output,
                mz_column=mz_column,
                peak_a_column=peak_a_column,
                peak_b_column=peak_b_column,
                write_metadata=False,
                preloaded_index=index,
                **settings,
            )
            counts: dict[str, int] = {}
            for row in rows:
                status = row["enantiomer_pair_status"]
                counts[status] = counts.get(status, 0) + 1
            summaries.append(
                {
                    "varied_parameter": parameter,
                    "varied_value": value,
                    "raw_load_seconds_shared": raw_load_seconds,
                    "analysis_seconds": time.perf_counter() - started,
                    "rows": len(rows),
                    "both_ms2_spectra": sum(
                        bool(row["peak_a_MS2"] and row["peak_b_MS2"]) for row in rows
                    ),
                    **counts,
                    **{f"parameter_{name}": setting for name, setting in settings.items()},
                }
            )
    summary_path = output_dir / "sensitivity_summary.csv"
    write_table(summary_path, summaries)
    return summaries


def _main(argv=None):
    """解析单因素敏感性扫描参数并报告完成的variant数量。"""

    parser = argparse.ArgumentParser(
        description="Run one-factor MSAI parameter sensitivity analysis."
    )
    parser.add_argument("--peaklist", required=True)
    parser.add_argument("--raw", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--mz-column", default="MZ")
    parser.add_argument("--peak-a-column", default="Peak1")
    parser.add_argument("--peak-b-column", default="Peak2")
    args = parser.parse_args(argv)
    rows = run_sensitivity_sweep(
        args.peaklist,
        args.raw,
        args.output_dir,
        mz_column=args.mz_column,
        peak_a_column=args.peak_a_column,
        peak_b_column=args.peak_b_column,
    )
    print(f"completed_variants={len(rows)}")


if __name__ == "__main__":
    _main()
