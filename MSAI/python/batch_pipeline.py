"""Manifest-driven, resumable multi-file MS1 + chiral MS2 processing."""

from __future__ import annotations

import argparse
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

try:
    from .get_chiral_frag import analyze_chiral_peak_pairs
    from .ms1_peak_picker import PeakPickingConfig, annotate_peaklist, load_peak_config
    from .ms2 import raw_files, read_table, write_table
except ImportError:
    from get_chiral_frag import analyze_chiral_peak_pairs  # type: ignore
    from ms1_peak_picker import PeakPickingConfig, annotate_peaklist, load_peak_config  # type: ignore
    from ms2 import raw_files, read_table, write_table  # type: ignore


REQUIRED_MANIFEST_COLUMNS = {"sample_id", "mode", "peaklist", "raw", "output"}


def prepare_well_manifest(
    combined_peaklist,
    raw_dir,
    split_peaklist_dir,
    result_dir,
    manifest_path,
    *,
    group_column: str = "Orginal Pooled Well",
    mz_column: str = "MZ",
    peak_config_json=None,
):
    """Split a multi-well table and require an unambiguous raw match per well."""

    combined_peaklist = Path(combined_peaklist).resolve()
    raw_dir = Path(raw_dir).resolve()
    split_peaklist_dir = Path(split_peaklist_dir).resolve()
    result_dir = Path(result_dir).resolve()
    manifest_path = Path(manifest_path).resolve()
    manifest_base = manifest_path.parent
    rows = read_table(combined_peaklist)
    if not rows or group_column not in rows[0]:
        raise ValueError(f"Combined peaklist must contain {group_column!r}.")
    groups: dict[str, list[dict]] = {}
    for row in rows:
        group = row.get(group_column, "").strip().upper()
        if group:
            groups.setdefault(group, []).append(row)

    available_raw = raw_files(raw_dir)
    manifest_rows: list[dict] = []
    for group in sorted(groups):
        matches = [path for path in available_raw if group.lower() in path.stem.lower()]
        if len(matches) != 1:
            names = ", ".join(path.name for path in matches) or "none"
            raise ValueError(f"Expected exactly one raw file for {group}; found: {names}")
        peaklist_path = split_peaklist_dir / f"{combined_peaklist.stem}_{group}.csv"
        output_path = result_dir / f"{combined_peaklist.stem}_{group}_full.csv"
        write_table(peaklist_path, groups[group])
        job = {
            "sample_id": group,
            "mode": "full",
            "peaklist": os.path.relpath(peaklist_path, manifest_base),
            "raw": os.path.relpath(matches[0], manifest_base),
            "output": os.path.relpath(output_path, manifest_base),
            "mz_column": mz_column,
        }
        if peak_config_json:
            job["peak_config_json"] = os.path.relpath(Path(peak_config_json).resolve(), manifest_base)
        manifest_rows.append(job)
    write_table(manifest_path, manifest_rows)
    return manifest_rows


def run_batch_manifest(manifest_path, *, jobs: int = 2, resume: bool = True):
    """Run independent raw files concurrently and isolate per-file failures."""

    manifest_path = Path(manifest_path).resolve()
    rows = read_table(manifest_path)
    if not rows:
        raise ValueError("Batch manifest is empty.")
    missing = REQUIRED_MANIFEST_COLUMNS - set(rows[0])
    if missing:
        raise ValueError(f"Batch manifest is missing columns: {', '.join(sorted(missing))}")
    for row in rows:
        for field in ("peaklist", "raw", "output", "peak_config_json"):
            if not row.get(field):
                continue
            path = Path(row[field])
            row[field] = str(path if path.is_absolute() else (manifest_path.parent / path).resolve())
    jobs = max(1, min(int(jobs), os.cpu_count() or 1))

    indexed_results: dict[int, dict] = {}
    if jobs == 1:
        for index, row in enumerate(rows):
            try:
                indexed_results[index] = _run_job(row, resume)
            except Exception as exc:
                indexed_results[index] = {
                    "sample_id": row.get("sample_id", ""),
                    "status": "failed",
                    "error": f"{type(exc).__name__}: {exc}",
                }
    else:
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            future_map = {
                executor.submit(_run_job, row, resume): index
                for index, row in enumerate(rows)
            }
            for future in as_completed(future_map):
                index = future_map[future]
                try:
                    indexed_results[index] = future.result()
                except Exception as exc:  # isolate one sample without hiding its error
                    indexed_results[index] = {
                        "sample_id": rows[index].get("sample_id", ""),
                        "status": "failed",
                        "error": f"{type(exc).__name__}: {exc}",
                    }

    results = [indexed_results[index] for index in range(len(rows))]
    output_parents = [str(Path(row["output"]).parent) for row in rows]
    summary_dir = Path(os.path.commonpath(output_parents))
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary_path = summary_dir / "batch_summary.json"
    summary_path.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    return results


def _run_job(job: dict, resume: bool) -> dict:
    started = time.perf_counter()
    sample_id = job["sample_id"]
    mode = job["mode"].strip().lower()
    peaklist = Path(job["peaklist"])
    raw = Path(job["raw"])
    output = Path(job["output"])
    mz_column = job.get("mz_column", "MZ") or "MZ"
    peak_config = _peak_config(job)
    metadata = output.with_suffix(output.suffix + ".metadata.json")
    if resume and output.exists() and (mode != "full" or metadata.exists()):
        if mode == "full" and metadata.exists():
            _augment_full_metadata(metadata, sample_id, peak_config)
        cached = json.loads(metadata.read_text(encoding="utf-8")) if metadata.exists() else {}
        return {
            "sample_id": sample_id,
            "status": "skipped_existing",
            "output": str(output),
            "row_count": cached.get("rows_processed"),
            "status_counts": cached.get("status_counts", {}),
            "elapsed_seconds": 0.0,
        }

    if mode == "ms1":
        rows = annotate_peaklist(peaklist, raw, output, mz_column=mz_column, config=peak_config)
    elif mode == "ms2":
        rows = analyze_chiral_peak_pairs(
            peaklist,
            raw,
            output,
            mz_column=mz_column,
            peak_a_column=job.get("peak_a_column", "Peak1") or "Peak1",
            peak_b_column=job.get("peak_b_column", "Peak2") or "Peak2",
        )
    elif mode == "full":
        intermediate = output.with_name(output.stem + "_ms1.csv")
        annotate_peaklist(
            peaklist,
            raw,
            intermediate,
            mz_column=mz_column,
            config=peak_config,
        )
        rows = analyze_chiral_peak_pairs(
            intermediate,
            raw,
            output,
            mz_column=mz_column,
            peak_a_column="auto_peak1_rt_min",
            peak_b_column="auto_peak2_rt_min",
        )
        _augment_full_metadata(metadata, sample_id, peak_config)
    else:
        raise ValueError(f"Unsupported batch mode {mode!r} for {sample_id}.")

    status_counts: dict[str, int] = {}
    for row in rows:
        status = row.get("enantiomer_pair_status") or row.get("auto_chromatographic_status", "")
        status_counts[status] = status_counts.get(status, 0) + 1
    return {
        "sample_id": sample_id,
        "status": "completed",
        "output": str(output),
        "row_count": len(rows),
        "status_counts": status_counts,
        "elapsed_seconds": time.perf_counter() - started,
    }


def _peak_config(job: dict) -> PeakPickingConfig:
    path = job.get("peak_config_json", "").strip()
    if not path:
        return PeakPickingConfig()
    return load_peak_config(path)


def _augment_full_metadata(metadata_path: Path, sample_id: str, config: PeakPickingConfig) -> None:
    payload = json.loads(metadata_path.read_text(encoding="utf-8"))
    payload["batch_sample_id"] = sample_id
    payload["upstream_ms1_peak_picking"] = asdict(config)
    metadata_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def _main(argv=None):
    parser = argparse.ArgumentParser(description="Prepare or run explicit multi-file MSAI manifests.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare", help="Split a combined well table and build a strict manifest.")
    prepare.add_argument("--combined-peaklist", required=True)
    prepare.add_argument("--raw-dir", required=True)
    prepare.add_argument("--split-peaklist-dir", required=True)
    prepare.add_argument("--result-dir", required=True)
    prepare.add_argument("--manifest", required=True)
    prepare.add_argument("--group-column", default="Orginal Pooled Well")
    prepare.add_argument("--mz-column", default="MZ")
    prepare.add_argument("--peak-config-json")
    run = subparsers.add_parser("run", help="Run a prepared manifest.")
    run.add_argument("--manifest", required=True)
    run.add_argument("--jobs", type=int, default=2)
    run.add_argument("--no-resume", action="store_true")
    args = parser.parse_args(argv)

    if args.command == "prepare":
        rows = prepare_well_manifest(
            args.combined_peaklist,
            args.raw_dir,
            args.split_peaklist_dir,
            args.result_dir,
            args.manifest,
            group_column=args.group_column,
            mz_column=args.mz_column,
            peak_config_json=args.peak_config_json,
        )
        print(f"prepared_jobs={len(rows)}")
    else:
        results = run_batch_manifest(args.manifest, jobs=args.jobs, resume=not args.no_resume)
        counts: dict[str, int] = {}
        for result in results:
            counts[result["status"]] = counts.get(result["status"], 0) + 1
        print("; ".join(f"{key}={counts[key]}" for key in sorted(counts)))


if __name__ == "__main__":
    _main()
