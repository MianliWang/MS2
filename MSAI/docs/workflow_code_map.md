# MS2 workflow: steps 4-8 and code ownership

This document maps the scientific workflow to one canonical implementation.
`IG`, `IA`, and `IB` are experiment/machine/batch identifiers.  A column named
`IG` contains the source label produced for the IG experiment; `GREEN`,
`YELLOW`, `RED`, and `CHECK` are values in that column, not machine IDs.

## Inputs before step 4

The core route expects one target row with at least `Compound_ID`, selected
precursor `MZ`, and supplied `Peak1`/`Peak2` RT coordinates.  A production
manifest should additionally carry experiment ID, pool, ion mode, DIA segment,
and replicate so that the row is routed to the correct raw file.  Automatic
MS1 peak picking is optional when reviewed RT coordinates already exist.

## Step 4 - raw-file routing and DIA coverage

Primary modules:

- `ms2/raw_io.py`: stream mzXML/mzML MS2 scans;
- `ms2/raw_metadata.py`: retain observed collision energy and isolation bounds;
- `ms2/indexing.py`: `build_ms2_index()` and `matching_dia_windows()`;
- `ms2/workspace.py`: strict file discovery helpers;
- `ms2/pipeline.py`: validate each row and select its real DIA window.

Output: a real DIA center/lower/upper bound or an explicit
`NO_ACQUIRED_DIA_WINDOW` reason.  Similarity thresholds are not involved here.

## Step 5 - reconstruct spectra at Peak1 and Peak2

Primary modules:

- `ms2/chromatograms.py`: raw XIC extraction and apex/area summaries;
- `ms2/extraction.py`: `extract_window_result()` performs candidate-fragment,
  intensity, chromatographic-correlation, and consensus-scan filtering;
- `ms2/pipeline.py`: `extract_fragments_for_rt_window()` runs the same operation
  independently at Peak1 and Peak2.

Output: two full-precision fragment arrays plus quality and scan-count fields.

## Step 6 - compare the two multi-fragment spectra

Primary module: `ms2/similarity.py`.

`compare_fragment_spectra()` aligns fragments one-to-one and reports cosine,
entropy similarity, matched-fragment count, and explained intensity on both
sides.  This step asks whether the spectra are mutually consistent; it does
not prescribe a theoretical fragmentation pathway or assign R/S.

## Step 7 - assign an independent MS2 diagnosis

Primary module: `ms2/diagnostics.py`.

`classify_ms2_diagnostic()` emits `supported_same_compound`,
`conflicting_spectra`, `insufficient_evidence`, `not_evaluable`, or
`not_attempted`, together with issue codes and review priority.  The source
experiment label and MS1 RT availability remain separate fields.

`ms2/pipeline.py::analyze_chiral_peak_pairs()` is the orchestration boundary
for steps 4-7 and writes the result CSV plus provenance sidecar.

## Step 8 - visualization and human review

Primary package: `ms2_review/`.

- `report.py`: compact technical HTML/CSV report;
- `exporter.py`: canonical per-target artifacts and review manifests;
- `svg.py` / `png.py`: static mirror-spectrum renderers;
- `classification.py`: independent folder views by status, issue, experiment,
  pool, and review priority;
- `acquisition.py`: displayed raw-observed versus reference-method provenance.

## Supported commands

```powershell
python -m MSAI.python.cli.ms2 analyze --help
python -m MSAI.python.cli.ms2 report --help
python -m MSAI.python.cli.ms2 export-review --help
```

Legacy script names remain available, but the module command above is the
documented interface.

## Result ownership

```text
results/
├── final/ig_easmsv1/       # current deliverable tables/report/review gallery
├── development/ms2/        # sensitivity runs and superseded MS2 tables
└── reference/ms1/          # optional MS1 calibration and review artifacts
```

Generated data remain gitignored; `results/README.md` documents the stable
locations and regeneration commands.
