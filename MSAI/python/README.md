# Python code map

The supported MS2 entry point is:

```powershell
python -m MSAI.python.cli.ms2 --help
```

Canonical code is grouped by workflow responsibility:

```text
python/
├── cli/ms2.py              # one router: analysis, review, and ML shadow commands
├── ms2/                    # steps 4-7: raw data through diagnosis
│   ├── raw_io.py           # mzXML/mzML MS2 streaming reader
│   ├── raw_metadata.py     # acquisition provenance
│   ├── indexing.py         # DIA-window index and coverage lookup
│   ├── chromatograms.py    # precursor/fragment XIC primitives
│   ├── extraction.py       # one RT-window spectrum reconstruction
│   ├── similarity.py       # Peak1/Peak2 spectral comparison
│   ├── diagnostics.py      # independent MS2 status vocabulary
│   └── pipeline.py         # row-level orchestration and result table
├── ms2_review/             # step 8: report, SVG/PNG, review folders
├── ms2_ml/                 # truth locks, grouped training, shadow inference/package
├── ms1_review/             # optional MS1 reference workflow
└── *_legacy_facade.py      # conceptually: old root names remain thin facades
```

The existing root modules `get_chiral_frag.py`, `chiral_similarity.py`,
`diagnostic_standards.py`, and `ms2_review_report.py` are compatibility facades.
New imports should target `MSAI.python.ms2` or `MSAI.python.ms2_review`.

See [`../docs/workflow_code_map.md`](../docs/workflow_code_map.md) for the
step-by-step call chain and output layout.

The independently labelled MS2 classifier design, frozen parameter boundary,
feature definitions, grouped evaluation contract, and current feasibility-only
status are documented in
[`../docs/ms2_shadow_classifier_plan.md`](../docs/ms2_shadow_classifier_plan.md).
