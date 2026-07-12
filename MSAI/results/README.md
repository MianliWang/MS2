# Results layout

Generated results are intentionally separated by scientific role:

```text
results/
├── final/
│   └── ig_easmsv1/
│       ├── tables/
│       │   ├── ms2_analysis.csv
│       │   └── ms2_analysis.csv.metadata.json
│       ├── report/             # compact technical report and review queue
│       └── spectrum_review/    # all per-target SVG/PNG evidence and views
├── development/
│   └── ms2/                    # sweeps, shadow profiles, superseded results
└── reference/
    └── ms1/                    # optional MS1 peak-picking/review artifacts
```

Only `final/ig_easmsv1/` is the current MS2 deliverable.  Development outputs
must not be quoted as final counts.  MS1 artifacts are method references and do
not gate the core MS2 analysis when reviewed Peak1/Peak2 coordinates exist.

Regenerate the three final layers from the repository root:

```powershell
python -m MSAI.python.cli.ms2 analyze `
  --peaklist MSAI\peaklist\merged_result_IG_EASMSV1.csv `
  --raw MSAI\data\IG_STD500nM_EASMSV1_2.mzXML `
  --output MSAI\results\final\ig_easmsv1\tables\ms2_analysis.csv

python -m MSAI.python.cli.ms2 report `
  --input MSAI\results\final\ig_easmsv1\tables\ms2_analysis.csv `
  --output-dir MSAI\results\final\ig_easmsv1\report

python -m MSAI.python.cli.ms2 export-review `
  --input MSAI\results\final\ig_easmsv1\tables\ms2_analysis.csv `
  --output-dir MSAI\results\final\ig_easmsv1\spectrum_review `
  --method-profile MSAI\config\acquisition_profiles\adductmlib_chemrxiv_20260129_v1.json `
  --source-machine-label-column IG
```
