# MS1 EIC visualization and multi-reviewer plan

This plan treats MS1 images as chromatographic/RT context only.  It does not
change `ms2_diagnostic_status`, and an MS1 double peak plus an MS2 conflict is a
cross-stage follow-up rather than a classification disagreement.

## Current pilot artifact

The EASMSV1 pilot writes 525 canonical SVG images plus matched 2x PNG copies
to:

```text
MSAI/results/reference/ms1/ms1_eic_review_easmsv1/
├── assets/ms1/cfg-<hash>/<raw-run>/
│   ├── svg/                               # canonical vector image per target
│   └── png/                               # matching 2x raster image per target
├── views/
│   ├── reference_ms1/supplied_rt/         # double/single/no supplied peak
│   ├── auto_ms1/baseline/                 # current automatic classes
│   ├── auto_ms1/consensus_stable/         # two parameter sets agree in class+RT
│   ├── reference_comparison/baseline/     # development audit only
│   └── review_queues/                     # instability, shoulders, disagreement
├── galleries/index.html                   # browser thumbnail entry point
├── metadata/
│   ├── target_manifest.csv
│   ├── view_index.csv
│   ├── run_manifest.json
│   └── generation_summary.json
└── annotations/_template/review_labels.csv
```

`views/` uses NTFS hardlinks to both canonical image formats.  Hardlink failure stops
the exporter instead of silently copying thousands of images.  `target_uid`
is a stable hash of dataset, run, compound, and m/z; every image also records a
configuration hash and SHA-256.

## Contents of each image

1. Full 30-minute raw EIC on a zero-based linear intensity scale.
2. A supplied-RT-centered detail panel that is not stretched by a mislocalized
   automatic peak elsewhere in the run.
3. Full `log10(1 + intensity)` EIC to expose sparse background and low signal.
4. The exact m/z window, ppm, baseline SG profile, experimental SG profile,
   supplied RTs, and high-recall candidate markers.
5. Peak RT, intensity, discrete FWHM, FWHM area, SNR, prominence, separation,
   FWHM-based Rs, valley ratio, height floor, background p99, continuous
   half-height support, parameter stability, and cross-stage MS2 status.

Color is not the only encoding: baseline calls use blue circles/solid lines,
experimental calls use orange triangles/dashes, high-recall candidates are
open purple squares, and supplied RTs are dark dashed vertical lines.

## Review workflow

1. Copy the annotation template to
   `annotations/<reviewer_id>/round_1/review_labels.csv`.
2. Review priority: parameter class flips, parameter RT shifts, baseline
   ambiguity, supplied close/shoulder doubles, baseline/reference mismatch,
   experimental regressions, then a stratified audit of stable calls.
3. Allowed human classes should include `single_peak`, `double_resolved`,
   `double_shoulder`, `no_peak`, `interference_or_coelution`, and `uncertain`.
4. Record confidence, reviewer RTs, QC flags, and notes; never move or edit the
   canonical SVG to encode the label.
5. Critical/uncertain rows should receive two blinded reviews.  A separate
   adjudicator resolves class disagreement or RT differences >0.1 min.
6. Freeze consensus labels before parameter selection.  Split future
   calibration by compound and acquisition batch, not by individual scan.

The supplied `Peak1`/`Peak2` values are experimental-design/legacy candidate
coordinates until their reviewer provenance and adjudication process are
confirmed. `IG` is a source machine/run identifier with multiple pools, not a
generic QC field. Its current colour-coded labels are preserved under
`source_machine/IG/`: GREEN = expected usable double peak, YELLOW = expected
usable single peak, RED = no peak or multiple peaks/unusable, and CHECK =
manual review. RED is deliberately not converted into an MS1 no-peak reference
label. The exporter exposes `--source-machine-label-column` so IA, IB, and
other machines can be exported independently without assuming their labels use
the IG vocabulary.
