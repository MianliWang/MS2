# MS2 static visualization and multi-reviewer plan

This workflow gives every result row an independent MS2 evidence image. It
never classifies a spectrum as single/double peak: `Peak1` and `Peak2` are
supplied extraction coordinates, while the MS2 status asks only whether the
two extracted fragment spectra support the same compound.

## Pilot artifact

The EASMSV1 pilot exports all 525 rows, including diagnostic placeholders for
records with no retained spectrum:

```text
MSAI/results/final/ig_easmsv1/spectrum_review/
├── assets/ms2/<standard>/cfg-<hash>/<run>/
│   ├── svg/                              # canonical vector images
│   └── png/                              # matching 2x raster images
├── views/
│   ├── diagnostic_status/               # supported/conflict/insufficient/not evaluable/not attempted
│   ├── review_priority/                  # P1 through P5
│   ├── spectrum_availability/            # both, Peak1-only, Peak2-only, none
│   ├── issue_codes/                      # one non-exclusive folder per issue code
│   ├── acquisition_coverage/             # acquired or outside DIA range
│   ├── source_machine/IG/                # raw IG colour result
│   ├── source_pool/                      # pool provenance
│   ├── review_queues/                    # focused parameter/failure audits
│   └── cross_stage_followup/             # MS1 context plus MS2 conflict, not a shared class
├── galleries/index.html
├── metadata/{target_manifest.csv,view_index.csv,run_manifest.json,generation_summary.json}
└── annotations/_template/review_labels.csv
```

`assets/` contains the only image entities. `views/` uses NTFS hardlinks and
hardlink failure stops the export instead of silently copying thousands of
files. The annotation directory is never deleted on rerun.

## Contents of each image

1. Compound, precursor m/z, supplied Peak1/Peak2 RT, source machine, pool, and
   well.
2. Independent MS2 diagnostic status, review priority, issue codes, and
   spectrum availability.
3. Threshold meters for square-root cosine, reported entropy similarity,
   matched fragments, and explained intensity on both sides.
4. A relative-intensity mirror spectrum. Each side is normalized to its own
   base peak; matched fragments use thick blue/orange sticks plus different
   endpoint shapes, while unmatched fragments are thin grey sticks.
5. The exact alignment parameters loaded from the result metadata sidecar:
   1% base-peak filter and 0.01 Da exclusive one-to-one matching for this run.
6. The actual acquired DIA-window strip, target precursor position, DIA
   bounds/count, and acquisition attributes aggregated over the complete raw
   XML. A separate line shows any unconfirmed reference-method profile and its
   provenance-discrepancy count; reference values never overwrite raw values.
7. Per-side RT extraction window (seconds), apex RT/intensity, area, MS2 and
   consensus scan counts, candidate/retained fragment counts, and quality
   flags.

The result table does not preserve precursor or fragment time-series arrays,
so these images do not pretend to contain raw fragment EIC curves. A future
coelution panel must reload the raw file or persist the extracted time series.

## Current pilot findings

- Diagnostic statuses: 21 supported, 4 conflicting, 28 insufficient, 351 not
  evaluable, and 121 not attempted.
- Spectrum availability: 53 two-sided, 23 Peak1-only, 9 Peak2-only, and 440
  with neither spectrum.
- All 306 `NO_ACQUIRED_DIA_WINDOW` records have precursor m/z below the actual
  acquired lower bound of 367.5. This is a coverage gap in the currently
  imported raw file, not a fragment-matching failure. The AdductMLib reference
  method reports three DIA files per sample, so matching EASMSV1 `_1`/`_3`
  files must be checked before concluding that the experiment never acquired
  these masses.
- 73 rows retain less than 5% of candidate fragments on at least one side.
  They are queued for checking whether the 0.9 coelution-correlation filter,
  precursor XIC quality, or single-scan consensus is too restrictive.
- 61 supplied RT pairs are separated by at most 20 seconds; 32 rows have only
  one retained spectrum; six have high cosine but fewer than six matches.

These are review queues, not new automatic decision rules.

## Review workflow

1. Copy the template to
   `annotations/<reviewer_id>/round_1/review_labels.csv`.
2. Review P1 conflict, P2 boundary/interference, and P3 insufficient rows
   first; then audit P4 supported rows and acquisition/attrition queues.
3. Use only `supports_same_compound`, `conflicts_same_compound`,
   `insufficient_evidence`, `not_evaluable`, or `uncertain`. Do not assign R/S
   from these spectra.
4. Record each side's spectrum quality, interference flags, diagnostic ions,
   confidence, and notes. Critical rows should receive two blinded reviews and
   adjudication.
5. Freeze the consensus labels before changing similarity, RT-window,
   correlation, fragment-intensity, or consensus-scan parameters.

The versioned annotation definition is
[`standards/ms2_manual_review_schema_v1.json`](../standards/ms2_manual_review_schema_v1.json).
