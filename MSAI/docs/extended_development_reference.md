# Extended development and validation reference

> This is the detailed parameter, validation, performance, and legacy-command
> reference preserved from the former long README. Start with
> [`../README.md`](../README.md) for the supported workflow and current paths.

MSAI is a dependency-light Python workflow whose primary goal is DIA-MS2
extraction and spectrum comparison for two supplied chiral LC retention-time
peaks.  It also contains experimental targeted-MS1 peak-picking utilities for
method development, but MS1 automation is a reference workflow rather than the
main deliverable.

The supported project path is therefore: reviewed or otherwise supplied
`Peak1`/`Peak2` values -> independent DIA-MS2 extraction -> cosine/entropy and
matched-fragment evidence -> candidate chiral-doublet screening.  Automatic
MS1 peak classification and E-ASMS enrichment statistics remain optional
research helpers and do not gate the core MS2 workflow.

## Scientific scope

The workflow deliberately separates four conclusions:

1. `chromatographic_status`: the supplied/manual peak context, or an optional
   experimental MS1 EIC call;
2. `compound_identity_status`: whether the two extracted MS2 spectra support
   the same compound identity;
3. `chiral_doublet_status`: whether a chromatographic doublet is supported as a
   candidate chiral pair;
4. `enantioselective_enrichment_status`: requires input/eluate peak areas,
   same-batch controls, and replicates, and is not inferred from the current
   standard-pool data.

MS2 similarity alone cannot assign R/S configuration or prove
enantioselective enrichment.  Authentic enantiopure standards and orthogonal
validation remain required.

## What changed in the upgraded workflow

- reads conventional ProteoWizard mzXML and mzML with the Python standard
  library; `pyopenms` is optional;
- retains the real isolation lower/upper bounds from every raw file instead of
  assuming one fixed DIA width;
- resolves overlapping DIA windows deterministically by nearest center and
  reports the number of matching windows;
- builds a reusable `DIA window -> sorted RT scans` index and uses binary search
  for RT/m/z extraction;
- compares in-memory full-precision fragment arrays, avoiding CSV rounding
  before scoring;
- reports square-root cosine and weighted spectral-entropy similarity;
- supports multi-scan candidate consensus and configurable fragment intensity,
  relative intensity, and coelution-correlation thresholds;
- removes fully blank spreadsheet rows while retaining all partial records;
- writes a sidecar metadata JSON with input SHA-256, parameters, runtime, status
  counts, and acquired DIA windows;
- adds automatic MS1 EIC smoothing, peak detection, raw-apex refinement, FWHM
  area, S/N, valley ratio, and chromatographic resolution;
- adds independent threshold calibration and one-factor sensitivity tools;
- adds strict manifest-based multi-file processing, process-level parallelism,
  checkpoint/resume, and per-file failure isolation;
- provides explicitly named target/control enrichment, input/eluate recovery,
  enantiomer-ratio-shift, and Benjamini-Hochberg helpers for future E-ASMS
  datasets.

## Project layout

```text
MSAI/
├── config/       # method-specific calibrated profiles
├── standards/    # versioned, independent MS1-reference and MS2 rules
├── docs/         # scientific definitions and current visual audit notes
├── data/         # mzXML/mzML raw files (gitignored)
├── peaklist/     # input tables and generated well splits (gitignored)
├── results/      # final, development, and MS1-reference outputs (gitignored)
├── python/
│   ├── cli/      # supported analyze/report/export-review command router
│   ├── ms2/      # steps 4-7: raw I/O, DIA routing, extraction, similarity, diagnosis
│   ├── ms2_review/ # step 8: reports, SVG/PNG, and human-review views
│   └── ms2_core.py # small backward-compatible import facade
└── batch_*.csv   # explicit reproducible batch manifests
```

The former 940-line `ms2_core.py` is now split under `python/ms2/`:

| Module | Responsibility |
| --- | --- |
| `models.py` | `Spectrum`, `DiaData`, `Ms2Index`, and `Eic` data structures |
| `xml_codec.py` / `raw_io.py` | mzXML/mzML binary decoding and streaming MS2 reads |
| `indexing.py` | one-time DIA-window grouping and nearest-window lookup |
| `chromatograms.py` | XIC extraction, apex/area summaries, and result schemas |
| `extraction.py` | precursor/fragment coelution and consensus fragment selection |
| `similarity.py` | one-to-one fragment alignment, cosine, entropy, and explained intensity |
| `diagnostics.py` | independent MS2 status and review-priority rules |
| `pipeline.py` | row-level orchestration and result/provenance writing |
| `tables.py` / `workspace.py` | tabular I/O, project paths, and strict raw matching |
| `utils.py` | m/z/RT bounds, finite-number checks, and correlation helpers |

Project code imports the package directly.  `ms2_core.py` re-exports the same
objects and signatures so existing callers do not break.

## 1. Analyze manually reviewed Peak1/Peak2 values

From the repository root:

```powershell
python -m MSAI.python.cli.ms2 analyze `
  --peaklist MSAI\peaklist\merged_result_IG_EASMSV1.csv `
  --raw MSAI\data\IG_STD500nM_EASMSV1_2.mzXML `
  --output MSAI\results\final\ig_easmsv1\tables\ms2_analysis.csv
```

The compatibility defaults remain a shared 10 ppm precursor/fragment EIC
tolerance, 8-second RT half-window, Pearson coelution correlation >0.9, one candidate
scan, 0.01 Da spectrum-alignment tolerance, cosine >=0.7, at least six matched
ions, and at least 50% explained intensity on both sides.  These defaults are
screening guardrails, not validated universal identity thresholds.

For calibrated or shadow runs, precursor and fragment coelution EIC tolerances
can now be separated with `--precursor-eic-mz-tol` and
`--fragment-eic-mz-tol` (plus their unit options). `--dia-window-fallback`
is the explicit name for the old `--dia-window` option: it is used only when
the raw file lacks real isolation bounds. The legacy option remains accepted.

`--fragment-correlation-mode active_support` enables an experimental
common-zero-resistant Pearson calculation. It uses only precursor active
support and exposes minimum support scans, apex-offset, and consecutive-scan
guards. The frozen shadow settings are documented under
`config/ms2_processing_profiles/`; `full_window` remains the production-
compatibility default until blinded review is available.

Entropy similarity is always reported.  A calibrated guardrail can be enabled
with `--min-entropy-similarity`; it is disabled by default because a published
library-search threshold is not automatically transferable to this DIA method.

### Independent MS2 diagnosis and manual visualization

The pipeline now writes `ms1_reference_status` and `ms2_diagnostic_status` as
separate columns.  MS1 says only whether two RT coordinates were supplied;
MS2 independently reports `supported_same_compound`, `conflicting_spectra`,
`insufficient_evidence`, `not_evaluable`, or `not_attempted`.  The legacy
integrated columns remain for compatibility.

Generate the MS2 technical report, mirror-spectrum gallery, and editable
manual-review queue from an existing result:

```powershell
python -m MSAI.python.cli.ms2 report `
  --input MSAI\results\final\ig_easmsv1\tables\ms2_analysis.csv `
  --output-dir MSAI\results\final\ig_easmsv1\report
```

The static report works without JavaScript.  A local Data Analytics runtime
helper can optionally be passed through `--runtime-helper` to upgrade the three
summary charts to live Recharts while preserving the same-data SVG fallbacks.
The review queue intentionally leaves `manual_ms2_label`, reviewer, date, issue
codes, and notes blank for blinded labeling.

Generate one standalone SVG and matching 2x PNG for every result row, plus
hard-linked filesystem views by status, priority, spectrum availability,
issue, DIA coverage, source machine, and pool:

```powershell
py -3.12 -m MSAI.python.cli.ms2 export-review `
  --input MSAI\results\final\ig_easmsv1\tables\ms2_analysis.csv `
  --output-dir MSAI\results\final\ig_easmsv1\spectrum_review `
  --method-profile MSAI\config\acquisition_profiles\adductmlib_chemrxiv_20260129_v1.json `
  --source-machine-label-column IG
```

The exporter automatically reads `<input>.metadata.json`, so image matching,
threshold, raw-file fingerprint, and acquired DIA-window labels agree with the
analysis that produced the CSV. It creates honest reason/coverage images for
empty spectra instead of excluding them. Human labels live in preserved
per-reviewer copies of `annotations/_template/review_labels.csv`.
The optional method profile is displayed as `reference method`, while all
values aggregated from the complete mzXML are displayed as `raw observed`.
The two layers are deliberately not merged. See
[`docs/acquisition_parameter_reconciliation.md`](docs/acquisition_parameter_reconciliation.md)
for the 5/4.5-vs-15 m/z and stepped-15/30/60-vs-35 provenance audit.

See `standards/ms2_diagnostic_standard_v1.json` and
`docs/ms2_diagnostic_standard.md` for the decision sequence and calibration
boundary. See `docs/ms2_static_visual_review_plan.md` and
`standards/ms2_manual_review_schema_v1.json` for static image review. MS1 remains separately documented in
`standards/ms1_reference_standard_v1.json` and
`docs/ms1_reference_standard.md`.

## 2. Automatically detect MS1 chiral peaks

Manuscript-starting parameters:

```powershell
python MSAI\python\ms1_peak_picker.py `
  --peaklist MSAI\peaklist\merged_result_IG_EASMSV1.csv `
  --raw MSAI\data\IG_STD500nM_EASMSV1_2.mzXML `
  --output MSAI\results\reference\ms1\merged_result_IG_EASMSV1_auto_peaks.csv
```

Validated IG profile for the current standardized method:

```powershell
python MSAI\python\ms1_peak_picker.py `
  --peaklist MSAI\peaklist\merged_result_IG_EASMSV1.csv `
  --raw MSAI\data\IG_STD500nM_EASMSV1_2.mzXML `
  --output MSAI\results\reference\ms1\merged_result_IG_EASMSV1_auto_peaks_calibrated.csv `
  --config-json MSAI\config\ig_peak_picker_calibrated.json
```

The profile is method-specific.  Recalibrate after changing column, gradient,
concentration, instrument, or pool design.

### MS1 implementation and peak metrics

For a target mass `m`, each MS1 scan contributes the maximum centroid
intensity inside `m +/- m * eic_ppm / 1e6`.  All target EICs are collected in
one streaming pass over the raw file.  Each trace then follows this sequence:

1. smooth the EIC with a zero-derivative Savitzky-Golay filter;
2. find local maxima above `min_height` and greedily enforce the minimum scan
   distance, strongest first;
3. refine every candidate to the raw-EIC maximum within +/-3 scans;
4. retain the two strongest distinct raw apices and order them by RT;
5. compute quality metrics on the unsmoothed EIC and apply valley, resolution,
   and second-peak-ratio guardrails.

The calibrated profile uses a seven-scan, second-order Savitzky-Golay filter.
The manuscript supplied the 3 ppm EIC, odd smoothing window, 20-scan distance,
top-two selection, and raw-apex refinement design.  Polynomial order 2,
refinement radius 3, minimum five scans, and top two are engineering choices.
The calibration grid evaluated SG windows 5/7/9; heights 100k/200k/500k;
distances 15/20/25/30; valley ratios 0.5/0.7/0.85/1.0; resolution thresholds
0/0.5/1.0; and second-peak ratios 0.05/0.1/0.2/0.3, for 1,728 total
combinations.

Metrics use the raw trace and a global median baseline `b`:

- robust noise is `1.4826 * median(|x-b|)` and
  `SNR = max(0, (apex-b)/noise)`;
- FWHM uses the first sampled points on each side below
  `b + (apex-b)/2`; no between-scan interpolation is applied;
- `area_fwhm` is the trapezoidal integral of `max(0, x-b)` between those
  FWHM bounds, in intensity-seconds, rather than whole-peak area;
- two-peak resolution from measured half-height widths is
  `Rs = 1.17741*(RT2-RT1)/(FWHM1+FWHM2)`; the factor-2 IUPAC form is for
  baseline widths rather than FWHM;
- valley ratio is the minimum raw signal between apices divided by the weaker
  raw-apex intensity.

For sparse centroided traces, the median absolute deviation is often zero.  In
that case SNR is deliberately left blank rather than reported as infinity;
quantitative SNR needs a prespecified off-peak or blank-derived noise model.
The reported `Rs` uses measured discrete FWHM, whereas the supplied manuscript
used fixed 0.5-minute widths for one labeling rule, so those values must not be
treated as identical definitions.

### Export MS1 EIC images for independent human review

MS1 image review remains a separate RT/chromatography layer and never changes
the MS2 same-compound status. Generate one canonical SVG and matching 2x PNG
per target, plus NTFS hard-linked folder views:

```powershell
python MSAI\python\export_ms1_eic_review.py `
  --peaklist MSAI\peaklist\merged_result_IG_EASMSV1.csv `
  --raw MSAI\data\IG_STD500nM_EASMSV1_2.mzXML `
  --output-dir MSAI\results\reference\ms1\ms1_eic_review_easmsv1 `
  --baseline-config MSAI\config\ig_peak_picker_calibrated.json `
  --experimental-config MSAI\config\ig_peak_picker_rt_aware_experimental.json `
  --candidate-config MSAI\config\ig_peak_picker_high_recall_candidates_v1.json `
  --ms2-results MSAI\results\final\ig_easmsv1\tables\ms2_analysis.csv `
  --source-machine-label-column IG
```

Open `galleries/index.html` to browse by supplied RT reference, automatic
status, parameter instability, close/shoulder peaks, or reference disagreement.
Each gallery previews PNG but links to both PNG and SVG. Use `--formats svg` or
`--formats png` to export only one format, and `--png-scale` to change the
default 2x raster size.
`IG` is the source machine/run identifier; its colour labels are retained as
source context, together with pool metadata, rather than converted into a
universal chromatographic truth. Use `--source-machine-label-column IA` or
`IB` for a separately supplied machine-specific table after confirming that
source's label vocabulary.
The generated manifest is under `metadata/`; human labels belong in separate
per-reviewer copies of `annotations/_template/review_labels.csv`, which the
exporter never deletes.  See `docs/ms1_eic_visual_review_plan.md` and
`docs/ms1_parameter_audit_EASMSV1.md`.

Calibration uses a deterministic compound-level 80/20 split and maximizes
balanced accuracy, then precision and recall.  It covers 404 reviewed doubles
and 99 reviewed singles.  The 22 no-peak records are excluded, so the reported
held-out metrics validate double-versus-single agreement, not no-signal
detection, chemical identity, or universal transferability.  The lightweight
SG/local-max implementation also need not reproduce SciPy's edge and plateau
behavior exactly.

### MS1 versus MS2

| Aspect | MS1 stage | MS2 stage |
| --- | --- | --- |
| Input signal | precursor-ion intensity over chromatographic time | fragment ions inside an acquired DIA window |
| Main operation | smooth and locate zero/one/two chromatographic peaks | extract coeluting fragments at each supplied/picked RT peak |
| Evidence produced | peak RT, intensity, FWHM area, SNR, valley ratio, Rs | fragment spectrum, cosine/entropy similarity, matched ions, explained intensity |
| Primary question | "Are two chromatographic peaks present and separated?" | "Do both peaks have sufficiently similar fragment evidence?" |
| Main failure modes | noise, shoulders, baseline drift, under-resolved peaks | window interference, sparse fragments, wrong DIA window, weak coelution |

MS1 therefore proposes the two RT regions; MS2 evaluates molecular-fragment
consistency between them.  Neither stage alone assigns R/S configuration or
demonstrates selective enrichment.

### Experimental adaptive MS1 mode

The historical absolute height and 20-scan distance rules are intentionally
preserved for reproducibility.  An opt-in adaptive detector now adds:

- baseline-corrected peak prominence and half-prominence scan support;
- an RT distance in seconds rather than assuming a fixed scan rate;
- a valley exception that can retain two nearby candidates;
- valley-bounded raw-apex refinement so adjacent candidates do not collapse to
  the same stronger apex;
- explicit `auto_manual_review_required` and `auto_review_reasons` columns for
  close pairs, shallow valleys, weak second peaks, and runs with no calibrated
  absolute height floor.

This follows the same broad motivation as centWave: chromatographic feature
detection should use local peak shape and support close or partially
overlapping features, rather than one global amplitude/distance rule.  Peak
prominence follows the topographic definition used by SciPy.  See the
[centWave paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC2639432/) and
[SciPy prominence documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.peak_prominences.html).

The EASMSV1 error audit found that 25/404 reviewed doubles were separated by
fewer than 20 MS1 scans, but only 2/25 legacy automatic peak pairs matched both
reviewed RTs within 0.1 min.  In the original holdout's 14 false negatives,
six were closer than 20 scans, ten had a reviewed raw valley ratio above 0.7,
and seven automatic pairs were at the wrong RTs.  Conversely, none of the 404
reviewed doubles had a weaker local apex below 200k.  Also, 401/404 traces had
global median and MAD equal to zero.  These observations justify changing the
close-peak and scoring logic, but they do not validate low-intensity calls.

Adaptive calibration is RT-aware and three-class: a double is correct only if
both predicted RTs match reviewed Peak1/Peak2 within 0.1 min; single and
explicitly reviewed no-peak rows are also scored; wrong-location pairs are
`mislocalized`.  Run it with:

```powershell
python MSAI\python\tune_peak_picker.py `
  --peaklist MSAI\peaklist\merged_result_IG_EASMSV1.csv `
  --raw MSAI\data\IG_STD500nM_EASMSV1_2.mzXML `
  --output MSAI\results\reference\ms1\peak_picker_adaptive_calibration.json `
  --adaptive
```

Apply the experimental profile while exporting only the uncertain rows:

```powershell
python MSAI\python\ms1_peak_picker.py `
  --peaklist MSAI\peaklist\merged_result_IG_EASMSV1.csv `
  --raw MSAI\data\IG_STD500nM_EASMSV1_2.mzXML `
  --output MSAI\results\reference\ms1\merged_result_IG_EASMSV1_auto_rt_aware.csv `
  --review-output MSAI\results\reference\ms1\merged_result_IG_EASMSV1_manual_review.csv `
  --config-json MSAI\config\ig_peak_picker_rt_aware_experimental.json
```

The current exploratory profile is
`config/ig_peak_picker_rt_aware_experimental.json`.  On the same strict metric,
accuracy increased from 0.761 to 0.826 and macro recall from 0.790 to 0.872,
mainly by improving single-peak and ordinary double-peak handling.  However,
the search still selected a 200k absolute height floor and close-double recall
was zero on the reused exploratory holdout.  It is therefore a shadow-mode
review profile, not a production replacement.  The current 500 nM labels have
no reviewed double with a weaker peak below 200k, so low-clean ChiralRTdb calls
cannot be validated without adding those manual labels.

## 3. Process A01-A11 as an explicit batch

The combined table must first be split by well.  This prevents the old unsafe
behavior of silently matching a multi-well table to the first raw file.

```powershell
python MSAI\python\batch_pipeline.py prepare `
  --combined-peaklist MSAI\peaklist\ChiralAlist_ID_IG_A01-A11.csv `
  --raw-dir MSAI\data `
  --split-peaklist-dir MSAI\peaklist\by_well `
  --result-dir MSAI\results\reference\ms1\batch_A01_A11_calibrated `
  --manifest MSAI\batch_A01_A11_calibrated.csv `
  --peak-config-json MSAI\config\ig_peak_picker_calibrated.json

python MSAI\python\batch_pipeline.py run `
  --manifest MSAI\batch_A01_A11_calibrated.csv `
  --jobs 2
```

Existing outputs with valid metadata are skipped.  Use `--no-resume` to force a
rerun.  The summary is written to the batch result directory as
`batch_summary.json`.

## 4. Parameter sensitivity and calibration

One-factor sensitivity analysis parses the raw file once and reuses its index:

```powershell
python MSAI\python\parameter_sweep.py `
  --peaklist MSAI\peaklist\merged_result_IG_EASMSV1.csv `
  --raw MSAI\data\IG_STD500nM_EASMSV1_2.mzXML `
  --output-dir MSAI\results\development\ms2\sensitivity_EASMSV1
```

The sweep now varies the legacy shared EIC tolerance as well as precursor-only
and fragment-EIC-only tolerances, and includes the manuscript-referenced
10-second RT window. To turn any candidate run into a reproducible human-review
queue, compare it with the frozen baseline:

```powershell
python MSAI\python\compare_ms2_runs.py `
  --baseline MSAI\results\final\ig_easmsv1\tables\ms2_analysis.csv `
  --candidate MSAI\results\development\ms2\sensitivity_EASMSV1_provenance\precursor10_fragment10_rt10.csv `
  --output MSAI\results\development\ms2\sensitivity_EASMSV1_provenance\rt10_status_change_review_queue.csv `
  --status-changes-only
```

The comparison preserves candidate spectra plus baseline statuses/scores, so
the standard MS2 exporter can make a focused SVG/PNG gallery without replacing
the main result.

Tune MS2 decision thresholds only on a held-out table containing real positive
and negative standard labels:

```powershell
python MSAI\python\calibrate_thresholds.py `
  --input labeled_standard_results.csv `
  --label-column truth `
  --objective balanced_accuracy `
  --output calibrated_ms2_thresholds.json
```

Reproduce the automatic peak-picker calibration:

```powershell
python MSAI\python\tune_peak_picker.py `
  --peaklist MSAI\peaklist\merged_result_IG_EASMSV1.csv `
  --raw MSAI\data\IG_STD500nM_EASMSV1_2.mzXML `
  --output MSAI\results\reference\ms1\peak_picker_calibration.json
```

Do not select parameters on the same rows later used to claim accuracy.
`tune_peak_picker.py` uses a deterministic compound-level 80/20 holdout.

## 5. E-ASMS effect sizes and multiple testing

`easms_statistics.py` keeps distinct biological questions under distinct
function names:

- target/control enrichment:
  `EF = (target_eluate + p) / (mean(background_control_eluates) + p)`;
- input/eluate recovery:
  `RF = (eluate + p) / (input + p)` after compatible normalization for
  injection amount, dilution, total volume, and preferably an internal standard;
- peak-1 fraction: `f1 = (A1+p)/(A1+A2+2p)`, with fraction shift
  `delta_f1 = f1_eluate - f1_input`;
- enantiomer log-ratio shift:
  `S = log2((E1+p)/(E2+p)) - log2((I1+p)/(I2+p))`.

Positive `S` means peak 1 became relatively enriched and negative `S` means
peak 2 did.  A magnitude of one is approximately a two-fold relative-ratio
change when the pseudocount is negligible.  Peak 1/2 must not be called R/S
without authentic standards, and `p` must be predeclared from blank/LOD/LOQ
behavior because it strongly affects low-signal ratios.

Benjamini-Hochberg takes already-calculated p-values, sorts the `m` tests, and
computes monotone adjusted values from `m*p(i)/i`.  It controls the expected
false-discovery rate across the declared test family; a q-value is not the
probability that one specific hit is false.  The family must be defined before
filtering (for example, every compound tested for one protein), and technical
injections must not be counted as independent biological replicates.  The
helper adjusts supplied p-values only; it does not choose a statistical model
or generate p-values from three intensities.

## Validation results

### Reviewed EASMSV1 table

The source CSV contains 525 non-empty records plus 10,255 fully blank trailing
rows.  Blank rows are now removed during input:

- 404 reviewed double peaks;
- 99 reviewed single peaks;
- 22 rows without a detected peak;
- 21 MS2-supported candidate chiral doublets;
- 28 doublets with insufficient matched-fragment evidence;
- 4 populated but MS2-dissimilar doublets;
- 351 doublets not evaluable from this raw MS2 acquisition.

The previous 20-candidate result changed to 21 because six targets occupy
overlapping acquired DIA windows.  The upgraded code selects the nearest real
window instead of whichever center appeared first in the file.

Automatic MS1 peak-picker held-out performance against reviewed double/single
labels:

- balanced accuracy: 0.919;
- precision: 1.000;
- recall: 0.837;
- specificity: 1.000;
- 68/72 predicted-positive test rows had both RTs within 0.1 min of review.

This validates agreement with reviewed chromatographic labels, not chemical
identity.

### A01-A11 full batch

- 11/11 raw files completed;
- 3,631 compounds processed;
- approximately 2.4 GB raw input;
- two-worker wall time: about 110 s;
- resume check: about 0.46 s;
- calibrated-profile screen: 1,044 candidate doublets, 1,289 insufficient-MS2,
  119 MS2-dissimilar, 962 not evaluable, and 217 not evaluated.

These are screening counts, not precision estimates, because A01-A11 do not
provide independent MS2 identity truth labels.

### Performance

On the 65.6 MB EASMSV1 mzXML, cProfile runtime decreased from 15.98 s to
6.10 s (2.6x including profiler and metadata work).  Ordinary end-to-end wall
time decreased from about 14.3 s to about 3.3 s (about 4.3x).  The largest gain
comes from eliminating per-row DIA-window reconstruction and using RT/m/z
binary search.

## Literature basis

- The supplied ChiralRTlib manuscript motivated 3 ppm targeted MS1 EICs,
  Savitzky-Golay smoothing, minimum 20-scan peak distance, raw-apex refinement,
  explicit negative data, held-out threshold selection, and balanced accuracy.
- [E-ASMS](https://doi.org/10.1038/s41467-025-67403-2) motivated strict
  MS/RT database matching, triplicates, batch controls, QC/blank monitoring,
  MS2 interference confirmation, and separating identity from enantioselective
  enrichment.
- The supplied AdductMLib ChemRxiv preprint documents a related Exploris 240
  design with three DIA files per sample, seven spectra per loop, RT alignment
  within 10 s, and chromatographic correlation above 0.9. It is retained as a
  reference method because its 6-min C18 acquisition is not the current 30-min
  EASMSV1 chiral run.
- [MS-DIAL](https://pmc.ncbi.nlm.nih.gov/articles/PMC4449330/) supports using
  precursor-fragment chromatographic coelution to deconvolve DIA spectra.
- [Spectral entropy](https://www.nature.com/articles/s41592-021-01331-z)
  motivates reporting an entropy-based score alongside cosine because of its
  improved noise robustness in large small-molecule benchmarks.
- [LC-MS QC guidance](https://pmc.ncbi.nlm.nih.gov/articles/PMC5960010/)
  supports periodic pooled QC and blanks; the supplied E-ASMS paper used both
  after every 12 samples and maintained QC response deviation below 10%.

## Tests

```powershell
python -m unittest discover -s tests -v
```

Tests cover binary mzXML decoding, isolation-window overlap, indexed
extraction, exclusive MS2 matching, entropy similarity, automatic MS1 peak
picking, blank-row removal, threshold calibration, evidence-layer semantics,
and E-ASMS statistics helpers.
