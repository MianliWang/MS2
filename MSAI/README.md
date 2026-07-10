# MSAI chiral LC-MS processing

MSAI is a dependency-light Python workflow whose primary goal is DIA-MS2
extraction and spectrum comparison for two supplied chiral LC retention-time
peaks.  It also contains experimental targeted-MS1 peak-picking utilities for
method development, but MS1 automation is a reference workflow rather than the
main deliverable.  The original R functions remain available for legacy use.

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
├── data/         # mzXML/mzML raw files (gitignored)
├── peaklist/     # input tables and generated well splits (gitignored)
├── results/      # generated CSV/JSON results (gitignored)
├── python/
│   ├── ms2/      # focused MS2 models, readers, index, XIC, extraction, and I/O
│   └── ms2_core.py # small backward-compatible import facade
├── R/            # legacy package implementation
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
| `tables.py` / `workspace.py` | tabular I/O, project paths, and strict raw matching |
| `utils.py` | m/z/RT bounds, finite-number checks, and correlation helpers |

Project code imports the package directly.  `ms2_core.py` re-exports the same
objects and signatures so existing callers do not break.

## 1. Analyze manually reviewed Peak1/Peak2 values

From the repository root:

```powershell
python MSAI\python\get_chiral_frag.py `
  --peaklist MSAI\peaklist\merged_result_IG_EASMSV1.csv `
  --raw MSAI\data\IG_STD500nM_EASMSV1_2.mzXML `
  --output MSAI\results\merged_result_IG_EASMSV1_chiral_MS2_similarity.csv
```

The compatibility defaults remain 10 ppm precursor/fragment EIC tolerance,
8-second RT half-window, Pearson coelution correlation >0.9, one candidate
scan, 0.01 Da spectrum-alignment tolerance, cosine >=0.7, at least six matched
ions, and at least 50% explained intensity on both sides.  These defaults are
screening guardrails, not validated universal identity thresholds.

Entropy similarity is always reported.  A calibrated guardrail can be enabled
with `--min-entropy-similarity`; it is disabled by default because a published
library-search threshold is not automatically transferable to this DIA method.

## 2. Automatically detect MS1 chiral peaks

Manuscript-starting parameters:

```powershell
python MSAI\python\ms1_peak_picker.py `
  --peaklist MSAI\peaklist\merged_result_IG_EASMSV1.csv `
  --raw MSAI\data\IG_STD500nM_EASMSV1_2.mzXML `
  --output MSAI\results\merged_result_IG_EASMSV1_auto_peaks.csv
```

Validated IG profile for the current standardized method:

```powershell
python MSAI\python\ms1_peak_picker.py `
  --peaklist MSAI\peaklist\merged_result_IG_EASMSV1.csv `
  --raw MSAI\data\IG_STD500nM_EASMSV1_2.mzXML `
  --output MSAI\results\merged_result_IG_EASMSV1_auto_peaks_calibrated.csv `
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
- two-peak resolution is `Rs = 2*(RT2-RT1)/(FWHM1+FWHM2)`;
- valley ratio is the minimum raw signal between apices divided by the weaker
  raw-apex intensity.

For sparse centroided traces, the median absolute deviation is often zero.  In
that case SNR is deliberately left blank rather than reported as infinity;
quantitative SNR needs a prespecified off-peak or blank-derived noise model.
The reported `Rs` uses measured discrete FWHM, whereas the supplied manuscript
used fixed 0.5-minute widths for one labeling rule, so those values must not be
treated as identical definitions.

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
  --output MSAI\results\peak_picker_adaptive_calibration.json `
  --adaptive
```

Apply the experimental profile while exporting only the uncertain rows:

```powershell
python MSAI\python\ms1_peak_picker.py `
  --peaklist MSAI\peaklist\merged_result_IG_EASMSV1.csv `
  --raw MSAI\data\IG_STD500nM_EASMSV1_2.mzXML `
  --output MSAI\results\merged_result_IG_EASMSV1_auto_rt_aware.csv `
  --review-output MSAI\results\merged_result_IG_EASMSV1_manual_review.csv `
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
  --result-dir MSAI\results\batch_A01_A11_calibrated `
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
  --output-dir MSAI\results\sensitivity_EASMSV1
```

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
  --output MSAI\results\peak_picker_calibration.json
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
