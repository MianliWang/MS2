# EASMSV1 MS1 parameter and image audit

This is an exploratory audit on the same 525-row development table.  It is not
an independent accuracy claim and MS1 remains reference-only for the MS2
project.

## Data and baseline

- 404 supplied double peaks, 99 supplied single peaks, and 22 rows with no
  supplied peak; 10,255 fully blank trailing rows are discarded.
- The raw file contains 3,064 MS1 scans over 30 minutes with a median 0.60 s
  interval.
- The current adaptive result agrees in class and RT (within 0.1 min) on
  462/525 rows.  Its 63 discrepancies comprise 31 double-peak RT
  mislocalizations, 18 double-to-single calls, three double-to-ambiguous calls,
  ten single-to-double calls, and one single-to-ambiguous call.
- All 31 mislocalized doubles and all 18 missed doubles still locate at least
  one supplied peak.  The dominant problem is second-candidate generation and
  pairing rather than total signal absence.

## Parameter findings

### Keep the 3 ppm EIC window for this run

With the same SG3 adaptive configuration, strict reference agreement was
396/525 at 1 ppm, 456/525 at 2 ppm, 462/525 at 3 ppm, 452/525 at 5 ppm, and
404/525 at 10 ppm.  All 907 supplied apices were covered by 3 ppm; their
absolute centroid error had p95 1.141 ppm, p99 2.423 ppm, and maximum
2.865 ppm.  Narrower windows fragment the trace under mass drift; wider windows
admit more same-mass interference.

### The close-double problem is not fixed by one global threshold

Thirty supplied doubles are separated by <=12 s; the current adaptive picker
strictly localizes only 3/30.  Every missed real second peak remains a local
maximum after SG3, but 29 fail `min_prominence_ratio=0.25` and 11 fail
`min_support_scans=3`.

The new shadow candidate profile lowers only these two candidate gates to 0.03
and two scans.  Before `top_k` truncation it covers 388/404 supplied doubles and
17/30 close doubles, but 57/99 supplied singles also contain two or more
candidates.  It therefore exports open-square review candidates and never
assigns a final class.

### Use agreement with abstention as the safe current iteration

The baseline SG7 and experimental SG3/adaptive models are considered stable
only when they return the same class and every picked RT agrees within 0.05 min.
This accepts 378/525 rows; 372/378 (98.4%) agree with the supplied reference on
this development set.  The other 147 rows are explicitly routed to manual
review.  This improves reliability by abstaining rather than pretending that a
single parameter set solved shoulders and global mislocalization.

### Image-derived failure modes

- `XS061098d` and `XS841653b` visibly contain close shoulder pairs; the
  high-recall candidates mark both, while stricter prominence/support or valley
  gates remove one.
- `XS840844b` contains several local shoulders.  SG3 selects an internal
  fluctuation and becomes RT-mislocalized, while SG7 matches the supplied pair.
- `XS843092a` has a dominant late asymmetric peak plus an unrelated earlier
  same-mass peak.  Global top-two ranking pairs the wrong features.
- `XS843055a` is labeled as a supplied single but visually has a pronounced
  second lobe; this may be a model false positive or a noisy legacy label and
  requires adjudication before training.
- The original seven “low-clean” flags were one-scan spikes.  Requiring at
  least three consecutive half-height scans reduces this false queue to zero.
  The current table therefore contains no validated low-clean example.

### Metric corrections

The EIC is extremely sparse, so global median/MAD noise is zero for nearly all
peaks and the current SNR is blank for most rows.  Future quantitative SNR must
use a local off-peak, off-mass, blank, or replicate-derived noise definition.

The code measures FWHM but previously used the factor-2 baseline-width
resolution equation.  It now reports the Gaussian half-height convention:

```text
Rs = 1.17741 * (RT2 - RT1) / (FWHM1 + FWHM2)
```

The IUPAC factor-2 definition explicitly uses widths at base.  For heavily
overlapped or asymmetric peaks, even the corrected FWHM metric is descriptive
rather than proof of two components.

## Is asymmetric Gaussian / EMG / EGH cascade fitting viable?

Yes, but as a second-stage local ROI model, not as a replacement for SG
candidate generation.

- The exponential-Gaussian hybrid (EGH; possibly the remembered “ESG”) was
  designed as a stable asymmetric chromatographic peak model and is close to
  EMG at modest asymmetry: [Lan and Jorgenson, 2001](https://doi.org/10.1016/S0021-9673(01)00594-5).
- A bi-Gaussian mixture plus statistical model selection has been used to
  quantify and deconvolve asymmetric, overlapping LC-MS peaks, while its paper
  cautions that BIC is heuristic for LC-MS traces:
  [Yu and Peng, 2010](https://doi.org/10.1186/1471-2105-11-559).
- centWave uses ROI construction and continuous wavelets, with optional
  Gaussian fitting, rather than asking one global nonlinear model to discover
  every component: [Tautenhahn et al., 2008](https://doi.org/10.1186/1471-2105-9-504).
- MZmine retains local-minimum, wavelet, and SG approaches; its original paper
  described Gaussian/EMG modeling as experimental and not thoroughly
  validated: [Pluskal et al., 2010](https://doi.org/10.1186/1471-2105-11-395).
- Iterative peak fitting becomes more identifiable when multiple detector
  profiles are fitted together and can remain ambiguous on one-dimensional
  traces: [Erny et al., 2021](https://doi.org/10.3390/separations8100178).

The recommended next experiment is therefore: SG3 high-recall candidates ->
crop only close/shoulder ROIs -> compare bounded one-component and two-component
bi-Gaussian models with multiple starts -> require AICc/BIC improvement,
structured-residual reduction, plausible widths/separation, and bootstrap
stability -> leave unresolved cases as `uncertain`.  MS2 fragment coelution or
replicate injections supply the extra dimension needed when one EIC cannot
distinguish a skewed single peak from two near-coincident peaks.
