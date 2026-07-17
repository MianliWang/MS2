# MS2 diagnostic standard (current primary workflow)

The current standard is an explicit, conservative screening rule rather than
a finished biological gold standard.  It asks one question: do independently
extracted fragment spectra at two supplied LC peaks support the same compound?
It does not call R/S configuration and does not infer enrichment.

The current versioned machine-readable definition is
[`standards/ms2_diagnostic_standard_v2.json`](../standards/ms2_diagnostic_standard_v2.json).
Version 2 changes the primary RT half-window from 8 to 10 seconds to follow the
AdductMLib Methods section; version 1 remains unchanged as the historical RT8
definition.

## Decision sequence

1. Confirm that two retention times, an acquired DIA window, precursor signal,
   and two non-empty fragment spectra exist.
2. Align fragments one-to-one within 0.01 Da after removing peaks below 1% of
   each spectrum's base peak.
3. Require at least six matched fragments.  A row failing here is
   `insufficient_evidence`, even if the cosine from those few ions is high.
4. On the full aligned spectra, require square-root cosine >= 0.7 and at least
   50% explained intensity in each spectrum.  Entropy is reported but is not a
   gate until it is calibrated on this method.
5. Passing rows are `supported_same_compound`; failing rows with enough matched
   ions are `conflicting_spectra`.

These defaults were chosen as defensible initial guardrails and to avoid a few
coincidental low-intensity matches dominating the result.  They have not yet
been tuned against an independent labeled MS2 truth set, so counts produced by
them are screening counts, not estimates of diagnostic accuracy.

The supplied AdductMLib Methods section directly specifies the current
10-second precursor-fragment RT alignment and corroborates the chromatographic-
correlation cutoff (>0.9), but it does not validate the cosine, six-match,
explained-intensity, 2,000-intensity, or 0.01 Da thresholds. The author's 3 ppm
note remains a shadow variant because its scope is not stated in the PDF.
Acquisition and processing provenance is reconciled in
[`acquisition_parameter_reconciliation.md`](acquisition_parameter_reconciliation.md).

The compatibility extractor computes Pearson across the full RT window. A new
experimental `active_support` mode excludes common zero/tail scans and adds
minimum-support plus apex-offset guards. It addresses a failure mode found in
the historical 8-vs-10-second visual audit, but its thresholds are not part of
diagnostic standard v2 and cannot silently change the primary classification.

## How a real standard should be established

Build compound-level positive controls from repeat injections of authentic
standards and negative controls from near-isobaric compounds, coeluting matrix,
wrong-window extractions, and deliberate interference.  Freeze the extraction
parameters first, label mirror spectra without seeing the automatic result,
split compounds (not scans) across training and held-out batches, then select
thresholds for a prespecified error target.  Report the held-out confusion
matrix, confidence intervals, and sensitivity to RT window, fragment tolerance,
coelution threshold, and consensus-scan count.  The current report/queue is the
first step for creating those blinded manual labels.

## Manual issue codes

`NO_ACQUIRED_DIA_WINDOW`, `EMPTY_ONE_OR_BOTH_SPECTRA`,
`TOO_FEW_MATCHED_FRAGMENTS`, `LOW_COSINE_SIMILARITY`,
`LOW_EXPLAINED_INTENSITY`, `FRAGMENT_COUNT_IMBALANCE`,
`OVERLAPPING_DIA_WINDOWS`, `ASYMMETRIC_UNMATCHED_INTENSITY`,
`SHARED_RT_DIA_TARGETS`, and `NEAR_DECISION_BOUNDARY` are diagnostic/audit
reasons.  Only the versioned decision rule changes the automatic status;
review-only flags do not silently change it.

The per-target SVG/PNG workflow, folder taxonomy, and preserved reviewer
annotations are specified in
[`ms2_static_visual_review_plan.md`](ms2_static_visual_review_plan.md) and
[`standards/ms2_manual_review_schema_v1.json`](../standards/ms2_manual_review_schema_v1.json).
