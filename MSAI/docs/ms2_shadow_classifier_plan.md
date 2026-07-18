# MS2 same-compound shadow classifier plan

## Purpose and current eligibility

The classifier asks one bounded question: given two spectra that were already
extracted at the supplied Peak1/Peak2 retention times, how strongly do those
spectra support the same-compound hypothesis?  It does not assign R/S, prove a
chiral doublet, or estimate enantioselective enrichment.

The current RT10 coverage result contains 129 rows and about 58 rows with two
non-empty spectra.  It has no independent identity truth and represents one raw
acquisition batch.  It is therefore **feasibility-only**: it may be used to test
feature extraction, abstention, manifests, and shadow-report generation, but it
is not eligible for fitting, threshold selection, validation, or accuracy
claims.  No v3 diagnostic standard is introduced by this work.

## Parameter ownership

Parameters are separated by scientific role.  Training code must not blur
these layers.

| Layer | Parameters | Source and rule |
| --- | --- | --- |
| Paper Method, locked | RT half-window = 10 s; fragment chromatographic Pearson correlation strictly > 0.9; acquisition/DIA conditions | Wang et al., *Orthogonal molecular annotation in mass spectrometry with AdductMLib*, Methods page 17 (manuscript lines 337-339 and 345-348 for the processing values). These describe an experimental method and are never ML features or tuning variables. |
| Extraction, frozen for this model family | precursor EIC = 10 ppm; fragment EIC = 10 ppm; fragment alignment = 0.01 Da; spectrum relative-intensity filter = 1%; minimum fragment intensity = 2,000; exclusive one-to-one matching; close-peak anti-overlap rule; consensus scan count = 1 | Current v2 implementation/default profile. Some are project empirical defaults rather than paper values, but changing them changes the input-generating process. They must be studied as a separately versioned extraction experiment, never silently optimized by this classifier. |
| Trainable decision layer | cosine threshold; matched-fragment threshold; explained-intensity threshold; optional entropy threshold; Logistic Regression coefficients, L2 `C`, class weighting, and probability abstention boundaries | Learned only from an independently labelled, batch-separated truth set. |
| Not allowed in this model family | RT/Pearson/acquisition tuning; L1; trees; interactions; neural networks; current v2 status fields as features | Out of scope by design. |

The source-sidecar analysis mapping is closed-world.  It may contain only the
25 frozen input-generation keys plus the four downstream v2 decision
guardrails (`min_cosine`, `min_matched_peaks`, `min_explained_intensity`, and
`min_entropy_similarity`).  The downstream four may vary but must appear as a
complete group; any unknown future extraction key requires an explicit
model-family contract update.  If both metadata parameter views are present,
they must be identical.

The paper profile is contextual rather than proof that it generated the current
run.  The cited method used a 6-minute C18 experiment, nominal 5 m/z DIA windows
with 1 m/z overlap (and separately reports 4.5 Da preprocessing windows),
stepped HCD 15/30/60, and three files per sample.  The current EASMSV1 raw run
is a 30-minute chiral acquisition with different observed settings.  The
classifier reports this provenance discrepancy and never edits either layer.

## Truth contract

The annotation column accepts exactly four values:

- `positive_same_compound`
- `negative_different_or_interference`
- `uncertain`
- `not_evaluable`

Only the first two are eligible for fitting.  `uncertain` and `not_evaluable`
are retained in the manifest and excluded from fitting.  A positive or negative
label on a row without two non-empty spectra is rejected as a contradiction;
missing bilateral spectra are never synthesized into negative controls.

Labels must be assigned independently of `ms2_diagnostic_status`,
`enantiomer_pair_status`, `ms2_issue_codes`, `IG`, or any existing
supported/conflicting rule output.  Those fields are denied as features and are
used only after inference to construct the shadow disagreement queue.

## Fixed features

All features are recomputed directly from `peak_a_MS2` and `peak_b_MS2` through
the existing exclusive one-to-one spectrum matcher, using the frozen 0.01 Da
alignment and 1% filter:

1. `cosine`: square-root-intensity cosine on the complete filtered spectra.
2. `entropy_similarity`: spectral-entropy similarity from the same matched spectra.
3. `log1p_matched_peaks`: `log(1 + matched fragment count)`.
4. `min_explained_intensity`: the smaller matched-intensity fraction across the two spectra.
5. `explained_intensity_asymmetry`: absolute difference between those two fractions.
6. `fragment_count_balance`: smaller filtered fragment count divided by the larger count.

No identifier, existing status, issue code, manually assigned label, RT-window
setting, or acquisition setting enters the feature vector.

## Development and locked evaluation design

One independently acquired batch is named in advance and locked for the final
evaluation.  Every compound appearing in that batch is purged from development,
even if it also appears in another batch.  Development uses leave-one-batch-out
outer folds; within every outer training partition, tuning uses another set of
leave-one-batch-out inner folds.  At every boundary both conditions must hold:

- validation/test batches do not occur in the corresponding training set;
- validation/test `Compound_ID` values do not occur in the corresponding training set.

There is no row-random fallback.  A dataset with too few independent batches,
one class in a required partition, missing batch/compound identifiers, or no
locked batch fails closed.  The artifact records seeds, row IDs, batch and
compound groups, fold roles and purged rows, absolute input paths, and SHA-256
hashes.

A multi-batch truth CSV must also provide a `batch_sources` manifest with
exactly one source-metadata sidecar for every eligible `batch_id`.  Each sidecar
is checked independently against the same frozen extraction/v2 preprocessing
contract; a single compliant run cannot stand in for the other batches.  The
model artifact records only each resolved path, SHA-256 and proof type, never a
copy of Method or extraction values.

## Prespecified models and selection rule

Baseline A extends the existing threshold grid:

- cosine: 0.50 to 0.95 in steps of 0.025;
- matched fragments: 2, 3, 4, 5, 6, 8, 10;
- minimum explained intensity: 0.3 to 0.8;
- optional entropy: disabled or 0.5 to 0.9.

Baseline B is `StandardScaler` followed by L2 Logistic Regression.  Its only
model choices are `C` in `1e-4 ... 1e2` and class weighting in
`none/balanced`.  The deterministic implementation serializes scaler means and
scales, coefficients, intercept, convergence information, and feature order;
it does not pickle executable objects.
All 14 configurations are attempted in each tuning run.  A configuration that
fails to converge in any inner fold is excluded and recorded; at least one
configuration must remain.  Outer and final refits fail closed on
non-convergence, and deployment rejects artifacts unless the fixed optimizer
records `converged=true` with 1--200 iterations.

None of those search ranges comes from the paper.  They are prespecified
project engineering/empirical defaults that must be revisited only with an
independent truth set.  The same applies to the minimum specificity target
(`0.95`), ten equal-width calibration bins, 1,000 compound-group bootstrap
replicates at 95% confidence, base seed `20260129`, Laplace `(+1)/(+2)` rule-
bucket smoothing, and the `0.05` near-boundary review flag.  The Logistic
optimizer is fixed rather than tuned: deterministic Newton/IRLS,
`max_iter=200`, `tolerance=1e-9`, at most 40 step-halving attempts, and a
`1e-4` sufficient-decrease constant.  Actual seeds and derived bootstrap seeds
are written into the artifact.

Within training data, the support boundary is selected to maximize recall
subject to specificity >= 95%.  Logistic Regression also learns a lower
conflict boundary; probabilities between the conflict and support boundaries
abstain as `insufficient_evidence`.  Missing or invalid feature inputs are
`not_evaluable`.  The threshold baseline remains selected on ties and whenever
Logistic Regression does not improve recall consistently across held-out
development batches while maintaining the specificity constraint.  The locked
batch is opened once for reporting, not for model selection.

Reported evaluation includes confusion counts, precision, recall, specificity,
average precision/PR-AUC, F0.5, Brier score, calibration bins, abstention
coverage, and deterministic group-bootstrap 95% confidence intervals for
all-row calls, emitted-probability quality, selective non-abstained performance,
status rates, and coverage.  A run
without an independent locked batch may not publish these as validation
metrics.

## Shadow contract and review artifacts

Inference preserves every v2 column and appends exactly:

- `ml_same_compound_probability`
- `ml_diagnostic_status`
- `ml_model_id`
- `ml_abstention_reason`
- `ml_review_flags`

Allowed ML statuses are `supported_same_compound`, `conflicting_spectra`,
`insufficient_evidence`, and `not_evaluable`.  The existing v2 result remains
authoritative.  A feasibility run without a trained artifact emits a blank
probability, `not_evaluable`, an untrained model sentinel, and explicit
`NO_INDEPENDENT_TRUTH_MODEL`/`FEASIBILITY_ONLY` flags.

The shadow review writes a combined CSV, a v2/ML-disagreement CSV and HTML,
summary/manifest JSON, and paired mirror-spectrum PNG/SVG for every
disagreement.  It also writes empty-but-valid comparison artifacts when no
disagreement exists.  Model, source table, source metadata, v2 standard, and
method profile fingerprints are carried as provenance; none is merged back
into the source analysis.

The production loader deeply validates the fixed grids, >=95% specificity
floor, nested row/batch/compound partitions, locked-final invariants, and exact
per-batch source hashes before accepting trained calls.  `ml_model_id` is a
canonical content identifier and tamper check, not a cryptographic assertion of
who trained the model.  A deployment whose threat model includes maliciously
fabricated but internally consistent JSON must additionally use an externally
controlled signed model registry or model-ID allowlist.

## Implementation and acceptance gates

1. Lock the truth vocabulary, frozen-source signature, feature allowlist,
   forbidden-feature denylist, and single-section artifact schema.
2. Implement deterministic features, grouped/nested splitting, metrics,
   threshold search, Logistic Regression, abstention, and conservative model
   selection.
3. Add shadow CSV inference while byte-for-byte preserving the logical v2 row
   content.
4. Reuse the existing MS2 mirror-spectrum preparation and static exporter for
   disagreement PNG/SVG and HTML review.
5. Test parameter locks, artifact shape, leakage resistance, label/missing-data
   handling, two-axis group isolation, determinism, known metric vectors,
   model-selection ties, shadow integration, and zero-disagreement output.
6. Run focused tests, the relevant full suite, Ruff, and Pyright.  Do not commit
   or push as part of this implementation gate.

The implementation becomes eligible for a real shadow evaluation only after a
separate truth set supplies authentic same-compound positives, deliberately
designed different/interference negatives, at least several independent
batches, and a prespecified locked final batch.
