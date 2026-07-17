# Acquisition and processing parameter reconciliation

This note separates three evidence layers that must not overwrite one another:

1. `raw_observed`: attributes and real DIA bounds read from the supplied raw
   file;
2. `reference_method`: the AdductMLib ChemRxiv Methods section and the author's
   supplementary lab notes;
3. `analysis_parameters`: the values that actually produced a result CSV.

The machine-readable reference profile is
[`config/acquisition_profiles/adductmlib_chemrxiv_20260129_v1.json`](../config/acquisition_profiles/adductmlib_chemrxiv_20260129_v1.json).

## What the reference method establishes

The manuscript describes an Orbitrap Exploris 240 method with positive full MS
from m/z 150–520 at 60,000 resolution, DIA at 22,500 resolution, stepped HCD
15/30/60%, one m/z overlap, seven sequential spectra per loop, and three DIA
files per sample. Its preprocessing associates precursor and fragment traces
within 10 seconds and requires chromatographic correlation greater than 0.9.

The manuscript is internally inconsistent about DIA width: page 16 says a
5 m/z isolation window with 1 m/z overlap, while page 17 says a fixed 4.5 Da
window. They are preserved as distinct fields. The PDF does not define a
3 ppm extraction tolerance; 3 ppm is recorded separately as an author/lab note
whose scope still needs confirmation.

## What the current EASMSV1 raw file establishes

`IG_STD500nM_EASMSV1_2.mzXML` is a different 30-minute chiral acquisition, not
the manuscript's 6-minute C18 AdductMLib method. Across the complete raw file:

- 3,064 MS1 and 21,441 MS2 scans are present (about seven MS2 scans per MS1);
- every MS2 scan reports HCD, `collisionEnergy=35.0`, and
  `windowWideness=15.0`;
- the 14 observed window centers are 375, 389, 403, …, 557;
- the centers are 14 m/z apart, so the 15 m/z windows also overlap by 1 m/z;
- real coverage is 367.5–564.5 m/z.

The common instrument, HCD, one-m/z-overlap, and seven-spectrum loop logic are
consistent at a design level. Window width, collision energy, chromatography,
and run time show that these are different acquisition versions. Visualizations
therefore label 35/15 as `raw observed` and 15/30/60 plus 5/4.5 as an unconfirmed
`reference method`.

## Why so many rows are currently not evaluable

Of the 351 not-evaluable rows, 306 have no window in the supplied raw file and
all 306 are below its 367.5 m/z lower bound. The manuscript adds an important
new possibility: it says that three DIA files were collected per sample, while
the workspace contains only one EASMSV1 file and its name ends in `_2`.

This does not prove that `_1` and `_3` are missing companion mass segments;
the number may instead denote a replicate or method version. It does mean the
306 rows should be described as **not covered by the currently imported raw
file**, not as failed spectrum matching. Before tuning cosine or fragment
thresholds, verify whether matching `_1` and `_3` files exist and whether their
window ranges should be combined for this peak list.

## Parameter consequences

| Parameter | Current analysis | New evidence | Action |
| --- | ---: | ---: | --- |
| MS1 EIC tolerance | 3 ppm | author lab note: 3 ppm | already aligned |
| MS2 shared EIC tolerance | 10 ppm | lab note: 3 ppm, scope unknown | shadow comparison; do not overwrite |
| precursor–fragment RT half-window | 10 s | reference alignment within 10 s | adopted as primary; close doublets still use the anti-overlap guard |
| chromatographic correlation | >0.9 | reference >0.9 | directly corroborated |
| DIA fallback width | 15 m/z | reference 5/4.5 | keep as fallback only; current real bounds come from raw |
| raw DIA width | 15 m/z | raw scan metadata | authoritative for this file |
| collision energy | raw field 35 | reference stepped 15/30/60 | retain both provenance layers |

The existing one-factor sweep shows that changing the shared MS2 EIC tolerance
from 10 to 3 ppm reduces rows with two spectra from 53 to 45 and changes 19
diagnostic statuses. After separating the two operations, changing only the
precursor EIC to 3 ppm while retaining 10 ppm fragment EICs reduces two-sided
spectra further to 43 and changes 17 statuses. This is not evidence for adopting
3 ppm in MS2; it is evidence that the note's scope must be confirmed.

The historical full-table audit that increased only the RT half-window from 8
to 10 seconds raised two-sided spectra
from 53 to 58 and changes 11 statuses: six not-evaluable rows become
insufficient, two insufficient rows and one conflict become supported, one
insufficient row becomes conflict, and one becomes not evaluable. The focused
queue and SVG/PNG gallery are under
`results/development/ms2/sensitivity_EASMSV1_provenance/rt10_status_change_review_*`. Visual
inspection shows that several promoted rows have MS1 FWHM values of roughly
14–34 seconds, so an 8-second half-window can truncate a substantial part of
their chromatographic profile. Following the instruction to use the paper
Method as the parameter source, 10 seconds is now the primary default; the
changed rows still require manual review and must not be treated as truth.

The 10-second run changes the retained spectrum in 71 rows even though only 11
statuses flip. A complete PNG/SVG queue is available at
`results/development/ms2/sensitivity_EASMSV1_provenance/rt10_all_evidence_change_review_gallery`;
its `shadow_comparison/` views separate 60 spectrum-only changes from 11 status
changes and rank two strong promotion candidates separately.

Mechanistic audit showed that nearly all RT-window flips retain the same apex
scan and candidate fragments; the full-window Pearson filter alone changes
which fragments pass. To test common-zero inflation, the code now offers an
`active_support` correlation mode. It calculates Pearson only where the
precursor is at least 5% of its apex and can require minimum scan support,
consecutive fragment support, and bounded apex displacement. A permissive
shadow (`max apex offset=3 scans`, `min consecutive=2`) retains the two strong
RT10 examples and removes the obvious one-fragment wide-tail artifact, but it
also changes 12 statuses relative to full-window RT10 and creates four new
supported candidates. These are queued under
`active_support_status_change_review_gallery`; the mode is not a production
replacement without blinded labels.

## Data request that unblocks the largest gain

Please confirm whether corresponding EASMSV1 `_1` and `_3` raw files exist,
their DIA centers/bounds, and whether all three files represent one sample.
If so, the next code change should support multi-file DIA routing/merging by
real isolation bounds. That can make currently uncovered rows evaluable;
similarity-threshold tuning cannot recover spectra that were never imported.
