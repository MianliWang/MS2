# EASMSV1 MS2 visual audit (current iteration)

This note records the first stratified visual inspection of the generated
mirror-spectrum gallery.  It is an engineering audit of the current dataset,
not a blinded truth set and not an accuracy claim.

## What was inspected

- all four `P1_conflict` mirror spectra;
- the first boundary/interference page, including overlap, empty-spectrum, and
  exact-six-match examples;
- representative `P3_insufficient` rows with zero, one, and two dominant
  matched fragments;
- representative high-confidence `P4_supported_audit` rows.

The complete sortable queue and all rows with at least one spectrum are in the
generated `MSAI/results/final/ig_easmsv1/report/` report and galleries.

## Observed failure modes

1. Acquisition coverage dominates the total result.  Most non-evaluable
   doublets have no acquired DIA window, so changing cosine cannot rescue them.
2. The four conflict spectra are not merely weak.  They contain shared fragment
   masses but substantial intensity inversion and strong unmatched peaks,
   consistent with co-isolation/background interference or window mislocalization.
3. `XS714599b` and `XS839010b` use the same RT pair and DIA center for adjacent
   precursor targets and display a strikingly similar unmatched background
   pattern.  This is a cross-target interference warning; the report now flags
   shared RT/DIA target groups.
4. Evidence-insufficient rows are heterogeneous.  Some have no credible mirror
   relationship, while `XS837856b` has only two matched peaks but both are
   dominant and visually concordant.  A fixed six-match floor is conservative
   but can discard sparse potentially useful spectra.
5. Exact-six-match supported rows such as `XS840008b` can have very high cosine
   while being sparse and base-peak dominated.  Their score is mathematically
   valid but their evidence breadth is narrow, so they remain boundary reviews.
6. Some supported spectra show asymmetric unmatched intensity despite passing
   cosine and the 50% explained-intensity floor (for example `XS841031b`).  A
   new audit flag captures this without changing the automatic pass/fail rule.
7. High-confidence supported examples (`XS635896b`, `XS836659b`,
   `XS840767b`) show the expected dense or clean mirror symmetry, demonstrating
   that the visualization separates convincing cases from interference cases.

## Next iteration priorities

1. Resolve acquisition/metadata failures before similarity tuning: verify ion
   mode and adduct, target m/z, real isolation bounds, and RT alignment.
2. Blind-label P1-P3 plus a random P4 audit sample.  Freeze these labels before
   any parameter search.
3. Compare multi-scan consensus against the current single-candidate spectrum;
   require stability across scans or replicate injections for sparse spectra.
4. Calibrate a compound-grouped rule that considers matched count, matched
   fraction, explained-intensity asymmetry, entropy, and recurrent background
   ions.  Do not lower the six-match floor globally from this visual review.
5. Test whether RT-window adjustment or fragment coelution/deconvolution removes
   the recurrent background seen for targets sharing an RT/DIA group.
