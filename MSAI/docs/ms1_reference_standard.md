# MS1 reference standard (kept separate)

MS1 is temporarily a context provider for the MS2 project.  The primary MS2
workflow consumes reviewed or otherwise supplied `Peak1`/`Peak2` retention
times.  Its MS1 reference status records only whether two, one, or zero RT
coordinates were supplied.

`Peak1`/`Peak2` are experimental-design/legacy candidate coordinates, not
confirmed chromatographic truth until reviewer provenance and adjudication are
available. `IG` is a laboratory machine/run identifier with multiple pools;
the current `IG` column stores that source's colour-coded result. For IG,
`GREEN` means expected usable double peak, `YELLOW` expected usable single
peak, `RED` no peak or multiple peaks/unusable (not a no-peak truth label), and
`CHECK` further human review required. The image exporter preserves the source
machine, raw label, pool ID, and pooled well in metadata and places images in
`source_machine/IG/` folders. IA, IB, and other source machines remain
separate: their raw labels are retained but their colour semantics are not
assumed without confirmation.

The versioned definition is
[`standards/ms1_reference_standard_v1.json`](../standards/ms1_reference_standard_v1.json).
Experimental automatic MS1 peak-picking retains its own configurations under
`MSAI/config/` and its own FWHM, area, SNR, valley, and resolution diagnostics.
Those thresholds do not pass or fail `ms2_diagnostic_status`.

This boundary is intentional: a clean low MS1 peak or a visibly close doublet
can still be submitted to MS2, while MS2 independently reports whether usable
and mutually consistent fragment evidence exists.
