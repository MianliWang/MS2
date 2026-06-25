source(file.path("MSAI", "R", "GetChiralFrag.R"))

stopifnot(identical(.normalize_rt_window(c(1, 2), "min"), c(60, 120)))
stopifnot(identical(.normalize_rt_window(c(10, 20), "sec"), c(10, 20)))

stopifnot(all.equal(.mz_bounds(100, 10, "ppm"), c(99.999, 100.001)))
stopifnot(all.equal(.mz_bounds(100, 0.001, "Da"), c(99.999, 100.001)))
stopifnot(all.equal(.mz_bounds(100, 0.00001, "legacy_fraction"), c(99.999, 100.001)))

peaks <- list(
  rt = c(0, 1, 2),
  scan = c(1, 2, 3),
  intensity = c(0, 10, 0)
)
summary <- .summarize_xic_peak(peaks)
stopifnot(identical(summary$apex_rt, 1))
stopifnot(identical(summary$apex_intensity, 10))
stopifnot(all.equal(summary$area, 10))
stopifnot(identical(summary$quality_flags, "ok"))

stopifnot(identical(.format_fragment_string(numeric(), numeric()), NA_character_))
stopifnot(identical(.format_fragment_string(c(50, 75), c(1000, 2000)), "50,1000;75,2000"))
