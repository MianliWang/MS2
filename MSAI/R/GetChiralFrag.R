#' function to extract fragment summaries for two chiral LC peaks
#'
#' @param mz_tol
#' @param DIAisowin
#' @param peak_a_rt_window
#' @param peak_b_rt_window
#' @param rt_unit
#' @param mz_tol_unit
#' @export
GetChiralFrag <- function(
  mz_tol,
  DIAisowin,
  peak_a_rt_window,
  peak_b_rt_window,
  rt_unit = "min",
  mz_tol_unit = "legacy_fraction"
) {
  path <- getwd()

  path.data <- paste0(path, "/data")
  path.peak <- paste0(path, "/peaklist")
  path.results <- paste0(path, "/results")

  peak_a_rt_window <- .normalize_rt_window(peak_a_rt_window, rt_unit)
  peak_b_rt_window <- .normalize_rt_window(peak_b_rt_window, rt_unit)
  .mz_bounds(1, mz_tol, mz_tol_unit)

  if (peak_a_rt_window[2] > peak_b_rt_window[1] &&
      peak_b_rt_window[2] > peak_a_rt_window[1]) {
    warning(
      "Chiral RT windows overlap; MVP does not perform deconvolution.",
      call. = FALSE
    )
  }

  setwd(path.peak)
  peakfiles <- list.files()
  if (length(peakfiles) < 1) {
    warning(
      "No peaklist files found in ", path.peak, "; returning NULL.",
      call. = FALSE
    )
    setwd(path)
    return(NULL)
  }

  setwd(path.data)
  msfiles <- list.files()

  for (k in seq_along(peakfiles)) {
    setwd(path.peak)
    print(c("getchiralfragment...", k))
    mycpd <- readxl::read_excel(peakfiles[k])
    mycpd <- .initialize_chiral_columns(mycpd)

    setwd(path.data)
    temp <- strsplit(peakfiles[k], ".xlsx")[[1]]
    temp <- temp[1]
    raw_index <- grep(temp, msfiles)

    if (length(raw_index) == 0) {
      warning(
        "No matching raw MS data found for peaklist ",
        peakfiles[k],
        "; skipping.",
        call. = FALSE
      )
      next
    }
    if (length(raw_index) > 1) {
      warning(
        "Multiple raw MS data files matched peaklist ",
        peakfiles[k],
        "; using first match: ",
        msfiles[raw_index[1]],
        call. = FALSE
      )
      raw_index <- raw_index[1]
    }

    xraw <- xcms::xcmsRaw(msfiles[raw_index], includeMSn = TRUE)
    precursor <- preclist(xraw)

    for (j in seq_len(nrow(mycpd))) {
      mz <- mycpd$mz[j]
      DIAwin <- which(abs(mz - precursor) <= (0.5 * DIAisowin))
      if (length(DIAwin) < 1) {
        mycpd <- .assign_chiral_result(
          mycpd,
          j,
          "peak_a",
          .empty_chiral_result("no_matching_DIA_window")
        )
        mycpd <- .assign_chiral_result(
          mycpd,
          j,
          "peak_b",
          .empty_chiral_result("no_matching_DIA_window")
        )
        next
      }
      DIAwin <- precursor[DIAwin[1]]

      peak_a <- .extract_fragments_for_rt_window(
        xraw,
        j,
        mycpd,
        mz_tol,
        mz_tol_unit,
        DIAwin,
        peak_a_rt_window
      )
      peak_b <- .extract_fragments_for_rt_window(
        xraw,
        j,
        mycpd,
        mz_tol,
        mz_tol_unit,
        DIAwin,
        peak_b_rt_window
      )

      mycpd <- .assign_chiral_result(mycpd, j, "peak_a", peak_a)
      mycpd <- .assign_chiral_result(mycpd, j, "peak_b", peak_b)
    }

    setwd(path.results)
    write.table(
      mycpd,
      file = paste(c(temp, "_chiral_MS2.csv"), collapse = ""),
      sep = ",",
      row.names = FALSE
    )
  }

  setwd(path)
  return(mycpd)
}

.normalize_rt_window <- function(rt_window, rt_unit) {
  if (!is.character(rt_unit) || length(rt_unit) != 1 ||
      !rt_unit %in% c("min", "sec")) {
    stop("rt_unit must be either 'min' or 'sec'.", call. = FALSE)
  }
  if (!is.numeric(rt_window) || length(rt_window) != 2) {
    stop("rt_window must be a numeric vector of length 2.", call. = FALSE)
  }
  if (any(is.na(rt_window)) || any(!is.finite(rt_window))) {
    stop("rt_window must contain finite non-NA values.", call. = FALSE)
  }
  if (rt_window[1] >= rt_window[2]) {
    stop("rt_window start must be smaller than rt_window end.", call. = FALSE)
  }
  if (rt_unit == "min") {
    return(rt_window * 60)
  }
  return(rt_window)
}

.mz_bounds <- function(mz, mz_tol, mz_tol_unit) {
  if (!is.character(mz_tol_unit) || length(mz_tol_unit) != 1) {
    stop(
      "mz_tol_unit must be one of 'legacy_fraction', 'ppm', or 'Da'.",
      call. = FALSE
    )
  }
  if (!is.numeric(mz) || length(mz) != 1 || is.na(mz) || !is.finite(mz)) {
    stop("mz must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(mz_tol) ||
      length(mz_tol) != 1 ||
      is.na(mz_tol) ||
      !is.finite(mz_tol) ||
      mz_tol < 0) {
    stop("mz_tol must be a finite non-negative numeric scalar.", call. = FALSE)
  }

  if (mz_tol_unit == "legacy_fraction") {
    delta <- mz * mz_tol
  } else if (mz_tol_unit == "ppm") {
    delta <- mz * mz_tol / 1e6
  } else if (mz_tol_unit == "Da") {
    delta <- mz_tol
  } else {
    stop(
      "mz_tol_unit must be one of 'legacy_fraction', 'ppm', or 'Da'.",
      call. = FALSE
    )
  }

  return(c(mz - delta, mz + delta))
}

.summarize_xic_peak <- function(peaks) {
  if (is.null(peaks) || is.null(peaks$intensity)) {
    return(.empty_xic_summary("empty_xic"))
  }

  intensity <- suppressWarnings(as.numeric(peaks$intensity))
  if (length(intensity) < 1 || all(is.na(intensity) | !is.finite(intensity))) {
    return(.empty_xic_summary("empty_xic"))
  }

  rt <- .extract_eic_rt(peaks, length(intensity))
  flags <- character()
  valid_intensity <- !is.na(intensity) & is.finite(intensity)
  apex_index <- which.max(replace(intensity, !valid_intensity, -Inf))

  if (length(rt) == length(intensity)) {
    valid_rt <- !is.na(rt) & is.finite(rt)
  } else {
    valid_rt <- rep(FALSE, length(intensity))
  }

  if (!any(valid_rt)) {
    flags <- c(flags, "missing_rt")
    apex_rt <- NA_real_
    rt_start <- NA_real_
    rt_end <- NA_real_
    area <- NA_real_
  } else {
    apex_rt <- rt[apex_index]
    rt_start <- min(rt[valid_rt])
    rt_end <- max(rt[valid_rt])
    integrate_index <- valid_intensity & valid_rt
    if (sum(integrate_index) < 2) {
      flags <- c(flags, "insufficient_points")
      area <- NA_real_
    } else {
      rt_integrate <- rt[integrate_index]
      intensity_integrate <- intensity[integrate_index]
      ord <- order(rt_integrate)
      rt_integrate <- rt_integrate[ord]
      intensity_integrate <- intensity_integrate[ord]
      area <- sum(
        diff(rt_integrate) *
          (head(intensity_integrate, -1) + tail(intensity_integrate, -1)) / 2
      )
    }
  }

  if (length(flags) < 1) {
    flags <- "ok"
  }

  return(list(
    apex_rt = apex_rt,
    apex_intensity = intensity[apex_index],
    area = area,
    rt_start = rt_start,
    rt_end = rt_end,
    quality_flags = paste(flags, collapse = ";")
  ))
}

.format_fragment_string <- function(fragment_mz, intensity) {
  if (length(fragment_mz) < 1 || length(intensity) < 1) {
    return(NA_character_)
  }

  valid <- !is.na(fragment_mz) & !is.na(intensity)
  if (!any(valid)) {
    return(NA_character_)
  }

  return(paste0(fragment_mz[valid], ",", intensity[valid], collapse = ";"))
}

.extract_fragments_for_rt_window <- function(
  xraw,
  index,
  mycpd,
  mz_tol,
  mz_tol_unit,
  DIAwin,
  rt_window_sec
) {
  precurmz <- mycpd$mz[index]
  DIAdata <- ms2copy(xraw, DIAwin)

  if (length(DIAdata@scantime) < 1) {
    return(.empty_chiral_result("no_ms2_scans"))
  }

  mzrange <- DIAdata@mzrange
  mz_bounds <- .mz_bounds(precurmz, mz_tol, mz_tol_unit)
  mzmin <- max(mzrange[1], mz_bounds[1])
  mzmax <- min(mzrange[2], mz_bounds[2])
  if (!is.finite(mzmin) || !is.finite(mzmax) || mzmin >= mzmax) {
    return(.empty_chiral_result("mz_out_of_range"))
  }

  rtmin <- max(min(DIAdata@scantime), rt_window_sec[1])
  rtmax <- min(max(DIAdata@scantime), rt_window_sec[2])
  if (!is.finite(rtmin) || !is.finite(rtmax) || rtmin >= rtmax) {
    return(.empty_chiral_result("rt_window_out_of_range"))
  }

  scan_in_window <- DIAdata@scantime >= rtmin & DIAdata@scantime <= rtmax
  scan_in_window[is.na(scan_in_window)] <- FALSE
  ms2_count <- as.integer(sum(scan_in_window))
  if (ms2_count < 1) {
    return(.empty_chiral_result("no_ms2_scans", ms2_count = 0))
  }

  peaks <- xcms::rawEIC(
    DIAdata,
    mzrange = cbind(mzmin, mzmax),
    rtrange = cbind(rtmin, rtmax)
  )
  peaks <- .attach_scan_rt(peaks, DIAdata@scantime)
  peak_summary <- .summarize_xic_peak(peaks)
  flags <- .split_quality_flags(peak_summary$quality_flags)

  if (is.null(peaks$scan) || is.null(peaks$intensity)) {
    return(.merge_chiral_result(
      peak_summary,
      MS2 = NA_character_,
      ms2_count = ms2_count,
      flags = c(flags, "no_fragment_scan")
    ))
  }

  intensity <- suppressWarnings(as.numeric(peaks$intensity))
  if (length(intensity) < 1 || all(is.na(intensity) | !is.finite(intensity))) {
    return(.merge_chiral_result(
      peak_summary,
      MS2 = NA_character_,
      ms2_count = ms2_count,
      flags = flags
    ))
  }

  valid_intensity <- !is.na(intensity) & is.finite(intensity)
  scan.max <- peaks$scan[which.max(replace(intensity, !valid_intensity, -Inf))]
  if (is.na(scan.max) || scan.max >= length(DIAdata@scanindex)) {
    return(.merge_chiral_result(
      peak_summary,
      MS2 = NA_character_,
      ms2_count = ms2_count,
      flags = c(flags, "invalid_apex_scan")
    ))
  }

  scanNum <- c(DIAdata@scanindex[scan.max], DIAdata@scanindex[scan.max + 1])
  correctindex <- (scanNum[1] + 1):scanNum[2]
  mz.frag <- DIAdata@env$mz[correctindex]
  if (length(mz.frag) < 1) {
    return(.merge_chiral_result(
      peak_summary,
      MS2 = NA_character_,
      ms2_count = ms2_count,
      flags = c(flags, "no_candidate_fragments")
    ))
  }

  fragment_index <- which(mz.frag < precurmz - 10)
  if (length(fragment_index) < 1) {
    return(.merge_chiral_result(
      peak_summary,
      MS2 = NA_character_,
      ms2_count = ms2_count,
      flags = c(flags, "no_candidate_fragments")
    ))
  }
  mz.frag <- mz.frag[fragment_index]

  native.peak <- intensity
  fragment_mz <- numeric()
  fragment_intensity <- numeric()
  for (k in seq_along(mz.frag)) {
    mz0 <- mz.frag[k]
    frag_bounds <- .mz_bounds(mz0, mz_tol, mz_tol_unit)
    frag_mzmin <- max(mzrange[1], frag_bounds[1])
    frag_mzmax <- min(mzrange[2], frag_bounds[2])
    if (!is.finite(frag_mzmin) ||
        !is.finite(frag_mzmax) ||
        frag_mzmin >= frag_mzmax) {
      next
    }

    frag.peak <- xcms::rawEIC(
      DIAdata,
      mzrange = cbind(frag_mzmin, frag_mzmax),
      rtrange = cbind(rtmin, rtmax)
    )
    frag.peak <- suppressWarnings(as.numeric(frag.peak$intensity))
    if (length(frag.peak) != length(native.peak)) {
      next
    }
    if (!any(is.finite(frag.peak))) {
      next
    }
    native_sd <- sd(native.peak, na.rm = TRUE)
    frag_sd <- sd(frag.peak, na.rm = TRUE)
    if (!is.finite(native_sd) ||
        !is.finite(frag_sd) ||
        native_sd == 0 ||
        frag_sd == 0) {
      next
    }
    frag_max <- max(frag.peak, na.rm = TRUE)
    if (frag_max < 2000) {
      next
    }

    corr <- suppressWarnings(cor(native.peak, frag.peak, use = "complete.obs"))
    if (!is.na(corr) && corr > 0.9) {
      fragment_mz <- c(fragment_mz, mz.frag[k])
      fragment_intensity <- c(fragment_intensity, frag_max)
    }
  }

  MS2 <- .format_fragment_string(fragment_mz, fragment_intensity)
  if (is.na(MS2)) {
    flags <- c(flags, "no_fragments")
  }

  return(.merge_chiral_result(
    peak_summary,
    MS2 = MS2,
    ms2_count = ms2_count,
    flags = flags
  ))
}

.initialize_chiral_columns <- function(mycpd) {
  character_columns <- c(
    "peak_a_MS2",
    "peak_a_quality_flags",
    "peak_b_MS2",
    "peak_b_quality_flags"
  )
  numeric_columns <- c(
    "peak_a_apex_rt",
    "peak_a_apex_intensity",
    "peak_a_area",
    "peak_a_ms2_count",
    "peak_b_apex_rt",
    "peak_b_apex_intensity",
    "peak_b_area",
    "peak_b_ms2_count"
  )

  for (column in character_columns) {
    mycpd[[column]] <- rep(NA_character_, nrow(mycpd))
  }
  for (column in numeric_columns) {
    mycpd[[column]] <- rep(NA_real_, nrow(mycpd))
  }
  return(mycpd)
}

.assign_chiral_result <- function(mycpd, row_index, prefix, result) {
  mycpd[[paste0(prefix, "_MS2")]][row_index] <- .as_scalar_character(result$MS2)
  mycpd[[paste0(prefix, "_apex_rt")]][row_index] <- .as_scalar_numeric(result$apex_rt)
  mycpd[[paste0(prefix, "_apex_intensity")]][row_index] <- .as_scalar_numeric(result$apex_intensity)
  mycpd[[paste0(prefix, "_area")]][row_index] <- .as_scalar_numeric(result$area)
  mycpd[[paste0(prefix, "_ms2_count")]][row_index] <- .as_scalar_numeric(result$ms2_count)
  mycpd[[paste0(prefix, "_quality_flags")]][row_index] <- .as_scalar_character(
    result$quality_flags
  )
  return(mycpd)
}

.as_scalar_character <- function(value) {
  if (length(value) < 1 || all(is.na(value))) {
    return(NA_character_)
  }
  value <- as.character(value)
  value <- value[!is.na(value) & value != ""]
  if (length(value) < 1) {
    return(NA_character_)
  }
  return(paste(value, collapse = ";"))
}

.as_scalar_numeric <- function(value) {
  if (length(value) < 1 || is.na(value[1])) {
    return(NA_real_)
  }
  return(as.numeric(value[1]))
}

.attach_scan_rt <- function(peaks, scantime) {
  if (!is.null(peaks$rt)) {
    return(peaks)
  }
  if (is.null(peaks$scan)) {
    return(peaks)
  }

  scan <- suppressWarnings(as.integer(peaks$scan))
  rt <- rep(NA_real_, length(scan))
  valid <- !is.na(scan) & scan >= 1 & scan <= length(scantime)
  rt[valid] <- scantime[scan[valid]]
  peaks$rt <- rt
  return(peaks)
}

.extract_eic_rt <- function(peaks, expected_length) {
  for (field in c("rt", "rtime", "scantime")) {
    if (!is.null(peaks[[field]]) && length(peaks[[field]]) == expected_length) {
      return(suppressWarnings(as.numeric(peaks[[field]])))
    }
  }
  if (!is.null(peaks$scan) && length(peaks$scan) == expected_length) {
    return(suppressWarnings(as.numeric(peaks$scan)))
  }
  return(rep(NA_real_, expected_length))
}

.empty_xic_summary <- function(flag) {
  return(list(
    apex_rt = NA_real_,
    apex_intensity = NA_real_,
    area = NA_real_,
    rt_start = NA_real_,
    rt_end = NA_real_,
    quality_flags = flag
  ))
}

.empty_chiral_result <- function(flag, ms2_count = NA_integer_) {
  return(list(
    MS2 = NA_character_,
    apex_rt = NA_real_,
    apex_intensity = NA_real_,
    area = NA_real_,
    ms2_count = ms2_count,
    quality_flags = flag
  ))
}

.merge_chiral_result <- function(peak_summary, MS2, ms2_count, flags) {
  flags <- unique(flags[!is.na(flags) & flags != ""])
  if (length(flags) < 1) {
    flags <- "ok"
  }

  return(list(
    MS2 = MS2,
    apex_rt = peak_summary$apex_rt,
    apex_intensity = peak_summary$apex_intensity,
    area = peak_summary$area,
    ms2_count = ms2_count,
    quality_flags = paste(flags, collapse = ";")
  ))
}

.split_quality_flags <- function(quality_flags) {
  if (is.null(quality_flags) || is.na(quality_flags) || quality_flags == "") {
    return(character())
  }
  if (quality_flags == "ok") {
    return(character())
  }
  return(unlist(strsplit(quality_flags, ";", fixed = TRUE)))
}
