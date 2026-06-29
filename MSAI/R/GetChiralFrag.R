#' Extract fragment summaries for two chiral LC peaks
#'
#' 这个函数是 `GetFrag()` 的手性峰 MVP 版本：同一个 precursor m/z 不再只看
#' 一个 RT 中心点，而是在用户给定的两个 retention time windows 中分别提取
#' XIC peak summary 和 MS2 fragments。
#'
#' 当前 MVP 的文件约定仍然比较接近旧代码：
#' - `peaklist/` 中放 Excel peaklist，且需要小写 `mz` 列。
#' - peaklist 文件名主体需要能匹配 `data/` 中的 raw/mzML/mzXML 文件名。
#' - 输出写入 `results/{sample}_chiral_MS2.csv`，不会覆盖旧 `{sample}_MS2.csv`。
#'
#' @param mz_tol m/z tolerance。默认单位由 `mz_tol_unit` 控制。
#' @param DIAisowin DIA isolation window 宽度，用于把 feature m/z 匹配到 DIA window。
#' @param peak_a_rt_window 第一个手性峰 RT window，长度为 2。
#' @param peak_b_rt_window 第二个手性峰 RT window，长度为 2。
#' @param rt_unit RT window 的单位，`"min"` 会转成秒，`"sec"` 原样使用。
#' @param mz_tol_unit m/z tolerance 单位；默认 `"legacy_fraction"` 保留旧逻辑。
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

  # 工作流总览：
  # 1. 把用户输入的两个 RT windows 统一成秒，因为 xcmsRaw/rawEIC 使用秒。
  # 2. 逐个读取 peaklist，并用文件名匹配对应 raw MS 文件。
  # 3. 对 peaklist 中每个 feature：先按 m/z 找 DIA precursor window。
  # 4. 在 peak_a / peak_b 两个 RT window 内分别提取 XIC、峰顶、面积和 fragments。
  # 5. 把两个峰的结果写成不同前缀列，例如 `peak_a_MS2` 和 `peak_b_MS2`。
  rt_windows <- list(
    peak_a = .normalize_rt_window(peak_a_rt_window, rt_unit),
    peak_b = .normalize_rt_window(peak_b_rt_window, rt_unit)
  )

  # 只做一次 dummy call 来校验 mz tolerance 参数是否合法。
  .mz_bounds(1, mz_tol, mz_tol_unit)

  if (rt_windows$peak_a[2] > rt_windows$peak_b[1] &&
      rt_windows$peak_b[2] > rt_windows$peak_a[1]) {
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

    # 当前 MVP 仍沿用 Excel peaklist；如果输入是 CSV 或列名是 `MZ`，
    # 这里会失败或后续 `mycpd$mz` 为空，需要先转换数据或后续增加读取分支。
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
        for (prefix in names(rt_windows)) {
          mycpd <- .assign_chiral_result(
            mycpd,
            j,
            prefix,
            .empty_chiral_result("no_matching_DIA_window")
          )
        }
        next
      }
      DIAwin <- precursor[DIAwin[1]]

      for (prefix in names(rt_windows)) {
        peak_result <- .extract_fragments_for_rt_window(
          xraw,
          j,
          mycpd,
          mz_tol,
          mz_tol_unit,
          DIAwin,
          rt_windows[[prefix]]
        )
        mycpd <- .assign_chiral_result(mycpd, j, prefix, peak_result)
      }
    }

    setwd(path.results)
    write.table(
      mycpd,
      file = paste0(temp, "_chiral_MS2.csv"),
      sep = ",",
      row.names = FALSE
    )
  }

  setwd(path)
  return(mycpd)
}

.normalize_rt_window <- function(rt_window, rt_unit) {
  # 用户通常用分钟描述 LC peak window；xcmsRaw 对象里的 RT 是秒。
  # 这里把所有合法输入统一成秒，后面的 rawEIC(rtrange=...) 只吃这一种单位。
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
  # 旧 GetFrag 使用 `mz ± mz * mz_tol`，这不是 ppm，也不是 Da。
  # 为了兼容历史调用，默认继续使用 legacy_fraction；ppm/Da 只是新增分支。
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
  # rawEIC() 返回的是一个 EIC/XIC trace：每个 scan 有一个 intensity。
  # 这个 helper 只做 MVP 级别的峰摘要：
  # - apex_rt: intensity 最大点的 RT，输出单位为秒。
  # - apex_intensity: apex 点强度。
  # - area: 对 RT-intensity 做 trapezoidal integration。
  # - quality_flags: 空 trace、缺 RT、点数不足等情况的可读标签。
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
  # 沿用旧输出风格：每个 fragment 写成 "mz,intensity"，多个 fragment 用分号拼接。
  # 这样结果仍然是普通 CSV scalar cell，不会变成 list column。
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
  # 这是核心 extraction helper。它只负责“一个 feature + 一个 RT window”：
  # 1. 用 ms2copy() 把目标 DIA window 的 MS2 scans 复制成 pseudo xcmsRaw object。
  # 2. 在 precursor m/z tolerance 内提取 native XIC，并总结 apex/area。
  # 3. 取 native XIC apex scan 里的候选 fragment m/z。
  # 4. 对每个候选 fragment 再提取 fragment EIC。
  # 5. 保留与 native XIC 高度相关、且强度过阈值的 fragments。
  #
  # 手性区分只发生在 `rt_window_sec`；MS2 本身不负责判断 R/S。
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
  native_sd <- sd(native.peak, na.rm = TRUE)
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
  # 在原 peaklist 后面预先加输出列。每一行 feature 会得到两组结果：
  # peak_a_* 和 peak_b_*。这里显式初始化类型，避免后续写 CSV 时出现 list column。
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
  # 把一个 helper 返回的 list 写回 data frame 的一行。
  # `.as_scalar_*()` 是防御层：即使 helper 意外返回向量，也压成单个 CSV cell。
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
  # CSV 的一个单元格只能安全承载一个字符串；多个 flag 或 fragment 在这里用分号合并。
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
  # 数值输出只取第一个元素，保证 apex/area/ms2_count 都是 scalar。
  if (length(value) < 1 || is.na(value[1])) {
    return(NA_real_)
  }
  return(as.numeric(value[1]))
}

.attach_scan_rt <- function(peaks, scantime) {
  # 不同 xcms/rawEIC 版本返回结构可能略有差异。
  # 如果 EIC 没有 RT，但有 scan index，就用 xcmsRaw@scantime 补出秒单位 RT。
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
  # 兼容不同字段名：优先用真实 RT；没有 RT 时最后才退回 scan number。
  # scan number 只是保守兜底，不等价于真实秒单位 RT；正常路径会由 .attach_scan_rt() 补 RT。
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
  # XIC 层失败时使用统一 NA 结构，让上层可以继续写完整 CSV。
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
  # 整个 peak/window 无法提取时的统一返回结构。
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
  # 把 XIC summary 和 MS2 fragment summary 合并成最终写入 peak_a/peak_b 列的结构。
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
  # quality_flags 在 CSV 中是分号字符串；内部处理时拆回 character vector。
  if (is.null(quality_flags) || is.na(quality_flags) || quality_flags == "") {
    return(character())
  }
  if (quality_flags == "ok") {
    return(character())
  }
  return(unlist(strsplit(quality_flags, ";", fixed = TRUE)))
}
