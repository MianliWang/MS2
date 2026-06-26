# Legacy GetFrag.R 总览：
# 这个文件保留原始的 single-peak DIA MS2 fragment extraction 逻辑。
# 它会读取 peaklist 中每个 feature 的 `mz` / `rt`，在匹配的 raw file 里
# 找对应的 DIA precursor window，并尝试提取与 precursor EIC 共洗脱的 fragments。
# NOTE: 这不是手性双峰提取逻辑；双 RT window 的手性峰 MVP 在 `GetChiralFrag()` 中。
# 本文件中的中文注释只解释 legacy 行为，不改变任何代码路径、阈值或输出格式。
#' function to extract fragment from DIA window
#'
#' @param Library
#' @param path
#' @param msfiles
#' @param mz_tol
#' @export
# GetFrag() 是旧流程的主入口。
# 参数：
# - `mz_tol`: legacy m/z tolerance，旧代码按 `mz +/- mz * mz_tol` 解释。
# - `DIAisowin`: DIA isolation window 宽度，用于匹配 precursor mz 所在窗口。
# - `RTwin`: legacy 参数；当前函数会传给 `getfrag()`，但旧核心逻辑没有真正使用它。
# 目录约定：
# - `peaklist/`: Excel peak list，至少依赖 `mz` 和 `rt` 两列。
# - `data/`: 与 peaklist 文件名前缀匹配的 raw MS/MS DIA 数据。
# - `results/`: 输出 `{sample}_MS2.csv`。
# side effect: 会切换 working directory，并写出结果 CSV。
# 返回值：最后一个处理过的 `mycpd` data frame / tibble；无 peaklist 时返回 NULL。
GetFrag <- function(mz_tol, DIAisowin, RTwin) {
  #' the path to save results
  # `getwd()` 取得当前工作目录；旧流程假设用户先 `setwd(path)` 到 package 项目目录。
  path <- getwd()

  #' the path to save raw data
  # `paste0()` 会把字符串直接拼接起来；这里拼出三个固定工作目录。
  path.data <- paste0(path, "/data")
  path.peak <- paste0(path, "/peaklist")
  path.results <- paste0(path, "/results")

  # `setwd()` 切换目录；`list.files()` 枚举当前目录下的文件名。
  # 这里先进入 `peaklist/`，收集所有 peaklist Excel 文件。
  setwd(path.peak) # the list of peak documents
  peakfiles <- list.files()
  if (length(peakfiles) < 1) {
    # `warning()` 给出可恢复告警；`call. = FALSE` 让信息更短。
    warning(
      "No peaklist files found in ", path.peak, "; returning NULL.",
      call. = FALSE
    )
    setwd(path)
    # `return(...)` 显式结束函数并返回值。
    return(NULL)
  }

  # 再进入 `data/`，收集 raw MS data 文件名。文件类型由 `xcmsRaw()` 支持情况决定。
  setwd(path.data) # the list of raw MS data
  msfiles <- list.files()

  # `fragments` 保存单个 feature 提取到的 fragment string；循环中会反复覆盖。
  fragments <- NULL
  # `seq_along(peakfiles)` 比 `1:length(peakfiles)` 更安全，空向量时不会生成 `1:0`。
  for (k in seq_along(peakfiles)) {

    # get the peak list
    setwd(path.peak)
    print(c("getfragment...", k))
    # `readxl::read_excel()` 读取 Excel peak list；后续代码假设存在 `mz` 和 `rt` 列。
    mycpd <- readxl::read_excel(peakfiles[k])
    # `$` 用于访问 data frame / list 的列或成员；这里新增/覆盖 `MS2` 输出列。
    mycpd$MS2 <- rep(0, nrow(mycpd))

    # get the raw data
    setwd(path.data)
    # `strsplit()` 把 peaklist 文件名按 ".xlsx" 分开，用第一段作为 sample name。
    temp <- strsplit(peakfiles[k], ".xlsx")[[1]]
    temp <- temp[1] # sample name
    # `grep()` 在 raw 文件名列表中找包含 sample name 的文件。
    index <- grep(temp, msfiles)

    if (length(index) == 0) {
      warning(
        "No matching raw MS data found for peaklist ",
        peakfiles[k],
        "; skipping.",
        call. = FALSE
      )
      # `next` 跳过当前 peaklist，继续处理下一个 peaklist 文件。
      next
    } # no MS data
    if (length(index) > 1) {
      # 多个 raw file 匹配时沿用稳定规则：使用 `list.files()` 返回顺序中的第一个。
      warning(
        "Multiple raw MS data files matched peaklist ",
        peakfiles[k],
        "; using first match: ",
        msfiles[index[1]],
        call. = FALSE
      )
      index <- index[1]
    }

    # `xcms::xcmsRaw(..., includeMSn = TRUE)` 读取 raw 数据，并保留 MS/MS 信息。
    # NOTE: 这里依赖 xcms 支持的 raw/mzXML/mzML 等格式；本文件不做格式判断。
    xraw <- xcms::xcmsRaw(msfiles[index], includeMSn = TRUE)

    #'extracting precursor DIA windows
    # 从 raw object 中提取 DIA precursor window 列表。
    precursor <- preclist(xraw)
    # `seq_len(nrow(mycpd))` 按 peaklist 每一行 feature 循环。
    for (j in seq_len(nrow(mycpd))) {
      # `mycpd$mz[j]` 是当前 feature 的 precursor m/z。
      mz <- mycpd$mz[j]
      # 找出与当前 m/z 距离不超过半个 `DIAisowin` 的 DIA precursor window。
      # `which()` 返回满足条件的位置索引。
      DIAwin <- which(abs(mz - precursor) <= (0.5 * DIAisowin))
      if (length(DIAwin) < 1) {
        # 当前 feature 没有匹配 DIA window 时，不写 MS2，直接处理下一行。
        next
      }
      # 如果多个窗口匹配，legacy 行为只取第一个。
      DIAwin <- precursor[DIAwin[1]]

      #' get fragments from each window based correlation >0.9
      # 核心 fragment extraction 在 `getfrag()` 中完成。
      fragments <- getfrag(xraw, j, mycpd, mz_tol, DIAwin, RTwin)
      if (length(fragments) == 0) {
        next
      }
      # 成功时把 fragment string 写入当前 feature 的 `MS2` 列。
      mycpd$MS2[j] <- fragments
    }
    setwd(path.results)
    # 输出 `{sample}_MS2.csv`。`paste(..., collapse = "")` 把 sample name 和后缀拼成文件名。
    write.table(
      mycpd,
      file = paste(c(temp, "_MS2.csv"), collapse = ""),
      sep = ",",
      row.names = FALSE
    )
  }
  setwd(path)
  return(mycpd)
}


#' -------------------
#' Finding DIA windows
#'--------------------
#' @param xmsn
#' @return
# preclist() 从 xcmsRaw S4 object 中提取 unique precursor mz / DIA window 列表。
# R 里 `@` 用于访问 S4 object slot，例如 `xmsn@msnPrecursorMz`；
# `$` 则常用于 list/data frame 的字段，例如 `peaks$intensity`。
# 空 precursor 列表时直接返回空对象，避免后续循环产生运行风险。
preclist <- function(xmsn) {
  #'extracting DIA windows
  # `@msnPrecursorMz` 是 xcmsRaw 中保存 MS/MS precursor mz 的 slot。
  precmz <- xmsn@msnPrecursorMz
  if (length(precmz) < 1) {
    return(precmz)
  }
  # 从第一个 precursor 开始，后面遇到未出现过的 mz 就追加到 `precur`。
  precur <- precmz[1]
  for (i in seq_along(precmz)[-1]) {
    if (length(which(precur == precmz[i])) == 0) {
      precur <- c(precur, precmz[i])
    }
  }
  return(precur)
}


#' -------------------------------
#' This function is used to extract fragmetns from each DIA window
#' -------------------------------
#' @param rawdata
#' @param prec_list
#' @param index2
#' @param mz_tol
#' @param Library
#' @param DIAmzwin
#' @param rtwindow
#'
#' @return
#'
# getfrag() 是旧代码的核心 fragment extraction 函数。
# 输入：
# - `xraw`: 由 `xcmsRaw(..., includeMSn = TRUE)` 读取的 raw object。
# - `index`: 当前 peaklist 行号。
# - `mycpd`: peaklist 表，至少需要 `mz` 和 `rt` 列。
# - `mz_tol`: legacy m/z tolerance，按 `mz +/- mz * mz_tol` 计算。
# - `DIAwin`: 当前 feature 匹配到的 DIA precursor window。
# - `RTwin`: legacy 参数；NOTE: 当前代码没有用它控制 RT window。
# 输出：
# - 成功时返回 `"fragment_mz,max_intensity;fragment_mz,max_intensity"` 格式字符串。
# - 没有 candidate fragment 或没有通过筛选时返回 `NA_character_`。
# - apex scan 落在最后一个 scan 时返回 `NULL`，保持旧行为。
getfrag <- function(xraw, index, mycpd, mz_tol, DIAwin, RTwin) {
  # legacy 假设 peaklist 里的 `rt` 单位是分钟；这里乘 60 转成秒。
  mycpd$rt <- mycpd$rt * 60
  # 当前 feature 的 precursor m/z。
  precurmz <- mycpd$mz[index]

  #' read rawdata for each DIA window
  # `ms2copy()` 把指定 DIA window 的 MS2 scans 复制为 pseudo xcmsRaw object。
  DIAdata <- ms2copy(xraw, DIAwin)
  # `@mzrange` 是 pseudo raw object 中所有 mz 的范围，用于裁剪 EIC mz window。
  mzrange <- DIAdata@mzrange
  minmz <- mzrange[1]
  maxmz <- mzrange[2]
  # precursor EIC 的 mz window：`mz +/- mz * mz_tol`，并限制在原始 mzrange 内。
  mzmin <- max(minmz, precurmz - precurmz * mz_tol)
  mzmax <- min(maxmz, precurmz + precurmz * mz_tol)
  # `@scantime` 是 scan 的 retention time；当前 pseudo object 里来自 MS2 RT。
  rtrange <- DIAdata@scantime
  # NOTE: 这是 legacy behavior：虽然函数有 `RTwin` 参数，实际仍写死为中心 RT +/- 10 秒。
  rtmin <- max(min(rtrange), mycpd$rt[index] - 10)
  rtmax <- min(max(rtrange), mycpd$rt[index] + 10)

  #'extracting precursor ions and peaks
  # `rawEIC()` 提取 extracted ion chromatogram (EIC/XIC)。
  # 这里提取 precursor m/z 在 RT window 内的 EIC，作为后续相关性比较的 native peak。
  peaks <- xcms::rawEIC(
    DIAdata,
    mzrange = cbind(mzmin, mzmax),
    rtrange = cbind(rtmin, rtmax)
  )

  #'finding the scan number of the peak top
  # `which.max()` 返回最大 intensity 的位置；这里用它找 precursor EIC apex。
  scan.max <- which.max(peaks$intensity)
  # `peaks$scan` 是 rawEIC 返回的 scan 编号；从 apex 位置取真正 scan id。
  scan.max <- peaks$scan[scan.max[1]]
  if (scan.max >= length(DIAdata@scanindex)) { # the last scanning point
    return(NULL)
  }
  # `@scanindex` 记录每个 scan 在 mz/intensity 向量里的起点。
  # 通过 apex scan 和下一个 scan 的 index 边界，取出 apex scan 内的所有 fragment mz。
  scanNum <- c(DIAdata@scanindex[scan.max], DIAdata@scanindex[scan.max + 1])
  correctindex <- (scanNum[1] + 1):scanNum[2]

  #'finding co-eluting ions
  # 这些 mz 是 apex scan 中出现的 fragment candidate。
  mz.frag <- DIAdata@env$mz[correctindex]
  if (length(mz.frag) < 1) {
    # `NA_character_` 是 R 中明确的 character 类型缺失值。
    return(NA_character_)
  }
  # legacy 规则：只保留比 precursor mz 至少小 10 Da 的 fragments。
  index <- which(mz.frag < precurmz - 10)
  if (length(index) < 1) {
    return(NA_character_)
  }
  mz.frag <- mz.frag[index]

  #' using correlations to find fragments
  # precursor EIC intensity 被作为 native peak trace。
  native.peak <- peaks$intensity
  fragment.list <- NULL
  for (k in seq_along(mz.frag)) {
    mz0 <- mz.frag[k]
    # 对每个 fragment candidate，用同样的 legacy m/z tolerance 提取 fragment EIC。
    mzmin <- max(minmz, mz0 - mz0 * mz_tol)
    mzmax <- min(maxmz, mz0 + mz0 * mz_tol)

    #'fragment peaks
    # fragment EIC 使用同一个 RT window，目的是检查它是否与 precursor 共洗脱。
    frag.peak <- xcms::rawEIC(
      DIAdata,
      mzrange = cbind(mzmin, mzmax),
      rtrange = cbind(rtmin, rtmax)
    )
    frag.peak <- frag.peak$intensity
    # 如果任一 trace 没有变化，correlation 没有意义，跳过。
    if (sd(native.peak) == 0 || sd(frag.peak) == 0) {
      next
    }
    # legacy intensity threshold：fragment EIC 最大强度低于 2000 时跳过。
    if (max(frag.peak) < 2000) {
      next
    }

    #'correlation calculations
    # `cor()` 计算 precursor EIC 和 fragment EIC 的相关性。
    corr <- cor(native.peak, frag.peak)

    #'save the fragments to list
    # legacy 阈值：只有 `corr > 0.9` 的 fragment 被保留。
    if (corr > 0.9) {
      # 保存格式为 `fragment_mz,max_intensity`。
      list.frag <- paste(c(mz.frag[k], ",", max(frag.peak)), collapse = "")
      if (length(fragment.list) < 1) {
        fragment.list <- list.frag
      } else {
        # 多个 fragments 用 `;` 串联。
        fragment.list <- paste(c(fragment.list, ";", list.frag), collapse = "")
      }
    }
  }
  if (length(fragment.list) < 1) {
    return(NA_character_)
  }
  return(fragment.list)
}

#' -------------------------------------------
#' function to copy MS2 data to MS1 matrix, but without precursor information
#'--------------------------------------------
#' @param xmsn
#' @param precursor
#'
#' @return
#'
# ms2copy() 把指定 `precursor` window 的 MS2 scans 拷贝成一个 pseudo `xcmsRaw` object。
# 这样后续可以继续调用 `rawEIC()`，在 MS2 fragment mz 维度上提取类似 chromatogram 的 trace。
# CAVEAT: 这是对 `xcmsRaw` 内部 S4 slot / env 结构的低层操作，依赖 xcms 对象结构。
# 当前代码实际写入 `@env$mz`, `@env$intensity`, `@mzrange`, `@scanindex`,
# `@scantime`, `@tic`, `@acquisitionNum`, `@polarity` 等字段；
# 没有显式填充 `@env$profile` 或 `@env$scanindex` 这类可能存在的其他内部字段。
ms2copy <- function(xmsn, precursor) {
  # `new("xcmsRaw")` 创建一个空的 xcmsRaw S4 object。
  x <- new("xcmsRaw")
  # `@env` 是 xcmsRaw 中保存大向量数据的 environment。
  x@env <- new.env(parent = .GlobalEnv)
  # 找出所有 precursor mz 等于目标 DIA window 的 MS2 scans。
  index <- which(xmsn@msnPrecursorMz == precursor)
  # 下面这些 slot 让 pseudo object 保留 scan 编号、RT 和 polarity 等基本信息。
  x@tic <- xmsn@msnAcquisitionNum[index]
  x@scantime <- xmsn@msnRt[index]
  x@acquisitionNum <- xmsn@msnAcquisitionNum[index]
  x@polarity <- xmsn@polarity[seq_along(index)]
  len2 <- length(xmsn@msnPrecursorMz)
  index_total <- 0
  index3 <- 0
  for (j in seq_along(index)) {
    # `@msnScanindex` 是 MS2 mz/intensity 大向量里的 scan 起点。
    # 如果当前 scan 是最后一个，就取到 `@env$msnMz` 末尾。
    if (index[j] == len2) {
      index2 <- (xmsn@msnScanindex[index[j]] + 1):length(xmsn@env$msnMz)
    } else {
      index2 <- (xmsn@msnScanindex[index[j]] + 1):xmsn@msnScanindex[index[j] + 1]
    }
    # `index3` 记录新 pseudo object 中每个 scan 在合并 mz/intensity 向量里的起点。
    index3 <- c(index3, index3[length(index3)] + length(index2))
    # `index_total` 收集所有目标 scans 对应的原始 mz/intensity 索引。
    index_total <- c(index_total, index2)
  }
  # 去掉初始化时放进去的占位元素。
  index_total <- index_total[-1]
  index3 <- index3[-length(index3)]
  # 把目标 precursor window 的 MS2 mz/intensity 复制到新 object 的 env 中。
  x@env$mz <- xmsn@env$msnMz[index_total]
  x@env$intensity <- xmsn@env$msnIntensity[index_total]
  # `rawEIC()` 需要 mz range 和 scan index 这类基本结构。
  x@mzrange <- c(min(x@env$mz), max(x@env$mz))
  x@scanindex <- as.integer(index3)
  return(x)
}
