#' function to extract fragment from DIA window
#'
#' @param Library
#' @param path
#' @param msfiles
#' @param mz_tol
#' @export
GetFrag <- function(mz_tol, DIAisowin, RTwin) {
  #' the path to save results
  path <- getwd()

  #' the path to save raw data
  path.data <- paste0(path, "/data")
  path.peak <- paste0(path, "/peaklist")
  path.results <- paste0(path, "/results")

  setwd(path.peak) # the list of peak documents
  peakfiles <- list.files()
  if (length(peakfiles) < 1) {
    warning(
      "No peaklist files found in ", path.peak, "; returning NULL.",
      call. = FALSE
    )
    setwd(path)
    return(NULL)
  }

  setwd(path.data) # the list of raw MS data
  msfiles <- list.files()

  fragments <- NULL
  for (k in seq_along(peakfiles)) {

    # get the peak list
    setwd(path.peak)
    print(c("getfragment...", k))
    mycpd <- readxl::read_excel(peakfiles[k])
    mycpd$MS2 <- rep(0, nrow(mycpd))

    # get the raw data
    setwd(path.data)
    temp <- strsplit(peakfiles[k], ".xlsx")[[1]]
    temp <- temp[1] # sample name
    index <- grep(temp, msfiles)

    if (length(index) == 0) {
      warning(
        "No matching raw MS data found for peaklist ",
        peakfiles[k],
        "; skipping.",
        call. = FALSE
      )
      next
    } # no MS data
    if (length(index) > 1) {
      warning(
        "Multiple raw MS data files matched peaklist ",
        peakfiles[k],
        "; using first match: ",
        msfiles[index[1]],
        call. = FALSE
      )
      index <- index[1]
    }

    xraw <- xcms::xcmsRaw(msfiles[index], includeMSn = TRUE)

    #'extracting precursor DIA windows
    precursor <- preclist(xraw)
    for (j in seq_len(nrow(mycpd))) {
      mz <- mycpd$mz[j]
      DIAwin <- which(abs(mz - precursor) <= (0.5 * DIAisowin))
      if (length(DIAwin) < 1) {
        next
      }
      DIAwin <- precursor[DIAwin[1]]

      #' get fragments from each window based correlation >0.9
      fragments <- getfrag(xraw, j, mycpd, mz_tol, DIAwin, RTwin)
      if (length(fragments) == 0) {
        next
      }
      mycpd$MS2[j] <- fragments
    }
    setwd(path.results)
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
preclist <- function(xmsn) {
  #'extracting DIA windows
  precmz <- xmsn@msnPrecursorMz
  if (length(precmz) < 1) {
    return(precmz)
  }
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
getfrag <- function(xraw, index, mycpd, mz_tol, DIAwin, RTwin) {
  mycpd$rt <- mycpd$rt * 60
  precurmz <- mycpd$mz[index]

  #' read rawdata for each DIA window
  DIAdata <- ms2copy(xraw, DIAwin)
  mzrange <- DIAdata@mzrange
  minmz <- mzrange[1]
  maxmz <- mzrange[2]
  mzmin <- max(minmz, precurmz - precurmz * mz_tol)
  mzmax <- min(maxmz, precurmz + precurmz * mz_tol)
  rtrange <- DIAdata@scantime
  rtmin <- max(min(rtrange), mycpd$rt[index] - 10)
  rtmax <- min(max(rtrange), mycpd$rt[index] + 10)

  #'extracting precursor ions and peaks
  peaks <- xcms::rawEIC(
    DIAdata,
    mzrange = cbind(mzmin, mzmax),
    rtrange = cbind(rtmin, rtmax)
  )

  #'finding the scan number of the peak top
  scan.max <- which.max(peaks$intensity)
  scan.max <- peaks$scan[scan.max[1]]
  if (scan.max >= length(DIAdata@scanindex)) { # the last scanning point
    return(NULL)
  }
  scanNum <- c(DIAdata@scanindex[scan.max], DIAdata@scanindex[scan.max + 1])
  correctindex <- (scanNum[1] + 1):scanNum[2]

  #'finding co-eluting ions
  mz.frag <- DIAdata@env$mz[correctindex]
  if (length(mz.frag) < 1) {
    return(NA_character_)
  }
  index <- which(mz.frag < precurmz - 10)
  if (length(index) < 1) {
    return(NA_character_)
  }
  mz.frag <- mz.frag[index]

  #' using correlations to find fragments
  native.peak <- peaks$intensity
  fragment.list <- NULL
  for (k in seq_along(mz.frag)) {
    mz0 <- mz.frag[k]
    mzmin <- max(minmz, mz0 - mz0 * mz_tol)
    mzmax <- min(maxmz, mz0 + mz0 * mz_tol)

    #'fragment peaks
    frag.peak <- xcms::rawEIC(
      DIAdata,
      mzrange = cbind(mzmin, mzmax),
      rtrange = cbind(rtmin, rtmax)
    )
    frag.peak <- frag.peak$intensity
    if (sd(native.peak) == 0 || sd(frag.peak) == 0) {
      next
    }
    if (max(frag.peak) < 2000) {
      next
    }

    #'correlation calculations
    corr <- cor(native.peak, frag.peak)

    #'save the fragments to list
    if (corr > 0.9) {
      list.frag <- paste(c(mz.frag[k], ",", max(frag.peak)), collapse = "")
      if (length(fragment.list) < 1) {
        fragment.list <- list.frag
      } else {
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
ms2copy <- function(xmsn, precursor) {
  x <- new("xcmsRaw")
  x@env <- new.env(parent = .GlobalEnv)
  index <- which(xmsn@msnPrecursorMz == precursor)
  x@tic <- xmsn@msnAcquisitionNum[index]
  x@scantime <- xmsn@msnRt[index]
  x@acquisitionNum <- xmsn@msnAcquisitionNum[index]
  x@polarity <- xmsn@polarity[seq_along(index)]
  len2 <- length(xmsn@msnPrecursorMz)
  index_total <- 0
  index3 <- 0
  for (j in seq_along(index)) {
    if (index[j] == len2) {
      index2 <- (xmsn@msnScanindex[index[j]] + 1):length(xmsn@env$msnMz)
    } else {
      index2 <- (xmsn@msnScanindex[index[j]] + 1):xmsn@msnScanindex[index[j] + 1]
    }
    index3 <- c(index3, index3[length(index3)] + length(index2))
    index_total <- c(index_total, index2)
  }
  index_total <- index_total[-1]
  index3 <- index3[-length(index3)]
  x@env$mz <- xmsn@env$msnMz[index_total]
  x@env$intensity <- xmsn@env$msnIntensity[index_total]
  x@mzrange <- c(min(x@env$mz), max(x@env$mz))
  x@scanindex <- as.integer(index3)
  return(x)
}
