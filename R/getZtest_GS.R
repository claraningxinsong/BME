#' Get test statistics for Group Sequential method
#'
#' @param Zstats A list of Z statistics returned from \link{getZstats_GS}
#' @param design Group sequential design parameters from \link{getDesignGroupSequential}
#'
#' @return A list with Z-statistics and rejection decision
#' @importFrom rpact
#'
#' @export
#'
#' @examples
#' d <- simu_enrich_trial(n = 200, prop_S = 0.5, duration = 10)
#' design <- rpact::getDesignGroupSequential(kMax = 2, alpha = 0.0125, informationRates = c(0.5, 1))
#' Zstats <- getZstats_GS(d, targetEvents = c(100,200))
#' getZtest_GS(Zstats, design)
getZtest_GS <- function(Zstats, design) {
  z <- Zstats$z
  cut <- Zstats$cutTime
  C <- design$criticalValues
  k <- design$kMax
  FA_time <- Zstats$FA_time

  reject <- FALSE
  timing <- 0

  for (i in seq_len(k)) {
    if (!is.na(z[i]) && z[i] > C[i]) {
      reject <- TRUE
      timing <- FA_time[i]
      break
    }
  }
  if (!reject) {
    timing <- FA_time[k]
  }
  return(list(
    z = z,
    reject = reject,
    timing = timing
  ))
}
