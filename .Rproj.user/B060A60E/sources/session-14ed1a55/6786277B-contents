#' Get test statistics for Group Sequential method
#'
#' @param dat_S Dataset for subgroup S (BM+), with columns: surtime, event, trt
#' @param design_S Group sequential design object from rpact::getDesignGroupSequential()
#' @param n_target Total number of events planned for final analysis
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
  criticalValues <- design$criticalValues
  k <- design$kMax

  reject <- FALSE

  for (i in seq_len(k)) {
    if (!is.na(z[i]) && z[i] > criticalValues[i]) {
      reject <- TRUE
      break
    }
  }

  return(list(
    z = z,
    reject = reject
  ))
}
