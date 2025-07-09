#' Get Z-statistics and test results for Group Sequential method
#'
#' @param dat_S Dataset for subgroup S (BM+), with columns: surtime, event, trt
#' @param design_S Group sequential design object from rpact::getDesignGroupSequential()
#' @param n_target Total number of events planned for final analysis
#'
#' @return A list with Z-statistics, observed events, and rejection decision
#'
#' @export
#' #' @examples
#' d <- simu_enrich_trial(n = 200, prop_S = 0.5, duration = 10)
#' design <- getDesignGroupSequential(kMax = 2, alpha = 0.0125, informationRates = c(0.5, 1))
#' getZstats_GS(d,design)
getZtest_GS <- function(dat, design) {
  infoRates <- design$informationRates
  criticalValues <- design$criticalValues
  k <- length(infoRates)
  n <- nrow(dat)

  z_stats <- c()
  events_obs <- c()
  reject_each <- c()
  stopped_at <- NA

  for (i in seq_len(k)) {
    targetEvents <- ceiling(infoRates[i] * n)

    # Cut data based on event count
    d <- cut_by_event(dat, targetEvents = targetEvents)

    # Perform one-sided log-rank test
    res <- logrank.one.sided(
      time = d$survTimeCut,
      event = d$eventCut,
      group = d$trt,
      STRATA = NULL
    )

    z_i <- res$z
    z_stats <- c(z_stats, z_i)
    events_obs <- c(events_obs, sum(d$eventCut))

    # Check for rejection
    reject_i <- !is.na(z_i) && z_i > criticalValues[i]
    reject_each <- c(reject_each, reject_i)

    if (reject_i) {
      stopped_at <- i
      break
    }
  }

  return(list(
    z = z_stats,
    obsEvents = events_obs,
    reject_each = reject_each,
    reject = any(reject_each),
    stopped_at = stopped_at
  ))
}
