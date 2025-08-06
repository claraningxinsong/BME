#' Get Z-statistics for a given trial at interim analysis.
#'
#' @param cutTime
#' @param dat A time-to-event dataset returned from \link{simu_enrich_trial}.
#'
#' @return It returns a list of test statistics used for later adjustments.
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' d <- simu_enrich_trial(n = 600, prop_S = 0.5, duration = 20)
#' getZstats_IA(d, targetEvents.Sc = 105)
getZstats_IA_bydate <- function(dat, cutTime) {
  ## Step 1: Subset for subgroup Sc
  dSc <- dat %>% filter(.data$subgroup == 0)

  ## Initialize outputs
  z.Sc <- NA
  p.Sc <- NA
  obsEvents.Sc <- NA
  hr.Sc.IA <- NA

  ## Step 2: Cut data by event count
  d <- cut_by_date(dSc, cut_time = cutTime)

  # Check for failure
  if (nrow(d) == 0 || !all(c("calendarCutoff", "survTimeCut", "eventCut", "trt") %in% names(d))) {
    warning("cut_by_event() returned invalid or empty result in getZstats_IA()")
    return(list(z.Sc = z.Sc, p.Sc = p.Sc, hr.Sc.IA = hr.Sc.IA, obsEvents.Sc = obsEvents.Sc))
  }

  ## Step 3: Perform one-sided log-rank test
  IA_time <- d$calendarCutoff[1]
  res <- logrank.one.sided(time = d$survTimeCut, event = d$eventCut, group = d$trt, STRATA = NULL)

  ## Step 4: Extract stats
  z.Sc <- res$z
  p.Sc <- res$p
  obsEvents.Sc <- sum(res$obs)
  hr.Sc.IA <- exp(-(res$z) / sqrt(1 / sum(1 / res$obs)))

  ## Step 5: Return
  return(list(
    z.Sc = z.Sc,
    p.Sc = p.Sc,
    hr.Sc.IA = hr.Sc.IA,
    obsEvents.Sc = obsEvents.Sc,
    timing.IA = IA_time
  ))
}
