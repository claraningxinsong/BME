#' Get test statistics for a given trial.
#'
#' @param dat A time-to-event dataset returned from \link{simu_enrich_trial}.
#' @param targetEvents.Sc The target number of events in Sc; length of 1,
#' futility or population selection will be given at IA without efficacy testing.
#'
#' @return It returns a list of test statistics used for later adjustments.
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' d <- simu_enrich_trial(n = 600, prop_S = 0.5, duration = 20)
#' getZstats_IA(d, targetEvents.Sc = 105)
getZstats_IA <- function(dat, targetEvents.Sc){

  ## split dataset by subgroup
  dSc <- dat %>% filter(.data$subgroup==0)

  ## initial results
  z.Sc <- NA # z-statistics
  p.Sc <- NA # p-values
  obsEvents.Sc <- NA # observed number of events

  ## cut data at IA
  # population Sc
  d <- cut_by_event(dSc, targetEvents = targetEvents.Sc)
  IA_time <- d$calendarCutoff[1]
  res <- logrank.one.sided(time = d$survTimeCut, event = d$eventCut,
                           group = d$trt, STRATA = NULL)
  z.Sc <- res$z # non-adjusted Z statistic
  p.Sc <- res$p
  obsEvents.Sc <- sum(res$obs)
  hr.Sc.IA <- exp(-(res$z)/sqrt(1/sum(1/res$obs)))

  # resulted test statistic
  return(list(z.Sc = z.Sc, p.Sc = p.Sc,
              hr.Sc.IA = hr.Sc.IA,
              obsEvents.Sc = obsEvents.Sc))
}
