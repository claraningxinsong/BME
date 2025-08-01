#' Get Z-statistics for a given trial at final analysis
#'
#' @param dat_initial Initial time-to-event dataset returned from \link{simu_enrich_trial}.
#' @param dat_additional Additional time-to-event dataset returned from \link{simu_enrich_trial}.
#' @param targetEvents A numeric vector of length 3, specifying the number of events for: FA_S_noexpand, FA_S_expand, FA_F.
#' @param expand 0=No, 1=Yes
#'
#' @return It returns a list of test statistics used for later adjustments.
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' d1 <- simu_enrich_trial(n = 600, prop_S = 0.5, duration = 20)
#' d2 <- simu_enrich_trial(n = 100, prop_S = 0.5, duration = 20)
#' getZstats_FA(d1, d2, targetEvents = c(210, 280, 420), expand=0)
getZstats_FA <- function(dat_initial, dat_additional, targetEvents, expand){

  ## split dataset by subgroup
  if (expand) {
    dat <- bind_rows(dat_initial, dat_additional) %>% arrange(subgroup)
    targetEvents_new <- c(targetEvents[2], targetEvents[3])
  } else {
    dat <- dat_initial
    targetEvents_new <- c(targetEvents[1], targetEvents[3])
  }

  dS <- dat %>% filter(.data$subgroup==1)
  dF <- dat_initial

  ## initial results
  z.S <- z.F <- NA # z-statistics
  p.S <- p.F <- NA # p-values
  obsEvents.S <- obsEvents.F <- NA # observed number of events


  ## cut data at FA when S reaches certain percent
  # population S
  d <- cut_by_event(dS, targetEvents = targetEvents_new[1])
  FA_time_S <- d$calendarCutoff[1]
  res <- logrank.one.sided(time = d$survTimeCut, event = d$eventCut,
                           group = d$trt, STRATA = NULL)
  z.S <- res$z # non-adjusted Z statistic
  p.S <- res$p
  obsEvents.S <- sum(d$eventCut)
  nS <- nrow(d)


  ## cut data at FA when F reaches certain percent
  # population F
  d <- cut_by_event(dF, targetEvents = targetEvents_new[2])
  FA_time_F <- d$calendarCutoff[1]
  res <- logrank.one.sided(time = d$survTimeCut, event = d$eventCut,
                           group = d$trt, STRATA = NULL)
  z.F <- res$z # non-adjusted Z statistic
  p.F <- res$p
  obsEvents.F <- sum(res$obs)
  nF <- nrow(d)

  # resulted test statistic
  return(list(z.S = z.S, z.F = z.F, p.S = p.S, p.F = p.F,
              obsEvents.S = obsEvents.S, obsEvents.F = obsEvents.F,
              nS.FA = nS, nF.FA = nF,
              timing.S = FA_time_S, timing.F = FA_time_F))
}
