#' Get operating characteristics via simulations for an enrichment design
#'
#' @param seed random seed for reproducibility
#' @param nsim number of replicates
#' @param nF total number of subjects for full population.
#' @param nS_additional number of subjects for subgroup (BM+)
#' @param prop_S proportion for sub-population group S.
#' @param duration enrollment duration in months for full population
#' @param duration_additional if expand, enrollment duration in months for additional patients in BM+
#' @param targetEvents.Sc The target number of events in Sc; length of 1,
#' futility or population selection will be given at IA without efficacy testing.
#' @param targetEvents  A numeric vector of length 2, specifying the number of events for:
#'        FA_S (final analysis in S), and  FA_F (final analysis in F).
#' @param HR.Sc.threshold Hazard ratio threshold for futility.
#' @param hazard_S hazard rates (h_control, h_treatment) for S
#' @param hazard_Sc hazard rates (h_control, h_treatment) for Sc
#' @param dropout_S dropout hazard rates (h_control, h_treatment) for S
#' @param dropout_Sc dropout hazard rates (h_control, h_treatment) for Sc
#' @param w weight parameter in cumulative enrollment pattern.
#' @param alpha1 Significance level for subgroup S
#' @param alpha2 Significance level for full population F
#' @param ratio randomization ratio r:1, where r refers to treatment arm; for equal randomization, r=1.
#'
#' @return It returns a list.
#' @export
#'
#' @examples
#' res <- getOC(seed = 24232, nsim=10)
#' lapply(res, function(x) mean(apply(x[,1:2], 1, any, na.rm=TRUE)))
getOC_BME <- function(seed = 2025, nsim = 1000, nF = 600, nS_additional = 100,  prop_S = 0.5,
                  duration = 25, duration_additional = 3,
                  targetEvents.Sc = 105,  targetEvents = c(210, 280, 420),
                  HR.Sc.threshold = 0.9,
                  hazard_S = c(0.1, 0.1), hazard_Sc = c(0.12, 0.12),
                  dropout_S = c(0, 0), dropout_Sc = c(0, 0), w = 1, ratio = 1,
                  alpha1 = 0.0125, alpha2 = 0.0125){

  ## Start simulation
  set.seed(seed)

  results <- vector("list", nsim)


  for(i in 1:nsim){
    # Step 1: Simulate initial dataset
    dat_initial <- simu_enrich_trial(n = nF, prop_S = prop_S, ratio = ratio, duration = duration,
                                     hazard_S = hazard_S, hazard_Sc = hazard_Sc,
                                     dropout_S = dropout_S, dropout_Sc = dropout_Sc, w = w)

    # Step 2: Interim analysis
    zstats_IA <- getZstats_IA(dat_initial, targetEvents.Sc = targetEvents.Sc)
    expand <- zstats_IA$hr.Sc.IA > HR.Sc.threshold

    # Step 3: Simulate additional data if expand
    if (expand){
      dat_additional <- simu_enrich_trial(n = nS_additional*2, prop_S = prop_S, ratio = ratio, duration = duration_additional,
                                          hazard_S = hazard_S, dropout_S = dropout_S, w = w) %>%
        filter(subgroup == 1) %>%
        mutate(enterTime = enterTime + duration, calendarTime = calendarTime + duration)
    } else {
      dat_additional <- NULL
    }

    # Step 4: Final analysis
    zstats_FA <- getZstats_FA(dat_initial, dat_additional, targetEvents = targetEvents, expand = expand)

    # Step 5: Hypothesis testing
    reject <- getZtests_BME(zstats_FA, alpha1 = alpha1, alpha2 = alpha2)

    # Store results
    results[[i]] <- list(
      reject = reject,
      expand = expand
    )
  }

  return(results)
}
