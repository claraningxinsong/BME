#' Get operating characteristics via simulations for an biomarker enrichment design with group sequential method
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
#' @param HR.Sc.threshold Hazard ratio threshold for futility.
#' @param targetEvents.S.noexpand The target number of events in S without expansion
#' @param targetEvents.S.expand The target number of events in S with expansion
#' @param targetEvents.F  The target number of events in F
#' @param hazard_S hazard rates (h_control, h_treatment) for S
#' @param hazard_Sc hazard rates (h_control, h_treatment) for Sc
#' @param dropout_S dropout hazard rates (h_control, h_treatment) for S
#' @param dropout_Sc dropout hazard rates (h_control, h_treatment) for Sc
#' @param w weight parameter in cumulative enrollment pattern.
#' @param ratio randomization ratio r:1, where r refers to treatment arm; for equal randomization, r=1.
#' @param orr_S Objective response (binary: 1 = response, 0 = non-response) for S
#' @param orr_Sc Objective response (binary: 1 = response, 0 = non-response) for Sc
#' @param rho_S Correlation between ORR and time to event endpoint for S
#' @param rho_Sc Correlation between ORR and time to event endpoint for Sc
#' @param alpha1 Significance level for subgroup S
#' @param alpha2 Significance level for full population F
#' @param kMax_S maximum number of stages K for subgroup S. Must be a positive integer of length 1
#' @param kMax_F maximum number of stages K for full population F. Must be a positive integer of length 1
#' @param informationRates_F  The information rates t_1, ..., t_kMax (that must be fixed prior to the trial), default is (1:kMax) / kMax for F
#' @param informationRates_S  The information rates t_1, ..., t_kMax (that must be fixed prior to the trial), default is (1:kMax) / kMax for S
#' @param selection "survival" or "orr:
#'
#' @return It returns a list.
#' @import dplyr rpact
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' getOC_BME(seed = 2025, nsim = 10, nF = 600, nS_additional = 100,  prop_S = 0.5,duration = 20)
getOC_BME_GS <- function(seed = 2025, nsim = 1000, nF = 600, nS_additional = 100,  prop_S = 0.5,
                        duration = 20, duration_additional = 3,
                        targetEvents.Sc = 105, HR.Sc.threshold = 0.9,
                        targetEvents.S.noexpand = c(150, 210),
                        targetEvents.S.expand = c(200, 280),
                        targetEvents.F = c(300, 420),
                        hazard_S = c(0.1, 0.1), hazard_Sc = c(0.12, 0.12),
                        dropout_S = c(0, 0), dropout_Sc = c(0, 0), w = 1, ratio = 1,
                        orr_S = 0.2, orr_Sc = 0.2, rho_S = 0.7, rho_Sc = 0.7,
                        alpha1 = 0.0125, alpha2 = 0.0125,
                        selection = "survival",
                        orr_thres = 0.1, orr_number = 200)
    {

  ## Start simulation
  set.seed(seed)
  results <- data.frame(
    F.reject = numeric(nsim),
    S.reject = numeric(nsim),
    expand = numeric(nsim)
  )


  for(i in 1:nsim){
    # Step 1: Simulate initial dataset
    dat_initial <- simu_enrich_trial(n = nF, prop_S = prop_S, ratio = ratio, duration = duration,
                                     hazard_S = hazard_S, hazard_Sc = hazard_Sc,
                                     dropout_S = dropout_S, dropout_Sc = dropout_Sc, w = w,
                                     orr_S = orr_S, orr_Sc = orr_Sc, rho_S = rho_S, rho_Sc = rho_Sc)
    dat_F <- dat_initial

    # Step 2: Interim analysis
    if (selection == "survival"){

      zstats_IA <- getZstats_IA(dat_initial, targetEvents.Sc = targetEvents.Sc)
      expand <- zstats_IA$hr.Sc.IA > HR.Sc.threshold

    } else {

      dSc <- dat_F %>% filter(.data$subgroup == 0) %>% arrange(enterTime) %>% slice_head(n = orr_number)
      #assume decision based on first 200 patients

      dSc_trt <- dSc %>% filter(trt == 1)
      dSc_con <- dSc %>% filter(trt == 0)
      orr_diff <- mean(dSc_trt$response) - mean(dSc_con$response)
      expand <- orr_diff > orr_thres
    }


    # Step 3: Simulate additional data if expand
    if (expand){
      dat_additional <- simu_enrich_trial(n = nS_additional*2, prop_S = prop_S, ratio = ratio, duration = duration_additional,
                                          hazard_S = hazard_S, dropout_S = dropout_S, w = w,
                                          orr_S = orr_S, orr_Sc = orr_Sc, rho_S = rho_S, rho_Sc = rho_Sc) %>%
        filter(subgroup == 1) %>%
        mutate(enterTime = enterTime + duration, calendarTime = calendarTime + duration)

      dat_S <- bind_rows(dat_initial, dat_additional) %>% arrange(subgroup) %>%  filter(subgroup == 1)
      targetEvents_S <- targetEvents.S.expand

    } else {
      dat_S <-dat_initial %>%  filter(subgroup == 1)
      targetEvents_S <- targetEvents.S.noexpand
    }

    # Step 4: Group sequential loop
    # get critical boundaries
    targetEvents_F <- targetEvents.F

    informationRates_S <- targetEvents_S / max(targetEvents_S)
    informationRates_F <- targetEvents_F / max(targetEvents_F)

    design_S <- getDesignGroupSequential(kMax = length(targetEvents_S), alpha = alpha1, informationRates = informationRates_S)
    design_F <- getDesignGroupSequential(kMax = length(targetEvents_F), alpha = alpha2, informationRates = informationRates_F)

    # get Z statistics
    Zstats_S <- getZstats_GS(dat_S, targetEvents_S)
    Zstats_F <- getZstats_GS(dat_F, targetEvents_F)

    # Step 5: Test: compare to critical boundaries
    reject_S <- getZtest_GS(Zstats_S, design_S)
    reject_F <- getZtest_GS(Zstats_F, design_F)

    # Store results
    results$F.reject[i] <- reject_F$reject
    results$S.reject[i] <- reject_S$reject
    results$expand[i] <- expand
  }

  return(results)
}
