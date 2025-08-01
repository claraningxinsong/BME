
#' Title
#'
#' @param seed
#' @param nF
#' @param prop_S
#' @param duration
#' @param hazard_S
#' @param hazard_Sc
#' @param dropout_S
#' @param dropout_Sc
#' @param w
#' @param ratio
#' @param orr_S
#' @param orr_Sc
#' @param rho_S
#' @param rho_Sc
#' @param nsim
#' @param enroll_month
#'
#' @return
#' @export
#'
#' @examples
getOC_data <- function(seed = 2025, nsim = 1000, nF = 600, prop_S = 0.5,
                       duration = 20, hazard_S = c(0.1, 0.1), hazard_Sc = c(0.12, 0.12),
                       dropout_S = c(0, 0), dropout_Sc = c(0, 0),
                       w = 1, ratio = 1,
                       orr_S = c(0.2, 0.3), orr_Sc = c(0.2, 0.4),
                       rho_S = 0.7, rho_Sc = 0.7,
                       enroll_month = 19) {

  set.seed(seed)
  enrolled_by_month <- numeric(nsim)

  for(i in 1:nsim){

  # Simulate initial data
  dat_initial <- simu_enrich_trial(
    n = nF,
    prop_S = prop_S,
    ratio = ratio,
    duration = duration,
    hazard_S = hazard_S,
    hazard_Sc = hazard_Sc,
    dropout_S = dropout_S,
    dropout_Sc = dropout_Sc,
    w = w,
    orr_S = orr_S,
    orr_Sc = orr_Sc,
    rho_S = rho_S,
    rho_Sc = rho_Sc
  ) %>%
    dplyr::arrange(enterTime)

  # Count how many enrolled by month
  enrolled_by_month[i] <- sum(dat_initial$enterTime <= enroll_month)
  }

  return(enrolled_by_month)
}
