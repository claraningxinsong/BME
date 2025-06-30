#' Get adaptive enrichment design test statistics for a given trial.
#'
#' @param zstatsAn object returned by getZstats_FA()
#' @param alpha1 Significance level for subgroup S
#' @param alpha2 Significance level for full population F
#'
#' @return It returns a list of test statistics used for later adjustments.
#' @importFrom stats pnorm qnorm
#' @export A list indicating whether the null hypothesis is rejected in S and/or F
#'
#' @examples
getZtests_BME <- function(zstats, alpha1 = 0.0125, alpha2 = 0.0125){
  z.alpha1 <- qnorm(1 - alpha1)
  z.alpha2 <- qnorm(1 - alpha2)

  S_reject <- zstats$z.S > z.alpha1
  F_reject <- zstats$z.F > z.alpha2

  return(list(p_S = zstats$p.S, S_reject = S_reject,
              p_F = zstats$p.F, F_reject = F_reject))

}
