#' Perform one-sided logrank test
#'
#' This function performs a one-sided logrank test. The standard logrank test in the survival package
#' only provides two-sided tests. This version facilitates one-sided testing where a larger z-statistic
#' indicates better treatment effect (i.e., lower hazard).
#'
#' @param time Survival time (numeric vector)
#' @param event Event status (0 = censor, 1 = event; integer or logical vector)
#' @param group Group indicator (0 = control, 1 = experimental; must be binary)
#' @param STRATA Optional stratification variable
#'
#' @return A list with:
#' \describe{
#'   \item{z}{Z-statistic from the one-sided logrank test}
#'   \item{p}{One-sided p-value}
#'   \item{obs}{Observed number of events per group}
#' }
#'
#' @examples
#' n <- 100
#' time <- c(rexp(n, rate=log(2)/12), rexp(n, rate=log(2)/12*1.2))
#' event <- rep(1, n*2)
#' group <- c(rep(0, n), rep(1, n))
#' STRATA <- rep(c(0,1), n)
#' logrank.one.sided(time, event, group, STRATA)
#'
#' @importFrom survival survdiff Surv strata
#' @export
logrank.one.sided <- function(time, event, group, STRATA = NULL) {
  # Ensure proper types
  time <- as.numeric(time)
  event <- as.integer(event)
  group <- as.integer(group)

  # Check group variability
  if (length(unique(group[!is.na(group)])) < 2) {
    return(list(z = NA, p = NA, obs = c(NA, NA)))
  }

  # Perform logrank test
  if (!is.null(STRATA)) {
    STRATA <- as.factor(STRATA)
    lr.test <- survdiff(Surv(time, event) ~ group + strata(STRATA))
  } else {
    lr.test <- survdiff(Surv(time, event) ~ group)
  }

  # Extract observed and expected
  if (is.matrix(lr.test$obs)) {
    otmp <- rowSums(lr.test$obs)
    etmp <- rowSums(lr.test$exp)
  } else {
    otmp <- lr.test$obs
    etmp <- lr.test$exp
  }

  # Compute z
  better <- as.numeric(otmp[2] < etmp[2])
  sign <- 2 * better - 1
  z <- sqrt(lr.test$chisq) * sign
  p <- 1 - pnorm(z)

  return(list(z = z, p = p, obs = otmp))
}
