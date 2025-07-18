#' Cut a dataset for analysis at a specified event count
#'
#' @param data A time-to-event dataset returned from \link{simu_enrich_trial}.
#' @param targetEvents Event count at which data cutoff is to be made.
#'
#' @return A data frame ready for survival analysis including columns \code{calendarCutoff},
#' \code{survTimeCut}, \code{eventCut}.
#' @export
#'
#' @examples
#' d <- simu_enrich_trial(n = 100, prop_S = 0.5, ratio = 1)
#' dcut <- cut_by_event(d, 10)
cut_by_event <- function(data, targetEvents) {
  # Basic check: are required columns present?
  required_cols <- c("calendarTime", "event", "survTime", "enterTime")
  if (!all(required_cols %in% names(data))) {
    warning("Missing required columns in data.")
    return(data.frame())  # return empty safely
  }

  # Order by event time
  data0 <- data
  data0.order <- data0[order(data0$calendarTime), ]

  # Get rows where events occurred
  data.event <- data0.order[data0.order$event == 1, ]

  # Check if enough events exist
  if (nrow(data.event) < targetEvents || targetEvents < 1) {
    warning(paste("Insufficient events: requested", targetEvents,
                  "but only", nrow(data.event), "available."))
    return(data.frame())  # return empty
  }

  # Assign event sequence and determine cutoff
  data.event$event.seq <- seq_len(nrow(data.event))
  cutoff_time <- data.event$calendarTime[data.event$event.seq == targetEvents]
  data0$calendarCutoff <- cutoff_time

  # Apply cutoff
  data0$survTimeCut <- ifelse(data0$calendarTime <= cutoff_time,
                              data0$survTime,
                              pmax(0, cutoff_time - data0$enterTime))
  data0$eventCut <- ifelse(data0$calendarTime <= cutoff_time,
                           data0$event, 0)

  # Return patients who entered before cutoff
  return(data0[data0$enterTime < cutoff_time, ])
}
