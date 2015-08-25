#'Calculates the remaining exposure at given time
#'
#'This function takes a set of exposures and exposures times (or a single exposure) and calculates the remaining level at the end of a certain time frame
#'The exposures decay exponentially (see \code{\link{calc_risk}}) and if end_time -> infty, the remaining exposure will be zero
#' @param exposure_series The values of repeated pathogen exposures
#' @param exposure_times The times of exposures
#' @param final_time The time to calculate remaining exposure at
#' @param clr The clearance rate of exposure
#' @return remaining exposure and end of time frame (given by final_times)
#' @examples
#' calc_current_exposure_cpp(rep(1e7, 5), 0:4, 5, 2)
#'
#' @export
calc_current_exposure = function(exposure_series, exposure_times, final_time, clr, na.zero = F){
  check_exposure_input(exposure_series, exposure_times)
  if(any(is.na(exposure_series)) & na.zero) exposure_series[which(is.na(exposure_series))] <- 0

  if(length(exposure_series) == 1) { #for a single exposure, method assumed time = 0 unless specified
    if(is.na(exposure_times)) exposure_times <- 0
    return(exposure_series * exp(-clr * (final_time - exposure_times)))
  }
  sum(exposure_series * exp(-clr * (final_time - exposure_times)))
}
