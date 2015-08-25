#' Takes a series of exposures and times and samples from them based on an exposure rate
#'
#' Effectively maps a continuum of exposure to a discrete set of actual exposures that could cause infection.
#'
#' @param raw_exposure_data The subsetted data.frame containing time (days) and exposure (count) for one subject
#' @param exposure_rate The rate of sampling from the raw_exposure_data
#' @return A data.frame containing only the sampled exposure times
#' @export

#this does it by rate
sample_exposure_data = function(raw_exposure_data, exposure_rate){
  if(exposure_rate > 5 & exposure_rate != 10) stop("Exposure rate must be an integer between 1-5 or 10")
  days = raw_exposure_data$days
  exposure_times = round(seq(min(days), max(days), by = 1/exposure_rate), 1)
  subset(raw_exposure_data, round(days, 1) %in% exposure_times)
}
