#' Looks at exposure series and exposure times for common errors
#'
#' @keywords internal
check_exposure_input = function(exposure_series, exposure_times, na.zero = F){
  if(any(is.na(exposure_series))){
    if(!na.zero) stop("Missing data in the exposure series; set na.zero = T to convert to 0")
  }
  if(any(is.na(exposure_times))) stop("Missing data in the exposure time series")

  return(NULL)

}
