#'Processes data for and sends data to survival likelihood function
#'
#' @param exposure_data A data frame that contains time and exposures
#' @param hazard_risk The per time risk of exposure during persistence
#' @param clr The clearance rate of exposure
#' @param censor_start Start of the interval censored outcome, should be passed in from the data
#' @return Returns a log survival likelihood (see \code{\link{calc_risk_mult}}) for uninfected and an interval censored likelihood
#' (see \code{\link{calc_censor_likelihood}}) for infected
likelihood_setup_call = function(exposure_data, betat, clr, lag, end_time, censor_start_time){

  }
