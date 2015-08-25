#'Calculates the likelihood of the exposure data
#'
#' @param exposure_series The vector of repeated pathogen exposure values, should be pre-processed
#' @param exposure_times The times of exposures
#' @param end_time Final time: either censor time (infected = 0) or infection time (infected = 1).  If there is a lag, process it that first.
#' @param hazard_risk The per time risk of exposure during persistence
#' @param clr The clearance rate of exposure
#' @param censor_start Start of the interval censored data (end_time is the end)
#' @return Returns a log survival likelihood (see \code{\link{calc_risk_mult}}) for uninfected and an interval censored likelihood
#' (see \code{\link{calc_censor_likelihood}}) for infected
#' @export
log_likelihood_fun = function(infected, exposure_series, exposure_times, end_time, hazard_risk, clr = 2, censor_start = NA){
  if(length(exposure_series) != length(exposure_times)) stop("Uneven match between exposures and exposure times")
  if(!infected) lik = calc_risk_mult(exposure_series, exposure_times, survival_time = end_time, hazard_risk, clr, log.surv = T)
  if(infected) lik = calc_censor_likelihood(exposure_series, exposure_times, censor_start, censor_end = end_time, hazard_risk, clr)
  -lik #return positive values
}


#'Calculates the interval censored likelihood when there is censoring related to the infection time
#'
#'For interval censoring with left censor time Tl and right censor time Tr, this returns S(Tl) - S(Tr).
#'
#' @param exposure_series The vector of repeated pathogen exposure values, should be pre-processed
#' @param exposure_times The times of exposures
#' @param censor_start Start of the censor period (exposure up to this time are survived)
#' @param censor_end End of the censor period (infection occured in exposures between censor_start and censor_end)
#' @param hazard_risk The per time risk of exposure during persistence
#' @param clr The clearance rate of exposure
#' @param loglik Whether to return log (default TRUE)
#' @return Returns a survival likelihood for interval censored outcome (see \code{\link{calc_risk_mult}} for calculations)
#' @export
calc_censor_likelihood = function(exposure_series, exposure_times, censor_start, censor_end, hazard_risk, clr, loglik = T){
  right_surv = exp(calc_risk_mult(exposure_series, exposure_times, survival_time = censor_end, hazard_risk, clr, log.surv = T))
  if(length(which(exposure_times < censor_start)) > 0){
    left_surv = exp(calc_risk_mult(exposure_series[which(exposure_times < censor_start)], exposure_times[which(exposure_times < censor_start)],
                                  survival_time = censor_start, hazard_risk, clr, log.surv = T))
  }else left_surv = 1 #if exposure time starts at day 0

  if(loglik) {return(log(left_surv - right_surv + 1e-12)) #add a correction when they are the same
    }else return(left_surv - right_surv)
}
