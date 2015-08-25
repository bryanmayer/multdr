#' Calculate risk from one pathogen exposure
#'
#' This function takes a pathogen exposure and calculates the risk of infection over its persistence within host
#'
#' The risk is derived using a surival function.  exposure = p, survival_time = T, hazard_risk = \eqn{\beta}, and clearance rate = c.
#' Note: every calculation assumes that time starts at 0
#'
#' The pathogen is assumed to be cleared exponentially (\eqn{P(0) = p})
#' \deqn{\frac{\mathrm{d}P}{\mathrm{d}t} = -cP}{dP/dt = -cP}
#'
#' The pathogens pose a risk of infection during their persistance that can be used to calculate the survival function
#' \deqn{\lambda(t) = \beta * P(t)}
#' \deqn{S(T) =  \mathrm{e}^{-\int_0^T{\lambda(t)\mathrm{d}t}}}{S(T) =  exp(-integral(0, T) \lambda(t)dt)}
#' \deqn{S(T) = \mathrm{e}^{-\frac{\beta}{c} * p * (1 - \mathrm{e}^{-c(T)})}}{S(t) = exp(-\beta/c * p * (1 - exp(-cT))}
#' which simplifies as \eqn{T \rightarrow \infty}{T -> Inf} to
#' \deqn{S(T) = \mathrm{e}^{-\frac{\beta}{c}p}}{S(T) = exp(-\beta/c * p)}
#'
#'The risk is \eqn{1 - S(T)}.
#'
#' @param exposure The initial level of the pathogen exposure
#' @param survival_time The range of time to calculate the risk over (default infinity)
#' @param hazard_risk The per time risk of exposure during persistence
#' @param clr The clearance rate of exposure
#' @return risk calculation over the time specified in survival time
#'
#' @seealso \code{\link{calc_risk_mult} for mulitple exposures.}
#'
#' @examples
#' calc_risk(1e8, 1, 1e-8, 1)
#' calc_risk(exposure = 1e8, hazard_risk = 1e-8, clr = 1)
#' @export
calc_risk = function(exposure, survival_time = Inf, hazard_risk, clr){
  1 - exp(- hazard_risk/clr * exposure * (1 - exp(-clr * survival_time)))
}


#' Calculate risk from mulitple pathogen exposures
#'
#' This function takes mulitple pathogen exposure and calculates the risk of infection over their persistence within host
#'
#' The risk for a single exposure is derived in \code{\link{calc_risk}}.  The accumulation of exposure assuming exponential clearance is independent.
#' \deqn{S(T) = \mathrm{e}^{-\frac{\beta}{c}\sum\limits_{i=1}^n{p_i(1 - \mathrm{e}^{-c(T - t_i)})}}}{S(t) = exp(-\beta/c *sum(i = 1..n) p_i *(1 - exp(-c(T - t_i))}
#' which simplifies as \eqn{T \rightarrow \infty}{T -> Inf} to
#' \deqn{S(T) = \mathrm{e}^{-\frac{\beta}{c}\sum\limits_{i=1}^n{p_i}}}{S(T) = exp(-\beta/c *sum(i = 1..n) p_i)}
#'
#' @param exposure_series The values of repeated pathogen exposures
#' @param exposure_times The times of exposures
#' @param survival_time the final time to calculate survival to (default infinity)
#' @param hazard_risk The per time risk of exposure during persistence
#' @param clr The clearance rate of exposure
#' @return risk calculation over the time specified in survival time
#'
#' @seealso \code{\link{calc_risk} for single exposures.}
#'
#' @examples
#' ##Single exposure examples
#' calc_risk_mult(1e8, 0, 1, 1e-8, 1) # should return same as calc_risk(1e8, 1, 1e-8, 1)
#' calc_risk_mult(1e8, 1, 2, 1e-8, 1) # should return same as calc_risk(1e8, 1, 1e-8, 1)
#' ## multiple exposures
#' exposures = c(1e7, 1e9, 1e5); exp_times = c(0, 1, 2)
#' calc_risk_mult(exposures, exp_times, 3, hazard_risk = 1e-9, clr = 1)
#' calc_risk_mult(exposures, exp_times, hazard_risk = 1e-9, clr = 1) #total risk
#' 1-exp(-1e-9/1 * sum(exposures))
#' @export
#'
calc_risk_mult = function(exposure_series, exposure_times, survival_time = Inf, hazard_risk, clr, na.zero = F){
  check_exposure_input(exposure_series, exposure_times)
  if(any(is.na(exposure_series)) & na.zero) exposure_series[which(is.na(exposure_series))] <- 0
  1 - exp(- hazard_risk/clr * sum(exposure_series * (1 - exp(-clr * (survival_time - exposure_times)))))
}


