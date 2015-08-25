#' Simulates infection time for one exposure time series based on given infection detection times
#'
#' This function takes a time series of exposures, exposure times, and infection detection times and returns a stochastic
#' infection time that is linked to the given detection time. For more flexibility on returning exact infection times,
#' use TBD.
#'
#'@param exposure_ts The data.frame time series of exposure (must include a day and count name)
#'@param sample_times The sampling days that infection would actually be observed
#'@param hazard_risk The per time hazard risk for the hazard function
#'@param clr The exponential pathogen clearance rate
#'@param lag The minimum amount of time between infection and when it could be detected (lag = 0 means infection is instantly detectable)
#'@return An infection time that matches a time given in sample_times
#'@export

simulate_infection_sample_times = function(exposure_ts, sample_times, hazard_risk, clr, lag = 0, debug = F){
  exposure_ts = dplyr::arrange(exposure_ts, days)
  sample_times = sort(sample_times)
  sample_times_lag = c(min(sample_times), tail(sample_times, - 1) - lag)
  if(any(sample_times_lag < 0) | sample_times_lag[1] >= sample_times_lag[2]) stop(paste("Check lag setting. First two lagged samples times:",
                                                                                        sample_times_lag[1:2]))
  final_time_index = length(sample_times_lag)
  infection = 0
  cur_exposure = 0
  t_index = 1
  while(t_index < final_time_index & !infection){
    start_time = sample_times_lag[t_index]
    end_time = sample_times_lag[t_index + 1]
    exposure_set = subset(exposure_ts, days >= start_time & days < end_time)$count
    time_set =  subset(exposure_ts, days >= start_time & days < end_time)$days

    exposure_set[1] = exposure_set[1] + cur_exposure

    risk = calc_risk_mult(exposure_set, time_set, survival_time = end_time, hazard_risk, clr)

    prob_sample = runif(1)
    if(is.na(prob_sample < risk)) stop(paste("Missing value in risk estimate at time:", sample_times_lag[t_index]))
    if(prob_sample < risk){
      infection = 1
      t_index = t_index + 1 #detected at the end of the time bounds
    } else{
      cur_exposure = calc_current_exposure(exposure_set, time_set, final_time = end_time, clr)
      t_index = t_index + 1
    }
  }
  if(sample_times_lag[t_index] + lag != sample_times[t_index]) print(paste("Lag not lining up, (lagged time+lag, detect time)", c(sample_times_lag[t_index] + lag,
                                                                                                                                  sample_times[t_index])))
  if(infection == 1){
    detect_time = sample_times[t_index] #
  } else detect_time = max(sample_times) #if censored then always on the final exposure day

  data.frame(
    detect_time = detect_time, #even if uninfected, infection cannot be observed in the lag period after the simulation finishes
    infected = infection
  )
}



