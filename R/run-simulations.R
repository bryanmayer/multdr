#' Runs multiple simulations for infection time by ID and given infection detection times
#'
#' This function takes data of multiple subjects with repeated exposure data and runs mulitple stochastic simulations estimating infection time
#' infection time that is linked to the given detection time. For more flexibility on returning exact infection times,
#' use TBD.
#'
#'@param nsim The amount of simulations to run
#'@param raw_exposure_data The data.frame time containing series of all possible exposures (must include a day and count name) for all subjects (FamilyID)
#'@param raw_sample_times The sampling days that infection would actually be observed for each subject
#'@param exposure_rate The rate of exposures to be sampled from raw_exposure_data
#'@param hazard_risk The per time hazard risk for the hazard function
#'@param clr The exponential pathogen clearance rate
#'@param lag The minimum amount of time between infection and when it could be detected (lag = 0 means infection is instantly detectable)
#'@param return_KM return the Kaplan-Meier estimates from all the runs (if false, just returns raw data)
#'@param parallel, Use parallel computation for simulations across nsim
#'@return Either a data.frame of infection time simulations for each family matched to sample times or a Kaplan-Meier object
#'@export

run_sample_times_simulations =  function(nsims = 100, raw_exposure_data, raw_sample_times, exposure_rate,
                                         hazard_risk, clr = 2, lag = 0, return_KM = F, parallel = F){

  exposure_data = plyr::ldply(unique(raw_exposure_data$FamilyID), function(fid){
    sample_exposure_data(raw_exposure_data = subset(raw_exposure_data, FamilyID == fid), exposure_rate = exposure_rate)
  })

  simulations = plyr::ldply(1:nsims, function(i){
    plyr::ldply(unique(exposure_data$FamilyID), function(fid){
      out = simulate_infection_sample_times(exposure_ts = subset(exposure_data, FamilyID == fid),
                                            sample_times = subset(raw_sample_times, FamilyID == fid)$days,
                                            hazard_risk = hazard_risk, clr = clr, lag = lag)
      out$FamilyID = fid
      out$run = i
      out
    })
  }, .parallel = parallel)

  if(return_KM) {return(survfit(Surv(time = detect_time, event = infected) ~ 1, data = simulations))
  } else return(simulations)
}


#' Runs multiple simulations for infection time by ID
#'
#' This function takes data of multiple subjects with repeated exposure data and runs mulitple stochastic simulations estimating infection time
#' The simulation function (written in cpp) slowly crawls through all time points, the infection times are mapped to sample times after
#' use TBD.
#'
#'@param nsim The amount of simulations to run
#'@param raw_exposure_data The data.frame time containing series of all possible exposures (must include a day and count name) for all subjects (FamilyID)
#'@param raw_sample_times The sampling days that infection would actually be observed for each subject
#'@param exposure_rate The rate of exposures to be sampled from raw_exposure_data
#'@param hazard_risk The per time hazard risk for the hazard function
#'@param clr The exponential pathogen clearance rate
#'@param lag The minimum amount of time between infection and when it could be detected (lag = 0 means infection is instantly detectable)
#'@param return_KM return the Kaplan-Meier estimates from all the runs (if false, just returns raw data)
#'@param parallel, Use parallel computation for simulations across nsim
#'@return Either a data.frame of infection time simulations for each family matched to sample times or a Kaplan-Meier object
#'@export
run_simulations = function(nsims = 100, raw_exposure_data, raw_sample_times, exposure_rate,
                                         hazard_risk, clr = 2, lag = 0, return_KM = F, parallel = F){

  exposure_data = plyr::ldply(unique(raw_exposure_data$FamilyID), function(fid){
    sample_exposure_data(raw_exposure_data = subset(raw_exposure_data, FamilyID == fid), exposure_rate = exposure_rate)
  })

  simulations = plyr::ldply(unique(exposure_data$FamilyID), function(fid){
    sub_data = subset(exposure_data, FamilyID == fid)
    max_time_lag = sub_data$max_time[1] - lag
    sub_data_lag = dplyr::arrange(subset(sub_data, days < max_time_lag), days)
    exposures_in = sub_data_lag$count
    time_in = sub_data_lag$days

    plyr::ldply(1:nsims, function(i){
      #this returns -1 if there isn't infection
      infection_time = simulate_infection_times_cpp(exposure_set = exposures_in, time_set =  time_in, lag = lag,
                                                    hazard_risk = hazard_risk, clr = clr)
      detect_time =  if(infection_time < 0) {
        max(subset(raw_sample_times, FamilyID == fid)$days)
        }else assign_detection_time(infection_time, subset(raw_sample_times, FamilyID == fid)$days)

      data.frame(
        FamilyID = fid,
        run = i,
        infected = 1 * (infection_time > 0),
        detect_time = detect_time,
        infection_time = if(infection_time > 0) infection_time else sub_data$max_time[1]
      )
    }, .parallel = parallel)
  })

  if(return_KM) {return(survival::survfit(survival::Surv(time = detect_time, event = infected) ~ 1, data = simulations))
  } else return(simulations)
}

