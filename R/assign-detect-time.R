#' Takes a simulated infection time and maps it to an actual detection time from the data
#'
#'@param sim_infection_time Infection time from a simulation
#'@param raw_data_times Sampling times from the data where detection would be detect
#'@return A time that matches the detection/observation times in the data
#'@examples
#' assign_detection_time(7.2, c(0,13,21,33))
#' assign_detection_time(7.2, c(0,13,21,33))
#'@export
assign_detection_time = function(sim_infection_time, raw_data_times){
  detect_times = sort(raw_data_times)
  detect_times[min(which(sim_infection_time <= detect_times))]
}
