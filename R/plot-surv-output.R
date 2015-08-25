#'Plots KM from simulations to compare to the data
#'
#'@param data_KM The KM object of the data
#'@param sim_data The simulated data: either a KM object or a data.frame from  \code{\link{run-simulations.R}} and \code{\link{simulate-infections.R}}
#'@param rate_str Legend text when a simulated KM is passed in.  Otherwise this is populated from the simulation results
#'@param leg.x,leg.y Legend position coordinates
#'
#'@return a KM plot with multiple outputs


plot_surv_output = function(data_KM, sim_data, rate_str = NULL, leg.x = 250, leg.y = 1){
  survival::plot.survfit(data_KM)
  if(any(is(sim_data) == "survfit")){
    survival::lines.survfit(sim_data, conf.int = F, col = "red", mark.time = F, lty = 1)
    legend_text = c("Data", if(is.null(rate_str)) "Simulation" else rate_str)
    print(legend_text)
    legend(250, 1, legend_text, col = c("black", "red"), lty = 1)
  } else{ # a set of survival data
    rates = unique(sim_data$rate)
    plot_col = rev(brewer.pal(8, "Accent"))
    rate_str = NULL
    for(i in 1:length(rates)){
      sub_sim_data = subset(sim_data, rate == rates[i])
      temp_surv = survival::survfit(Surv(time = detect_time, event = infected) ~ 1, data = sub_sim_data)
      survival::lines.survfit(temp_surv, conf.int = F, col = plot_col[i], mark.time = F)
      rate_str = c(rate_str, paste("Rate = ", rates[i], "/day", sep = ""))
    }
    legend(leg.x, leg.y, c("Data", rate_str), col = c(1, plot_col[1:length(rates)]), lty = 1)
  }
}
