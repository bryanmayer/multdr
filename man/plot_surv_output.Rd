% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/plot-surv-output.R
\name{plot_surv_output}
\alias{plot_surv_output}
\title{Plots KM from simulations to compare to the data}
\usage{
plot_surv_output(data_KM, sim_data, rate_str = NULL, leg.x = 250,
  leg.y = 1)
}
\arguments{
\item{data_KM}{The KM object of the data}

\item{sim_data}{The simulated data: either a KM object or a data.frame from  \code{\link{run-simulations.R}} and \code{\link{simulate-infections.R}}}

\item{rate_str}{Legend text when a simulated KM is passed in.  Otherwise this is populated from the simulation results}

\item{leg.x,leg.y}{Legend position coordinates}
}
\value{
a KM plot with multiple outputs
}
\description{
Plots KM from simulations to compare to the data
}

