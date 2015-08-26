#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;
//'  This is a C++ version of the current remaining exposure function
// [[Rcpp::export]]
double calc_current_exposure_cpp(double exposure, double start_time, double final_day, double clr){
  return exposure * exp(-clr * (final_day - start_time));
}

//'  This is a C++ version of the risk calculation function
// [[Rcpp::export]]
double calc_risk_mult_cpp(double exposure, double start_time, double final_day, double hazard_risk, double clr){
  return 1 - exp(-hazard_risk/clr * exposure * (1 -  exp(-clr * (final_day - start_time))));
}


//'  This is a C++ version of the simulation function, that will replace the while loop in the simulation call
//'  Assigning infection time to sample time is done later
//'  @export
// [[Rcpp::export]]
double simulate_infection_times_cpp(NumericVector exposure_set, NumericVector time_set,
                       double lag, double hazard_risk, double clr) {
  double start_day;
  double end_day;
  double infection_time;
  double risk;
  NumericVector prob_sample(1);

  int t_index = 0;
  int infection = 0;
  double cur_exposure = 0;
  int final_time_index = (exposure_set.size() - 1);

  while(t_index < final_time_index && infection == 0){

    start_day = time_set[t_index];
    end_day = time_set[t_index + 1];
    risk = calc_risk_mult_cpp(exposure_set[t_index] + cur_exposure, 0, end_day - start_day, hazard_risk, clr);

    prob_sample = runif(1);
    if(prob_sample[0] < risk){
      infection = 1;
      t_index++;
    } else{
      cur_exposure = calc_current_exposure_cpp(exposure_set[t_index] + cur_exposure, 0, end_day - start_day, clr);
      t_index++;
    }
  }
  if(infection == 0){
    infection_time = -1;
  } else{
    infection_time = time_set[t_index] + lag;
  }
  return infection_time;
}



/*** R
calc_current_exposure(1e7, 0, 5, 2)
calc_current_exposure_cpp(1e7, 0, 5, clr = 2)
calc_risk_mult(1e7, 0, 3, hazard_risk = 1e-9, clr = 1)
calc_risk_mult_cpp(1e7, 0, 3, hazard_risk = 1e-9, clr = 1)
simulate_infection_times_cpp(c(1e7, 1e9, 0), 0:3, lag = 1, hazard_risk = 1e-9, clr = 1)
*/
