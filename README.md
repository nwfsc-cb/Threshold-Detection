# Threshold-Detection
Repo for simulation studies to explore detectability of ecological thresholds

Details on repo files:

run_simulations.Rmd = code for running all of the theoretical simulation studies (everything except the empirical hake case study)

sim_plots.Rmd = code for plotting the results of the theoretical simulation studies

conceptual_plots.Rmd = code for mkaing the conceptual figures (threshold definitions figure, simulation steps figures)

cciea_data.Rmd = code for downloading the CCIEA indicator data and analyzing it to inform choice of time series lenghts and observation errors in the generic simulations

hake_example.Rmd = code for the empirical case study on Pacific hake (analyzing relationships between ocean temperature and the proportion of hake in Canada, simulating datasets similar to the observed data)

simulation ouputs folder = outputs from the simulations run in run_simulations.Rmd (simulations using focal and nonfocal driver-response relationships, simulations with a linear and exponential covariate that is vs isn't included in the estimation models) and hake_example.Rmd

hake data folder = empirical hake and temperature data used in "hake_example.Rmd"

functions folder: all of the functions used in the simulation analyses
- sim_data.R = function for simulating data with a specified underlying driver-response relationship
- sim_emp_data.R = same as sim_data.R but with more options for the driver-response relationship (to better match potential empirical relationships)
- lin_check.R = function that fits nonlinear (GAM) and linear models to simulated driver-response data and uses model selection (AIC) to choose which is a better fit
- lin_check_cov.R = same as lin_check.R function but for cases with simulated driver-covariate-response data (candidate models include those with a covariate in addition 
  to the driver)
- jack_thresh.R = functions for estimating threshold locations in simulated data sets with leave-one-out jackknife resampling
- jack_thresh_cov.R = same as jack_thresh.R but fits GAMs that include a driver and covariate rather than just a driver
- jack_results.R = functions that calculates and returns results of the jackknifing for a simulated dataset(s) (jack_thresh.R just returns the threshold summary statistics, 
  jack_results.R returns the GAM and and derivative predictions for the full data set and each jackknife iteration)
- jack_results_cov.R = same as jack_results.R but fits GAMs that include both a driver and a covariate
- thresh_values.R = function that returns the values of the thresholds detected in each jackknife iteration in a simulated data set(s)
- emp_functions.R = functions for cleaning up empirical datasets and applying the threshold calculation functions to them

load_functions.R = loads all of the functions in the "functions" folder into the workspace

figurepds folder = pdfs of all the figures in the manuscript
