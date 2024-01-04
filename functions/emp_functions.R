# functions for threshold calculations on empirical data sets

# linearity test
#' @param dt the empirical data frame
#' @param xvar the name of driver variable
#' @param yvar the name of the response variable
#' @param knots knots for the gam() function
#' @param smooth_type the type of smoother for the gam() function

lin_test <- function(dt, xvar, yvar, knots = 4, smooth_type = "tp"){

  dtIn <- dt

  # first rename the xvar and yvar columns driver and response
  x_ind <- which(colnames(dtIn)==xvar)
  y_ind <- which(colnames(dtIn)==yvar)

  colnames(dtIn)[x_ind] <- "driver"
  colnames(dtIn)[y_ind] <- "obs_response"

  # fit the gam
  dt.lm <- lm(obs_response ~ driver, data = dtIn)
  dt.gm <- gam(obs_response~s(driver,k=knots,bs=smooth_type),data = dtIn)

  # get the AICs
  AICs <- AIC(dt.lm, dt.gm)

  return(AICs)
}


# function for fitting
#' @param dt the empirical data frame
#' @param xvar the name of driver variable
#' @param yvar the name of the response variable
#' @param xlength for generating predictions; make set of xvals with length equal to max(100, xlength*range(driver))
all_fits <- function(dt, xvar, yvar, xlength){

  dtIn <- dt

  # first rename the xvar and yvar columns driver and response
  x_ind <- which(colnames(dtIn)==xvar)
  y_ind <- which(colnames(dtIn)==yvar)

  colnames(dtIn)[x_ind] <- "driver"
  colnames(dtIn)[y_ind] <- "obs_response"

  # add "sim" column because jack_results expects this
  dtIn$sim <- rep(1, length(dtIn$driver))

  # filter out the NAs
  dtIn <- dtIn[which(is.na(dtIn$driver*dtIn$obs_response)==F),]

  # get the gam and derivative predictions for full data and each jackknifing iteration
  all_results <- jack_results(simdt = as.data.frame(dtIn), sim_choice = c(1), x_pred = NA, xlength, span = 0.1)
  #ind_results <- all_results$ind_dfs # results for each individual jackknife iteration
  #full_results <- all_results$full_dfs # results for detecting thresholds in full data set
  #summ_results <- all_results$summ_dfs # summary of jackknifed results

  return(all_results)

}


# function for threshold calculations
#' @param dt the empirical data frame
#' @param xvar the name of driver variable
#' @param yvar the name of the response variable
#' @param xlength for generating predictions; make set of xvals with length equal to max(100, xlength*range(driver))
#' @param all_results output of the all_fits function above
thresh_calc <- function(dt, xvar, yvar,xlength, all_results){

  ind_df <- all_results$ind_dfs # gam and deriv fits for each jackknife iteration
  full_df <- all_results$full_df # gam and deriv fits for full dataset
  summ_df <- all_results$summ_df # median and 95% quantiles of jackknifed gam and deriv fits

  dtIn <- dt

  # first rename the xvar and yvar columns driver and response
  x_ind <- which(colnames(dtIn)==xvar)
  y_ind <- which(colnames(dtIn)==yvar)

  colnames(dtIn)[x_ind] <- "driver"
  colnames(dtIn)[y_ind] <- "obs_response"

  # add "sim" and "thresh_loc" columns because jack_thresh expects this
  dtIn$sim <- rep(1, length(dtIn$driver))
  dtIn$thresh_loc <- rep(NA, length(dtIn$driver))
  dtIn$thresh_loc_z <- rep(NA, length(dtIn$driver))

  # filter out the NAs
  dtIn <- dtIn[which(is.na(dtIn$driver*dtIn$obs_response)==F),]

  # calculate the summary statistics across jackknifed thresholds
  jack_summ <- jack_thresh(as.data.frame(dtIn), x_pred = NA, xlength)

  # calculate values of thresholds for each jackknife iteration
  all_thresh <- ind_thresh(ind_df, summ_df)

  # calculate values of thresholds in full dataset
  tot_thresh <- full_thresh(full_df)

  return(list(jack_summ = jack_summ, all_thresh = all_thresh, tot_thresh = tot_thresh))

}



# function to. back-tranform standardized data to raw scale
#' @param xz value on standardized scale
#' @param x_raw vector with the raw values that were used for standardization
invzfun <- function(xz, x_raw){
  xz*sd(x_raw, na.rm = T) + mean(x_raw, na.rm = T)
}



