#' function for evaluating thresholds using jackknifing
#' required packages:
#' mgcv (for fitting the gam)
#' gratia (for derivative calculations),
#' pracma (for the findpeaks function to find max/mins in the derivatives),
#' data.table (for between function to determine whether CIs contain 0)

#' update from v5: add eps argument, change xvals argument, and standardize the data

#' function inputs:
#' @param simdt data frame with the simulated data (output of simfun2)
#' @param x_pred optional vector of driver values to use for generating predictions from the gams and
#' calculating their derivatives
#' @param xlength if x_pred = NA, make set of xvals for each dataset with length equal to max(100, xlength*range(driver))
#' @param z_scale if TRUE, standardize the data before fitting (default), if FALSE use raw data
#' @param thresh_methods vector of methods to use for estimating the threshold location.
#' Can include `abs_max_d2` (location where 2nd deriv of gam is furthest from zero), `min_d2`
#' (location where 2nd deriv is at its minimum- i.e., most negative- value), `zero_d2`
#' (location where 2nd deriv of gam crosses zero), or `zero_d1` (location where first deriv of gam
#' crosses zero/is zero). Default is to do all of the above
#' @param sig_criteria the type of criteria to use in evaluating whether a threshold is significant
#' Can be one or more of `none` (no significance criteria), `jack_quant` (use the quantiles from
#' across all jackknife iterations as upper/lower bounds), or `sim_int` (use the simultaneous
#' intervals calculated by gratia::derivatives for each individual jackknife iteration)
#' default is to use all 3
#' @param alpha quantile for calculating the CIs of the threshold estimate (defaults to 0.05 for 95% CIs)
#' @param knots number of knots to use in the gam fitting, defaults to 4
#' @param smooth_type type of smooth used with gams, defaults to "tp"
#' @param span smoothing step/span for thresholds, values closer to 1 = more smoothing (defaults to
#' 0.1, used by Holsman et al. 2020)
#' @param eps_val eps value (the finite difference) for the gratia::derivatives function, defaults to 5*10^-6
#' @param pk_height minimum absolute height of a peak in the derivative (relative to median value of the deriv)
#' required for that peak to be returned by the findpeaks function
#' @param pk_ups minimum number of increasing steps before a peak is reached (and decreasing steps after the peak)
#' required for a peak to be returned by the findpeaks function


jack_thresh <- function(simdt, x_pred, xlength = 50, z_scale = TRUE, thresh_methods = c("abs_max_d2", "min_d2", "zero_d2", "zero_d1"),
                        sig_criteria = c("none", "jack_quant", "sim_int"), alpha = 0.05, knots = 4,
                        smooth_type = "tp", eps_val = 5*10^-6, span = 0.1, pk_height = 0.05, pk_ups = 3){

  # get the number of simulations that were run
  nsim <- length(unique(simdt$sim))

  # get the true value of the threshold (for calculating rmse)
  if(z_scale==TRUE){
  thresh_loc <- simdt$thresh_loc_z[1]
  } else{
    thresh_loc <- simdt$thresh_loc[1]
  }

  # turn the driver values into a dataframe
  if(is.na(x_pred[1])==F){

    xvals <- x_pred

    xdt <- data.frame(
      driver = x_pred
    )
  }

  # findpeaks parameters to put into the functions for finding maxes and mins
  pk_height1 <- pk_height
  pk_ups1 <- pk_ups

  # length of sig_criteria vector
  lsc <- length(sig_criteria)

  z_crit <- qnorm(1 - alpha / 2) # for calculating CIs of threshold estimates based on assumption of normal distribution

  # make the holding matrices to store mean, SEs, and RMSEs of the thresholds, number of jackknife
  # iterations that detected a threshold for each simulation replicate, and average number of
  # thresholds detected per jackknife iteration for each threshold calculation method
  thresh_mean <- matrix(NA, nrow = nsim*lsc, ncol = length(thresh_methods))# mean threshold estimate
  thresh_se <- matrix(NA, nrow = nsim*lsc, ncol = length(thresh_methods))# SE of threshold estimates
  ci_n_up <- matrix(NA, nrow = nsim*lsc, ncol = length(thresh_methods)) # upper CI of threshold estimate calculated using normal dist. approx
  ci_n_low <- matrix(NA, nrow = nsim*lsc, ncol = length(thresh_methods)) # lower CI of threshold estimate calculated using normal dist. approx
  ci_em_up <- matrix(NA, nrow = nsim*lsc, ncol = length(thresh_methods)) # upper CI of threshold estimate calculated using empirical dist.
  ci_em_low <- matrix(NA, nrow = nsim*lsc, ncol = length(thresh_methods)) # lower CI of threshold estimate calculated using empirical dist.
  thresh_n <- matrix(NA, nrow = nsim*lsc, ncol = length(thresh_methods))# number of jackknife iterations that detected at least one threshold (goes into calculating the se, also criteria for whether threshold was detected for the simulation)
  full_thresh <- matrix(NA, nrow = nsim*lsc, ncol = length(thresh_methods)) # whether or not a threshold was detected in the full data set


  # I think the easiest thing to do is to first do the methods check on the full data set for all the versions
  # then do the jackknifing and store all the results, and then go through the results for each sig criteria type

  for(i in 1:nsim){ # for each simulation

    # subset the data for the ith simulation and standardize
    datIN <- simdt[which(simdt$sim == i), ]

    if(z_scale==TRUE){
    datIN$driver <- zfun(datIN$driver)
    datIN$obs_response <- zfun(datIN$obs_response)
    }


    # get the xvals
    if(is.na(x_pred[1])==T){

      rangei <- max(datIN$driver, na.rm = T) - min(datIN$driver, na.rm = T)
      xvals <- seq(from = min(datIN$driver, na.rm = T), to = max(datIN$driver, na.rm = T), length.out =max(100, rangei*xlength))# make sure length is at least 100

      xdt <- data.frame(
        driver = xvals
      )

    }

    # store the means and sd's of the driver
    #x_means[i] <- mean(datIN$driver)
    #x_sds[i] <- sd(datIN$driver)

    # first fit a gam to the whole data set and see whether a threshold was detected
    gam_full <- gam(obs_response~s(driver,k=knots,bs=smooth_type),data = datIN) # fit a gam to the jackd data set

    #if(z_scale==TRUE){
    #gam_full <- gam(obs_response~ -1 + s(driver,k=knots,bs=smooth_type),data = datIN) # fit a gam to the jackd data set
    #} else{
      #gam_full <- gam(obs_response~s(driver,k=knots,bs=smooth_type),data = datIN) # fit a gam to the jackd data set
   # }



    # calculate the smoothed derivatives of the gam (smooth them because the second deriv can be really noisy)
    # first deriv
    D1_full <- gratia::derivatives(gam_full, data = xdt, order = 1, eps = eps_val)
    D1_full_up <- predict(loess(d1 ~ driver, data=data.frame(d1 = D1_full$upper, driver = xvals), span=span))
    D1_full_mn <- predict(loess(d1 ~ driver, data=data.frame(d1 = D1_full$derivative, driver = xvals), span=span))
    D1_full_low <- predict(loess(d1 ~ driver, data=data.frame(d1 = D1_full$lower, driver = xvals), span=span))

    # second deriv
    D2_full <- gratia::derivatives(gam_full, data = xdt, order = 2, eps = eps_val)
    D2_full_up <- predict(loess(d2 ~ driver, data=data.frame(d2 = D2_full$upper, driver = xvals), span=span))
    D2_full_mn <- predict(loess(d2 ~ driver, data=data.frame(d2 = D2_full$derivative, driver = xvals), span=span))
    D2_full_low <- predict(loess(d2 ~ driver, data=data.frame(d2 = D2_full$lower, driver = xvals), span=span))


    # check with methods (if any) found a threshold on the full data set
    methods_check <- matrix(NA, nrow = length(sig_criteria), ncol = length(thresh_methods)) # max number of sig criteria options is 3
    #methods_check <- rep(NA, length(thresh_methods))

    for(mm in 1:length(thresh_methods)){ # for each threshold method

      thresh_method <- thresh_methods[mm]

      if(thresh_method=="abs_max_d2"){
        # find peaks (max and mins) in the second deriv

        for(ss in 1:length(sig_criteria)){

          sig_criteria_s <- sig_criteria[ss]

          if(sig_criteria_s=="none"){
            d2pks_full1 <- abs_max_threshF(D2_full_mn, sig_type = "none", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low, pkht = pk_height1, pkup = pk_ups1)
            methods_check[ss, mm] <- ifelse(is.na(d2pks_full1)==F, 1, 0)
          }

          if(sig_criteria_s == "jack_quant"){
            methods_check[ss, mm] <- NA # can't evaluate w/o jackknifing
          }

          if(sig_criteria_s=="sim_int"){
            d2pks_full3 <- abs_max_threshF(D2_full_mn, sig_type = "sim_int", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low, pkht = pk_height1, pkup = pk_ups1)
            methods_check[ss, mm] <- ifelse(is.na(d2pks_full3)==F, 1, 0)
          }

          add_index <- ss-1
          full_index <- lsc*i-(lsc-1) + add_index # get the correct index for the ith simulation and ssth sig criteria method

          full_thresh[full_index, mm] <- methods_check[ss, mm] # whether or not a threshold was detected in the full data set


        }

        #abs_max_threshF <- function(deriv_i, sig_type, xvals, deriv_up, deriv_low, pkht, pkup){
       # d2pks_full1 <- abs_max_threshF(D2_full_mn, sig_type = "none", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low, pkht = pk_height1, pkup = pk_ups1)
       # d2pks_full3 <- abs_max_threshF(D2_full_mn, sig_type = "sim_int", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low, pkht = pk_height1, pkup = pk_ups1)

        # 1 if threshold was detected in full data for sig_criteria = "none", also put 1 for sig_criteria = "jack_int", since
        #this step is just to decide whether to do the full jackknifing, which needs to be done to evaluate significance using the "jack_int" method
        #methods_check[mm, 1] <- ifelse(length(d2pks_full1)>0, 1, 0)
        #methods_check[mm, 2] <- methods_check[mm, 1]
        #methods_check[mm, 3] <- ifelse(length(d2pks_full3)>0, 1, 0) # 1 if threshold was detected in full data for sig_criteria = "sim_int"

     }

      if(thresh_method=="min_d2"){
        # find the minima in the second deriv

        for(ss in 1:length(sig_criteria)){

          sig_criteria_s <- sig_criteria[ss]

          if(sig_criteria_s=="none"){
            d2mins_full1 <- min_threshF(D2_full_mn, sig_type = "none", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low, pkht = pk_height1, pkup = pk_ups1)
            methods_check[ss, mm] <- ifelse(is.na(d2mins_full1)==F, 1, 0)
          }

          if(sig_criteria_s == "jack_quant"){
            methods_check[ss, mm] <- NA # can't evaluate w/o jackknifing
          }

          if(sig_criteria_s=="sim_int"){
            d2mins_full3 <- min_threshF(D2_full_mn, sig_type = "sim_int", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low, pkht = pk_height1, pkup = pk_ups1)
            methods_check[ss, mm] <- ifelse(is.na(d2mins_full3)==F, 1, 0)
          }


          add_index <- ss-1
          full_index <- lsc*i-(lsc-1) + add_index # get the correct index for the ith simulation and ssth sig criteria method

          full_thresh[full_index, mm] <- methods_check[ss, mm] # whether or not a threshold was detected in the full data set

        }

        #d2mins_full1 <- min_threshF(D2_full_mn, sig_type = "none", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low, pkht = pk_height1, pkup = pk_ups1)
        #d2mins_full3 <- min_threshF(D2_full_mn, sig_type = "sim_int", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low, pkht = pk_height1, pkup = pk_ups1)
        #methods_check[mm, 1] <- ifelse(length(d2mins_full1)>0, 1, 0)
        #methods_check[mm, 2] <- methods_check[mm, 1]
        #methods_check[mm, 3] <- ifelse(length(d2mins_full3)>0, 1, 0) # 1 if threshold was detected in full data for sig_criteria = "sim_int"

      }

      if(thresh_method=="zero_d2"){

        # find the points where the second deriv crosses zero (=changes sign)

        for(ss in 1:length(sig_criteria)){

          sig_criteria_s <- sig_criteria[ss]

          if(sig_criteria_s=="none"){
            d2roots_full1 <- root_threshF(D2_full_mn, sig_type = "none", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low)
            methods_check[ss, mm] <- ifelse(is.na(d2roots_full1)==F, 1, 0)
          }

          if(sig_criteria_s == "jack_quant"){
            methods_check[ss, mm] <- NA # can't evaluate w/o jackknifing
          }

          if(sig_criteria_s=="sim_int"){
            d2roots_full3 <- root_threshF(D2_full_mn, sig_type = "sim_int", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low)
            methods_check[ss, mm] <- ifelse(is.na(d2roots_full3)==F, 1, 0)
          }

          add_index <- ss-1
          full_index <- lsc*i-(lsc-1) + add_index # get the correct index for the ith simulation and ssth sig criteria method

          full_thresh[full_index, mm] <- methods_check[ss, mm] # whether or not a threshold was detected in the full data set


        }

        #root_threshF <- function(deriv_i, sig_type, xvals, deriv_up, deriv_low){
        #d2roots_full1 <- root_threshF(D2_full_mn, sig_type = "none", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low)
        #d2roots_full3 <- root_threshF(D2_full_mn, sig_type = "sim_int", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low)
        #methods_check[mm, 1] <- ifelse(length(d2roots_full1)>0, 1, 0)
        #methods_check[mm, 2] <- methods_check[mm, 1]
        #methods_check[mm, 3] <- ifelse(length(d2roots_full3)>0, 1, 0) # 1 if threshold was detected in full data for sig_criteria = "sim_int"

      }

      if(thresh_method=="zero_d1"){

        for(ss in 1:length(sig_criteria)){

          sig_criteria_s <- sig_criteria[ss]

          if(sig_criteria_s=="none"){
            d1roots_full1 <- root_threshF(D1_full_mn, sig_type = "none", xvals, deriv_up = D1_full_up, deriv_low = D1_full_low)
            methods_check[ss, mm] <- ifelse(is.na(d1roots_full1)==F, 1, 0)
          }

          if(sig_criteria_s == "jack_quant"){
            methods_check[ss, mm] <- NA # can't evaluate w/o jackknifing
          }

          if(sig_criteria_s=="sim_int"){
            d1roots_full3 <- root_threshF(D1_full_mn, sig_type = "sim_int", xvals, deriv_up = D1_full_up, deriv_low = D1_full_low)
            methods_check[ss, mm] <- ifelse(is.na(d1roots_full3)==F, 1, 0)
          }

          add_index <- ss-1
          full_index <- lsc*i-(lsc-1) + add_index # get the correct index for the ith simulation and ssth sig criteria method

          full_thresh[full_index, mm] <- methods_check[ss, mm] # whether or not a threshold was detected in the full data set


        }

        # find the points where the first deriv crosses zero (=changes sign)
        #d1roots_full1 <- root_threshF(D1_full_mn, sig_type = "none", xvals, deriv_up = D1_full_up, deriv_low = D1_full_low)
        #d1roots_full3 <- root_threshF(D1_full_mn, sig_type = "sim_int", xvals, deriv_up = D1_full_up, deriv_low = D1_full_low)
        #methods_check[mm, 1] <- ifelse(length(d1roots_full1)>0, 1, 0)
        #methods_check[mm, 2] <- methods_check[mm, 1]
        #methods_check[mm, 3] <- ifelse(length(d1roots_full3)>0, 1, 0) # 1 if threshold was detected in full data for sig_criteria = "sim_int"

      }

    } # end of threshold methods check

    # uncomment below to skip the jackknifing if none of the methods detected a threshold in the full data set
   # if(sum(c(methods_check))==0){ # if none of the methods detected a threshold in the full data set for the ith simulation
     # for(ss in 1:length(sig_criteria))
     # add_index <- ss-1
     # full_index <- lsc*i-(lsc-1) + add_index # get the correct index for the ith simulation and ssth sig criteria method
     # thresh_mean[full_index, ] <- rep(NA, length(thresh_methods))# mean threshold estimate
     # thresh_se[full_index, ] <- rep(NA, length(thresh_methods))# SE of threshold estimates
     # thresh_n[full_index, ] <- rep(NA, length(thresh_methods))# number of jackknife iterations that detected at least one threshold

   # } else{

      # make data set for the jackknife resampling
    dd <- datIN
    #dd$driver <- round(dd$driver, 3)
    dd <- dd[ ,c("driver", "obs_response")]

    # make holding matrices for the derivatives
    D1 <- matrix(NA,length(dd[,1]),length(xvals)) # first derivative
    D1_up <- matrix(NA,length(dd[,1]),length(xvals)) # first derivative, upper simultaneous interval
    D1_low <- matrix(NA,length(dd[,1]),length(xvals)) # first derivative, lower simultaneous interval
    D2 <- matrix(NA,length(dd[,1]),length(xvals)) # second derivative
    D2_up <- matrix(NA,length(dd[,1]),length(xvals)) # second derivative, upper simultaneous interval
    D2_low <- matrix(NA,length(dd[,1]),length(xvals)) # second derivative, lower simultaneous interval

    # perform jackknife
    for(int in 1:length(dd[,1])){ # for each observation
      # remove this observation from the dataset
      jackd <- dd[-int, ]

      gami <- gam(obs_response~s(driver,k=knots,bs=smooth_type),data = jackd) # fit a gam to the jackd data set

      #if(z_scale==TRUE){
      #gami <- gam(obs_response~ -1+ s(driver,k=knots,bs=smooth_type),data = jackd) # fit a gam to the jackd data set
      #} else{
       # gami <- gam(obs_response~s(driver,k=knots,bs=smooth_type),data = jackd) # fit a gam to the jackd data set
      #}

      D1i <- gratia::derivatives(gami, data = xdt, order = 1, eps = eps_val)
      D1_up[int, ] <- predict(loess(d1 ~ driver, data=data.frame(d1 = D1i$upper, driver = xvals), span=span))
      D1[int, ] <- predict(loess(d1 ~ driver, data=data.frame(d1 = D1i$derivative, driver = xvals), span=span))
      D1_low[int, ] <- predict(loess(d1 ~ driver, data=data.frame(d1 = D1i$lower, driver = xvals), span=span))

      # second deriv
      D2i <- gratia::derivatives(gami, data = xdt, order = 2, eps = eps_val)
      D2_up[int, ] <- predict(loess(d2 ~ driver, data=data.frame(d2 = D2i$upper, driver = xvals), span=span))
      D2[int, ] <- predict(loess(d2 ~ driver, data=data.frame(d2 = D2i$derivative, driver = xvals), span=span))
      D2_low[int, ] <- predict(loess(d2 ~ driver, data=data.frame(d2 = D2i$lower, driver = xvals), span=span))

      } # end of jackknife iterations


      # get the 2.5% and 97.5% quantiles of the jackknifed first and second derivatives
      d2_low <- apply(D2,2,quantile,probs=0.025) # 2 in apply() means apply quantile() across columns of D2 matrix
      d2_up <- apply(D2,2,quantile,probs=0.975)

      d1_low <- apply(D1,2,quantile,probs=0.025)
      d1_up <- apply(D1,2,quantile,probs=0.975)

      # then apply the significant criteria to all the thresholds that were detected

      # make holding lists for results for the different threshold methods
      # threshold values for each jackknife iteration, methods, and sig criteria
      thresh_l <- vector(mode = "list", length = length(thresh_methods))
      # put place holder matrixes in the list
      for(ll in 1:length(thresh_l)){
        thresh_l[[ll]] <- matrix(NA, nrow = length(dd[,1]), ncol = length(sig_criteria))
      }

      for(int in 1:length(dd[,1])){ # for each of the jackknife iterations
        # get the derivatives for this iteration
        D2int <- D2[int, ] # second deriv for int^th jackknife iteration
        D2_up_int <- D2_up[int, ] # upper simultaneous interval of 2nd deriv
        D2_low_int <- D2_low[int, ] # lower simultaneous interval of 2nd deriv

        D1int <- D1[int, ] # first deriv for int^th jackknife iteration
        D1_up_int <- D1_up[int, ] # upper simultaneous interval of 1st deriv
        D1_low_int <- D1_low[int, ] # lower simultaneous interval of 1st deriv


      for(mm in 1:length(thresh_methods)){
      #for(mm in thresh_found){ # for each threshold method that found a threshold on the full dataset
        thresh_method <- thresh_methods[mm]

      if(thresh_method =="abs_max_d2"){

        #for(int in 1:length(dd[,1])){ # for each of the jackknife iterations

          for(ss in 1:length(sig_criteria)){

            sig_criteria_s <- sig_criteria[ss]

            if(sig_criteria_s=="none"){
              thresh_int <- abs_max_threshF(D2int, sig_type = "none", xvals, deriv_up = d2_up, deriv_low = d2_low, pkht = pk_height1, pkup = pk_ups1)
              thresh_l[[mm]][int, ss] <- thresh_int #ifelse(length(thresh_int)>0, thresh_int, NA)
            }

            if(sig_criteria_s == "jack_quant"){
              thresh_int <- abs_max_threshF(D2int,sig_type = "jack_quant", xvals, deriv_up = d2_up, deriv_low = d2_low, pkht = pk_height1, pkup = pk_ups1)
              thresh_l[[mm]][int, ss] <- thresh_int #ifelse(length(thresh_int)>0, thresh_int, NA)

            }

            if(sig_criteria_s=="sim_int"){
              thresh_int <- abs_max_threshF(D2int, sig_type = "sim_int", xvals, deriv_up = D2_up_int, deriv_low = D2_low_int, pkht = pk_height1, pkup = pk_ups1)
              thresh_l[[mm]][int, ss] <- thresh_int #ifelse(length(thresh_int)>0, thresh_int, NA)
            }

          }
        #}

      }

      if(thresh_method =="min_d2"){

        for(ss in 1:length(sig_criteria)){

          sig_criteria_s <- sig_criteria[ss]

          if(sig_criteria_s=="none"){
            thresh_int <- min_threshF(D2int, sig_type = "none", xvals, deriv_up = d2_up, deriv_low = d2_low, pkht = pk_height1, pkup = pk_ups1)
            thresh_l[[mm]][int, ss] <- thresh_int #ifelse(length(thresh_int)>0, thresh_int, NA)
          }

          if(sig_criteria_s == "jack_quant"){
            thresh_int <- min_threshF(D2int, sig_type = "jack_quant", xvals, deriv_up = d2_up, deriv_low = d2_low, pkht = pk_height1, pkup = pk_ups1)
            thresh_l[[mm]][int, ss] <- thresh_int #ifelse(length(thresh_int)>0, thresh_int, NA)

          }

          if(sig_criteria_s=="sim_int"){
            thresh_int <- min_threshF(D2int, sig_type = "sim_int", xvals, deriv_up = D2_up_int, deriv_low = D2_low_int, pkht = pk_height1, pkup = pk_ups1)
            thresh_l[[mm]][int, ss] <- thresh_int #ifelse(length(thresh_int)>0, thresh_int, NA)
          }

        }
      }

      if(thresh_method =="zero_d2"){

        for(ss in 1:length(sig_criteria)){

          sig_criteria_s <- sig_criteria[ss]

          if(sig_criteria_s=="none"){
            thresh_int <- root_threshF(D2int, sig_type = "none", xvals, deriv_up = d2_up, deriv_low = d2_low)
            thresh_l[[mm]][int, ss] <- thresh_int# ifelse(length(thresh_int)>0, thresh_int, NA)
          }

          if(sig_criteria_s == "jack_quant"){
            thresh_int <- root_threshF(D2int, sig_type = "jack_quant", xvals, deriv_up = d2_up, deriv_low = d2_low)
            thresh_l[[mm]][int, ss] <- thresh_int #ifelse(length(thresh_int)>0, thresh_int, NA)

          }

          if(sig_criteria_s=="sim_int"){
            thresh_int <- root_threshF(D2int, sig_type = "sim_int", xvals, deriv_up = D2_up_int, deriv_low = D2_low_int)
            thresh_l[[mm]][int, ss] <- thresh_int #ifelse(length(thresh_int)>0, thresh_int, NA)
          }

        }

      }


      if(thresh_method =="zero_d1"){

        for(ss in 1:length(sig_criteria)){

          sig_criteria_s <- sig_criteria[ss]

          if(sig_criteria_s=="none"){
            thresh_int <- root_threshF(D1int, sig_type = "none", xvals, deriv_up = d1_up, deriv_low = d1_low)
            thresh_l[[mm]][int, ss] <- thresh_int #ifelse(length(thresh_int)>0, thresh_int, NA)
          }

          if(sig_criteria_s == "jack_quant"){
            thresh_int <- root_threshF(D1int, sig_type = "jack_quant", xvals, deriv_up = d1_up, deriv_low = d1_low)
            thresh_l[[mm]][int, ss] <- thresh_int #ifelse(length(thresh_int)>0, thresh_int, NA)

          }

          if(sig_criteria_s=="sim_int"){
            thresh_int <- root_threshF(D1int, sig_type = "sim_int", xvals, deriv_up = D1_up_int, deriv_low = D1_low_int)
            thresh_l[[mm]][int, ss] <- thresh_int #ifelse(length(thresh_int)>0, thresh_int, NA)
          }

        }


      }

     } # end of threshold calculations for each calculation method

  } # end of iterations across jackknife iterations


    # final threshold calculations
    for(mm in 1:length(thresh_methods)){
    #for(mm in thresh_found){ # for each threshold method that found a threshold on the full dataset

      #thresh_method <- thresh_methods[mm]

      #thresh <- unlist(thresh_l[[mm]][[ss]]) # just unlist thresh_l to get the vector of thresholds

     for(ss in 1:length(sig_criteria)){

       thresh <- thresh_l[[mm]][, ss]
       thresh_real <- thresh[which(is.na(thresh)==FALSE)] # thresholds that weren't NA

       add_index <- ss-1
       full_index <- lsc*i-(lsc-1) + add_index # get the correct index for the ith simulation and ssth sig criteria method
      # (lsc = length(sig_criteria))

    # get the mean and SEs of the chosen thresholds for the ith simulation
    #if(mm %in% thresh_found){ # if thresh method mm detected a threshold on the full data set
    thresh_mean[full_index,mm] <- mean(thresh, na.rm = T) # mean of thresholds detected across the jackknife iterations
    thresh_n[full_index,mm] <- length(thresh_real)# number of jackknife iterations that detected at least one threshold
    thresh_se[full_index,mm] <- ifelse(length(thresh_real) > 2, sd(thresh_real)/sqrt(length(thresh_real)-1), NA)
    #thresh_rmse[full_index, mm] <- ifelse(length(thresh_real) > 1, sqrt(sum((thresh_real-thresh_loc)^2)/length(thresh_real)), NA)
    #thresh_detected[full_index, mm] <- ifelse(thresh_n[full_index,mm]/length(dd[,1]) >= sensitivity, "y", "n")

    ci_n_up[full_index, mm] <- ifelse(length(thresh_real) > 2, thresh_mean[full_index,mm] + z_crit * thresh_se[full_index,mm], NA) # upper CI of threshold estimate calculated using normal dist. approx
    ci_n_low[full_index, mm] <- ifelse(length(thresh_real) > 2, thresh_mean[full_index,mm] - z_crit * thresh_se[full_index,mm], NA) # lower CI of threshold estimate calculated using normal dist. approx
    ci_em_up[full_index, mm] <- ifelse(length(thresh_real) > 1, quantile(thresh_real, 1-alpha/2), NA) # upper CI of threshold estimate calculated using empirical dist.
    ci_em_low[full_index, mm] <- ifelse(length(thresh_real) > 1,quantile(thresh_real, alpha/2), NA) # lower CI of threshold estimate calculated using empirical dist.

     } # end of calculations for each significance criteria

    #} else{ # if threshold wasn't detected in full dataset

     # thresh_mean[i,mm] <- NA # mean of thresholds detected across the jackknife iterations
     # thresh_n[i,mm] <- 0# number of jackknife iterations that detected at least one threshold
     # thresh_se[i,mm] <- NA
     # thresh_detected[i, mm] <- "n"

     # ci_n_up[i,mm] <- NA # upper CI of threshold estimate calculated using normal dist. approx
     # ci_n_low[i,mm] <- NA # lower CI of threshold estimate calculated using normal dist. approx
     # ci_em_up[i,mm] <- NA # upper CI of threshold estimate calculated using empirical dist.
     # ci_em_low[i,mm] <- NA # lower CI of threshold estimate calculated using empirical dist.

   # }


  } # end of calculations for each threshold method

 # } # end of if else for whether methods_check found a threshold

  } # end of iterations for each simulation i

  thresh_methods_ss <- rep(thresh_methods, each = length(sig_criteria))
  nsim_ss <- rep(c(1:nsim), each = length(sig_criteria))

  df <- data.frame(
    #sim = rep(c(1:nsim), length(thresh_methods)*length(sig_criteria)),
    sim = rep(nsim_ss, length(thresh_methods)),
    thresh_method = rep(thresh_methods_ss, each = nsim),
    thresh_mean = c(thresh_mean),# need c() to turn from matrix into a vector
    thresh_n = c(thresh_n), # number of jackknife iterations that found a threshold
    thresh_n_full = c(full_thresh), # whether a threshold was detected in full data
    thresh_se = c(thresh_se),
    sig_type = rep(sig_criteria, nsim*length(thresh_methods)),
    #thresh_detected = c(thresh_detected),
    ci_n_up = c(ci_n_up),
    ci_n_low = c(ci_n_low),
    ci_em_up = c(ci_em_up),
    ci_em_low = c(ci_em_low),
    thresh_diff = c(c(thresh_mean)-thresh_loc)#,

  )




  # return results
  return(df)
  # return(list(thresh_mean = thresh_mean, thresh_n = thresh_n, thresh_se = thresh_se, thresh_n_mean = thresh_n_mean))

}


# functions for finding thresholds and determining their significance
# common inputs:
#'@param deriv_i values of derivative for which to find a threshold (either the ith jackknifing
#'iteration or the full data)
#'@param xvals values of the driver used to generate the derivative predictions
#'@param sig_type which type of significance criteria to use for determining whether a threshold is
#'significant. Can be `none` (no significance criteria), `jack_quant` (use the quantiles from across
#'all jackknife iterations as upper/lower bounds), or `sim_int` (use the simultaneous intervals
#'calculated by gratia::derivatives for each individual jackknife iteration)
#'@param deriv_up either the 97.5% quantile of the deriv across  all jackknifing iterations
#'(if `sig_type` = `jack_quant`) or the upper simultaneous interval for the ith jackknife iteration
#'(if `sig_type` = `sim_int`) or `NULL` (if `sig_type` = `none`)
#'@param deriv_low either the 2.5% quantile of the deriv across  all jackknifing iterations
#'(if `sig_type` = `jack_quant` or the upper simultaneous interval for the ith jackknife iteration
#'(if `sig_type` = `sim_int`) or `NULL` (if `sig_type` = `none`)


# ##### #### absolute max in second deriv
abs_max_threshF <- function(deriv_i, sig_type, xvals, deriv_up, deriv_low, pkht, pkup){

 # first get the max and min points
  # get median value (shift deriv to relative to median so minpeakheight is relative to median rather than 0)
  med <- median(deriv_i)
  min_all <- findpeaks(-(deriv_i-med), nups =pkup, minpeakheight = pkht)[,2]
  max_all <- findpeaks((deriv_i-med), nups =pkup, minpeakheight = pkht)[,2]

  #min_all <- findpeaks(-deriv_i, nups =pkup, minpeakheight = pkht)[,2]
  #max_all <- findpeaks(deriv_i, nups =pkup, minpeakheight = pkht)[,2]


  if(is.null(min_all)==T & is.null(max_all)==T){ # findpeaks returns NULL if none are found
    abs_max_thresh <- NA#which(1==10)#length(min_all) # abs_max_thresh = integer(0)
  } else{

    if(sig_type=="none"){

      # get the min corresponding to the smallest deriv value
      min1 <- min_all
      if(is.null(min_all)==F){
        min1 <- min_all[which(deriv_i[min_all] == min(deriv_i[min_all]))]
      }

      # get the max corresponding to the largest deriv value
      max1 <- max_all
      if(is.null(max_all)==F){
        max1 <- max_all[which(deriv_i[max_all] == max(deriv_i[max_all]))]
      }

      pks1 <- c(min1, max1)
      pk_choice <- which(abs(deriv_i[pks1])==max(abs(deriv_i[pks1])))[1]

      abs_max_thresh <- xvals[pks1[pk_choice]]

    } # end of if sig_type = none

    if(sig_type%in% c("jack_quant", "sim_int")){

      min2 <- which(1==10) # holding place
      if(is.null(min_all)==F){

        # significant mins: upper interval at the min is less than the lower interval at the lower interval's highest point
        min1 <- min_all[which(deriv_up[min_all] < max(deriv_low))]

        if(length(min1)>0){
          min2 <- min1[which(deriv_i[min1] == min(deriv_i[min1]))]
        }
      } # end of if length(min_all) > 0

      max2 <- which(1==10) # empty holding place
      if(is.null(max_all)==F){

        # significant maxs: lower interval at the max is greater than the upper interval at the upper interval's lowest point
        max1 <- max_all[which(deriv_low[max_all] > min(deriv_up))]

        if(length(max1)>0){
          max2 <- max1[which(deriv_i[max1] == max(deriv_i[max1]))]
        }
      } # end of if length(max_all) > 0

      pks2 <- c(min2, max2)

      if(length(pks2)==0){
        abs_max_thresh <- NA#length(pks2) # abs_max_thresh = integer(0)
      } else{

        pk_choice2 <- which(abs(deriv_i[pks2])==max(abs(deriv_i[pks2])))[1]

        abs_max_thresh <- xvals[pks2[pk_choice2]]

      }

    } # end of if sig_type%in% c("jack_quant", "sim_int")

  } # end of else for length(min_all)>0 or length(max_all)>0

return(abs_max_thresh)

} # end of function

# ##### #### min in second deriv

min_threshF <- function(deriv_i, sig_type, xvals, deriv_up, deriv_low, pkht, pkup){

  # first get the max and min points
  # get median value (shift deriv to relative to median so minpeakheight is relative to median rather than 0)
  med <- median(deriv_i)
  min_all <- findpeaks(-(deriv_i-med), nups =pkup, minpeakheight = pkht)[,2]

 # min_all <- findpeaks(-deriv_i, nups =pkup, minpeakheight = pkht)[,2]


  if(is.null(min_all)==T){
    min_thresh <- NA#length(min_all) # min_thresh = integer(0)
  } else{

    if(sig_type=="none"){

      min_choice <- which(deriv_i[min_all] == min(deriv_i[min_all]))[1]

      min_thresh <- xvals[min_all[min_choice]]

    } # end of if sig_type = none

    if(sig_type%in% c("jack_quant", "sim_int")){

        # significant mins: upper interval at the min is less than the lower interval at the lower interval's highest point
        min1 <- min_all[which(deriv_up[min_all] < max(deriv_low))]

        if(length(min1)>0){
          min2 <- min1[which(deriv_i[min1] == min(deriv_i[min1]))][1]

          min_thresh <- xvals[min2]

        } else{
          min_thresh <- NA#length(min1) # min_thresh = integer(0)
        }


    } # end of if sig_type%in% c("jack_quant", "sim_int")

  } # end of else for length(min_all)>0 or length(max_all)>0

  return(min_thresh)

} # end of function

# ##### #### roots
root_threshF <- function(deriv_i, sig_type, xvals, deriv_up, deriv_low){

  # find the points where the second deriv crosses zero (=changes sign)
  # sign(-10) = -1, sign(10) = 1, sign(0) = 0
  droots <- which(diff(sign(deriv_i))!=0)

  if(length(droots)==0){
    root_thresh <- NA#which(1==10)#droots # root_thresh = integer(0)
  } else{

    if(sig_type=="none"){

      if(length(droots) ==1){

        root_thresh <- (xvals[droots] + xvals[droots + 1])/2

      } else{

        # for each root, calculate the min of the number of indeces to the left and right of it, before
        # the next root (or beginning/end of sequence) is reached

        min_dist <- rep(NA, length(droots))
        max_dist <- rep(NA, length(droots))

        for(rr in 1:length(droots)){

          n_left <- ifelse(rr==1, droots[rr]-1, droots[rr]-droots[rr-1])
          n_right <- ifelse(rr==length(droots), length(deriv_i)-droots[length(droots)], droots[rr+1]-droots[rr])
          min_dist[rr] <- min(n_left, n_right)
          max_dist[rr] <- max(n_left, n_right)

        }

        root_index <- droots[which(min_dist==max(min_dist))]

        if(length(root_index)>1){ # if the max min distance was the same for multiple roots
          min_dist_sub <- which(min_dist==max(min_dist))
          root_index <- droots[which(max_dist[min_dist_sub]==max(max_dist[min_dist_sub]))][1] # pick the one with the largest distance on the other side
        }

        root_thresh <- (xvals[root_index] + xvals[root_index + 1])/2


      } # end of else for if length(droots)>1

    } # end of if sig_type == "none"

    if(sig_type %in% c("jack_quant", "sim_int")){
# first get which roots are significant (those that are closest to where the upper/lower intervals cross zero)
      d_low_roots <- which(diff(sign(deriv_low))!=0)
      d_up_roots <- which(diff(sign(deriv_up))!=0)

      if(length(d_low_roots) > 0 & length(d_up_roots) > 0){ # if both intervals had roots
        # get all the roots that fall between roots in the upper and lower intervals
        # start with the shortest, and pair each of these roots with the closest one in the other interval on the correct side

        if(length(d_low_roots) >= length(d_up_roots)){ # if there were fewer roots in upper CI

          int_mat <- matrix(NA, nrow = length(d_up_roots), ncol = 2) # first col = lower root, 2nd col = corresponding upper root

          for(ii in 1:length(d_up_roots)){

            # if d_up is decreasing at root, d_low root is on the left
            if(deriv_up[d_up_roots[ii]] > deriv_up[d_up_roots[ii] + 1]){
              rootsL <- d_low_roots[which(d_low_roots <= d_up_roots[ii])] # d_low roots to the left

              int_mat[ii, 1] <- ifelse(length(rootsL)==0, NA, rootsL[which(abs(rootsL - d_up_roots[ii])==min(abs(rootsL - d_up_roots[ii])))])[1]
              int_mat[ii, 2] <- d_up_roots[ii]
            }

            # if d_up is increasing at root, d_low root is on the right
            if(deriv_up[d_up_roots[ii]] < deriv_up[d_up_roots[ii] + 1]){
              rootsR <- d_low_roots[which(d_low_roots >= d_up_roots[ii])] # d_low roots to the right

              int_mat[ii, 1] <- ifelse(length(rootsR)==0, NA, rootsR[which(abs(rootsR - d_up_roots[ii])==min(abs(rootsR - d_up_roots[ii])))])[1]
              int_mat[ii, 2] <- d_up_roots[ii]

            }

          }

        } # end of calculating int_mat if length(d_low_roots) >= length(d_up_roots)

        if(length(d_low_roots) < length(d_up_roots)){ # if there were fewer roots in lower CI

          int_mat <- matrix(NA, nrow = length(d_low_roots), ncol = 2) # first col = lower root, 2nd col = corresponding upper root

          for(ii in 1:length(d_low_roots)){

            # if d_low is increasing at root, d_up root is on the left
            if(deriv_low[d_low_roots[ii]] < deriv_low[d_low_roots[ii] + 1]){
              rootsL <- d_up_roots[which(d_up_roots <= d_low_roots[ii])] # d_up roots to the left

              int_mat[ii, 1] <- d_low_roots[ii]
              int_mat[ii, 2] <- ifelse(length(rootsL)==0, NA, rootsL[which(abs(rootsL - d_low_roots[ii])==min(abs(rootsL - d_low_roots[ii])))])[1]

            }

            # if d_low is increasing at root, d_up root is on the right
            if(deriv_low[d_low_roots[ii]] > deriv_low[d_low_roots[ii] + 1]){
              rootsR <- d_up_roots[which(d_up_roots >= d_low_roots[ii])] # d_up roots to the right

              int_mat[ii, 1] <- d_low_roots[ii]
              int_mat[ii, 2] <- ifelse(length(rootsR)==0, NA, rootsR[which(abs(rootsR - d_low_roots[ii])==min(abs(rootsR - d_low_roots[ii])))])[1]


            }
          }

        } # end of calculating int_mat if length(d_low_roots) < length(d_up_roots)

        # get all the non-NA intervals
        na_chk <- which(is.na(int_mat[, 1] + int_mat[, 2])==F)

        # if there weren't any
        if(length(na_chk)==0){
          root_thresh <- NA#which(1==10)#length(na_chk)# root_thresh = integer(0)
        } else{

        # for the non-NA intervals, get the roots in the deriv that fall between the intervals
       # int_mat <- int_mat[na_chk, ]
        int_mat <- matrix(int_mat[na_chk, ], nrow = length(na_chk), ncol = 2) # if na_chk is 1, this will
        #turn into a vector and get incorrect number of dimensions error, so need to force it to
        #be a matrix

        sig_roots <- list(NA)

        for(ij in 1:length(int_mat[,1])){
          # for each interval, get the roots within that interval
          sig_ind <- which(data.table::between(droots, min(int_mat[ij,1], int_mat[ij,2]), max(int_mat[ij,1], int_mat[ij,2]), incbounds=TRUE))
          sig_roots[[ij]] <- ifelse(length(sig_ind)==0, NA, droots[sig_ind]) #droots[sig_ind]


        }

        sig_roots2 <- unlist(sig_roots)

        # get all the locations where the intervals cross zero, then get the roots in these intervals
        # if multiple in an interval, pick the one with the largest consecutive negatives/positives on each side, to be consistent

        na_chk2 <- which(is.na(sig_roots2)==F)

        if(length(na_chk2)==0){
          root_thresh <- NA#which(1==10)#length(na_chk2)# root_thresh = integer(0)
        } else {

          sig_roots2 <- sig_roots2[which(is.na(sig_roots2)==F)]

          if(length(sig_roots2)==1){
            root_thresh <- (xvals[sig_roots2] + xvals[sig_roots2 + 1])/2
          }

          if(length(sig_roots2) > 1){
            # pick the root for which the deriv was consecutively pos and neg on either side for the longest amount
            min_dist <- rep(NA, length(sig_roots2))
            max_dist <- rep(NA, length(sig_roots2))

            for(rr in 1:length(sig_roots2)){

              n_left <- ifelse(rr==1, sig_roots2[rr]-1, sig_roots2[rr]-sig_roots2[rr-1])
              n_right <- ifelse(rr==length(sig_roots2), length(deriv_i)-sig_roots2[length(sig_roots2)], sig_roots2[rr+1]-sig_roots2[rr])
              min_dist[rr] <- min(n_left, n_right)
              max_dist[rr] <- max(n_left, n_right)

            }

            root_index <- sig_roots2[which(min_dist==max(min_dist))]

            if(length(root_index)>1){ # if the max min distance was the same for multiple roots
              min_dist_sub <- which(min_dist==max(min_dist))
              root_index <- sig_roots2[which(max_dist[min_dist_sub]==max(max_dist[min_dist_sub]))][1] # pick the one with the largest distance on the other side
            }

            root_thresh <- (xvals[root_index] + xvals[root_index + 1])/2

          } # end of if length(sig_roots2) > 1

        } # end of else for length(na_chk2) > 0

      } # end of else for length(na_chk) > 0

        # end of if length(d_low_roots) > 0 & length(d_up_roots) > 0
      } else{ # if the upper and lower intervals never crossed 0

        root_thresh <- NA#which(1==10)#length(d_low_roots)# root_thresh = integer(0)
      }


} # end of if sig_type %in% c("jack_quant", "sim_int")


  } # end of else{} for length(droots) >0

 return(root_thresh)

} # end of function



# function for standardizing data
zfun <- function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
}










