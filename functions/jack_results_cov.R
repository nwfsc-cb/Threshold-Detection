
#' same as jack_results function but for GAMs that include a covariate
#' function to calculate and return the results (gam predictions, first and second derivs)
#' of the jackknifing done by the jack_thresh_cov function for gams that include a covariate
#'
#' required packages:
#' mgcv (for fitting the gam)
#' gratia (for derivative calculations),


#' function inputs:
#' @param simdt data frame with the simulated data (output of simfun2)
#' @param x_pred optional vector of driver values to use for generating predictions from the gams and
#' calculating their derivatives
#' @param cov_val value of the covariate (standardized scale) to use for generating gam predictions
#' @param xlength if x_pred = NA, make set of xvals for each dataset with length equal to max(100, xlength*range(driver))
#' @param z_scale if TRUE, standardize the data before fitting (default), if FALSE use raw data
#' @param sim_choice the simulation(s) in simdt for which to return the jackknifing results
#' @param knots number of knots to use in the gam fitting, defaults to 4
#' @param smooth_type type of smooth used with gams, defaults to "tp"
#' @param span smoothing step/span for thresholds, values closer to 1 = more
#' smoothing (Holsman et al. 2020 used 0.1 as default)
#' @param eps_val eps value (the finite difference) for the gratia::derivatives function, defaults to 5*10^-6


#' note: make sure knots, smooth_type, and span all have the same values as used in jack_thresh


jack_results_cov <- function(simdt, sim_choice, x_pred, cov_val = 0, xlength = 50, z_scale = TRUE, knots = 4,
                         smooth_type = "tp", eps_val = 5*10^-6, span = 0.1){

  # get the simulations for which to do the jackknifing
  nsim <- sim_choice

  # turn the driver values into a dataframe
  if(is.na(x_pred[1])==F){

    xvals <- x_pred

    xdt <- data.frame(
      driver = x_pred,
      cov = rep(cov_val, length(x_pred))
    )
  }


  for(i in 1:length(nsim)){ # for each simulation


    # subset the data for the ith simulation and standardize
    datIN <- simdt[which(simdt$sim == nsim[i]), ]

    if(z_scale==TRUE){
    datIN$driver <- zfun(datIN$driver)
    datIN$obs_response <- zfun(datIN$obs_response)
    datIN$cov <- zfun(datIN$cov)
    }

    # get the xvals
    if(is.na(x_pred[1])==T){

      rangei <- max(datIN$driver, na.rm = T) - min(datIN$driver, na.rm = T)
      xvals <- seq(from = min(datIN$driver, na.rm = T), to = max(datIN$driver, na.rm = T), length.out =max(100, rangei*xlength))# make sure length is at least 100

      xdt <- data.frame(
        driver = xvals,
        cov = rep(cov_val, length(xvals))
      )

    }


    # first fit a gam to the whole data set

    gam_full <- gam(obs_response~s(driver,k=knots,bs=smooth_type) + s(cov,k=knots,bs=smooth_type),data = datIN) # fit a gam to the jackd data set


     #if(z_scale==TRUE){
      #gam_full <- gam(obs_response~ -1 + s(driver,k=knots,bs=smooth_type),data = datIN) # fit a gam to the jackd data set
    #} else{
     # gam_full <- gam(obs_response~s(driver,k=knots,bs=smooth_type),data = datIN) # fit a gam to the jackd data set
    #}

    # gam predictions and simultaneous intervals
    pred_full <- gratia::fitted_values(gam_full, data = xdt)
    pred_full_up <- predict(loess(pred ~ driver, data=data.frame(pred = pred_full$upper, driver = xvals), span=span))
    pred_full_mn <- predict(loess(pred ~ driver, data=data.frame(pred = pred_full$fitted, driver = xvals), span=span))
    pred_full_low <- predict(loess(pred ~ driver, data=data.frame(pred = pred_full$lower, driver = xvals), span=span))

    # first deriv
    D1_full <- gratia::derivatives(gam_full,term = "s(driver)", data = xdt, order = 1, eps = eps_val)
    D1_full_up <- predict(loess(d1 ~ driver, data=data.frame(d1 = D1_full$upper, driver = xvals), span=span))
    D1_full_mn <- predict(loess(d1 ~ driver, data=data.frame(d1 = D1_full$derivative, driver = xvals), span=span))
    D1_full_low <- predict(loess(d1 ~ driver, data=data.frame(d1 = D1_full$lower, driver = xvals), span=span))

    # second deriv
    D2_full <- gratia::derivatives(gam_full,term = "s(driver)", data = xdt, order = 2, eps = eps_val)
    D2_full_up <- predict(loess(d2 ~ driver, data=data.frame(d2 = D2_full$upper, driver = xvals), span=span))
    D2_full_mn <- predict(loess(d2 ~ driver, data=data.frame(d2 = D2_full$derivative, driver = xvals), span=span))
    D2_full_low <- predict(loess(d2 ~ driver, data=data.frame(d2 = D2_full$lower, driver = xvals), span=span))


    # save results
    response_full_df <- data.frame(
      sim = rep(nsim[i], length(xvals)),
      driver = xvals,
      output = rep("response", length(xvals)),
      low = pred_full_low,
      mn = pred_full_mn,
      up = pred_full_up
    )

    d2_full_df <- data.frame(
      sim = rep(nsim[i], length(xvals)),
      driver = xvals,
      output = rep("d2", length(xvals)),
      low = D2_full_low,
      mn = D2_full_mn,
      up = D2_full_up
    )


    d1_full_df <- data.frame(
      sim = rep(nsim[i], length(xvals)),
      driver = xvals,
      output = rep("d1", length(xvals)),
      low = D1_full_low,
      mn = D1_full_mn,
      up = D1_full_up
    )

    full_df <- rbind(response_full_df, d2_full_df, d1_full_df)


    # make data set for the jackknife resampling
    dd <- datIN
    #dd$driver <- round(dd$driver, 3)
    dd <- dd[ ,c("driver", "obs_response", "cov")]

    # make holding matrix for the gam predictions
    pred <- matrix(NA,length(dd[,1]),length(xvals))
    pred_lowi <- matrix(NA,length(dd[,1]),length(xvals))
    pred_upi <- matrix(NA,length(dd[,1]),length(xvals))

    # make holding matrices for the derivatives
    D1 <- matrix(NA,length(dd[,1]),length(xvals)) # first derivative
    D1_lowi <- matrix(NA,length(dd[,1]),length(xvals))
    D1_upi <- matrix(NA,length(dd[,1]),length(xvals))

    D2 <- matrix(NA,length(dd[,1]),length(xvals)) # second derivative
    D2_lowi <- matrix(NA,length(dd[,1]),length(xvals))
    D2_upi <- matrix(NA,length(dd[,1]),length(xvals))

    # make matrix to keep track of the jackknife iterations
    jack_index <- matrix(NA,length(dd[,1]),length(xvals))

    # perform jackknife
    for(int in 1:length(dd[,1])){ # for each observation
      # remove this observation from the dataset
      jackd <- dd[-int, ]

      gami <- gam(obs_response~s(driver,k=knots,bs=smooth_type) + s(cov,k=knots,bs=smooth_type),data = jackd) # fit a gam to the jackd data set

       #if(z_scale==TRUE){
        #gami <- gam(obs_response~ -1+ s(driver,k=knots,bs=smooth_type),data = jackd) # fit a gam to the jackd data set
     # } else{
        #gami <- gam(obs_response~s(driver,k=knots,bs=smooth_type),data = jackd) # fit a gam to the jackd data set
     # }

      #predi <- stats::predict(gami,se.fit=TRUE, newdata =xdt)$fit # get the gam's predictions
      #predis <- predict(loess(predi ~ driver, data=data.frame(predi = predi, driver = xvals), span=span))

      #D1i <- gratia::derivatives(gami, data = xdt, order = 1)$derivative # use the derivatives function from gratia to approximate 1st derivative for the values of the driver in xvals
      #D2i <- gratia::derivatives(gami, data = xdt, order = 2)$derivative # use the derivatives function from gratia to approximate 2nd derivative for the values of the driver in xvals
      #D1is <- predict(loess(D1i ~ driver, data=data.frame(D1i = D1i, driver = xvals), span=span))
      #D2is <- predict(loess(D2i ~ driver, data=data.frame(D2i = D2i, driver = xvals), span=span))

      # gam predictions and simultaneous intervals
      predi <- gratia::fitted_values(gami, data = xdt)
      predi_up <- predict(loess(pred ~ driver, data=data.frame(pred = predi$upper, driver = xvals), span=span))
      predi_mn <- predict(loess(pred ~ driver, data=data.frame(pred = predi$fitted, driver = xvals), span=span))
      predi_low <- predict(loess(pred ~ driver, data=data.frame(pred = predi$lower, driver = xvals), span=span))

      # first deriv
      D1i <- gratia::derivatives(gami, term = "s(driver)", data = xdt, order = 1, eps = eps_val)
      D1i_up <- predict(loess(d1 ~ driver, data=data.frame(d1 = D1i$upper, driver = xvals), span=span))
      D1i_mn <- predict(loess(d1 ~ driver, data=data.frame(d1 = D1i$derivative, driver = xvals), span=span))
      D1i_low <- predict(loess(d1 ~ driver, data=data.frame(d1 = D1i$lower, driver = xvals), span=span))


      D2i <- gratia::derivatives(gami, term = "s(driver)", data = xdt, order = 2, eps = eps_val)
      D2i_up <- predict(loess(d2 ~ driver, data=data.frame(d2 = D2i$upper, driver = xvals), span=span))
      D2i_mn <- predict(loess(d2 ~ driver, data=data.frame(d2 = D2i$derivative, driver = xvals), span=span))
      D2i_low <- predict(loess(d2 ~ driver, data=data.frame(d2 = D2i$lower, driver = xvals), span=span))


      # store these
      pred[int, ] <- predi_mn
      pred_lowi[int, ] <- predi_low
      pred_upi[int, ] <- predi_up

      D1[int, ] <- D1i_mn
      D1_lowi[int, ] <- D1i_low
      D1_upi[int, ] <- D1i_up

      D2[int, ] <- D2i_mn
      D2_lowi[int, ] <- D2i_low
      D2_upi[int, ] <- D2i_up

      # store the jackknife iteration
      jack_index[int, ] <- rep(int, length(xvals))

    }


    # get the 2.5%, 50%, and 97.5% quantiles of the jackknifed gam predictions and first and second derivatives
    response_df <- data.frame(
      sim = rep(nsim[i], length(xvals)),
      driver = xvals,
      output = rep("response", length(xvals)),
      low = apply(pred,2,quantile,probs=0.025),
      mn = apply(pred,2,quantile,probs=0.5),
      up = apply(pred,2,quantile,probs=0.975)
    )

    d2_df <- data.frame(
      sim = rep(nsim[i], length(xvals)),
      driver = xvals,
      output = rep("d2", length(xvals)),
      low = apply(D2,2,quantile,probs=0.025),
      mn = apply(D2,2,quantile,probs=0.5),
      up = apply(D2,2,quantile,probs=0.975)
    )


      d1_df <- data.frame(
      sim = rep(nsim[i], length(xvals)),
      driver = xvals,
      output = rep("d1", length(xvals)),
      low = apply(D1,2,quantile,probs=0.025),
      mn = apply(D1,2,quantile,probs=0.5),
      up = apply(D1,2,quantile,probs=0.975)
      )

      summ_df <- rbind(response_df, d2_df, d1_df)

      # store results of each jackknifing iteration as well as just the quantiles
      # convert the storage matrixes to data frames (first need to transpose them)
      r_all <- c(t(pred)) # gam predictions (response)
      r_all_up <- c(t(pred_upi)) # upper simultaneous interval
      r_all_low <- c(t(pred_lowi)) # lower simultaneous interval
      D1_all <- c(t(D1)) # first deriv
      D1_all_up <- c(t(D1_upi)) # upper simultaneous interval
      D1_all_low <- c(t(D1_lowi)) # lower simultaneous interval
      D2_all <- c(t(D2)) # second deriv
      D2_all_up <- c(t(D2_upi)) # upper simultaneous interval
      D2_all_low <- c(t(D2_lowi)) # lower simultaneous interval

      int_all <- c(t(jack_index)) # indeces for jackknifing iterations

      driver_all <- rep(xvals, length(dd[,1])) # driver values

      # turn into a dataframe
      ind_df <- data.frame(
        sim = rep(nsim[i], length(int_all)),
        jack_int = int_all,
        response = r_all,
        response_up = r_all_up,
        response_low = r_all_low,
        d1 = D1_all,
        d1_up = D1_all_up,
        d1_low = D1_all_low,
        d2= D2_all,
        d2_up = D2_all_up,
        d2_low = D2_all_low
      )

  # store results for all simulations
  if(i == 1){ # if this was the first simulation
  full_dfs <- full_df
  summ_dfs <- summ_df
  ind_dfs <- ind_df
  } else {
    full_dfs <- rbind(full_df, full_dfs)
    summ_dfs <- rbind(summ_df, summ_dfs)
    ind_dfs <- rbind(ind_df, ind_dfs)
  }

  }

  # return results
  return(list(summ_dfs = summ_dfs, ind_dfs = ind_dfs, full_dfs = full_dfs))

}



