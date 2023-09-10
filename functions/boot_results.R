#' function to calculate and return the results (gam predictions, first and second derivs)
#' of the bootstrapping done by the boot_thresh function (boot_thresh just returns results
#' of the threshold calculations; thresh_results returns the gam fits and derivatives that went into
#' these calculations)
#'
#' required packages:
#' mgcv (for fitting the gam)
#' gratia (for derivative calculations),
#' quantmod (for the findPeaks function to find max/mins in the derivatives),
#' data.table (for between function to determine whether CIs contain 0)


#' function inputs:
#' @param simdt data frame with the simulated data (output of simfun2)
#' @param xvals vector of driver values to use for generating predictions from the gams and
#' calculating their derivatives
#' @param sim_choice the simulation(s) in simdt for which to return the bootstrapping results
#' @param boot_n number of bootstrap iterations to run, default is 500
#' @param boot_nobs size of subsample taken each bootstrap iteration
#' @param boot_seed initial seed for the random sampling in the bootstrapping
#' @param knots number of knots to use in the gam fitting, defaults to 4
#' @param smooth_type type of smooth used with gams, defaults to "tp"
#' @param span smoothing step/span for thresholds, values closer to 1 = more
#' smoothing (Holsman et al. 2020 used 0.1 as default)

#' note: make sure boot_n, boot_nobs, boot_seed, knots, smooth_type, and span all have the same values as used in boot_thresh


boot_results <- function(simdt, xvals, sim_choice, boot_nobs, boot_n = 500, boot_seed = 204, knots = 4, smooth_type = "tp", span = 0.1){

  # get the number of simulations that were run
  #nsim <- length(unique(simdt$sim))

  # get the simulations for which to do the bootstrapping
  nsim <- sim_choice

  # turn the driver values into a dataframe
  xdt <- data.frame(
    driver = xvals
  )


  for(i in 1:length(nsim)){ # for each simulation

    # subset the data for the ith simulation
    datIN <- simdt[which(simdt$sim == nsim[i]), ]

    # make data set for the jackknife resampling
    dd <- datIN
    #dd$driver <- round(dd$driver, 3)
    dd <- dd[ ,c("driver", "obs_response")]
    dd$num <- c(1:length(dd$driver))

    # holding matrix for gam predictions
    pred <- matrix(NA,boot_n,length(xvals))

    # make holding matrices for the derivatives
    D1 <- matrix(NA,boot_n,length(xvals)) # first derivative
    D2 <- matrix(NA,boot_n,length(xvals)) # second derivative

    # set the seed
    match_i <- which(simdt$sim == nsim[i])
    set.seed(boot_seed+match_i^2)

    # perform bootstrapping
    for(int in 1:boot_n){ # for each bootstrap iteration

      # get bootstrapped sub-sample
      nobs <- length(dd$num) # number of observations in the data set
      if(boot_nobs > nobs){ # if number of observations to sample for bootstrap is greater than nobs
        boot_nobs   <- nobs} # set number of bootstrap samples equal to nobs

      bootd <- dd[sample(dd$num, boot_nobs, replace = TRUE),]# sample boot_nobs rows (with replacement) from the dd data set and name this subsetted data frame bootd


      gami <- gam(obs_response~s(driver,k=knots,bs=smooth_type),data = bootd) # fit a gam to the bootd data set

      prediunsm <- stats::predict(gami,se.fit=TRUE, newdata =xdt)$fit # get the gam's predictions

      D1iunsm <- gratia::derivatives(gami, data = xdt, order = 1)$derivative # use the derivatives function from gratia to approximate 1st derivative for the values of the driver in xvals
      D2iunsm <- gratia::derivatives(gami, data = xdt, order = 2)$derivative # use the derivatives function from gratia to approximate 2nd derivative for the values of the driver in xvals

      # smooth these (without smoothing, the 2nd deriv seems to be super noisy)
      #D1i <- predict(loess(D1iunsm ~ driver, data=data.frame(D1iunsm = D1iunsm, driver = xvals), span=span))
      #D2i <- predict(loess(D2iunsm ~ driver, data=data.frame(D2iunsm = D2iunsm, driver = xvals), span=span))

      # store these
      pred[int, ] <- prediunsm
      D1[int, ] <- D1iunsm
      D2[int, ] <- D2iunsm



    } # end of bootstrap iterations


    # get the 2.5%, 50%, and 97.5% quantiles of the jackknifed first and second derivatives and smootht them
    d2unsm_low <- apply(D2,2,quantile,probs=0.025) # 2 in apply() means apply quantile() across columns of D2 matrix
    d2unsm_mn <- apply(D2,2,quantile,probs=0.5)
    d2unsm_up <- apply(D2,2,quantile,probs=0.975)

    d1unsm_low <- apply(D1,2,quantile,probs=0.025)
    d1unsm_mn <- apply(D1,2,quantile,probs=0.5)
    d1unsm_up <- apply(D1,2,quantile,probs=0.975)

    predunsm_low <- apply(pred,2,quantile,probs=0.025)
    predunsm_mn <- apply(pred,2,quantile,probs=0.5)
    predunsm_up <- apply(pred,2,quantile,probs=0.975)

    # smooth these
    d2_low <- predict(loess(d2unsm_low ~ driver, data=data.frame(d2unsm_low = d2unsm_low, driver = xvals), span=span))
    d2_mn <- predict(loess(d2unsm_mn ~ driver, data=data.frame(d2unsm_mn = d2unsm_mn, driver = xvals), span=span))
    d2_up <- predict(loess(d2unsm_up ~ driver, data=data.frame(d2unsm_up = d2unsm_up, driver = xvals), span=span))

    d1_low <- predict(loess(d1unsm_low ~ driver, data=data.frame(d1unsm_low = d1unsm_low, driver = xvals), span=span))
    d1_mn <- predict(loess(d1unsm_mn ~ driver, data=data.frame(d1unsm_mn = d1unsm_mn, driver = xvals), span=span))
    d1_up <- predict(loess(d1unsm_up ~ driver, data=data.frame(d1unsm_up = d1unsm_up, driver = xvals), span=span))

    pred_low <- predict(loess(predunsm_low ~ driver, data=data.frame(predunsm_low = predunsm_low, driver = xvals), span=span))
    pred_mn <- predict(loess(predunsm_mn ~ driver, data=data.frame(predunsm_mn = predunsm_mn, driver = xvals), span=span))
    pred_up <- predict(loess(predunsm_up ~ driver, data=data.frame(predunsm_up = predunsm_up, driver = xvals), span=span))


    # save these all as data frames
    response_df <- data.frame(
      sim = rep(nsim[i], length(xvals)),
      driver = xvals,
      output = rep("response", length(xvals)),
      low = pred_low,
      mn = pred_mn,
      up = pred_up
    )

    d2_df <- data.frame(
      sim = rep(nsim[i], length(xvals)),
      driver = xvals,
      output = rep("d2", length(xvals)),
      low = d2_low,
      mn = d2_mn,
      up = d2_up
    )


    d1_df <- data.frame(
      sim = rep(nsim[i], length(xvals)),
      driver = xvals,
      output = rep("d1", length(xvals)),
      low = d1_low,
      mn = d1_mn,
      up = d1_up
    )

    summ_df <- rbind(response_df, d2_df, d1_df)

    # store results for all simulations
    if(i == 1){ # if this was the first simulation
      summ_dfs <- summ_df
    } else {
      summ_dfs <- rbind(summ_df, summ_dfs)
    }

  } # end of iterations for each simulation i


  # return results
  return(summ_dfs)
  # return(list(thresh_mean = thresh_mean, thresh_n = thresh_n, thresh_se = thresh_se, thresh_n_mean = thresh_n_mean))

}

