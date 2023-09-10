#' function to calculate and return the results (gam predictions, first and second derivs)
#' of the jackknifing done by the jack_thresh function (jack_thresh just returns summary statistics
#' of the threshold calculations; jack_results returns the gam fits and derivatives that went into
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
#' @param sim_choice the simulation(s) in simdt for which to return the jackknifing results
#' @param knots number of knots to use in the gam fitting, defaults to 4
#' @param smooth_type type of smooth used with gams, defaults to "tp"
#' @param span smoothing step/span for thresholds, values closer to 1 = more
#' smoothing (Holsman et al. 2020 used 0.1 as default)

#' note: make sure knots, smooth_type, and span all have the same values as used in jack_thresh


jack_results <- function(simdt, xvals, sim_choice, knots = 4, smooth_type = "tp", span = 0.1){

  # get the number of simulations that were run
  #nsim <- length(unique(simdt$sim))

  # get the simulations for which to do the jackknifing
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

    # make holding matrix for the gam predictions
    pred <- matrix(NA,length(dd[,1]),length(xvals))

    # make holding matrices for the derivatives
    D1 <- matrix(NA,length(dd[,1]),length(xvals)) # first derivative
    D2 <- matrix(NA,length(dd[,1]),length(xvals)) # second derivative

    # make matrix to keep track of the jackknife iterations
    jack_index <- matrix(NA,length(dd[,1]),length(xvals))

    # perform jackknife
    for(int in 1:length(dd[,1])){ # for each observation
      # remove this observation from the dataset
      jackd <- dd[-int, ]

      gami <- gam(obs_response~s(driver,k=knots,bs=smooth_type),data = jackd) # fit a gam to the bootd data set
      predi <- stats::predict(gami,se.fit=TRUE, newdata =xdt)$fit # get the gam's predictions

      predis <- predict(loess(predi ~ driver, data=data.frame(predi = predi, driver = xvals), span=span))

      D1i <- gratia::derivatives(gami, data = xdt, order = 1)$derivative # use the derivatives function from gratia to approximate 1st derivative for the values of the driver in xvals
      D2i <- gratia::derivatives(gami, data = xdt, order = 2)$derivative # use the derivatives function from gratia to approximate 2nd derivative for the values of the driver in xvals

      D1is <- predict(loess(D1i ~ driver, data=data.frame(D1i = D1i, driver = xvals), span=span))

      D2is <- predict(loess(D2i ~ driver, data=data.frame(D2i = D2i, driver = xvals), span=span))

      # store these
      pred[int, ] <- predis
      D1[int, ] <- D1is
      D2[int, ] <- D2is

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
      D1_all <- c(t(D1)) # first deriv
      D2_all <- c(t(D2)) # second deriv

      int_all <- c(t(jack_index)) # indeces for jackknifing iterations

      driver_all <- rep(xvals, length(dd[,1])) # driver values

      # turn into a dataframe
      ind_df <- data.frame(
        sim = rep(nsim[i], length(int_all)),
        jack_int = int_all,
        response = r_all,
        d1 = D1_all,
        d2= D2_all
      )

  # store results for all simulations
  if(i == 1){ # if this was the first simulation
  summ_dfs <- summ_df
  ind_dfs <- ind_df
  } else {
    summ_dfs <- rbind(summ_df, summ_dfs)
    ind_dfs <- rbind(ind_df, ind_dfs)
  }

  }

  # return results
  return(list(summ_dfs = summ_dfs, ind_dfs = ind_dfs))

}
