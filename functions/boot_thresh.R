#' function for evaluating thresholds using bootrapping methods from the Holsman et al. thresholds
#' workshop (only difference is using the gratia package to calculate derivatives):
#' https://github.com/kholsman/EBM_Holsman_NatComm/blob/master/R/sub_fun/threshold.R
#' From Holsman etal 2020
#' Holsman, K. K., A. Haynie, A. Hollowed, J.C.P. Reum, K. Aydin, A. J. Hermann, W. Cheng,
#' A. Faig, J. N. Ianelli, K. Kearney, A. E. Punt, 2020.Ecosystem-based fisheries
#'  management forestalls climate-driven collapse, Nature Communications.

#' required packages:
#' mgcv (for fitting the gam)
#' gratia (for derivative calculations),
#' quantmod (for the findPeaks function to find max/mins in the derivatives),
#' data.table (for between function to determine whether CIs contain 0)


#' function inputs:
#' @param simdt data frame with the simulated data (output of simfun2)
#' @param xvals vector of driver values to use for generating predictions from the gams and
#' calculating their derivatives
#' @param thresh_methods vector of methods to use for estimating the threshold location.
#' Can include `abs_max_d2` (location where 2nd deriv of gam is furthest from zero), `min_d2`
#' (location where 2nd deriv is at its minimum- i.e., most negative- value), `zero_d2`
#' (location where 2nd deriv of gam crosses zero), or `zero_d1` (location where first deriv of gam
#' crosses zero/is zero). Default is to do all of the above
#' @param boot_n number of bootstrap iterations to run, default is 500
#' @param boot_nobs size of subsample taken each bootstrap iteration
#' @param boot_seed initial seed for the random sampling in the bootstrapping
#' @param knots number of knots to use in the gam fitting, defaults to 4
#' @param smooth_type type of smooth used with gams, defaults to "tp"
#' @param span smoothing step/span for thresholds, values closer to 1 = more smoothing
#' (Holsman et al. 2020 used 0.1 as default)


boot_thresh <- function(simdt, xvals, thresh_methods = c("abs_max_d2", "min_d2", "zero_d2", "zero_d1"), boot_nobs, boot_n = 500, boot_seed = 204, knots = 4, smooth_type = "tp", span = 0.1){

  # get the number of simulations that were run
  nsim <- length(unique(simdt$sim))

  # get the true value of the threshold
  thresh_loc <- simdt$thresh_loc[1]

  # turn the driver values into a dataframe
  xdt <- data.frame(
    driver = xvals
  )

  #for(mm in 1:length(thresh_methods)){

  #thresh_method <- thresh_methods[mm]

  # make the holding matrices to the value of the threshold and whether a threshold was detected (0 or 1)
  # for each simulation
  thresh_val <- matrix(NA, nrow = nsim, ncol = length(thresh_methods))# estimated threshold value
  thresh_n <- matrix(NA, nrow = nsim, ncol = length(thresh_methods))# whether a threshold was detected

  for(i in 1:nsim){ # for each simulation

    # subset the data for the ith simulation
    datIN <- simdt[which(simdt$sim == i), ]

    # store the means and sd's of the driver
    #x_means[i] <- mean(datIN$driver)
    #x_sds[i] <- sd(datIN$driver)

    # make data set for the jackknife resampling
    dd <- datIN
    #dd$driver <- round(dd$driver, 3)
    dd <- dd[ ,c("driver", "obs_response")]
    dd$num <- c(1:length(dd$driver))


    # make holding matrices for the derivatives
    D1 <- matrix(NA,boot_n,length(xvals)) # first derivative
    D2 <- matrix(NA,boot_n,length(xvals)) # second derivative

    # set the seed
    set.seed(boot_seed+i^2)

    # perform bootstrapping
    for(int in 1:boot_n){ # for each bootstrap iteration

      # get bootstrapped sub-sample
      nobs <- length(dd$num) # number of observations in the data set
      if(boot_nobs > nobs){ # if number of observations to sample for bootstrap is greater than nobs
        boot_nobs   <- nobs} # set number of bootstrap samples equal to nobs

      bootd <- dd[sample(dd$num, boot_nobs, replace = TRUE),]# sample boot_nobs rows (with replacement) from the dd data set and name this subsetted data frame bootd


      gami <- gam(obs_response~s(driver,k=knots,bs=smooth_type),data = bootd) # fit a gam to the bootd data set

      D1iunsm <- gratia::derivatives(gami, data = xdt, order = 1)$derivative # use the derivatives function from gratia to approximate 1st derivative for the values of the driver in xvals
      D2iunsm <- gratia::derivatives(gami, data = xdt, order = 2)$derivative # use the derivatives function from gratia to approximate 2nd derivative for the values of the driver in xvals

      # smooth these (without smoothing, the 2nd deriv seems to be super noisy)
      #D1i <- predict(loess(D1iunsm ~ driver, data=data.frame(D1iunsm = D1iunsm, driver = xvals), span=span))
      #D2i <- predict(loess(D2iunsm ~ driver, data=data.frame(D2iunsm = D2iunsm, driver = xvals), span=span))

      # store these
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

      # smooth these
      d2_low <- predict(loess(d2unsm_low ~ driver, data=data.frame(d2unsm_low = d2unsm_low, driver = xvals), span=span))
      d2_mn <- predict(loess(d2unsm_mn ~ driver, data=data.frame(d2unsm_mn = d2unsm_mn, driver = xvals), span=span))
      d2_up <- predict(loess(d2unsm_up ~ driver, data=data.frame(d2unsm_up = d2unsm_up, driver = xvals), span=span))

      d1_low <- predict(loess(d1unsm_low ~ driver, data=data.frame(d1unsm_low = d1unsm_low, driver = xvals), span=span))
      d1_mn <- predict(loess(d1unsm_mn ~ driver, data=data.frame(d1unsm_mn = d1unsm_mn, driver = xvals), span=span))
      d1_up <- predict(loess(d1unsm_up ~ driver, data=data.frame(d1unsm_up = d1unsm_up, driver = xvals), span=span))


      # then use these to calculate the thresholds for all the different methods in thresh_methods
      for(mm in 1:length(thresh_methods)){

        thresh_method <- thresh_methods[mm]

        if(thresh_method =="abs_max_d2"){

            # peaks in 2nd deriv
            d2pks <- c(findPeaks(d2_mn), findPeaks(-d2_mn))

            # locations where 2nd deriv is significantly different from 0
            signifD2 <- which(!data.table::between(0, d2_low, d2_up, incbounds=TRUE))

            sigd2pks <- intersect(signifD2, d2pks)

            thresh_val[i, mm] <- NA # holding place

            if(length(sigd2pks)>0){
              # choose the detected threshold(s) for which the derivative is significantly different from zero
              thresh_val[i, mm] <- xvals[sigd2pks[which(abs(d2_mn[sigd2pks])==max(abs(d2_mn[sigd2pks])))]]
            }

            # store the number of thresholds detected
            thresh_n[i, mm] <- ifelse(is.na(thresh_val[i, mm])==F, 1, 0)



        }

        if(thresh_method =="min_d2"){

            # minima in second deriv
            d2mins <- findPeaks(-d2_mn)

            # locations where 2nd deriv is significantly different from 0
            signifD2 <- which(!data.table::between(0, d2_low, d2_up, incbounds=TRUE))

            sigd2mins <- intersect(signifD2, d2mins)

            thresh_val[i, mm] <- NA # holding place

            if(length(sigd2mins)>0){
              # choose the detected threshold(s) for which the derivative is significantly different from zero
              thresh_val[i, mm] <- xvals[sigd2mins[which(abs(d2_mn[sigd2mins])==max(abs(d2_mn[sigd2mins])))]]
            }

            # store the number of thresholds detected
            thresh_n[i, mm] <- ifelse(is.na(thresh_val[i, mm])==F, 1, 0)

        }

        if(thresh_method =="zero_d2"){

          # roots of second deriv
          d2roots <- which(diff(sign(d2_mn))!=0)

            thresh_val[i, mm] <- NA # holding place


            if(length(d2roots)>0){
              # locations where 2nd deriv is NOT significantly different from 0
              nsignifD2 <- which(data.table::between(0, d2_low, d2_up, incbounds=TRUE))

              if(length(nsignifD2)>0){
                # select the roots where the deriv was not significantly different from zero at at least one of the two points where the sign change occurred
                # (e.g., if d2root = 3, it means d2[3] and d2[4] had different signs so 2nd deriv crossed zero between them, and say
                # this root is significant if either d2[3] or d2[4] (or both) are not significantly different from zero
                nsigd2rootsL <- intersect(nsignifD2, d2roots) # left of crossing is not sig diff from 0
                nsigd2rootsR <- intersect(nsignifD2, d2roots + 1) # right of crossing is not sig diff from 0
                nsigd2roots <- unique(nsigd2rootsL, nsigd2rootsR-1)

                if(length(nsigd2roots)>0){
                  # get the sum of the differences between zero and the upper quantile, lower quantile, and mean of second deriv on both sides of the root
                  summdiff <- abs(d2_low[nsigd2roots]-0) + abs(d2_up[nsigd2roots]-0) + abs(d2_mn[nsigd2roots]-0) + abs(d2_low[nsigd2roots+1]-0) + abs(d2_up[nsigd2roots+1]-0) + abs(d2_mn[nsigd2roots+1]-0)
                  index<-nsigd2roots[which(summdiff==min(summdiff))] # get value where this sum is minimized

                  # save the threshold value
                  thresh_val[i, mm] <- (xvals[index] + xvals[index + 1])/2
                }
              }
            }

            # store the number of thresholds detected
            thresh_n[i, mm] <- ifelse(is.na(thresh_val[i, mm])==F, 1, 0)


        }


        if(thresh_method =="zero_d1"){

          # get roots in first derivative
          # find the points where the first deriv crosses zero (=changes sign)
          d1roots <- which(diff(sign(d1_mn))!=0)

          thresh_val[i, mm] <- NA # holding place

            if(length(d1roots)>0){
              # locations where 1st deriv is NOT significantly different from 0
              nsignifD1 <- which(data.table::between(0, d1_low, d1_up, incbounds=TRUE))

              if(length(nsignifD1)>0){

                # select the roots where the deriv was not significantly different from zero at at least one of the two points where the sign change occurred
                nsigd1rootsL <- intersect(nsignifD1, d1roots) # left of crossing is not sig diff from 0
                nsigd1rootsR <- intersect(nsignifD1, d1roots + 1) # right of crossing is not sig diff from 0
                nsigd1roots <- unique(nsigd1rootsL, nsigd1rootsR-1)

                if(length(nsigd1roots)>0){
                  # get the sum of the differences between zero and the upper quantile, lower quantile, and int'th second deriv on both sides of the root
                  summdiff <- abs(d1_low[nsigd1roots]-0) + abs(d1_up[nsigd1roots]-0) + abs(d1_mn[nsigd1roots]-0) + abs(d1_low[nsigd1roots+1]-0) + abs(d1_up[nsigd1roots+1]-0) + abs(d1_mn[nsigd1roots+1]-0)
                  index<-nsigd1roots[which(summdiff==min(summdiff))] # get value where this sum is minimized

                  # save the threshold value
                  thresh_val[i, mm] <- (xvals[index] + xvals[index + 1])/2

                }
              }
            }

            # store the number of thresholds detected
            thresh_n[i, mm] <- ifelse(is.na(thresh_val[i, mm])==F, 1, 0)


        }

      } # end of threshold calculations



  } # end of iterations for each simulation i





  df <- data.frame(
    sim = rep(c(1:nsim), length(thresh_methods)),
    thresh_method = rep(thresh_methods, each = nsim),
    thresh_vals = c(thresh_val),
    thresh_n = c(thresh_n),
    thresh_diff = c(c(thresh_val)-thresh_loc)#,
    #thresh_diffN = c((thresh_mean-thresh_loc)/ rep(x_sds, length(thresh_methods)))
  )



  # return results
  return(df)
  # return(list(thresh_mean = thresh_mean, thresh_n = thresh_n, thresh_se = thresh_se, thresh_n_mean = thresh_n_mean))

}
