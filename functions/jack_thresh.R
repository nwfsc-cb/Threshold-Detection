#' function for evaluating thresholds using jackknifing
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
#' @param sig_criteria whether to use jackknifed CIs as significance criteria when evaluating whether
#' a threshold is detected (TRUE or FALSE)
#' @param thresh_choice if multiple thresholds are detected in a single jackknife iteration, which one
#' should be selected? Options are  or `mean` (default; closest to mean threshold across all
#' iterations) or `median` (closest to median threshold across all iterations)
#' @param knots number of knots to use in the gam fitting, defaults to 4
#' @param smooth_type type of smooth used with gams, defaults to "tp"
#' @param span smoothing step/span for thresholds, values closer to 1 = more smoothing
#' (Holsman et al. 2020 used 0.1 as default)

# what to do about potential for multiple thresholds? Maybe sort them?
# I think its good to report if multiple thresholds were being detected, but not sure how to group them
# maybe simplest is to have another vector that records the number of thresholds detected each jackknife iteration,
# and then output the mean number detected per iteration for each simulation replicate
# which means I need a method for consistently selecting a threshold
# I would want the one closest to the mode, but would need to round to get a mode

jack_thresh <- function(simdt, xvals, thresh_methods = c("abs_max_d2", "min_d2", "zero_d2", "zero_d1"), sig_criteria, thresh_choice = "mean", knots = 4, smooth_type = "tp", span = 0.1){

  # get the number of simulations that were run
  nsim <- length(unique(simdt$sim))

  # get the true value of the threshold (for calculating rmse)
  thresh_loc <- simdt$thresh_loc[1]

  # turn the driver values into a dataframe
  xdt <- data.frame(
    driver = xvals
  )

  #for(mm in 1:length(thresh_methods)){

  #thresh_method <- thresh_methods[mm]

  # make the holding matrices to store mean, SEs, and RMSEs of the thresholds, number of jackknife
  # iterations that detected a threshold for each simulation replicate, and average number of
  # thresholds detected per jackknife iteration for each threshold calculation method
  thresh_mean <- matrix(NA, nrow = nsim, ncol = length(thresh_methods))# mean threshold estimate
  thresh_se <- matrix(NA, nrow = nsim, ncol = length(thresh_methods))# SE of threshold estimates
  thresh_rmse <- matrix(NA, nrow = nsim, ncol = length(thresh_methods))# RMSE of threshold estimates
  thresh_n <- matrix(NA, nrow = nsim, ncol = length(thresh_methods))# number of jackknife iterations that detected at least one threshold (goes into calculating the se)
  thresh_n_mean <- matrix(NA, nrow = nsim, ncol = length(thresh_methods))# average number of thresholds that were detected per jackknife iteration

  # make holding vectors for the means and sd's of the driver values for each simulation
 # x_means <- rep(NA, nsim)
 # x_sds <- rep(NA, nsim)

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

    # make the holding list for the threshold value(s) for each jackknife iteration
    #thresh_l <- vector(mode = "list", length = length(dd[,1]))

    # make temporary holding list for thresholds before significance criteria are applied (only needed if sig_criteria = TRUE)
    #thresh_temp_l <- vector(mode = "list", length = length(dd[,1]))

    # holding vector for total number of thresholds detected each iteration
    #thresh_tot <- rep(NA,length(dd[,1]))

    # make holding lists for results for the different threshold methods
    # threshold value(s) for each jackknife iteration
    thresh_l <- vector(mode = "list", length = length(thresh_methods))
    # thresholds before significance criteria are applied (only really needed if sig_criteria = TRUE)
    thresh_tot <- vector(mode = "list", length = length(thresh_methods))
    # total number of thresholds detected each iteration
    thresh_temp_l <- vector(mode = "list", length = length(thresh_methods))


    # make holding matrices for the derivatives
    D1 <- matrix(NA,length(dd[,1]),length(xvals)) # first derivative
    D2 <- matrix(NA,length(dd[,1]),length(xvals)) # second derivative

    # perform jackknife
    for(int in 1:length(dd[,1])){ # for each observation
      # remove this observation from the dataset
      jackd <- dd[-int, ]

      gami <- gam(obs_response~s(driver,k=knots,bs=smooth_type),data = jackd) # fit a gam to the bootd data set

      D1iunsm <- gratia::derivatives(gami, data = xdt, order = 1)$derivative # use the derivatives function from gratia to approximate 1st derivative for the values of the driver in xvals
      D2iunsm <- gratia::derivatives(gami, data = xdt, order = 2)$derivative # use the derivatives function from gratia to approximate 2nd derivative for the values of the driver in xvals

      # smooth these (without smoothing, the 2nd deriv seems to be super noisy)
      D1i <- predict(loess(D1iunsm ~ driver, data=data.frame(D1iunsm = D1iunsm, driver = xvals), span=span))
      D2i <- predict(loess(D2iunsm ~ driver, data=data.frame(D2iunsm = D2iunsm, driver = xvals), span=span))

      # store these
      D1[int, ] <- D1i
      D2[int, ] <- D2i


      # threshold calculations
      for(mm in 1:length(thresh_methods)){ # for each threshold method

        thresh_method <- thresh_methods[mm]

      if(thresh_method=="abs_max_d2"){
        # find peaks (max and mins) in the second deriv
        d2pks <- c(findPeaks(D2i), findPeaks(-D2i))

        if(sig_criteria==FALSE){

          thresh_l[[mm]][[int]] <- NA # holding place

        # if there are no significance criteria, pick the peak with max abs value and store the corresponding xvalue
        # note: need to store the x val of the threshold not the index, because later take the mean across jackknife iterations (and mean of indeces doesn't make sense)
        if(length(d2pks)>0){
        thresh_l[[mm]][[int]] <- xvals[d2pks[which(abs(D2i[d2pks])==max(abs(D2i[d2pks])))]]
        }

        # store the number of thresholds detected
        thresh_tot[[mm]][int] <- length(which(is.na(thresh_l[[mm]][[int]])==FALSE))

        } else { # otherwise store them all to process after jackknifing

          thresh_temp_l[[mm]][[int]] <- d2pks
        }


      }

      if(thresh_method=="min_d2"){
        # find the minima in the second deriv
        d2mins <- findPeaks(-D2i)

        if(sig_criteria==FALSE){

          thresh_l[[mm]][[int]] <- NA # holding place

          # if there are no significance criteria, pick the min with max abs value
          if(length(d2mins)>0){
          thresh_l[[mm]][[int]] <- xvals[d2mins[which(abs(D2i[d2mins])==max(abs(D2i[d2mins])))]]
          }

          # store the number of thresholds detected
          thresh_tot[[mm]][int] <- length(which(is.na(thresh_l[[mm]][[int]])==FALSE))

        } else { # otherwise store them all to process after jackknifing

          thresh_temp_l[[mm]][[int]] <- d2mins
        }

      }

      if(thresh_method=="zero_d2"){

        # find the points where the second deriv crosses zero (=changes sign)
        # sign(-10) = -1, sign(10) = 1, sign(0) = 0
        d2roots <- which(diff(sign(D2i))!=0)

        if(sig_criteria==FALSE){

          thresh_l[[mm]][[int]] <- NA # holding place

          # if there are no significance criteria, keep all roots and choose the values of x halfway between the sign change locations
          if(length(d2roots)>0){
          thresh_l[[mm]][[int]] <- (xvals[d2roots] + xvals[d2roots + 1])/2
          }

          # store the number of thresholds detected
          thresh_tot[[mm]][int] <- length(which(is.na(thresh_l[[mm]][[int]])==FALSE))

        } else { # otherwise store them all to process after jackknifing

          thresh_temp_l[[mm]][[int]] <- d2roots
        }

      }

      if(thresh_method=="zero_d1"){

        # find the points where the first deriv crosses zero (=changes sign)
        d1roots <- which(diff(sign(D1i))!=0)

        if(sig_criteria==FALSE){

          thresh_l[[mm]][[int]] <- NA # holding place

          # if there are no significance criteria, keep all roots and choose the values of x halfway between the sign change locations
          if(length(d1roots)>0){
          thresh_l[[mm]][[int]] <- (xvals[d1roots] + xvals[d1roots + 1])/2
          }

          # store the number of thresholds detected
          thresh_tot[[mm]][int] <- length(which(is.na(thresh_l[[mm]][[int]])==FALSE))

        } else { # otherwise store them all to process after jackknifing

          thresh_temp_l[[mm]][[int]] <- d1roots
        }

      }

     } # end of threshold calculations

      } # end of jackknife iterations



    # if significance criteria = true, get the CIs of the derivatives and update the thresh_l
    if(sig_criteria == TRUE){
      # get the 2.5% and 97.5% quantiles of the jackknifed first and second derivatives
      d2_low <- apply(D2,2,quantile,probs=0.025) # 2 in apply() means apply quantile() across columns of D2 matrix
      d2_up <- apply(D2,2,quantile,probs=0.975)

      d1_low <- apply(D1,2,quantile,probs=0.025)
      d1_up <- apply(D1,2,quantile,probs=0.975)

      # then apply the significant criteria to all the thresholds that were detected
      for(mm in 1:length(thresh_methods)){

        thresh_method <- thresh_methods[mm]

      if(thresh_method =="abs_max_d2"){

        for(int in 1:length(dd[,1])){ # for each of the jackknife iterations

          # locations where 2nd deriv is significantly different from 0
          signifD2 <- which(!data.table::between(0, d2_low, d2_up, incbounds=TRUE))

          sigd2pks <- intersect(signifD2, thresh_temp_l[[mm]][[int]])

          thresh_l[[mm]][[int]] <- NA # holding place

          if(length(sigd2pks)>0){
          # choose the detected threshold(s) for which the derivative is significantly different from zero
          thresh_l[[mm]][[int]] <- xvals[sigd2pks[which(abs(D2i[sigd2pks])==max(abs(D2i[sigd2pks])))]]
          }

          # store the number of thresholds detected
          thresh_tot[[mm]][int] <- length(which(is.na(thresh_l[[mm]][[int]])==FALSE))

        }

      }

      if(thresh_method =="min_d2"){

        for(int in 1:length(dd[,1])){ # for each of the jackknife iterations

          # locations where 2nd deriv is significantly different from 0
          signifD2 <- which(!data.table::between(0, d2_low, d2_up, incbounds=TRUE))

          sigd2mins <- intersect(signifD2, thresh_temp_l[[mm]][[int]])

          thresh_l[[mm]][[int]] <- NA # holding place

          if(length(sigd2mins)>0){
          # choose the detected threshold(s) for which the derivative is significantly different from zero
          thresh_l[[mm]][[int]] <- xvals[sigd2mins[which(abs(D2i[sigd2mins])==max(abs(D2i[sigd2mins])))]]
          }

          # store the number of thresholds detected
          thresh_tot[[mm]][int] <- length(which(is.na(thresh_l[[mm]][[int]])==FALSE))

        }

      }

      if(thresh_method =="zero_d2"){

        for(int in 1:length(dd[,1])){ # for each of the jackknife iterations

          thresh_l[[mm]][[int]] <- NA # holding place


          if(length(thresh_temp_l[[mm]][[int]])>0){
          # locations where 2nd deriv is NOT significantly different from 0
          nsignifD2 <- which(data.table::between(0, d2_low, d2_up, incbounds=TRUE))

          if(length(nsignifD2)>0){
          # select the roots where the deriv was not significantly different from zero at at least one of the two points where the sign change occurred
          # (e.g., if d2root = 3, it means d2[3] and d2[4] had different signs so 2nd deriv crossed zero between them, and say
          # this root is significant if either d2[3] or d2[4] (or both) are not significantly different from zero
          nsigd2rootsL <- intersect(nsignifD2, thresh_temp_l[[mm]][[int]]) # left of crossing is not sig diff from 0
          nsigd2rootsR <- intersect(nsignifD2, thresh_temp_l[[mm]][[int]] + 1) # right of crossing is not sig diff from 0
          nsigd2roots <- unique(nsigd2rootsL, nsigd2rootsR-1)

          if(length(nsigd2roots)>0){
          # get the sum of the differences between zero and the upper quantile, lower quantile, and int'th second deriv on both sides of the root
          summdiff <- abs(d2_low[nsigd2roots]-0) + abs(d2_up[nsigd2roots]-0) + abs(D2[int,][nsigd2roots]-0) + abs(d2_low[nsigd2roots+1]-0) + abs(d2_up[nsigd2roots+1]-0) + abs(D2[int,][nsigd2roots+1]-0)
          index<-nsigd2roots[which(summdiff==min(summdiff))] # get value where this sum is minimized

          # save the threshold value
          thresh_l[[mm]][[int]] <- (xvals[index] + xvals[index + 1])/2
          }
          }
          }

          # store the number of thresholds detected
          thresh_tot[[mm]][int] <- length(which(is.na(thresh_l[[mm]][[int]])==FALSE))

        }

      }


      if(thresh_method =="zero_d1"){

        for(int in 1:length(dd[,1])){ # for each of the jackknife iterations

          thresh_l[[mm]][[int]] <- NA # holding place

          if(length(thresh_temp_l[[mm]][[int]])>0){
            # locations where 1st deriv is NOT significantly different from 0
            nsignifD1 <- which(data.table::between(0, d1_low, d1_up, incbounds=TRUE))

            if(length(nsignifD1)>0){

              # select the roots where the deriv was not significantly different from zero at at least one of the two points where the sign change occurred
              nsigd1rootsL <- intersect(nsignifD1, thresh_temp_l[[mm]][[int]]) # left of crossing is not sig diff from 0
              nsigd1rootsR <- intersect(nsignifD1, thresh_temp_l[[mm]][[int]] + 1) # right of crossing is not sig diff from 0
              nsigd1roots <- unique(nsigd1rootsL, nsigd1rootsR-1)

              if(length(nsigd1roots)>0){
                # get the sum of the differences between zero and the upper quantile, lower quantile, and int'th second deriv on both sides of the root
                summdiff <- abs(d1_low[nsigd1roots]-0) + abs(d1_up[nsigd1roots]-0) + abs(D1[int,][nsigd1roots]-0) + abs(d1_low[nsigd1roots+1]-0) + abs(d1_up[nsigd1roots+1]-0) + abs(D1[int,][nsigd1roots+1]-0)
                index<-nsigd1roots[which(summdiff==min(summdiff))] # get value where this sum is minimized

                # save the threshold value
                thresh_l[[mm]][[int]] <- (xvals[index] + xvals[index + 1])/2

              }
            }
          }

          # store the number of thresholds detected
          thresh_tot[[mm]][int] <- length(which(is.na(thresh_l[[mm]][[int]])==FALSE))

        }

      }

     } # end of threshold calculations


    } # end of if(sig_criteria = TRUE)


    # final threshold calculations
    for(mm in 1:length(thresh_methods)){

      thresh_method <- thresh_methods[mm]
    # if more than one threshold was detected (e.g., if sig_criteria = FALSE and there were multiple roots), choose only one of them based on thresh_choice
    if(max(thresh_tot[[mm]], na.rm = T) > 1){ # if at least one of the jackknife iterations found more than one threshold

      if(thresh_choice=="mean"){
        # get the mean of all thresholds detected
        mean_all <- mean(unlist(thresh_l[[mm]]), na.rm = T)

        thresh_pick <- rep(NA, length(dd[,1])) # holding vector for selected thresholds

        for(ii in 1:length(dd[,1])){ # for each jackknife iteration

         # if more than one threshold was detected in the ith iteration, pick the one closest to the mean of all thresholds detected across all simulations
          thresh_pick[ii] <- ifelse(length(thresh_l[[mm]][[ii]])>1, thresh_l[[mm]][[ii]][which(abs(thresh_l[[mm]][[ii]]-mean_all)==min(abs(thresh_l[[mm]][[ii]]-mean_all)))], thresh_l[[mm]][[ii]][1])

        }

      }

      if(thresh_choice=="median"){
        # get the mean of all thresholds detected
        med_all <- median(unlist(thresh_l[[mm]]), na.rm = T)

        thresh_pick <- rep(NA, length(dd[,1])) # holding vector for selected thresholds

        for(ii in 1:length(dd[,1])){ # for each jackknife iteration

          # if more than one threshold was detected in the ith iteration, pick the one closest to the median of all thresholds detected across all simulations (if 2 are equidistant, just pick the first one)
          thresh_pick[ii] <- ifelse(length(thresh_l[[mm]][[ii]])>1, thresh_l[[mm]][[ii]][which(abs(thresh_l[[mm]][[ii]]-med_all)==min(abs(thresh_l[[mm]][[ii]]-med_all)))][1], thresh_l[[mm]][[ii]][1])

        }

      }

      thresh <- thresh_pick

    } else { # if every iteration detected at most one threshold for mm'th method
      thresh <- unlist(thresh_l[[mm]]) # just unlist thresh_l to get the vector of thresholds
    }

    # get the mean and SEs of the chosen thresholds for the ith simulation
    thresh_mean[i,mm] <- mean(thresh, na.rm = T) # mean of thresholds detected across the jackknife iterations
    thresh_n[i,mm] <- length(which(thresh_tot[[mm]]>0))# number of jackknife iterations that detected at least one threshold
    thresh_real <- thresh[which(is.na(thresh)==FALSE)] # thresholds that weren't NA
    thresh_se[i,mm] <- ifelse(length(thresh_real) > 2, sd(thresh_real)/sqrt(length(thresh_real)-1), NA)
    thresh_rmse[i, mm] <- ifelse(length(thresh_real) > 1, sqrt(sum((thresh_real-thresh_loc)^2)/length(thresh_real)), NA)

    # get the mean number of thresholds that were detected each iteration
    thresh_n_mean[i,mm] <- mean(thresh_tot[[mm]], na.rm = T) # average number of thresholds that were detected each jackknife iteration



  } # end of calculations for each threshold method



  } # end of iterations for each simulation i


  df <- data.frame(
    sim = rep(c(1:nsim), length(thresh_methods)),
    thresh_method = rep(thresh_methods, each = nsim),
    thresh_mean = c(thresh_mean),
    thresh_n = c(thresh_n),
    thresh_se = c(thresh_se),
    thresh_rmse = c(thresh_rmse),
    thresh_n_mean = c(thresh_n_mean),
    #x_mean = rep(x_means, length(thresh_methods)),
    #x_sd = rep(x_sds, length(thresh_methods)),
    thresh_diff = c(thresh_mean-thresh_loc)#,
    #thresh_diffN = c((thresh_mean-thresh_loc)/ rep(x_sds, length(thresh_methods)))
  )

  #df <- data.frame(
    #sim = c(1: nsim),
    #thresh_mean = thresh_mean,
    #thresh_n = thresh_n,
    #thresh_se = thresh_se,
    #thresh_n_mean = thresh_n_mean,
    #thresh_method = thresh_method
  #)

  #}

  # return results
  return(df)
  # return(list(thresh_mean = thresh_mean, thresh_n = thresh_n, thresh_se = thresh_se, thresh_n_mean = thresh_n_mean))

}
