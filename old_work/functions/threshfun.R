

# function for calculating threshold values based on results of gam fitting and bootstrapping
# of simulated data

# adapted from the threshold function in:
# https://github.com/kholsman/EBM_Holsman_NatComm/blob/master/R/sub_fun/threshold.R
# From Holsman etal 2020
# Holsman, K. K., A. Haynie, A. Hollowed, J.C.P. Reum, K. Aydin, A. J. Hermann, W. Cheng,
# A. Faig, J. N. Ianelli, K. Kearney, A. E. Punt, 2020.Ecosystem-based fisheries
#  management forestalls climate-driven collapse, Nature Communications.


# inputs
# nsim.f: number of simulations used in simfun1()
# gamfunOUT.f: output from gamfun()
# bootfunOUT.f: output from bootfun()
# xvals.f: sequence of driver values used in gamfun() and bootfun() for generating predictions
#sdmult.f: involved in calculating threshold significance in method 2 and method 3 (default from Holsman et al. 2020 is 1)
#method.f: method for getting the threshold significance (Holsman et al. 2020 used method 2)
#prob.f: probabilities for the quantile ranges
#span.f: smoothing step/span for thresholds (Holsman et al. 2020 used 0.1 as default)

# NOTE this function will return at most one threshold per threshold calculation method,
# so will need to adjust calculations to allow for the possibility of there being multiple thresholds

threshfun <- function(nsim.f, gamfunOUT.f, bootfunOUT.f, xvals.f, sdmult.f, method.f, prob.f, span.f){

  # make the holding objects (indices of threshold estimates, actual estimated threshold values, as well as hat_qnt, df1_qnt, df2_qnt data frames for each simulation)
  thresh_ind1 <- vector(mode = "list", length = nsim.f) # index for thresh based on max/min in first deriv
  thresh_ind2 <- vector(mode = "list", length = nsim.f) # index for thresh based on second deriv being furthest from zero

  thresh_ind20 <- vector(mode = "list", length = nsim.f) # index for thresh based on second deriv being closest to zero


  thresh_1 <- vector(mode = "list", length = nsim.f) # value of driver for thresh based on max/min in first deriv
  thresh_2 <- vector(mode = "list", length = nsim.f) # value of driver for thresh based on second deriv being furthest from zero

  thresh_20 <- vector(mode = "list", length = nsim.f) # value of driver for thresh based onsecond deriv being closest to zero

  hat_qnts <- vector(mode = "list", length = nsim.f) # data frame with values of driver, mean, and upper and lower quantiles for predictions of gams fit to data from each simulation

  df1_qnts <- vector(mode = "list", length = nsim.f) # data frame with values of driver, mean, and upper and lower quantiles for first derivative of predictions of gams fit to data from each simulation

  df2_qnts <- vector(mode = "list", length = nsim.f) # data frame with values of driver, mean, and upper and lower quantiles for second derivative of predictions of gams fit to data from each simulation

  for(i in 1:nsim.f){ # for each simulation

    # get the outputs from the gamfun() and bootfun() functions for the ith simulation:
    hat <- gamfunOUT.f$hats[[i]] # predictions from the gam model fit to data for ith simulation

    Deriv1 <- bootfunOUT.f$Deriv1s[[i]] # bootstrapped first deriv
    Deriv2 <- bootfunOUT.f$Deriv2s[[i]] # bootstrapped second deriv
    hatFit <- bootfunOUT.f$hatFits[[i]] # bootstrapped predictions
    hatse <- bootfunOUT.f$hatses[[i]] # bootstrapped se of predictions

    # apply quantiles to bootstrap replicates

    # apply the quantile() function to each column of Deriv1, Deriv2, hatFit, and hatse (with the quantiles specified by prob.f)
    D1_se  <- apply(Deriv1,2,quantile,probs=prob.f)
    D2_se  <- apply(Deriv2,2,quantile,probs=prob.f)
    qnt    <- apply(hatFit,2,quantile,probs=prob.f)
    qntse  <- apply(hatse,2,quantile,probs=prob.f)

    # get the gam quantiles using one of 3 methods for calculating confidence intervals

    if(method.f==1)
      hat_qnt<-data.frame(driver=xvals.f,
                          up=hat$fit+qnt[3,]-qnt[2,],
                          mn=hat$fit,
                          dwn=hat$fit+qnt[1,]-qnt[2,])
    if(method.f==2)
      hat_qnt<-data.frame(driver=xvals.f,
                          up=hat$fit+sdmult.f*qntse[2,],
                          mn=hat$fit,
                          dwn=hat$fit-sdmult.f*qntse[2,])

    if(method.f==3)
      hat_qnt<-data.frame(driver=xvals.f,
                          up=hat$fit+sdmult.f*hat$se,
                          mn=hat$fit,
                          dwn=hat$fit-sdmult.f*hat$se)

    # smooth the means and upper/lower quantiles using loess, with the span parameter (controlling degree of smoothing) set to span.f
    hat_qnt$smoothed_mn  <- predict(loess(mn ~ driver, data=hat_qnt, span=span.f))
    hat_qnt$smoothed_dwn <- predict(loess(dwn ~ driver, data=hat_qnt, span=span.f))
    hat_qnt$smoothed_up  <- predict(loess(up ~ driver, data=hat_qnt, span=span.f))

    # repeat for first derivative quantiles
    df1_qnt<-data.frame(driver =xvals.f,
                        up  = D1_se[3,],
                        mn  = D1_se[2,],
                        dwn = D1_se[1,])

    # smooth these
    df1_qnt$smoothed_mn  <- predict(loess(mn  ~ driver, data=df1_qnt, span=span.f))
    df1_qnt$smoothed_dwn <- predict(loess(dwn ~ driver, data=df1_qnt, span=span.f))
    df1_qnt$smoothed_up  <- predict(loess(up  ~ driver, data=df1_qnt, span=span.f))


    # repepat for second derivative quantiles
    df2_qnt<-data.frame(driver = xvals.f,
                        up  = D2_se[3,],
                        mn  = D2_se[2,],
                        dwn = D2_se[1,])

    # smooth these
    df2_qnt$smoothed_mn  <- predict(loess(mn ~ driver, data=df2_qnt, span=span.f))
    df2_qnt$smoothed_dwn <- predict(loess(dwn ~ driver, data=df2_qnt, span=span.f))
    df2_qnt$smoothed_up  <- predict(loess(up ~ driver, data=df2_qnt, span=span.f))

    # determine peaks and valleys
    #findPeaks() is a function from quantmod package that finds the peaks and valleys in a series
    # findPeaks(y) gives index of y corresponding to max, findPeaks(-y) gives index of y corresponding to min; if there isn't a max or min it returns 0
    # first derivative
    pks1    <- sort(c(findPeaks(df1_qnt$smoothed_mn),findPeaks(-df1_qnt$smoothed_mn)))
    #pks1_up  <- sort(c(findPeaks(df1_qnt$smoothed_up),findPeaks(-df1_qnt$smoothed_up)))
    #pks1_dwn <- sort(c(findPeaks(df1_qnt$smoothed_dwn),findPeaks(-df1_qnt$smoothed_dwn)))

    # get region (i.e., rows of df1_qnt) where the first derivative is significantly different from 0
    # between(x, xmin, xmax) is a function from the data.table package that asks if x is between xmin and xmax (inclusive if incbounds = TRUE), if so returns TRUE otherwise returns FALSE
    # the ! means only want values where between returns FALSE (i.e, 0 is not between the lower and upper confidence intervals)
    signif1 <- which(!data.table::between(0, df1_qnt$dwn, df1_qnt$up, incbounds=TRUE))

    # threshold based on first derivative is the intersection of the rows where the first deriv is significantly different from 0 and it is at a max or min
    thrsh1_index_all  <- intersect(signif1,pks1)

    # NOTE: I don't think the first derivative necessarily has to be significantly different from zero, since it could equal zero at its max or min...


    #thrsh1_index <- thrsh1_index_all

    thrsh1_index <- which(1==10) # make empty integer

    # if there is at least one value in thrsh1_index_all, say thresh1 is the value for which the first deriv is largest in magnitude
    if(length(thrsh1_index_all)>0)
      thrsh1_index<-thrsh1_index_all[which(abs(df1_qnt$smoothed_mn[thrsh1_index_all])==max(abs(df1_qnt$smoothed_mn[thrsh1_index_all])))]

    #thrsh1 <- df1_qnt$driver[thrsh1_index] # value(s) of the driver at the threshold(s) defined by first derivative


    # second derivative, most different from zero approach (used in Holsman et al. 2020)
    pks2     <- sort(c(findPeaks(df2_qnt$smoothed_mn),findPeaks(-df2_qnt$smoothed_mn)))
    #pks2_up  <- sort(c(findPeaks(df2_qnt$smoothed_up),findPeaks(-df2_qnt$smoothed_up)))
    #pks2_dwn <- sort(c(findPeaks(df2_qnt$smoothed_dwn),findPeaks(-df2_qnt$smoothed_dwn)))



    # get region (rows of df2_qnt) where the second derivative is significantly different from zero
    signif2 <- which(!between(0, df2_qnt$dwn, df2_qnt$up, incbounds=TRUE))

    # get the region for which the second deriv is furthest from zero and significantly different from zero
    thrsh2_index_all <- intersect(signif2,pks2)

    #thrsh2 <- df2_qnt$driver[thrsh2_index_all] # value(s) of the driver at the threshold(s) defined by second derivative being furthest from zero

    thrsh2_index <- which(1==10) # make empty integer

    # if there is at least one value in thrsh2_index_all, say thresh2 is the value for which the second deriv is largest in magnitude
    if(length(thrsh2_index_all)>0) # if there is at least one value in thrsh2_index_all
      thrsh2_index<-thrsh2_index_all[which(abs(df2_qnt$smoothed_mn[thrsh2_index_all])==max(abs(df2_qnt$smoothed_mn[thrsh2_index_all])))] # get which one of these corresponds to the value for which the smoothed mean of the second deriv is at its max value. Original code from the thresholds workshop (below) took the mean of these, but taking the mean of an index didn't seem to make sense
    #thrsh2_index<-mean(thrsh2_index_all[which(abs(df2_qnt$smoothed_mn[thrsh2_index_all])==max(abs(df2_qnt$smoothed_mn[thrsh2_index_all])))],na.rm=T)


    # add columns to df1_qnt and df2_qnt indicating whether first or second deriv is significantly different from zero, respectively
    # add new column "sig" to df1_qnt, and df2_qnt that is just NA
    df1_qnt$sig <- NA
    df2_qnt$sig <-NA
    #hat_qnt <- NA

    # if upper and lower CIs for second deriv don't include 0, sig is TRUE, otherwise it is FALSE
    df2_qnt$sig <- !between(0, df2_qnt$dwn, df2_qnt$up, incbounds=TRUE)
    # if sig is FALSE, replace it with NA
    df2_qnt$sig[!df2_qnt$sig] <-NA

    # do same thing for first derivative
    df1_qnt$sig <- !between(0, df1_qnt$dwn, df1_qnt$up, incbounds=TRUE)
    df1_qnt$sig[!df1_qnt$sig] <-NA

    # second derivative method 2: second deriv quantiles include and are closest to zero (since s''(x)=0 at inflection point)
    # get the region (rows) where the second deriv is NOT significantly different from zero
    nsignif20 <- which(between(0, df2_qnt$dwn, df2_qnt$up, incbounds=TRUE))

    # get the location where the confidence intervals are closest to zero
    thrsh20_index <- which(1==10) # make empty integer

    if(length(nsignif20)>0) # if there is at least one value in nsignif20
      # get the sum of the differences between zero and the upper quantile, lower quantile, and mean of the second deriv
      summdiff <- abs(df2_qnt$smoothed_up[nsignif20]-0) + abs(df2_qnt$smoothed_dwn[nsignif20]-0) + abs(df2_qnt$smoothed_mn[nsignif20]-0)
    thrsh20_index<-nsignif20[which(summdiff==min(summdiff))] # thrsh20_index is the value where this sum is minimized

    #store the results
    thresh_ind1[[i]] <- thrsh1_index # index for thresh based on max/min in first deriv for gam fit to data from ith simulation
    thresh_ind2[[i]] <- thrsh2_index# index for thresh based on second deriv being furthest from zero for gam fit to data from ith simulation

    thresh_ind20[[i]] <- thrsh20_index # index for thresh based on second deriv being close to zero

    thresh_1[[i]] <- df1_qnt$driver[thrsh1_index] # value(s) of the driver at the threshold(s) defined by max/min of first derivative for the ith simulation

    thresh_2[[i]] <- df2_qnt$driver[thrsh2_index]# value(s) of the driver at the threshold(s) defined by max/min of second derivative for the ith simulation

    thresh_20[[i]] <- df2_qnt$driver[thrsh20_index] # value(s) of the driver at the threshold(s) defined by second derivative close to 0 for the ith simulation

    hat_qnts[[i]] <- hat_qnt # data frame with values of driver, mean, and upper and lower quantiles for predictions of gams fit to data from ith simulation

    df1_qnts[[i]] <- df1_qnt # data frame with values of driver, mean, and upper and lower quantiles for first derivative of predictions of gams fit to data from ith simulation

    df2_qnts[[i]] <- df2_qnt # data frame with values of driver, mean, and upper and lower quantiles for second derivative of predictions of gams fit to data from ith simulation

  }

  # return the results
  return(list(thresh_ind1 = thresh_ind1, thresh_ind2 = thresh_ind2, thresh_ind20 = thresh_ind20, thresh_1 = thresh_1, thresh_2 = thresh_2, thresh_20 = thresh_20, hat_qnts = hat_qnts, df1_qnts = df1_qnts, df2_qnts = df2_qnts))

}


