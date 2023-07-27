
#' function for bootstrapping the gam predictions and their first and second derivatives
#' requires the deriv2 function for calculating derivatives

#' adapted from the threshold function in:
#' https://github.com/kholsman/EBM_Holsman_NatComm/blob/master/R/sub_fun/threshold.R
#' From Holsman etal 2020
#' Holsman, K. K., A. Haynie, A. Hollowed, J.C.P. Reum, K. Aydin, A. J. Hermann, W. Cheng,
#' A. Faig, J. N. Ianelli, K. Kearney, A. E. Punt, 2020.Ecosystem-based fisheries
#;  management forestalls climate-driven collapse, Nature Communications.

#' function inputs:
#' nsim.f: number of simulations run
#' simdt.f: output from simfun1 (outputs from the simulations)
#' knots.f: number of knots to use in the gam fitting
#' boot_n.f: number of bootstrap interations to run
#' boot_nobs.f: number of samples to take each bootstrap run
#' xvals.f: vector of driver values (xvals) to use for generating predictions from the gams

bootfun <- function(nsim.f, simdt.f, knots.f, boot_n.f, boot_nobs.f, xvals.f){

  # make the holding lists to store results for each simulation in simdt
  Deriv1s <- vector(mode = "list", length =nsim.f) # first derivative
  Deriv2s <- vector(mode = "list", length =nsim.f) # second derivative
  hatFits <- vector(mode = "list", length =nsim.f) # predicted values
  hatses <- vector(mode = "list", length =nsim.f) # se of fitted values
  gmlists <- vector(mode = "list", length =nsim.f) # gams fit in each bootstrap run (this will be a nested list, where outer list corresponds to each simulation and inner list to each bootstrap run for that simulation)

  for(i in 1:nsim.f){

    # subset the data for the ith simulation
    datIN <- simdt.f[which(simdt.f$sim == i), ]

    # make data set for the bootstrapping
    dd <- datIN%>%mutate(driver = round(driver,3))%>%select(driver, obs_response)
    dd$num <- 1:length(dd[,1])

    # make the holding matrices
    Deriv1 <- matrix(NA,boot_n.f,length(xvals.f)) # first derivative
    Deriv2 <- matrix(NA,boot_n.f,length(xvals.f)) # second derivative
    hatFit <- matrix(NA,boot_n.f,length(xvals.f)) # predicted values
    hatse  <- matrix(NA,boot_n.f,length(xvals.f)) # se of fitted values
    gmlist <- vector(mode = "list", length =boot_n.f) # empty list to store the gam model fit to the data subset each bootstrap run

    # Run the boot strap

    for(int in 1:boot_n.f){ # for each bootstrap interation
      # get bootstrapped sub-sample
      nobs <- length(dd$num) # number of observations in the data set
      if(boot_nobs.f > nobs) # if number of observations to sample for bootstrap is greater than nobs
        boot_nobs.f   <- nobs # set number of bootstrap samples equal to nobs

      bootd <- sample_n(dd,boot_nobs.f,replace = TRUE) # sample boot_nobs.f rows (with replacement) from the dd data set and name this subsetted data frame bootd
      tmpgam <- gam(obs_response~s(driver,k=knots.f,bs="tp"),data = bootd) # fit a gam to the bootd data set
      tmpd <- deriv2(tmpgam,xvals.f) # use the deriv2 function to approximate the first and second derivatives for the values of the driver in xvals.f
      gmlist[[int]] <- tmpgam # store the gam model for the int data subset in gmlist
      Deriv1[int,]  <- tmpd$fd_d1 # store the values of the first derivative across the values of the driver in xvals.f for the gam fit to the int data subset
      Deriv2[int,]  <- tmpd$fd_d2 # store the values of the second derivative across the values of the driver in xvals.f for the gam fit to the int data subset
      hatFit[int,]  <- predict(tmpgam,se.fit=TRUE,newdata=data.frame(driver=xvals.f))$fit # store the predictions for each value of the driver in xvals.f for the gam fit to the int data subset
      hatse[int,]   <- predict(tmpgam,se.fit=TRUE,newdata=data.frame(driver=xvals.f))$se # store the se of the predictions

    }

    # store results
    Deriv1s[[i]] <- Deriv1
    Deriv2s[[i]] <- Deriv2
    hatFits[[i]] <- hatFit
    hatses[[i]] <- hatse
    gmlists[[i]] <- gmlist

  }



  # return results
  return(list(Deriv1s = Deriv1s, Deriv2s = Deriv2s, hatFits = hatFits, hatses = hatses, gmlists = gmlists))

}
