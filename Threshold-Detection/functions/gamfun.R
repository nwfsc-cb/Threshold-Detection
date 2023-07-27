
#' function for fitting gams to simulated data


#' #' adapted from the threshold function in:
#' https://github.com/kholsman/EBM_Holsman_NatComm/blob/master/R/sub_fun/threshold.R
#' From Holsman etal 2020
#' Holsman, K. K., A. Haynie, A. Hollowed, J.C.P. Reum, K. Aydin, A. J. Hermann, W. Cheng,
#' A. Faig, J. N. Ianelli, K. Kearney, A. E. Punt, 2020.Ecosystem-based fisheries
#;  management forestalls climate-driven collapse, Nature Communications.

#' function inputs:
#' nsim.f: number of simulations run
#' simdt.f: output from simfun1 (outputs from the simulations)
#' knots.f: number of knots to use in the gam fitting
#' xvals.f: vector of driver values (xvals) to use for generating predictions from the gams


gamfun <- function(nsim.f, simdt.f, knots.f, xvals.f){

  # holding lists for the fitted gams
  gams <- vector(mode = "list", length = nsim.f)

  # holding list for the predicted values
  hats <- vector(mode = "list", length = nsim.f)

  for(i in 1:nsim.f){
    # subset the data for the ith simulation
    datIN <- simdt.f[which(simdt.f$sim == i), ]

    # fit a gam to the data (observed response as a function of driver)
    gam_i <- gam(obs_response ~ s(driver,k=knots.f,bs="tp"),data = datIN) # s() defines smoother in the GAM: takes list of covariate variables that the smoother will be a function of (here, "driver"); k = number of knots (wiggliness of smoother); bs indicates the (penalized) smoothing basic to use ("tp" means thin plate regression spline)

    # save gam_i as the ith element of the gams list
    gams[[i]] <- gam_i

    # get the predictions from the fitted gam
    hat_i <- predict(gam_i,se.fit=TRUE, newdata = data.frame(driver=xvals.f))

    # save hat_i
    hats[[i]] <- hat_i


  }

  # return list of the gam model for each simulation in simdt and a list of the hats for each

  return(list(gams = gams, hats = hats))

}
