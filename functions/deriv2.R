#' Get first and second derivatives based on Large et al. 2014
#' copied from:
#' https://github.com/kholsman/EBM_Holsman_NatComm/blob/master/R/sub_fun/threshold.R
#' From Holsman etal 2020
#' Holsman, K. K., A. Haynie, A. Hollowed, J.C.P. Reum, K. Aydin, A. J. Hermann, W. Cheng,
#' A. Faig, J. N. Ianelli, K. Kearney, A. E. Punt, 2020.Ecosystem-based fisheries
#'  management forestalls climate-driven collapse, Nature Communications.

#' inputs:
#' gam_mod = the gam model
#' simdata= values of the driver used for generating predictions

# function inputs: gam fit to simulated data and vector of driver values (xvals) to use for generating predictions from the gam
deriv2 <- function(gam_mod,xvals.f){
  # finite difference approach to derivatives following
  # example from ?predict.gam
  DF1  <- gam_mod$model
  eps  <- 1e-7
  # new data for prediction
  newDF <- with(DF1, data.frame(driver = xvals.f))

  # prediction of smoothed estimates at each unique year value
  # with standard error
  B <- predict.gam(gam_mod,  newDF, type="response", se.fit=TRUE)

  # lines(scale(newDF$S_driver),scale(B$fit),pch=16,lwd=2)
  X0 <- predict(gam_mod, newDF, type = 'lpmatrix')
  newDFeps_p <- newDF + eps
  X1 <- predict(gam_mod, newDFeps_p, type = 'lpmatrix')
  X1sim <- simulate(gam_mod,nsim=200,newdata=newDFeps_p)

  # finite difference approximation of first derivative
  # the design matrix
  Xp <- (X0 - X1) / eps

  # first derivative
  fd_d1  <- -1*Xp %*% coef(gam_mod)
  D1_dat <- data.frame(driver=newDF$driver,y=fd_d1)[order(newDF$driver),]

  # plot(D1_dat[,1],D1_dat[,2],type="l",lwd=2, xlim=c(-3,3),ylim=c(-2,2))

  # second derivative
  newDFeps_m <- newDF - eps
  X_1  <- predict(gam_mod, newDFeps_m, type = 'lpmatrix')

  # design matrix for second derivative
  Xpp <- (X1 + X_1 - 2*X0)  / eps^2
  # second derivative
  fd_d2 <- Xpp %*% coef(gam_mod)

  D2_dat <- data.frame(driver=newDF$driver,y=fd_d2)[order(newDF$driver),]
  # plot(D2_dat[,1],D2_dat[,2],type="l",lwd=2, xlim=c(-3,3),ylim=c(-2,2))

  return(list(newDF=newDF,fd_d1=as.numeric(fd_d1),fd_d2=as.numeric(fd_d2)))#Returns list with 1) newDF (new data used for prediction), 2) fd_d1 (vector of first deriv values at each value of the driver in x), and 3) fd_d2 (vector of second deriv values at each value of the driver in x)
}
