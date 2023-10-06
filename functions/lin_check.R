#' function to check whether there is support for a nonlinear model in the simulation data
#' required packages:
#' mgcv (for fitting the gam)
#' gratia (for derivative calculations),
#' quantmod (for the findPeaks function to find max/mins in the derivatives),
#' data.table (for between function to determine whether CIs contain 0)


#' function inputs:
#' @param simdt data frame with the simulated data (output of simfun2)
#' @param thresh_methods vector of methods to use for estimating the threshold location. (this is just
#' to later join the df produced by this function w/ the df output by the jack_thresh function)
#' @param sig_criteria vector of sig criteria to use for evaluating threshol significance (this is just
#' to later join the df produced by this function w/ the df output by the jack_thresh function)
#' @param knots number of knots to use in the gam fitting, defaults to 4
#' @param smooth_type type of smooth used with gams, defaults to "tp"


lin_check <- function(simdt, thresh_methods = c("abs_max_d2", "min_d2", "zero_d2", "zero_d1"),
                      sig_criteria = c("none", "jack_quant", "sim_int"), knots = 4, smooth_type = "tp", span = 0.1){

  # get the number of simulations that were run
  nsim <- length(unique(simdt$sim))

  # holding vector for which model is best (lowest AIC)
  best_mod <- rep(NA, nsim)

  for(i in 1:nsim){ # for each simulation replicate

    # get the data
    datIN <- simdt[which(simdt$sim == i), ]

    # fit the gam
    gam_full <- gam(obs_response~s(driver,k=knots,bs=smooth_type),data = datIN) # fit a gam to the jackd data set

    # fit the linear model
    lm_full <- lm(obs_response ~ driver, data = datIN)

    # get the AIC values
    gam_aic <- AIC(gam_full)
    lm_aic <- AIC(lm_full)

    if(gam_aic < lm_aic){
      best_mod[i] <- "gam"
    } else {
      best_mod[i] <- "lm"
    }


  }

  # save the results as a dataframe


  #df <- data.frame(
   # sim = rep(1:nsim, length(thresh_methods)),
   # best_mod = rep(best_mod, length(thresh_methods)),
   # thresh_method = rep(thresh_methods, each = nsim)
  #)

  thresh_methods_ss <- rep(thresh_methods, each = length(sig_criteria))
  nsim_ss <- rep(c(1:nsim), each = length(sig_criteria))

  df <- data.frame(
    sim = rep(nsim_ss, length(thresh_methods)),
    best_mod = rep(best_mod, length(thresh_methods)*length(sig_criteria)),
    thresh_method = rep(thresh_methods_ss, each = nsim),
    sig_type = rep(sig_criteria, nsim*length(thresh_methods))
  )


  return(df)

}


