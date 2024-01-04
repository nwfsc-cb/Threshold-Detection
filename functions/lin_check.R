#' function to check whether there is support for a nonlinear model in the simulation data
#' required packages:
#' mgcv (for fitting the gam)
#' gratia (for derivative calculations),
#' quantmod (for the findPeaks function to find max/mins in the derivatives),
#' data.table (for between function to determine whether CIs contain 0)


#' function inputs:
#' @param simdt data frame with the simulated data (output of sim_data function)
#' @param thresh_methods vector of threshold definitions to use for estimating the threshold location. (this is just
#' to later join the df produced by this function w/ the df output by the jack_thresh function)
#' @param sig_criteria vector of sig criteria to use for evaluating threshold significance (this is just
#' to later join the df produced by this function w/ the df output by the jack_thresh function)
#' @param z_scale if TRUE, standardize the data before fitting (default), if FALSE use raw data
#' @param knots number of knots to use in the gam fitting, defaults to 4
#' @param smooth_type type of smooth used with gams, defaults to "tp"


lin_check <- function(simdt, thresh_methods = c("abs_max_d2", "min_d2", "zero_d2", "zero_d1"),
                      sig_criteria = c("none", "jack_quant", "sim_int"), z_scale = TRUE, knots = 4, smooth_type = "tp", span = 0.1){

  # get the number of simulations that were run
  nsim <- length(unique(simdt$sim))

  # holding vector for which model is best (lowest AIC)
  best_mod <- rep(NA, nsim)

  # holding vector for differences in AIC between the best and next-best model
  aic_diff <- rep(NA, nsim)

  for(i in 1:nsim){ # for each simulation replicate

    # get the data
    datIN <- simdt[which(simdt$sim == i), ]

    if(z_scale==TRUE){ # standardize it
      datIN$driver <- zfun(datIN$driver)
      datIN$obs_response <- zfun(datIN$obs_response)
    }

    # fit the gam
    #if(z_scale==TRUE){
    #gam_full <- gam(obs_response~ -1 + s(driver,k=knots,bs=smooth_type),data = datIN) # fit a gam to the data set
    #lm_full <- lm(obs_response ~ 0 + driver, data = datIN)

    #} else{
     # gam_full <- gam(obs_response~s(driver,k=knots,bs=smooth_type),data = datIN) # fit a gam to the jackd data set
     # lm_full <- lm(obs_response ~ driver, data = datIN)
   # }

    # fit the GAM
    gam_full <- gam(obs_response~s(driver,k=knots,bs=smooth_type),data = datIN) # fit a gam to the data set

    # fit the linear model
    lm_full <- lm(obs_response ~ driver, data = datIN)

    # get the AIC values
    gam_aic <- AIC(gam_full)
    lm_aic <- AIC(lm_full)

    if(gam_aic < lm_aic){ # if gam has smaller AIC than the lm

      if(abs(lm_aic - gam_aic) >=2){ # if difference is at least 2
        best_mod[i] <- "gam" # best model is gam
      } else{
        best_mod[i] <- "gam_ns" # best model is gam but not significant
      }

    } else { # if lm has AIC that is less than or equal to gam AIC

      if(abs(lm_aic - gam_aic) >=2){ # if difference is at least 2
        best_mod[i] <- "lm" # best model is lm
      } else{
        best_mod[i] <- "lm_ns" # best model is lm but not significant
      }

    }

    # record the difference in AIC values
    aic_diff[i] <- lm_aic - gam_aic # values > +2 mean gam is better, values < -2 mean lm is better

  }

  # save the results as a dataframe


  #df <- data.frame(
   # sim = rep(1:nsim, length(thresh_methods)),
   # best_mod = rep(best_mod, length(thresh_methods)),
   # thresh_method = rep(thresh_methods, each = nsim)
  #)

  thresh_methods_ss <- rep(thresh_methods, each = length(sig_criteria))
  nsim_ss <- rep(c(1:nsim), each = length(sig_criteria))
  best_mod_ss <- rep(best_mod, each = length(sig_criteria)) # best model and AIC differences are the same for all the threshold definitions and significance criteria
  aic_diff_ss <- rep(aic_diff, each = length(sig_criteria))


  df <- data.frame(
    sim = rep(nsim_ss, length(thresh_methods)),
    best_mod = rep(best_mod_ss, length(thresh_methods)),
    aic_diff = rep(aic_diff_ss, length(thresh_methods)),
    thresh_method = rep(thresh_methods_ss, each = nsim),
    sig_type = rep(sig_criteria, nsim*length(thresh_methods))
  )


  return(df)

}


