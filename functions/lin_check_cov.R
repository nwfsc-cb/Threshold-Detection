#' function to check which model has the best support (linear or gam), with the covariate included
#' in candidate models
#' required packages:
#' mgcv (for fitting the gam)
#' gratia (for derivative calculations),



#' function inputs:
#' @param simdt data frame with the simulated data (output of sim_data)
#' @param thresh_methods vector of methods to use for estimating the threshold location. (this is just
#' to later join the df produced by this function w/ the df output by the jack_thresh function)
#' @param sig_criteria vector of sig criteria to use for evaluating threshol significance (this is just
#' to later join the df produced by this function w/ the df output by the jack_thresh function)
#' @param z_scale if TRUE, standardize the data before fitting (default), if FALSE use raw data
#' @param knots number of knots to use in the gam fitting, defaults to 4
#' @param smooth_type type of smooth used with gams, defaults to "tp"
#' @param aicc_size sample size below which AICc is used instead of AICc


lin_check_cov <- function(simdt, thresh_methods = c("abs_max_d2", "min_d2", "zero_d2", "zero_d1"),
                      sig_criteria = c("none", "jack_quant", "sim_int"), z_scale = TRUE, knots = 4,
                      smooth_type = "tp", aicc_size = 40){ #span = 0.1

  # get the number of simulations that were run
  nsim <- length(unique(simdt$sim))

  # holding vector for which model is best (lowest AIC)
  best_mod <- rep(NA, nsim)

  # holding vector for next best model
  next_mod <- rep(NA, nsim)

  # holding vector for differences in AIC
  aic_diff <- rep(NA, nsim)

  for(i in 1:nsim){ # for each simulation replicate

    # get the data
    datIN <- simdt[which(simdt$sim == i), ]

    if(z_scale==TRUE){
      datIN$driver <- zfun(datIN$driver)
      datIN$cov <- zfun(datIN$cov)
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

    # fit the gams
    gam_D <- gam(obs_response~s(driver,k=knots,bs=smooth_type),data = datIN) # no cov
    gam_DC <- gam(obs_response~s(driver,k=knots,bs=smooth_type) + s(cov,k=knots,bs=smooth_type),data = datIN) # nonlinear driver and cov
    gam_DlC <- gam(obs_response~driver + s(cov,k=knots,bs=smooth_type),data = datIN) # linear driver, nonlinear cov
    gam_DCl <- gam(obs_response~s(driver,k=knots,bs=smooth_type) + cov,data = datIN) # nonlinear driver, linear cov


    # fit the linear models
    lm_D <- lm(obs_response ~ driver, data = datIN) # no cov
    lm_DC <- lm(obs_response ~ driver + cov, data = datIN) # driver and cov

    # get the AIC values
    if(length(datIN$driver) >= aicc_size){
      all_aic <- AIC(gam_D, gam_DC, gam_DCl, gam_DlC, lm_D, lm_DC)$AIC
    } else{
      all_aic <- AICc(gam_D, gam_DC, gam_DCl, gam_DlC, lm_D, lm_DC)$AICc
    }



    mod_names <- c("gam_D", "gam_DC", "gam_DCl", "gam_DlC", "lm_D", "lm_DC") # names of the different models

    best_mod[i] <- mod_names[which(all_aic==min(all_aic))] # model with lowest AIC

    next_aic <- which(all_aic==min(all_aic[-which(all_aic==min(all_aic))])) # next-lowest AIC

    next_mod[i] <- mod_names[next_aic] # model with next-lowest AIC

    # aic_diff = difference btw best and next-best model
    aic_diff[i] <- all_aic[next_aic] - all_aic[which(all_aic==min(all_aic))]

    min_gam <- min(all_aic[1:3]) # min AIC of all models with a nonlinear driver term
    min_lm <- min(all_aic[4:6]) # min AIC of all models with a linear driver term (includes gam w/ linear driver term)

    # add ns if difference btw best model w/ nonlinear driver and best model with linear driver isn't >=2
    if(abs(min_lm - min_gam) <2){
        best_mod[i] <- paste(best_mod[i], "ns", sep = "_")
    }

  }

  # save the results as a dataframe

  thresh_methods_ss <- rep(thresh_methods, each = length(sig_criteria))
  nsim_ss <- rep(c(1:nsim), each = length(sig_criteria))
  best_mod_ss <- rep(best_mod, each = length(sig_criteria)) # best model and AIC differences are the same for all the threshold definitions and significance criteria
  next_mod_ss <- rep(next_mod, each = length(sig_criteria))
  aic_diff_ss <- rep(aic_diff, each = length(sig_criteria))

  df <- data.frame(
    sim = rep(nsim_ss, length(thresh_methods)),
    best_mod = rep(best_mod_ss, length(thresh_methods)),
    next_mod = rep(next_mod_ss, length(thresh_methods)),
    aic_diff = rep(aic_diff_ss, length(thresh_methods)),
    thresh_method = rep(thresh_methods_ss, each = nsim),
    sig_type = rep(sig_criteria, nsim*length(thresh_methods))
  )


  return(df)

}


