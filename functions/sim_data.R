#' function for simulating driver-response data with threshold relationship


#' arguments:
#' @param nsim number of replicate simulations to do
#' @param tmax length of each simulated time series
#' @param driver_pars A named list controlling the driver (predictor) values simulated, including
#' `thresh_quant` (which quantile of the driver's distribution the threshold falls in; used for
#' calculating the mean of the x values), `x_min` (the minimum), `x_max` (the maximum), `x_sd`
#' (the standard deviation, defaults to 1), and `x_df` (the Student-t df, defaults to 100 -- wider
#' tailed distributions are modeled by making this small ~ 3). Finally, the `uniform` boolean
#' argument allows uniform distributions to be used instead (defaults to FALSE)
#' @param fun function for driver-response relationship. Can be `skew` (values are simulated from a Gamma density, described by the
#' mode and variance), `hockeystick` , `sigmoidal`, or `linear`. Parameters for each can be modified
#' @param obs_sd = observation error sd (noise in observation of y values), defaults to 0.1
#' @param acf ACF of residual error, defaults to 0
#' @param seed initial seed for random number generation
#' @param control_pars A named list controlling the threshold related parameters. There are four functions
#' allowed (specified by the `fun` paramter above). All functions have a specified threshold location
#' (`thresh_loc`, defaults to 3) and simulations can have unique values of this (controlled by
#' `thresh_loc_sd`, defaults to 0). The functions also share a maximum value (`y_max`, defaults to 3). The hockey stick parameter `hs_a` represents an optional
#' lower x bound below which the function value becomes constant. The `skew_cv` and `skew_conc`
#' arguments are only used for the skewed distribution, and represents the CV of the Gamma
#' distribution and whether the curve is concave up or concave down, respectively. The `sig_k`
#' argument is only used for the sigmoidal function, and represents the steepness of the curve;
#' positive values mean y decreases w/ increasing x values, while negative values mean y increases
#' with increasing x values. The `lin_m` and `lin_b` arguments are used only for the linear function
#' and represent the slope and intercept, respectively
#' @param cov_pars A named list controlling the covariate parameters; `inc_cov` boolean argument for whether
#' to include a missing covariate (defaults to FALSE), `beta_mean` (mean of normal distribution for the slope
#' of the relationship between covariate and response, defaults to 0) and `beta_sd` (sd of this
#' distribution, controls amount of variability in beta values between simulations; defaults to 0). The absolute
#' values of the betas are drawn from this distribution, and `beta_sign` then controls whether to keep them
#' positive (`beta_sign` = 1) or make them negative (`beta_sign` = -1). The covariate is normally distributed
#' with a mean of `cov_mean` and an sd of `cov_sd`. For an interactive covariate, `cov_int` is the location
#' of the threshold when the cov is at its mean value
#' @param cov_type the type covariate-response relationship, can be `linear`(default), `exponential`,
#' or `interactive` (cov affects the relationship between the driver and the response)
sim_data <- function(nsim = 100, tmax = 30,
                   driver_pars = list(x_min = NULL, x_max = NULL, thresh_quant = 0.5, x_df = 10, x_sd = 1, uniform = FALSE),
                   fun = "skew",
                   obs_sd = 0.1,
                   acf = 0,
                   seed = 1234,
                   control_pars = list(thresh_loc = 3,
                                       thresh_loc_sd = 0,
                                       y_max = 3,
                                       y_min = 0,
                                       hs_type = 1,
                                       hs_a = NULL,
                                       skew_cv = 1,
                                       skew_conc = "up",
                                       sig_k = 1.5,
                                       lin_m = 1,
                                       lin_b = 0),
                   cov_pars = list(inc_cov = TRUE,
                                   beta_mean = 1,
                                   beta_sd = 0,
                                   beta_sign = 1,
                                   cov_mean = 0,
                                   cov_sd = 0,
                                   cov_int = NA),# location of threshold when cov is at mean value (if NA, use thresh_loc)
                   cov_type = "linear"){

  for(s in 1:nsim){ # for each simulation

    # set the seed
    set.seed(seed+s^2)

    # get the true threshold location
    thresh_loc_true <- rnorm(1, control_pars$thresh_loc, control_pars$thresh_loc_sd)

    # simulate driver values
    # NOTE: when cov_type = "interactive", thresh_quant now represents the distance between the center
    # of the driver distribution and the threshold location when the cov is at its mean value (if that
    # mean value is equal to thresh_loc_true); if that mean value differs from thresh_loc_true then
    # on average the distance between the driver center and the threshold at a given time point will
    # differ from thresh_quant

    if(driver_pars$uniform == FALSE) { # if the driver distribution is not uniform
      xset <- rt(tmax, df = driver_pars$x_df) # get the initial sequence of driver values at each time point
      x_mean <- thresh_loc_true - qt(driver_pars$thresh_quant, df = driver_pars$x_df) # calculate the mean of the driver distribution based on the threshold quantile
      xset <- x_mean + xset * driver_pars$x_sd # adjust the driver values based on the mean of the distribution
    } else { # if the driver is uniformly distributed
      x_mean <- thresh_loc_true - qunif(driver_pars$thresh_quant, df = driver_pars$x_df)
      a <- x_mean - (driver_pars$x_sd * sqrt(12)) / 2
      b <- x_mean + (driver_pars$x_sd * sqrt(12)) / 2
      xset <- runif(tmax, min = a, max = b)
    }


    # simulate covariate values
    if(cov_pars$inc_cov == TRUE){ # if including a covariate (x1)

      # get the covariate values
      x1set <- rnorm(tmax, mean = cov_pars$cov_mean, sd = cov_pars$cov_sd)

      # get the value of the slope of relationship btw response and covariate (or steepness, if exponential)
      beta <- cov_pars$beta_sign*abs(rnorm(1, mean = cov_pars$beta_mean, sd = cov_pars$beta_sd))

      if(cov_type == "linear"){
        cov_vals <- beta*x1set
      }

      if(cov_type == "exponential"){
        cov_vals <- exp(beta*x1set) - 1
      }

      if(cov_type == "interactive"){

        # intercept
        cov_int <- ifelse(is.na(cov_pars$cov_int)==T, thresh_loc_true, cov_pars$cov_int)

        cov_vals <- beta*(x1set-cov_pars$cov_mean) + cov_int # cov_vals = threshold; cov_int is the location of
        # the threshold in the driver-response relationship when cov = cov_mean

        thresh_loc_true <- cov_vals # these are the threshold locations at each time point
      }


    } else { # otherwise just input zero for the covariate values
      cov_vals <- rep(0, tmax)
      x1set <- rep(NA, tmax)
    }


    # get the threshold location on the standardized scale
    thresh_loc_true_z <- (thresh_loc_true - mean(xset))/sd(xset)

    # get the normal observation error -- acf controls mean reversion
    # using sd
    dev <- rnorm(1, mean = 0, sd = obs_sd)
    for(t in 2:tmax) {
      dev[t] <- rnorm(1, mean = acf*dev[t-1], sqrt(1 - acf*acf) * obs_sd)
    }

    # using cv
    #dev <- rnorm(1, mean = 0, sd = obs_cv*mean(yset))
    #for(t in 2:tmax) {
    # dev[t] <- rnorm(1, mean = acf*dev[t-1], sqrt(1 - acf*acf) * obs_cv*mean(yset))
    #}

    # relationship between driver and response
    if(cov_type %in% c("linear", "exponential")){


    if(fun == "skew") {
      # generate with a lognormal density
      # parameterized in terms of the mode, exp(u-sigma2) and sd
      skew_var <-  (control_pars$skew_cv * control_pars$thresh_loc) ^ 2
      pars <- solve_gamma(thresh_loc_true, skew_var) # get the parameters of the gamma dist (see end of this .R script)
      max_dens <- dgamma(thresh_loc_true, shape = pars$shape, scale = pars$scale)
      if(control_pars$skew_conc == "down"){
        yset <- dgamma(xset, shape = pars$shape, scale = pars$scale)
        yset <- yset * (control_pars$y_max / max_dens) # scale curve so max = y_max
      }
      if(control_pars$skew_conc == "up"){ # flip the curve and shift it up
        yset <- (-1*dgamma(xset, shape = pars$shape, scale = pars$scale))* (control_pars$y_max/max_dens) + max_dens*(control_pars$y_max/max_dens)
      }

    }
    if(fun == "hockeystick") {
      hs <- function(x, y_max, y_min, thresh_loc, type, a=NULL) {
        # x is vector of predictors, type is the type of breakpt relationship (1 = Beverton Holt,
        # 2= decrease then flat), a is value below which the relationship = 0 [optional, for type = 1]
        # thresh_loc is value of breakpt on x-axis,  y_max is maximum y value,  y_min is min y value
        if(type == 1){
        if(is.null(a)) a <- min(x)
        y <- ifelse(x < a, y_min,
                    ifelse(x < thresh_loc, (y_max-y_min) / (thresh_loc - a) * (x - thresh_loc) + y_max, y_max))
        }

        if(type==2){
          if(is.null(a)) a <- min(x)
          y <- ifelse(x < a, y_max,
                      ifelse(x < thresh_loc, (y_min-y_max) / (thresh_loc - a) * (x - thresh_loc) + y_min, y_min))
          }
        return(y)
      }
      #thresh_loc_true <- rnorm(1, control_pars$thresh_loc, control_pars$thresh_loc_sd)
      yset <- hs(xset, y_max = control_pars$y_max, y_min = control_pars$y_min, thresh_loc = thresh_loc_true, type = control_pars$hs_type, a = control_pars$hs_a)
    }

    if(fun == "sigmoidal"){

      # calculate the inflection point based on the specified threshold location (where the threshold
      # is defined as the location where the 2nd deriv is minimized) -- see end of this .R script
      infl_point <- sig_infl_pt(control_pars$y_max, control_pars$sig_k, thresh_loc_true)

      sig <- function(x, y_max, infl_pt, k){
        y <- y_max/(1 + exp(-k*(x-infl_pt)))
        return(y)
      }

      yset <- sig(xset, y_max = control_pars$y_max, infl_pt = infl_point, k = control_pars$sig_k)

    }

    if(fun == "linear"){
      thresh_loc_true <- NA
      yset <-  control_pars$lin_m*xset + control_pars$lin_b

    }

      # add effect of covariate and residual/observation error
      yobset <- yset + cov_vals + dev

      # make data frame of all the results for this simulation
      df <- data.frame(t = 1:tmax, # time steps
                       driver = xset, # value of driver at each time step
                       cov = x1set, # value of covariate at each time step
                       response = yset, # true value of response at each time step
                       obs_response = yobset, # observed value of response at each time step,
                       thresh_loc = rep(thresh_loc_true, tmax),
                       thresh_loc_z = rep(thresh_loc_true_z, tmax),
                       #beta = rep(beta, tmax), # slope of relationship btw covariate and response
                       sim = rep(s, tmax) # which simulation this was
      )

    }


    if(cov_type == "interactive"){

      yset <- rep(NA, length(thresh_loc_true))

      if(fun == "skew") {
        # generate with a lognormal density
        # parameterized in terms of the mode, exp(u-sigma2) and sd
        skew_var <-  (control_pars$skew_cv * control_pars$thresh_loc) ^ 2

        for(ii in 1:length(thresh_loc_true)){

        pars <- solve_gamma(thresh_loc_true[ii], skew_var) # get the parameters of the gamma dist
        max_dens <- dgamma(thresh_loc_true[ii], shape = pars$shape, scale = pars$scale)
        if(control_pars$skew_conc == "down"){
          yset[ii] <- dgamma(xset[ii], shape = pars$shape, scale = pars$scale)
          yset[ii] <- yset[ii] * (control_pars$y_max / max_dens) # scale curve so max = y_max
        }
        if(control_pars$skew_conc == "up"){ # flip the curve and shift it up
          yset[ii] <- (-1*dgamma(xset[ii], shape = pars$shape, scale = pars$scale))* (control_pars$y_max/max_dens) + max_dens*(control_pars$y_max/max_dens)
        }

      }

      }
      if(fun == "hockeystick") {
        hs <- function(x, y_max, y_min, thresh_loc, type, a=NULL) {
          # x is vector of predictors, type is the type of breakpt relationship (1 = Beverton Holt,
          # 2= decrease then flat), a is value below which the relationship = 0 [optional, for type = 1]
          # thresh_loc is value of breakpt on x-axis,  y_max is maximum y value,  y_min is min y value
          if(type == 1){
            if(is.null(a)) a <- min(x)
            y <- ifelse(x < a, y_min,
                        ifelse(x < thresh_loc, (y_max-y_min) / (thresh_loc - a) * (x - thresh_loc) + y_max, y_max))
          }

          if(type==2){
            if(is.null(a)) a <- min(x)
            y <- ifelse(x < a, y_max,
                        ifelse(x < thresh_loc, (y_min-y_max) / (thresh_loc - a) * (x - thresh_loc) + y_min, y_min))
          }
          return(y)
        }

        for(ii in 1:length(thresh_loc_true)){
        yset[ii] <- hs(xset[ii], y_max = control_pars$y_max, y_min = control_pars$y_min, thresh_loc = thresh_loc_true[ii], type = control_pars$hs_type, a = ifelse(is.null(control_pars$hs_a)==T, min(xset), control_pars$hs_a))
        }
      }

      if(fun == "sigmoidal"){

        # calculate the inflection point based on the specified threshold location (where the threshold
        # is defined as the location where the 2nd deriv is minimized)

        sig <- function(x, y_max, infl_pt, k){
          y <- y_max/(1 + exp(-k*(x-infl_pt)))
          return(y)
        }

        for(ii in 1:length(thresh_loc_true)){
        infl_point <- sig_infl_pt(control_pars$y_max, control_pars$sig_k, thresh_loc_true[ii])
        yset[ii] <- sig(xset[ii], y_max = control_pars$y_max, infl_pt = infl_point, k = control_pars$sig_k)
        }

      }

      if(fun == "linear"){ # no independent effect of cov so linear stays the same
        thresh_loc_true <- NA*thresh_loc_true
        yset <-  control_pars$lin_m*xset + control_pars$lin_b

      }

      # add residual/observation error
      yobset <- yset + dev

      # make data frame of all the results for this simulation
      df <- data.frame(t = 1:tmax, # time steps
                       driver = xset, # value of driver at each time step
                       cov = x1set, # value of covariate at each time step
                       response = yset, # true value of response at each time step
                       obs_response = yobset, # observed value of response at each time step,
                       thresh_loc = thresh_loc_true,
                       thresh_loc_z = thresh_loc_true_z,
                       #beta = rep(beta, tmax), # slope of relationship btw covariate and response
                       sim = rep(s, tmax) # which simulation this was
      )

    }


    if(s==1) { # if this was the first simulation
      all_df <- df # save this df as all_df
    } else { # otherwise
      all_df <- rbind(df,all_df) # rbind the results for simulation s to all_df (so all_df will have the results for all simulations, which simulation identified by the "sim" column)
    }

  }

  # return the data frame with the time series for all the simulations
  return(all_df)

}

# This function just returns the shape/scale of the gamma in terms of the mode and variance
#' @param mode Mode (location where maximum occurs)
#' @param variance variance of Gamma distribution
solve_gamma <- function(mode, variance) {

  # Find k using the variance and mode equations
  find_k <- function(k, m, v) {
    theta_estimate = m / (k - 1)
    v - k * theta_estimate^2
  }
  k_sol = uniroot(find_k, c(1.001, 100), m = mode, v = variance)$root

  # Calculate theta using the solution for k
  theta_sol = mode / (k_sol - 1)
  return(list(shape = k_sol, scale = theta_sol))
}

# This function returns the location of the inflection point of the sigmoidal function in terms of the
# location of the minimum of the 2nd deriv (the threshold location, thresh_loc)

#' @param y_max maximum value of the response
#' @param k steepness of logistic curve
#' @param thresh_loc location of the threshold, defined as the point where the 2nd deriv has a local min

sig_infl_pt <- function(y_max, k, thresh_loc){
 # 2nd deriv has max/min where 3rd deriv = 0, which occurs in 2 places
 # x = (k*ip - log(2 +- sqrt(3)))/k -> solve each of these for the inflection point (ip), then determine which one corresponds to the min
  ip1 <- (k*thresh_loc + log(2 + sqrt(3)))/k
  ip2 <- (k*thresh_loc + log(2 - sqrt(3)))/k

  # see which of these produces a minima in the second deriv at x = thresh_loc
  # minima means 4th deriv at thresh_loc is positive
  d41 <- (24*exp(-4*k*(-ip1 + thresh_loc))*k^4*y_max)/(1 + exp(-k*(-ip1 + thresh_loc)))^5 - (36*exp(-3*k*(-ip1 + thresh_loc))*k^4*y_max)/(1 + exp(-k*(-ip1 + thresh_loc)))^4 + (
    14*exp(-2*k*(-ip1 + thresh_loc))*k^4*y_max)/(1 + exp(-k*(-ip1 + thresh_loc)))^3 - (
      exp(-k*(-ip1 + thresh_loc))*k^4*y_max)/(1 + exp(-k*(-ip1 + thresh_loc)))^2

  d42 <- (24*exp(-4*k*(-ip2 + thresh_loc))*k^4*y_max)/(1 + exp(-k*(-ip2 + thresh_loc)))^5 - (36*exp(-3*k*(-ip2 + thresh_loc))*k^4*y_max)/(1 + exp(-k*(-ip2 + thresh_loc)))^4 + (
    14*exp(-2*k*(-ip2 + thresh_loc))*k^4*y_max)/(1 + exp(-k*(-ip2 + thresh_loc)))^3 - (
      exp(-k*(-ip2 + thresh_loc))*k^4*y_max)/(1 + exp(-k*(-ip2 + thresh_loc)))^2

  if(d41 > 0){
    ip <- ip1
  } else (
    ip <- ip2
  )

  return(ip)

}




