#' function for simulating driver-response data with threshold relationship
#' changes from simfun.R: added option sigmoidal function and modified the hockey stick function so
#' it could be either Beverton-Holt or decrease then flat; also modified the gamma function so the curve
#' can be concave up or concave down

#' should also add option for missing covariate

#' @param nsim number of simulations
#' @param tmax length of simulated time series
#' @param driver_pars A named list controlling the driver (predictor) values simulated, including
#' `thresh_quant` (which quantile of the driver's distribution the threshold falls in; used for
#' calculating the mean of the x values), `x_min` (the minimum), `x_max` (the maximum), `x_sd`
#' (the standard deviation, defaults to 1), and `x_df` (the Student-t df, defaults to 100 -- wider
#' tailed distributions are modeled by making this small ~ 3). Finally, the `uniform` boolean
#' argument allows uniform distributions to be used instead (defaults to FALSE)
#' @param acf = temporal autocorrelation of x variable
#' @param fun function for driver-response relationship. Can be `skew` (values are simulated from a Gamma density, described by the
#' mode and variance), `hockeystick` , `sigmoidal`, or `linear`. Parameters for each can be modified
#' @param obs_sd = observation error sd (noise in observation of y values), defaults to 0.1
#' @param acf ACF of residual error, defaults to 0
#' @param seed initial seed for random number generation
#' @param control_pars A named list controlling the threshold related parameters. There are four functions
#' allowed (`skew` for quadratic / skewed distributions, `hockeystick` for breakpoint models, `sigmoidal`, and `linear`).
#' All functions have a specified threshold location (`thresh_loc`, defaults to 3) and simulations can
#' have unique values of this (controlled by `thresh_loc_sd`, defaults to 0). The functions also share
#' a maximum value (`y_max`, defaults to 3). The hockey stick parameter `hs_a` represents an optional
#' lower x bound below which the function value becomes constant. The `skew_cv` and `skew_conc` arguments are only used for
#' the skewed distribution, and represents the CV of the Gamma distribution and whether the curve is concave up or concave down, respectively.
#' The `sig_k` argument is only used for the sigmoidal function, and represents the steepness of the curve; positive values mean y decreases
#' w/ increasing x values, while negative values mean y increases with increasing x values.
#' The `lin_m` and `lin_b` arguments are used only for the linear function and represent the slope and intercept, respectively
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
                                       lin_b = 0)){

  for(s in 1:nsim){ # for each simulation

    # set the seed
    set.seed(seed+s^2)

    # simulate driver values
    if(driver_pars$uniform == FALSE) {
      xset <- rt(tmax, df = driver_pars$x_df)
      x_mean <- control_pars$thresh_loc - qt(driver_pars$thresh_quant, df = driver_pars$x_df)
      xset <- x_mean + xset * driver_pars$x_sd
      #xset <- driver_pars$x_mean + xset * driver_pars$x_sd
    } else {
      x_mean <- control_pars$thresh_loc - qunif(driver_pars$thresh_quant, df = driver_pars$x_df)
      a <- x_mean - (driver_pars$x_sd * sqrt(12)) / 2
      b <- x_mean + (driver_pars$x_sd * sqrt(12)) / 2
      #a <- driver_pars$x_mean - (driver_pars$x_sd * sqrt(12)) / 2
      #b <- driver_pars$x_mean + (driver_pars$x_sd * sqrt(12)) / 2
      xset <- runif(tmax, min = a, max = b)
    }

    if(fun == "skew") {
      # generate with a lognormal density
      # parameterized in terms of the mode, exp(u-sigma2) and sd
      thresh_loc_true <- rnorm(1, control_pars$thresh_loc, control_pars$thresh_loc_sd)
      skew_var = (control_pars$skew_cv * control_pars$thresh_loc) ^ 2
      pars <- solve_gamma(thresh_loc_true, skew_var) # get the parameters of the gamma dist
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
        # x is vector of predictors
        # type is the type of breakpt relationship (1 = Beverton Holt, 2= decrease then flat)
        # a is value below which the relationship = 0 [optional, for type = 1]
        # thresh_loc is value of breakpt on x-axis
        # y_max is maximum y value
        # y_min is minimum y value
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
      thresh_loc_true <- rnorm(1, control_pars$thresh_loc, control_pars$thresh_loc_sd)
      yset <- hs(xset, y_max = control_pars$y_max, y_min = control_pars$y_min, thresh_loc = thresh_loc_true, type = control_pars$hs_type, a = control_pars$hs_a)
    }

    if(fun == "sigmoidal"){

      sig <- function(x, y_max, thresh_loc, k){
        y <- y_max/(1 + exp(-k*(x-thresh_loc)))
        return(y)
      }

      thresh_loc_true <- rnorm(1, control_pars$thresh_loc, control_pars$thresh_loc_sd)
      yset <- sig(xset, y_max = control_pars$y_max, thresh_loc = thresh_loc_true, k = control_pars$sig_k)

    }

    if(fun == "linear"){

      thresh_loc_true <- NA
      yset <-  control_pars$lin_m*xset + control_pars$lin_b

    }

    # add normal observation error here -- acf controls mean reversion
    dev <- rnorm(1, mean = 0, sd = obs_sd)
    for(t in 2:tmax) {
      dev[t] <- rnorm(1, mean = acf*dev[t-1], sqrt(1 - acf*acf) * obs_sd)
    }

    # add missing covariate

    # add residual / observation error
    yobset <- yset + dev

    # make data frame of all the results for this simulation
    df <- data.frame(t = 1:tmax, # time steps
                     driver = xset, # value of driver at each time step
                     response = yset, # true value of response at each time step
                     obs_response = yobset, # observed value of response at each time step,
                     thresh_loc = thresh_loc_true,
                     sim = s # which simulation this was
    )
    if(s==1) { # if this was the first simulation
      all_df <- df # save this df as all_df
    } else { # otherwise
      all_df <- rbind(df,all_df) # rbind the results for simulation s to all_df (so all_df will have the results for all simulations, which simulation identified by the "sim" column)
    }

  }

  # return the data frame with the results for all the simulations
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

