#' function for simulating driver-response data with threshold relationship

#' @param nsim number of simulations
#' @param tmax length of simulated time series
#' @param driver_pars A named list controlling the driver (predictor) values simulated, including
#' `x_mean` (the mean of x, defaults to 0), `x_min` (the minimum), `x_max` (the maximum), `x_sd` (the standard deviation, defaults to 1),
#' and `x_df` (the Student-t df, defaults to 100 -- wider tailed distributions are modeled by making this small ~ 3). Finally,
#' the `uniform` boolean argument allows uniform distributions to be used instead (defaults to FALSE)
#' @param acf = temporal autocorrelation of x variable
#' @param fun function for driver-response relationship. Can be `skew` (values are simulated from a Gamma density, described by the
#' mode and variance) or `hockeystick`. Parameters for each can be modified
#' @param obs_sd = observation error sd (noise in observation of y values), defaults to 0.1
#' @param acf ACF of residual error, defaults to 0
#' @param seed initial seed for random number generation
#' @param control_pars A named list controlling the threshold related parameters. There are two functions
#' allowed (`skew` for quadratic / skewed distributions, and `hockeystick` for breakpoint models). Both
#' functions have a specified threshold location (`thresh_loc`, defaults to 3) and simulations can
#' have unique values of this (controlled by `thresh_loc_sd`, defaults to 0). Both functions also share
#' a maximum value (`y_max`, defaults to 3). The hockey stick parameter `hs_a` represents an optional
#' lower bound below which the function value becomes 0. The `skew_cv` argument is only used for
#' the skewed distribution, and represents the CV of the Gamma distribution
simfun <- function(nsim = 100, tmax = 30,
                    driver_pars = list(x_min = NULL, x_max = NULL, x_mean = 3, x_df = 10, x_sd = 1, uniform = FALSE),
                    fun = "skew",
                    obs_sd = 0.1,
                    acf = 0,
                    seed = 1234,
                    control_pars = list(thresh_loc = 3,
                                        thresh_loc_sd = 0,
                                        y_max = 3,
                                        hs_a = NULL,
                                        skew_cv = 1)){

  for(s in 1:nsim){ # for each simulation

    # set the seed
    set.seed(seed+s^2)

    # simulate driver values
    if(driver_pars$uniform == FALSE) {
      xset <- rt(tmax, df = driver_pars$x_df)
      xset <- driver_pars$x_mean + xset * driver_pars$x_sd
    } else {
      a <- driver_pars$x_mean - (driver_pars$x_sd * sqrt(12)) / 2
      b <- driver_pars$x_mean + (driver_pars$x_sd * sqrt(12)) / 2
      xset <- runif(tmax, min = a, max = b)
    }

    if(fun == "skew") {
      # generate with a lognormal density
      # parameterized in terms of the mode, exp(u-sigma2) and sd
      thresh_loc_true <- rnorm(1, control_pars$thresh_loc, control_pars$thresh_loc_sd)
      skew_var = (control_pars$skew_cv * control_pars$thresh_loc) ^ 2
      pars <- solve_gamma(thresh_loc_true, skew_var)
      yset <- dgamma(xset, shape = pars$shape, scale = pars$scale)
      max_dens <- dgamma(thresh_loc_true, shape = pars$shape, scale = pars$scale)
      yset <- yset * (control_pars$y_max / max_dens)
    }
    if(fun == "hockeystick") {
      hs <- function(x, y_max, thresh_loc, a=NULL) {
        # x is vector of predictors
        # a is value below which the relationship = 0 [optional]
        # x_max is value of breakpt on x-axis
        # y_max is maximum y value
        if(is.null(a)) a <- min(x)
        y <- ifelse(x < a, 0,
                    ifelse(x < thresh_loc, y_max / (thresh_loc - a) * (x - a), y_max))
        return(y)
      }
      thresh_loc_true <- rnorm(1, control_pars$thresh_loc, control_pars$thresh_loc_sd)
      yset <- hs(xset, a = control_pars$hs_a, thresh_loc = thresh_loc_true, y_max = control_pars$y_max)
    }

    # add normal observation error here -- acf controls mean reversion
    dev <- rnorm(1, mean = 0, sd = obs_sd)
    for(t in 2:tmax) {
      dev[t] <- rnorm(1, mean = acf*dev[t-1], sqrt(1 - acf*acf) * obs_sd)
    }

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

