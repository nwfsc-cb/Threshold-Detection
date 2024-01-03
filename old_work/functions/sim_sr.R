

#' function for simulating driver-response data with threshold relationship

#' function inputs:
#' nsim = number of simulations
#' tmax = length of simulated time series
#' x_vec = vector of x-variable (driver, e.g. spawners)
#' rec_sd = sd of x variable
#' rec_acf = temporal autocorrelation of x variable
#' xy_fun = function for relationship between driver (x) and response (y), either "ricker" or "beverton-holt"
#' p_error = process error (noise in relationship between x and y)
#' obs_error = observation error (noise in observation of y values)
#' seed = initial seed for random number generation
#' sr_pars = stock recruit parameters for

sim_sr <- function(nsim.f, tmax.f, x_vec, rec_sd, rec_acf, xy_fun, obs_error.f, seed.f, sr_pars = c(1,1)){

  for(s in 1:nsim.f){ # for each simulation

    # set the seed
    set.seed(seed.f+s^2)

    # holding vector of x values -- these are recruitment deviations
    xset <- rep(NA, tmax.f)
    xset[1] <- rnorm(1, mean = 0, sd = rec_sd) # first value of xset: single value drawn from normal dist with mean = x_mean.f and sd = x_sd.f,
    for(t in 2:tmax.f){ # for each time point starting at t = 2
      # get the value of x (driver)
      xset[t] <- rec_acf*xset[t-1] + sqrt(1 - rec_acf*rec_acf) * rnorm(1, mean = 0, sd = rec_sd)
    }

    # generate predicted values w/o stochasticity
    if(xy_fun == "ricker") {
      pred_rec <- ricker(x_vec, sr_pars[1], sr_pars[2])
    } else {
      pred_rec <- beverton_holt(x_vec, sr_pars[1], sr_pars[2])
    }

    # now add stochastic recruitment -- lognormal process error
    # the rec_sd^2 / 2 is correcting for bias
    pred_rec <- pred_rec * exp(xset - rec_sd*rec_sd/2)

    # add obs error
    yobset <- pred_rec * exp(rnorm(length(pred_rec), 0, obs_error.f) - obs_error.f*obs_error.f/2)

    # make data frame of all the results for this simulation
    df <- data.frame(t = 1:tmax.f, # time steps
                     driver = xset, # value of driver at each time step
                     response = pred_rec, # true value of response at each time step
                     obs_response = yobset, # observed value of response at each time step
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


# Define the Ricker function
ricker <- function(S, a, b) {
  R = S * exp(a - b * S)
  return(R)
}

beverton_holt <- function(S, alpha, beta) {
  R = (alpha * S) / (1 + beta * S)
  return(R)
}

