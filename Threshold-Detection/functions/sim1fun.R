

#' function for simulating driver-response data with threshold relationship

#' function inputs:
#' nsim = number of simulations
#' tmax = length of simulated time series
#' x_mean = mean of x variable (driver)
#' x_sd = sd of x variable
#' x_acf = temporal autocorrelation of x variable
#' xy_fun = function for relationship between driver (x) and response (y)
#' p_error = process error (noise in relationship between x and y)
#' obs_error = observation error (noise in observation of y values)
#' seed = initial seed for random number generation
#' y_pos = should y be restricted to positive values? (TRUE = restrict to pos values, FALSE = y can be negative)
#' x_max = max x value allowed (to prevent y from being negative)
# 'x_min = min x value allowed (to prevent y from being negative)

simfun1 <- function(nsim.f, tmax.f, x_mean.f, x_sd.f, x_acf.f, xy_fun, p_error.f, obs_error.f, seed.f, y_pos.f, x_max.f, x_min.f){

  for(s in 1:nsim.f){ # for each simulation

    # set the seed
    set.seed(seed.f+s^2)

    # holding vector of x values
    xset <- rep(NA, tmax.f)

    if(y_pos.f==FALSE){
      xset[1] <-rnorm(1, mean = x_mean.f, sd = x_sd.f) # first value of xset: single value drawn from normal dist with mean = x_mean.f and sd = x_sd.f,
    } else{ # if y is restricted to positive values than restrict this to the range between x intercept(s)
      xset[1] <- min(max(x_min.f, rnorm(1, mean = x_mean.f, sd = x_sd.f)), x_max.f)
    }

    # holding vector of y values
    yset <- rep(NA, tmax.f)
    yset[1] <- rnorm(1, mean = xy_fun(xset[1]), sd = p_error.f) # first value of yset: single value drawn from normal dist with mean = the value of y corresponding to x = xset[1], as given by the relationship in xy_fun, and sd = process error (the noise in this relationship between x and y)

    # holding vector of observed y values
    yobset <- rep(NA, tmax.f)
    yobset[1] <- rnorm(1, mean = yset[1], sd = obs_error.f) # first value of yobset: single value drawn from normal dist with mean = yset[1] (true value of y) and sd = observation error

    for(t in 2:tmax.f){ # for each time point starting at t = 2


      # get the value of x (driver)
      if(y_pos.f==FALSE){
        xset[t] <-x_mean.f + x_acf.f*xset[t-1] + rnorm(1, mean = 0, sd = x_sd.f)
      } else{ # if y is restricted to positive values than restrict this to the range between x intercept(s)
        xset[t] <- min(max(x_min.f, x_mean.f + x_acf.f*xset[t-1] + rnorm(1, mean = 0, sd = x_sd.f)), x_max.f)
      }

      # get the true value of y (response)
      yset[t] <- rnorm(1, mean = xy_fun(xset[t]), sd = p_error.f)

      # get the observed value of y
      if(y_pos.f==FALSE){

        yobset[t] <- rnorm(1, mean = yset[t], sd = obs_error.f)# observed y is drawn from normal dist with mean = true value of y at time t and sd = observation error

      } else{ # if yobs should be restricted to positive values

        yobset[t] <- max(0, rnorm(1, mean = yset[t], sd = obs_error.f)) # restricted to be greater than or equal to zero
      }

      #}

    }


    # make data frame of all the results for this simulation
    df <- data.frame(t = 1:tmax.f, # time steps
                     driver = xset, # value of driver at each time step
                     response = yset, # true value of response at each time step
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
