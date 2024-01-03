

#' function for calculating the number of thresholds detected and the difference between those
#' thresholds and the true value of the threshold
#'
#' inputs:
#' nsim.f: number of simulations that were run
#' thresh_true.f: true value of the threshold
#' threshfunOUT.f: output from threshfun()
#' thresh_choice.f: in case there are multiple thresholds detected, specify which of these (where 1st =
#' smallest and last = largest) should be used for calculating difference from true value

#' NOTE: the current version of the function for calculating thresholds, threshfun(), should only
#' return a single value of a threshold for each calculation method (since the code specifies
#' that if there are more than one, only return the max/min), but the function below allows for
#' the possibility that there were multiple thresholds detected in a given simulation

# function of number of simulations, the true threshold value,  and the outputs of threshfun()
threshdiff_fun <- function(nsim.f, thresh_true.f, threshfunOUT.f, thresh_choice.f){

  # make holding vectors for number of thresholds detected (for each calculation method)
  nthresh1 <- rep(NA, nsim.f)
  nthresh2 <- rep(NA, nsim.f)
  nthresh20 <- rep(NA, nsim.f)

  # make holding vectors for value of threshold (for each calculation method)
  thresh1 <- rep(NA, nsim.f)
  thresh2 <- rep(NA, nsim.f)
  thresh20 <- rep(NA, nsim.f)

  for(i in 1:nsim.f){

    # get the threshold values
    thresh_1 <- threshfunOUT.f$thresh_1[[i]] # list of threshold values based on max/min (and significantly nonzero) in first deriv for ith simulation
    thresh_2 <- threshfunOUT.f$thresh_2[[i]] # list of threshold values based on second deriv being furthest from zero (and significantly nonzero) for ith simulation
    thresh_20 <- threshfunOUT.f$thresh_20[[i]] # list of threshold values based on second deriv CIs being closest to zero for ith simulation

    # start with the threshold calculated based on first derivative
    if(length(thresh_1)==0){ # if no threshold was detected

      nthresh1[i] <- 0

    } else{ # otherwise

      nthresh1[i] <- length(thresh_1) # record number of thresholds detected

      if(nthresh1[i] >= thresh_choice.f){ # if there are at least thresh_choice.f thresholds detected

        # sort all the thresholds from smallest to largest
        thresh_1s <- sort(threshfunOUT.f$thresh_1[[i]])

        thresh1[i] <- thresh_1s[thresh_choice.f] # save the threshold corresponding to thresh_choice.f (e.g., if thresh_choice.f = 1, then thresh1[i] would be the smallest of all the thresholds detected, etc.)

      } else{ # if there aren't at least thresh_choice.f thresholds detected,

        thresh1[i] <- NA # the thresh_choice.fth threshold is NA

      }
    }

    # repeat for the other two threshold calculations
    # based on max/min in second deriv
    if(length(thresh_2)==0){ # if no threshold was detected

      nthresh2[i] <- 0

    } else{ # otherwise

      nthresh2[i] <- length(thresh_2) # record number of thresholds detected

      if(nthresh2[i] >= thresh_choice.f){ # if there are at least thresh_choice.f thresholds detected

        # sort all the thresholds from smallest to largest
        thresh_2s <- sort(threshfunOUT.f$thresh_2[[i]])

        thresh2[i] <- thresh_2s[thresh_choice.f] # save the threshold corresponding to thresh_choice.f

      } else{ # if there aren't at least thresh_choice.f thresholds detected,

        thresh2[i] <- NA # the thresh_choice.fth threshold is NA

      }
    }

    # based on second deriv being closest to 0
    if(length(thresh_20)==0){ # if no threshold was detected

      nthresh20[i] <- 0

    } else{ # otherwise

      nthresh20[i] <- length(thresh_20) # record number of thresholds detected

      if(nthresh20[i] >= thresh_choice.f){ # if there are at least thresh_choice.f thresholds detected

        # sort all the thresholds from smallest to largest
        thresh_20s <- sort(threshfunOUT.f$thresh_20[[i]])

        thresh20[i] <- thresh_20s[thresh_choice.f] # save the threshold corresponding to thresh_choice.f

      } else{ # if there aren't at least thresh_choice.f thresholds detected,

        thresh20[i] <- NA # the thresh_choice.fth threshold is NA

      }
    }

  }

  # make data frames where first col = simulation number, second col = number of thresholds detected, third col = true value of threshold, 4th col = estimated threshold value, 5th col = difference

  thresh1df <- data.frame(
    sim = c(1:nsim.f), # simulation number
    type = rep("max_d1", nsim.f), # type of threshold calculation used
    nthresh = nthresh1, # number of thresholds detected
    threshT = rep(thresh_true.f, nsim.f), # true threshold value
    threshE = thresh1, # estimated value of threshold (if nthresh>1, this is the one closest to true value)
    threshdiff = rep(thresh_true.f, nsim.f) - thresh1 # difference btw true and estimated threshold

  )


  thresh2df <- data.frame(
    sim = c(1:nsim.f), # simulation number
    type = rep("max_d2", nsim.f), # type of threshold calculation used
    nthresh = nthresh2, # number of thresholds detected
    threshT = rep(thresh_true.f, nsim.f), # true threshold value
    threshE = thresh2, # estimated value of threshold (if nthresh>1, this is the one closest to true value)
    threshdiff = rep(thresh_true.f, nsim.f) - thresh2 # difference btw true and estimated threshold

  )

  thresh20df <- data.frame(
    sim = c(1:nsim.f), # simulation number
    type = rep("zero_d2", nsim.f), # type of threshold calculation used
    nthresh = nthresh20, # number of thresholds detected
    threshT = rep(thresh_true.f, nsim.f), # true threshold value
    threshE = thresh20, # estimated value of threshold (if nthresh>1, this is the one closest to true value)
    threshdiff = rep(thresh_true.f, nsim.f) - thresh20 # difference btw true and estimated threshold

  )


  return(list(thresh1df = thresh1df, thresh2df = thresh2df, thresh20df = thresh20df))


  # return data frames with results for each type of threshold calculation

}
