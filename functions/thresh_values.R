#' functions to return the values of the thresholds found in either the full data set (full_thresh)
#' or each individual jackknife iteration in a specified simulation (ind_thresh)

#' ind_thresh returns the threshold values calculated for each jackknife iteration in a specified simulation
#'@param ind_dfs data frame with the jackknifed derivative values and their simultaneous intervals
#'for each jackknife iteration (an output of jack_results function)
#'@param summ_dfs data frame with the 95% quantiles of the jackknifed derivatives for each jackknife
#'iteration (an output of jack_results function)
#'@param simx simulation replicate for which to return the jackknifed thresholds
#' other parameters are as defined for the jack_thresh function

ind_thresh <- function(ind_dfs, summ_dfs, simx = 1,
                       thresh_methods = c("abs_max_d2", "min_d2", "zero_d2", "zero_d1"),
                       sig_criteria = c("none", "jack_quant", "sim_int"),
                       pk_height = 0.01, pk_ups = 3){ #pk_height = 0.05, -Inf

xvals <- summ_dfs$driver[which(summ_dfs$sim==simx)]

# findpeaks parameters to put into the functions for finding maxes and mins
pk_height1 <- pk_height
pk_ups1 <- pk_ups


# holding list for each sig criteria
all_mats <- vector(mode  = "list", length = length(sig_criteria))

for(ss in 1:length(sig_criteria)){

  sig_criteria_s <- sig_criteria[ss]

  # holding matrix where col = calc method and row = jackknife iteration
  thresh_mat <- matrix(NA, nrow = length(unique(ind_dfs$jack_int)), ncol = length(thresh_methods))

  for(int in 1:length(unique(ind_dfs$jack_int))){

    D2int <- ind_dfs$d2[which(ind_dfs$sim==simx & ind_dfs$jack_int==int)]
    D1int <- ind_dfs$d1[which(ind_dfs$sim==simx & ind_dfs$jack_int==int)]

    # individual simultanenous intervals to use for sig_criteria = sim_int
    D2_low_int <- ind_dfs$d2_low[which(ind_dfs$sim==simx & ind_dfs$jack_int==int)]
    D2_up_int <- ind_dfs$d2_up[which(ind_dfs$sim==simx & ind_dfs$jack_int==int)]

    D1_low_int <- ind_dfs$d1_low[which(ind_dfs$sim==simx & ind_dfs$jack_int==int)]
    D1_up_int <- ind_dfs$d1_up[which(ind_dfs$sim==simx & ind_dfs$jack_int==int)]


    # jackknifed CIs to use for sig_criteria = jack_quant
    d2_low <- summ_dfs$low[which(summ_dfs$sim==simx & summ_dfs$output=="d2")]
    d2_up <- summ_dfs$up[which(summ_dfs$sim==simx & summ_dfs$output=="d2")]

    d1_low <- summ_dfs$low[which(summ_dfs$sim==simx & summ_dfs$output=="d1")]
    d1_up <- summ_dfs$up[which(summ_dfs$sim==simx & summ_dfs$output=="d1")]


    # then apply the significant criteria to all the thresholds that were detected
    for(mm in 1:length(thresh_methods)){

      thresh_method <- thresh_methods[mm]

      if(thresh_method =="abs_max_d2"){

        if(sig_criteria_s=="none"){
          thresh_int <- abs_max_threshF(D2int, sig_type = "none", xvals, deriv_up = d2_up, deriv_low = d2_low, pkht = pk_height1, pkup = pk_ups1)
          thresh_mat[int, mm] <- ifelse(is.na(thresh_int)==F, thresh_int, NA)
        }

        if(sig_criteria_s == "jack_quant"){
          thresh_int <- abs_max_threshF(D2int, sig_type = "jack_quant", xvals, deriv_up = d2_up, deriv_low = d2_low, pkht = pk_height1, pkup = pk_ups1)
          thresh_mat[int, mm] <- ifelse(is.na(thresh_int)==F, thresh_int, NA) #ifelse(length(thresh_int)>0, thresh_int, NA)
        }

        if(sig_criteria_s=="sim_int"){
          thresh_int <- abs_max_threshF(D2int, sig_type = "sim_int", xvals, deriv_up = D2_up_int, deriv_low = D2_low_int, pkht = pk_height1, pkup = pk_ups1)
          thresh_mat[int, mm] <- ifelse(is.na(thresh_int)==F, thresh_int, NA)
        }

      }

      if(thresh_method =="min_d2"){

        if(sig_criteria_s=="none"){
          thresh_int <- min_threshF(D2int, sig_type = "none", xvals, deriv_up = d2_up, deriv_low = d2_low, pkht = pk_height1, pkup = pk_ups1)
          thresh_mat[int, mm] <- ifelse(is.na(thresh_int)==F, thresh_int, NA)
        }

        if(sig_criteria_s == "jack_quant"){
          thresh_int <- min_threshF(D2int, sig_type = "jack_quant", xvals, deriv_up = d2_up, deriv_low = d2_low, pkht = pk_height1, pkup = pk_ups1)
          thresh_mat[int, mm] <- ifelse(is.na(thresh_int)==F, thresh_int, NA)

        }

        if(sig_criteria_s=="sim_int"){
          thresh_int <- min_threshF(D2int, sig_type = "sim_int", xvals, deriv_up = D2_up_int, deriv_low = D2_low_int, pkht = pk_height1, pkup = pk_ups1)
          thresh_mat[int, mm] <- ifelse(is.na(thresh_int)==F, thresh_int, NA)
        }

      }

      if(thresh_method =="zero_d2"){

        if(sig_criteria_s=="none"){
          thresh_int <- root_threshF(D2int, sig_type = "none", xvals, deriv_up = d2_up, deriv_low = d2_low)
          thresh_mat[int, mm] <- ifelse(is.na(thresh_int)==F, thresh_int, NA)
        }

        if(sig_criteria_s == "jack_quant"){
          thresh_int <- root_threshF(D2int, sig_type = "jack_quant", xvals, deriv_up = d2_up, deriv_low = d2_low)
          thresh_mat[int, mm] <- ifelse(is.na(thresh_int)==F, thresh_int, NA)

        }

        if(sig_criteria_s=="sim_int"){
          thresh_int <- root_threshF(D2int, sig_type = "sim_int", xvals, deriv_up = D2_up_int, deriv_low = D2_low_int)
          thresh_mat[int, mm] <- ifelse(is.na(thresh_int)==F, thresh_int, NA)
        }


      }


      if(thresh_method =="zero_d1"){

        if(sig_criteria_s=="none"){
          thresh_int <- root_threshF(D1int, sig_type = "none", xvals, deriv_up = d1_up, deriv_low = d1_low)
          thresh_mat[int, mm] <- ifelse(is.na(thresh_int)==F, thresh_int, NA)
        }

        if(sig_criteria_s == "jack_quant"){
          thresh_int <- root_threshF(D1int, sig_type = "jack_quant", xvals, deriv_up = d1_up, deriv_low = d1_low)
          thresh_mat[int, mm] <- ifelse(is.na(thresh_int)==F, thresh_int, NA)

        }

        if(sig_criteria_s=="sim_int"){
          thresh_int <- root_threshF(D1int, sig_type = "sim_int", xvals, deriv_up = D1_up_int, deriv_low = D1_low_int)
          thresh_mat[int, mm] <- ifelse(is.na(thresh_int)==F, thresh_int, NA)
        }
      }

    } # end of threshold calculations mm^th method

  } # end of calculations for int^th jackknife iteration

  all_mats[[ss]] <- thresh_mat

} # end of calculations for ss^th significance criteria

return(all_mats)

} # end of function


#' full_thresh returns the values of the thresholds calculated on the full data set for the specified
#' simulation replicate
#'@param full_dfs data frame with the derivative values and their simultaneous intervals
#'for the full simulated data set
#'@param simx simulation replicate for which to return the thresholds
#' other parameters are as defined for the jack_thresh function

full_thresh <- function(full_dfs, simx = 1,
                        thresh_methods = c("abs_max_d2", "min_d2", "zero_d2", "zero_d1"),
                        sig_criteria = c("none", "jack_quant", "sim_int"),
                        pk_height = 0.01, pk_ups = 3){ #pk_height = 0.05

  xvals <- full_dfs$driver[which(full_dfs$sim==simx)]

  # findpeaks parameters to put into the functions for finding maxes and mins
  pk_height1 <- pk_height
  pk_ups1 <- pk_ups


  # holding list for each sig criteria
  all_vals <- vector(mode  = "list", length = length(sig_criteria))

  for(ss in 1:length(sig_criteria)){

    sig_criteria_s <- sig_criteria[ss]

    # holding vector for thresholds for each calculation method
    thresh_vals <- rep(NA, length(thresh_methods))

    D2_full_mn <- full_dfs$mn[which(full_dfs$sim==simx & full_dfs$output=="d2")]
    D1_full_mn <- full_dfs$mn[which(full_dfs$sim==simx & full_dfs$output=="d1")]

    D2_full_low <- full_dfs$low[which(full_dfs$sim==simx & full_dfs$output=="d2")]
    D2_full_up <- full_dfs$up[which(full_dfs$sim==simx & full_dfs$output=="d2")]

    D1_full_low <- full_dfs$low[which(full_dfs$sim==simx & full_dfs$output=="d1")]
    D1_full_up <- full_dfs$up[which(full_dfs$sim==simx & full_dfs$output=="d1")]

    # then for each threshold method, calculate the threshold and if there are more than one, pick the one closest to the mean value in from the jackknifing threshold calculations
    for(mm in 1:length(thresh_methods)){

      thresh_method <- thresh_methods[mm]

      if(thresh_method =="abs_max_d2"){

        if(sig_criteria_s=="none"){
          d2pks_full1 <- abs_max_threshF(D2_full_mn, sig_type = "none", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low, pkht = pk_height1, pkup = pk_ups1)
          thresh_vals[mm] <- ifelse(is.na(d2pks_full1)==F, d2pks_full1, NA)
        }

        if(sig_criteria_s == "jack_quant"){
          thresh_vals[mm] <- NA # can't evaluate w/o jackknifing
        }

        if(sig_criteria_s=="sim_int"){
          d2pks_full3 <- abs_max_threshF(D2_full_mn, sig_type = "sim_int", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low, pkht = pk_height1, pkup = pk_ups1)
          thresh_vals[mm] <- ifelse(is.na(d2pks_full3)==F, d2pks_full3, NA)
        }
      }

      if(thresh_method =="min_d2"){

        # calculations
        if(sig_criteria_s=="none"){
          d2mins_full1 <- min_threshF(D2_full_mn, sig_type = "none", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low, pkht = pk_height1, pkup = pk_ups1)
          thresh_vals[mm] <- ifelse(is.na(d2mins_full1)==F, d2mins_full1, NA)
        }

        if(sig_criteria_s == "jack_quant"){
          thresh_vals[mm] <- NA # can't evaluate w/o jackknifing
        }

        if(sig_criteria_s=="sim_int"){
          d2mins_full3 <- min_threshF(D2_full_mn, sig_type = "sim_int", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low, pkht = pk_height1, pkup = pk_ups1)
          thresh_vals[mm] <- ifelse(is.na(d2mins_full3)==F, d2mins_full3, NA)
        }

      }

      if(thresh_method =="zero_d2"){

        if(sig_criteria_s=="none"){
          d2roots_full1 <- root_threshF(D2_full_mn, sig_type = "none", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low)
          thresh_vals[mm] <- ifelse(is.na(d2roots_full1)==F, d2roots_full1, NA)
        }

        if(sig_criteria_s == "jack_quant"){
          thresh_vals[mm] <- NA # can't evaluate w/o jackknifing
        }

        if(sig_criteria_s=="sim_int"){
          d2roots_full3 <- root_threshF(D2_full_mn, sig_type = "sim_int", xvals, deriv_up = D2_full_up, deriv_low = D2_full_low)
          thresh_vals[mm] <- ifelse(is.na(d2roots_full3)==F, d2roots_full3, NA)
        }
      }


      if(thresh_method =="zero_d1"){

        if(sig_criteria_s=="none"){
          d1roots_full1 <- root_threshF(D1_full_mn, sig_type = "none", xvals, deriv_up = D1_full_up, deriv_low = D1_full_low)
          thresh_vals[mm] <- ifelse(is.na(d1roots_full1)==F, d1roots_full1, NA)
        }

        if(sig_criteria_s == "jack_quant"){
          thresh_vals[mm] <- NA # can't evaluate w/o jackknifing
        }

        if(sig_criteria_s=="sim_int"){
          d1roots_full3 <- root_threshF(D1_full_mn, sig_type = "sim_int", xvals, deriv_up = D1_full_up, deriv_low = D1_full_low)
          thresh_vals[mm] <- ifelse(is.na(d1roots_full3)==F, d1roots_full3, NA)
        }
      }

    } # end of threshold calculations for mm^th method

    all_vals[[ss]] <- thresh_vals

  } # end of calculations for ss^th sig_criteria

  return(all_vals)

} # end of function

