#### 11 - Miscellaneous analyses 

source("R/07_model_analysis/07-01_load_stan_models.R")


###### Lower Granite Spillway impacts ######

# Do we see a year effect for the change at Lower Granite Dam starting following
# when that site came online on 4/3/2020?

# So - there is no random effect of year for the 9 -> 8 movement for Tucannon River fish,
# because this is a post-overshoot fallback movement (and we therefore want to allow
# the spill differences to be captured here)
# So to see the effect of the LGR spillway on this movement, we'd need to just
# look at the raw data

#### Define year movement probability functions ####
estimate_year_move_prob_UCW <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- rep(0,3)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = wild_origin_params[1]
  origin2 = wild_origin_params[2]
  origin3 = wild_origin_params[3]
  
  # set up years to predict across
  year_predict <- seq(1,20)
  
  # get an array to store probabilities of movements in different years
  niter <- 4000 # this is the number of draws we have
  
  # set this up where rownames = years, columns = iter, slices = movements
  year_move_prob_array <- array(dim = c(length(year_predict), niter, nrow(movements)))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    UCW_states_dates_years <- data.frame(state = as.vector(UCW_envir$data$y),
                                         date = as.vector(UCW_envir$data$transition_dates))
    
    UCW_states_dates_years$year <- ceiling(UCW_states_dates_years$date/365.25)+1
    
    # get median spill window for this state, by year - this will be indexed the same as year, using j
    spillwindow_data <- UCW_envir$data$spill_window_data
    
    med_spillwindow <- vector(length = 20)
    
    # for spill and temp, if there were no fish in year/state combination,
    # there will be NAs in the vector - and that's ok
    
    for (x in 1:length(med_spillwindow)){
      med_spillwindow[x] <- median(spillwindow_data[subset(UCW_states_dates_years, state == from & year == x)$date,from])
    }
    
    # get median winter spill for this state (note that this is just one value, since each year already only has one winter spill value), 
    # by year - this will be indexed the same as year, using j
    winterspill_data <- UCW_envir$data$winter_spill_days_data
    med_winterspill <- vector(length = 20)
    
    for (y in 1:length(med_winterspill)){
      med_winterspill[y] <- winterspill_data[y,from]
    }
    
    # get median temperature for this state, by year - this will be indexed the same as year, using j
    temp_data <- UCW_envir$data$temperature_data
    
    med_temp <- vector(length = 20)
    
    for (z in 1:length(med_temp)){
      med_temp[z] <- median(temp_data[subset(UCW_states_dates_years, state == from & year == z)$date,from])
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(year_predict)){
        
        # if there aren't any movements in that year, skip it
        if (!(is.na(med_temp[j]))){
          # evaluation movement 
          year_move_prob_array[j, iter, i] <- exp(b0_array_UCW[from,to,iter] +
                                                    # btemp0_array_UCW[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                    btemp1_array_UCW[from,to,iter]*med_temp[j] + 
                                                    bspillwindow_array_UCW[from,to,iter]*med_spillwindow[j] +
                                                    bwinterspill_array_UCW[from,to,iter]*med_winterspill[j] +
                                                    # btemp0xorigin1_array_UCW[from,to,iter]*origin1 +
                                                    btemp1xorigin1_array_UCW[from,to,iter]*med_temp[j]*origin1 +
                                                    # btemp0xorigin2_array_UCW[from,to,iter]*origin2 + 
                                                    btemp1xorigin2_array_UCW[from,to,iter]*med_temp[j]*origin2 + 
                                                    # btemp0xorigin3_array_UCW[from,to,iter]*origin3 +
                                                    btemp1xorigin3_array_UCW[from,to,iter]*med_temp[j]*origin3 +
                                                    borigin1_array_UCW[from,to,iter]*origin1 +
                                                    borigin2_array_UCW[from,to,iter]*origin2 +
                                                    borigin3_array_UCW[from,to,iter]*origin3 +
                                                    
                                                    # year effects
                                                    origin1_year_param_array_UCW[paste0(from, "_", to), j, iter]*origin1 +
                                                    origin2_year_param_array_UCW[paste0(from, "_", to), j, iter]*origin2 +
                                                    origin3_year_param_array_UCW[paste0(from, "_", to), j, iter]*origin3)/
            sum(exp(b0_array_UCW[from,possible_movements,iter] +
                      # btemp0_array_UCW[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                      btemp1_array_UCW[from,possible_movements,iter]*med_temp[j] + 
                      bspillwindow_array_UCW[from,possible_movements,iter]*med_spillwindow[j] +
                      bwinterspill_array_UCW[from,possible_movements,iter]*med_winterspill[j] +
                      # btemp0xorigin1_array_UCW[from,possible_movements,iter]*origin1 +
                      btemp1xorigin1_array_UCW[from,possible_movements,iter]*med_temp[j]*origin1 +
                      # btemp0xorigin2_array_UCW[from,possible_movements,iter]*origin2 + 
                      btemp1xorigin2_array_UCW[from,possible_movements,iter]*med_temp[j]*origin2 + 
                      # btemp0xorigin3_array_UCW[from,possible_movements,iter]*origin3 +
                      btemp1xorigin3_array_UCW[from,possible_movements,iter]*med_temp[j]*origin3 +
                      borigin1_array_UCW[from,possible_movements,iter]*origin1 +
                      borigin2_array_UCW[from,possible_movements,iter]*origin2 +
                      borigin3_array_UCW[from,possible_movements,iter]*origin3+
                      # year effects
                      c(origin1_year_param_array_UCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin1 +
                      c(origin2_year_param_array_UCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin2 +
                      c(origin3_year_param_array_UCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin3))
          
        } else {
          
        }
        
        
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(year_move_prob_array)
  
}

estimate_year_move_prob_UCH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,3)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  origin3 = hatchery_origin_params[3]
  
  # set up years to predict across
  year_predict <- seq(1,20)
  
  # get an array to store probabilities of movements in different years
  niter <- 4000 # this is the number of draws we have
  
  # set this up where rownames = years, columns = iter, slices = movements
  year_move_prob_array <- array(dim = c(length(year_predict), niter, nrow(movements)))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    UCH_states_dates_years <- data.frame(state = as.vector(UCH_envir$data$y),
                                         date = as.vector(UCH_envir$data$transition_dates))
    
    UCH_states_dates_years$year <- ceiling(UCH_states_dates_years$date/365.25)+1
    
    # get median spill window for this state, by year - this will be indexed the same as year, using j
    spillwindow_data <- UCH_envir$data$spill_window_data
    
    med_spillwindow <- vector(length = 20)
    
    # for spill and temp, if there were no fish in year/state combination,
    # there will be NAs in the vector - and that's ok
    
    for (x in 1:length(med_spillwindow)){
      med_spillwindow[x] <- median(spillwindow_data[subset(UCH_states_dates_years, state == from & year == x)$date,from])
    }
    
    # get median winter spill for this state (note that this is just one value, since each year already only has one winter spill value), 
    # by year - this will be indexed the same as year, using j
    winterspill_data <- UCH_envir$data$winter_spill_days_data
    med_winterspill <- vector(length = 20)
    
    for (y in 1:length(med_winterspill)){
      med_winterspill[y] <- winterspill_data[y,from]
    }
    
    # get median temperature for this state, by year - this will be indexed the same as year, using j
    temp_data <- UCH_envir$data$temperature_data
    
    med_temp <- vector(length = 20)
    
    for (z in 1:length(med_temp)){
      med_temp[z] <- median(temp_data[subset(UCH_states_dates_years, state == from & year == z)$date,from])
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(year_predict)){
        
        # if there aren't any movements in that year, skip it
        if (!(is.na(med_temp[j]))){
          # evaluation movement 
          year_move_prob_array[j, iter, i] <- exp(b0_array_UCH[from,to,iter] +
                                                    # btemp0_array_UCH[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                    btemp1_array_UCH[from,to,iter]*med_temp[j] + 
                                                    bspillwindow_array_UCH[from,to,iter]*med_spillwindow[j] +
                                                    bwinterspill_array_UCH[from,to,iter]*med_winterspill[j] +
                                                    # btemp0xorigin1_array_UCH[from,to,iter]*origin1 +
                                                    btemp1xorigin1_array_UCH[from,to,iter]*med_temp[j]*origin1 +
                                                    # btemp0xorigin2_array_UCH[from,to,iter]*origin2 + 
                                                    btemp1xorigin2_array_UCH[from,to,iter]*med_temp[j]*origin2 + 
                                                    # btemp0xorigin3_array_UCH[from,to,iter]*origin3 +
                                                    btemp1xorigin3_array_UCH[from,to,iter]*med_temp[j]*origin3 +
                                                    borigin1_array_UCH[from,to,iter]*origin1 +
                                                    borigin2_array_UCH[from,to,iter]*origin2 +
                                                    borigin3_array_UCH[from,to,iter]*origin3 +
                                                    
                                                    # year effects
                                                    origin1_year_param_array_UCH[paste0(from, "_", to), j, iter]*origin1 +
                                                    origin2_year_param_array_UCH[paste0(from, "_", to), j, iter]*origin2 +
                                                    origin3_year_param_array_UCH[paste0(from, "_", to), j, iter]*origin3)/
            sum(exp(b0_array_UCH[from,possible_movements,iter] +
                      # btemp0_array_UCH[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                      btemp1_array_UCH[from,possible_movements,iter]*med_temp[j] + 
                      bspillwindow_array_UCH[from,possible_movements,iter]*med_spillwindow[j] +
                      bwinterspill_array_UCH[from,possible_movements,iter]*med_winterspill[j] +
                      # btemp0xorigin1_array_UCH[from,possible_movements,iter]*origin1 +
                      btemp1xorigin1_array_UCH[from,possible_movements,iter]*med_temp[j]*origin1 +
                      # btemp0xorigin2_array_UCH[from,possible_movements,iter]*origin2 + 
                      btemp1xorigin2_array_UCH[from,possible_movements,iter]*med_temp[j]*origin2 + 
                      # btemp0xorigin3_array_UCH[from,possible_movements,iter]*origin3 +
                      btemp1xorigin3_array_UCH[from,possible_movements,iter]*med_temp[j]*origin3 +
                      borigin1_array_UCH[from,possible_movements,iter]*origin1 +
                      borigin2_array_UCH[from,possible_movements,iter]*origin2 +
                      borigin3_array_UCH[from,possible_movements,iter]*origin3+
                      # year effects
                      c(origin1_year_param_array_UCH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin1 +
                      c(origin2_year_param_array_UCH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin2 +
                      c(origin3_year_param_array_UCH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin3))
          
        } else {
          
        }
        
        
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(year_move_prob_array)
  
}

estimate_year_move_prob_MCW <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- rep(0,6)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = wild_origin_params[1]
  origin2 = wild_origin_params[2]
  origin3 = wild_origin_params[3]
  origin4 = wild_origin_params[4]
  origin5 = wild_origin_params[5]
  origin6 = wild_origin_params[6]
  
  # set up years to predict across
  year_predict <- seq(1,20)
  
  # get an array to store probabilities of movements in different years
  niter <- 4000 # this is the number of draws we have
  
  # set this up where rownames = years, columns = iter, slices = movements
  year_move_prob_array <- array(dim = c(length(year_predict), niter, nrow(movements)))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    MCW_states_dates_years <- data.frame(state = as.vector(MCW_envir$data$y),
                                         date = as.vector(MCW_envir$data$transition_dates))
    
    MCW_states_dates_years$year <- ceiling(MCW_states_dates_years$date/365.25)+1
    
    # get median spill window for this state, by year - this will be indexed the same as year, using j
    spillwindow_data <- MCW_envir$data$spill_window_data
    
    med_spillwindow <- vector(length = 20)
    
    # for spill and temp, if there were no fish in year/state combination,
    # there will be NAs in the vector - and that's ok
    
    for (x in 1:length(med_spillwindow)){
      med_spillwindow[x] <- median(spillwindow_data[subset(MCW_states_dates_years, state == from & year == x)$date,from])
    }
    
    # get median winter spill for this state (note that this is just one value, since each year already only has one winter spill value), 
    # by year - this will be indexed the same as year, using j
    winterspill_data <- MCW_envir$data$winter_spill_days_data
    med_winterspill <- vector(length = 20)
    
    for (y in 1:length(med_winterspill)){
      med_winterspill[y] <- winterspill_data[y,from]
    }
    
    # get median temperature for this state, by year - this will be indexed the same as year, using j
    temp_data <- MCW_envir$data$temperature_data
    
    med_temp <- vector(length = 20)
    
    for (z in 1:length(med_temp)){
      med_temp[z] <- median(temp_data[subset(MCW_states_dates_years, state == from & year == z)$date,from])
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(year_predict)){
        
        # if there aren't any movements in that year, skip it
        if (!(is.na(med_temp[j]))){
          # evaluation movement 
          year_move_prob_array[j, iter, i] <- exp(b0_array_MCW[from,to,iter] +
                                                    # btemp0_array_MCW[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                    btemp1_array_MCW[from,to,iter]*med_temp[j] + 
                                                    bspillwindow_array_MCW[from,to,iter]*med_spillwindow[j] +
                                                    bwinterspill_array_MCW[from,to,iter]*med_winterspill[j] +
                                                    # btemp0xorigin1_array_MCW[from,to,iter]*origin1 +
                                                    btemp1xorigin1_array_MCW[from,to,iter]*med_temp[j]*origin1 +
                                                    # btemp0xorigin2_array_MCW[from,to,iter]*origin2 + 
                                                    btemp1xorigin2_array_MCW[from,to,iter]*med_temp[j]*origin2 + 
                                                    # btemp0xorigin3_array_MCW[from,to,iter]*origin3 +
                                                    btemp1xorigin3_array_MCW[from,to,iter]*med_temp[j]*origin3 +
                                                    # btemp0xorigin4_array_MCW[from,to,iter]*origin4 +
                                                    btemp1xorigin4_array_MCW[from,to,iter]*med_temp[j]*origin4 +
                                                    # btemp0xorigin5_array_MCW[from,to,iter]*origin5 +
                                                    btemp1xorigin5_array_MCW[from,to,iter]*med_temp[j]*origin5 +
                                                    # btemp0xorigin6_array_MCW[from,to,iter]*origin6 +
                                                    btemp1xorigin6_array_MCW[from,to,iter]*med_temp[j]*origin6 +
                                                    borigin1_array_MCW[from,to,iter]*origin1 +
                                                    borigin2_array_MCW[from,to,iter]*origin2 +
                                                    borigin3_array_MCW[from,to,iter]*origin3 +
                                                    borigin4_array_MCW[from,to,iter]*origin4 +
                                                    borigin5_array_MCW[from,to,iter]*origin5 +
                                                    borigin6_array_MCW[from,to,iter]*origin6 +
                                                    
                                                    # year effects
                                                    origin1_year_param_array_MCW[paste0(from, "_", to), j, iter]*origin1 +
                                                    origin2_year_param_array_MCW[paste0(from, "_", to), j, iter]*origin2 +
                                                    origin3_year_param_array_MCW[paste0(from, "_", to), j, iter]*origin3 +
                                                    origin4_year_param_array_MCW[paste0(from, "_", to), j, iter]*origin4 +
                                                    origin5_year_param_array_MCW[paste0(from, "_", to), j, iter]*origin5 +
                                                    origin6_year_param_array_MCW[paste0(from, "_", to), j, iter]*origin6)/
            sum(exp(b0_array_MCW[from,possible_movements,iter] +
                      # btemp0_array_MCW[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                      btemp1_array_MCW[from,possible_movements,iter]*med_temp[j] + 
                      bspillwindow_array_MCW[from,possible_movements,iter]*med_spillwindow[j] +
                      bwinterspill_array_MCW[from,possible_movements,iter]*med_winterspill[j] +
                      # btemp0xorigin1_array_MCW[from,possible_movements,iter]*origin1 +
                      btemp1xorigin1_array_MCW[from,possible_movements,iter]*med_temp[j]*origin1 +
                      # btemp0xorigin2_array_MCW[from,possible_movements,iter]*origin2 + 
                      btemp1xorigin2_array_MCW[from,possible_movements,iter]*med_temp[j]*origin2 + 
                      # btemp0xorigin3_array_MCW[from,possible_movements,iter]*origin3 +
                      btemp1xorigin3_array_MCW[from,possible_movements,iter]*med_temp[j]*origin3 +
                      # btemp0xorigin4_array_MCW[from,possible_movements,iter]*origin4 +
                      btemp1xorigin4_array_MCW[from,possible_movements,iter]*med_temp[j]*origin4 +
                      # btemp0xorigin5_array_MCW[from,possible_movements,iter]*origin5 +
                      btemp1xorigin5_array_MCW[from,possible_movements,iter]*med_temp[j]*origin5 +
                      # btemp0xorigin6_array_MCW[from,possible_movements,iter]*origin6 +
                      btemp1xorigin6_array_MCW[from,possible_movements,iter]*med_temp[j]*origin6 +
                      borigin1_array_MCW[from,possible_movements,iter]*origin1 +
                      borigin2_array_MCW[from,possible_movements,iter]*origin2 +
                      borigin3_array_MCW[from,possible_movements,iter]*origin3 +
                      borigin4_array_MCW[from,possible_movements,iter]*origin4 +
                      borigin5_array_MCW[from,possible_movements,iter]*origin5 +
                      borigin6_array_MCW[from,possible_movements,iter]*origin6 +
                      # year effects
                      c(origin1_year_param_array_MCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin1 +
                      c(origin2_year_param_array_MCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin2 +
                      c(origin3_year_param_array_MCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin3 +
                      c(origin4_year_param_array_MCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin4 +
                      c(origin5_year_param_array_MCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin5 +
                      c(origin6_year_param_array_MCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin6))
          
        } else {
          
        }
        
        
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(year_move_prob_array)
  
}

estimate_year_move_prob_MCH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,2)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  
  # set up years to predict across
  year_predict <- seq(1,20)
  
  # get an array to store probabilities of movements in different years
  niter <- 4000 # this is the number of draws we have
  
  # set this up where rownames = years, columns = iter, slices = movements
  year_move_prob_array <- array(dim = c(length(year_predict), niter, nrow(movements)))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    MCH_states_dates_years <- data.frame(state = as.vector(MCH_envir$data$y),
                                         date = as.vector(MCH_envir$data$transition_dates))
    
    MCH_states_dates_years$year <- ceiling(MCH_states_dates_years$date/365.25)+1
    
    # get median spill window for this state, by year - this will be indexed the same as year, using j
    spillwindow_data <- MCH_envir$data$spill_window_data
    
    med_spillwindow <- vector(length = 20)
    
    # for spill and temp, if there were no fish in year/state combination,
    # there will be NAs in the vector - and that's ok
    
    for (x in 1:length(med_spillwindow)){
      med_spillwindow[x] <- median(spillwindow_data[subset(MCH_states_dates_years, state == from & year == x)$date,from])
    }
    
    # get median winter spill for this state (note that this is just one value, since each year already only has one winter spill value), 
    # by year - this will be indexed the same as year, using j
    winterspill_data <- MCH_envir$data$winter_spill_days_data
    med_winterspill <- vector(length = 20)
    
    for (y in 1:length(med_winterspill)){
      med_winterspill[y] <- winterspill_data[y,from]
    }
    
    # get median temperature for this state, by year - this will be indexed the same as year, using j
    temp_data <- MCH_envir$data$temperature_data
    
    med_temp <- vector(length = 20)
    
    for (z in 1:length(med_temp)){
      med_temp[z] <- median(temp_data[subset(MCH_states_dates_years, state == from & year == z)$date,from])
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(year_predict)){
        
        # if there aren't any movements in that year, skip it
        if (!(is.na(med_temp[j]))){
          # evaluation movement 
          year_move_prob_array[j, iter, i] <- exp(b0_array_MCH[from,to,iter] +
                                                    # btemp0_array_MCH[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                    btemp1_array_MCH[from,to,iter]*med_temp[j] + 
                                                    bspillwindow_array_MCH[from,to,iter]*med_spillwindow[j] +
                                                    bwinterspill_array_MCH[from,to,iter]*med_winterspill[j] +
                                                    # btemp0xorigin1_array_MCH[from,to,iter]*origin1 +
                                                    btemp1xorigin1_array_MCH[from,to,iter]*med_temp[j]*origin1 +
                                                    # btemp0xorigin2_array_MCH[from,to,iter]*origin2 + 
                                                    btemp1xorigin2_array_MCH[from,to,iter]*med_temp[j]*origin2 + 
                                                    borigin1_array_MCH[from,to,iter]*origin1 +
                                                    borigin2_array_MCH[from,to,iter]*origin2 +
                                                    
                                                    # year effects
                                                    origin1_year_param_array_MCH[paste0(from, "_", to), j, iter]*origin1 +
                                                    origin2_year_param_array_MCH[paste0(from, "_", to), j, iter]*origin2)/
            sum(exp(b0_array_MCH[from,possible_movements,iter] +
                      # btemp0_array_MCH[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                      btemp1_array_MCH[from,possible_movements,iter]*med_temp[j] + 
                      bspillwindow_array_MCH[from,possible_movements,iter]*med_spillwindow[j] +
                      bwinterspill_array_MCH[from,possible_movements,iter]*med_winterspill[j] +
                      # btemp0xorigin1_array_MCH[from,possible_movements,iter]*origin1 +
                      btemp1xorigin1_array_MCH[from,possible_movements,iter]*med_temp[j]*origin1 +
                      # btemp0xorigin2_array_MCH[from,possible_movements,iter]*origin2 + 
                      btemp1xorigin2_array_MCH[from,possible_movements,iter]*med_temp[j]*origin2 + 
                      borigin1_array_MCH[from,possible_movements,iter]*origin1 +
                      borigin2_array_MCH[from,possible_movements,iter]*origin2 +
                      # year effects
                      c(origin1_year_param_array_MCH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin1 +
                      c(origin2_year_param_array_MCH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin2))
          
        } else {
          
        }
        
        
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(year_move_prob_array)
  
}

estimate_year_move_prob_SRW <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- rep(0,6)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = wild_origin_params[1]
  origin2 = wild_origin_params[2]
  origin3 = wild_origin_params[3]
  origin4 = wild_origin_params[4]
  origin5 = wild_origin_params[5]
  origin6 = wild_origin_params[6]
  
  # set up years to predict across
  year_predict <- seq(1,20)
  
  # get an array to store probabilities of movements in different years
  niter <- 4000 # this is the number of draws we have
  
  # set this up where rownames = years, columns = iter, slices = movements
  year_move_prob_array <- array(dim = c(length(year_predict), niter, nrow(movements)))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    SRW_states_dates_years <- data.frame(state = as.vector(SRW_envir$data$y),
                                         date = as.vector(SRW_envir$data$transition_dates))
    
    SRW_states_dates_years$year <- ceiling(SRW_states_dates_years$date/365.25)+1
    
    # get median spill window for this state, by year - this will be indexed the same as year, using j
    spillwindow_data <- SRW_envir$data$spill_window_data
    
    med_spillwindow <- vector(length = 20)
    
    # for spill and temp, if there were no fish in year/state combination,
    # there will be NAs in the vector - and that's ok
    
    for (x in 1:length(med_spillwindow)){
      med_spillwindow[x] <- median(spillwindow_data[subset(SRW_states_dates_years, state == from & year == x)$date,from])
    }
    
    # get median winter spill for this state (note that this is just one value, since each year already only has one winter spill value), 
    # by year - this will be indexed the same as year, using j
    winterspill_data <- SRW_envir$data$winter_spill_days_data
    med_winterspill <- vector(length = 20)
    
    for (y in 1:length(med_winterspill)){
      med_winterspill[y] <- winterspill_data[y,from]
    }
    
    # get median temperature for this state, by year - this will be indexed the same as year, using j
    temp_data <- SRW_envir$data$temperature_data
    
    med_temp <- vector(length = 20)
    
    for (z in 1:length(med_temp)){
      med_temp[z] <- median(temp_data[subset(SRW_states_dates_years, state == from & year == z)$date,from])
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(year_predict)){
        
        # if there aren't any movements in that year, skip it
        if (!(is.na(med_temp[j]))){
          # evaluation movement 
          year_move_prob_array[j, iter, i] <- exp(b0_array_SRW[from,to,iter] +
                                                    # btemp0_array_SRW[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                    btemp1_array_SRW[from,to,iter]*med_temp[j] + 
                                                    bspillwindow_array_SRW[from,to,iter]*med_spillwindow[j] +
                                                    bwinterspill_array_SRW[from,to,iter]*med_winterspill[j] +
                                                    # btemp0xorigin1_array_SRW[from,to,iter]*origin1 +
                                                    btemp1xorigin1_array_SRW[from,to,iter]*med_temp[j]*origin1 +
                                                    # btemp0xorigin2_array_SRW[from,to,iter]*origin2 + 
                                                    btemp1xorigin2_array_SRW[from,to,iter]*med_temp[j]*origin2 + 
                                                    # btemp0xorigin3_array_SRW[from,to,iter]*origin3 +
                                                    btemp1xorigin3_array_SRW[from,to,iter]*med_temp[j]*origin3 +
                                                    # btemp0xorigin4_array_SRW[from,to,iter]*origin4 +
                                                    btemp1xorigin4_array_SRW[from,to,iter]*med_temp[j]*origin4 +
                                                    # btemp0xorigin5_array_SRW[from,to,iter]*origin5 +
                                                    btemp1xorigin5_array_SRW[from,to,iter]*med_temp[j]*origin5 +
                                                    # btemp0xorigin6_array_SRW[from,to,iter]*origin6 +
                                                    btemp1xorigin6_array_SRW[from,to,iter]*med_temp[j]*origin6 +
                                                    borigin1_array_SRW[from,to,iter]*origin1 +
                                                    borigin2_array_SRW[from,to,iter]*origin2 +
                                                    borigin3_array_SRW[from,to,iter]*origin3 +
                                                    borigin4_array_SRW[from,to,iter]*origin4 +
                                                    borigin5_array_SRW[from,to,iter]*origin5 +
                                                    borigin6_array_SRW[from,to,iter]*origin6 +
                                                    
                                                    # year effects
                                                    origin1_year_param_array_SRW[paste0(from, "_", to), j, iter]*origin1 +
                                                    origin2_year_param_array_SRW[paste0(from, "_", to), j, iter]*origin2 +
                                                    origin3_year_param_array_SRW[paste0(from, "_", to), j, iter]*origin3 +
                                                    origin4_year_param_array_SRW[paste0(from, "_", to), j, iter]*origin4 +
                                                    origin5_year_param_array_SRW[paste0(from, "_", to), j, iter]*origin5 +
                                                    origin6_year_param_array_SRW[paste0(from, "_", to), j, iter]*origin6)/
            sum(exp(b0_array_SRW[from,possible_movements,iter] +
                      # btemp0_array_SRW[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                      btemp1_array_SRW[from,possible_movements,iter]*med_temp[j] + 
                      bspillwindow_array_SRW[from,possible_movements,iter]*med_spillwindow[j] +
                      bwinterspill_array_SRW[from,possible_movements,iter]*med_winterspill[j] +
                      # btemp0xorigin1_array_SRW[from,possible_movements,iter]*origin1 +
                      btemp1xorigin1_array_SRW[from,possible_movements,iter]*med_temp[j]*origin1 +
                      # btemp0xorigin2_array_SRW[from,possible_movements,iter]*origin2 + 
                      btemp1xorigin2_array_SRW[from,possible_movements,iter]*med_temp[j]*origin2 + 
                      # btemp0xorigin3_array_SRW[from,possible_movements,iter]*origin3 +
                      btemp1xorigin3_array_SRW[from,possible_movements,iter]*med_temp[j]*origin3 +
                      # btemp0xorigin4_array_SRW[from,possible_movements,iter]*origin4 +
                      btemp1xorigin4_array_SRW[from,possible_movements,iter]*med_temp[j]*origin4 +
                      # btemp0xorigin5_array_SRW[from,possible_movements,iter]*origin5 +
                      btemp1xorigin5_array_SRW[from,possible_movements,iter]*med_temp[j]*origin5 +
                      # btemp0xorigin6_array_SRW[from,possible_movements,iter]*origin6 +
                      btemp1xorigin6_array_SRW[from,possible_movements,iter]*med_temp[j]*origin6 +
                      borigin1_array_SRW[from,possible_movements,iter]*origin1 +
                      borigin2_array_SRW[from,possible_movements,iter]*origin2 +
                      borigin3_array_SRW[from,possible_movements,iter]*origin3 +
                      borigin4_array_SRW[from,possible_movements,iter]*origin4 +
                      borigin5_array_SRW[from,possible_movements,iter]*origin5 +
                      borigin6_array_SRW[from,possible_movements,iter]*origin6 +
                      # year effects
                      c(origin1_year_param_array_SRW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin1 +
                      c(origin2_year_param_array_SRW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin2 +
                      c(origin3_year_param_array_SRW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin3 +
                      c(origin4_year_param_array_SRW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin4 +
                      c(origin5_year_param_array_SRW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin5 +
                      c(origin6_year_param_array_SRW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin6))
          
        } else {
          
        }
        
        
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(year_move_prob_array)
  
}

estimate_year_move_prob_SRH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,5)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  origin3 = hatchery_origin_params[3]
  origin4 = hatchery_origin_params[4]
  origin5 = hatchery_origin_params[5]
  
  # set up years to predict across
  year_predict <- seq(1,20)
  
  # get an array to store probabilities of movements in different years
  niter <- 4000 # this is the number of draws we have
  
  # set this up where rownames = years, columns = iter, slices = movements
  year_move_prob_array <- array(dim = c(length(year_predict), niter, nrow(movements)))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    SRH_states_dates_years <- data.frame(state = as.vector(SRH_envir$data$y),
                                         date = as.vector(SRH_envir$data$transition_dates))
    
    SRH_states_dates_years$year <- ceiling(SRH_states_dates_years$date/365.25)+1
    
    # get median spill window for this state, by year - this will be indexed the same as year, using j
    spillwindow_data <- SRH_envir$data$spill_window_data
    
    med_spillwindow <- vector(length = 20)
    
    # for spill and temp, if there were no fish in year/state combination,
    # there will be NAs in the vector - and that's ok
    
    for (x in 1:length(med_spillwindow)){
      med_spillwindow[x] <- median(spillwindow_data[subset(SRH_states_dates_years, state == from & year == x)$date,from])
    }
    
    # get median winter spill for this state (note that this is just one value, since each year already only has one winter spill value), 
    # by year - this will be indexed the same as year, using j
    winterspill_data <- SRH_envir$data$winter_spill_days_data
    med_winterspill <- vector(length = 20)
    
    for (y in 1:length(med_winterspill)){
      med_winterspill[y] <- winterspill_data[y,from]
    }
    
    # get median temperature for this state, by year - this will be indexed the same as year, using j
    temp_data <- SRH_envir$data$temperature_data
    
    med_temp <- vector(length = 20)
    
    for (z in 1:length(med_temp)){
      med_temp[z] <- median(temp_data[subset(SRH_states_dates_years, state == from & year == z)$date,from])
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(year_predict)){
        
        # if there aren't any movements in that year, skip it
        if (!(is.na(med_temp[j]))){
          # evaluation movement 
          year_move_prob_array[j, iter, i] <- exp(b0_array_SRH[from,to,iter] +
                                                    # btemp0_array_SRH[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                    btemp1_array_SRH[from,to,iter]*med_temp[j] + 
                                                    bspillwindow_array_SRH[from,to,iter]*med_spillwindow[j] +
                                                    bwinterspill_array_SRH[from,to,iter]*med_winterspill[j] +
                                                    # btemp0xorigin1_array_SRH[from,to,iter]*origin1 +
                                                    btemp1xorigin1_array_SRH[from,to,iter]*med_temp[j]*origin1 +
                                                    # btemp0xorigin2_array_SRH[from,to,iter]*origin2 + 
                                                    btemp1xorigin2_array_SRH[from,to,iter]*med_temp[j]*origin2 + 
                                                    # btemp0xorigin3_array_SRH[from,to,iter]*origin3 +
                                                    btemp1xorigin3_array_SRH[from,to,iter]*med_temp[j]*origin3 +
                                                    # btemp0xorigin4_array_SRH[from,to,iter]*origin4 +
                                                    btemp1xorigin4_array_SRH[from,to,iter]*med_temp[j]*origin4 +
                                                    # btemp0xorigin5_array_SRH[from,to,iter]*origin5 +
                                                    btemp1xorigin5_array_SRH[from,to,iter]*med_temp[j]*origin5 +
                                                    borigin1_array_SRH[from,to,iter]*origin1 +
                                                    borigin2_array_SRH[from,to,iter]*origin2 +
                                                    borigin3_array_SRH[from,to,iter]*origin3 +
                                                    borigin4_array_SRH[from,to,iter]*origin4 +
                                                    borigin5_array_SRH[from,to,iter]*origin5 +
                                                    
                                                    # year effects
                                                    origin1_year_param_array_SRH[paste0(from, "_", to), j, iter]*origin1 +
                                                    origin2_year_param_array_SRH[paste0(from, "_", to), j, iter]*origin2 +
                                                    origin3_year_param_array_SRH[paste0(from, "_", to), j, iter]*origin3 +
                                                    origin4_year_param_array_SRH[paste0(from, "_", to), j, iter]*origin4 +
                                                    origin5_year_param_array_SRH[paste0(from, "_", to), j, iter]*origin5)/
            sum(exp(b0_array_SRH[from,possible_movements,iter] +
                      # btemp0_array_SRH[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                      btemp1_array_SRH[from,possible_movements,iter]*med_temp[j] + 
                      bspillwindow_array_SRH[from,possible_movements,iter]*med_spillwindow[j] +
                      bwinterspill_array_SRH[from,possible_movements,iter]*med_winterspill[j] +
                      # btemp0xorigin1_array_SRH[from,possible_movements,iter]*origin1 +
                      btemp1xorigin1_array_SRH[from,possible_movements,iter]*med_temp[j]*origin1 +
                      # btemp0xorigin2_array_SRH[from,possible_movements,iter]*origin2 + 
                      btemp1xorigin2_array_SRH[from,possible_movements,iter]*med_temp[j]*origin2 + 
                      # btemp0xorigin3_array_SRH[from,possible_movements,iter]*origin3 +
                      btemp1xorigin3_array_SRH[from,possible_movements,iter]*med_temp[j]*origin3 +
                      # btemp0xorigin4_array_SRH[from,possible_movements,iter]*origin4 +
                      btemp1xorigin4_array_SRH[from,possible_movements,iter]*med_temp[j]*origin4 +
                      # btemp0xorigin5_array_SRH[from,possible_movements,iter]*origin5 +
                      btemp1xorigin5_array_SRH[from,possible_movements,iter]*med_temp[j]*origin5 +
                      borigin1_array_SRH[from,possible_movements,iter]*origin1 +
                      borigin2_array_SRH[from,possible_movements,iter]*origin2 +
                      borigin3_array_SRH[from,possible_movements,iter]*origin3 +
                      borigin4_array_SRH[from,possible_movements,iter]*origin4 +
                      borigin5_array_SRH[from,possible_movements,iter]*origin5 +
                      # year effects
                      c(origin1_year_param_array_SRH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin1 +
                      c(origin2_year_param_array_SRH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin2 +
                      c(origin3_year_param_array_SRH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin3 +
                      c(origin4_year_param_array_SRH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin4 +
                      c(origin5_year_param_array_SRH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin5))
          
        } else {
          
        }
        
        
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(year_move_prob_array)
  
}


plot_RE_year <- function(wild_year_param_array = NULL,
                         hatchery_year_param_array = NULL,
                         from, to, plot_title = NULL){
  rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")
  year_predict <- seq(1,20)
  
  if(is.null(hatchery_year_param_array)){
    wild_RE_df <- as.data.frame(wild_year_param_array[paste0(from, "_", to), ,])
    
    niter <- 4000 # for the number of draws
    
    colnames(wild_RE_df) <- paste0("iter", 1:niter) 
    wild_year_predict <- 1:20
    wild_RE_df$year <- year_predict
    
    # Add a column with the actual temperatures
    wild_RE_df$year_actual <- 2004:2023
    
    # drop years without observations
    wild_RE_df <- na.omit(wild_RE_df)
    
    
    wild_RE_df %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(year_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) -> wild_RE_df_quantiles
    
    wild_RE_df_quantiles %>% 
      mutate(rear = "wild") -> RE_df_quantiles
    
  } else if(is.null(wild_year_param_array)){
    hatchery_RE_df <- as.data.frame(hatchery_year_param_array[paste0(from, "_", to), ,])
    
    niter <- 4000 # for the number of draws
    
    colnames(hatchery_RE_df) <- paste0("iter", 1:niter) 
    hatchery_year_predict <- 1:20
    hatchery_RE_df$year <- year_predict
    
    # Add a column with the actual temperatures
    hatchery_RE_df$year_actual <- 2004:2023
    
    # drop years without observations
    hatchery_RE_df <- na.omit(hatchery_RE_df)
    
    
    hatchery_RE_df %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(year_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) -> hatchery_RE_df_quantiles
    
    hatchery_RE_df_quantiles %>% 
      mutate(rear = "hatchery") -> RE_df_quantiles
    
  } else {
    wild_RE_df <- as.data.frame(wild_year_param_array[paste0(from, "_", to), ,])
    
    niter <- 4000 # for the number of draws
    
    colnames(wild_RE_df) <- paste0("iter", 1:niter) 
    wild_year_predict <- 1:20
    wild_RE_df$year <- year_predict
    
    # Add a column with the actual temperatures
    wild_RE_df$year_actual <- 2004:2023
    
    # drop years without observations
    wild_RE_df <- na.omit(wild_RE_df)
    
    
    wild_RE_df %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(year_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) -> wild_RE_df_quantiles
    
    wild_RE_df_quantiles %>% 
      mutate(rear = "wild") -> wild_RE_df_quantiles
    
    hatchery_RE_df <- as.data.frame(hatchery_year_param_array[paste0(from, "_", to), ,])
    
    niter <- 4000 # for the number of draws
    
    colnames(hatchery_RE_df) <- paste0("iter", 1:niter) 
    hatchery_year_predict <- 1:20
    hatchery_RE_df$year <- year_predict
    
    # Add a column with the actual temperatures
    hatchery_RE_df$year_actual <- 2004:2023
    
    # drop years without observations
    hatchery_RE_df <- na.omit(hatchery_RE_df)
    
    
    hatchery_RE_df %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(year_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) -> hatchery_RE_df_quantiles
    
    hatchery_RE_df_quantiles %>% 
      mutate(rear = "hatchery") -> hatchery_RE_df_quantiles
    
    wild_RE_df_quantiles %>% 
      bind_rows(., hatchery_RE_df_quantiles) -> RE_df_quantiles
    
  }
  
  
  
  # start here!
  RE_year_plot <- ggplot(RE_df_quantiles, aes(x = year_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`,
                                              color = rear, fill = rear)) +
    geom_line() +
    geom_ribbon(alpha = 0.2, color = NA) +
    scale_color_manual(values = rear_colors) +
    scale_fill_manual(values = rear_colors) +
    # scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
    scale_x_continuous(lim = c(2003, 2024), breaks = seq(2005, 2020, 5), expand = c(0,0)) +
    xlab("Run Year") +
    ylab("Parameter estimate") +
    ggtitle(plot_title)
  
  return(RE_year_plot)
}



#### Extract year parameter values from the model fit objects ####

# Here it's a bit weird with NDE vs DE - because those change by year
# We're only going to look at DE, because we trust those more than NDE

# function to store all year effects in an array
make_year_origin_parameter_draws_array <- function(fit, fit_summary, origin_select_numeric, envir){
  
  parameters <- fit_summary$variable
  # Figure out which year effects we're actually estimating
  parameters[grepl(paste0("raw_vector"), parameters)] -> year_raw_params
  year_raw_params_from = as.numeric(sub("[^_]*_[^_]*_[^_]*_", "", str_extract(year_raw_params, "[^_]*_[^_]*_[^_]*_[^_]*")))
  year_raw_params_to = as.numeric(str_extract(sub("[^_]*_[^_]*_[^_]*_[^_]*_", "", year_raw_params), "\\d+"))
  year_raw_params_indices <- data.frame(parameter = year_raw_params, from = year_raw_params_from, to = year_raw_params_to)
  year_raw_params_indices %>% 
    dplyr::select(-parameter) %>% 
    mutate(from_to = paste0(from, "_", to)) %>% 
    filter(!(duplicated(from_to)))-> year_raw_params_indices
  
  
  # get the parameter indexing
  param_indices_matrix <- envir$data$parameter_indices_matrix
  param_indices <- data.frame(row = vector(length = max(param_indices_matrix)-1),
                              col = vector(length = max(param_indices_matrix)-1))
  for (i in 1:(max(param_indices_matrix)-1)){
    param_indices[i,] <- which(param_indices_matrix == i, arr.ind = TRUE)
  }
  param_indices$index <- as.numeric(rownames(param_indices))
  param_indices$from_to <- paste0(param_indices$row, "_", param_indices$col)
  
  # Now, include only the indices that actually have year effects
  year_movement_indices <- filter(param_indices, from_to %in% year_raw_params_indices$from_to)
  year_movement_indices %>% 
    dplyr::rename(from = row, to = col) -> year_movement_indices
  # Loop through this (?) to create the param array below
  
  
  
  
  # extract year effect as an 4-d array with rows=from, columns=to, years = slices, and iter=4th dimension
  parameters <- fit_summary$variable
  parameters[grepl(paste0("actual"), parameters)] -> year_params
  # drop the NDE parameters
  year_params <- year_params[!(grepl("_NDE", year_params))]
  # select only one origin
  year_origin_params <- year_params[grepl(origin_select_numeric, year_params)]
  
  
  year_origin_params_index <- as.numeric(str_split_i(str_extract(year_origin_params, "(?<=\\[).*(?=\\])"), "[,]", 1))
  year_origin_params_year <-  as.numeric(str_split_i(str_extract(year_origin_params, "(?<=\\[).*(?=\\])"), "[,]", 2))
  
  year_origin_params_indices <- data.frame(parameter = year_origin_params, index = year_origin_params_index, year = year_origin_params_year)
  year_origin_params_indices %>% 
    left_join(., year_movement_indices, by = "index") %>% 
    filter(!(is.na(from))) -> year_origin_params_indices
  
  # get a repeating index for indexing the 3d array below
  year_origin_params_indices$array_index <- rep(1:length(unique(year_origin_params_indices$index)), length(unique(year_origin_params_indices$year)))
  
  
  # update: extract year effect as 3-d array, with rows = movements, columns = years, iter = slices
  # include movements that don't have a year effect, to make our indexing easier later on
  # these will all have zero (for no year effect)
  year_origin_param_array <- array(data = 0, dim = c(nrow(param_indices),
                                                     length(unique(year_origin_params_indices$year)),
                                                     length(as.matrix(fit[,,1]))))
  
  # give them names that indicate movements
  rownames(year_origin_param_array) <- param_indices$from_to
  
  # now populate the array with the iteration draws
  for(i in 1:nrow(year_origin_params_indices)){
    year_origin_param_array[year_origin_params_indices[i, "index"], year_origin_params_indices[i, "year"], ] <- as.matrix(fit[,,year_origin_params_indices[i, "parameter"]])
  }
  
  
  
  
  return(year_origin_param_array)
}



### UCW ###
origin1_year_param_array_UCW <- make_year_origin_parameter_draws_array(fit = UCW_fit, fit_summary = UCW_fit_summary, 
                                                                       origin_select_numeric = "origin1", envir = UCW_envir)

origin2_year_param_array_UCW <- make_year_origin_parameter_draws_array(fit = UCW_fit, fit_summary = UCW_fit_summary, 
                                                                       origin_select_numeric = "origin2", envir = UCW_envir)

origin3_year_param_array_UCW <- make_year_origin_parameter_draws_array(fit = UCW_fit, fit_summary = UCW_fit_summary, 
                                                                       origin_select_numeric = "origin3", envir = UCW_envir)

### UCH ###
origin1_year_param_array_UCH <- make_year_origin_parameter_draws_array(fit = UCH_fit, fit_summary = UCH_fit_summary, 
                                                                       origin_select_numeric = "origin1", envir = UCH_envir)

origin2_year_param_array_UCH <- make_year_origin_parameter_draws_array(fit = UCH_fit, fit_summary = UCH_fit_summary, 
                                                                       origin_select_numeric = "origin2", envir = UCH_envir)

origin3_year_param_array_UCH <- make_year_origin_parameter_draws_array(fit = UCH_fit, fit_summary = UCH_fit_summary, 
                                                                       origin_select_numeric = "origin3", envir = UCH_envir)

### MCW ###
origin1_year_param_array_MCW <- make_year_origin_parameter_draws_array(fit = MCW_fit, fit_summary = MCW_fit_summary, 
                                                                       origin_select_numeric = "origin1", envir = MCW_envir)

origin2_year_param_array_MCW <- make_year_origin_parameter_draws_array(fit = MCW_fit, fit_summary = MCW_fit_summary, 
                                                                       origin_select_numeric = "origin2", envir = MCW_envir)

origin3_year_param_array_MCW <- make_year_origin_parameter_draws_array(fit = MCW_fit, fit_summary = MCW_fit_summary, 
                                                                       origin_select_numeric = "origin3", envir = MCW_envir)

origin4_year_param_array_MCW <- make_year_origin_parameter_draws_array(fit = MCW_fit, fit_summary = MCW_fit_summary, 
                                                                       origin_select_numeric = "origin4", envir = MCW_envir)

origin5_year_param_array_MCW <- make_year_origin_parameter_draws_array(fit = MCW_fit, fit_summary = MCW_fit_summary, 
                                                                       origin_select_numeric = "origin5", envir = MCW_envir)

origin6_year_param_array_MCW <- make_year_origin_parameter_draws_array(fit = MCW_fit, fit_summary = MCW_fit_summary, 
                                                                       origin_select_numeric = "origin6", envir = MCW_envir)

### MCH ###
origin1_year_param_array_MCH <- make_year_origin_parameter_draws_array(fit = MCH_fit, fit_summary = MCH_fit_summary, 
                                                                       origin_select_numeric = "origin1", envir = MCH_envir)

origin2_year_param_array_MCH <- make_year_origin_parameter_draws_array(fit = MCH_fit, fit_summary = MCH_fit_summary, 
                                                                       origin_select_numeric = "origin2", envir = MCH_envir)

### SRW ###
origin1_year_param_array_SRW <- make_year_origin_parameter_draws_array(fit = SRW_fit, fit_summary = SRW_fit_summary, 
                                                                       origin_select_numeric = "origin1", envir = SRW_envir)

origin2_year_param_array_SRW <- make_year_origin_parameter_draws_array(fit = SRW_fit, fit_summary = SRW_fit_summary, 
                                                                       origin_select_numeric = "origin2", envir = SRW_envir)

origin3_year_param_array_SRW <- make_year_origin_parameter_draws_array(fit = SRW_fit, fit_summary = SRW_fit_summary, 
                                                                       origin_select_numeric = "origin3", envir = SRW_envir)

origin4_year_param_array_SRW <- make_year_origin_parameter_draws_array(fit = SRW_fit, fit_summary = SRW_fit_summary, 
                                                                       origin_select_numeric = "origin4", envir = SRW_envir)

origin5_year_param_array_SRW <- make_year_origin_parameter_draws_array(fit = SRW_fit, fit_summary = SRW_fit_summary, 
                                                                       origin_select_numeric = "origin5", envir = SRW_envir)

origin6_year_param_array_SRW <- make_year_origin_parameter_draws_array(fit = SRW_fit, fit_summary = SRW_fit_summary, 
                                                                       origin_select_numeric = "origin6", envir = SRW_envir)

### SRH ###
origin1_year_param_array_SRH <- make_year_origin_parameter_draws_array(fit = SRH_fit, fit_summary = SRH_fit_summary, 
                                                                       origin_select_numeric = "origin1", envir = SRH_envir)

origin2_year_param_array_SRH <- make_year_origin_parameter_draws_array(fit = SRH_fit, fit_summary = SRH_fit_summary, 
                                                                       origin_select_numeric = "origin2", envir = SRH_envir)

origin3_year_param_array_SRH <- make_year_origin_parameter_draws_array(fit = SRH_fit, fit_summary = SRH_fit_summary, 
                                                                       origin_select_numeric = "origin3", envir = SRH_envir)

origin4_year_param_array_SRH <- make_year_origin_parameter_draws_array(fit = SRH_fit, fit_summary = SRH_fit_summary, 
                                                                       origin_select_numeric = "origin4", envir = SRH_envir)

origin5_year_param_array_SRH <- make_year_origin_parameter_draws_array(fit = SRH_fit, fit_summary = SRH_fit_summary, 
                                                                       origin_select_numeric = "origin5", envir = SRH_envir)



#### Year effect for 9 -> 8 movement ####

## Natural origin

fallback_LGR_ASO_wild_RE_year_plot <- plot_RE_year(wild_year_param_array = origin1_year_param_array_SRW, from = 9, to = 8, plot_title = "Year effect for wild Asotin: fallback at LGR")
ggsave(here::here("figures", "year_effects", "fallback_LGR_ASO_wild_RE_year_plot.png"), fallback_LGR_ASO_wild_RE_year_plot, height = 5, width = 8)

fallback_LGR_CLE_wild_RE_year_plot <- plot_RE_year(wild_year_param_array = origin2_year_param_array_SRW, from = 9, to = 8, plot_title = "Year effect for wild Clearwater: fallback at LGR")
ggsave(here::here("figures", "year_effects", "fallback_LGR_CLE_wild_RE_year_plot.png"), fallback_LGR_CLE_wild_RE_year_plot, height = 5, width = 8)

fallback_LGR_GR_wild_RE_year_plot <- plot_RE_year(wild_year_param_array = origin3_year_param_array_SRW, from = 9, to = 8, plot_title = "Year effect for wild Grande Ronde: fallback at LGR")
ggsave(here::here("figures", "year_effects", "fallback_LGR_GR_wild_RE_year_plot.png"), fallback_LGR_GR_wild_RE_year_plot, height = 5, width = 8)

fallback_LGR_IMN_wild_RE_year_plot <- plot_RE_year(wild_year_param_array = origin4_year_param_array_SRW, from = 9, to = 8, plot_title = "Year effect for wild Imnaha: fallback at LGR")
ggsave(here::here("figures", "year_effects", "fallback_LGR_IMN_wild_RE_year_plot.png"), fallback_LGR_IMN_wild_RE_year_plot, height = 5, width = 8)

fallback_LGR_SAL_wild_RE_year_plot <- plot_RE_year(wild_year_param_array = origin5_year_param_array_SRW, from = 9, to = 8, plot_title = "Year effect for wild Salmon: fallback at LGR")
ggsave(here::here("figures", "year_effects", "fallback_LGR_SAL_wild_RE_year_plot.png"), fallback_LGR_SAL_wild_RE_year_plot, height = 5, width = 8)

## Hatchery origin

fallback_LGR_CLE_hatchery_RE_year_plot <- plot_RE_year(hatchery_year_param_array = origin1_year_param_array_SRH, from = 9, to = 8, plot_title = "Year effect for hatchery Clearwater: fallback at LGR")
ggsave(here::here("figures", "year_effects", "fallback_LGR_CLE_hatchery_RE_year_plot.png"), fallback_LGR_CLE_hatchery_RE_year_plot, height = 5, width = 8)

fallback_LGR_GR_hatchery_RE_year_plot <- plot_RE_year(hatchery_year_param_array = origin2_year_param_array_SRH, from = 9, to = 8, plot_title = "Year effect for hatchery Grande Ronde: fallback at LGR")
ggsave(here::here("figures", "year_effects", "fallback_LGR_GR_hatchery_RE_year_plot.png"), fallback_LGR_GR_hatchery_RE_year_plot, height = 5, width = 8)

fallback_LGR_IMN_hatchery_RE_year_plot <- plot_RE_year(hatchery_year_param_array = origin3_year_param_array_SRH, from = 9, to = 8, plot_title = "Year effect for hatchery Imnaha: fallback at LGR")
ggsave(here::here("figures", "year_effects", "fallback_LGR_IMN_hatchery_RE_year_plot.png"), fallback_LGR_IMN_hatchery_RE_year_plot, height = 5, width = 8)

fallback_LGR_SAL_hatchery_RE_year_plot <- plot_RE_year(hatchery_year_param_array = origin4_year_param_array_SRH, from = 9, to = 8, plot_title = "Year effect for hatchery Salmon: fallback at LGR")
ggsave(here::here("figures", "year_effects", "fallback_LGR_SAL_hatchery_RE_year_plot.png"), fallback_LGR_SAL_hatchery_RE_year_plot, height = 5, width = 8)


#### Inspect data ####

plot(x = 2005:2024, y = SRW_envir$data$winter_spill_days_data[,"LGR"])
# much higher spill in the last two years, and the last four years
# have had occupied 4/5 top spill days.
# so perhaps there is some confounding here with the spillway detector having
# come online around the same time that there is more spill to help fish.

#### Plot fit to data ####

estimate_spilldays_effect_SRW <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- rep(0,6)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = wild_origin_params[1]
  origin2 = wild_origin_params[2]
  origin3 = wild_origin_params[3]
  origin4 = wild_origin_params[4]
  origin5 = wild_origin_params[5]
  origin6 = wild_origin_params[6]
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  spill_move_prob_array <- array(dim = c(100, niter, nrow(movements)))
  dimnames(spill_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  # spill_move_prob_array <- array(dim = c(length(model_states), length(model_states), 100, niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # set up spill window values to predict across
    # extract the data
    winterspill_data <- SRW_envir$data$winter_spill_days_data
    # create a sequence of values from 1 to max
    spill_days_predict <- seq(0, max(winterspill_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    SRW_states_dates <- data.frame(state = as.vector(SRW_envir$data$y),
                                   date_numeric = as.vector(SRW_envir$data$transition_dates))
    
    # add month to states dates so that we can subset
    SRW_states_dates %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> SRW_states_dates
    
    # get median temperature for this state - only in Jan/Feb/Mar
    temp_data <- as.data.frame(SRW_envir$data$temperature)
    temp_data %>% 
      rownames_to_column("date_numeric") %>% 
      pivot_longer(cols = -c(date_numeric)) -> temp_data_long
    
    temp_data_long %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> temp_data_long
    
    # only keep actual observed fish
    med_temp <- median(temp_data[subset(SRW_states_dates, state == from & month %in% c(1,2,3))$date_numeric,from])
    
    # do not include spill window, because code is set up to index so that it's only one or the other
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_days_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_SRW[movements$from[i],movements$to[i],iter] +
                                                 btemp0_array_SRW[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # btemp1_array_SRW[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # bspillwindow_array_SRW[movements$from[i],movements$to[i],iter]*med_spillwindow + # no spill window if spill days is indexed
                                                 bwinterspill_array_SRW[movements$from[i],movements$to[i],iter]*spill_days_predict[j] +
                                                 btemp0xorigin1_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp1xorigin1_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 btemp0xorigin2_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin2 +
                                                 # btemp1xorigin2_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 btemp0xorigin3_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 # btemp1xorigin3_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 btemp0xorigin4_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 # btemp1xorigin4_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 btemp0xorigin5_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 # btemp1xorigin5_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 btemp0xorigin6_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin6 +
                                                 # btemp1xorigin6_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin6 +
                                                 borigin1_array_SRW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_SRW[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_SRW[movements$from[i],movements$to[i],iter]*origin3 +
                                                 borigin4_array_SRW[movements$from[i],movements$to[i],iter]*origin4 +
                                                 borigin5_array_SRW[movements$from[i],movements$to[i],iter]*origin5 +
                                                 borigin6_array_SRW[movements$from[i],movements$to[i],iter]*origin6)/
          sum(exp(b0_array_SRW[movements$from[i],possible_movements,iter] +
                    btemp0_array_SRW[movements$from[i],possible_movements,iter]*med_temp + # select temp0 for any winterspill estimates
                    # btemp1_array_SRW[movements$from[i],possible_movements,iter]*med_temp + 
                    # bspillwindow_array_SRW[movements$from[i],possible_movements,iter]*med_spillwindow + # no spill window if spill days is indexed
                    bwinterspill_array_SRW[movements$from[i],possible_movements,iter]*spill_days_predict[j] +
                    btemp0xorigin1_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp1xorigin1_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    btemp0xorigin2_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin2 +
                    # btemp1xorigin2_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    btemp0xorigin3_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    # btemp1xorigin3_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    btemp0xorigin4_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    # btemp1xorigin4_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    btemp0xorigin5_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin5 +
                    # btemp1xorigin5_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin5 +
                    btemp0xorigin6_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin6 +
                    # btemp1xorigin6_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin6 +
                    borigin1_array_SRW[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_SRW[movements$from[i],possible_movements,iter]*origin2 +
                    borigin3_array_SRW[movements$from[i],possible_movements,iter]*origin3 +
                    borigin4_array_SRW[movements$from[i],possible_movements,iter]*origin4 +
                    borigin5_array_SRW[movements$from[i],possible_movements,iter]*origin5 +
                    borigin6_array_SRW[movements$from[i],possible_movements,iter]*origin6))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}

estimate_spilldays_effect_SRH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,5)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  origin3 = hatchery_origin_params[3]
  origin4 = hatchery_origin_params[4]
  origin5 = hatchery_origin_params[5]
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  spill_move_prob_array <- array(dim = c(100, niter, nrow(movements)))
  dimnames(spill_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  # spill_move_prob_array <- array(dim = c(length(model_states), length(model_states), 100, niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # set up spill window values to predict across
    # extract the data
    winterspill_data <- SRH_envir$data$winter_spill_days_data
    # create a sequence of values from 1 to max
    spill_days_predict <- seq(0, max(winterspill_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    SRH_states_dates <- data.frame(state = as.vector(SRH_envir$data$y),
                                   date_numeric = as.vector(SRH_envir$data$transition_dates))
    
    # add month to states dates so that we can subset
    SRH_states_dates %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> SRH_states_dates
    
    # get median temperature for this state - only in Jan/Feb/Mar
    temp_data <- as.data.frame(SRH_envir$data$temperature)
    temp_data %>% 
      rownames_to_column("date_numeric") %>% 
      pivot_longer(cols = -c(date_numeric)) -> temp_data_long
    
    temp_data_long %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> temp_data_long
    
    # only keep actual observed fish
    med_temp <- median(temp_data[subset(SRH_states_dates, state == from & month %in% c(1,2,3))$date_numeric,from])
    
    # do not include spill window, because code is set up to index so that it's only one or the other
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_days_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_SRH[movements$from[i],movements$to[i],iter] +
                                                 btemp0_array_SRH[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # btemp1_array_SRH[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # bspillwindow_array_SRH[movements$from[i],movements$to[i],iter]*med_spillwindow + # no spill window if spill days is indexed
                                                 bwinterspill_array_SRH[movements$from[i],movements$to[i],iter]*spill_days_predict[j] +
                                                 btemp0xorigin1_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp1xorigin1_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 btemp0xorigin2_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin2 +
                                                 # btemp1xorigin2_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 btemp0xorigin3_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 # btemp1xorigin3_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 btemp0xorigin4_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 # btemp1xorigin4_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 btemp0xorigin5_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 # btemp1xorigin5_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 borigin1_array_SRH[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_SRH[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_SRH[movements$from[i],movements$to[i],iter]*origin3 +
                                                 borigin4_array_SRH[movements$from[i],movements$to[i],iter]*origin4 +
                                                 borigin5_array_SRH[movements$from[i],movements$to[i],iter]*origin5)/
          sum(exp(b0_array_SRH[movements$from[i],possible_movements,iter] +
                    btemp0_array_SRH[movements$from[i],possible_movements,iter]*med_temp + # select temp0 for any winterspill estimates
                    # btemp1_array_SRH[movements$from[i],possible_movements,iter]*med_temp + 
                    # bspillwindow_array_SRH[movements$from[i],possible_movements,iter]*med_spillwindow + # no spill window if spill days is indexed
                    bwinterspill_array_SRH[movements$from[i],possible_movements,iter]*spill_days_predict[j] +
                    btemp0xorigin1_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp1xorigin1_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    btemp0xorigin2_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin2 +
                    # btemp1xorigin2_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    btemp0xorigin3_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    # btemp1xorigin3_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    btemp0xorigin4_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    # btemp1xorigin4_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    btemp0xorigin5_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin5 +
                    # btemp1xorigin5_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin5 +
                    borigin1_array_SRH[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_SRH[movements$from[i],possible_movements,iter]*origin2 +
                    borigin3_array_SRH[movements$from[i],possible_movements,iter]*origin3 +
                    borigin4_array_SRH[movements$from[i],possible_movements,iter]*origin4 +
                    borigin5_array_SRH[movements$from[i],possible_movements,iter]*origin5))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}


plot_compare_rear_spilldays_effect_fit_to_data <- function(origin_select,
                                                           wild_move_prob_array = NULL, hatchery_move_prob_array = NULL,
                                                           wild_covariate_experiences = NULL, hatchery_covariate_experiences = NULL,
                                                           movements_evaluated, spill_predict,
                                                           from, to, plot_title = NULL){
  array_index <- which(movements_evaluated$from == from & movements_evaluated$to == to)
  
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    # modify covariate experiences to show what the next state is
    wild_covariate_experiences %>% 
      mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> wild_covariate_experiences
    
    wild_spill_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(wild_spill_move_prob) <- paste0("iter", 1:niter) 
    wild_spill_move_prob$spill <- spill_predict
    
    wild_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_spill_move_prob_quantiles
    
    wild_spill_move_prob_quantiles -> rear_spill_move_prob_quantiles
    
    # convert back to days
    rear_spill_move_prob_quantiles %>% 
      mutate(spill_actual = spill*100) -> rear_spill_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    
    wild_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
  } else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    # modify covariate experiences to show what the next state is
    hatchery_covariate_experiences %>% 
      mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> hatchery_covariate_experiences
    
    # If wild is NA, run hatchery only
    hatchery_spill_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(hatchery_spill_move_prob) <- paste0("iter", 1:niter) 
    hatchery_spill_move_prob$spill <- spill_predict
    
    hatchery_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_spill_move_prob_quantiles
    
    hatchery_spill_move_prob_quantiles -> rear_spill_move_prob_quantiles
    
    # convert back to days
    rear_spill_move_prob_quantiles %>% 
      mutate(spill_actual = spill*100) -> rear_spill_move_prob_quantiles
    
    # get data organized for rug plot
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    
    hatchery_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from)  -> covariate_experiences
    
  } else {
    # else run both
    
    # modify covariate experiences to show what the next state is
    wild_covariate_experiences %>% 
      mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> wild_covariate_experiences
    hatchery_covariate_experiences %>% 
      mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> hatchery_covariate_experiences
    
    wild_spill_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(wild_spill_move_prob) <- paste0("iter", 1:niter) 
    wild_spill_move_prob$spill <- spill_predict
    
    wild_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_spill_move_prob_quantiles
    
    hatchery_spill_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(hatchery_spill_move_prob) <- paste0("iter", 1:niter) 
    hatchery_spill_move_prob$spill <- spill_predict
    
    hatchery_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_spill_move_prob_quantiles
    
    wild_spill_move_prob_quantiles %>% 
      bind_rows(., hatchery_spill_move_prob_quantiles) -> rear_spill_move_prob_quantiles
    
    # convert back to days
    rear_spill_move_prob_quantiles %>% 
      mutate(spill_actual = spill*100) -> rear_spill_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    wild_covariate_experiences %>% 
      bind_rows(., hatchery_covariate_experiences) %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
    
    
  }
  
  # create a spill actual column
  covariate_experiences %>% 
    mutate(spill_actual = winter_spill*100) -> covariate_experiences
  
  # Keep only fish where they could have experienced spill conditions for rug plot
  # this will keep jan/feb/march, unless it's the march only model
  covariate_experiences %>% 
    filter(winter_post_overshoot_vector == 1) -> covariate_experiences
  
  # now, create DF to show fit to data
  # run this only if the fish actually visited that state
  if (nrow(covariate_experiences) > 0){
    
    # create a new DF to show bubbles of points
    covariate_experiences %>% 
      group_by(spill_actual, year, rear) %>% 
      summarise(total = n()) -> spill_rear_total_counts
    
    covariate_experiences %>% 
      group_by(spill_actual, year, rear, next_state) %>% 
      summarise(count = n()) %>% 
      ungroup() %>% 
      complete(spill_actual, nesting(rear, next_state), fill = list(count = 0)) -> spill_rear_move_counts
    
    spill_rear_move_counts %>% 
      left_join(spill_rear_total_counts, by = c("spill_actual", "year", "rear")) %>% 
      mutate(prop = count/total) -> spill_move_props
    
    # get the actual years
    years_forjoin <- data.frame(year = 1:20, year_actual = 2005:2024)
    
    spill_move_props %>% 
      left_join(., years_forjoin, by = "year") -> spill_move_props
    
  } else {
    spill_move_props <- data.frame(spill_actual = as.numeric(NA), rear = as.character(NA), 
                                   next_state = as.character(NA), count  = as.numeric(NA), 
                                   total = as.numeric(NA), prop = as.numeric(NA))
  }
  
  rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")
  
  # Calculate the credible interval coverage
  # also calculate the average bias (from model estimated median)
  rear_spill_move_prob_quantiles %>% 
    left_join(subset(spill_move_props, next_state == to), by = c("spill_actual", "rear")) %>% 
    # round the model estimate prop to the nearest 0.001 - otherwise we are getting points that are
    # falling outside of CI because points are zero and model estimates a probability that is basically
    # zero but just above
    mutate(CI_coverage = ifelse(prop >= round(`0.025`, 3) & prop <= round(`0.975`,3), 1, 0)) %>% 
    mutate(CI_coverage_weight = CI_coverage*total) %>% 
    # calculate bias
    mutate(bias = `0.5` - prop) %>% 
    mutate(bias_weight = bias*total) %>% 
    mutate(expected_count = `0.5` * total) -> rear_spill_move_prob_quantiles_CIcov
  
  rear_spill_move_prob_quantiles_CIcov %>% 
    group_by(rear) %>% 
    summarise(points_within_CI = sum(CI_coverage_weight, na.rm = T),
              total_bias = sum(bias_weight, na.rm = T),
              total = sum(total, na.rm = T),
              actual_movements = sum(count, na.rm = T),
              total_expected_count = sum(expected_count, na.rm = T)) %>% 
    mutate(CI_coverage = points_within_CI/total) %>% 
    mutate(overall_bias = total_bias/total) %>% 
    mutate(model_v_data = total_expected_count/actual_movements) %>% 
    # add information on origin and movement
    mutate(origin = origin_select, from = from, to = to) %>% 
    relocate(origin) %>% 
    relocate(from, .after = rear) %>% 
    relocate(to, .after = from) -> goodness_of_fit
  
  rear_spill_move_prob_plot <- ggplot(rear_spill_move_prob_quantiles, aes(x = spill_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                          color = rear, fill = rear)) +
    geom_line() +
    geom_ribbon(alpha = 0.2, color = NA) +
    geom_point(data = subset(spill_move_props, next_state == to & rear == "wild"),
               aes(x = spill_actual, y = prop, color = rear, size = total), alpha = 0.5,
               inherit.aes = FALSE) +
    geom_point(data = subset(spill_move_props, next_state == to & rear == "hatchery"),
               aes(x = spill_actual, y = prop, color = rear, size = total), alpha = 0.5,
               inherit.aes = FALSE) +
    geom_text(data = subset(spill_move_props, next_state == to & rear == "wild"),
               aes(x = spill_actual, y = prop, color = rear, label = year_actual),
               inherit.aes = FALSE, size = 2) +
    geom_text(data = subset(spill_move_props, next_state == to & rear == "hatchery"),
               aes(x = spill_actual, y = prop, color = rear, label = year_actual),
               inherit.aes = FALSE, size = 2) +
    scale_y_continuous(lim = c(-0.01, 1.01), expand = c(0,0)) +
    scale_x_continuous(lim = c(-1,max(rear_spill_move_prob_quantiles$spill_actual)+1), expand = c(0,0)) +
    scale_color_manual(values = rear_colors) +
    scale_fill_manual(values = rear_colors) +
    xlab("Days of winter spill") +
    ylab("Movement probability") +
    ggtitle(plot_title)
  
  return(list(plot = rear_spill_move_prob_plot,
              goodness_of_fit = goodness_of_fit))
}


SRH_spill_days_data <- SRH_envir$data$winter_spill_days_data
SRW_spill_days_data <- SRW_envir$data$winter_spill_days_data
SR_spill_days_data <- bind_rows(as.data.frame(SRH_spill_days_data), as.data.frame(SRW_spill_days_data))

# Tucannon River - evaluate all fallback move probs
TUC_wild_spilldays_move_prob_array <- estimate_spilldays_effect_SRW(origin_select= "Tucannon River", movements = SR_postovershoot_fallback_movements)
TUC_hatchery_spilldays_move_prob_array <- estimate_spilldays_effect_SRH(origin_select= "Tucannon River", movements = SR_postovershoot_fallback_movements)


# Tucannon River - get covariate experiences
TUC_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Tucannon River")
TUC_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Tucannon River")

# Tucannon River - use a loop to plot all of them
SR_postovershoot_fallback_movements <- data.frame(dam = c("LGR"),
                                                  from = c(9), to = c(8))
TUC_fallback_spilldays_plots <- vector(mode = "list", length = nrow(SR_postovershoot_fallback_movements))
for (i in 1:nrow(SR_postovershoot_fallback_movements)){
  # if a movement doesn't have any data, don't plot it
  filter(TUC_wild_covariate_experiences, state == SR_postovershoot_fallback_movements$from[i]) -> TUC_wild_data_to_plot
  filter(TUC_hatchery_covariate_experiences, state == SR_postovershoot_fallback_movements$from[i]) -> TUC_hatchery_data_to_plot
  
  if(nrow(TUC_wild_data_to_plot)>0 | nrow(TUC_hatchery_data_to_plot)>0){
    TUC_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect_fit_to_data(origin_select = "Tucannon River",
                                                                                        wild_move_prob_array = TUC_wild_spilldays_move_prob_array,
                                                                                        hatchery_move_prob_array = TUC_hatchery_spilldays_move_prob_array,
                                                                                        wild_covariate_experiences = TUC_wild_covariate_experiences, 
                                                                                        hatchery_covariate_experiences = TUC_hatchery_covariate_experiences,
                                                                                        movements_evaluated = SR_postovershoot_fallback_movements,
                                                                                        spill_predict = seq(0, max(SR_spill_days_data[,SR_postovershoot_fallback_movements$from[i]]),length = 100),
                                                                                        from = SR_postovershoot_fallback_movements$from[i], to = SR_postovershoot_fallback_movements$to[i], 
                                                                                        plot_title = paste0("Tucannon River - post-overshoot fallback at ", SR_postovershoot_fallback_movements$dam[i]))
    
    ggsave(here::here("figures", "spilldays_effects", paste0("TUC_compare_fallback_", 
                                                                                  SR_postovershoot_fallback_movements$dam[i],
                                                                                  "_spilldays.png")),
           TUC_fallback_spilldays_plots[[i]]$plot, height = 8, width = 8)
  } else {
    # if it is zero, do nothing
    
  }
}

# Okay, so this plot shows that 2024 has the most spill days

# can we just plot fallback events at TUC? do we see incidence of observations changing over time? pathway would be informative
# did this locally



