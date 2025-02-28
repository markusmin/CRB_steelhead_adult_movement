# 07-08_Fig7_spillvolume_effects_plot

# This script takes the output from the stan model runs and
# plots the effects of spill volume

# First, need to load in all of the model runs and all of the packages.
source("R/07_model_analysis/07-01_load_stan_models.R")


#### Define functions to estimate effect of spill window on functions of interest ####

# Tell it the movements for which you want to estimate spill volume effects
# movements are formatted as matrix, with column for from and column for to
estimate_spillwindow_effect_UCW <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- rep(0,3)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = wild_origin_params[1]
  origin2 = wild_origin_params[2]
  origin3 = wild_origin_params[3]
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  spill_move_prob_array <- array(dim = c(100, niter, nrow(movements)))
  dimnames(spill_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  # spill_move_prob_array <- array(dim = c(length(model_states), length(model_states), 100, niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # set up spill values to predict across
    # extract the data
    spill_window_data <- UCW_envir$data$spill_window_data
    # create a sequence of values from 1 to max
    spill_predict <- seq(0, max(spill_window_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    UCW_states_dates <- data.frame(state = as.vector(UCW_envir$data$y),
                                   date = as.vector(UCW_envir$data$transition_dates))
    
    # get median temperature for this state
    temp_data <- UCW_envir$data$temperature
    med_temp <- median(temp_data[subset(UCW_states_dates, state == from)$date,from])
    
    # get median winter spill for this state
    winterspill_data <- UCW_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of spill volume values to get predicted response
      # spill volume was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_UCW[movements$from[i],movements$to[i],iter] +
                                                 # btemp0_array_UCW[movements$from[i],movements$to[i],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                 btemp1_array_UCW[movements$from[i],movements$to[i],iter]*med_temp + 
                                                 bspillwindow_array_UCW[movements$from[i],movements$to[i],iter]*spill_predict[j] +
                                                 # bwinterspill_array_UCW[movements$from[i],movements$to[i],iter]*med_winterspill + # indexing is set up so that if spill window is selected, spill days is ignored
                                                 # btemp0xorigin1_array_UCW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 btemp1xorigin1_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp0xorigin2_array_UCW[movements$from[i],movements$to[i],iter]*origin2 + 
                                                 btemp1xorigin2_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 # btemp0xorigin3_array_UCW[movements$from[i],movements$to[i],iter]*origin3 +
                                                 btemp1xorigin3_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 borigin1_array_UCW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_UCW[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_UCW[movements$from[i],movements$to[i],iter]*origin3)/
          sum(exp(b0_array_UCW[movements$from[i],possible_movements,iter] +
                    # btemp0_array_UCW[movements$from[i],possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_UCW[movements$from[i],possible_movements,iter]*med_temp + 
                    bspillwindow_array_UCW[movements$from[i],possible_movements,iter]*spill_predict[j] +
                    # bwinterspill_array_UCW[movements$from[i],possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_UCW[movements$from[i],possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp0xorigin2_array_UCW[movements$from[i],possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    # btemp0xorigin3_array_UCW[movements$from[i],possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    borigin1_array_UCW[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_UCW[movements$from[i],possible_movements,iter]*origin2 +
                    borigin3_array_UCW[movements$from[i],possible_movements,iter]*origin3))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}

estimate_spillwindow_effect_UCH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,3)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  origin3 = hatchery_origin_params[3]
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  spill_move_prob_array <- array(dim = c(100, niter, nrow(movements)))
  dimnames(spill_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  # spill_move_prob_array <- array(dim = c(length(model_states), length(model_states), 100, niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # set up spill values to predict across
    # extract the data
    spill_window_data <- UCH_envir$data$spill_window_data
    # create a sequence of values from 1 to max
    spill_predict <- seq(0, max(spill_window_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    UCH_states_dates <- data.frame(state = as.vector(UCH_envir$data$y),
                                   date = as.vector(UCH_envir$data$transition_dates))
    
    # get median temperature for this state
    temp_data <- UCH_envir$data$temperature
    med_temp <- median(temp_data[subset(UCH_states_dates, state == from)$date,from])
    
    # get median winter spill for this state
    winterspill_data <- UCH_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of spill volume values to get predicted response
      # spill volume was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_UCH[movements$from[i],movements$to[i],iter] +
                                                 # btemp0_array_UCH[movements$from[i],movements$to[i],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                 btemp1_array_UCH[movements$from[i],movements$to[i],iter]*med_temp + 
                                                 bspillwindow_array_UCH[movements$from[i],movements$to[i],iter]*spill_predict[j] +
                                                 # bwinterspill_array_UCH[movements$from[i],movements$to[i],iter]*med_winterspill + # indexing is set up so that if spill window is selected, spill days is ignored
                                                 # btemp0xorigin1_array_UCH[movements$from[i],movements$to[i],iter]*origin1 +
                                                 btemp1xorigin1_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp0xorigin2_array_UCH[movements$from[i],movements$to[i],iter]*origin2 + 
                                                 btemp1xorigin2_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 # btemp0xorigin3_array_UCH[movements$from[i],movements$to[i],iter]*origin3 +
                                                 btemp1xorigin3_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 borigin1_array_UCH[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_UCH[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_UCH[movements$from[i],movements$to[i],iter]*origin3)/
          sum(exp(b0_array_UCH[movements$from[i],possible_movements,iter] +
                    # btemp0_array_UCH[movements$from[i],possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_UCH[movements$from[i],possible_movements,iter]*med_temp + 
                    bspillwindow_array_UCH[movements$from[i],possible_movements,iter]*spill_predict[j] +
                    # bwinterspill_array_UCH[movements$from[i],possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_UCH[movements$from[i],possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp0xorigin2_array_UCH[movements$from[i],possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    # btemp0xorigin3_array_UCH[movements$from[i],possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    borigin1_array_UCH[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_UCH[movements$from[i],possible_movements,iter]*origin2 +
                    borigin3_array_UCH[movements$from[i],possible_movements,iter]*origin3))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}
estimate_spillwindow_effect_MCW <- function(origin_select, movements){
  
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
    
    # set up spill values to predict across
    # extract the data
    spill_window_data <- MCW_envir$data$spill_window_data
    # create a sequence of values from 1 to max
    spill_predict <- seq(0, max(spill_window_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    MCW_states_dates <- data.frame(state = as.vector(MCW_envir$data$y),
                                   date = as.vector(MCW_envir$data$transition_dates))
    
    # get median temperature for this state
    temp_data <- MCW_envir$data$temperature
    med_temp <- median(temp_data[subset(MCW_states_dates, state == from)$date,from])
    
    # get median winter spill for this state
    winterspill_data <- MCW_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of spill volume values to get predicted response
      # spill volume was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_MCW[movements$from[i],movements$to[i],iter] +
                                                 # btemp0_array_MCW[movements$from[i],movements$to[i],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                 btemp1_array_MCW[movements$from[i],movements$to[i],iter]*med_temp + 
                                                 bspillwindow_array_MCW[movements$from[i],movements$to[i],iter]*spill_predict[j] +
                                                 # bwinterspill_array_MCW[movements$from[i],movements$to[i],iter]*med_winterspill + # indexing is set up so that if spill window is selected, spill days is ignored
                                                 # btemp0xorigin1_array_MCW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 btemp1xorigin1_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp0xorigin2_array_MCW[movements$from[i],movements$to[i],iter]*origin2 + 
                                                 btemp1xorigin2_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 # btemp0xorigin3_array_MCW[movements$from[i],movements$to[i],iter]*origin3 +
                                                 btemp1xorigin3_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 # btemp0xorigin4_array_MCW[movements$from[i],movements$to[i],iter]*origin4 +
                                                 btemp1xorigin4_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 # btemp0xorigin5_array_MCW[movements$from[i],movements$to[i],iter]*origin5 +
                                                 btemp1xorigin5_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 # btemp0xorigin6_array_MCW[movements$from[i],movements$to[i],iter]*origin6 +
                                                 btemp1xorigin6_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin6 +
                                                 borigin1_array_MCW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_MCW[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_MCW[movements$from[i],movements$to[i],iter]*origin3 +
                                                 borigin4_array_MCW[movements$from[i],movements$to[i],iter]*origin4 +
                                                 borigin5_array_MCW[movements$from[i],movements$to[i],iter]*origin5 +
                                                 borigin6_array_MCW[movements$from[i],movements$to[i],iter]*origin6)/
          sum(exp(b0_array_MCW[movements$from[i],possible_movements,iter] +
                    # btemp0_array_MCW[movements$from[i],possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_MCW[movements$from[i],possible_movements,iter]*med_temp + 
                    bspillwindow_array_MCW[movements$from[i],possible_movements,iter]*spill_predict[j] +
                    # bwinterspill_array_MCW[movements$from[i],possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_MCW[movements$from[i],possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp0xorigin2_array_MCW[movements$from[i],possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    # btemp0xorigin3_array_MCW[movements$from[i],possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    # btemp0xorigin4_array_MCW[movements$from[i],possible_movements,iter]*origin4 +
                    btemp1xorigin4_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    # btemp0xorigin5_array_MCW[movements$from[i],possible_movements,iter]*origin5 +
                    btemp1xorigin5_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin5 +
                    # btemp0xorigin6_array_MCW[movements$from[i],possible_movements,iter]*origin6 +
                    btemp1xorigin6_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin6 +
                    borigin1_array_MCW[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_MCW[movements$from[i],possible_movements,iter]*origin2 +
                    borigin3_array_MCW[movements$from[i],possible_movements,iter]*origin3 +
                    borigin4_array_MCW[movements$from[i],possible_movements,iter]*origin4 +
                    borigin5_array_MCW[movements$from[i],possible_movements,iter]*origin5 +
                    borigin6_array_MCW[movements$from[i],possible_movements,iter]*origin6))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}

estimate_spillwindow_effect_MCH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,2)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  spill_move_prob_array <- array(dim = c(100, niter, nrow(movements)))
  dimnames(spill_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  # spill_move_prob_array <- array(dim = c(length(model_states), length(model_states), 100, niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # set up spill values to predict across
    # extract the data
    spill_window_data <- MCH_envir$data$spill_window_data
    # create a sequence of values from 1 to max
    spill_predict <- seq(0, max(spill_window_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    MCH_states_dates <- data.frame(state = as.vector(MCH_envir$data$y),
                                   date = as.vector(MCH_envir$data$transition_dates))
    
    # get median temperature for this state
    temp_data <- MCH_envir$data$temperature
    med_temp <- median(temp_data[subset(MCH_states_dates, state == from)$date,from])
    
    # get median winter spill for this state
    winterspill_data <- MCH_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of spill volume values to get predicted response
      # spill volume was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_MCH[movements$from[i],movements$to[i],iter] +
                                                 # btemp0_array_MCH[movements$from[i],movements$to[i],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                 btemp1_array_MCH[movements$from[i],movements$to[i],iter]*med_temp + 
                                                 bspillwindow_array_MCH[movements$from[i],movements$to[i],iter]*spill_predict[j] +
                                                 # bwinterspill_array_MCH[movements$from[i],movements$to[i],iter]*med_winterspill + # indexing is set up so that if spill window is selected, spill days is ignored
                                                 # btemp0xorigin1_array_MCH[movements$from[i],movements$to[i],iter]*origin1 +
                                                 btemp1xorigin1_array_MCH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp0xorigin2_array_MCH[movements$from[i],movements$to[i],iter]*origin2 + 
                                                 btemp1xorigin2_array_MCH[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 borigin1_array_MCH[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_MCH[movements$from[i],movements$to[i],iter]*origin2)/
          sum(exp(b0_array_MCH[movements$from[i],possible_movements,iter] +
                    # btemp0_array_MCH[movements$from[i],possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_MCH[movements$from[i],possible_movements,iter]*med_temp + 
                    bspillwindow_array_MCH[movements$from[i],possible_movements,iter]*spill_predict[j] +
                    # bwinterspill_array_MCH[movements$from[i],possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_MCH[movements$from[i],possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_MCH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp0xorigin2_array_MCH[movements$from[i],possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_MCH[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    borigin1_array_MCH[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_MCH[movements$from[i],possible_movements,iter]*origin2))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}

estimate_spillwindow_effect_SRW <- function(origin_select, movements){
  
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
    
    # set up spill values to predict across
    # extract the data
    spill_window_data <- SRW_envir$data$spill_window_data
    # create a sequence of values from 1 to max
    spill_predict <- seq(0, max(spill_window_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    SRW_states_dates <- data.frame(state = as.vector(SRW_envir$data$y),
                                   date = as.vector(SRW_envir$data$transition_dates))
    
    # get median temperature for this state
    temp_data <- SRW_envir$data$temperature
    med_temp <- median(temp_data[subset(SRW_states_dates, state == from)$date,from])
    
    # get median winter spill for this state
    winterspill_data <- SRW_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of spill volume values to get predicted response
      # spill volume was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_SRW[movements$from[i],movements$to[i],iter] +
                                                 # btemp0_array_SRW[movements$from[i],movements$to[i],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                 btemp1_array_SRW[movements$from[i],movements$to[i],iter]*med_temp + 
                                                 bspillwindow_array_SRW[movements$from[i],movements$to[i],iter]*spill_predict[j] +
                                                 # bwinterspill_array_SRW[movements$from[i],movements$to[i],iter]*med_winterspill + # indexing is set up so that if spill window is selected, spill days is ignored
                                                 # btemp0xorigin1_array_SRW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 btemp1xorigin1_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp0xorigin2_array_SRW[movements$from[i],movements$to[i],iter]*origin2 + 
                                                 btemp1xorigin2_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 # btemp0xorigin3_array_SRW[movements$from[i],movements$to[i],iter]*origin3 +
                                                 btemp1xorigin3_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 # btemp0xorigin4_array_SRW[movements$from[i],movements$to[i],iter]*origin4 +
                                                 btemp1xorigin4_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 # btemp0xorigin5_array_SRW[movements$from[i],movements$to[i],iter]*origin5 +
                                                 btemp1xorigin5_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 # btemp0xorigin6_array_SRW[movements$from[i],movements$to[i],iter]*origin6 +
                                                 btemp1xorigin6_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin6 +
                                                 borigin1_array_SRW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_SRW[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_SRW[movements$from[i],movements$to[i],iter]*origin3 +
                                                 borigin4_array_SRW[movements$from[i],movements$to[i],iter]*origin4 +
                                                 borigin5_array_SRW[movements$from[i],movements$to[i],iter]*origin5 +
                                                 borigin6_array_SRW[movements$from[i],movements$to[i],iter]*origin6)/
          sum(exp(b0_array_SRW[movements$from[i],possible_movements,iter] +
                    # btemp0_array_SRW[movements$from[i],possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_SRW[movements$from[i],possible_movements,iter]*med_temp + 
                    bspillwindow_array_SRW[movements$from[i],possible_movements,iter]*spill_predict[j] +
                    # bwinterspill_array_SRW[movements$from[i],possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_SRW[movements$from[i],possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp0xorigin2_array_SRW[movements$from[i],possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    # btemp0xorigin3_array_SRW[movements$from[i],possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    # btemp0xorigin4_array_SRW[movements$from[i],possible_movements,iter]*origin4 +
                    btemp1xorigin4_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    # btemp0xorigin5_array_SRW[movements$from[i],possible_movements,iter]*origin5 +
                    btemp1xorigin5_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin5 +
                    # btemp0xorigin6_array_SRW[movements$from[i],possible_movements,iter]*origin6 +
                    btemp1xorigin6_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin6 +
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

estimate_spillwindow_effect_SRH <- function(origin_select, movements){
  
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
    
    # set up spill values to predict across
    # extract the data
    spill_window_data <- SRH_envir$data$spill_window_data
    # create a sequence of values from 1 to max
    spill_predict <- seq(0, max(spill_window_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    SRH_states_dates <- data.frame(state = as.vector(SRH_envir$data$y),
                                   date = as.vector(SRH_envir$data$transition_dates))
    
    # get median temperature for this state
    temp_data <- SRH_envir$data$temperature
    med_temp <- median(temp_data[subset(SRH_states_dates, state == from)$date,from])
    
    # get median winter spill for this state
    winterspill_data <- SRH_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of spill volume values to get predicted response
      # spill volume was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_SRH[movements$from[i],movements$to[i],iter] +
                                                 # btemp0_array_SRH[movements$from[i],movements$to[i],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                 btemp1_array_SRH[movements$from[i],movements$to[i],iter]*med_temp + 
                                                 bspillwindow_array_SRH[movements$from[i],movements$to[i],iter]*spill_predict[j] +
                                                 # bwinterspill_array_SRH[movements$from[i],movements$to[i],iter]*med_winterspill + # indexing is set up so that if spill window is selected, spill days is ignored
                                                 # btemp0xorigin1_array_SRH[movements$from[i],movements$to[i],iter]*origin1 +
                                                 btemp1xorigin1_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp0xorigin2_array_SRH[movements$from[i],movements$to[i],iter]*origin2 + 
                                                 btemp1xorigin2_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 # btemp0xorigin3_array_SRH[movements$from[i],movements$to[i],iter]*origin3 +
                                                 btemp1xorigin3_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 # btemp0xorigin4_array_SRH[movements$from[i],movements$to[i],iter]*origin4 +
                                                 btemp1xorigin4_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 # btemp0xorigin5_array_SRH[movements$from[i],movements$to[i],iter]*origin5 +
                                                 btemp1xorigin5_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 borigin1_array_SRH[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_SRH[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_SRH[movements$from[i],movements$to[i],iter]*origin3 +
                                                 borigin4_array_SRH[movements$from[i],movements$to[i],iter]*origin4 +
                                                 borigin5_array_SRH[movements$from[i],movements$to[i],iter]*origin5)/
          sum(exp(b0_array_SRH[movements$from[i],possible_movements,iter] +
                    # btemp0_array_SRH[movements$from[i],possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_SRH[movements$from[i],possible_movements,iter]*med_temp + 
                    bspillwindow_array_SRH[movements$from[i],possible_movements,iter]*spill_predict[j] +
                    # bwinterspill_array_SRH[movements$from[i],possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_SRH[movements$from[i],possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp0xorigin2_array_SRH[movements$from[i],possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    # btemp0xorigin3_array_SRH[movements$from[i],possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    # btemp0xorigin4_array_SRH[movements$from[i],possible_movements,iter]*origin4 +
                    btemp1xorigin4_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    # btemp0xorigin5_array_SRH[movements$from[i],possible_movements,iter]*origin5 +
                    btemp1xorigin5_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin5 +
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


#### Estimate the effect of spill window on fallback at each dam for each population ####

## First, create data frames that store all of the movements that we care about
# Here, we need to only list dams that are en-route fallback (downstream
# of natal tributaries)

## Middle Columbia
JDR_spillvolume_movements <- data.frame(from = c(2), to = c(1),
                                        Dam = c("BON"))
DES_spillvolume_movements <- data.frame(from = c(2), to = c(1),
                                        Dam = c("BON"))
UMA_spillvolume_movements <- data.frame(from = c(2), to = c(1),
                                        Dam = c("BON"))
FIF_spillvolume_movements <- data.frame(from = c(2), to = c(1),
                                        Dam = c("BON"))
YAK_spillvolume_movements <- data.frame(from = c(3, 2), to = c(2, 1),
                                        Dam = c("MCN", "BON"))
WAWA_spillvolume_movements <- data.frame(from = c(3, 2), to = c(2, 1),
                                        Dam = c("MCN", "BON"))
middle_columbia_spillvolume_movements <- data.frame(from = c(3, 2), to = c(2, 1),
                                                    Dam = c("MCN", "BON"))

# Upper Columbia
WEN_spillvolume_movements <- data.frame(from = c(5, 4, 3, 2), to = c(4, 3, 2, 1),
                                         Dam = c("RIS", "PRA", "MCN", "BON"))
ENT_spillvolume_movements <- data.frame(from = c(6, 5, 4, 3, 2), to = c(5, 4, 3, 2, 1),
                                        Dam = c("RRE", "RIS", "PRA", "MCN", "BON"))
MET_spillvolume_movements <- data.frame(from = c(7, 6, 5, 4, 3, 2), to = c(6, 5, 4, 3, 2, 1),
                                        Dam = c("WEL", "RRE", "RIS", "PRA", "MCN", "BON"))
OKA_spillvolume_movements <- data.frame(from = c(7, 6, 5, 4, 3, 2), to = c(6, 5, 4, 3, 2, 1),
                                        Dam = c("WEL", "RRE", "RIS", "PRA", "MCN", "BON"))
upper_columbia_spillvolume_movements <- data.frame(from = c(7, 6, 5, 4, 3, 2), to = c(6, 5, 4, 3, 2, 1),
                                        Dam = c("WEL", "RRE", "RIS", "PRA", "MCN", "BON"))

# Snake River
TUC_spillvolume_movements <- data.frame(from = c(8, 3, 2), to = c(3, 2, 1),
                  Dam = c("ICH", "MCN", "BON"))
ASO_spillvolume_movements <- data.frame(from = c(9, 8, 3, 2), to = c(8, 3, 2, 1),
                                        Dam = c("LGR", "ICH", "MCN", "BON"))
CLE_spillvolume_movements <- data.frame(from = c(9, 8, 3, 2), to = c(8, 3, 2, 1),
                                        Dam = c("LGR", "ICH", "MCN", "BON"))
IMN_spillvolume_movements <- data.frame(from = c(9, 8, 3, 2), to = c(8, 3, 2, 1),
                                        Dam = c("LGR", "ICH", "MCN", "BON"))
GR_spillvolume_movements <- data.frame(from = c(9, 8, 3, 2), to = c(8, 3, 2, 1),
                                        Dam = c("LGR", "ICH", "MCN", "BON"))
SAL_spillvolume_movements <- data.frame(from = c(9, 8, 3, 2), to = c(8, 3, 2, 1),
                                        Dam = c("LGR", "ICH", "MCN", "BON"))
snake_river_spillvolume_movements <- data.frame(from = c(9, 8, 3, 2), to = c(8, 3, 2, 1),
                                        Dam = c("LGR", "ICH", "MCN", "BON"))

## Now, evaluate the effect of spill for each of these movements

# Middle Columbia, wild
JDR_wild_spillvolume_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select = "John Day River", 
                                                                movements = JDR_spillvolume_movements)
JDR_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, 
                                                                rear = "wild", origin_select = "John Day River")

DES_wild_spillvolume_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select = "Deschutes River", 
                                                                        movements = DES_spillvolume_movements)
DES_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, 
                                                                rear = "wild", origin_select = "Deschutes River")

UMA_wild_spillvolume_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select = "Umatilla River", 
                                                                        movements = UMA_spillvolume_movements)
UMA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, 
                                                                rear = "wild", origin_select = "Umatilla River")

FIF_wild_spillvolume_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select = "Fifteenmile Creek", 
                                                                        movements = FIF_spillvolume_movements)
FIF_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, 
                                                                rear = "wild", origin_select = "Fifteenmile Creek")

YAK_wild_spillvolume_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select = "Yakima River", 
                                                                        movements = YAK_spillvolume_movements)
YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, 
                                                                rear = "wild", origin_select = "Yakima River")

WAWA_wild_spillvolume_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select = "Walla Walla River", 
                                                                        movements = WAWA_spillvolume_movements)
WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, 
                                                                rear = "wild", origin_select = "Walla Walla River")


# Middle Columbia, hatchery
UMA_hatchery_spillvolume_move_prob_array <- estimate_spillwindow_effect_MCH(origin_select = "Umatilla River", 
                                                                        movements = UMA_spillvolume_movements)
UMA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, 
                                                                rear = "hatchery", origin_select = "Umatilla River")

WAWA_hatchery_spillvolume_move_prob_array <- estimate_spillwindow_effect_MCH(origin_select = "Walla Walla River", 
                                                                        movements = WAWA_spillvolume_movements)
WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, 
                                                                rear = "hatchery", origin_select = "Walla Walla River")

# Upper Columbia, wild
WEN_wild_spillvolume_move_prob_array <- estimate_spillwindow_effect_UCW(origin_select = "Wenatchee River", 
                                                                        movements = WEN_spillvolume_movements)
WEN_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, 
                                                                rear = "wild", origin_select = "Wenatchee River")

ENT_wild_spillvolume_move_prob_array <- estimate_spillwindow_effect_UCW(origin_select = "Entiat River", 
                                                                        movements = ENT_spillvolume_movements)
ENT_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, 
                                                                rear = "wild", origin_select = "Entiat River")

MET_wild_spillvolume_move_prob_array <- estimate_spillwindow_effect_UCW(origin_select = "Methow River", 
                                                                        movements = MET_spillvolume_movements)
MET_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, 
                                                                rear = "wild", origin_select = "Methow River")

# Upper Columbia, hatchery
WEN_hatchery_spillvolume_move_prob_array <- estimate_spillwindow_effect_UCH(origin_select = "Wenatchee River", 
                                                                        movements = WEN_spillvolume_movements)
WEN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, 
                                                                rear = "hatchery", origin_select = "Wenatchee River")

OKA_hatchery_spillvolume_move_prob_array <- estimate_spillwindow_effect_UCH(origin_select = "Okanogan River", 
                                                                            movements = OKA_spillvolume_movements)
OKA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, 
                                                                    rear = "hatchery", origin_select = "Okanogan River")

MET_hatchery_spillvolume_move_prob_array <- estimate_spillwindow_effect_UCH(origin_select = "Methow River", 
                                                                            movements = MET_spillvolume_movements)
MET_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, 
                                                                    rear = "hatchery", origin_select = "Methow River")

WEN_hatchery_spillvolume_move_prob_array <- estimate_spillwindow_effect_UCH(origin_select = "Wenatchee River", 
                                                                        movements = WEN_spillvolume_movements)
WEN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, 
                                                                rear = "hatchery", origin_select = "Wenatchee River")# Snake River, wild

# Snake River, wild
TUC_wild_spillvolume_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select = "Tucannon River", 
                                                                            movements = TUC_spillvolume_movements)
TUC_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, 
                                                                    rear = "wild", origin_select = "Tucannon River")

ASO_wild_spillvolume_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select = "Asotin Creek", 
                                                                        movements = ASO_spillvolume_movements)
ASO_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, 
                                                                rear = "wild", origin_select = "Asotin Creek")

CLE_wild_spillvolume_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select = "Clearwater River", 
                                                                        movements = CLE_spillvolume_movements)
CLE_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, 
                                                                rear = "wild", origin_select = "Clearwater River")

SAL_wild_spillvolume_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select = "Salmon River", 
                                                                        movements = SAL_spillvolume_movements)
SAL_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, 
                                                                rear = "wild", origin_select = "Salmon River")

GR_wild_spillvolume_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select = "Grande Ronde River", 
                                                                        movements = GR_spillvolume_movements)
GR_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, 
                                                                rear = "wild", origin_select = "Grande Ronde River")

IMN_wild_spillvolume_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select = "Imnaha River", 
                                                                        movements = IMN_spillvolume_movements)
IMN_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, 
                                                                rear = "wild", origin_select = "Imnaha River")

# Snake River, hatchery
TUC_hatchery_spillvolume_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select = "Tucannon River", 
                                                                        movements = TUC_spillvolume_movements)
TUC_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, 
                                                                rear = "hatchery", origin_select = "Tucannon River")

CLE_hatchery_spillvolume_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select = "Clearwater River", 
                                                                        movements = CLE_spillvolume_movements)
CLE_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, 
                                                                rear = "hatchery", origin_select = "Clearwater River")

SAL_hatchery_spillvolume_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select = "Salmon River", 
                                                                        movements = SAL_spillvolume_movements)
SAL_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, 
                                                                rear = "hatchery", origin_select = "Salmon River")

GR_hatchery_spillvolume_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select = "Grande Ronde River", 
                                                                       movements = GR_spillvolume_movements)
GR_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, 
                                                               rear = "hatchery", origin_select = "Grande Ronde River")

IMN_hatchery_spillvolume_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select = "Imnaha River", 
                                                                        movements = IMN_spillvolume_movements)
IMN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, 
                                                                rear = "hatchery", origin_select = "Imnaha River")


#### Reformat spill window effects for plot ####
# First, get the spillwindow data. It doesn't matter which envir you pull from, they all have the same data
spill_window_data <- SRH_envir$data$spill_window_data
# Function to reformat movement prob arrays
reformat_spillwindow_prob_array <- function(move_prob_array, movements, origin_name, rear_type){
  # loop through each of the movements that we evaluated
  for (i in 1:dim(move_prob_array)[3]){
  
      movement_df <- as.data.frame(move_prob_array[,,i])
      colnames(movement_df) <- 1:ncol(movement_df)
      # get the spill data for that state
      
      # create a sequence of values from 1 to max
      spillwindow_predict <- seq(0, max(spill_window_data[,movements$from[i]]),length = 100)
      # convert this to actual spill by multiplying by 100
      spillwindow_actual <- spillwindow_predict * 100
      
      movement_df$spillwindow_actual <- spillwindow_actual
      
      movement_df %>% 
        pivot_longer(!spillwindow_actual, names_to = "iter", values_to = "prob") %>% 
        mutate(origin = origin_name,
               rear = rear_type,
               from = movements$from[i],
               to = movements$to[i]) -> movement_df
    # if it's the first movement, store it as a new array
    if(i == 1){
      all_movements_df <- movement_df
    } 
      # if it's not, bind it to the existing one
      else {
      all_movements_df %>% 
          bind_rows(movement_df) -> all_movements_df
    }
    # if it's not, add it to the existing
  }
  return(all_movements_df)
}


# reformat spill window movement prob estimates for each population
# Spill window parameters are shared across all origins within a DPS/rear type combination
# BUT - the movement probabilities won't be shared across all populations, because 
# the other parameters are often origin-specific.

### Middle Columbia, wild ###
JDR_natural_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = JDR_wild_spillvolume_move_prob_array,
                                                                                 movements = JDR_spillvolume_movements,
                                                                                 origin_name = "John Day River",
                                                                                 rear_type = "wild")

DES_natural_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = DES_wild_spillvolume_move_prob_array,
                                                                     movements = DES_spillvolume_movements,
                                                                     origin_name = "Deschutes River",
                                                                     rear_type = "wild")

FIF_natural_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = FIF_wild_spillvolume_move_prob_array,
                                                                     movements = FIF_spillvolume_movements,
                                                                     origin_name = "Fifteenmile Creek",
                                                                     rear_type = "wild")

UMA_natural_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = UMA_wild_spillvolume_move_prob_array,
                                                                     movements = UMA_spillvolume_movements,
                                                                     origin_name = "Umatilla River",
                                                                     rear_type = "wild")

YAK_natural_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = YAK_wild_spillvolume_move_prob_array,
                                                                     movements = YAK_spillvolume_movements,
                                                                     origin_name = "Yakima River",
                                                                     rear_type = "wild")

WAWA_natural_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = WAWA_wild_spillvolume_move_prob_array,
                                                                     movements = WAWA_spillvolume_movements,
                                                                     origin_name = "Walla Walla River",
                                                                     rear_type = "wild")





### Middle Columbia, hatchery ###






IMN_hatchery_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = IMN_hatchery_spillvolume_move_prob_array,
                                                                      movements = IMN_spillvolume_movements,
                                                                      origin_name = "Imnaha River",
                                                                      rear_type = "hatchery")


### combine all populations into one df for plotting ###

JDR_natural_spillvolume_move_long %>% 
  bind_rows(., UMA_natural_spillvolume_move_long, YAK_natural_spillvolume_move_long, WAWA_natural_spillvolume_move_long, 
            FIF_natural_spillvolume_move_long, DES_natural_spillvolume_move_long) -> combined_MCW_spillvolume_move_long

# keep only BON (for now)
combined_MCW_spillvolume_move_long %>% 
  filter(from == 2) -> combined_MCW_spillvolume_move_long

# summarise into median and 95% CI
combined_MCW_spillvolume_move_long %>% 
  group_by(spillwindow_actual, origin, rear) %>% 
  summarise(prob = quantile(prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
  pivot_wider(names_from = q, values_from = prob) %>% 
  mutate(origin_rear = paste0(origin, ", ", rear)) -> combined_MCW_spillvolume_move_summary

#### Generate plot ####

MCW_spillvolume_move_prob_plot <- ggplot(combined_MCW_spillvolume_move_summary, aes(x = spillwindow_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                 color = origin_rear, fill = origin_rear)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, color = NA) +
  # scale_y_continuous(lim = c(0,1)) +
  # common x-axis scale across all populations
  coord_cartesian(xlim = c(0, 250), expand = FALSE) +
  # scale_x_continuous(lim = c(floor(min(covariate_experiences$temp_actual)),ceiling(max(covariate_experiences$temp_actual))), expand = c(0,0)) +
  # scale_color_manual(values = movement_colors) +
  # scale_fill_manual(values =  movement_colors) +
  xlab("Spill volume (kcfs)") +
  ylab("Movement probability") +
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.4),
        # turn off the axis titles on each individual plot and just show one for whole plot
        axis.title = element_blank(),
        # these plot margins are to leave space for the population name on the big figure
        plot.margin = unit(c(0, 0.2, 0.2, 0.2),"cm"))


ggsave(here::here("figures", "paper_figures", "Fig7_MCW_spillwindow_movements.png"), MCW_spillvolume_move_prob_plot, height = 8, width = 10)







