# 07-08_Fig7_spillvolume_effects_plot

# This script takes the output from the stan model runs and
# plots the effects of spill volume

# First, need to load in all of the model runs and all of the packages.
source("R/07_model_analysis/07-01_load_stan_models.R")


#### Define functions to estimate effect of spill window on functions of interest ####

extract_covariate_experiences <- function(envir, rear, origin_select){
  origin_vector <- vector(length = nrow(envir$data$cat_X_mat))
  for(i in 1:nrow(envir$data$cat_X_mat)){
    origin_vector[i] <- which(envir$data$cat_X_mat[i,]==1)
  }
  
  
  # for spill days - include the winter post-overshoot vector, which contains
  # info on whether they could have experienced winter spill conditions or not
  # add a new fish_ID column, which is not the same as tag code but will allow us to differentiate between fish
  pop_states_dates <- data.frame(fish_ID = rep(1:length(origin_vector), each = ncol(envir$data$y)),
                                 state = as.vector(t(envir$data$y)),
                                 date = as.vector(t(envir$data$transition_dates)),
                                 year = ceiling(as.vector(t(envir$data$transition_dates))/365.25)+1,
                                 origin = rep(origin_vector, each = ncol(envir$data$y)))
  
  
  
  # add mainstem dam for each state
  dam_index <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL", "ICH", "LGR"),
                          state = seq(2,9))
  
  pop_states_dates %>% 
    left_join(., dam_index, by = "state") -> pop_states_dates
  
  
  # reformat covariates so that they can be joined
  as.data.frame(envir$data$spill_window_data) %>% 
    rownames_to_column("date") %>% 
    dplyr::rename(LGR = LGR_quick, WEL = WEL_quick) %>% 
    mutate(date = as.numeric(date)) %>% 
    pivot_longer(cols = -c(date), names_to = "dam", values_to = "spill_window") -> spill_window_long
  
  as.data.frame(envir$data$temperature_data) %>% 
    rownames_to_column("date") %>% 
    dplyr::rename(LGR = LGR_quick, WEL = WEL_quick) %>% 
    mutate(date = as.numeric(date)) %>% 
    pivot_longer(cols = -c(date), names_to = "dam", values_to = "temperature") -> temp_long
  
  as.data.frame(envir$data$winter_spill_days_data) %>% 
    rownames_to_column("year") %>% 
    mutate(year = as.numeric(year)) %>% 
    pivot_longer(cols = -c(year), names_to = "dam", values_to = "winter_spill") -> spill_days_long
  
  
  # add temperature, spill window, winter spill days
  pop_states_dates %>% 
    left_join(., spill_window_long, by = c("date", "dam")) %>% 
    left_join(., temp_long, by = c("date", "dam")) %>% 
    left_join(., spill_days_long, by = c("year", "dam")) -> pop_states_dates_covariates
  
  # drop observations in the loss state and with index 0
  pop_states_dates_covariates %>% 
    filter(!(state %in% c(0,41))) -> pop_states_dates_covariates
  
  # now add winter post-overshoot vector
  pop_states_dates_covariates$winter_post_overshoot_vector <- as.vector(envir$data$winter_post_overshoot_vector)
  
  
  # Now, keep only the origin selected
  if(rear == "wild"){
    origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$wild
  } else {
    origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$hatchery
  }
  
  subset(pop_states_dates_covariates, origin == origin_numeric) -> origin_states_dates_covariates
  
  return(origin_states_dates_covariates)
  
}

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
UMA_hatchery_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = UMA_hatchery_spillvolume_move_prob_array,
                                                                     movements = UMA_spillvolume_movements,
                                                                     origin_name = "Umatilla River",
                                                                     rear_type = "hatchery")

WAWA_hatchery_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = WAWA_hatchery_spillvolume_move_prob_array,
                                                                      movements = WAWA_spillvolume_movements,
                                                                      origin_name = "Walla Walla River",
                                                                      rear_type = "hatchery")

### Upper Columbia, wild ###
WEN_natural_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = WEN_wild_spillvolume_move_prob_array,
                                                                      movements = WEN_spillvolume_movements,
                                                                      origin_name = "Wenatchee River",
                                                                      rear_type = "wild")

ENT_natural_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = ENT_wild_spillvolume_move_prob_array,
                                                                     movements = ENT_spillvolume_movements,
                                                                     origin_name = "Entiat River",
                                                                     rear_type = "wild")

MET_natural_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = MET_wild_spillvolume_move_prob_array,
                                                                     movements = MET_spillvolume_movements,
                                                                     origin_name = "Methow River",
                                                                     rear_type = "wild")

### Upper Columbia, hatchery ###
WEN_hatchery_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = WEN_hatchery_spillvolume_move_prob_array,
                                                                     movements = WEN_spillvolume_movements,
                                                                     origin_name = "Wenatchee River",
                                                                     rear_type = "hatchery")

OKA_hatchery_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = OKA_hatchery_spillvolume_move_prob_array,
                                                                     movements = OKA_spillvolume_movements,
                                                                     origin_name = "Okanogan River",
                                                                     rear_type = "hatchery")

MET_hatchery_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = MET_hatchery_spillvolume_move_prob_array,
                                                                     movements = MET_spillvolume_movements,
                                                                     origin_name = "Methow River",
                                                                     rear_type = "hatchery")

### Snake River, wild ###
TUC_natural_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = TUC_wild_spillvolume_move_prob_array,
                                                                     movements = TUC_spillvolume_movements,
                                                                     origin_name = "Tucannon River",
                                                                     rear_type = "wild")

ASO_natural_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = ASO_wild_spillvolume_move_prob_array,
                                                                     movements = ASO_spillvolume_movements,
                                                                     origin_name = "Asotin Creek",
                                                                     rear_type = "wild")

CLE_natural_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = CLE_wild_spillvolume_move_prob_array,
                                                                     movements = CLE_spillvolume_movements,
                                                                     origin_name = "Clearwater River",
                                                                     rear_type = "wild")

SAL_natural_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = SAL_wild_spillvolume_move_prob_array,
                                                                     movements = SAL_spillvolume_movements,
                                                                     origin_name = "Salmon River",
                                                                     rear_type = "wild")

GR_natural_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = GR_wild_spillvolume_move_prob_array,
                                                                     movements = GR_spillvolume_movements,
                                                                     origin_name = "Grande Ronde River",
                                                                     rear_type = "wild")

IMN_natural_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = IMN_wild_spillvolume_move_prob_array,
                                                                     movements = IMN_spillvolume_movements,
                                                                     origin_name = "Imnaha River",
                                                                     rear_type = "wild")

### Snake River, hatchery
TUC_hatchery_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = TUC_hatchery_spillvolume_move_prob_array,
                                                                     movements = TUC_spillvolume_movements,
                                                                     origin_name = "Tucannon River",
                                                                     rear_type = "hatchery")

CLE_hatchery_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = CLE_hatchery_spillvolume_move_prob_array,
                                                                     movements = CLE_spillvolume_movements,
                                                                     origin_name = "Clearwater River",
                                                                     rear_type = "hatchery")

SAL_hatchery_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = SAL_hatchery_spillvolume_move_prob_array,
                                                                     movements = SAL_spillvolume_movements,
                                                                     origin_name = "Salmon River",
                                                                     rear_type = "hatchery")

GR_hatchery_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = GR_hatchery_spillvolume_move_prob_array,
                                                                    movements = GR_spillvolume_movements,
                                                                    origin_name = "Grande Ronde River",
                                                                    rear_type = "hatchery")

IMN_hatchery_spillvolume_move_long <- reformat_spillwindow_prob_array(move_prob_array = IMN_hatchery_spillvolume_move_prob_array,
                                                                     movements = IMN_spillvolume_movements,
                                                                     origin_name = "Imnaha River",
                                                                     rear_type = "hatchery")


### combine all populations into one df for plotting ###

## combine movement probabilities vs spill actual
JDR_natural_spillvolume_move_long %>% 
  bind_rows(., UMA_natural_spillvolume_move_long, YAK_natural_spillvolume_move_long, 
            WAWA_natural_spillvolume_move_long, 
            FIF_natural_spillvolume_move_long, DES_natural_spillvolume_move_long) -> combined_MCW_spillvolume_move_long

JDR_natural_spillvolume_move_long %>% 
  bind_rows( # MCW
    UMA_natural_spillvolume_move_long, YAK_natural_spillvolume_move_long, 
            WAWA_natural_spillvolume_move_long, FIF_natural_spillvolume_move_long, 
            DES_natural_spillvolume_move_long, 
            # MCH
            UMA_hatchery_spillvolume_move_long, WAWA_hatchery_spillvolume_move_long,
    # UCW
    WEN_natural_spillvolume_move_long, ENT_natural_spillvolume_move_long,
    MET_natural_spillvolume_move_long,
    # UCH
    WEN_hatchery_spillvolume_move_long, OKA_hatchery_spillvolume_move_long,
    MET_hatchery_spillvolume_move_long,
    # SRW
    TUC_natural_spillvolume_move_long, ASO_natural_spillvolume_move_long,
    CLE_natural_spillvolume_move_long, SAL_natural_spillvolume_move_long,
    GR_natural_spillvolume_move_long, IMN_natural_spillvolume_move_long,
    # SRH
    TUC_hatchery_spillvolume_move_long,
    CLE_hatchery_spillvolume_move_long, SAL_hatchery_spillvolume_move_long,
    GR_hatchery_spillvolume_move_long, IMN_hatchery_spillvolume_move_long
    ) -> combined_spillvolume_move_long

combined_spillvolume_move_long %>% 
  group_by(spillwindow_actual, origin, rear, from) %>% 
  summarise(prob = quantile(prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
  pivot_wider(names_from = q, values_from = prob) %>% 
  mutate(origin_rear = paste0(origin, ", ", rear)) -> combined_spillvolume_move_summary

# Note which populations belong to which DPS
populations_DPS <- data.frame(origin_rear = unique(combined_spillvolume_move_summary$origin_rear),
                              DPS = c("Snake River, natural",
                                      "Snake River, hatchery",
                                      "Snake River, natural",
                                      "Middle Columbia, natural",
                                      "Upper Columbia, natural",
                                      "Middle Columbia, natural",
                                      "Snake River, hatchery",
                                      "Snake River, natural",
                                      "Snake River, hatchery",
                                      "Snake River, natural",
                                      "Middle Columbia, natural",
                                      "Upper Columbia, hatchery",
                                      "Upper Columbia, natural",
                                      "Upper Columbia, hatchery",
                                      "Snake River, hatchery",
                                      "Snake River, natural",
                                      "Snake River, hatchery",
                                      "Snake River, natural",
                                      "Middle Columbia, natural",
                                      "Middle Columbia, hatchery",
                                      "Middle Columbia, natural",
                                      "Middle Columbia, hatchery",
                                      "Upper Columbia, hatchery",
                                      "Upper Columbia, natural",
                                      "Middle Columbia, natural"))

# Create the population groups for fallback at Bonneville
populations_DPS %>% 
  mutate(BON_groups = ifelse(grepl("Middle Columbia", DPS), origin_rear, DPS)) -> populations_DPS


combined_spillvolume_move_summary %>% 
  left_join(populations_DPS, by = "origin_rear") -> combined_spillvolume_move_summary



## combine covariate experiences
mutate(JDR_wild_covariate_experiences, origin_rear = "John Day River, wild") %>% 
  bind_rows( # MCW
    mutate(UMA_wild_covariate_experiences, origin_rear = "Umatilla River, wild"),
    mutate(YAK_wild_covariate_experiences, origin_rear = "Yakima River, wild"), 
    mutate(WAWA_wild_covariate_experiences, origin_rear = "Walla Walla River, wild"), 
    mutate(FIF_wild_covariate_experiences, origin_rear = "Fifteenmile Creek, wild"), 
    mutate(DES_wild_covariate_experiences, origin_rear = "Deschutes River, wild"), 
    # MCH
    mutate(UMA_hatchery_covariate_experiences, origin_rear = "Umatilla River, hatchery"), 
    mutate(WAWA_hatchery_covariate_experiences, origin_rear = "Walla Walla River, hatchery"),
    # UCW
    mutate(WEN_wild_covariate_experiences, origin_rear = "Wenatchee River, wild"), 
    mutate(ENT_wild_covariate_experiences, origin_rear = "Entiat River, wild"),
    mutate(MET_wild_covariate_experiences, origin_rear = "Methow River, wild"),
    # UCH
    mutate(WEN_hatchery_covariate_experiences, origin_rear = "Wenatchee River, hatchery"), 
    mutate(OKA_hatchery_covariate_experiences, origin_rear = "Okanogan River, hatchery"),
    mutate(MET_hatchery_covariate_experiences, origin_rear = "Methow River, hatchery"),
    # SRW
    mutate(TUC_wild_covariate_experiences, origin_rear = "Tucannon River, wild"), 
    mutate(ASO_wild_covariate_experiences, origin_rear = "Asotin Creek, wild"),
    mutate(CLE_wild_covariate_experiences, origin_rear = "Clearwater River, wild"), 
    mutate(SAL_wild_covariate_experiences, origin_rear = "Salmon River, wild"),
    mutate(GR_wild_covariate_experiences, origin_rear = "Grande Ronde River, wild"), 
    mutate(IMN_wild_covariate_experiences, origin_rear = "Imnaha River, wild"),
    # SRH
    mutate(TUC_hatchery_covariate_experiences, origin_rear = "Tucannon River, hatchery"),
    mutate(CLE_hatchery_covariate_experiences, origin_rear = "Clearwater River, hatchery"), 
    mutate(SAL_hatchery_covariate_experiences, origin_rear = "Salmon River, hatchery"),
    mutate(GR_hatchery_covariate_experiences, origin_rear = "Grande Ronde River, hatchery"), 
    mutate(IMN_hatchery_covariate_experiences, origin_rear = "Imnaha River, hatchery")
  ) %>% 
  mutate(spill_window_actual = spill_window * 100) -> combined_covariate_experiences

combined_covariate_experiences %>% 
  left_join(populations_DPS, by = "origin_rear") -> combined_covariate_experiences


#### Define function for multipanel plots ####

# function to create a plot, with a histogram, for one population/DPS

plot_spill_fallback_onepop_onedam <- function(from_state, origin_rear_select = NA, DPS_select = NA,
                                              lower_spill_lim, upper_spill_lim){
  
  # can either subset a single origin/rear or a whole DPS
  if(!is.na(origin_rear_select)){
    move_prob_data <- subset(combined_spillvolume_move_summary, from == from_state & origin_rear == origin_rear_select)
    covariate_data <- subset(combined_covariate_experiences, state == from_state & origin_rear == origin_rear_select)
  } else {
    move_prob_data <- subset(combined_spillvolume_move_summary, from == from_state & DPS == DPS_select)
    covariate_data <- subset(combined_covariate_experiences, state == from_state & DPS == DPS_select)
  }
  
  
  spill_move_prob_plot_onepop_onedam <- ggplot(move_prob_data, aes(x = spillwindow_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`)) +
    geom_line() +
    geom_ribbon(alpha = 0.2, color = NA) +
    # scale_y_continuous(lim = c(0,0.23), expand = c(0,0), breaks = seq(0, 0.2, 0.05)) +
    # xlim(0,150) +
    # common x-axis scale across all populations
    coord_cartesian(xlim = c(lower_spill_lim, upper_spill_lim), ylim = c(0,0.23), expand = FALSE) +
    xlab("Winter spill days") +
    ylab("Movement probability") +
    theme(# legend.position = "none",
      panel.grid.major = element_line(color = "gray90"),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=0.4),
      # turn off the axis titles on each individual plot and just show one for whole plot
      axis.title = element_blank(),
      # these plot margins are to leave space for the population name on the big figure
      plot.margin = unit(c(0, 0.3, 0.2, 0.2),"cm"))
  
  density_plot_onepop_onedam <- ggplot(covariate_data, 
                                  aes(spill_window_actual))+
    geom_histogram(aes(y=..count../sum(..count..)), alpha = 0.5, bins = 60, boundary = 0) +
    ylab("Density") +
    # these labels aren't real, they just need to be the same size as the labels from the temp effect plot
    scale_y_continuous(n.breaks = 2, labels = c("0.00", "1.00")) +
    # xlim(0,150) +
    # common x-axis scale across all populations
    coord_cartesian(xlim = c(lower_spill_lim, upper_spill_lim), expand = FALSE) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(color = "white"),
          axis.ticks.x = element_blank(),
          # axis.ticks.length.x=unit(.1, "cm"),
          axis.ticks.length.x=unit(0, "cm"),
          axis.ticks.y = element_line(color = "white"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          # these plot margins are to leave space for the population name on the big figure
          plot.margin = unit(c(1.0, 0.3, 0, 0.2),"cm"))
  
  combined_plot <- ggarrange(density_plot_onepop_onedam, spill_move_prob_plot_onepop_onedam, nrow = 2, ncol = 1,
                                      heights = c(2,6))
  
  return(combined_plot)
}

#### Fallback at BON plot ####

DES_N_BON_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 2, origin_rear_select = "Deschutes River, wild", 
                                  DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

FIF_N_BON_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 2, origin_rear_select = "Fifteenmile Creek, wild", 
                                                            DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

JDR_N_BON_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 2, origin_rear_select = "John Day River, wild", 
                                                            DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

UMA_H_BON_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 2, origin_rear_select = "Umatilla River, hatchery", 
                                                            DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

UMA_N_BON_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 2, origin_rear_select = "Umatilla River, wild", 
                                                            DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

WAWA_H_BON_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 2, origin_rear_select = "Walla Walla River, hatchery", 
                                                            DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

WAWA_N_BON_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 2, origin_rear_select = "Walla Walla River, wild", 
                                                            DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

YAK_N_BON_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 2, origin_rear_select = "Yakima River, wild", 
                                                            DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

Snake_River_N_BON_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 2, origin_rear_select = NA, 
                                                              DPS_select = "Snake River, natural", lower_spill_lim = 0, upper_spill_lim = 150)

Snake_River_H_BON_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 2, origin_rear_select = NA, 
                                                                      DPS_select = "Snake River, hatchery", lower_spill_lim = 0, upper_spill_lim = 150)

Upper_Columbia_N_BON_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 2, origin_rear_select = NA, 
                                                                      DPS_select = "Upper Columbia, natural", lower_spill_lim = 0, upper_spill_lim = 150)

Upper_Columbia_H_BON_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 2, origin_rear_select = NA, 
                                                                      DPS_select = "Upper Columbia, hatchery", lower_spill_lim = 0, upper_spill_lim = 150)



combined_BON_spill_fallback <- ggarrange(DES_BON_spill_fallback, FIF_N_BON_spill_fallback, JDR_N_BON_spill_fallback,
                                         UMA_N_BON_spill_fallback, UMA_H_BON_spill_fallback, WAWA_N_BON_spill_fallback, 
                                         WAWA_H_BON_spill_fallback, YAK_N_BON_spill_fallback, Snake_River_N_BON_spill_fallback,
                                         Snake_River_H_BON_spill_fallback, Upper_Columbia_N_BON_spill_fallback, Upper_Columbia_H_BON_spill_fallback, 
                                         nrow = 3, ncol = 4,
                                         labels = c("(A) DES, Natural", "(B) FIF, Natural",
                                                    "(C) JDR, Natural", "(D) UMA, Natural",
                                                    "(E) UMA, Hatchery", "(F) WAWA, Natural",
                                                    "(G) WAWA, Hatchery", "(H) YAK, Natural",
                                                    "(I) Snake R., Natural", "(J) Snake R., Hatchery",
                                                    "(K) Upper Col., Natural", "(L) Upper Col., Hatchery"),
                                         label.x = 0.05, label.y = 0.925, font.label = list(size = 14, face = "plain"),
                                         hjust = 0, vjust = 0)


combined_BON_spill_fallback <- cowplot::ggdraw(annotate_figure(combined_BON_spill_fallback,
                                                               bottom = textGrob("Daily Spill at Bonneville Dam (kcfs)", gp = gpar(cex = 1.3)),
                                                               left = textGrob("Fallback probability at Bonneville Dam", rot = 90, gp = gpar(cex = 1.3)))) +
  theme(plot.background = element_rect(fill="white", color = NA))

ggsave(here::here("figures", "paper_figures", "spill_volume", "BON_spill_fallback.png"), 
       combined_BON_spill_fallback, height = 8, width = 10)

#### Fallback at MCN plot ####


WAWA_N_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Walla Walla River, wild", 
                                                               DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

WAWA_H_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Walla Walla River, hatchery", 
                                                               DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

YAK_N_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Yakima River, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

WEN_N_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Wenatchee River, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

WEN_H_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Wenatchee River, hatchery", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

ENT_N_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Entiat River, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

OKA_H_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Okanogan River, hatchery", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

MET_N_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Methow River, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

MET_H_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Methow River, hatchery", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

TUC_N_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Tucannon River, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

TUC_H_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Tucannon River, hatchery", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

ASO_N_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Asotin Creek, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

CLE_N_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Clearwater River, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

CLE_H_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Clearwater River, hatchery", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

IMN_N_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Imnaha River, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

IMN_H_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Imnaha River, hatchery", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

GR_N_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Grande Ronde River, wild", 
                                                             DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

GR_H_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Grande Ronde River, hatchery", 
                                                             DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

SAL_N_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Salmon River, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)

SAL_H_MCN_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 3, origin_rear_select = "Salmon River, hatchery", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 150)



combined_MCN_spill_fallback <- ggarrange(WAWA_N_MCN_spill_fallback,
                                         WAWA_H_MCN_spill_fallback,
                                         YAK_N_MCN_spill_fallback,
                                         WEN_N_MCN_spill_fallback,
                                         WEN_H_MCN_spill_fallback,
                                         ENT_N_MCN_spill_fallback,
                                         OKA_H_MCN_spill_fallback,
                                         MET_N_MCN_spill_fallback,
                                         MET_H_MCN_spill_fallback,
                                         TUC_N_MCN_spill_fallback,
                                         TUC_H_MCN_spill_fallback,
                                         ASO_N_MCN_spill_fallback,
                                         CLE_N_MCN_spill_fallback,
                                         CLE_H_MCN_spill_fallback,
                                         IMN_N_MCN_spill_fallback,
                                         IMN_H_MCN_spill_fallback,
                                         GR_N_MCN_spill_fallback,
                                         GR_H_MCN_spill_fallback,
                                         SAL_N_MCN_spill_fallback,
                                         SAL_H_MCN_spill_fallback, 
                                         nrow = 5, ncol = 4,
                                         labels = c("(A) WAWA, Natural", "(B) WAWA, Hatchery",
                                                    "(C) YAK, Natural", "(D) WEN, Natural",
                                                    "(E) WEN, Hatchery", "(F) ENT, Natural",
                                                    "(G) OKA, Hatchery", "(H) MET, Natural",
                                                    "(I) MET, Hatchery", "(J) TUC, Natural",
                                                    "(K) TUC, Hatchery", "(L) ASO, Natural",
                                                    "(M) CLE, Natural", "(N) CLE, Hatchery",
                                                    "(O) IMN, Natural", "(P) IMN, Hatchery",
                                                    "(Q) GR, Natural", "(R) GR, Hatchery",
                                                    "(S) SAL, Natural", "(T) SAL, Hatchery"),
                                         label.x = 0.05, label.y = 0.925, font.label = list(size = 14, face = "plain"),
                                         hjust = 0, vjust = 0)


combined_MCN_spill_fallback <- cowplot::ggdraw(annotate_figure(combined_MCN_spill_fallback,
                                                               bottom = textGrob("Daily Spill at McNary Dam (kcfs)", gp = gpar(cex = 1.3)),
                                                               left = textGrob("Fallback probability at McNary Dam", rot = 90, gp = gpar(cex = 1.3)))) +
  theme(plot.background = element_rect(fill="white", color = NA))

ggsave(here::here("figures", "paper_figures", "spill_volume", "MCN_spill_fallback.png"), 
       combined_MCN_spill_fallback, height = 40/3, width = 10)



#### Fallback at ICH plot ####


TUC_N_ICH_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 8, origin_rear_select = "Tucannon River, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 100)

TUC_H_ICH_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 8, origin_rear_select = "Tucannon River, hatchery", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 100)

ASO_N_ICH_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 8, origin_rear_select = "Asotin Creek, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 100)

CLE_N_ICH_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 8, origin_rear_select = "Clearwater River, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 100)

CLE_H_ICH_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 8, origin_rear_select = "Clearwater River, hatchery", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 100)

IMN_N_ICH_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 8, origin_rear_select = "Imnaha River, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 100)

IMN_H_ICH_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 8, origin_rear_select = "Imnaha River, hatchery", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 100)

GR_N_ICH_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 8, origin_rear_select = "Grande Ronde River, wild", 
                                                             DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 100)

GR_H_ICH_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 8, origin_rear_select = "Grande Ronde River, hatchery", 
                                                             DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 100)

SAL_N_ICH_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 8, origin_rear_select = "Salmon River, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 100)

SAL_H_ICH_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 8, origin_rear_select = "Salmon River, hatchery", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 100)



combined_ICH_spill_fallback <- ggarrange(TUC_N_ICH_spill_fallback,
                                         TUC_H_ICH_spill_fallback,
                                         ASO_N_ICH_spill_fallback,
                                         CLE_N_ICH_spill_fallback,
                                         CLE_H_ICH_spill_fallback,
                                         IMN_N_ICH_spill_fallback,
                                         IMN_H_ICH_spill_fallback,
                                         GR_N_ICH_spill_fallback,
                                         GR_H_ICH_spill_fallback,
                                         SAL_N_ICH_spill_fallback,
                                         SAL_H_ICH_spill_fallback, 
                                         nrow = 3, ncol = 4,
                                         labels = c("(A) TUC, Natural",
                                                    "(B) TUC, Hatchery", "(C) ASO, Natural",
                                                    "(D) CLE, Natural", "(E) CLE, Hatchery",
                                                    "(F) IMN, Natural", "(G) IMN, Hatchery",
                                                    "(H) GR, Natural", "(I) GR, Hatchery",
                                                    "(J) SAL, Natural", "(K) SAL, Hatchery"),
                                         label.x = 0.05, label.y = 0.925, font.label = list(size = 14, face = "plain"),
                                         hjust = 0, vjust = 0)


combined_ICH_spill_fallback <- cowplot::ggdraw(annotate_figure(combined_ICH_spill_fallback,
                                                               bottom = textGrob("Daily Spill at Ice Harbor Dam (kcfs)", gp = gpar(cex = 1.3)),
                                                               left = textGrob("Fallback probability at Ice Harbor Dam", rot = 90, gp = gpar(cex = 1.3)))) +
  theme(plot.background = element_rect(fill="white", color = NA))

ggsave(here::here("figures", "paper_figures", "spill_volume", "ICH_spill_fallback.png"), 
       combined_ICH_spill_fallback, height = 8, width = 10)


#### Fallback at LGR plot ####

ASO_N_LGR_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 9, origin_rear_select = "Asotin Creek, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 60)

CLE_N_LGR_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 9, origin_rear_select = "Clearwater River, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 60)

CLE_H_LGR_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 9, origin_rear_select = "Clearwater River, hatchery", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 60)

IMN_N_LGR_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 9, origin_rear_select = "Imnaha River, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 60)

IMN_H_LGR_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 9, origin_rear_select = "Imnaha River, hatchery", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 60)

GR_N_LGR_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 9, origin_rear_select = "Grande Ronde River, wild", 
                                                             DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 60)

GR_H_LGR_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 9, origin_rear_select = "Grande Ronde River, hatchery", 
                                                             DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 60)

SAL_N_LGR_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 9, origin_rear_select = "Salmon River, wild", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 60)

SAL_H_LGR_spill_fallback <- plot_spill_fallback_onepop_onedam(from_state = 9, origin_rear_select = "Salmon River, hatchery", 
                                                              DPS_select = NA, lower_spill_lim = 0, upper_spill_lim = 60)



combined_LGR_spill_fallback <- ggarrange(ASO_N_LGR_spill_fallback,
                                         CLE_N_LGR_spill_fallback,
                                         CLE_H_LGR_spill_fallback,
                                         IMN_N_LGR_spill_fallback,
                                         IMN_H_LGR_spill_fallback,
                                         GR_N_LGR_spill_fallback,
                                         GR_H_LGR_spill_fallback,
                                         SAL_N_LGR_spill_fallback,
                                         SAL_H_LGR_spill_fallback, 
                                         nrow = 3, ncol = 4,
                                         labels = c("(A) ASO, Natural",
                                                    "(B) CLE, Natural", "(C) CLE, Hatchery",
                                                    "(D) IMN, Natural", "(E) IMN, Hatchery",
                                                    "(F) GR, Natural", "(G) GR, Hatchery",
                                                    "(H) SAL, Natural", "(I) SAL, Hatchery"),
                                         label.x = 0.05, label.y = 0.925, font.label = list(size = 14, face = "plain"),
                                         hjust = 0, vjust = 0)


combined_LGR_spill_fallback <- cowplot::ggdraw(annotate_figure(combined_LGR_spill_fallback,
                                                               bottom = textGrob("Daily Spill at Lower Granite Dam (kcfs)", gp = gpar(cex = 1.3)),
                                                               left = textGrob("Fallback probability at Lower Granite Dam", rot = 90, gp = gpar(cex = 1.3)))) +
  theme(plot.background = element_rect(fill="white", color = NA))

ggsave(here::here("figures", "paper_figures", "spill_volume", "LGR_spill_fallback.png"), 
       combined_LGR_spill_fallback, height = 8, width = 10)


#### OLD CODE BELOW ####
#### Generate BON plot ####

combined_rear_spill_move_prob_plot <- ggplot(subset(combined_spillvolume_move_summary, from == 2), aes(x = spillwindow_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                              color = BON_groups, fill = BON_groups)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, color = NA) +
  # geom_rug(data = covariate_experiences, aes(x = spill_actual), inherit.aes = FALSE,
  #          sides = "t", length = unit(0.3, "cm"), outside = TRUE) +
  scale_y_continuous(lim = c(0,0.23), expand = c(0,0), breaks = seq(0, 0.2, 0.05)) +
  xlim(0,150) +
  # common x-axis scale across all populations
  coord_cartesian(xlim = c(0, 150), expand = FALSE) +
  # scale_color_manual(values = movement_colors) +
  # scale_fill_manual(values =  movement_colors) +
  scale_color_tableau(palette = "Tableau 20") +
  scale_fill_tableau(palette = "Tableau 20") +
  xlab("Winter spill days") +
  ylab("Movement probability") +
  theme(# legend.position = "none",
        panel.grid.major = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.4),
        # turn off the axis titles on each individual plot and just show one for whole plot
        axis.title = element_blank(),
        # these plot margins are to leave space for the population name on the big figure
        plot.margin = unit(c(0, 0.3, 0.2, 0.2),"cm"))

combined_density_plot <- ggplot(data = subset(combined_covariate_experiences, state == 2), 
                                aes(spill_window_actual, fill = BON_groups))+
  # geom_density(alpha = 0.1) +
  # geom_rug(sides = "t", length = unit(0.2, "cm"), outside = FALSE) +
  # geom_histogram(alpha = 0.5, bins = 60) +
  geom_histogram(aes(y=..count../sum(..count..)), alpha = 0.5, bins = 60, boundary = 0) +
  ylab("Density") +
  # scale_x_continuous(lim = c(floor(min(covariate_experiences$temp_actual)),ceiling(max(covariate_experiences$temp_actual))), expand = c(0,0)) +
  # scale_x_continuous(lim = c(0, 23), expand = c(0,0)) +
  # scale_y_continuous(lim = c(0,0.75), expand = c(0,0),
  #                    breaks = c(0, 0.25, 0.50)) +
  # these labels aren't real, they just need to be the same size as the labels from the temp effect plot
  scale_y_continuous(n.breaks = 2, labels = c("0.00", "1.00")) +
  xlim(0,150) +
  # common x-axis scale across all populations
  coord_cartesian(xlim = c(0, 150), expand = FALSE) +
  scale_fill_tableau(palette = "Tableau 20") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "white"),
        axis.ticks.x = element_blank(),
        # axis.ticks.length.x=unit(.1, "cm"),
        axis.ticks.length.x=unit(0, "cm"),
        axis.ticks.y = element_line(color = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        # these plot margins are to leave space for the population name on the big figure
        plot.margin = unit(c(1.0, 0.3, 0, 0.2),"cm"))
# testing theme
# theme(
#       # these plot margins are to leave space for the population name on the big figure
#       plot.margin = unit(c(0.2, 0.2, 0, 0.2),"cm"))

ggsave(here::here("figures", "paper_figures", "test_combined_density_plot.png"), combined_density_plot, height = 8, width = 10)
ggsave(here::here("figures", "paper_figures", "test_combined_spillvolume_plot.png"), combined_rear_spill_move_prob_plot, height = 8, width = 10)

combined_combined_plot <- ggarrange(combined_density_plot, combined_rear_spill_move_prob_plot, nrow = 2, ncol = 1,
                           heights = c(2,6))

ggsave(here::here("figures", "paper_figures", "test_spill_combined_combined_plot.png"), 
       combined_combined_plot, height = 8, width = 10)


#### Generate MCN plot ####

combined_rear_spill_move_prob_plot <- ggplot(subset(combined_spillvolume_move_summary, from == 3), aes(x = spillwindow_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                                                       color = origin_rear, fill = origin_rear)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, color = NA) +
  # geom_rug(data = covariate_experiences, aes(x = spill_actual), inherit.aes = FALSE,
  #          sides = "t", length = unit(0.3, "cm"), outside = TRUE) +
  scale_y_continuous(lim = c(0,0.23), expand = c(0,0), breaks = seq(0, 0.2, 0.05)) +
  xlim(0,150) +
  # common x-axis scale across all populations
  coord_cartesian(xlim = c(0, 150), expand = FALSE) +
  # scale_color_manual(values = movement_colors) +
  # scale_fill_manual(values =  movement_colors) +
  scale_color_tableau(palette = "Tableau 20") +
  scale_fill_tableau(palette = "Tableau 20") +
  xlab("Winter spill days") +
  ylab("Movement probability") +
  theme(# legend.position = "none",
    panel.grid.major = element_line(color = "gray90"),
    panel.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.4),
    # turn off the axis titles on each individual plot and just show one for whole plot
    axis.title = element_blank(),
    # these plot margins are to leave space for the population name on the big figure
    plot.margin = unit(c(0, 0.3, 0.2, 0.2),"cm"))

combined_density_plot <- ggplot(data = subset(combined_covariate_experiences, state == 3 & origin_rear %in% unique(subset(combined_spillvolume_move_summary, from == 3)$origin_rear)), 
                                aes(spill_window_actual, fill = origin_rear))+
  # geom_density(alpha = 0.1) +
  # geom_rug(sides = "t", length = unit(0.2, "cm"), outside = FALSE) +
  # geom_histogram(alpha = 0.5, bins = 60) +
  geom_histogram(aes(y=..count../sum(..count..)), alpha = 0.5, bins = 60, boundary = 0) +
  ylab("Density") +
  # scale_x_continuous(lim = c(floor(min(covariate_experiences$temp_actual)),ceiling(max(covariate_experiences$temp_actual))), expand = c(0,0)) +
  # scale_x_continuous(lim = c(0, 23), expand = c(0,0)) +
  # scale_y_continuous(lim = c(0,0.75), expand = c(0,0),
  #                    breaks = c(0, 0.25, 0.50)) +
  # these labels aren't real, they just need to be the same size as the labels from the temp effect plot
  scale_y_continuous(n.breaks = 2, labels = c("0.00", "1.00")) +
  xlim(0,150) +
  # common x-axis scale across all populations
  coord_cartesian(xlim = c(0, 150), expand = FALSE) +
  scale_fill_tableau(palette = "Tableau 20") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "white"),
        axis.ticks.x = element_blank(),
        # axis.ticks.length.x=unit(.1, "cm"),
        axis.ticks.length.x=unit(0, "cm"),
        axis.ticks.y = element_line(color = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        # these plot margins are to leave space for the population name on the big figure
        plot.margin = unit(c(1.0, 0.3, 0, 0.2),"cm"))
# testing theme
# theme(
#       # these plot margins are to leave space for the population name on the big figure
#       plot.margin = unit(c(0.2, 0.2, 0, 0.2),"cm"))

combined_combined_plot_MCN <- ggarrange(combined_density_plot, combined_rear_spill_move_prob_plot, nrow = 2, ncol = 1,
                                    heights = c(2,6))

ggsave(here::here("figures", "paper_figures", "MCN_test_spill_combined_combined_plot.png"), 
       combined_combined_plot_MCN, height = 8, width = 10)


#### Generate ICH plot ####

combined_rear_spill_move_prob_plot <- ggplot(subset(combined_spillvolume_move_summary, from == 8), aes(x = spillwindow_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                                                       color = origin_rear, fill = origin_rear)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, color = NA) +
  # geom_rug(data = covariate_experiences, aes(x = spill_actual), inherit.aes = FALSE,
  #          sides = "t", length = unit(0.3, "cm"), outside = TRUE) +
  scale_y_continuous(lim = c(0,0.23), expand = c(0,0), breaks = seq(0, 0.2, 0.05)) +
  xlim(0,150) +
  # common x-axis scale across all populations
  coord_cartesian(xlim = c(0, 150), expand = FALSE) +
  # scale_color_manual(values = movement_colors) +
  # scale_fill_manual(values =  movement_colors) +
  scale_color_tableau(palette = "Tableau 20") +
  scale_fill_tableau(palette = "Tableau 20") +
  xlab("Winter spill days") +
  ylab("Movement probability") +
  theme(# legend.position = "none",
    panel.grid.major = element_line(color = "gray90"),
    panel.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.4),
    # turn off the axis titles on each individual plot and just show one for whole plot
    axis.title = element_blank(),
    # these plot margins are to leave space for the population name on the big figure
    plot.margin = unit(c(0, 0.3, 0.2, 0.2),"cm"))

combined_density_plot <- ggplot(data = subset(combined_covariate_experiences, state == 8 & origin_rear %in% unique(subset(combined_spillvolume_move_summary, from == 8)$origin_rear)), 
                                aes(spill_window_actual, fill = origin_rear))+
  # geom_density(alpha = 0.1) +
  # geom_rug(sides = "t", length = unit(0.2, "cm"), outside = FALSE) +
  # geom_histogram(alpha = 0.5, bins = 60) +
  geom_histogram(aes(y=..count../sum(..count..)), alpha = 0.5, bins = 60, boundary = 0) +
  ylab("Density") +
  # scale_x_continuous(lim = c(floor(min(covariate_experiences$temp_actual)),ceiling(max(covariate_experiences$temp_actual))), expand = c(0,0)) +
  # scale_x_continuous(lim = c(0, 23), expand = c(0,0)) +
  # scale_y_continuous(lim = c(0,0.75), expand = c(0,0),
  #                    breaks = c(0, 0.25, 0.50)) +
  # these labels aren't real, they just need to be the same size as the labels from the temp effect plot
  scale_y_continuous(n.breaks = 2, labels = c("0.00", "1.00")) +
  xlim(0,150) +
  # common x-axis scale across all populations
  coord_cartesian(xlim = c(0, 150), expand = FALSE) +
  scale_fill_tableau(palette = "Tableau 20") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "white"),
        axis.ticks.x = element_blank(),
        # axis.ticks.length.x=unit(.1, "cm"),
        axis.ticks.length.x=unit(0, "cm"),
        axis.ticks.y = element_line(color = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        # these plot margins are to leave space for the population name on the big figure
        plot.margin = unit(c(1.0, 0.3, 0, 0.2),"cm"))
# testing theme
# theme(
#       # these plot margins are to leave space for the population name on the big figure
#       plot.margin = unit(c(0.2, 0.2, 0, 0.2),"cm"))

combined_combined_plot_ICH <- ggarrange(combined_density_plot, combined_rear_spill_move_prob_plot, nrow = 2, ncol = 1,
                                        heights = c(2,6))

ggsave(here::here("figures", "paper_figures", "ICH_test_spill_combined_combined_plot.png"), 
       combined_combined_plot_ICH, height = 8, width = 10)


#### Generate LGR plot ####

combined_rear_spill_move_prob_plot <- ggplot(subset(combined_spillvolume_move_summary, from == 9), aes(x = spillwindow_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                                                       color = origin_rear, fill = origin_rear)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, color = NA) +
  # geom_rug(data = covariate_experiences, aes(x = spill_actual), inherit.aes = FALSE,
  #          sides = "t", length = unit(0.3, "cm"), outside = TRUE) +
  scale_y_continuous(lim = c(0,0.23), expand = c(0,0), breaks = seq(0, 0.2, 0.05)) +
  xlim(0,150) +
  # common x-axis scale across all populations
  coord_cartesian(xlim = c(0, 150), expand = FALSE) +
  # scale_color_manual(values = movement_colors) +
  # scale_fill_manual(values =  movement_colors) +
  scale_color_tableau(palette = "Tableau 20") +
  scale_fill_tableau(palette = "Tableau 20") +
  xlab("Winter spill days") +
  ylab("Movement probability") +
  theme(# legend.position = "none",
    panel.grid.major = element_line(color = "gray90"),
    panel.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.4),
    # turn off the axis titles on each individual plot and just show one for whole plot
    axis.title = element_blank(),
    # these plot margins are to leave space for the population name on the big figure
    plot.margin = unit(c(0, 0.3, 0.2, 0.2),"cm"))

combined_density_plot <- ggplot(data = subset(combined_covariate_experiences, state == 9 & origin_rear %in% unique(subset(combined_spillvolume_move_summary, from == 9)$origin_rear)), 
                                aes(spill_window_actual, fill = origin_rear))+
  # geom_density(alpha = 0.1) +
  # geom_rug(sides = "t", length = unit(0.2, "cm"), outside = FALSE) +
  # geom_histogram(alpha = 0.5, bins = 60) +
  geom_histogram(aes(y=..count../sum(..count..)), alpha = 0.5, bins = 60, boundary = 0) +
  ylab("Density") +
  # scale_x_continuous(lim = c(floor(min(covariate_experiences$temp_actual)),ceiling(max(covariate_experiences$temp_actual))), expand = c(0,0)) +
  # scale_x_continuous(lim = c(0, 23), expand = c(0,0)) +
  # scale_y_continuous(lim = c(0,0.75), expand = c(0,0),
  #                    breaks = c(0, 0.25, 0.50)) +
  # these labels aren't real, they just need to be the same size as the labels from the temp effect plot
  scale_y_continuous(n.breaks = 2, labels = c("0.00", "1.00")) +
  xlim(0,150) +
  # common x-axis scale across all populations
  coord_cartesian(xlim = c(0, 150), expand = FALSE) +
  scale_fill_tableau(palette = "Tableau 20") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "white"),
        axis.ticks.x = element_blank(),
        # axis.ticks.length.x=unit(.1, "cm"),
        axis.ticks.length.x=unit(0, "cm"),
        axis.ticks.y = element_line(color = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        # these plot margins are to leave space for the population name on the big figure
        plot.margin = unit(c(1.0, 0.3, 0, 0.2),"cm"))
# testing theme
# theme(
#       # these plot margins are to leave space for the population name on the big figure
#       plot.margin = unit(c(0.2, 0.2, 0, 0.2),"cm"))

combined_combined_plot_LGR <- ggarrange(combined_density_plot, combined_rear_spill_move_prob_plot, nrow = 2, ncol = 1,
                                        heights = c(2,6))

ggsave(here::here("figures", "paper_figures", "LGR_test_spill_combined_combined_plot.png"), 
       combined_combined_plot_LGR, height = 8, width = 10)


