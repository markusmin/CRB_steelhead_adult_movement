# 07-06_Fig5_spilldays_effects_plot

# This script takes the output from the stan model runs and
# plots the effects of temperature

# First, need to load in all of the model runs and all of the packages.
source("R/07_model_analysis/07-01_load_stan_models.R")

#### Extract transition data by spill for rug plot ####


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
    filter(!(state %in% c(0,43))) -> pop_states_dates_covariates
  
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





#### Spill window plot ####

# We will plot all of these - they affect fallback through mainstem states
# Note that spill window parameter, like spill days, is shared for the DPS. So we 
# will only notice differences in overall movement probability in from
# states that have origin effects for other movements (e.g., fallback over BON
# will have same prob for different UC origins, fallback over RIS will not)

# to capture posterior correlations, we'll need to select that same iteration
# for all parameters that affect movement

# Tell it the movements for which you want to estimate temperature effects
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
    possible_movements <- c(possible_movements, 43)
    
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
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
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
    possible_movements <- c(possible_movements, 43)
    
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
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
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
    possible_movements <- c(possible_movements, 43)
    
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
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
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
    possible_movements <- c(possible_movements, 43)
    
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
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
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
    possible_movements <- c(possible_movements, 43)
    
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
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
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
    possible_movements <- c(possible_movements, 43)
    
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
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
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




#### plot spill window function ####

# create df to index to right dam temp for plotting
dam_index <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL", "ICH", "LGR"),
                        state = seq(2,9))

# another plot option to compare between hatchery and wild
rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")

# set up post-overshoot info so that we don't select fish that got the spill days covariate instead
post_overshoot_combos <- list(c(NA), # Asotin Creek
                              c(NA), # Clearwater River
                              c(3,4,5,6,7,8,9), # Deschutes River
                              c(7), # Entiat River
                              c(3,4,5,6,7,8,9), # Fifteenmile Creek
                              c(NA), # Grande Ronde River
                              c(3,4,5,6,7,8,9), # Hood River
                              c(NA), # Imnaha River
                              c(3,4,5,6,7,8,9), # John Day River
                              c(NA), # Methow River
                              c(NA), # Okanogan River
                              c(NA), # Salmon River
                              c(9), # Tucannon River
                              c(3,4,5,6,7,8,9), # Umatilla River
                              c(4,5,6,7,8,9), # Walla Walla River
                              c(6,7), # Wenatchee River
                              c(4,5,6,7,8,9)) # Yakima River

names(post_overshoot_combos) <- natal_origins[order(natal_origins)]


plot_compare_rear_spillwindow_effect <- function(origin_select,
                                           wild_move_prob_array = NULL, hatchery_move_prob_array = NULL,
                                           wild_covariate_experiences = NULL, hatchery_covariate_experiences = NULL,
                                          movements_evaluated, spill_predict,
                                          from, to, plot_title = NULL){
  array_index <- which(movements_evaluated$from == from & movements_evaluated$to == to)
  
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    wild_spill_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(wild_spill_move_prob) <- paste0("iter", 1:niter) 
    wild_spill_move_prob$spill <- spill_predict
    
    wild_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_spill_move_prob_quantiles
    
    wild_spill_move_prob_quantiles -> rear_spill_move_prob_quantiles
    
    # convert back to kcfs
    rear_spill_move_prob_quantiles %>% 
      mutate(spill_actual = spill*100) -> rear_spill_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences

    wild_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    hatchery_spill_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(hatchery_spill_move_prob) <- paste0("iter", 1:niter) 
    hatchery_spill_move_prob$spill <- spill_predict
    
    hatchery_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_spill_move_prob_quantiles
    
    hatchery_spill_move_prob_quantiles -> rear_spill_move_prob_quantiles
    
    # convert back to kcfs
    rear_spill_move_prob_quantiles %>% 
      mutate(spill_actual = spill*100) -> rear_spill_move_prob_quantiles
    
    # get data organized for rug plot
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    
    hatchery_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
  } 
  # else run both
  else {
    wild_spill_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(wild_spill_move_prob) <- paste0("iter", 1:niter) 
    wild_spill_move_prob$spill <- spill_predict
    
    wild_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_spill_move_prob_quantiles
    
    hatchery_spill_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(hatchery_spill_move_prob) <- paste0("iter", 1:niter) 
    hatchery_spill_move_prob$spill <- spill_predict
    
    hatchery_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_spill_move_prob_quantiles
    
    wild_spill_move_prob_quantiles %>% 
      bind_rows(., hatchery_spill_move_prob_quantiles) -> rear_spill_move_prob_quantiles
    
    # convert back to kcfs
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
  
  # Remove any instances for rug plot of fish that would have received winter spill days covariate
  # if there aren't any for this combination, don't try to manipulate it
  if(nrow(covariate_experiences) != 0) {
    covariate_experiences %>% 
      dplyr::rename(date_numeric = date) %>% 
      # keep only jan/feb/mar 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) %>% 
      filter(!(state %in% post_overshoot_combos[[origin_select]] & month %in% c(1,2,3))) -> covariate_experiences
    
  } else {
    
    
  }

  
  rear_spill_move_prob_plot <- ggplot(rear_spill_move_prob_quantiles, aes(x = spill_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                          color = rear, fill = rear)) +
    geom_line() +
    geom_ribbon(alpha = 0.2, color = NA) +
    geom_rug(data = subset(covariate_experiences, rear == "wild"), aes(x = spill_window*100, color = rear), inherit.aes = FALSE,
             sides = "t", length = unit(0.3, "cm")) +
    geom_rug(data = subset(covariate_experiences, rear == "hatchery"), aes(x = spill_window*100, color = rear), inherit.aes = FALSE,
             sides = "b", length = unit(0.3, "cm")) +
    scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
    scale_x_continuous(lim = c(-5,NA), expand = c(0,0)) +
    scale_color_manual(values = rear_colors) +
    scale_fill_manual(values = rear_colors) +
    xlab("Spill (kcfs)") +
    ylab("Movement probability") +
    ggtitle(plot_title)
    
  return(rear_spill_move_prob_plot)
}


# #### Plot rear type comparison fallback spill window for all ####
# mainstem_fallback_movements <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL", "ICH", "LGR"),
#   from = c(2,3,4,5,6,7,8,9), to = c(1,2,3,4,5,6,3,8))
# 
# #### Upper Columbia ####
# UCH_spill_window_data <- UCH_envir$data$spill_window_data
# UCW_spill_window_data <- UCW_envir$data$spill_window_data
# UC_spill_window_data <- bind_rows(as.data.frame(UCH_spill_window_data), as.data.frame(UCW_spill_window_data))
# 
# # Wenatchee - evaluate all fallback move probs
# WEN_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCW(origin_select= "Wenatchee River", movements = mainstem_fallback_movements)
# WEN_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCH(origin_select= "Wenatchee River", movements = mainstem_fallback_movements)
# 
# # Wenatchee - get covariate experiences
# WEN_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Wenatchee River")
# WEN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Wenatchee River")
# 
# # Wenatchee - use a loop to plot all of them
# WEN_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   WEN_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect(origin_select = "Wenatchee River",
#                                                                         wild_move_prob_array = WEN_wild_spillwindow_move_prob_array,
#                                                                         hatchery_move_prob_array = WEN_hatchery_spillwindow_move_prob_array,
#                                                                         wild_covariate_experiences = WEN_wild_covariate_experiences, 
#                                                                         hatchery_covariate_experiences = WEN_hatchery_covariate_experiences,
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(UC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Wenatchee - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spillwindow_v3", paste0("WEN_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          WEN_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # Entiat - evaluate all fallback move probs
# # Only wild for entiat
# ENT_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCW(origin_select= "Entiat River", movements = mainstem_fallback_movements)
# 
# # Entiat - get covariate experiences
# ENT_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Entiat River")
# 
# # Entiat - use a loop to plot all of them
# ENT_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   ENT_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect(origin_select = "Entiat River",
#                                                                         wild_move_prob_array = ENT_wild_spillwindow_move_prob_array,
#                                                                         wild_covariate_experiences = ENT_wild_covariate_experiences, 
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(UC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Entiat - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spillwindow_v3", paste0("ENT_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          ENT_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # Methow - evaluate all fallback move probs
# MET_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCW(origin_select= "Methow River", movements = mainstem_fallback_movements)
# MET_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCH(origin_select= "Methow River", movements = mainstem_fallback_movements)
# 
# # Methow - get covariate experiences
# MET_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Methow River")
# MET_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Methow River")
# 
# # Methow - use a loop to plot all of them
# MET_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   MET_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect(origin_select = "Methow River",
#                                                                         wild_move_prob_array = MET_wild_spillwindow_move_prob_array,
#                                                                         hatchery_move_prob_array = MET_hatchery_spillwindow_move_prob_array,
#                                                                         wild_covariate_experiences = MET_wild_covariate_experiences, 
#                                                                         hatchery_covariate_experiences = MET_hatchery_covariate_experiences,
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(UC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Methow - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spillwindow_v3", paste0("MET_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          MET_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # Okanoagn - evaluate all fallback move probs
# # Only hatchery for Okanogan
# OKA_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCH(origin_select= "Okanogan River", movements = mainstem_fallback_movements)
# 
# # Okanogan - get covariate experiences
# OKA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Okanogan River")
# 
# # Okanogan - use a loop to plot all of them
# OKA_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   OKA_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect(origin_select = "Okanogan River",
#                                                                         hatchery_move_prob_array = OKA_hatchery_spillwindow_move_prob_array,
#                                                                         hatchery_covariate_experiences = OKA_hatchery_covariate_experiences,
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(UC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Okanogan - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spillwindow_v3", paste0("OKA_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          OKA_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# #### Middle Columbia ####
# MCH_spill_window_data <- MCH_envir$data$spill_window_data
# MCW_spill_window_data <- MCW_envir$data$spill_window_data
# MC_spill_window_data <- bind_rows(as.data.frame(MCH_spill_window_data), as.data.frame(MCW_spill_window_data))
# 
# # Deschutes - evaluate all fallback move probs
# # Deschutes is wild only
# DES_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "Deschutes River", movements = mainstem_fallback_movements)
# 
# # Deschutes - get covariate experiences
# DES_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Deschutes River")
# 
# # Deschutes - use a loop to plot all of them
# DES_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   DES_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect(origin_select = "Deschutes River",
#                                                                         wild_move_prob_array = DES_wild_spillwindow_move_prob_array,
#                                                                         
#                                                                         wild_covariate_experiences = DES_wild_covariate_experiences, 
#                                                                         
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Deschutes - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spillwindow_v3", paste0("DES_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          DES_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # John Day - evaluate all fallback move probs
# # JDR is wild only
# JDR_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "John Day River", movements = mainstem_fallback_movements)
# 
# # John Day - get covariate experiences
# JDR_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "John Day River")
# 
# # John Day - use a loop to plot all of them
# JDR_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   JDR_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect(origin_select = "John Day River",
#                                                                         wild_move_prob_array = JDR_wild_spillwindow_move_prob_array,
#                                                                         
#                                                                         wild_covariate_experiences = JDR_wild_covariate_experiences, 
#                                                                         
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("John Day - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spillwindow_v3", paste0("JDR_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          JDR_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# 
# # Fifteenmile - evaluate all fallback move probs
# # Fifteenmile is wild only
# FIF_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "Fifteenmile Creek", movements = mainstem_fallback_movements)
# 
# # Fifteenmile - get covariate experiences
# FIF_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Fifteenmile Creek")
# 
# # Fifteenmile - use a loop to plot all of them
# FIF_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   FIF_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect(origin_select = "Fifteenmile Creek",
#                                                                         wild_move_prob_array = FIF_wild_spillwindow_move_prob_array,
#                                                                         
#                                                                         wild_covariate_experiences = FIF_wild_covariate_experiences, 
#                                                                         
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Fifteenmile - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spillwindow_v3", paste0("FIF_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          FIF_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # Umatilla - evaluate all fallback move probs
# UMA_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "Umatilla River", movements = mainstem_fallback_movements)
# UMA_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCH(origin_select= "Umatilla River", movements = mainstem_fallback_movements)
# 
# # Umatilla - get covariate experiences
# UMA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Umatilla River")
# UMA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Umatilla River")
# 
# # Umatilla - use a loop to plot all of them
# UMA_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   UMA_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect(origin_select = "Umatilla River",
#                                                                         wild_move_prob_array = UMA_wild_spillwindow_move_prob_array,
#                                                                         hatchery_move_prob_array = UMA_hatchery_spillwindow_move_prob_array,
#                                                                         wild_covariate_experiences = UMA_wild_covariate_experiences, 
#                                                                         hatchery_covariate_experiences = UMA_hatchery_covariate_experiences,
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Umatilla - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spillwindow_v3", paste0("UMA_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          UMA_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # Yakima - evaluate all fallback move probs
# # Yakima is wild only
# YAK_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "Yakima River", movements = mainstem_fallback_movements)
# 
# # Yakima - get covariate experiences
# YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Yakima River")
# 
# # Yakima - use a loop to plot all of them
# YAK_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   YAK_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect(origin_select = "Yakima River",
#                                                                         wild_move_prob_array = YAK_wild_spillwindow_move_prob_array,
#                                                                         
#                                                                         wild_covariate_experiences = YAK_wild_covariate_experiences, 
#                                                                         
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Yakima - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spillwindow_v3", paste0("YAK_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          YAK_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # Walla Walla - evaluate all fallback move probs
# WAWA_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "Walla Walla River", movements = mainstem_fallback_movements)
# WAWA_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCH(origin_select= "Walla Walla River", movements = mainstem_fallback_movements)
# 
# # Walla Walla - get covariate experiences
# WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Walla Walla River")
# WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Walla Walla River")
# 
# # Walla Walla - use a loop to plot all of them
# WAWA_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   WAWA_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect(origin_select = "Walla Walla River",
#                                                                          wild_move_prob_array = WAWA_wild_spillwindow_move_prob_array,
#                                                                          hatchery_move_prob_array = WAWA_hatchery_spillwindow_move_prob_array,
#                                                                          wild_covariate_experiences = WAWA_wild_covariate_experiences, 
#                                                                          hatchery_covariate_experiences = WAWA_hatchery_covariate_experiences,
#                                                                          movements_evaluated = mainstem_fallback_movements,
#                                                                          spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                          from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                          plot_title = paste0("Walla Walla - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spillwindow_v3", paste0("WAWA_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          WAWA_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# #### Snake River ####
# SRH_spill_window_data <- SRH_envir$data$spill_window_data
# SRW_spill_window_data <- SRW_envir$data$spill_window_data
# SR_spill_window_data <- bind_rows(as.data.frame(SRH_spill_window_data), as.data.frame(SRW_spill_window_data))
# 
# # Asotin Creek - evaluate all fallback move probs
# # Asotin is wild only
# ASO_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Asotin Creek", movements = mainstem_fallback_movements)
# 
# # Asotin Creek - get covariate experiences
# ASO_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Asotin Creek")
# 
# # Asotin Creek - use a loop to plot all of them
# ASO_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   ASO_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect(origin_select = "Asotin Creek",
#                                                                         wild_move_prob_array = ASO_wild_spillwindow_move_prob_array,
#                                                                         
#                                                                         wild_covariate_experiences = ASO_wild_covariate_experiences, 
#                                                                         
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Asotin Creek - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spillwindow_v3", paste0("ASO_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          ASO_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # Clearwater River - evaluate all fallback move probs
# CLE_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Clearwater River", movements = mainstem_fallback_movements)
# CLE_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select= "Clearwater River", movements = mainstem_fallback_movements)
# 
# # Clearwater River - get covariate experiences
# CLE_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Clearwater River")
# CLE_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Clearwater River")
# 
# # Clearwater River - use a loop to plot all of them
# CLE_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   CLE_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect(origin_select = "Clearwater River",
#                                                                         wild_move_prob_array = CLE_wild_spillwindow_move_prob_array,
#                                                                         hatchery_move_prob_array = CLE_hatchery_spillwindow_move_prob_array,
#                                                                         wild_covariate_experiences = CLE_wild_covariate_experiences, 
#                                                                         hatchery_covariate_experiences = CLE_hatchery_covariate_experiences,
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Clearwater River - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spillwindow_v3", paste0("CLE_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          CLE_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # Grande Ronde River - evaluate all fallback move probs
# GR_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Grande Ronde River", movements = mainstem_fallback_movements)
# GR_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select= "Grande Ronde River", movements = mainstem_fallback_movements)
# 
# # Grande Ronde River - get covariate experiences
# GR_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Grande Ronde River")
# GR_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Grande Ronde River")
# 
# # Grande Ronde River - use a loop to plot all of them
# GR_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   GR_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect(origin_select = "Grande Ronde River",
#                                                                        wild_move_prob_array = GR_wild_spillwindow_move_prob_array,
#                                                                        hatchery_move_prob_array = GR_hatchery_spillwindow_move_prob_array,
#                                                                        wild_covariate_experiences = GR_wild_covariate_experiences, 
#                                                                        hatchery_covariate_experiences = GR_hatchery_covariate_experiences,
#                                                                        movements_evaluated = mainstem_fallback_movements,
#                                                                        spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                        from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                        plot_title = paste0("Grande Ronde River - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spillwindow_v3", paste0("GR_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          GR_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # Imnaha River - evaluate all fallback move probs
# IMN_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Imnaha River", movements = mainstem_fallback_movements)
# IMN_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select= "Imnaha River", movements = mainstem_fallback_movements)
# 
# # Imnaha River - get covariate experiences
# IMN_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Imnaha River")
# IMN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Imnaha River")
# 
# # Imnaha River - use a loop to plot all of them
# IMN_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   IMN_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect(origin_select = "Imnaha River",
#                                                                         wild_move_prob_array = IMN_wild_spillwindow_move_prob_array,
#                                                                         hatchery_move_prob_array = IMN_hatchery_spillwindow_move_prob_array,
#                                                                         wild_covariate_experiences = IMN_wild_covariate_experiences, 
#                                                                         hatchery_covariate_experiences = IMN_hatchery_covariate_experiences,
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Imnaha River - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spillwindow_v3", paste0("IMN_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          IMN_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # Salmon River - evaluate all fallback move probs
# SAL_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Salmon River", movements = mainstem_fallback_movements)
# SAL_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select= "Salmon River", movements = mainstem_fallback_movements)
# 
# # Salmon River - get covariate experiences
# SAL_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Salmon River")
# SAL_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Salmon River")
# 
# # Salmon River - use a loop to plot all of them
# SAL_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   SAL_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect(origin_select = "Salmon River",
#                                                                         wild_move_prob_array = SAL_wild_spillwindow_move_prob_array,
#                                                                         hatchery_move_prob_array = SAL_hatchery_spillwindow_move_prob_array,
#                                                                         wild_covariate_experiences = SAL_wild_covariate_experiences, 
#                                                                         hatchery_covariate_experiences = SAL_hatchery_covariate_experiences,
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Salmon River - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spillwindow_v3", paste0("SAL_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          SAL_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# # Tucannon River - evaluate all fallback move probs
# TUC_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Tucannon River", movements = mainstem_fallback_movements)
# TUC_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select= "Tucannon River", movements = mainstem_fallback_movements)
# 
# # Tucannon River - get covariate experiences
# TUC_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Tucannon River")
# TUC_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Tucannon River")
# 
# # Tucannon River - use a loop to plot all of them
# TUC_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   TUC_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect(origin_select = "Tucannon River",
#                                                                         wild_move_prob_array = TUC_wild_spillwindow_move_prob_array,
#                                                                         hatchery_move_prob_array = TUC_hatchery_spillwindow_move_prob_array,
#                                                                         wild_covariate_experiences = TUC_wild_covariate_experiences, 
#                                                                         hatchery_covariate_experiences = TUC_hatchery_covariate_experiences,
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Tucannon River - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spillwindow_v3", paste0("TUC_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          TUC_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# 

#### Spill days movement probability function ####
estimate_spilldays_effect_UCW <- function(origin_select, movements){
  
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
    
    # set up spill window values to predict across
    # extract the data
    winterspill_data <- UCW_envir$data$winter_spill_days_data
    # create a sequence of values from 1 to max
    spill_days_predict <- seq(0, max(winterspill_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    UCW_states_dates <- data.frame(state = as.vector(UCW_envir$data$y),
                                   date_numeric = as.vector(UCW_envir$data$transition_dates))
    
    # add month to states dates so that we can subset
    UCW_states_dates %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> UCW_states_dates
    
    # get median temperature for this state - only in Jan/Feb/Mar
    temp_data <- as.data.frame(UCW_envir$data$temperature)
    temp_data %>% 
      rownames_to_column("date_numeric") %>% 
      pivot_longer(cols = -c(date_numeric)) -> temp_data_long
    
    temp_data_long %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> temp_data_long
    
    # only keep actual observed fish
    med_temp <- median(temp_data[subset(UCW_states_dates, state == from & month %in% c(1,2,3))$date_numeric,from])
    
    # do not include spill window, because code is set up to index so that it's only one or the other
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_days_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_UCW[movements$from[i],movements$to[i],iter] +
                                                 btemp0_array_UCW[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # btemp1_array_UCW[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # bspillwindow_array_UCW[movements$from[i],movements$to[i],iter]*med_spillwindow + # no spill window if spill days is indexed
                                                 bwinterspill_array_UCW[movements$from[i],movements$to[i],iter]*spill_days_predict[j] +
                                                 btemp0xorigin1_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp1xorigin1_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 btemp0xorigin2_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin2 +
                                                 # btemp1xorigin2_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 btemp0xorigin3_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 # btemp1xorigin3_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 borigin1_array_UCW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_UCW[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_UCW[movements$from[i],movements$to[i],iter]*origin3)/
          sum(exp(b0_array_UCW[movements$from[i],possible_movements,iter] +
                    btemp0_array_UCW[movements$from[i],possible_movements,iter]*med_temp + # select temp0 for any winterspill estimates
                    # btemp1_array_UCW[movements$from[i],possible_movements,iter]*med_temp + 
                    # bspillwindow_array_UCW[movements$from[i],possible_movements,iter]*med_spillwindow + # no spill window if spill days is indexed
                    bwinterspill_array_UCW[movements$from[i],possible_movements,iter]*spill_days_predict[j] +
                    btemp0xorigin1_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp1xorigin1_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    btemp0xorigin2_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin2 +
                    # btemp1xorigin2_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    btemp0xorigin3_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    # btemp1xorigin3_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    borigin1_array_UCW[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_UCW[movements$from[i],possible_movements,iter]*origin2 +
                    borigin3_array_UCW[movements$from[i],possible_movements,iter]*origin3))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}

estimate_spilldays_effect_UCH <- function(origin_select, movements){
  
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
    
    # set up spill window values to predict across
    # extract the data
    winterspill_data <- UCH_envir$data$winter_spill_days_data
    # create a sequence of values from 1 to max
    spill_days_predict <- seq(0, max(winterspill_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    UCH_states_dates <- data.frame(state = as.vector(UCH_envir$data$y),
                                   date_numeric = as.vector(UCH_envir$data$transition_dates))
    
    # add month to states dates so that we can subset
    UCH_states_dates %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> UCH_states_dates
    
    # get median temperature for this state - only in Jan/Feb/Mar
    temp_data <- as.data.frame(UCH_envir$data$temperature)
    temp_data %>% 
      rownames_to_column("date_numeric") %>% 
      pivot_longer(cols = -c(date_numeric)) -> temp_data_long
    
    temp_data_long %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> temp_data_long
    
    # only keep actual observed fish
    med_temp <- median(temp_data[subset(UCH_states_dates, state == from & month %in% c(1,2,3))$date_numeric,from])
    
    # do not include spill window, because code is set up to index so that it's only one or the other
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_days_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_UCH[movements$from[i],movements$to[i],iter] +
                                                 btemp0_array_UCH[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # btemp1_array_UCH[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # bspillwindow_array_UCH[movements$from[i],movements$to[i],iter]*med_spillwindow + # no spill window if spill days is indexed
                                                 bwinterspill_array_UCH[movements$from[i],movements$to[i],iter]*spill_days_predict[j] +
                                                 btemp0xorigin1_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp1xorigin1_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 btemp0xorigin2_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin2 +
                                                 # btemp1xorigin2_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 btemp0xorigin3_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 # btemp1xorigin3_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 borigin1_array_UCH[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_UCH[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_UCH[movements$from[i],movements$to[i],iter]*origin3)/
          sum(exp(b0_array_UCH[movements$from[i],possible_movements,iter] +
                    btemp0_array_UCH[movements$from[i],possible_movements,iter]*med_temp + # select temp0 for any winterspill estimates
                    # btemp1_array_UCH[movements$from[i],possible_movements,iter]*med_temp + 
                    # bspillwindow_array_UCH[movements$from[i],possible_movements,iter]*med_spillwindow + # no spill window if spill days is indexed
                    bwinterspill_array_UCH[movements$from[i],possible_movements,iter]*spill_days_predict[j] +
                    btemp0xorigin1_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp1xorigin1_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    btemp0xorigin2_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin2 +
                    # btemp1xorigin2_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    btemp0xorigin3_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    # btemp1xorigin3_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    borigin1_array_UCH[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_UCH[movements$from[i],possible_movements,iter]*origin2 +
                    borigin3_array_UCH[movements$from[i],possible_movements,iter]*origin3))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}


estimate_spilldays_effect_MCW <- function(origin_select, movements){
  
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
    winterspill_data <- MCW_envir$data$winter_spill_days_data
    # create a sequence of values from 1 to max
    spill_days_predict <- seq(0, max(winterspill_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    MCW_states_dates <- data.frame(state = as.vector(MCW_envir$data$y),
                                   date_numeric = as.vector(MCW_envir$data$transition_dates))
    
    # add month to states dates so that we can subset
    MCW_states_dates %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> MCW_states_dates
    
    # get median temperature for this state - only in Jan/Feb/Mar
    temp_data <- as.data.frame(MCW_envir$data$temperature)
    temp_data %>% 
      rownames_to_column("date_numeric") %>% 
      pivot_longer(cols = -c(date_numeric)) -> temp_data_long
    
    temp_data_long %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> temp_data_long
    
    # only keep actual observed fish
    med_temp <- median(temp_data[subset(MCW_states_dates, state == from & month %in% c(1,2,3))$date_numeric,from])
    
    # do not include spill window, because code is set up to index so that it's only one or the other
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_days_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_MCW[movements$from[i],movements$to[i],iter] +
                                                 btemp0_array_MCW[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # btemp1_array_MCW[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # bspillwindow_array_MCW[movements$from[i],movements$to[i],iter]*med_spillwindow + # no spill window if spill days is indexed
                                                 bwinterspill_array_MCW[movements$from[i],movements$to[i],iter]*spill_days_predict[j] +
                                                 btemp0xorigin1_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp1xorigin1_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 btemp0xorigin2_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin2 +
                                                 # btemp1xorigin2_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 btemp0xorigin3_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 # btemp1xorigin3_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 btemp0xorigin4_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 # btemp1xorigin4_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 btemp0xorigin5_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 # btemp1xorigin5_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 btemp0xorigin6_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin6 +
                                                 # btemp1xorigin6_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin6 +
                                                 borigin1_array_MCW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_MCW[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_MCW[movements$from[i],movements$to[i],iter]*origin3 +
                                                 borigin4_array_MCW[movements$from[i],movements$to[i],iter]*origin4 +
                                                 borigin5_array_MCW[movements$from[i],movements$to[i],iter]*origin5 +
                                                 borigin6_array_MCW[movements$from[i],movements$to[i],iter]*origin6)/
          sum(exp(b0_array_MCW[movements$from[i],possible_movements,iter] +
                    btemp0_array_MCW[movements$from[i],possible_movements,iter]*med_temp + # select temp0 for any winterspill estimates
                    # btemp1_array_MCW[movements$from[i],possible_movements,iter]*med_temp + 
                    # bspillwindow_array_MCW[movements$from[i],possible_movements,iter]*med_spillwindow + # no spill window if spill days is indexed
                    bwinterspill_array_MCW[movements$from[i],possible_movements,iter]*spill_days_predict[j] +
                    btemp0xorigin1_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp1xorigin1_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    btemp0xorigin2_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin2 +
                    # btemp1xorigin2_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    btemp0xorigin3_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    # btemp1xorigin3_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    btemp0xorigin4_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    # btemp1xorigin4_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    btemp0xorigin5_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin5 +
                    # btemp1xorigin5_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin5 +
                    btemp0xorigin6_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin6 +
                    # btemp1xorigin6_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin6 +
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

estimate_spilldays_effect_MCH <- function(origin_select, movements){
  
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
    
    # set up spill window values to predict across
    # extract the data
    winterspill_data <- MCH_envir$data$winter_spill_days_data
    # create a sequence of values from 1 to max
    spill_days_predict <- seq(0, max(winterspill_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    MCH_states_dates <- data.frame(state = as.vector(MCH_envir$data$y),
                                   date_numeric = as.vector(MCH_envir$data$transition_dates))
    
    # add month to states dates so that we can subset
    MCH_states_dates %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> MCH_states_dates
    
    # get median temperature for this state - only in Jan/Feb/Mar
    temp_data <- as.data.frame(MCH_envir$data$temperature)
    temp_data %>% 
      rownames_to_column("date_numeric") %>% 
      pivot_longer(cols = -c(date_numeric)) -> temp_data_long
    
    temp_data_long %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> temp_data_long
    
    # only keep actual observed fish
    med_temp <- median(temp_data[subset(MCH_states_dates, state == from & month %in% c(1,2,3))$date_numeric,from])
    
    # do not include spill window, because code is set up to index so that it's only one or the other
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_days_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_MCH[movements$from[i],movements$to[i],iter] +
                                                 btemp0_array_MCH[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # btemp1_array_MCH[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # bspillwindow_array_MCH[movements$from[i],movements$to[i],iter]*med_spillwindow + # no spill window if spill days is indexed
                                                 bwinterspill_array_MCH[movements$from[i],movements$to[i],iter]*spill_days_predict[j] +
                                                 btemp0xorigin1_array_MCH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp1xorigin1_array_MCH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 btemp0xorigin2_array_MCH[movements$from[i],movements$to[i],iter]*med_temp*origin2 +
                                                 # btemp1xorigin2_array_MCH[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 borigin1_array_MCH[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_MCH[movements$from[i],movements$to[i],iter]*origin2)/
          sum(exp(b0_array_MCH[movements$from[i],possible_movements,iter] +
                    btemp0_array_MCH[movements$from[i],possible_movements,iter]*med_temp + # select temp0 for any winterspill estimates
                    # btemp1_array_MCH[movements$from[i],possible_movements,iter]*med_temp + 
                    # bspillwindow_array_MCH[movements$from[i],possible_movements,iter]*med_spillwindow + # no spill window if spill days is indexed
                    bwinterspill_array_MCH[movements$from[i],possible_movements,iter]*spill_days_predict[j] +
                    btemp0xorigin1_array_MCH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp1xorigin1_array_MCH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    btemp0xorigin2_array_MCH[movements$from[i],possible_movements,iter]*med_temp*origin2 +
                    # btemp1xorigin2_array_MCH[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    borigin1_array_MCH[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_MCH[movements$from[i],possible_movements,iter]*origin2))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}

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
    possible_movements <- c(possible_movements, 43)
    
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
    possible_movements <- c(possible_movements, 43)
    
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



#### plot spill days function ####

# create df to index to right dam temp for plotting
dam_index <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL", "ICH", "LGR"),
                        state = seq(2,9))

# another plot option to compare between hatchery and wild
rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")

plot_compare_rear_spilldays_effect <- function(origin_select,
                                           wild_move_prob_array = NULL, hatchery_move_prob_array = NULL,
                                           wild_covariate_experiences = NULL, hatchery_covariate_experiences = NULL,
                                           movements_evaluated, spill_predict,
                                           from, to, plot_title = NULL){
  array_index <- which(movements_evaluated$from == from & movements_evaluated$to == to)
  
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    wild_spill_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(wild_spill_move_prob) <- paste0("iter", 1:niter) 
    wild_spill_move_prob$spill <- spill_predict
    
    wild_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_spill_move_prob_quantiles
    
    wild_spill_move_prob_quantiles -> rear_spill_move_prob_quantiles
    
    # convert back to kcfs
    rear_spill_move_prob_quantiles %>% 
      mutate(spill_actual = spill*100) -> rear_spill_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    
    wild_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    hatchery_spill_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(hatchery_spill_move_prob) <- paste0("iter", 1:niter) 
    hatchery_spill_move_prob$spill <- spill_predict
    
    hatchery_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_spill_move_prob_quantiles
    
    hatchery_spill_move_prob_quantiles -> rear_spill_move_prob_quantiles
    
    # convert back to kcfs
    rear_spill_move_prob_quantiles %>% 
      mutate(spill_actual = spill*100) -> rear_spill_move_prob_quantiles
    
    # get data organized for rug plot
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    
    hatchery_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from)  -> covariate_experiences
    
  } 
  # else run both
  else {
    wild_spill_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(wild_spill_move_prob) <- paste0("iter", 1:niter) 
    wild_spill_move_prob$spill <- spill_predict
    
    wild_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_spill_move_prob_quantiles
    
    hatchery_spill_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(hatchery_spill_move_prob) <- paste0("iter", 1:niter) 
    hatchery_spill_move_prob$spill <- spill_predict
    
    hatchery_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_spill_move_prob_quantiles
    
    wild_spill_move_prob_quantiles %>% 
      bind_rows(., hatchery_spill_move_prob_quantiles) -> rear_spill_move_prob_quantiles
    
    # convert back to days
    rear_spill_move_prob_quantiles %>% 
      mutate(spill_actual = spill*100) -> rear_spill_move_prob_quantiles
    
    # get data organized for rug plot
    # filter to only keep jan/feb/mar
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    wild_covariate_experiences %>% 
      bind_rows(., hatchery_covariate_experiences) %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
    
    
  }
  
  # Keep only fish where they could have experienced spill conditions for rug plot
  covariate_experiences %>% 
    filter(winter_post_overshoot_vector == 1) -> covariate_experiences
  
  rear_spill_move_prob_plot <- ggplot(rear_spill_move_prob_quantiles, aes(x = spill_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                          color = rear, fill = rear)) +
    geom_line() +
    geom_ribbon(alpha = 0.2, color = NA) +
    geom_rug(data = subset(covariate_experiences, rear == "wild"), aes(x = winter_spill*100, color = rear, y = 0), inherit.aes = FALSE,
             sides = "t", length = unit(0.3, "cm"), position = "jitter") +
    geom_rug(data = subset(covariate_experiences, rear == "hatchery"), aes(x = winter_spill*100, color = rear, y = 0), inherit.aes = FALSE,
             sides = "b", length = unit(0.3, "cm"), position = "jitter") +
    scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
    scale_x_continuous(lim = c(-1,max(rear_spill_move_prob_quantiles$spill_actual)+1), expand = c(0,0)) +
    scale_color_manual(values = rear_colors) +
    scale_fill_manual(values = rear_colors) +
    xlab("Days of winter spill") +
    ylab("Movement probability") +
    ggtitle(plot_title)
  
  return(rear_spill_move_prob_plot)
}



# #### Plot rear type comparison fallback spill days for post-overshot fallback movements ####
# UC_postovershoot_fallback_movements <- data.frame(dam = c("RRE", "WEL"),
#                                           from = c(6,7), to = c(5,6))
# 
# MC_postovershoot_fallback_movements <- data.frame(dam = c("MCN", "PRA", "ICH"),
#                                                   from = c(3,4,8), to = c(2,3,3))
# 
# SR_postovershoot_fallback_movements <- data.frame(dam = c("LGR"),
#                                                   from = c(9), to = c(8))
# 
# #### Upper Columbia ####
# # only Entiat and Wenatchee
# 
# UCH_spill_days_data <- UCH_envir$data$winter_spill_days_data
# UCW_spill_days_data <- UCW_envir$data$winter_spill_days_data
# UC_spill_days_data <- bind_rows(as.data.frame(UCH_spill_days_data), as.data.frame(UCW_spill_days_data))
# 
# # Wenatchee - evaluate all fallback move probs
# WEN_wild_spilldays_move_prob_array <- estimate_spilldays_effect_UCW(origin_select= "Wenatchee River", movements = UC_postovershoot_fallback_movements)
# WEN_hatchery_spilldays_move_prob_array <- estimate_spilldays_effect_UCH(origin_select= "Wenatchee River", movements = UC_postovershoot_fallback_movements)
# 
# # Wenatchee - get covariate experiences
# WEN_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Wenatchee River")
# WEN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Wenatchee River")
# 
# # Wenatchee - use a loop to plot all of them
# WEN_fallback_spilldays_plots <- vector(mode = "list", length = nrow(UC_postovershoot_fallback_movements))
# for (i in 1:nrow(UC_postovershoot_fallback_movements)){
#   # if a movement doesn't have any data, don't plot it
#   filter(WEN_wild_covariate_experiences, state == UC_postovershoot_fallback_movements$from[i]) -> WEN_data_to_plot
#   
#   if(nrow(WEN_data_to_plot)>0){
#     WEN_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "Wenatchee River",
#                                                                           wild_move_prob_array = WEN_wild_spilldays_move_prob_array,
#                                                                           hatchery_move_prob_array = WEN_hatchery_spilldays_move_prob_array,
#                                                                           wild_covariate_experiences = WEN_wild_covariate_experiences, 
#                                                                           hatchery_covariate_experiences = WEN_hatchery_covariate_experiences,
#                                                                           movements_evaluated = UC_postovershoot_fallback_movements,
#                                                                           spill_predict = seq(0, max(UC_spill_days_data[,UC_postovershoot_fallback_movements$from[i]]),length = 100),
#                                                                           from = UC_postovershoot_fallback_movements$from[i], to = UC_postovershoot_fallback_movements$to[i], 
#                                                                           plot_title = paste0("Wenatchee - post-overshoot fallback at ", UC_postovershoot_fallback_movements$dam[i]))
#     
#     ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays_v3", paste0("WEN_compare_fallback_", 
#                                                                                         UC_postovershoot_fallback_movements$dam[i],
#                                                                                           "_spilldays.png")),
#            WEN_fallback_spilldays_plots[[i]], height = 8, width = 8)
#   } else {
#     # if it is zero, do nothing
#     
#   }
#   
# }
# 
# # Entiat - evaluate all fallback move probs
# # Only wild for entiat
# ENT_wild_spilldays_move_prob_array <- estimate_spilldays_effect_UCW(origin_select= "Entiat River", movements = UC_postovershoot_fallback_movements)
# 
# # Entiat - get covariate experiences
# ENT_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Entiat River")
# 
# # Entiat - use a loop to plot all of them
# ENT_fallback_spilldays_plots <- vector(mode = "list", length = nrow(UC_postovershoot_fallback_movements))
# for (i in 1:nrow(UC_postovershoot_fallback_movements)){
#   # if a movement doesn't have any data, don't plot it
#   filter(ENT_wild_covariate_experiences, state == UC_postovershoot_fallback_movements$from[i]) -> ENT_data_to_plot
#   
#   if(nrow(ENT_data_to_plot)>0){
#     
#     ENT_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "Entiat River",
#                                                                             wild_move_prob_array = ENT_wild_spilldays_move_prob_array,
#                                                                             wild_covariate_experiences = ENT_wild_covariate_experiences, 
#                                                                             movements_evaluated = UC_postovershoot_fallback_movements,
#                                                                             spill_predict = seq(0, max(UC_spill_days_data[,UC_postovershoot_fallback_movements$from[i]]),length = 100),
#                                                                             from = UC_postovershoot_fallback_movements$from[i], to = UC_postovershoot_fallback_movements$to[i], 
#                                                                             plot_title = paste0("Entiat - post-overshoot fallback at ", UC_postovershoot_fallback_movements$dam[i]))
#     
#     ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays_v3", paste0("ENT_compare_fallback_", 
#                                                                                         UC_postovershoot_fallback_movements$dam[i],
#                                                                                         "_spilldays.png")),
#            ENT_fallback_spilldays_plots[[i]], height = 8, width = 8)
#   } else {
#     # if it is zero, do nothing
#     
#   }
# }
# 
# 
# 
# #### Middle Columbia ####
# MCH_spill_days_data <- MCH_envir$data$winter_spill_days_data
# MCW_spill_days_data <- MCW_envir$data$winter_spill_days_data
# MC_spill_days_data <- bind_rows(as.data.frame(MCH_spill_days_data), as.data.frame(MCW_spill_days_data))
# 
# # Deschutes - evaluate all fallback move probs
# # Deschutes is wild only
# DES_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "Deschutes River", movements = MC_postovershoot_fallback_movements)
# 
# # Deschutes - get covariate experiences
# DES_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Deschutes River")
# 
# # Deschutes - use a loop to plot all of them
# DES_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
# for (i in 1:nrow(MC_postovershoot_fallback_movements)){
#   # if a movement doesn't have any data, don't plot it
#   filter(DES_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> DES_data_to_plot
#   
#   if(nrow(DES_data_to_plot)>0){
#     DES_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "Deschutes River",
#                                                                             wild_move_prob_array = DES_wild_spilldays_move_prob_array,
#                                                                             
#                                                                             wild_covariate_experiences = DES_wild_covariate_experiences, 
#                                                                             
#                                                                             movements_evaluated = MC_postovershoot_fallback_movements,
#                                                                             spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
#                                                                             from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
#                                                                             plot_title = paste0("Deschutes - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
#     
#     ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays_v3", paste0("DES_compare_fallback_", 
#                                                                                         MC_postovershoot_fallback_movements$dam[i],
#                                                                                         "_spilldays.png")),
#            DES_fallback_spilldays_plots[[i]], height = 8, width = 8)
#     
#   } else {
#     # if it is zero, do nothing
#     
#   }
#   
# 
#   
# }
# 
# # John Day - evaluate all fallback move probs
# # JDR is wild only
# JDR_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "John Day River", movements = MC_postovershoot_fallback_movements)
# 
# # John Day - get covariate experiences
# JDR_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "John Day River")
# 
# # John Day - use a loop to plot all of them
# JDR_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
# for (i in 1:nrow(MC_postovershoot_fallback_movements)){
#   # if a movement doesn't have any data, don't plot it
#   filter(JDR_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> JDR_data_to_plot
#   
#   if(nrow(JDR_data_to_plot)>0){
#   JDR_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "John Day River",
#                                                                           wild_move_prob_array = JDR_wild_spilldays_move_prob_array,
#                                                                           
#                                                                           wild_covariate_experiences = JDR_wild_covariate_experiences, 
#                                                                           
#                                                                           movements_evaluated = MC_postovershoot_fallback_movements,
#                                                                           spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
#                                                                           from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
#                                                                           plot_title = paste0("John Day - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays_v3", paste0("JDR_compare_fallback_", 
#                                                                                       MC_postovershoot_fallback_movements$dam[i],
#                                                                                       "_spilldays.png")),
#          JDR_fallback_spilldays_plots[[i]], height = 8, width = 8)
#   } else {
#     # if it is zero, do nothing
#     
#   }
# }
# 
# 
# # Fifteenmile - evaluate all fallback move probs
# # Fifteenmile is wild only
# FIF_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "Fifteenmile Creek", movements = MC_postovershoot_fallback_movements)
# 
# # Fifteenmile - get covariate experiences
# FIF_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Fifteenmile Creek")
# 
# # Fifteenmile - use a loop to plot all of them
# FIF_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
# for (i in 1:nrow(MC_postovershoot_fallback_movements)){
#   # if a movement doesn't have any data, don't plot it
#   filter(FIF_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> FIF_data_to_plot
#   
#   if(nrow(FIF_data_to_plot)>0){
#   FIF_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "Fifteenmile Creek",
#                                                                           wild_move_prob_array = FIF_wild_spilldays_move_prob_array,
#                                                                           
#                                                                           wild_covariate_experiences = FIF_wild_covariate_experiences, 
#                                                                           
#                                                                           movements_evaluated = MC_postovershoot_fallback_movements,
#                                                                           spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
#                                                                           from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
#                                                                           plot_title = paste0("Fifteenmile - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays_v3", paste0("FIF_compare_fallback_", 
#                                                                                       MC_postovershoot_fallback_movements$dam[i],
#                                                                                       "_spilldays.png")),
#          FIF_fallback_spilldays_plots[[i]], height = 8, width = 8)
#   } else {
#     # if it is zero, do nothing
#     
#   }
# }
# 
# # Umatilla - evaluate all fallback move probs
# UMA_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "Umatilla River", movements = MC_postovershoot_fallback_movements)
# UMA_hatchery_spilldays_move_prob_array <- estimate_spilldays_effect_MCH(origin_select= "Umatilla River", movements = MC_postovershoot_fallback_movements)
# 
# # Umatilla - get covariate experiences
# UMA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Umatilla River")
# UMA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Umatilla River")
# 
# # Umatilla - use a loop to plot all of them
# UMA_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
# for (i in 1:nrow(MC_postovershoot_fallback_movements)){
#   # if a movement doesn't have any data, don't plot it
#   filter(UMA_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> UMA_wild_data_to_plot
#   filter(UMA_hatchery_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> UMA_hatchery_data_to_plot
#   
#   if(nrow(UMA_wild_data_to_plot)>0 | nrow(UMA_hatchery_data_to_plot)>0){
#   UMA_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "Umatilla River",
#                                                                           wild_move_prob_array = UMA_wild_spilldays_move_prob_array,
#                                                                           hatchery_move_prob_array = UMA_hatchery_spilldays_move_prob_array,
#                                                                           wild_covariate_experiences = UMA_wild_covariate_experiences, 
#                                                                           hatchery_covariate_experiences = UMA_hatchery_covariate_experiences,
#                                                                           movements_evaluated = MC_postovershoot_fallback_movements,
#                                                                           spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
#                                                                           from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
#                                                                           plot_title = paste0("Umatilla - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays_v3", paste0("UMA_compare_fallback_", 
#                                                                                       MC_postovershoot_fallback_movements$dam[i],
#                                                                                       "_spilldays.png")),
#          UMA_fallback_spilldays_plots[[i]], height = 8, width = 8)
#   } else {
#     # if it is zero, do nothing
#     
#   }
# }
# 
# # Yakima - evaluate all fallback move probs
# # Yakima is wild only
# YAK_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "Yakima River", movements = MC_postovershoot_fallback_movements)
# 
# # Yakima - get covariate experiences
# YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Yakima River")
# 
# # Yakima - use a loop to plot all of them
# YAK_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
# for (i in 1:nrow(MC_postovershoot_fallback_movements)){
#   # if a movement doesn't have any data, don't plot it
#   filter(YAK_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> YAK_data_to_plot
#   
#   if(nrow(YAK_data_to_plot)>0){
#   YAK_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "Yakima River",
#                                                                           wild_move_prob_array = YAK_wild_spilldays_move_prob_array,
#                                                                           
#                                                                           wild_covariate_experiences = YAK_wild_covariate_experiences, 
#                                                                           
#                                                                           movements_evaluated = MC_postovershoot_fallback_movements,
#                                                                           spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
#                                                                           from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
#                                                                           plot_title = paste0("Yakima - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays_v3", paste0("YAK_compare_fallback_", 
#                                                                                       MC_postovershoot_fallback_movements$dam[i],
#                                                                                       "_spilldays.png")),
#          YAK_fallback_spilldays_plots[[i]], height = 8, width = 8)
#   } else {
#     # if it is zero, do nothing
#     
#   }
# }
# 
# # Walla Walla - evaluate all fallback move probs
# WAWA_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "Walla Walla River", movements = MC_postovershoot_fallback_movements)
# WAWA_hatchery_spilldays_move_prob_array <- estimate_spilldays_effect_MCH(origin_select= "Walla Walla River", movements = MC_postovershoot_fallback_movements)
# 
# # Walla Walla - get covariate experiences
# WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Walla Walla River")
# WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Walla Walla River")
# 
# # Walla Walla - use a loop to plot all of them
# WAWA_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
# for (i in 1:nrow(MC_postovershoot_fallback_movements)){
#   # if a movement doesn't have any data, don't plot it
#   filter(WAWA_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> WAWA_wild_data_to_plot
#   filter(WAWA_hatchery_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> WAWA_hatchery_data_to_plot
#   
#   if(nrow(WAWA_wild_data_to_plot)>0 | nrow(WAWA_hatchery_data_to_plot)>0){
#   WAWA_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "Walla Walla River",
#                                                                            wild_move_prob_array = WAWA_wild_spilldays_move_prob_array,
#                                                                            hatchery_move_prob_array = WAWA_hatchery_spilldays_move_prob_array,
#                                                                            wild_covariate_experiences = WAWA_wild_covariate_experiences, 
#                                                                            hatchery_covariate_experiences = WAWA_hatchery_covariate_experiences,
#                                                                            movements_evaluated = MC_postovershoot_fallback_movements,
#                                                                            spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
#                                                                            from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
#                                                                            plot_title = paste0("Walla Walla - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays_v3", paste0("WAWA_compare_fallback_", 
#                                                                                       MC_postovershoot_fallback_movements$dam[i],
#                                                                                       "_spilldays.png")),
#          WAWA_fallback_spilldays_plots[[i]], height = 8, width = 8)
#   } else {
#     # if it is zero, do nothing
#     
#   }
# }
# 
# 
# #### Snake River ####
# # only one post-overshoot fallback: Tucannon
# 
# SRH_spill_days_data <- SRH_envir$data$winter_spill_days_data
# SRW_spill_days_data <- SRW_envir$data$winter_spill_days_data
# SR_spill_days_data <- bind_rows(as.data.frame(SRH_spill_days_data), as.data.frame(SRW_spill_days_data))
# 
# # Tucannon River - evaluate all fallback move probs
# TUC_wild_spilldays_move_prob_array <- estimate_spilldays_effect_SRW(origin_select= "Tucannon River", movements = SR_postovershoot_fallback_movements)
# TUC_hatchery_spilldays_move_prob_array <- estimate_spilldays_effect_SRH(origin_select= "Tucannon River", movements = SR_postovershoot_fallback_movements)
# 
# # Tucannon River - get covariate experiences
# TUC_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Tucannon River")
# TUC_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Tucannon River")
# 
# # Tucannon River - use a loop to plot all of them
# TUC_fallback_spilldays_plots <- vector(mode = "list", length = nrow(SR_postovershoot_fallback_movements))
# for (i in 1:nrow(SR_postovershoot_fallback_movements)){
#   # if a movement doesn't have any data, don't plot it
#   filter(TUC_wild_covariate_experiences, state == SR_postovershoot_fallback_movements$from[i]) -> TUC_wild_data_to_plot
#   filter(TUC_hatchery_covariate_experiences, state == SR_postovershoot_fallback_movements$from[i]) -> TUC_hatchery_data_to_plot
#   
#   if(nrow(TUC_wild_data_to_plot)>0 | nrow(TUC_hatchery_data_to_plot)>0){
#   TUC_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "Tucannon River",
#                                                                           wild_move_prob_array = TUC_wild_spilldays_move_prob_array,
#                                                                           hatchery_move_prob_array = TUC_hatchery_spilldays_move_prob_array,
#                                                                           wild_covariate_experiences = TUC_wild_covariate_experiences, 
#                                                                           hatchery_covariate_experiences = TUC_hatchery_covariate_experiences,
#                                                                           movements_evaluated = SR_postovershoot_fallback_movements,
#                                                                           spill_predict = seq(0, max(SR_spill_days_data[,SR_postovershoot_fallback_movements$from[i]]),length = 100),
#                                                                           from = SR_postovershoot_fallback_movements$from[i], to = SR_postovershoot_fallback_movements$to[i], 
#                                                                           plot_title = paste0("Tucannon River - post-overshoot fallback at ", SR_postovershoot_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays_v3", paste0("TUC_compare_fallback_", 
#                                                                                       SR_postovershoot_fallback_movements$dam[i],
#                                                                                       "_spilldays.png")),
#          TUC_fallback_spilldays_plots[[i]], height = 8, width = 8)
#   } else {
#     # if it is zero, do nothing
#     
#   }
# }
# 
#### Plot fit to data for spill window, get goodness of fit ####

plot_compare_rear_spillwindow_effect_fit_to_data <- function(origin_select,
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
      dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_spill_move_prob_quantiles
    
    wild_spill_move_prob_quantiles -> rear_spill_move_prob_quantiles
    
    # convert back to kcfs
    rear_spill_move_prob_quantiles %>% 
      mutate(spill_actual = spill*100) -> rear_spill_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    
    wild_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
  } else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    # If wild is NA, run hatchery only
    # modify covariate experiences to show what the next state is
    hatchery_covariate_experiences %>% 
      mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> hatchery_covariate_experiences
    
    hatchery_spill_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(hatchery_spill_move_prob) <- paste0("iter", 1:niter) 
    hatchery_spill_move_prob$spill <- spill_predict
    
    hatchery_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_spill_move_prob_quantiles
    
    hatchery_spill_move_prob_quantiles -> rear_spill_move_prob_quantiles
    
    # convert back to kcfs
    rear_spill_move_prob_quantiles %>% 
      mutate(spill_actual = spill*100) -> rear_spill_move_prob_quantiles
    
    # get data organized for rug plot
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    
    hatchery_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
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
      dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_spill_move_prob_quantiles
    
    hatchery_spill_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(hatchery_spill_move_prob) <- paste0("iter", 1:niter) 
    hatchery_spill_move_prob$spill <- spill_predict
    
    hatchery_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_spill_move_prob_quantiles
    
    wild_spill_move_prob_quantiles %>% 
      bind_rows(., hatchery_spill_move_prob_quantiles) -> rear_spill_move_prob_quantiles
    
    # convert back to kcfs
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
    mutate(spill_actual = spill_window*100) -> covariate_experiences
  
  # Remove any instances for rug plot of fish that would have received winter spill days covariate
  covariate_experiences %>% 
    filter(winter_post_overshoot_vector != 1) -> covariate_experiences
  
  # now, create DF to show fit to data
  # run this only if the fish actually visited that state
  if(nrow(covariate_experiences) != 0) {
    # create a new DF to show bubbles of points
    covariate_experiences %>% 
      mutate(spill_round = round_any(spill_actual, 5)) %>% 
      group_by(spill_round, rear) %>% 
      dplyr::summarise(total = n()) -> spill_rear_total_counts
    
    covariate_experiences %>% 
      mutate(spill_round = round_any(spill_actual, 5)) %>% 
      group_by(spill_round, rear, next_state) %>% 
      dplyr::summarise(count = n()) %>% 
      ungroup() %>% 
      complete(spill_round, nesting(rear, next_state), fill = list(count = 0)) -> spill_rear_move_counts
    
    spill_rear_move_counts %>% 
      left_join(spill_rear_total_counts, by = c("spill_round", "rear")) %>% 
      mutate(prop = count/total) -> spill_move_props
    
  } else {
    spill_move_props <- data.frame(spill_round = as.numeric(NA), rear = as.character(NA), 
                                   next_state = as.character(NA), count  = as.numeric(NA), 
                                   total = as.numeric(NA), prop = as.numeric(NA))
    
  }
  
  # Calculate the credible interval coverage
  # also calculate the average bias (from model estimated median)
  rear_spill_move_prob_quantiles %>% 
    mutate(spill_round = round(spill_actual, 0)) %>% 
    # only keep the closest prediction per rounded spill
    mutate(spill_diff = abs(spill_round-spill_actual)) %>% 
    group_by(spill_round) %>% 
    top_n(-1, spill_diff) %>% 
    left_join(subset(spill_move_props, next_state == to), by = c("spill_round", "rear")) %>% 
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
    dplyr::summarise(points_within_CI = sum(CI_coverage_weight, na.rm = T),
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
               aes(x = spill_round, y = prop, color = rear, size = total), alpha = 0.5,
               inherit.aes = FALSE) +
    geom_point(data = subset(spill_move_props, next_state == to & rear == "hatchery"),
               aes(x = spill_round, y = prop, color = rear, size = total), alpha = 0.5,
               inherit.aes = FALSE) +
    scale_y_continuous(lim = c(-0.01, 1.01), expand = c(0,0)) +
    scale_x_continuous(lim = c(-1,max(rear_spill_move_prob_quantiles$spill_actual)+1), expand = c(0,0)) +
    scale_color_manual(values = rear_colors) +
    scale_fill_manual(values = rear_colors) +
    xlab("Spill (kcfs)") +
    ylab("Movement probability") +
    ggtitle(plot_title)
  
  return(list(plot = rear_spill_move_prob_plot,
              goodness_of_fit = goodness_of_fit))
}


# #### Plot rear type comparison fallback spill window for all ####
# mainstem_fallback_movements <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL", "ICH", "LGR"),
#                                           from = c(2,3,4,5,6,7,8,9), to = c(1,2,3,4,5,6,3,8))
# 
# #### Upper Columbia ####
# UCH_spill_window_data <- UCH_envir$data$spill_window_data
# UCW_spill_window_data <- UCW_envir$data$spill_window_data
# UC_spill_window_data <- bind_rows(as.data.frame(UCH_spill_window_data), as.data.frame(UCW_spill_window_data))
# 
# # Wenatchee - evaluate all fallback move probs
# WEN_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCW(origin_select= "Wenatchee River", movements = mainstem_fallback_movements)
# WEN_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCH(origin_select= "Wenatchee River", movements = mainstem_fallback_movements)
# 
# # Wenatchee - get covariate experiences
# WEN_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Wenatchee River")
# WEN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Wenatchee River")
# 
# # Wenatchee - use a loop to plot all of them
# WEN_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   WEN_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect_fit_to_data(origin_select = "Wenatchee River",
#                                                                               wild_move_prob_array = WEN_wild_spillwindow_move_prob_array,
#                                                                               hatchery_move_prob_array = WEN_hatchery_spillwindow_move_prob_array,
#                                                                               wild_covariate_experiences = WEN_wild_covariate_experiences, 
#                                                                               hatchery_covariate_experiences = WEN_hatchery_covariate_experiences,
#                                                                               movements_evaluated = mainstem_fallback_movements,
#                                                                               spill_predict = seq(0, max(UC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                               from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                               plot_title = paste0("Wenatchee - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "fit_to_data", "spillwindow_v3", paste0("WEN_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          WEN_fallback_spillwindow_plots[[i]]$plot, height = 8, width = 8)
#   
# }
# 
# # Entiat - evaluate all fallback move probs
# # Only wild for entiat
# ENT_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCW(origin_select= "Entiat River", movements = mainstem_fallback_movements)
# 
# # Entiat - get covariate experiences
# ENT_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Entiat River")
# 
# # Entiat - use a loop to plot all of them
# ENT_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   ENT_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect_fit_to_data(origin_select = "Entiat River",
#                                                                               wild_move_prob_array = ENT_wild_spillwindow_move_prob_array,
#                                                                               wild_covariate_experiences = ENT_wild_covariate_experiences, 
#                                                                               movements_evaluated = mainstem_fallback_movements,
#                                                                               spill_predict = seq(0, max(UC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                               from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                               plot_title = paste0("Entiat - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "fit_to_data", "spillwindow_v3", paste0("ENT_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          ENT_fallback_spillwindow_plots[[i]]$plot, height = 8, width = 8)
#   
# }
# 
# # Methow - evaluate all fallback move probs
# MET_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCW(origin_select= "Methow River", movements = mainstem_fallback_movements)
# MET_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCH(origin_select= "Methow River", movements = mainstem_fallback_movements)
# 
# # Methow - get covariate experiences
# MET_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Methow River")
# MET_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Methow River")
# 
# # Methow - use a loop to plot all of them
# MET_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   MET_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect_fit_to_data(origin_select = "Methow River",
#                                                                               wild_move_prob_array = MET_wild_spillwindow_move_prob_array,
#                                                                               hatchery_move_prob_array = MET_hatchery_spillwindow_move_prob_array,
#                                                                               wild_covariate_experiences = MET_wild_covariate_experiences, 
#                                                                               hatchery_covariate_experiences = MET_hatchery_covariate_experiences,
#                                                                               movements_evaluated = mainstem_fallback_movements,
#                                                                               spill_predict = seq(0, max(UC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                               from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                               plot_title = paste0("Methow - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "fit_to_data", "spillwindow_v3", paste0("MET_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          MET_fallback_spillwindow_plots[[i]]$plot, height = 8, width = 8)
#   
# }
# 
# # Okanoagn - evaluate all fallback move probs
# # Only hatchery for Okanogan
# OKA_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCH(origin_select= "Okanogan River", movements = mainstem_fallback_movements)
# 
# # Okanogan - get covariate experiences
# OKA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Okanogan River")
# 
# # Okanogan - use a loop to plot all of them
# OKA_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   OKA_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect_fit_to_data(origin_select = "Okanogan River",
#                                                                               hatchery_move_prob_array = OKA_hatchery_spillwindow_move_prob_array,
#                                                                               hatchery_covariate_experiences = OKA_hatchery_covariate_experiences,
#                                                                               movements_evaluated = mainstem_fallback_movements,
#                                                                               spill_predict = seq(0, max(UC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                               from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                               plot_title = paste0("Okanogan - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "fit_to_data", "spillwindow_v3", paste0("OKA_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          OKA_fallback_spillwindow_plots[[i]]$plot, height = 8, width = 8)
#   
# }
# 
# #### Middle Columbia ####
# MCH_spill_window_data <- MCH_envir$data$spill_window_data
# MCW_spill_window_data <- MCW_envir$data$spill_window_data
# MC_spill_window_data <- bind_rows(as.data.frame(MCH_spill_window_data), as.data.frame(MCW_spill_window_data))
# 
# # Deschutes - evaluate all fallback move probs
# # Deschutes is wild only
# DES_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "Deschutes River", movements = mainstem_fallback_movements)
# 
# # Deschutes - get covariate experiences
# DES_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Deschutes River")
# 
# # Deschutes - use a loop to plot all of them
# DES_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   DES_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect_fit_to_data(origin_select = "Deschutes River",
#                                                                               wild_move_prob_array = DES_wild_spillwindow_move_prob_array,
#                                                                               
#                                                                               wild_covariate_experiences = DES_wild_covariate_experiences, 
#                                                                               
#                                                                               movements_evaluated = mainstem_fallback_movements,
#                                                                               spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                               from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                               plot_title = paste0("Deschutes - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "fit_to_data", "spillwindow_v3", paste0("DES_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          DES_fallback_spillwindow_plots[[i]]$plot, height = 8, width = 8)
#   
# }
# 
# # John Day - evaluate all fallback move probs
# # JDR is wild only
# JDR_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "John Day River", movements = mainstem_fallback_movements)
# 
# # John Day - get covariate experiences
# JDR_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "John Day River")
# 
# # John Day - use a loop to plot all of them
# JDR_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   JDR_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect_fit_to_data(origin_select = "John Day River",
#                                                                               wild_move_prob_array = JDR_wild_spillwindow_move_prob_array,
#                                                                               
#                                                                               wild_covariate_experiences = JDR_wild_covariate_experiences, 
#                                                                               
#                                                                               movements_evaluated = mainstem_fallback_movements,
#                                                                               spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                               from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                               plot_title = paste0("John Day - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "fit_to_data", "spillwindow_v3", paste0("JDR_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          JDR_fallback_spillwindow_plots[[i]]$plot, height = 8, width = 8)
#   
# }
# 
# 
# # Fifteenmile - evaluate all fallback move probs
# # Fifteenmile is wild only
# FIF_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "Fifteenmile Creek", movements = mainstem_fallback_movements)
# 
# # Fifteenmile - get covariate experiences
# FIF_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Fifteenmile Creek")
# 
# # Fifteenmile - use a loop to plot all of them
# FIF_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   FIF_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect_fit_to_data(origin_select = "Fifteenmile Creek",
#                                                                               wild_move_prob_array = FIF_wild_spillwindow_move_prob_array,
#                                                                               
#                                                                               wild_covariate_experiences = FIF_wild_covariate_experiences, 
#                                                                               
#                                                                               movements_evaluated = mainstem_fallback_movements,
#                                                                               spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                               from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                               plot_title = paste0("Fifteenmile - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "fit_to_data", "spillwindow_v3", paste0("FIF_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          FIF_fallback_spillwindow_plots[[i]]$plot, height = 8, width = 8)
#   
# }
# 
# # Umatilla - evaluate all fallback move probs
# UMA_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "Umatilla River", movements = mainstem_fallback_movements)
# UMA_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCH(origin_select= "Umatilla River", movements = mainstem_fallback_movements)
# 
# # Umatilla - get covariate experiences
# UMA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Umatilla River")
# UMA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Umatilla River")
# 
# # Umatilla - use a loop to plot all of them
# UMA_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   UMA_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect_fit_to_data(origin_select = "Umatilla River",
#                                                                               wild_move_prob_array = UMA_wild_spillwindow_move_prob_array,
#                                                                               hatchery_move_prob_array = UMA_hatchery_spillwindow_move_prob_array,
#                                                                               wild_covariate_experiences = UMA_wild_covariate_experiences, 
#                                                                               hatchery_covariate_experiences = UMA_hatchery_covariate_experiences,
#                                                                               movements_evaluated = mainstem_fallback_movements,
#                                                                               spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                               from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                               plot_title = paste0("Umatilla - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "fit_to_data", "spillwindow_v3", paste0("UMA_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          UMA_fallback_spillwindow_plots[[i]]$plot, height = 8, width = 8)
#   
# }
# 
# # Yakima - evaluate all fallback move probs
# # Yakima is wild only
# YAK_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "Yakima River", movements = mainstem_fallback_movements)
# 
# # Yakima - get covariate experiences
# YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Yakima River")
# 
# # Yakima - use a loop to plot all of them
# YAK_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   YAK_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect_fit_to_data(origin_select = "Yakima River",
#                                                                               wild_move_prob_array = YAK_wild_spillwindow_move_prob_array,
#                                                                               
#                                                                               wild_covariate_experiences = YAK_wild_covariate_experiences, 
#                                                                               
#                                                                               movements_evaluated = mainstem_fallback_movements,
#                                                                               spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                               from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                               plot_title = paste0("Yakima - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "fit_to_data", "spillwindow_v3", paste0("YAK_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          YAK_fallback_spillwindow_plots[[i]]$plot, height = 8, width = 8)
#   
# }
# 
# # Walla Walla - evaluate all fallback move probs
# WAWA_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "Walla Walla River", movements = mainstem_fallback_movements)
# WAWA_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCH(origin_select= "Walla Walla River", movements = mainstem_fallback_movements)
# 
# # Walla Walla - get covariate experiences
# WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Walla Walla River")
# WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Walla Walla River")
# 
# # Walla Walla - use a loop to plot all of them
# WAWA_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   WAWA_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect_fit_to_data(origin_select = "Walla Walla River",
#                                                                                wild_move_prob_array = WAWA_wild_spillwindow_move_prob_array,
#                                                                                hatchery_move_prob_array = WAWA_hatchery_spillwindow_move_prob_array,
#                                                                                wild_covariate_experiences = WAWA_wild_covariate_experiences, 
#                                                                                hatchery_covariate_experiences = WAWA_hatchery_covariate_experiences,
#                                                                                movements_evaluated = mainstem_fallback_movements,
#                                                                                spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                                from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                                plot_title = paste0("Walla Walla - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "fit_to_data", "spillwindow_v3", paste0("WAWA_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          WAWA_fallback_spillwindow_plots[[i]]$plot, height = 8, width = 8)
#   
# }
# 
# #### Snake River ####
# SRH_spill_window_data <- SRH_envir$data$spill_window_data
# SRW_spill_window_data <- SRW_envir$data$spill_window_data
# SR_spill_window_data <- bind_rows(as.data.frame(SRH_spill_window_data), as.data.frame(SRW_spill_window_data))
# 
# # Asotin Creek - evaluate all fallback move probs
# # Asotin is wild only
# ASO_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Asotin Creek", movements = mainstem_fallback_movements)
# 
# # Asotin Creek - get covariate experiences
# ASO_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Asotin Creek")
# 
# # Asotin Creek - use a loop to plot all of them
# ASO_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   ASO_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect_fit_to_data(origin_select = "Asotin Creek",
#                                                                               wild_move_prob_array = ASO_wild_spillwindow_move_prob_array,
#                                                                               
#                                                                               wild_covariate_experiences = ASO_wild_covariate_experiences, 
#                                                                               
#                                                                               movements_evaluated = mainstem_fallback_movements,
#                                                                               spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                               from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                               plot_title = paste0("Asotin Creek - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "fit_to_data", "spillwindow_v3", paste0("ASO_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          ASO_fallback_spillwindow_plots[[i]]$plot, height = 8, width = 8)
#   
# }
# 
# # Clearwater River - evaluate all fallback move probs
# CLE_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Clearwater River", movements = mainstem_fallback_movements)
# CLE_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select= "Clearwater River", movements = mainstem_fallback_movements)
# 
# # Clearwater River - get covariate experiences
# CLE_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Clearwater River")
# CLE_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Clearwater River")
# 
# # Clearwater River - use a loop to plot all of them
# CLE_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   CLE_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect_fit_to_data(origin_select = "Clearwater River",
#                                                                               wild_move_prob_array = CLE_wild_spillwindow_move_prob_array,
#                                                                               hatchery_move_prob_array = CLE_hatchery_spillwindow_move_prob_array,
#                                                                               wild_covariate_experiences = CLE_wild_covariate_experiences, 
#                                                                               hatchery_covariate_experiences = CLE_hatchery_covariate_experiences,
#                                                                               movements_evaluated = mainstem_fallback_movements,
#                                                                               spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                               from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                               plot_title = paste0("Clearwater River - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "fit_to_data", "spillwindow_v3", paste0("CLE_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          CLE_fallback_spillwindow_plots[[i]]$plot, height = 8, width = 8)
#   
# }
# 
# # Grande Ronde River - evaluate all fallback move probs
# GR_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Grande Ronde River", movements = mainstem_fallback_movements)
# GR_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select= "Grande Ronde River", movements = mainstem_fallback_movements)
# 
# # Grande Ronde River - get covariate experiences
# GR_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Grande Ronde River")
# GR_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Grande Ronde River")
# 
# # Grande Ronde River - use a loop to plot all of them
# GR_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   GR_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect_fit_to_data(origin_select = "Grande Ronde River",
#                                                                              wild_move_prob_array = GR_wild_spillwindow_move_prob_array,
#                                                                              hatchery_move_prob_array = GR_hatchery_spillwindow_move_prob_array,
#                                                                              wild_covariate_experiences = GR_wild_covariate_experiences, 
#                                                                              hatchery_covariate_experiences = GR_hatchery_covariate_experiences,
#                                                                              movements_evaluated = mainstem_fallback_movements,
#                                                                              spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                              from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                              plot_title = paste0("Grande Ronde River - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "fit_to_data", "spillwindow_v3", paste0("GR_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          GR_fallback_spillwindow_plots[[i]]$plot, height = 8, width = 8)
#   
# }
# 
# # Imnaha River - evaluate all fallback move probs
# IMN_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Imnaha River", movements = mainstem_fallback_movements)
# IMN_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select= "Imnaha River", movements = mainstem_fallback_movements)
# 
# # Imnaha River - get covariate experiences
# IMN_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Imnaha River")
# IMN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Imnaha River")
# 
# # Imnaha River - use a loop to plot all of them
# IMN_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   IMN_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect_fit_to_data(origin_select = "Imnaha River",
#                                                                               wild_move_prob_array = IMN_wild_spillwindow_move_prob_array,
#                                                                               hatchery_move_prob_array = IMN_hatchery_spillwindow_move_prob_array,
#                                                                               wild_covariate_experiences = IMN_wild_covariate_experiences, 
#                                                                               hatchery_covariate_experiences = IMN_hatchery_covariate_experiences,
#                                                                               movements_evaluated = mainstem_fallback_movements,
#                                                                               spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                               from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                               plot_title = paste0("Imnaha River - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "fit_to_data", "spillwindow_v3", paste0("IMN_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          IMN_fallback_spillwindow_plots[[i]]$plot, height = 8, width = 8)
#   
# }
# 
# # Salmon River - evaluate all fallback move probs
# SAL_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Salmon River", movements = mainstem_fallback_movements)
# SAL_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select= "Salmon River", movements = mainstem_fallback_movements)
# 
# # Salmon River - get covariate experiences
# SAL_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Salmon River")
# SAL_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Salmon River")
# 
# # Salmon River - use a loop to plot all of them
# SAL_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   SAL_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect_fit_to_data(origin_select = "Salmon River",
#                                                                               wild_move_prob_array = SAL_wild_spillwindow_move_prob_array,
#                                                                               hatchery_move_prob_array = SAL_hatchery_spillwindow_move_prob_array,
#                                                                               wild_covariate_experiences = SAL_wild_covariate_experiences, 
#                                                                               hatchery_covariate_experiences = SAL_hatchery_covariate_experiences,
#                                                                               movements_evaluated = mainstem_fallback_movements,
#                                                                               spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                               from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                               plot_title = paste0("Salmon River - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "fit_to_data", "spillwindow_v3", paste0("SAL_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          SAL_fallback_spillwindow_plots[[i]]$plot, height = 8, width = 8)
#   
# }
# # Tucannon River - evaluate all fallback move probs
# TUC_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Tucannon River", movements = mainstem_fallback_movements)
# TUC_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select= "Tucannon River", movements = mainstem_fallback_movements)
# 
# # Tucannon River - get covariate experiences
# TUC_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Tucannon River")
# TUC_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Tucannon River")
# 
# # Tucannon River - use a loop to plot all of them
# TUC_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   TUC_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spillwindow_effect_fit_to_data(origin_select = "Tucannon River",
#                                                                               wild_move_prob_array = TUC_wild_spillwindow_move_prob_array,
#                                                                               hatchery_move_prob_array = TUC_hatchery_spillwindow_move_prob_array,
#                                                                               wild_covariate_experiences = TUC_wild_covariate_experiences, 
#                                                                               hatchery_covariate_experiences = TUC_hatchery_covariate_experiences,
#                                                                               movements_evaluated = mainstem_fallback_movements,
#                                                                               spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                               from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                               plot_title = paste0("Tucannon River - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output", "fit_to_data", "spillwindow_v3", paste0("TUC_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          TUC_fallback_spillwindow_plots[[i]]$plot, height = 8, width = 8)
#   
# }
# 

#### Plot fit to data for spill days, get goodness of fit #### 

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
      dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
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
      dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
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
      dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_spill_move_prob_quantiles
    
    hatchery_spill_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(hatchery_spill_move_prob) <- paste0("iter", 1:niter) 
    hatchery_spill_move_prob$spill <- spill_predict
    
    hatchery_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
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
      group_by(spill_actual, rear) %>% 
      dplyr::summarise(total = n()) -> spill_rear_total_counts
    
    covariate_experiences %>% 
      group_by(spill_actual, rear, next_state) %>% 
      dplyr::summarise(count = n()) %>% 
      ungroup() %>% 
      complete(spill_actual, nesting(rear, next_state), fill = list(count = 0)) -> spill_rear_move_counts
    
    spill_rear_move_counts %>% 
      left_join(spill_rear_total_counts, by = c("spill_actual", "rear")) %>% 
      mutate(prop = count/total) -> spill_move_props
    
  } else {
    spill_move_props <- data.frame(spill_actual = as.numeric(NA), rear = as.character(NA), 
                                  next_state = as.character(NA), count  = as.numeric(NA), 
                                  total = as.numeric(NA), prop = as.numeric(NA))
  }
  
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
    dplyr::summarise(points_within_CI = sum(CI_coverage_weight, na.rm = T),
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


#### Plot rear type comparison fallback spill days for post-overshot fallback movements ####
UC_postovershoot_fallback_movements <- data.frame(dam = c("RRE", "WEL"),
                                                  from = c(6,7), to = c(5,6))

MC_postovershoot_fallback_movements <- data.frame(dam = c("MCN", "PRA", "ICH"),
                                                  from = c(3,4,8), to = c(2,3,3))

SR_postovershoot_fallback_movements <- data.frame(dam = c("LGR"),
                                                  from = c(9), to = c(8))

#### Upper Columbia ####
# only Entiat and Wenatchee

UCH_spill_days_data <- UCH_envir$data$winter_spill_days_data
UCW_spill_days_data <- UCW_envir$data$winter_spill_days_data
UC_spill_days_data <- bind_rows(as.data.frame(UCH_spill_days_data), as.data.frame(UCW_spill_days_data))

# Wenatchee - evaluate all fallback move probs
WEN_wild_spilldays_move_prob_array <- estimate_spilldays_effect_UCW(origin_select= "Wenatchee River", movements = UC_postovershoot_fallback_movements)
WEN_hatchery_spilldays_move_prob_array <- estimate_spilldays_effect_UCH(origin_select= "Wenatchee River", movements = UC_postovershoot_fallback_movements)

# Wenatchee - get covariate experiences
WEN_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Wenatchee River")
WEN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Wenatchee River")

# Wenatchee - use a loop to plot all of them
WEN_fallback_spilldays_plots <- vector(mode = "list", length = nrow(UC_postovershoot_fallback_movements))
for (i in 1:nrow(UC_postovershoot_fallback_movements)){
  # if a movement doesn't have any data, don't plot it
  filter(WEN_wild_covariate_experiences, state == UC_postovershoot_fallback_movements$from[i]) -> WEN_data_to_plot
  
  if(nrow(WEN_data_to_plot)>0){
    WEN_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect_fit_to_data(origin_select = "Wenatchee River",
                                                                            wild_move_prob_array = WEN_wild_spilldays_move_prob_array,
                                                                            hatchery_move_prob_array = WEN_hatchery_spilldays_move_prob_array,
                                                                            wild_covariate_experiences = WEN_wild_covariate_experiences, 
                                                                            hatchery_covariate_experiences = WEN_hatchery_covariate_experiences,
                                                                            movements_evaluated = UC_postovershoot_fallback_movements,
                                                                            spill_predict = seq(0, max(UC_spill_days_data[,UC_postovershoot_fallback_movements$from[i]]),length = 100),
                                                                            from = UC_postovershoot_fallback_movements$from[i], to = UC_postovershoot_fallback_movements$to[i], 
                                                                            plot_title = paste0("Wenatchee - post-overshoot fallback at ", UC_postovershoot_fallback_movements$dam[i]))
    
    ggsave(here::here("stan_actual", "output", "fit_to_data", "spilldays_v3", paste0("WEN_compare_fallback_", 
                                                                                        UC_postovershoot_fallback_movements$dam[i],
                                                                                        "_spilldays.png")),
           WEN_fallback_spilldays_plots[[i]]$plot, height = 8, width = 8)
  } else {
    # if it is zero, do nothing
    
  }
  
}

# Entiat - evaluate all fallback move probs
# Only wild for entiat
ENT_wild_spilldays_move_prob_array <- estimate_spilldays_effect_UCW(origin_select= "Entiat River", movements = UC_postovershoot_fallback_movements)

# Entiat - get covariate experiences
ENT_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Entiat River")

# Entiat - use a loop to plot all of them
ENT_fallback_spilldays_plots <- vector(mode = "list", length = nrow(UC_postovershoot_fallback_movements))
for (i in 1:nrow(UC_postovershoot_fallback_movements)){
  # if a movement doesn't have any data, don't plot it
  filter(ENT_wild_covariate_experiences, state == UC_postovershoot_fallback_movements$from[i]) -> ENT_data_to_plot
  
  if(nrow(ENT_data_to_plot)>0){
    
    ENT_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect_fit_to_data(origin_select = "Entiat River",
                                                                            wild_move_prob_array = ENT_wild_spilldays_move_prob_array,
                                                                            wild_covariate_experiences = ENT_wild_covariate_experiences, 
                                                                            movements_evaluated = UC_postovershoot_fallback_movements,
                                                                            spill_predict = seq(0, max(UC_spill_days_data[,UC_postovershoot_fallback_movements$from[i]]),length = 100),
                                                                            from = UC_postovershoot_fallback_movements$from[i], to = UC_postovershoot_fallback_movements$to[i], 
                                                                            plot_title = paste0("Entiat - post-overshoot fallback at ", UC_postovershoot_fallback_movements$dam[i]))
    
    ggsave(here::here("stan_actual", "output", "fit_to_data", "spilldays_v3", paste0("ENT_compare_fallback_", 
                                                                                        UC_postovershoot_fallback_movements$dam[i],
                                                                                        "_spilldays_v3.png")),
           ENT_fallback_spilldays_plots[[i]]$plot, height = 8, width = 8)
  } else {
    # if it is zero, do nothing
    
  }
}



#### Middle Columbia ####
MCH_spill_days_data <- MCH_envir$data$winter_spill_days_data
MCW_spill_days_data <- MCW_envir$data$winter_spill_days_data
MC_spill_days_data <- bind_rows(as.data.frame(MCH_spill_days_data), as.data.frame(MCW_spill_days_data))

# Deschutes - evaluate all fallback move probs
# Deschutes is wild only
DES_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "Deschutes River", movements = MC_postovershoot_fallback_movements)

# Deschutes - get covariate experiences
DES_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Deschutes River")

# Deschutes - use a loop to plot all of them
DES_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
for (i in 1:nrow(MC_postovershoot_fallback_movements)){
  # if a movement doesn't have any data, don't plot it
  filter(DES_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> DES_data_to_plot
  
  if(nrow(DES_data_to_plot)>0){
    DES_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect_fit_to_data(origin_select = "Deschutes River",
                                                                            wild_move_prob_array = DES_wild_spilldays_move_prob_array,
                                                                            
                                                                            wild_covariate_experiences = DES_wild_covariate_experiences, 
                                                                            
                                                                            movements_evaluated = MC_postovershoot_fallback_movements,
                                                                            spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
                                                                            from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
                                                                            plot_title = paste0("Deschutes - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
    
    ggsave(here::here("stan_actual", "output", "fit_to_data", "spilldays_v3", paste0("DES_compare_fallback_", 
                                                                                        MC_postovershoot_fallback_movements$dam[i],
                                                                                        "_spilldays.png")),
           DES_fallback_spilldays_plots[[i]]$plot, height = 8, width = 8)
    
  } else {
    # if it is zero, do nothing
    
  }
  
  
  
}

# John Day - evaluate all fallback move probs
# JDR is wild only
JDR_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "John Day River", movements = MC_postovershoot_fallback_movements)

# John Day - get covariate experiences
JDR_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "John Day River")

# John Day - use a loop to plot all of them
JDR_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
for (i in 1:nrow(MC_postovershoot_fallback_movements)){
  # if a movement doesn't have any data, don't plot it
  filter(JDR_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> JDR_data_to_plot
  
  if(nrow(JDR_data_to_plot)>0){
    JDR_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect_fit_to_data(origin_select = "John Day River",
                                                                            wild_move_prob_array = JDR_wild_spilldays_move_prob_array,
                                                                            
                                                                            wild_covariate_experiences = JDR_wild_covariate_experiences, 
                                                                            
                                                                            movements_evaluated = MC_postovershoot_fallback_movements,
                                                                            spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
                                                                            from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
                                                                            plot_title = paste0("John Day - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
    
    ggsave(here::here("stan_actual", "output", "fit_to_data", "spilldays_v3", paste0("JDR_compare_fallback_", 
                                                                                        MC_postovershoot_fallback_movements$dam[i],
                                                                                        "_spilldays.png")),
           JDR_fallback_spilldays_plots[[i]]$plot, height = 8, width = 8)
  } else {
    # if it is zero, do nothing
    
  }
}


# Fifteenmile - evaluate all fallback move probs
# Fifteenmile is wild only
FIF_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "Fifteenmile Creek", movements = MC_postovershoot_fallback_movements)

# Fifteenmile - get covariate experiences
FIF_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Fifteenmile Creek")

# Fifteenmile - use a loop to plot all of them
FIF_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
for (i in 1:nrow(MC_postovershoot_fallback_movements)){
  # if a movement doesn't have any data, don't plot it
  filter(FIF_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> FIF_data_to_plot
  
  if(nrow(FIF_data_to_plot)>0){
    FIF_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect_fit_to_data(origin_select = "Fifteenmile Creek",
                                                                            wild_move_prob_array = FIF_wild_spilldays_move_prob_array,
                                                                            
                                                                            wild_covariate_experiences = FIF_wild_covariate_experiences, 
                                                                            
                                                                            movements_evaluated = MC_postovershoot_fallback_movements,
                                                                            spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
                                                                            from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
                                                                            plot_title = paste0("Fifteenmile - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
    
    ggsave(here::here("stan_actual", "output", "fit_to_data", "spilldays_v3", paste0("FIF_compare_fallback_", 
                                                                                        MC_postovershoot_fallback_movements$dam[i],
                                                                                        "_spilldays.png")),
           FIF_fallback_spilldays_plots[[i]]$plot, height = 8, width = 8)
  } else {
    # if it is zero, do nothing
    
  }
}

# Umatilla - evaluate all fallback move probs
UMA_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "Umatilla River", movements = MC_postovershoot_fallback_movements)
UMA_hatchery_spilldays_move_prob_array <- estimate_spilldays_effect_MCH(origin_select= "Umatilla River", movements = MC_postovershoot_fallback_movements)

# Umatilla - get covariate experiences
UMA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Umatilla River")
UMA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Umatilla River")

# Umatilla - use a loop to plot all of them
UMA_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
for (i in 1:nrow(MC_postovershoot_fallback_movements)){
  # if a movement doesn't have any data, don't plot it
  filter(UMA_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> UMA_wild_data_to_plot
  filter(UMA_hatchery_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> UMA_hatchery_data_to_plot
  
  if(nrow(UMA_wild_data_to_plot)>0 | nrow(UMA_hatchery_data_to_plot)>0){
    UMA_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect_fit_to_data(origin_select = "Umatilla River",
                                                                            wild_move_prob_array = UMA_wild_spilldays_move_prob_array,
                                                                            hatchery_move_prob_array = UMA_hatchery_spilldays_move_prob_array,
                                                                            wild_covariate_experiences = UMA_wild_covariate_experiences, 
                                                                            hatchery_covariate_experiences = UMA_hatchery_covariate_experiences,
                                                                            movements_evaluated = MC_postovershoot_fallback_movements,
                                                                            spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
                                                                            from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
                                                                            plot_title = paste0("Umatilla - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
    
    ggsave(here::here("stan_actual", "output", "fit_to_data", "spilldays_v3", paste0("UMA_compare_fallback_", 
                                                                                        MC_postovershoot_fallback_movements$dam[i],
                                                                                        "_spilldays.png")),
           UMA_fallback_spilldays_plots[[i]]$plot, height = 8, width = 8)
  } else {
    # if it is zero, do nothing
    
  }
}

# Yakima - evaluate all fallback move probs
# Yakima is wild only
YAK_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "Yakima River", movements = MC_postovershoot_fallback_movements)

# Yakima - get covariate experiences
YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Yakima River")

# Yakima - use a loop to plot all of them
YAK_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
for (i in 1:nrow(MC_postovershoot_fallback_movements)){
  # if a movement doesn't have any data, don't plot it
  filter(YAK_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> YAK_data_to_plot
  
  if(nrow(YAK_data_to_plot)>0){
    YAK_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect_fit_to_data(origin_select = "Yakima River",
                                                                            wild_move_prob_array = YAK_wild_spilldays_move_prob_array,
                                                                            
                                                                            wild_covariate_experiences = YAK_wild_covariate_experiences, 
                                                                            
                                                                            movements_evaluated = MC_postovershoot_fallback_movements,
                                                                            spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
                                                                            from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
                                                                            plot_title = paste0("Yakima - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
    
    ggsave(here::here("stan_actual", "output", "fit_to_data", "spilldays_v3", paste0("YAK_compare_fallback_", 
                                                                                        MC_postovershoot_fallback_movements$dam[i],
                                                                                        "_spilldays.png")),
           YAK_fallback_spilldays_plots[[i]]$plot, height = 8, width = 8)
  } else {
    # if it is zero, do nothing
    
  }
}

# Walla Walla - evaluate all fallback move probs
WAWA_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "Walla Walla River", movements = MC_postovershoot_fallback_movements)
WAWA_hatchery_spilldays_move_prob_array <- estimate_spilldays_effect_MCH(origin_select= "Walla Walla River", movements = MC_postovershoot_fallback_movements)

# Walla Walla - get covariate experiences
WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Walla Walla River")
WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Walla Walla River")

# Walla Walla - use a loop to plot all of them
WAWA_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
for (i in 1:nrow(MC_postovershoot_fallback_movements)){
  # if a movement doesn't have any data, don't plot it
  filter(WAWA_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> WAWA_wild_data_to_plot
  filter(WAWA_hatchery_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> WAWA_hatchery_data_to_plot
  
  if(nrow(WAWA_wild_data_to_plot)>0 | nrow(WAWA_hatchery_data_to_plot)>0){
    WAWA_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect_fit_to_data(origin_select = "Walla Walla River",
                                                                             wild_move_prob_array = WAWA_wild_spilldays_move_prob_array,
                                                                             hatchery_move_prob_array = WAWA_hatchery_spilldays_move_prob_array,
                                                                             wild_covariate_experiences = WAWA_wild_covariate_experiences, 
                                                                             hatchery_covariate_experiences = WAWA_hatchery_covariate_experiences,
                                                                             movements_evaluated = MC_postovershoot_fallback_movements,
                                                                             spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
                                                                             from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
                                                                             plot_title = paste0("Walla Walla - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
    
    ggsave(here::here("stan_actual", "output", "fit_to_data", "spilldays_v3", paste0("WAWA_compare_fallback_", 
                                                                                        MC_postovershoot_fallback_movements$dam[i],
                                                                                        "_spilldays.png")),
           WAWA_fallback_spilldays_plots[[i]]$plot, height = 8, width = 8)
  } else {
    # if it is zero, do nothing
    
  }
}



#### Snake River ####
# only one post-overshoot fallback: Tucannon

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
    
    ggsave(here::here("stan_actual", "output", "fit_to_data", "spilldays_v3", paste0("TUC_compare_fallback_", 
                                                                                        SR_postovershoot_fallback_movements$dam[i],
                                                                                        "_spilldays.png")),
           TUC_fallback_spilldays_plots[[i]]$plot, height = 8, width = 8)
  } else {
    # if it is zero, do nothing
    
  }
}

# Tucannon fit is terrible
# Let's look at other origins with Priest Rapids - do we have a consistent issue across SR populations?
# They're all single digit percentages, but we can look at Asotin Creek or Imnaha as some examples

SR_columbia_postovershoot_fallback_movements <- data.frame(dam = c("PRA"),
                                                           from = c(4), to = c(3))
# Imnaha River - evaluate all fallback move probs
IMN_wild_spilldays_move_prob_array <- estimate_spilldays_effect_SRW(origin_select= "Imnaha River", movements = SR_columbia_postovershoot_fallback_movements)
IMN_hatchery_spilldays_move_prob_array <- estimate_spilldays_effect_SRH(origin_select= "Imnaha River", movements = SR_columbia_postovershoot_fallback_movements)

# Imnaha River - get covariate experiences
IMN_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Imnaha River")
IMN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Imnaha River")

# Imnaha River - use a loop to plot all of them
IMN_fallback_spilldays_plots <- vector(mode = "list", length = nrow(SR_columbia_postovershoot_fallback_movements))
for (i in 1:nrow(SR_columbia_postovershoot_fallback_movements)){
  # if a movement doesn't have any data, don't plot it
  filter(IMN_wild_covariate_experiences, state == SR_columbia_postovershoot_fallback_movements$from[i]) -> IMN_wild_data_to_plot
  filter(IMN_hatchery_covariate_experiences, state == SR_columbia_postovershoot_fallback_movements$from[i]) -> IMN_hatchery_data_to_plot
  
  if(nrow(IMN_wild_data_to_plot)>0 | nrow(IMN_hatchery_data_to_plot)>0){
    IMN_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect_fit_to_data(origin_select = "Imnaha River",
                                                                                        wild_move_prob_array = IMN_wild_spilldays_move_prob_array,
                                                                                        hatchery_move_prob_array = IMN_hatchery_spilldays_move_prob_array,
                                                                                        wild_covariate_experiences = IMN_wild_covariate_experiences, 
                                                                                        hatchery_covariate_experiences = IMN_hatchery_covariate_experiences,
                                                                                        movements_evaluated = SR_columbia_postovershoot_fallback_movements,
                                                                                        spill_predict = seq(0, max(SR_spill_days_data[,SR_columbia_postovershoot_fallback_movements$from[i]]),length = 100),
                                                                                        from = SR_columbia_postovershoot_fallback_movements$from[i], to = SR_columbia_postovershoot_fallback_movements$to[i], 
                                                                                        plot_title = paste0("Imnaha River - post-overshoot fallback at ", SR_columbia_postovershoot_fallback_movements$dam[i]))
    
    ggsave(here::here("stan_actual", "output", "fit_to_data", "spilldays_v3", paste0("IMN_compare_fallback_", 
                                                                                  SR_columbia_postovershoot_fallback_movements$dam[i],
                                                                                  "_spilldays.png")),
           IMN_fallback_spilldays_plots[[i]]$plot, height = 8, width = 8)
  } else {
    # if it is zero, do nothing
    
  }
}

# No - Imnaha doesn't look bad - it seems to be a Tucannon-specific issue. See code below for year effect issue

#### year effects plot - this might explain the really poor fit ####

# get the winterspill effects
SRH_winterspill_variables <- SRH_fit_summary$variable[grep("bwinterspill_matrix_9_8", SRH_fit_summary$variable)]
subset(SRH_fit_summary, variable %in% SRH_winterspill_variables)

# get the year effects - raw and matt trick scaling
SRH_9_8_raw_year_variables <- SRH_fit_summary$variable[grep("byearxorigin5_raw_vector_9_8", SRH_fit_summary$variable)]
subset(SRH_fit_summary, variable %in% SRH_9_8_raw_year_variables) -> SRH_9_8_raw_year_variables_summary
SRH_9_8_scale_year_variables <- SRH_fit_summary$variable[grep("sigma_yearxorigin5_vector_9_8", SRH_fit_summary$variable)]
subset(SRH_fit_summary, variable %in% SRH_9_8_scale_year_variables)

# how much of an outlier are these year effects?
## check out raw
SRH_all_year_raw_variables <- SRH_fit_summary$variable[grep("byearxorigin", SRH_fit_summary$variable)]
SRH_all_year_raw_variables <- SRH_all_year_raw_variables[grep("raw_vector", SRH_all_year_raw_variables)]
# drop any NDE parameters
SRH_all_year_raw_variables <- SRH_all_year_raw_variables[!(SRH_all_year_raw_variables %in% SRH_all_year_raw_variables[grep("_NDE", SRH_all_year_raw_variables)])]
subset(SRH_fit_summary, variable %in% SRH_all_year_raw_variables) -> SRH_all_year_raw_variables_summary
# highlight this movement & note if 90% CI contains 0
SRH_all_year_raw_variables_summary %>% 
  mutate(movement = ifelse(variable %in% SRH_9_8_raw_year_variables, "LGR fallback", 
                           ifelse(grepl("_36|_37|_38", variable), "Salmon/Grande Ronde/Clearwater", "other"))) %>% 
  mutate(significant = ifelse(sign(q5) == sign(q95), "yes", "no")) -> SRH_all_year_raw_variables_summary

SRH_year_raw_plot <- ggplot(SRH_all_year_raw_variables_summary, aes(x = variable, y = median, ymin = q5, ymax = q95, color = movement, alpha = significant)) +
  geom_point() + 
  geom_errorbar() +
  scale_alpha_manual(values = c(0.1, 1)) +
  ggtitle("Snake River Hatchery Model: Year effect, unscaled") +
  ylab("Parameter Estimate") +
  xlab("Parameter") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(here::here("stan_actual", "output", "fit_to_data", "spilldays_v3", "SRH_year_raw_plot.png"),
       SRH_year_raw_plot, 
       height = 8, width = 8)

## check out sigma (scaling)
SRH_all_year_sigma_variables <- SRH_fit_summary$variable[grep("sigma_yearxorigin", SRH_fit_summary$variable)]
SRH_all_year_sigma_variables <- SRH_all_year_sigma_variables[!(SRH_all_year_sigma_variables %in% SRH_all_year_sigma_variables[grep("\\[", SRH_all_year_sigma_variables)])]
# drop any NDE parameters
SRH_all_year_sigma_variables <- SRH_all_year_sigma_variables[!(SRH_all_year_sigma_variables %in% SRH_all_year_sigma_variables[grep("_NDE", SRH_all_year_sigma_variables)])]
subset(SRH_fit_summary, variable %in% SRH_all_year_sigma_variables) -> SRH_all_year_sigma_variables_summary
# highlight this movement & note if 90% CI contains 0
SRH_all_year_sigma_variables_summary %>% 
  mutate(movement = ifelse(variable %in% SRH_9_8_scale_year_variables, "LGR fallback", 
                           ifelse(grepl("_36|_37|_38", variable), "Salmon/Grande Ronde/Clearwater", "other")))  -> SRH_all_year_sigma_variables_summary

SRH_year_scale_plot <- ggplot(SRH_all_year_sigma_variables_summary, aes(x = variable, y = median, ymin = q5, ymax = q95, color = movement)) +
  geom_point() + 
  geom_errorbar() +
  ggtitle("Snake River Hatchery Model: Year effect, scaling values") +
  ylab("Parameter Estimate") +
  xlab("Parameter") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(here::here("stan_actual", "output", "fit_to_data", "spilldays_v3", "SRH_year_scale_plot.png"),
       SRH_year_scale_plot, 
       height = 8, width = 8)

print(subset(SRH_all_year_raw_variables_summary, significant == "yes"), n = 100)
# the other ones make a lot of sense - 9 - 36 is Clearwater, 9-37 is salmon, 9-38 is Grande Ronde
print(subset(SRH_all_year_sigma_variables_summary, q5 > 0.5), n = 100)

# One more plot: Show the raw parameter estimates vs. the winter spill values for each year
as.data.frame(SRH_envir$data$winter_spill_days_data) %>% 
  rownames_to_column("year") %>% 
  mutate(year = as.numeric(year))-> winter_spill_days_data
SRH_9_8_raw_year_variables_summary %>% 
  mutate(year = row_number()) -> SRH_9_8_raw_year_variables_summary

SRH_9_8_raw_year_variables_summary %>% 
  left_join(., winter_spill_days_data, by = "year") -> SRH_9_8_raw_year_variables_summary_spilldata

TUC_year_v_winterspill_9_8_plot <- ggplot(SRH_9_8_raw_year_variables_summary_spilldata, aes(x = LGR*100, y= median)) +
  geom_point() +
  xlab("Winter spill days at Lower Granite Dam") +
  ylab("Year effect, median parameter estimate") +
  ggtitle("Winter spill days vs. year effect for that year")

ggsave(here::here("stan_actual", "output", "fit_to_data", "spilldays_v3", "TUC_year_v_winterspill_9_8_plot.png"),
       TUC_year_v_winterspill_9_8_plot, 
       height = 8, width = 8)

### Let's look at Upper Columbia - Wenatchee hatchery is fitting fallback at WEL poorly

# get the winterspill effects
UCH_winterspill_variables_7_6 <- UCH_fit_summary$variable[grep("bwinterspill_matrix_7_6", UCH_fit_summary$variable)]
subset(UCH_fit_summary, variable %in% UCH_winterspill_variables_7_6)
UCH_winterspill_variables_7_6 <- UCH_fit_summary$variable[grep("bwinterspill_matrix_7_6", UCH_fit_summary$variable)]
subset(UCH_fit_summary, variable %in% UCH_winterspill_variables_7_6)

# get the year effects - raw and matt trick scaling
UCH_7_6_raw_year_variables <- UCH_fit_summary$variable[grep("byearxorigin1_raw_vector_7_6", UCH_fit_summary$variable)]
subset(UCH_fit_summary, variable %in% UCH_7_6_raw_year_variables) -> UCH_7_6_raw_year_variables_summary
UCH_7_6_scale_year_variables <- UCH_fit_summary$variable[grep("sigma_yearxorigin1_vector_7_6", UCH_fit_summary$variable)]
subset(UCH_fit_summary, variable %in% UCH_7_6_scale_year_variables)
UCH_7_6_raw_year_variables <- UCH_fit_summary$variable[grep("byearxorigin1_raw_vector_7_6", UCH_fit_summary$variable)]
subset(UCH_fit_summary, variable %in% UCH_7_6_raw_year_variables)
UCH_7_6_scale_year_variables <- UCH_fit_summary$variable[grep("sigma_yearxorigin1_vector_7_6", UCH_fit_summary$variable)]
subset(UCH_fit_summary, variable %in% UCH_7_6_scale_year_variables)

# how much of an outlier are these year effects?
## check out raw
UCH_all_year_raw_variables <- UCH_fit_summary$variable[grep("byearxorigin", UCH_fit_summary$variable)]
UCH_all_year_raw_variables <- UCH_all_year_raw_variables[grep("raw_vector", UCH_all_year_raw_variables)]
# drop any NDE parameters
UCH_all_year_raw_variables <- UCH_all_year_raw_variables[!(UCH_all_year_raw_variables %in% UCH_all_year_raw_variables[grep("_NDE", UCH_all_year_raw_variables)])]
subset(UCH_fit_summary, variable %in% UCH_all_year_raw_variables) -> UCH_all_year_raw_variables_summary
# highlight this movement & note if 90% CI contains 0
UCH_all_year_raw_variables_summary %>% 
  # mutate(movement = ifelse(variable %in% UCH_7_6_raw_year_variables, "LGR fallback", 
  #                          ifelse(grepl("_36|_37|_38", variable), "Salmon/Grande Ronde/Clearwater", "other"))) %>% 
  mutate(movement = ifelse(variable %in% UCH_7_6_raw_year_variables, "WEL fallback", "other")) %>% 
  mutate(significant = ifelse(sign(q5) == sign(q95), "yes", "no")) -> UCH_all_year_raw_variables_summary

UCH_year_raw_plot <- ggplot(UCH_all_year_raw_variables_summary, aes(x = variable, y = median, ymin = q5, ymax = q95, color = movement, alpha = significant)) +
  geom_point() + 
  geom_errorbar() +
  scale_alpha_manual(values = c(0.1, 1)) +
  ggtitle("Upper Columbia Hatchery Model: Year effect, unscaled") +
  ylab("Parameter Estimate") +
  xlab("Parameter") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(here::here("stan_actual", "output", "fit_to_data", "spilldays_v3", "UCH_year_raw_plot.png"),
       UCH_year_raw_plot, 
       height = 8, width = 8)

## check out sigma (scaling)
UCH_all_year_sigma_variables <- UCH_fit_summary$variable[grep("sigma_yearxorigin", UCH_fit_summary$variable)]
UCH_all_year_sigma_variables <- UCH_all_year_sigma_variables[!(UCH_all_year_sigma_variables %in% UCH_all_year_sigma_variables[grep("\\[", UCH_all_year_sigma_variables)])]
# drop any NDE parameters
UCH_all_year_sigma_variables <- UCH_all_year_sigma_variables[!(UCH_all_year_sigma_variables %in% UCH_all_year_sigma_variables[grep("_NDE", UCH_all_year_sigma_variables)])]
subset(UCH_fit_summary, variable %in% UCH_all_year_sigma_variables) -> UCH_all_year_sigma_variables_summary
# highlight this movement & note if 90% CI contains 0
UCH_all_year_sigma_variables_summary %>% 
  mutate(movement = ifelse(variable %in% UCH_7_6_scale_year_variables, "WEL fallback", "other"))  -> UCH_all_year_sigma_variables_summary

UCH_year_scale_plot <- ggplot(UCH_all_year_sigma_variables_summary, aes(x = variable, y = median, ymin = q5, ymax = q95, color = movement)) +
  geom_point() + 
  geom_errorbar() +
  ggtitle("Upper Columbia Hatchery Model: Year effect, scaling values") +
  ylab("Parameter Estimate") +
  xlab("Parameter") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(here::here("stan_actual", "output", "fit_to_data", "spilldays_v3", "UCH_year_scale_plot.png"),
       UCH_year_scale_plot, 
       height = 8, width = 8)

print(arrange(subset(UCH_all_year_raw_variables_summary, significant == "yes"), desc(median)), n = 100)
# the other ones make a lot of sense - 9 - 36 is Clearwater, 9-37 is salmon, 9-38 is Grande Ronde
print(arrange(subset(UCH_all_year_sigma_variables_summary, q5 > 0.5), desc(median)), n = 100)

# One more plot: Show the raw parameter estimates vs. the winter spill values for each year
as.data.frame(UCH_envir$data$winter_spill_days_data) %>% 
  rownames_to_column("year") %>% 
  mutate(year = as.numeric(year))-> winter_spill_days_data
UCH_7_6_raw_year_variables_summary %>% 
  mutate(year = row_number()) -> UCH_7_6_raw_year_variables_summary

UCH_7_6_raw_year_variables_summary %>% 
  left_join(., winter_spill_days_data, by = "year") -> UCH_7_6_raw_year_variables_summary_spilldata

WEN_year_v_winterspill_7_6_plot <- ggplot(UCH_7_6_raw_year_variables_summary_spilldata, aes(x = WEL*100, y= median)) +
  geom_point() +
  xlab("Winter spill days at Wells Dam") +
  ylab("Year effect, median parameter estimate") +
  ggtitle("Winter spill days vs. year effect for that year")

ggsave(here::here("stan_actual", "output", "fit_to_data", "spilldays_v3", "WEN_year_v_winterspill_7_6_plot.png"),
       WEN_year_v_winterspill_7_6_plot, 
       height = 8, width = 8)


#### Multiple movements on same plot ####

# John Day River (W), Umatilla River (H, W), Yakima River (W), Walla Walla River (H, W),
# Wenatchee River (H, W), Entiat River (W), Tucannon River (H, W), Imnaha River (H, W)

# For some origins, there are multiple interesting and deleterious movements - 
# for example, Middle Columbia fish could go to loss, Deschutes, or overshoot

# plotting function
plot_compare_rear_spill_effect_multiple_movements <- function(origin_select,
                                                              wild_move_prob_array = NULL, hatchery_move_prob_array = NULL,
                                                              wild_covariate_experiences = NULL, hatchery_covariate_experiences = NULL,
                                                              movements_evaluated, spill_predict,
                                                              plot_title = NULL, plot_legend = FALSE){
  
  # create df to index to right dam spill joining; here we need to change names of 
  # WEL and LGR because we have slow v fast there
  dam_index <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL_quick", "ICH", "LGR_quick"),
                          state = seq(2,9))
  
  niter <- 4000 # this is the number of draws we have
  
  # Initialize a df to store every movement, then loop through each movement and bind
  # rear_spill_move_prob_quantiles <- data.frame(spill_actual = NA, `0.025` = NA, `0.5` = NA, `0.975` = NA,
  #                                             rear = NA, from = NA, to = NA)
  rear_spill_move_prob_quantiles <- data.frame()
  for (i in 1:nrow(movements_evaluated)){
    # First, determine if this origin has both a hatchery and a wild population
    
    # If hatchery is NA, run wild only
    if (is.null(hatchery_covariate_experiences)){
      wild_spill_move_prob <- as.data.frame(wild_move_prob_array[,,i])
      
      colnames(wild_spill_move_prob) <- paste0("iter", 1:niter) 
      wild_spill_move_prob$spill <- spill_predict
      
      # Add a column with the actual spilldays
      if (movements_evaluated$from[i] %in% c(1:9)){
        wild_spill_move_prob %>% 
          mutate(spill_actual = spill*100)  -> wild_spill_move_prob
      } else {
        wild_spill_move_prob %>% 
          mutate(spill_actual = 1) -> wild_spill_move_prob
      }
      
      
      
      
      wild_spill_move_prob %>% 
        pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
        group_by(spill_actual) %>% 
        dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
        pivot_wider(names_from = q, values_from = prob) %>% 
        mutate(rear = "wild") %>% 
        mutate(from = movements_evaluated$from[i],
               to = movements_evaluated$to[i]) -> wild_spill_move_prob_quantiles
      
      # get data organized for rug plot
      wild_covariate_experiences %>% 
        mutate(rear = "wild") -> wild_covariate_experiences
      
      wild_covariate_experiences %>% 
        # keep only the state that is our from state
        filter(state == movements_evaluated$from[i]) -> covariate_experiences
      
      rear_spill_move_prob_quantiles %>% 
        bind_rows(., wild_spill_move_prob_quantiles) -> rear_spill_move_prob_quantiles
      
    } 
    # If wild is NA, run hatchery only
    else if (is.null(wild_covariate_experiences)){
      hatchery_spill_move_prob <- as.data.frame(hatchery_move_prob_array[,,i])
      
      colnames(hatchery_spill_move_prob) <- paste0("iter", 1:niter)
      hatchery_spill_move_prob$spill <- spill_predict
      
      # Add a column with the actual spilldays
      if (movements_evaluated$from[i] %in% c(1:9)){
        hatchery_spill_move_prob %>% 
          mutate(spill_actual = spill*100) -> hatchery_spill_move_prob
      } else {
        hatchery_spill_move_prob %>% 
          mutate(spill_actual = 1) -> hatchery_spill_move_prob
      }
      
      
      
      
      hatchery_spill_move_prob %>% 
        pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
        group_by(spill_actual) %>% 
        summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
        pivot_wider(names_from = q, values_from = prob) %>% 
        mutate(rear = "hatchery") %>% 
        mutate(from = movements_evaluated$from[i],
               to = movements_evaluated$to[i]) -> hatchery_spill_move_prob_quantiles
      
      # get data organized for rug plot
      hatchery_covariate_experiences %>% 
        mutate(rear = "hatchery") -> hatchery_covariate_experiences
      
      hatchery_covariate_experiences %>% 
        # keep only the state that is our from state
        filter(state == movements_evaluated$from[i]) -> covariate_experiences
      
      # combine wild and hatchery
      rear_spill_move_prob_quantiles %>% 
        bind_rows(., hatchery_spill_move_prob_quantiles) -> rear_spill_move_prob_quantiles
      
    } 
    # else run both
    else {
      wild_spill_move_prob <- as.data.frame(wild_move_prob_array[,,i])
      
      colnames(wild_spill_move_prob) <- paste0("iter", 1:niter) 
      wild_spill_move_prob$spill <- spill_predict
      
      # Add a column with the actual spilldays
      if (movements_evaluated$from[i] %in% c(1:9)){
        wild_spill_move_prob %>% 
          mutate(spill_actual = spill*100) -> wild_spill_move_prob
      } else {
        wild_spill_move_prob %>% 
          mutate(spill_actual = 1) -> wild_spill_move_prob
      }
      
      
      wild_spill_move_prob %>% 
        pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
        group_by(spill_actual) %>% 
        summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
        pivot_wider(names_from = q, values_from = prob) %>% 
        mutate(rear = "wild") %>% 
        mutate(from = movements_evaluated$from[i],
               to = movements_evaluated$to[i]) -> wild_spill_move_prob_quantiles
      
      # get data organized for rug plot
      wild_covariate_experiences %>% 
        mutate(rear = "wild") -> wild_covariate_experiences
      
      wild_covariate_experiences %>% 
        # keep only the state that is our from state
        filter(state == movements_evaluated$from[i]) -> wild_covariate_experiences
      
      hatchery_spill_move_prob <- as.data.frame(hatchery_move_prob_array[,,i])
      
      colnames(hatchery_spill_move_prob) <- paste0("iter", 1:niter) 
      hatchery_spill_move_prob$spill <- spill_predict
      
      # Add a column with the actual spilldays
      if (movements_evaluated$from[i] %in% c(1:9)){
        hatchery_spill_move_prob %>% 
          mutate(spill_actual = spill*100) -> hatchery_spill_move_prob
      } else {
        hatchery_spill_move_prob %>% 
          mutate(spill_actual = 1) -> hatchery_spill_move_prob
      }
      
      hatchery_spill_move_prob %>% 
        pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
        group_by(spill_actual) %>% 
        summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
        pivot_wider(names_from = q, values_from = prob) %>% 
        mutate(rear = "hatchery") %>% 
        mutate(from = movements_evaluated$from[i],
               to = movements_evaluated$to[i]) -> hatchery_spill_move_prob_quantiles
      
      # get data organized for rug plot
      hatchery_covariate_experiences %>% 
        mutate(rear = "hatchery") -> hatchery_covariate_experiences
      
      hatchery_covariate_experiences %>% 
        # keep only the state that is our from state
        filter(state == movements_evaluated$from[i]) -> hatchery_covariate_experiences
      
      # combine wild and hatchery
      rear_spill_move_prob_quantiles %>% 
        bind_rows(., wild_spill_move_prob_quantiles) %>% 
        bind_rows(., hatchery_spill_move_prob_quantiles) -> rear_spill_move_prob_quantiles
      
      wild_covariate_experiences %>% 
        bind_rows(., hatchery_covariate_experiences) -> covariate_experiences
      
    }
    
  }
  
  
  
  
  # convert spill to spill_actual in covariate experiences
  if (movements_evaluated$from[i] %in% c(1:9)){
    covariate_experiences %>% 
      mutate(spill_actual = winter_spill*100) -> covariate_experiences
  } else {
    covariate_experiences %>% 
      mutate(spill_actual = 1) -> covariate_experiences
  }
  
  
  # Remove any instances for rug plot of fish that would have received spill0 covariate
  covariate_experiences %>% 
    dplyr::rename(date_numeric = date) %>% 
    # keep only jan/feb/mar 
    mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
    mutate(month = month(date)) %>% 
    filter(!(month %in% c(1,2,3,4,5))) -> covariate_experiences
  
  # rear_spill_move_prob_quantiles$to <- as.character(rear_spill_move_prob_quantiles$to)
  
  
  movement_colors <- c("Fallback" = "#6a3d9a", "Overshoot - PRA" = "#ff7f00",
                       "Overshoot - RIS" = "#ff7f00", "Overshoot - LGR" = "#ff7f00",
                       "Overshoot - RRE" = "#ff7f00", "Overshoot - WEL" = "#ff7f00",
                       "Overshoot - ICH" = "#ff7f00", "Overshoot" = "#ff7f00",
                       "Additional Overshoot" = "#ff7f00",
                       "Loss" = "#e31a1c")
  
  rear_spill_move_prob_quantiles %>% 
    left_join(., movements_evaluated, by = "to") -> rear_spill_move_prob_quantiles
  
  # if (plot_legend == TRUE){
  #   rear_spill_move_prob_plot <- ggplot(rear_spill_move_prob_quantiles, aes(x = spill_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`,
  #                                                                           color = Movement, fill = Movement)) +
  #     geom_line() +
  #     geom_ribbon(alpha = 0.2, color = NA) +
  #     geom_rug(data = covariate_experiences, aes(x = spill_actual), inherit.aes = FALSE,
  #              sides = "t", length = unit(0.3, "cm"), outside = TRUE) +
  #     scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
  #     scale_x_continuous(lim = c(0,NA), expand = c(0,0)) +
  #     scale_color_manual(values = movement_colors) +
  #     scale_fill_manual(values =  movement_colors) +
  #     xlab("Winter spill days") +
  #     ylab("Movement probability") +
  #     coord_cartesian(clip = "off")
  # } else {
  #   # suppress common legend - for combined plot
  #   rear_spill_move_prob_plot <- ggplot(rear_spill_move_prob_quantiles, aes(x = spill_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
  #                                                                           color = Movement, fill = Movement)) +
  #     geom_line() +
  #     geom_ribbon(alpha = 0.2, color = NA) +
  #     geom_rug(data = covariate_experiences, aes(x = spill_actual), inherit.aes = FALSE,
  #              sides = "t", length = unit(0.3, "cm"), outside = TRUE) +
  #     scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
  #     scale_x_continuous(lim = c(0,NA), expand = c(0,0)) +
  #     scale_color_manual(values = movement_colors) +
  #     scale_fill_manual(values =  movement_colors) +
  #     xlab("Winter spill days") +
  #     ylab("Movement probability") +
  #     coord_cartesian(clip = "off") +
  #     theme(legend.position = "none")
  # }
  
  if (plot_legend == TRUE){
    
    combined_plot <- ggplot(rear_spill_move_prob_quantiles, aes(x = spill_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                               color = Movement, fill = Movement)) +
      geom_line(linewidth = 2.5) +
      geom_ribbon(alpha = 0.2, color = NA) +
      # geom_rug(data = covariate_experiences, aes(x = spill_actual), inherit.aes = FALSE,
      #          sides = "t", length = unit(0.3, "cm"), outside = TRUE) +
      scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
      scale_x_continuous(lim = c(0,ceiling(max(covariate_experiences$spill_actual))), expand = c(0,0)) +
      scale_color_manual(values = movement_colors) +
      scale_fill_manual(values =  movement_colors) +
          xlab("Winter spill days") +
      ylab("Movement probability") +
      coord_cartesian(clip = "off") +
      theme(panel.grid.major = element_line(color = "gray90"),
            panel.background = element_rect(fill = "white", color = NA),
            panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
            legend.key.height = unit(1.25, "cm"),
            legend.key.width = unit(1.25, "cm"),
            legend.title = element_text(size = 25),
            legend.text = element_text(size = 15),
            # these plot margins are to leave space for the population name on the big figure
            plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))
    # for testing
    spill_legend <- ggpubr::get_legend(combined_plot)

    spill_plot_legend_gg <- as_ggplot(spill_legend)
    
    # for testing
    # ggsave(here::here("stan_actual", "output", "paper_figures", "01_legend_test.png"), spill_plot_legend_gg, height = 4, width = 4)
    
  } else {
    # suppress common legend - for combined plot
    rear_spill_move_prob_plot <- ggplot(rear_spill_move_prob_quantiles, aes(x = spill_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                          color = Movement, fill = Movement)) +
      geom_line() +
      geom_ribbon(alpha = 0.2, color = NA) +
      # geom_rug(data = covariate_experiences, aes(x = spill_actual), inherit.aes = FALSE,
      #          sides = "t", length = unit(0.3, "cm"), outside = TRUE) +
      scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
      # common x-axis scale across all populations
      coord_cartesian(xlim = c(0, 80), expand = FALSE) +
      # scale_x_continuous(lim = c(0,ceiling(max(covariate_experiences$spill_actual))), expand = c(0,0)) +
      scale_color_manual(values = movement_colors) +
      scale_fill_manual(values =  movement_colors) +
          xlab("Winter spill days") +
      ylab("Movement probability") +
      theme(legend.position = "none",
            panel.grid.major = element_line(color = "gray90"),
            panel.background = element_rect(fill = "white", color = "black"),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=0.4),
            # turn off the axis titles on each individual plot and just show one for whole plot
            axis.title = element_blank(),
            # these plot margins are to leave space for the population name on the big figure
            plot.margin = unit(c(0, 0.2, 0.2, 0.2),"cm"))
    
    density_plot <- ggplot(data = covariate_experiences, aes(spill_actual))+
      # geom_density(alpha = 0.1) +
      # geom_rug(sides = "t", length = unit(0.2, "cm"), outside = FALSE) +
      # geom_histogram(alpha = 0.5, bins = 60) +
      geom_histogram(aes(y=..count../sum(..count..)), alpha = 0.5, bins = 60) +
      ylab("Density") +
      # scale_x_continuous(lim = c(floor(min(covariate_experiences$temp_actual)),ceiling(max(covariate_experiences$temp_actual))), expand = c(0,0)) +
      # scale_x_continuous(lim = c(0, 23), expand = c(0,0)) +
      # scale_y_continuous(lim = c(0,0.75), expand = c(0,0),
      #                    breaks = c(0, 0.25, 0.50)) +
      # these labels aren't real, they just need to be the same size as the labels from the temp effect plot
      scale_y_continuous(n.breaks = 2, labels = c("0.00", "1.00")) +
      xlim(0,80)+
      # coord_cartesian(xlim = c(floor(min(covariate_experiences$temp_actual)),ceiling(max(covariate_experiences$temp_actual)))) +
      # common x-axis scale across all populations
      coord_cartesian(xlim = c(0, 80), expand = FALSE) +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_text(color = "white"),
            axis.ticks.x = element_blank(),
            # axis.ticks.length.x=unit(.1, "cm"),
            axis.ticks.length.x=unit(0, "cm"),
            axis.ticks.y = element_line(color = "white"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(color = "white"),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            # these plot margins are to leave space for the population name on the big figure
            plot.margin = unit(c(1.0, 0.2, 0, 0.2),"cm"))
    # testing theme
    # theme(
    #       # these plot margins are to leave space for the population name on the big figure
    #       plot.margin = unit(c(0.2, 0.2, 0, 0.2),"cm"))
    
    combined_plot <- ggarrange(density_plot, rear_spill_move_prob_plot, nrow = 2, ncol = 1,
                               heights = c(2,6))
    
    # for testing
    # ggsave(here::here("stan_actual", "output", "paper_figures", "01_test.png"), combined_plot, height = 6, width = 6)
    
    
  }
  
  
  
  
  
  
  return(combined_plot)
}

#### Prepare data to create figures for individual origins: ####

### John Day River ###
JDR_spilldays_movements <- data.frame(from = c(3, 3, 3, 3), to = c(2, 4, 8, 43),
                                      Movement = c("Fallback", "Overshoot - ICH", "Overshoot - PRA", 
                                                   "Loss"))
JDR_wild_spill_move_prob_array <- estimate_spilldays_effect_MCW(origin_select = "John Day River", movements = JDR_spilldays_movements)
JDR_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "John Day River")

### Umatilla River ###
UMA_spilldays_movements <- data.frame(from = c(3, 3, 3, 3), to = c(2, 4, 8, 43),
                                      Movement = c("Fallback", "Overshoot - ICH", "Overshoot - PRA", 
                                                   "Loss"))
UMA_wild_spill_move_prob_array <- estimate_spilldays_effect_MCW(origin_select = "Umatilla River", movements = UMA_spilldays_movements)
UMA_hatchery_spill_move_prob_array <- estimate_spilldays_effect_MCH(origin_select = "Umatilla River", movements = UMA_spilldays_movements)
UMA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Umatilla River")
UMA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Umatilla River")

### Yakima River ###
YAK_spilldays_movements <- data.frame(from = c(4,4,4), to = c(3,5,43),
                                      Movement = c("Fallback",
                                                   "Additional Overshoot", "Loss"))

YAK_wild_spill_move_prob_array <- estimate_spilldays_effect_MCW(origin_select = "Yakima River", movements = YAK_spilldays_movements)
YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Yakima River")

### Walla Walla River ###
WAWA_spilldays_movements <- data.frame(from = c(8,8,8), to = c(3, 9, 43),
                                       Movement = c("Fallback", "Overshoot - LGR", 
                                                    "Loss"))
WAWA_wild_spill_move_prob_array <- estimate_spilldays_effect_MCW(origin_select = "Walla Walla River", movements = WAWA_spilldays_movements)
WAWA_hatchery_spill_move_prob_array <- estimate_spilldays_effect_MCH(origin_select = "Walla Walla River", movements = WAWA_spilldays_movements)
WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Walla Walla River")
WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Walla Walla River")

### Wenatchee River ###
WEN_spilldays_movements <- data.frame(from = c(6,6,6), to = c(5, 7, 43),
                                      Movement = c("Fallback", "Overshoot - WEL", 
                                                   "Loss"))
WEN_wild_spill_move_prob_array <- estimate_spilldays_effect_UCW(origin_select = "Wenatchee River", movements = WEN_spilldays_movements)
WEN_hatchery_spill_move_prob_array <- estimate_spilldays_effect_UCH(origin_select = "Wenatchee River", movements = WEN_spilldays_movements)
WEN_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Wenatchee River")
WEN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Wenatchee River")

### Entiat River ###
ENT_spilldays_movements <- data.frame(from = c(7,7), to = c(6, 43),
                                      Movement = c("Fallback",
                                                   "Loss"))
ENT_wild_spill_move_prob_array <- estimate_spilldays_effect_UCW(origin_select = "Entiat River", movements = ENT_spilldays_movements)
ENT_hatchery_spill_move_prob_array <- estimate_spilldays_effect_UCH(origin_select = "Entiat River", movements = ENT_spilldays_movements)
ENT_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Entiat River")
ENT_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Entiat River")

### Tucannon River ###
TUC_spilldays_movements <- data.frame(from = c(9,9), to = c(8, 43),
                                      Movement = c("Fallback",
                                                   "Loss"))
TUC_wild_spill_move_prob_array <- estimate_spilldays_effect_SRW(origin_select = "Tucannon River", movements = TUC_spilldays_movements)
TUC_hatchery_spill_move_prob_array <- estimate_spilldays_effect_SRH(origin_select = "Tucannon River", movements = TUC_spilldays_movements)
TUC_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Tucannon River")
TUC_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Tucannon River")





#### Create figures for individual origins: ####




### John Day River ###

JDR_wild_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "John Day River",
                                                                                     wild_move_prob_array = JDR_wild_spill_move_prob_array,
                                                                                     wild_covariate_experiences = JDR_wild_covariate_experiences,
                                                                                     spill_predict = seq(0, 0.90, length = 100), 
                                                                                     movements_evaluated = JDR_spilldays_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays", "JDR_wild_compare_movement_spill.png"), JDR_wild_compare_movement_spill, height = 8, width = 8)

### Umatilla River ###


# UMA wild plot
UMA_wild_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Umatilla River",
                                                                                     wild_move_prob_array = UMA_wild_spill_move_prob_array,
                                                                                     wild_covariate_experiences = UMA_wild_covariate_experiences,
                                                                                     spill_predict = seq(0, 0.90, length = 100), 
                                                                                     movements_evaluated = UMA_spilldays_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays", "UMA_wild_compare_movement_spill.png"), UMA_wild_compare_movement_spill, height = 8, width = 8)

# UMA hatchery plot
UMA_hatchery_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Umatilla River",
                                                                                         hatchery_move_prob_array = UMA_hatchery_spill_move_prob_array,
                                                                                         hatchery_covariate_experiences = UMA_hatchery_covariate_experiences,
                                                                                         spill_predict = seq(0, 0.90, length = 100), 
                                                                                         movements_evaluated = UMA_spilldays_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays", "UMA_hatchery_compare_movement_spill.png"), UMA_hatchery_compare_movement_spill, height = 8, width = 8)

### Yakima River ###

YAK_wild_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Yakima River",
                                                                                     wild_move_prob_array = YAK_wild_spill_move_prob_array,
                                                                                     
                                                                                     wild_covariate_experiences = YAK_wild_covariate_experiences,
                                                                                     spill_predict = seq(0, 0.90, length = 100), 
                                                                                     movements_evaluated = YAK_spilldays_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays", "YAK_wild_compare_movement_spill.png"), YAK_wild_compare_movement_spill, height = 8, width = 8)



### Walla Walla River ###


# WAWA wild plot
WAWA_wild_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Walla Walla River",
                                                                                      wild_move_prob_array = WAWA_wild_spill_move_prob_array,
                                                                                      wild_covariate_experiences = WAWA_wild_covariate_experiences,
                                                                                      spill_predict = seq(0, 0.90, length = 100), 
                                                                                      movements_evaluated = WAWA_spilldays_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays", "WAWA_wild_compare_movement_spill.png"), WAWA_wild_compare_movement_spill, height = 8, width = 8)

# WAWA hatchery plot
WAWA_hatchery_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Walla Walla River",
                                                                                          hatchery_move_prob_array = WAWA_hatchery_spill_move_prob_array,
                                                                                          hatchery_covariate_experiences = WAWA_hatchery_covariate_experiences,
                                                                                          spill_predict = seq(0, 0.90, length = 100), 
                                                                                          movements_evaluated = WAWA_spilldays_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays", "WAWA_hatchery_compare_movement_spill.png"), WAWA_hatchery_compare_movement_spill, height = 8, width = 8)

### Wenatchee River ###


# WEN wild plot
WEN_wild_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Wenatchee River",
                                                                                     wild_move_prob_array = WEN_wild_spill_move_prob_array,
                                                                                     wild_covariate_experiences = WEN_wild_covariate_experiences,
                                                                                     spill_predict = seq(0, 0.90, length = 100), 
                                                                                     movements_evaluated = WEN_spilldays_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays", "WEN_wild_compare_movement_spill.png"), WEN_wild_compare_movement_spill, height = 8, width = 8)

# WEN hatchery plot
WEN_hatchery_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Wenatchee River",
                                                                                         hatchery_move_prob_array = WEN_hatchery_spill_move_prob_array,
                                                                                         hatchery_covariate_experiences = WEN_hatchery_covariate_experiences,
                                                                                         spill_predict = seq(0, 0.90, length = 100), 
                                                                                         movements_evaluated = WEN_spilldays_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays", "WEN_hatchery_compare_movement_spill.png"), WEN_hatchery_compare_movement_spill, height = 8, width = 8)

### Entiat River ###


# ENT wild plot
ENT_wild_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Entiat River",
                                                                                     wild_move_prob_array = ENT_wild_spill_move_prob_array,
                                                                                     wild_covariate_experiences = ENT_wild_covariate_experiences,
                                                                                     spill_predict = seq(0, 0.90, length = 100), 
                                                                                     movements_evaluated = ENT_spilldays_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays", "ENT_wild_compare_movement_spill.png"), ENT_wild_compare_movement_spill, height = 8, width = 8)


### Tucannon River ###


# TUC wild plot
TUC_wild_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Tucannon River",
                                                                                     wild_move_prob_array = TUC_wild_spill_move_prob_array,
                                                                                     wild_covariate_experiences = TUC_wild_covariate_experiences,
                                                                                     spill_predict = seq(0, 0.90, length = 100), 
                                                                                     movements_evaluated = TUC_spilldays_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays", "TUC_wild_compare_movement_spill.png"), TUC_wild_compare_movement_spill, height = 8, width = 8)

# TUC hatchery plot
TUC_hatchery_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Tucannon River",
                                                                                         hatchery_move_prob_array = TUC_hatchery_spill_move_prob_array,
                                                                                         hatchery_covariate_experiences = TUC_hatchery_covariate_experiences,
                                                                                         spill_predict = seq(0, 0.90, length = 100), 
                                                                                         movements_evaluated = TUC_spilldays_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays", "TUC_hatchery_compare_movement_spill.png"), TUC_hatchery_compare_movement_spill, height = 8, width = 8)

# ### Imnaha River ###
# IMN_spilldays_movements <- data.frame(from = c(9,9), to = c(8, 43),
#                             Movement = c("Fallback",
#                                          "Loss"))
# IMN_wild_spill_move_prob_array <- estimate_spilldays_effect_SRW(origin_select = "Imnaha River", movements = IMN_spilldays_movements)
# IMN_hatchery_spill_move_prob_array <- estimate_spilldays_effect_SRH(origin_select = "Imnaha River", movements = IMN_spilldays_movements)
# IMN_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Imnaha River")
# IMN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Imnaha River")
# 
# 
# # IMN wild plot
# IMN_wild_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Imnaha River",
#                                                                                      wild_move_prob_array = IMN_wild_spill_move_prob_array,
#                                                                                      wild_covariate_experiences = IMN_wild_covariate_experiences,
#                                                                                      spill_predict = seq(0, 0.90, length = 100), 
#                                                                                      movements_evaluated = IMN_spilldays_movements)
# 
# ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays", "IMN_wild_compare_movement_spill.png"), IMN_wild_compare_movement_spill, height = 8, width = 8)
# 
# # IMN hatchery plot
# IMN_hatchery_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Imnaha River",
#                                                                                          hatchery_move_prob_array = IMN_hatchery_spill_move_prob_array,
#                                                                                          hatchery_covariate_experiences = IMN_hatchery_covariate_experiences,
#                                                                                          spill_predict = seq(0, 0.90, length = 100), 
#                                                                                          movements_evaluated = IMN_spilldays_movements)
# 
# ggsave(here::here("stan_actual", "output", "covariate_effects", "spilldays", "IMN_hatchery_compare_movement_spill.png"), IMN_hatchery_compare_movement_spill, height = 8, width = 8)


# Create the legend figure by itself
spill_plot_for_legend <-  plot_compare_rear_spill_effect_multiple_movements(origin_select = "Yakima River",
                                                                      wild_move_prob_array = YAK_wild_spill_move_prob_array,
                                                                      
                                                                      wild_covariate_experiences = YAK_wild_covariate_experiences,
                                                                      spill_predict = seq(0, 0.90, length = 100), 
                                                                      movements_evaluated = YAK_spilldays_movements,
                                                                      plot_legend= TRUE)
spill_legend <- ggpubr::get_legend(spill_plot_for_legend)
spill_plot_legend_gg <- as_ggplot(spill_legend) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), 
                                                        panel.background = element_rect(fill = "white", color = "white"))


#### Generate the figure for the paper ####

# combined_movement_spill_plot <- ggarrange(JDR_wild_compare_movement_spill, UMA_wild_compare_movement_spill, UMA_hatchery_compare_movement_spill,
#                                           YAK_wild_compare_movement_spill, WAWA_wild_compare_movement_spill, WAWA_hatchery_compare_movement_spill,
#                                           WEN_wild_compare_movement_spill, WEN_hatchery_compare_movement_spill, ENT_wild_compare_movement_spill,
#                                           TUC_wild_compare_movement_spill, TUC_hatchery_compare_movement_spill, IMN_wild_compare_movement_spill, 
#                                           IMN_hatchery_compare_movement_spill, plot_legend_gg, nrow = 4, ncol = 4,
#                                           labels = c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)", "(G)", "(H)", "(I)",
#                                                      "(J)", "(K)", "(L)", "(M)"),
#                                           label.x = 0.00, label.y = 0.95, font.label = list(size = 10, face = "plain"),
#                                           hjust = 0, vjust = 0)
# 
# ggsave(here::here("stan_actual", "output", "paper_figures", "combined_movement_spill_plot_v1.png"), combined_movement_spill_plot, height = 12, width = 16)

# drop imnaha
combined_movement_spill_plot <- ggarrange(JDR_wild_compare_movement_spill, UMA_wild_compare_movement_spill, UMA_hatchery_compare_movement_spill,
                                          YAK_wild_compare_movement_spill, WAWA_wild_compare_movement_spill, WAWA_hatchery_compare_movement_spill,
                                          WEN_wild_compare_movement_spill, WEN_hatchery_compare_movement_spill, ENT_wild_compare_movement_spill,
                                          TUC_wild_compare_movement_spill, TUC_hatchery_compare_movement_spill, 
                                          spill_plot_legend_gg, nrow = 3, ncol = 4,
                                          labels = c("(A) JDR, Natural", "(B) UMA, Natural", 
                                                     "(C) UMA, Hatchery", "(D) YAK, Natural", 
                                                     "(E) WAWA, Natural", "(F) WAWA, Hatchery", 
                                                     "(G) WEN, Natural", "(H) WEN, Hatchery", 
                                                     "(I) ENT, Natural",
                                                     "(J) TUC, Natural", "(K) TUC, Hatchery"),
                                          label.x = 0.05, label.y = 0.925, font.label = list(size = 14, face = "plain"),
                                          hjust = 0, vjust = 0)

# combined_movement_spill_plot <- annotate_figure(combined_movement_spill_plot,
#                                                bottom = textGrob("Days of winter spill", gp = gpar(cex = 1.3)),
#                                                left = textGrob("Movement probability", rot = 90, gp = gpar(cex = 1.3))) + bgcolor("white")

# Let's try this again, this time using a cowplot solution since ggpubr is
# struggling with the background color
combined_movement_spill_plot <- cowplot::ggdraw(annotate_figure(combined_movement_spill_plot,
                                                bottom = textGrob("Days of winter spill", gp = gpar(cex = 1.3)),
                                                left = textGrob("Movement probability", rot = 90, gp = gpar(cex = 1.3)))) +
  theme(plot.background = element_rect(fill="white", color = NA))

ggsave(here::here("stan_actual", "output", "paper_figures", "combined_movement_spill_plot_v2.png"), combined_movement_spill_plot, height = 12, width = 16)
