# 07-06_Fig5_spilldays_effects_plot

# This script takes the output from the stan model runs and
# plots the effects of days of winter spill

# First, need to load in all of the model runs and all of the packages.
source("R/07_model_analysis/alternative_spill_treatment/07-01_load_stan_models_marchspill.R")

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
    possible_movements <- c(possible_movements, 41)
    
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
    possible_movements <- c(possible_movements, 41)
    
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
    possible_movements <- c(possible_movements, 41)
    
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
    possible_movements <- c(possible_movements, 41)
    
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
    # ggsave(here::here("figures", "alternative_spill_treatment", "paper_figures",  "01_legend_test.png"), spill_plot_legend_gg, height = 4, width = 4)
    
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
      coord_cartesian(xlim = c(0, 32), expand = FALSE) +
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
      xlim(0,32)+
      # coord_cartesian(xlim = c(floor(min(covariate_experiences$temp_actual)),ceiling(max(covariate_experiences$temp_actual)))) +
      # common x-axis scale across all populations
      coord_cartesian(xlim = c(0, 32), expand = FALSE) +
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
    # ggsave(here::here("figures", "alternative_spill_treatment", "paper_figures",  "01_test.png"), combined_plot, height = 6, width = 6)
    
    
  }
  
  
  
  
  
  
  return(combined_plot)
}

#### Prepare data to create figures for individual origins: ####

### John Day River ###
JDR_spilldays_movements <- data.frame(from = c(3, 3, 3, 3), to = c(2, 4, 8, 41),
                                      Movement = c("Fallback", "Overshoot - ICH", "Overshoot - PRA", 
                                                   "Loss"))
JDR_wild_spill_move_prob_array <- estimate_spilldays_effect_MCW(origin_select = "John Day River", movements = JDR_spilldays_movements)
JDR_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "John Day River")

### Umatilla River ###
UMA_spilldays_movements <- data.frame(from = c(3, 3, 3, 3), to = c(2, 4, 8, 41),
                                      Movement = c("Fallback", "Overshoot - ICH", "Overshoot - PRA", 
                                                   "Loss"))
UMA_wild_spill_move_prob_array <- estimate_spilldays_effect_MCW(origin_select = "Umatilla River", movements = UMA_spilldays_movements)
UMA_hatchery_spill_move_prob_array <- estimate_spilldays_effect_MCH(origin_select = "Umatilla River", movements = UMA_spilldays_movements)
UMA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Umatilla River")
UMA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Umatilla River")

### Yakima River ###
YAK_spilldays_movements <- data.frame(from = c(4,4,4), to = c(3,5,41),
                                      Movement = c("Fallback",
                                                   "Additional Overshoot", "Loss"))

YAK_wild_spill_move_prob_array <- estimate_spilldays_effect_MCW(origin_select = "Yakima River", movements = YAK_spilldays_movements)
YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Yakima River")

### Walla Walla River ###
WAWA_spilldays_movements <- data.frame(from = c(8,8,8), to = c(3, 9, 41),
                                       Movement = c("Fallback", "Overshoot - LGR", 
                                                    "Loss"))
WAWA_wild_spill_move_prob_array <- estimate_spilldays_effect_MCW(origin_select = "Walla Walla River", movements = WAWA_spilldays_movements)
WAWA_hatchery_spill_move_prob_array <- estimate_spilldays_effect_MCH(origin_select = "Walla Walla River", movements = WAWA_spilldays_movements)
WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Walla Walla River")
WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Walla Walla River")

### Wenatchee River ###
WEN_spilldays_movements <- data.frame(from = c(6,6,6), to = c(5, 7, 41),
                                      Movement = c("Fallback", "Overshoot - WEL", 
                                                   "Loss"))
WEN_wild_spill_move_prob_array <- estimate_spilldays_effect_UCW(origin_select = "Wenatchee River", movements = WEN_spilldays_movements)
WEN_hatchery_spill_move_prob_array <- estimate_spilldays_effect_UCH(origin_select = "Wenatchee River", movements = WEN_spilldays_movements)
WEN_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Wenatchee River")
WEN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Wenatchee River")

### Entiat River ###
ENT_spilldays_movements <- data.frame(from = c(7,7), to = c(6, 41),
                                      Movement = c("Fallback",
                                                   "Loss"))
ENT_wild_spill_move_prob_array <- estimate_spilldays_effect_UCW(origin_select = "Entiat River", movements = ENT_spilldays_movements)
ENT_hatchery_spill_move_prob_array <- estimate_spilldays_effect_UCH(origin_select = "Entiat River", movements = ENT_spilldays_movements)
ENT_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Entiat River")
ENT_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Entiat River")

### Tucannon River ###
TUC_spilldays_movements <- data.frame(from = c(9,9), to = c(8, 41),
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

ggsave(here::here("figures", "alternative_spill_treatment", "spilldays", "JDR_wild_compare_movement_spill.png"), JDR_wild_compare_movement_spill, height = 8, width = 8)

### Umatilla River ###


# UMA wild plot
UMA_wild_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Umatilla River",
                                                                                     wild_move_prob_array = UMA_wild_spill_move_prob_array,
                                                                                     wild_covariate_experiences = UMA_wild_covariate_experiences,
                                                                                     spill_predict = seq(0, 0.90, length = 100), 
                                                                                     movements_evaluated = UMA_spilldays_movements)

ggsave(here::here("figures", "alternative_spill_treatment", "spilldays", "UMA_wild_compare_movement_spill.png"), UMA_wild_compare_movement_spill, height = 8, width = 8)

# UMA hatchery plot
UMA_hatchery_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Umatilla River",
                                                                                         hatchery_move_prob_array = UMA_hatchery_spill_move_prob_array,
                                                                                         hatchery_covariate_experiences = UMA_hatchery_covariate_experiences,
                                                                                         spill_predict = seq(0, 0.90, length = 100), 
                                                                                         movements_evaluated = UMA_spilldays_movements)

ggsave(here::here("figures", "alternative_spill_treatment", "spilldays", "UMA_hatchery_compare_movement_spill.png"), UMA_hatchery_compare_movement_spill, height = 8, width = 8)

### Yakima River ###

YAK_wild_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Yakima River",
                                                                                     wild_move_prob_array = YAK_wild_spill_move_prob_array,
                                                                                     
                                                                                     wild_covariate_experiences = YAK_wild_covariate_experiences,
                                                                                     spill_predict = seq(0, 0.90, length = 100), 
                                                                                     movements_evaluated = YAK_spilldays_movements)

ggsave(here::here("figures", "alternative_spill_treatment", "spilldays", "YAK_wild_compare_movement_spill.png"), YAK_wild_compare_movement_spill, height = 8, width = 8)



### Walla Walla River ###


# WAWA wild plot
WAWA_wild_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Walla Walla River",
                                                                                      wild_move_prob_array = WAWA_wild_spill_move_prob_array,
                                                                                      wild_covariate_experiences = WAWA_wild_covariate_experiences,
                                                                                      spill_predict = seq(0, 0.90, length = 100), 
                                                                                      movements_evaluated = WAWA_spilldays_movements)

ggsave(here::here("figures", "alternative_spill_treatment", "spilldays", "WAWA_wild_compare_movement_spill.png"), WAWA_wild_compare_movement_spill, height = 8, width = 8)

# WAWA hatchery plot
WAWA_hatchery_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Walla Walla River",
                                                                                          hatchery_move_prob_array = WAWA_hatchery_spill_move_prob_array,
                                                                                          hatchery_covariate_experiences = WAWA_hatchery_covariate_experiences,
                                                                                          spill_predict = seq(0, 0.90, length = 100), 
                                                                                          movements_evaluated = WAWA_spilldays_movements)

ggsave(here::here("figures", "alternative_spill_treatment", "spilldays", "WAWA_hatchery_compare_movement_spill.png"), WAWA_hatchery_compare_movement_spill, height = 8, width = 8)

### Wenatchee River ###


# WEN wild plot
WEN_wild_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Wenatchee River",
                                                                                     wild_move_prob_array = WEN_wild_spill_move_prob_array,
                                                                                     wild_covariate_experiences = WEN_wild_covariate_experiences,
                                                                                     spill_predict = seq(0, 0.90, length = 100), 
                                                                                     movements_evaluated = WEN_spilldays_movements)

ggsave(here::here("figures", "alternative_spill_treatment", "spilldays", "WEN_wild_compare_movement_spill.png"), WEN_wild_compare_movement_spill, height = 8, width = 8)

# WEN hatchery plot
WEN_hatchery_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Wenatchee River",
                                                                                         hatchery_move_prob_array = WEN_hatchery_spill_move_prob_array,
                                                                                         hatchery_covariate_experiences = WEN_hatchery_covariate_experiences,
                                                                                         spill_predict = seq(0, 0.90, length = 100), 
                                                                                         movements_evaluated = WEN_spilldays_movements)

ggsave(here::here("figures", "alternative_spill_treatment", "spilldays", "WEN_hatchery_compare_movement_spill.png"), WEN_hatchery_compare_movement_spill, height = 8, width = 8)

### Entiat River ###


# ENT wild plot
ENT_wild_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Entiat River",
                                                                                     wild_move_prob_array = ENT_wild_spill_move_prob_array,
                                                                                     wild_covariate_experiences = ENT_wild_covariate_experiences,
                                                                                     spill_predict = seq(0, 0.90, length = 100), 
                                                                                     movements_evaluated = ENT_spilldays_movements)

ggsave(here::here("figures", "alternative_spill_treatment", "spilldays", "ENT_wild_compare_movement_spill.png"), ENT_wild_compare_movement_spill, height = 8, width = 8)


### Tucannon River ###


# TUC wild plot
TUC_wild_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Tucannon River",
                                                                                     wild_move_prob_array = TUC_wild_spill_move_prob_array,
                                                                                     wild_covariate_experiences = TUC_wild_covariate_experiences,
                                                                                     spill_predict = seq(0, 0.90, length = 100), 
                                                                                     movements_evaluated = TUC_spilldays_movements)

ggsave(here::here("figures", "alternative_spill_treatment", "spilldays", "TUC_wild_compare_movement_spill.png"), TUC_wild_compare_movement_spill, height = 8, width = 8)

# TUC hatchery plot
TUC_hatchery_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Tucannon River",
                                                                                         hatchery_move_prob_array = TUC_hatchery_spill_move_prob_array,
                                                                                         hatchery_covariate_experiences = TUC_hatchery_covariate_experiences,
                                                                                         spill_predict = seq(0, 0.90, length = 100), 
                                                                                         movements_evaluated = TUC_spilldays_movements)

ggsave(here::here("figures", "alternative_spill_treatment", "spilldays", "TUC_hatchery_compare_movement_spill.png"), TUC_hatchery_compare_movement_spill, height = 8, width = 8)

# ### Imnaha River ###
# IMN_spilldays_movements <- data.frame(from = c(9,9), to = c(8, 41),
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
# ggsave(here::here("figures", "alternative_spill_treatment", "spilldays", "IMN_wild_compare_movement_spill.png"), IMN_wild_compare_movement_spill, height = 8, width = 8)
# 
# # IMN hatchery plot
# IMN_hatchery_compare_movement_spill <- plot_compare_rear_spill_effect_multiple_movements(origin_select = "Imnaha River",
#                                                                                          hatchery_move_prob_array = IMN_hatchery_spill_move_prob_array,
#                                                                                          hatchery_covariate_experiences = IMN_hatchery_covariate_experiences,
#                                                                                          spill_predict = seq(0, 0.90, length = 100), 
#                                                                                          movements_evaluated = IMN_spilldays_movements)
# 
# ggsave(here::here("figures", "alternative_spill_treatment", "spilldays", "IMN_hatchery_compare_movement_spill.png"), IMN_hatchery_compare_movement_spill, height = 8, width = 8)


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
# ggsave(here::here("figures", "alternative_spill_treatment", "paper_figures",  "combined_movement_spill_plot_v1.png"), combined_movement_spill_plot, height = 12, width = 16)

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
                                                                bottom = textGrob("Days of March spill", gp = gpar(cex = 1.3)),
                                                                left = textGrob("Movement probability", rot = 90, gp = gpar(cex = 1.3)))) +
  theme(plot.background = element_rect(fill="white", color = NA))

ggsave(here::here("figures", "alternative_spill_treatment", "paper_figures",  "Fig6_spilldays_plot.png"), combined_movement_spill_plot, height = 12, width = 16)
