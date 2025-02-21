### 07-04_Fig3_final_fate_overshoot

# Figure 3: Compare homing final fate vs. overshoot probability for different origins
# The homing final fate is going to be based on the 08 script that estimates
# final fates in median conditions and takes into account differences for 
# the direction that fish are moving.
# the overshoot probability is going to be based on their overshoot probability
# when they first encounter a dam, so in season1 (this is summer/fall), and 
# when the fish is coming from a downstream state (moving upstream)
# This is the primary time that fish are actually overshooting

# First, need to load in all of the model runs and all of the packages.
source("R/07_model_analysis/07-01_load_stan_models.R")

# 08-model-final-fates.R also needs to be run first, because this will compute
# the final fates based on median conditions across all runs.
load(here::here("figures", "final_fates", "FF_comp_data.rda"))


#### few helpful functions from other scripts ####

get_origin_states_dates <- function(envir, origin_select, rear){
  
  transitions <- envir$data$y
  transition_dates <- envir$data$transition_dates
  # convert the seasons vector to the same shape as the state transitions
  transition_seasons <- envir$data$transition_seasons_vector
  transition_seasons_matrix <- matrix(nrow = nrow(transitions), ncol = ncol(transitions))
  
  movements_counter <- 0
  for (i in 1:nrow(transition_seasons_matrix)){
    # first, count how many observed transitions for that fish
    ntransitions <- length(transitions[i,][which(!(transitions[i,] %in% c(0,41)))])
    # then, take that many of the transitions seasons vector to populate
    transition_seasons_matrix[i,1:ntransitions] <- transition_seasons[(movements_counter+1):(movements_counter+ntransitions)]
    
    # increase the counter
    movements_counter <- movements_counter + ntransitions
  }
  
  
  state_history <- as.vector(t(transitions))
  transition_date_history <- as.vector(t(transition_dates))
  transition_seasons_history <- as.vector(t(transition_seasons_matrix))
  
  # get origin info
  origin_vector <- vector(length = nrow(envir$data$cat_X_mat))
  for(i in 1:nrow(envir$data$cat_X_mat)){
    origin_vector[i] <- which(envir$data$cat_X_mat[i,]==1)
  }
  
  origin_vector_matched <- rep(origin_vector, each = ncol(envir$data$y))
  
  # Now, keep only the origin selected
  if(rear == "wild"){
    origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$wild
  } else {
    origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$hatchery
  }
  
  origin_to_keep <- which(origin_vector_matched == origin_numeric)
  
  state_history <- state_history[origin_to_keep]
  transition_date_history <- transition_date_history[origin_to_keep]
  transition_seasons_history <- transition_seasons_history[origin_to_keep]
  
  # drop the non-data
  nostate <- which(state_history == 0)
  state_history <- state_history[-nostate]
  transition_date_history <- transition_date_history[-nostate]
  transition_seasons_history <- transition_seasons_history[-nostate]
  
  transitions_df <- data.frame(state = rep(NA, length(state_history)-1),
                               previous_state = rep(NA, length(state_history)-1),
                               date = rep(NA, length(state_history)-1),
                               season = rep(NA, length(state_history)-1))
  
  for (i in 1:nrow(transitions_df)){
    if (i == 1){ # first observation has to come from the ether
      transitions_df$state[i] <- state_history[i]
      transitions_df$previous_state[i] <- 0
      transitions_df$date[i] <- transition_date_history[i]
      transitions_df$season[i] <- transition_seasons_history[i]
    } else {
      transitions_df$state[i] <- state_history[i]
      transitions_df$previous_state[i] <- ifelse(state_history[i-1]==41, 0, state_history[i-1]) # all first detections come from state 0 (fish purgatory)
      transitions_df$date[i] <- transition_date_history[i]
      transitions_df$season[i] <- transition_seasons_history[i]
    }
  }
  
  # drop the non-data
  transitions_df %>% 
    filter(!(state %in% c(0,41))) -> transitions_df
  
  # determine if the fish is going upstream or downstream
  transitions_df %>% 
    mutate(direction = ifelse(state > previous_state, "upstream", "downstream")) -> transitions_df
  
  return(transitions_df)
  
}


get_states_dates_direction <- function(envir){
  transitions <- envir$data$y
  transition_dates <- envir$data$transition_dates
  # convert the seasons vector to the same shape as the state transitions
  transition_seasons <- envir$data$transition_seasons_vector
  transition_seasons_matrix <- matrix(nrow = nrow(transitions), ncol = ncol(transitions))
  
  movements_counter <- 0
  for (i in 1:nrow(transition_seasons_matrix)){
    # first, count how many observed transitions for that fish
    ntransitions <- length(transitions[i,][which(!(transitions[i,] %in% c(0,41)))])
    # then, take that many of the transitions seasons vector to populate
    transition_seasons_matrix[i,1:ntransitions] <- transition_seasons[(movements_counter+1):(movements_counter+ntransitions)]
    
    # increase the counter
    movements_counter <- movements_counter + ntransitions
  }
  
  state_history <- as.vector(t(transitions))
  transition_date_history <- as.vector(t(transition_dates))
  transition_seasons_history <- as.vector(t(transition_seasons_matrix))
  
  # drop the non-data
  nostate <- which(state_history == 0)
  state_history <- state_history[-nostate]
  transition_date_history <- transition_date_history[-nostate]
  transition_seasons_history <- transition_seasons_history[-nostate]
  
  transitions_df <- data.frame(state = rep(NA, length(state_history)-1),
                               previous_state = rep(NA, length(state_history)-1),
                               date = rep(NA, length(state_history)-1),
                               season = rep(NA, length(state_history)-1))
  
  for (i in 1:nrow(transitions_df)){
    if (i == 1){ # first observation has to come from the ether
      transitions_df$state[i] <- state_history[i]
      transitions_df$previous_state[i] <- 0
      transitions_df$date[i] <- transition_date_history[i]
      transitions_df$season[i] <- transition_seasons_history[i]
    } else {
      transitions_df$state[i] <- state_history[i]
      transitions_df$previous_state[i] <- ifelse(state_history[i-1]==41, 0, state_history[i-1]) # all first detections come from state 0 (fish purgatory)
      transitions_df$date[i] <- transition_date_history[i]
      transitions_df$season[i] <- transition_seasons_history[i]
    }
  }
  
  # drop the non-data
  transitions_df %>% 
    filter(!(state %in% c(0,41))) -> transitions_df
  
  # determine if the fish is going upstream or downstream
  transitions_df %>% 
    mutate(direction = ifelse(state > previous_state, "upstream", "downstream")) -> transitions_df
  
  return(transitions_df)
}

# create a list that maps origin numbers (params) to what they actually are
natal_origins <- gsub(" Mouth| Upstream", "", model_states)
natal_origins <- natal_origins[!(duplicated(natal_origins))]
natal_origins <- natal_origins[!(grepl("mainstem", natal_origins))]
natal_origins <- natal_origins[!(grepl("other tributaries", natal_origins))]
natal_origins <- natal_origins[!(natal_origins == "loss")]

# Use the parameter map to index the right effects
origin_param_map <- data.frame(
  natal_origin = natal_origins,
  hatchery = c(NA, NA, NA, NA, 1, NA, 2, # MC
               1,NA,2,3, # UC
               5,NA,1,4,2,3), # SR,
  wild = c(1,3,NA,2,4,6,5, # MC
           1,2,NA,3, # UC
           6,1,2,5,3,4)) # SR


#### Functions to compute overshoot probabilities based on median conditions ####

# Write a function that takes an origin, a from state, and computes the probability of movements
# out of that state under median conditions
# This is based primarily on the final fates script: We are going to resample the data
# to 


one_state_movement_probs_UCW <- function(niter = 1000, from_state, states_dates,
                                         origin_select,
                                         origin_param_map = origin_param_map,
                                         temp_data = UCW_envir$data$temperature_data, 
                                         spillwindow_data = UCW_envir$data$spill_window_data, 
                                         winterspill_data = UCW_envir$data$winter_spill_days_data){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- c(0,0,0)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  origin1 <- wild_origin_params[1]
  origin2 <- wild_origin_params[2]
  origin3 <- wild_origin_params[3]
  
  
  movement_prob_data_frame <- data.frame(state = model_states[1:length(model_states)])
  
  ## Step 1: Determine average conditions
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[from_state] <- sample(subset(UCW_states_dates, state == i)$date, 1)
    } else{
      sample_date[from_state] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # Take the median conditions for each state across all years
    winterspill[from_state] <- median(winterspill_data[,i])
  }
  
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[from_state] <- 0
      } else {
        season_upstream[from_state] <- sum(subset(UCW_states_dates, state == i  & direction == "upstream")$season)/length(subset(UCW_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[from_state] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[from_state] <- 0
      } else {
        season_downstream[from_state] <- sum(subset(UCW_states_dates, state == i  & direction == "downstream")$season)/length(subset(UCW_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[from_state] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[from_state] <- 0
      } else {
        temp_upstream_season0[from_state] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[from_state] <- 0
      } else {
        temp_downstream_season0[from_state] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[from_state] <- 0
      } else {
        temp_upstream_season1[from_state] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[from_state] <- 0
      } else {
        temp_downstream_season1[from_state] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[from_state] <- 0
      } else {
        spillwindow_upstream_season0[from_state] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[from_state] <- 0
      } else {
        spillwindow_downstream_season0[from_state] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[from_state] <- 0
      } else {
        spillwindow_upstream_season1[from_state] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[from_state] <- 0
      } else {
        spillwindow_downstream_season1[from_state] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  
  
  ## Step 2: Run through iterations and select different iterations of the model
  
  for(iter in 1:niter) { # iter isn't used to index anything, we're just running through for niter
    # select the iteration you will use
    iter <- sample(1:4000, 1)
    
    
    # make an upstream and a downstream prob matrix
    
    move_prob_vector <- rep(0, length(model_states))
    
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == from_state]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # change from final fates: Only look at upstream season 1 movements
    
    for (j in 1:length(possible_movements)){
      
      move_prob_vector[possible_movements[j]] <- exp(b0_array_UCW[from_state,possible_movements[j],iter] +
                                                       btemp1_array_UCW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state] + 
                                                       bspillwindow_array_UCW[from_state,possible_movements[j],iter]*spillwindow_upstream_season1[from_state] + 
                                                       bwinterspill_array_UCW[from_state,possible_movements[j],iter]*winterspill[from_state] +
                                                       btemp1xorigin1_array_UCW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin1 +
                                                       btemp1xorigin2_array_UCW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin2 + 
                                                       btemp1xorigin3_array_UCW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin3 +
                                                       borigin1_array_UCW[from_state,possible_movements[j],iter]*origin1 +
                                                       borigin2_array_UCW[from_state,possible_movements[j],iter]*origin2 +
                                                       borigin3_array_UCW[from_state,possible_movements[j],iter]*origin3)/
        sum(exp(b0_array_UCW[from_state,possible_movements,iter] +
                  btemp1_array_UCW[from_state,possible_movements,iter]*temp_upstream_season1[from_state] + 
                  bspillwindow_array_UCW[from_state,possible_movements,iter]*spillwindow_upstream_season1[from_state] + 
                  bwinterspill_array_UCW[from_state,possible_movements,iter]*winterspill[from_state] +
                  btemp1xorigin1_array_UCW[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin1 +
                  btemp1xorigin2_array_UCW[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin2 + 
                  btemp1xorigin3_array_UCW[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin3 +
                  borigin1_array_UCW[from_state,possible_movements,iter]*origin1 +
                  borigin2_array_UCW[from_state,possible_movements,iter]*origin2 +
                  borigin3_array_UCW[from_state,possible_movements,iter]*origin3))
      
      
      
    }
    
    movement_prob_data_frame %>% 
      bind_cols(move_prob_vector) -> movement_prob_data_frame
  }
  
  
  
  
  ## BREAK ##
  
  # reformat our output
  rownames(movement_prob_data_frame) <- NULL
  column_to_rownames(movement_prob_data_frame, "state") -> movement_prob_data_frame
  colnames(movement_prob_data_frame) <- paste0("iter", 1:niter)
  
  rownames_to_column(movement_prob_data_frame, "state") -> movement_prob_data_frame
  
  movement_prob_data_frame %>% 
    mutate(state = gsub(" Upstream", "", state)) %>% 
    mutate(state = gsub(" Mouth", "", state)) -> movement_prob_data_frame
  
  movement_prob_data_frame %>% 
    group_by(state) %>% 
    summarise(across(where(is.numeric), sum)) -> movement_prob_data_frame
  
  movement_prob_data_frame %>%
    pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "prob") -> movement_prob_data_frame_long
  
  movement_prob_data_frame_long %>% 
    group_by(state) %>%
    summarise(prob = quantile(prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, values_from = prob) -> movement_prob_data_frame_quantiles
  
  # Return the df that contains movement probabilities under average conditions out of that state
  
  return(movement_prob_data_frame_quantiles)
  
}

one_state_movement_probs_UCH <- function(niter = 1000, from_state, states_dates,
                                         origin_select,
                                         origin_param_map = origin_param_map,
                                         temp_data = UCH_envir$data$temperature_data, 
                                         spillwindow_data = UCH_envir$data$spill_window_data, 
                                         winterspill_data = UCH_envir$data$winter_spill_days_data){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- c(0,0,0)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  origin1 <- hatchery_origin_params[1]
  origin2 <- hatchery_origin_params[2]
  origin3 <- hatchery_origin_params[3]
  
  
  movement_prob_data_frame <- data.frame(state = model_states[1:length(model_states)])
  
  ## Step 1: Determine average conditions
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[from_state] <- sample(subset(UCH_states_dates, state == i)$date, 1)
    } else{
      sample_date[from_state] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # Take the median conditions for each state across all years
    winterspill[from_state] <- median(winterspill_data[,i])
  }
  
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[from_state] <- 0
      } else {
        season_upstream[from_state] <- sum(subset(UCH_states_dates, state == i  & direction == "upstream")$season)/length(subset(UCH_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[from_state] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[from_state] <- 0
      } else {
        season_downstream[from_state] <- sum(subset(UCH_states_dates, state == i  & direction == "downstream")$season)/length(subset(UCH_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[from_state] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[from_state] <- 0
      } else {
        temp_upstream_season0[from_state] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[from_state] <- 0
      } else {
        temp_downstream_season0[from_state] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[from_state] <- 0
      } else {
        temp_upstream_season1[from_state] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[from_state] <- 0
      } else {
        temp_downstream_season1[from_state] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[from_state] <- 0
      } else {
        spillwindow_upstream_season0[from_state] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[from_state] <- 0
      } else {
        spillwindow_downstream_season0[from_state] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[from_state] <- 0
      } else {
        spillwindow_upstream_season1[from_state] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[from_state] <- 0
      } else {
        spillwindow_downstream_season1[from_state] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  
  
  ## Step 2: Run through iterations and select different iterations of the model
  
  for(iter in 1:niter) { # iter isn't used to index anything, we're just running through for niter
    # select the iteration you will use
    iter <- sample(1:4000, 1)
    
    
    # make an upstream and a downstream prob matrix
    
    move_prob_vector <- rep(0, length(model_states))
    
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == from_state]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # change from final fates: Only look at upstream season 1 movements
    
    for (j in 1:length(possible_movements)){
      
      move_prob_vector[possible_movements[j]] <- exp(b0_array_UCH[from_state,possible_movements[j],iter] +
                                                       btemp1_array_UCH[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state] + 
                                                       bspillwindow_array_UCH[from_state,possible_movements[j],iter]*spillwindow_upstream_season1[from_state] + 
                                                       bwinterspill_array_UCH[from_state,possible_movements[j],iter]*winterspill[from_state] +
                                                       btemp1xorigin1_array_UCH[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin1 +
                                                       btemp1xorigin2_array_UCH[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin2 + 
                                                       btemp1xorigin3_array_UCH[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin3 +
                                                       borigin1_array_UCH[from_state,possible_movements[j],iter]*origin1 +
                                                       borigin2_array_UCH[from_state,possible_movements[j],iter]*origin2 +
                                                       borigin3_array_UCH[from_state,possible_movements[j],iter]*origin3)/
        sum(exp(b0_array_UCH[from_state,possible_movements,iter] +
                  btemp1_array_UCH[from_state,possible_movements,iter]*temp_upstream_season1[from_state] + 
                  bspillwindow_array_UCH[from_state,possible_movements,iter]*spillwindow_upstream_season1[from_state] + 
                  bwinterspill_array_UCH[from_state,possible_movements,iter]*winterspill[from_state] +
                  btemp1xorigin1_array_UCH[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin1 +
                  btemp1xorigin2_array_UCH[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin2 + 
                  btemp1xorigin3_array_UCH[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin3 +
                  borigin1_array_UCH[from_state,possible_movements,iter]*origin1 +
                  borigin2_array_UCH[from_state,possible_movements,iter]*origin2 +
                  borigin3_array_UCH[from_state,possible_movements,iter]*origin3))
      
      
      
    }
    
    movement_prob_data_frame %>% 
      bind_cols(move_prob_vector) -> movement_prob_data_frame
  }
  
  
  
  
  ## BREAK ##
  
  # reformat our output
  rownames(movement_prob_data_frame) <- NULL
  column_to_rownames(movement_prob_data_frame, "state") -> movement_prob_data_frame
  colnames(movement_prob_data_frame) <- paste0("iter", 1:niter)
  
  rownames_to_column(movement_prob_data_frame, "state") -> movement_prob_data_frame
  
  movement_prob_data_frame %>% 
    mutate(state = gsub(" Upstream", "", state)) %>% 
    mutate(state = gsub(" Mouth", "", state)) -> movement_prob_data_frame
  
  movement_prob_data_frame %>% 
    group_by(state) %>% 
    summarise(across(where(is.numeric), sum)) -> movement_prob_data_frame
  
  movement_prob_data_frame %>%
    pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "prob") -> movement_prob_data_frame_long
  
  movement_prob_data_frame_long %>% 
    group_by(state) %>%
    summarise(prob = quantile(prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, values_from = prob) -> movement_prob_data_frame_quantiles
  
  # Return the df that contains movement probabilities under average conditions out of that state
  
  return(movement_prob_data_frame_quantiles)
  
}

one_state_movement_probs_MCW <- function(niter = 1000, from_state, states_dates,
                                         origin_select,
                                         origin_param_map = origin_param_map,
                                         temp_data = MCW_envir$data$temperature_data, 
                                         spillwindow_data = MCW_envir$data$spill_window_data, 
                                         winterspill_data = MCW_envir$data$winter_spill_days_data){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- rep(0,6)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  origin1 <- wild_origin_params[1]
  origin2 <- wild_origin_params[2]
  origin3 <- wild_origin_params[3]
  origin4 <- wild_origin_params[4]
  origin5 <- wild_origin_params[5]
  origin6 <- wild_origin_params[6]
  
  
  movement_prob_data_frame <- data.frame(state = model_states[1:length(model_states)])
  
  ## Step 1: Determine average conditions
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[from_state] <- sample(subset(MCW_states_dates, state == i)$date, 1)
    } else{
      sample_date[from_state] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # Take the median conditions for each state across all years
    winterspill[from_state] <- median(winterspill_data[,i])
  }
  
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[from_state] <- 0
      } else {
        season_upstream[from_state] <- sum(subset(MCW_states_dates, state == i  & direction == "upstream")$season)/length(subset(MCW_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[from_state] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[from_state] <- 0
      } else {
        season_downstream[from_state] <- sum(subset(MCW_states_dates, state == i  & direction == "downstream")$season)/length(subset(MCW_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[from_state] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[from_state] <- 0
      } else {
        temp_upstream_season0[from_state] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[from_state] <- 0
      } else {
        temp_downstream_season0[from_state] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[from_state] <- 0
      } else {
        temp_upstream_season1[from_state] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[from_state] <- 0
      } else {
        temp_downstream_season1[from_state] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[from_state] <- 0
      } else {
        spillwindow_upstream_season0[from_state] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[from_state] <- 0
      } else {
        spillwindow_downstream_season0[from_state] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[from_state] <- 0
      } else {
        spillwindow_upstream_season1[from_state] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[from_state] <- 0
      } else {
        spillwindow_downstream_season1[from_state] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  
  
  ## Step 2: Run through iterations and select different iterations of the model
  
  for(iter in 1:niter) { # iter isn't used to index anything, we're just running through for niter
    # select the iteration you will use
    iter <- sample(1:4000, 1)
    
    
    # make an upstream and a downstream prob matrix
    
    move_prob_vector <- rep(0, length(model_states))
    
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == from_state]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # change from final fates: Only look at upstream season 1 movements
    
    for (j in 1:length(possible_movements)){
      
      move_prob_vector[possible_movements[j]] <- exp(b0_array_MCW[from_state,possible_movements[j],iter] +
                                                       btemp1_array_MCW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state] + 
                                                       bspillwindow_array_MCW[from_state,possible_movements[j],iter]*spillwindow_upstream_season1[from_state] + 
                                                       bwinterspill_array_MCW[from_state,possible_movements[j],iter]*winterspill[from_state] +
                                                       btemp1xorigin1_array_MCW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin1 +
                                                       btemp1xorigin2_array_MCW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin2 + 
                                                       btemp1xorigin3_array_MCW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin3 +
                                                       btemp1xorigin4_array_MCW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin4 +
                                                       btemp1xorigin5_array_MCW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin5 + 
                                                       btemp1xorigin6_array_MCW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin6 +
                                                       borigin1_array_MCW[from_state,possible_movements[j],iter]*origin1 +
                                                       borigin2_array_MCW[from_state,possible_movements[j],iter]*origin2 +
                                                       borigin3_array_MCW[from_state,possible_movements[j],iter]*origin3 +
                                                       borigin4_array_MCW[from_state,possible_movements[j],iter]*origin4 +
                                                       borigin5_array_MCW[from_state,possible_movements[j],iter]*origin5 +
                                                       borigin6_array_MCW[from_state,possible_movements[j],iter]*origin6)/
        sum(exp(b0_array_MCW[from_state,possible_movements,iter] +
                  btemp1_array_MCW[from_state,possible_movements,iter]*temp_upstream_season1[from_state] + 
                  bspillwindow_array_MCW[from_state,possible_movements,iter]*spillwindow_upstream_season1[from_state] + 
                  bwinterspill_array_MCW[from_state,possible_movements,iter]*winterspill[from_state] +
                  btemp1xorigin1_array_MCW[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin1 +
                  btemp1xorigin2_array_MCW[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin2 + 
                  btemp1xorigin3_array_MCW[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin3 +
                  btemp1xorigin4_array_MCW[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin4 +
                  btemp1xorigin5_array_MCW[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin5 + 
                  btemp1xorigin6_array_MCW[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin6 +
                  borigin1_array_MCW[from_state,possible_movements,iter]*origin1 +
                  borigin2_array_MCW[from_state,possible_movements,iter]*origin2 +
                  borigin3_array_MCW[from_state,possible_movements,iter]*origin3 +
                  borigin4_array_MCW[from_state,possible_movements,iter]*origin4 +
                  borigin5_array_MCW[from_state,possible_movements,iter]*origin5 +
                  borigin6_array_MCW[from_state,possible_movements,iter]*origin6))
      
      
      
    }
    
    movement_prob_data_frame %>% 
      bind_cols(move_prob_vector) -> movement_prob_data_frame
  }
  
  
  
  
  ## BREAK ##
  
  # reformat our output
  rownames(movement_prob_data_frame) <- NULL
  column_to_rownames(movement_prob_data_frame, "state") -> movement_prob_data_frame
  colnames(movement_prob_data_frame) <- paste0("iter", 1:niter)
  
  rownames_to_column(movement_prob_data_frame, "state") -> movement_prob_data_frame
  
  movement_prob_data_frame %>% 
    mutate(state = gsub(" Upstream", "", state)) %>% 
    mutate(state = gsub(" Mouth", "", state)) -> movement_prob_data_frame
  
  movement_prob_data_frame %>% 
    group_by(state) %>% 
    summarise(across(where(is.numeric), sum)) -> movement_prob_data_frame
  
  movement_prob_data_frame %>%
    pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "prob") -> movement_prob_data_frame_long
  
  movement_prob_data_frame_long %>% 
    group_by(state) %>%
    summarise(prob = quantile(prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, values_from = prob) -> movement_prob_data_frame_quantiles
  
  # Return the df that contains movement probabilities under average conditions out of that state
  
  return(movement_prob_data_frame_quantiles)
  
}

one_state_movement_probs_MCH <- function(niter = 1000, from_state, states_dates,
                                         origin_select,
                                         origin_param_map = origin_param_map,
                                         temp_data = MCH_envir$data$temperature_data, 
                                         spillwindow_data = MCH_envir$data$spill_window_data, 
                                         winterspill_data = MCH_envir$data$winter_spill_days_data){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,2)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  origin1 <- hatchery_origin_params[1]
  origin2 <- hatchery_origin_params[2]
  
  
  movement_prob_data_frame <- data.frame(state = model_states[1:length(model_states)])
  
  ## Step 1: Determine average conditions
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[from_state] <- sample(subset(MCH_states_dates, state == i)$date, 1)
    } else{
      sample_date[from_state] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # Take the median conditions for each state across all years
    winterspill[from_state] <- median(winterspill_data[,i])
  }
  
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[from_state] <- 0
      } else {
        season_upstream[from_state] <- sum(subset(MCH_states_dates, state == i  & direction == "upstream")$season)/length(subset(MCH_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[from_state] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[from_state] <- 0
      } else {
        season_downstream[from_state] <- sum(subset(MCH_states_dates, state == i  & direction == "downstream")$season)/length(subset(MCH_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[from_state] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[from_state] <- 0
      } else {
        temp_upstream_season0[from_state] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[from_state] <- 0
      } else {
        temp_downstream_season0[from_state] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[from_state] <- 0
      } else {
        temp_upstream_season1[from_state] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[from_state] <- 0
      } else {
        temp_downstream_season1[from_state] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[from_state] <- 0
      } else {
        spillwindow_upstream_season0[from_state] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[from_state] <- 0
      } else {
        spillwindow_downstream_season0[from_state] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[from_state] <- 0
      } else {
        spillwindow_upstream_season1[from_state] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[from_state] <- 0
      } else {
        spillwindow_downstream_season1[from_state] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  
  
  ## Step 2: Run through iterations and select different iterations of the model
  
  for(iter in 1:niter) { # iter isn't used to index anything, we're just running through for niter
    # select the iteration you will use
    iter <- sample(1:4000, 1)
    
    
    # make an upstream and a downstream prob matrix
    
    move_prob_vector <- rep(0, length(model_states))
    
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == from_state]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # change from final fates: Only look at upstream season 1 movements
    
    for (j in 1:length(possible_movements)){
      
      move_prob_vector[possible_movements[j]] <- exp(b0_array_MCH[from_state,possible_movements[j],iter] +
                                                       btemp1_array_MCH[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state] + 
                                                       bspillwindow_array_MCH[from_state,possible_movements[j],iter]*spillwindow_upstream_season1[from_state] + 
                                                       bwinterspill_array_MCH[from_state,possible_movements[j],iter]*winterspill[from_state] +
                                                       btemp1xorigin1_array_MCH[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin1 +
                                                       btemp1xorigin2_array_MCH[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin2 + 
                                                       borigin1_array_MCH[from_state,possible_movements[j],iter]*origin1 +
                                                       borigin2_array_MCH[from_state,possible_movements[j],iter]*origin2)/
        sum(exp(b0_array_MCH[from_state,possible_movements,iter] +
                  btemp1_array_MCH[from_state,possible_movements,iter]*temp_upstream_season1[from_state] + 
                  bspillwindow_array_MCH[from_state,possible_movements,iter]*spillwindow_upstream_season1[from_state] + 
                  bwinterspill_array_MCH[from_state,possible_movements,iter]*winterspill[from_state] +
                  btemp1xorigin1_array_MCH[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin1 +
                  btemp1xorigin2_array_MCH[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin2 + 
                  borigin1_array_MCH[from_state,possible_movements,iter]*origin1 +
                  borigin2_array_MCH[from_state,possible_movements,iter]*origin2))
      
      
      
    }
    
    movement_prob_data_frame %>% 
      bind_cols(move_prob_vector) -> movement_prob_data_frame
  }
  
  
  
  
  ## BREAK ##
  
  # reformat our output
  rownames(movement_prob_data_frame) <- NULL
  column_to_rownames(movement_prob_data_frame, "state") -> movement_prob_data_frame
  colnames(movement_prob_data_frame) <- paste0("iter", 1:niter)
  
  rownames_to_column(movement_prob_data_frame, "state") -> movement_prob_data_frame
  
  movement_prob_data_frame %>% 
    mutate(state = gsub(" Upstream", "", state)) %>% 
    mutate(state = gsub(" Mouth", "", state)) -> movement_prob_data_frame
  
  movement_prob_data_frame %>% 
    group_by(state) %>% 
    summarise(across(where(is.numeric), sum)) -> movement_prob_data_frame
  
  movement_prob_data_frame %>%
    pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "prob") -> movement_prob_data_frame_long
  
  movement_prob_data_frame_long %>% 
    group_by(state) %>%
    summarise(prob = quantile(prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, values_from = prob) -> movement_prob_data_frame_quantiles
  
  # Return the df that contains movement probabilities under average conditions out of that state
  
  return(movement_prob_data_frame_quantiles)
  
}

one_state_movement_probs_SRW <- function(niter = 1000, from_state, states_dates,
                                         origin_select,
                                         origin_param_map = origin_param_map,
                                         temp_data = SRW_envir$data$temperature_data, 
                                         spillwindow_data = SRW_envir$data$spill_window_data, 
                                         winterspill_data = SRW_envir$data$winter_spill_days_data){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- rep(0,6)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  origin1 <- wild_origin_params[1]
  origin2 <- wild_origin_params[2]
  origin3 <- wild_origin_params[3]
  origin4 <- wild_origin_params[4]
  origin5 <- wild_origin_params[5]
  origin6 <- wild_origin_params[6]
  
  
  movement_prob_data_frame <- data.frame(state = model_states[1:length(model_states)])
  
  ## Step 1: Determine average conditions
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[from_state] <- sample(subset(SRW_states_dates, state == i)$date, 1)
    } else{
      sample_date[from_state] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # Take the median conditions for each state across all years
    winterspill[from_state] <- median(winterspill_data[,i])
  }
  
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[from_state] <- 0
      } else {
        season_upstream[from_state] <- sum(subset(SRW_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRW_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[from_state] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[from_state] <- 0
      } else {
        season_downstream[from_state] <- sum(subset(SRW_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRW_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[from_state] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[from_state] <- 0
      } else {
        temp_upstream_season0[from_state] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[from_state] <- 0
      } else {
        temp_downstream_season0[from_state] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[from_state] <- 0
      } else {
        temp_upstream_season1[from_state] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[from_state] <- 0
      } else {
        temp_downstream_season1[from_state] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[from_state] <- 0
      } else {
        spillwindow_upstream_season0[from_state] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[from_state] <- 0
      } else {
        spillwindow_downstream_season0[from_state] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[from_state] <- 0
      } else {
        spillwindow_upstream_season1[from_state] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[from_state] <- 0
      } else {
        spillwindow_downstream_season1[from_state] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  
  
  ## Step 2: Run through iterations and select different iterations of the model
  
  for(iter in 1:niter) { # iter isn't used to index anything, we're just running through for niter
    # select the iteration you will use
    iter <- sample(1:4000, 1)
    
    
    # make an upstream and a downstream prob matrix
    
    move_prob_vector <- rep(0, length(model_states))
    
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == from_state]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # change from final fates: Only look at upstream season 1 movements
    
    for (j in 1:length(possible_movements)){
      
      move_prob_vector[possible_movements[j]] <- exp(b0_array_SRW[from_state,possible_movements[j],iter] +
                                                       btemp1_array_SRW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state] + 
                                                       bspillwindow_array_SRW[from_state,possible_movements[j],iter]*spillwindow_upstream_season1[from_state] + 
                                                       bwinterspill_array_SRW[from_state,possible_movements[j],iter]*winterspill[from_state] +
                                                       btemp1xorigin1_array_SRW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin1 +
                                                       btemp1xorigin2_array_SRW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin2 + 
                                                       btemp1xorigin3_array_SRW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin3 +
                                                       btemp1xorigin4_array_SRW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin4 +
                                                       btemp1xorigin5_array_SRW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin5 + 
                                                       btemp1xorigin6_array_SRW[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin6 +
                                                       borigin1_array_SRW[from_state,possible_movements[j],iter]*origin1 +
                                                       borigin2_array_SRW[from_state,possible_movements[j],iter]*origin2 +
                                                       borigin3_array_SRW[from_state,possible_movements[j],iter]*origin3 +
                                                       borigin4_array_SRW[from_state,possible_movements[j],iter]*origin4 +
                                                       borigin5_array_SRW[from_state,possible_movements[j],iter]*origin5 +
                                                       borigin6_array_SRW[from_state,possible_movements[j],iter]*origin6)/
        sum(exp(b0_array_SRW[from_state,possible_movements,iter] +
                  btemp1_array_SRW[from_state,possible_movements,iter]*temp_upstream_season1[from_state] + 
                  bspillwindow_array_SRW[from_state,possible_movements,iter]*spillwindow_upstream_season1[from_state] + 
                  bwinterspill_array_SRW[from_state,possible_movements,iter]*winterspill[from_state] +
                  btemp1xorigin1_array_SRW[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin1 +
                  btemp1xorigin2_array_SRW[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin2 + 
                  btemp1xorigin3_array_SRW[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin3 +
                  btemp1xorigin4_array_SRW[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin4 +
                  btemp1xorigin5_array_SRW[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin5 + 
                  btemp1xorigin6_array_SRW[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin6 +
                  borigin1_array_SRW[from_state,possible_movements,iter]*origin1 +
                  borigin2_array_SRW[from_state,possible_movements,iter]*origin2 +
                  borigin3_array_SRW[from_state,possible_movements,iter]*origin3 +
                  borigin4_array_SRW[from_state,possible_movements,iter]*origin4 +
                  borigin5_array_SRW[from_state,possible_movements,iter]*origin5 +
                  borigin6_array_SRW[from_state,possible_movements,iter]*origin6))
      
      
      
    }
    
    movement_prob_data_frame %>% 
      bind_cols(move_prob_vector) -> movement_prob_data_frame
  }
  
  
  
  
  ## BREAK ##
  
  # reformat our output
  rownames(movement_prob_data_frame) <- NULL
  column_to_rownames(movement_prob_data_frame, "state") -> movement_prob_data_frame
  colnames(movement_prob_data_frame) <- paste0("iter", 1:niter)
  
  rownames_to_column(movement_prob_data_frame, "state") -> movement_prob_data_frame
  
  movement_prob_data_frame %>% 
    mutate(state = gsub(" Upstream", "", state)) %>% 
    mutate(state = gsub(" Mouth", "", state)) -> movement_prob_data_frame
  
  movement_prob_data_frame %>% 
    group_by(state) %>% 
    summarise(across(where(is.numeric), sum)) -> movement_prob_data_frame
  
  movement_prob_data_frame %>%
    pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "prob") -> movement_prob_data_frame_long
  
  movement_prob_data_frame_long %>% 
    group_by(state) %>%
    summarise(prob = quantile(prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, values_from = prob) -> movement_prob_data_frame_quantiles
  
  # Return the df that contains movement probabilities under average conditions out of that state
  
  return(movement_prob_data_frame_quantiles)
  
}

one_state_movement_probs_SRH <- function(niter = 1000, from_state, states_dates,
                                         origin_select,
                                         origin_param_map = origin_param_map,
                                         temp_data = SRH_envir$data$temperature_data, 
                                         spillwindow_data = SRH_envir$data$spill_window_data, 
                                         winterspill_data = SRH_envir$data$winter_spill_days_data){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,5)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  origin1 <- hatchery_origin_params[1]
  origin2 <- hatchery_origin_params[2]
  origin3 <- hatchery_origin_params[3]
  origin4 <- hatchery_origin_params[4]
  origin5 <- hatchery_origin_params[5]
  
  
  movement_prob_data_frame <- data.frame(state = model_states[1:length(model_states)])
  
  ## Step 1: Determine average conditions
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[from_state] <- sample(subset(SRH_states_dates, state == i)$date, 1)
    } else{
      sample_date[from_state] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # Take the median conditions for each state across all years
    winterspill[from_state] <- median(winterspill_data[,i])
  }
  
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[from_state] <- 0
      } else {
        season_upstream[from_state] <- sum(subset(SRH_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRH_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[from_state] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[from_state] <- 0
      } else {
        season_downstream[from_state] <- sum(subset(SRH_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRH_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[from_state] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[from_state] <- 0
      } else {
        temp_upstream_season0[from_state] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[from_state] <- 0
      } else {
        temp_downstream_season0[from_state] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[from_state] <- 0
      } else {
        temp_upstream_season1[from_state] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[from_state] <- 0
      } else {
        temp_downstream_season1[from_state] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[from_state] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[from_state] <- 0
      } else {
        spillwindow_upstream_season0[from_state] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[from_state] <- 0
      } else {
        spillwindow_downstream_season0[from_state] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[from_state] <- 0
      } else {
        spillwindow_upstream_season1[from_state] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[from_state] <- 0
      } else {
        spillwindow_downstream_season1[from_state] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[from_state] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  
  
  ## Step 2: Run through iterations and select different iterations of the model
  
  for(iter in 1:niter) { # iter isn't used to index anything, we're just running through for niter
    # select the iteration you will use
    iter <- sample(1:4000, 1)
    
    
    # make an upstream and a downstream prob matrix
    
    move_prob_vector <- rep(0, length(model_states))
    
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == from_state]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # change from final fates: Only look at upstream season 1 movements
    
    for (j in 1:length(possible_movements)){
      
      move_prob_vector[possible_movements[j]] <- exp(b0_array_SRH[from_state,possible_movements[j],iter] +
                                                       btemp1_array_SRH[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state] + 
                                                       bspillwindow_array_SRH[from_state,possible_movements[j],iter]*spillwindow_upstream_season1[from_state] + 
                                                       bwinterspill_array_SRH[from_state,possible_movements[j],iter]*winterspill[from_state] +
                                                       btemp1xorigin1_array_SRH[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin1 +
                                                       btemp1xorigin2_array_SRH[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin2 + 
                                                       btemp1xorigin3_array_SRH[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin3 +
                                                       btemp1xorigin4_array_SRH[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin4 +
                                                       btemp1xorigin5_array_SRH[from_state,possible_movements[j],iter]*temp_upstream_season1[from_state]*origin5 + 
                                                       borigin1_array_SRH[from_state,possible_movements[j],iter]*origin1 +
                                                       borigin2_array_SRH[from_state,possible_movements[j],iter]*origin2 +
                                                       borigin3_array_SRH[from_state,possible_movements[j],iter]*origin3 +
                                                       borigin4_array_SRH[from_state,possible_movements[j],iter]*origin4 +
                                                       borigin5_array_SRH[from_state,possible_movements[j],iter]*origin5)/
        sum(exp(b0_array_SRH[from_state,possible_movements,iter] +
                  btemp1_array_SRH[from_state,possible_movements,iter]*temp_upstream_season1[from_state] + 
                  bspillwindow_array_SRH[from_state,possible_movements,iter]*spillwindow_upstream_season1[from_state] + 
                  bwinterspill_array_SRH[from_state,possible_movements,iter]*winterspill[from_state] +
                  btemp1xorigin1_array_SRH[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin1 +
                  btemp1xorigin2_array_SRH[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin2 + 
                  btemp1xorigin3_array_SRH[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin3 +
                  btemp1xorigin4_array_SRH[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin4 +
                  btemp1xorigin5_array_SRH[from_state,possible_movements,iter]*temp_upstream_season1[from_state]*origin5 + 
                  borigin1_array_SRH[from_state,possible_movements,iter]*origin1 +
                  borigin2_array_SRH[from_state,possible_movements,iter]*origin2 +
                  borigin3_array_SRH[from_state,possible_movements,iter]*origin3 +
                  borigin4_array_SRH[from_state,possible_movements,iter]*origin4 +
                  borigin5_array_SRH[from_state,possible_movements,iter]*origin5))
      
      
      
    }
    
    movement_prob_data_frame %>% 
      bind_cols(move_prob_vector) -> movement_prob_data_frame
  }
  
  
  
  
  ## BREAK ##
  
  # reformat our output
  rownames(movement_prob_data_frame) <- NULL
  column_to_rownames(movement_prob_data_frame, "state") -> movement_prob_data_frame
  colnames(movement_prob_data_frame) <- paste0("iter", 1:niter)
  
  rownames_to_column(movement_prob_data_frame, "state") -> movement_prob_data_frame
  
  movement_prob_data_frame %>% 
    mutate(state = gsub(" Upstream", "", state)) %>% 
    mutate(state = gsub(" Mouth", "", state)) -> movement_prob_data_frame
  
  movement_prob_data_frame %>% 
    group_by(state) %>% 
    summarise(across(where(is.numeric), sum)) -> movement_prob_data_frame
  
  movement_prob_data_frame %>%
    pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "prob") -> movement_prob_data_frame_long
  
  movement_prob_data_frame_long %>% 
    group_by(state) %>%
    summarise(prob = quantile(prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, values_from = prob) -> movement_prob_data_frame_quantiles
  
  # Return the df that contains movement probabilities under average conditions out of that state
  
  return(movement_prob_data_frame_quantiles)
  
}

# Function to take that output and reformat it only as the overshoot
reformat_overshoot_prob <- function(movement_probs, overshoot_state, origin_select,
                                    rear_type){
  subset(movement_probs, state == overshoot_state) -> overshoot_only
  overshoot_only %>% 
    mutate(origin = origin_select,
           rear = rear_type) -> overshoot_only
  return(overshoot_only)
}

# Code below keeps crashing, maybe because we're out of memory. Let's remove a bunch of things and try again.
# rm(UCW_fit_summary)
# rm(UCW_fit_raw)
# rm(UCH_fit_summary)
# rm(UCH_fit_raw)
# rm(MCW_fit_summary)
# rm(MCW_fit_raw)
# rm(MCH_fit_summary)
# rm(MCH_fit_raw)
# rm(SRW_fit_summary)
# rm(SRW_fit_raw)
# rm(SRH_fit_summary)
# rm(SRH_fit_raw)



#### Run functions to calculate overshoot probabilities ####

# Run for Upper Columbia, wild populations
UCW_states_dates <- get_states_dates_direction(envir = UCW_envir)

WEN_wild_home_state_probs <- one_state_movement_probs_UCW(niter = 1000, from_state = 5, states_dates = get_origin_states_dates(envir = UCW_envir, origin_select = "Wenatchee River", rear = "wild"),
                                                          origin_select = "Wenatchee River",
                                                          origin_param_map = origin_param_map,
                                                          temp_data = UCW_envir$data$temperature_data, 
                                                          spillwindow_data = UCW_envir$data$spill_window_data, 
                                                          winterspill_data = UCW_envir$data$winter_spill_days_data)

WEN_wild_overshoot <- reformat_overshoot_prob(WEN_wild_home_state_probs, 
                                              overshoot_state = "mainstem, RRE to WEL",
                                              origin = "Wenatchee River",
                                              rear_type = "wild")

ENT_wild_home_state_probs <- one_state_movement_probs_UCW(niter = 1000, from_state = 6, states_dates = get_origin_states_dates(envir = UCW_envir, origin_select = "Entiat River", rear = "wild"),
                                                          origin_select = "Entiat River",
                                                          origin_param_map = origin_param_map,
                                                          temp_data = UCW_envir$data$temperature_data, 
                                                          spillwindow_data = UCW_envir$data$spill_window_data, 
                                                          winterspill_data = UCW_envir$data$winter_spill_days_data)

ENT_wild_overshoot <- reformat_overshoot_prob(ENT_wild_home_state_probs, 
                                              overshoot_state = "mainstem, upstream of WEL",
                                              origin = "Entiat River",
                                              rear_type = "wild")

# rm(UCW_states_dates)
# Run for Upper Columbia, hatchery populations
UCH_states_dates <- get_states_dates_direction(envir = UCH_envir)

WEN_hatchery_home_state_probs <- one_state_movement_probs_UCH(niter = 1000, from_state = 5, states_dates = get_origin_states_dates(envir = UCH_envir, origin_select = "Wenatchee River", rear = "hatchery"),
                                                              origin_select = "Wenatchee River",
                                                              origin_param_map = origin_param_map,
                                                              temp_data = UCH_envir$data$temperature_data, 
                                                              spillwindow_data = UCH_envir$data$spill_window_data, 
                                                              winterspill_data = UCH_envir$data$winter_spill_days_data)

WEN_hatchery_overshoot <- reformat_overshoot_prob(WEN_hatchery_home_state_probs, 
                                                  overshoot_state = "mainstem, RRE to WEL",
                                                  origin = "Wenatchee River",
                                                  rear_type = "hatchery")

# rm(UCH_states_dates)
# Run for Middle Columbia, wild populations
MCW_states_dates <- get_states_dates_direction(envir = MCW_envir)

DES_wild_home_state_probs <- one_state_movement_probs_MCW(niter = 1000, from_state = 2, states_dates = get_origin_states_dates(envir = MCW_envir, origin_select = "Deschutes River", rear = "wild"),
                                                          origin_select = "Deschutes River",
                                                          origin_param_map = origin_param_map,
                                                          temp_data = MCW_envir$data$temperature_data, 
                                                          spillwindow_data = MCW_envir$data$spill_window_data, 
                                                          winterspill_data = MCW_envir$data$winter_spill_days_data)

DES_wild_overshoot <- reformat_overshoot_prob(DES_wild_home_state_probs, 
                                              overshoot_state = "mainstem, MCN to ICH or PRA",
                                              origin = "Deschutes River",
                                              rear_type = "wild")

FIF_wild_home_state_probs <- one_state_movement_probs_MCW(niter = 1000, from_state = 2, states_dates = get_origin_states_dates(envir = MCW_envir, origin_select = "Fifteenmile Creek", rear = "wild"),
                                                          origin_select = "Fifteenmile Creek",
                                                          origin_param_map = origin_param_map,
                                                          temp_data = MCW_envir$data$temperature_data, 
                                                          spillwindow_data = MCW_envir$data$spill_window_data, 
                                                          winterspill_data = MCW_envir$data$winter_spill_days_data)

FIF_wild_overshoot <- reformat_overshoot_prob(FIF_wild_home_state_probs, 
                                              overshoot_state = "mainstem, MCN to ICH or PRA",
                                              origin = "Fifteenmile Creek",
                                              rear_type = "wild")

JDR_wild_home_state_probs <- one_state_movement_probs_MCW(niter = 1000, from_state = 2, states_dates = get_origin_states_dates(envir = MCW_envir, origin_select = "John Day River", rear = "wild"),
                                                          origin_select = "John Day River",
                                                          origin_param_map = origin_param_map,
                                                          temp_data = MCW_envir$data$temperature_data, 
                                                          spillwindow_data = MCW_envir$data$spill_window_data, 
                                                          winterspill_data = MCW_envir$data$winter_spill_days_data)

JDR_wild_overshoot <- reformat_overshoot_prob(JDR_wild_home_state_probs, 
                                              overshoot_state = "mainstem, MCN to ICH or PRA",
                                              origin = "John Day River",
                                              rear_type = "wild")

UMA_wild_home_state_probs <- one_state_movement_probs_MCW(niter = 1000, from_state = 2, states_dates = get_origin_states_dates(envir = MCW_envir, origin_select = "Umatilla River", rear = "wild"),
                                                          origin_select = "Umatilla River",
                                                          origin_param_map = origin_param_map,
                                                          temp_data = MCW_envir$data$temperature_data, 
                                                          spillwindow_data = MCW_envir$data$spill_window_data, 
                                                          winterspill_data = MCW_envir$data$winter_spill_days_data)

UMA_wild_overshoot <- reformat_overshoot_prob(UMA_wild_home_state_probs, 
                                              overshoot_state = "mainstem, MCN to ICH or PRA",
                                              origin = "Umatilla River",
                                              rear_type = "wild")

WAWA_wild_home_state_probs <- one_state_movement_probs_MCW(niter = 1000, from_state = 3, states_dates = get_origin_states_dates(envir = MCW_envir, origin_select = "Walla Walla River", rear = "wild"),
                                                           origin_select = "Walla Walla River",
                                                           origin_param_map = origin_param_map,
                                                           temp_data = MCW_envir$data$temperature_data, 
                                                           spillwindow_data = MCW_envir$data$spill_window_data, 
                                                           winterspill_data = MCW_envir$data$winter_spill_days_data)

WAWA_wild_overshoot_PRA <- reformat_overshoot_prob(WAWA_wild_home_state_probs, 
                                                   overshoot_state = "mainstem, PRA to RIS",
                                                   origin = "Walla Walla River",
                                                   rear_type = "wild")

WAWA_wild_overshoot_ICH <- reformat_overshoot_prob(WAWA_wild_home_state_probs, 
                                                   overshoot_state = "mainstem, ICH to LGR",
                                                   origin = "Walla Walla River",
                                                   rear_type = "wild")

YAK_wild_home_state_probs <- one_state_movement_probs_MCW(niter = 1000, from_state = 3, states_dates = get_origin_states_dates(envir = MCW_envir, origin_select = "Yakima River", rear = "wild"),
                                                          origin_select = "Yakima River",
                                                          origin_param_map = origin_param_map,
                                                          temp_data = MCW_envir$data$temperature_data, 
                                                          spillwindow_data = MCW_envir$data$spill_window_data, 
                                                          winterspill_data = MCW_envir$data$winter_spill_days_data)

YAK_wild_overshoot_PRA <- reformat_overshoot_prob(YAK_wild_home_state_probs, 
                                                  overshoot_state = "mainstem, PRA to RIS",
                                                  origin = "Yakima River",
                                                  rear_type = "wild")

YAK_wild_overshoot_ICH <- reformat_overshoot_prob(YAK_wild_home_state_probs, 
                                                  overshoot_state = "mainstem, ICH to LGR",
                                                  origin = "Yakima River",
                                                  rear_type = "wild")
# rm(MCW_states_dates)

# Run for Middle Columbia, hatchery populations
MCH_states_dates <- get_states_dates_direction(envir = MCH_envir)

UMA_hatchery_home_state_probs <- one_state_movement_probs_MCH(niter = 1000, from_state = 2, states_dates = get_origin_states_dates(envir = MCH_envir, origin_select = "Umatilla River", rear = "hatchery"),
                                                              origin_select = "Umatilla River",
                                                              origin_param_map = origin_param_map,
                                                              temp_data = MCH_envir$data$temperature_data, 
                                                              spillwindow_data = MCH_envir$data$spill_window_data, 
                                                              winterspill_data = MCH_envir$data$winter_spill_days_data)

UMA_hatchery_overshoot <- reformat_overshoot_prob(UMA_hatchery_home_state_probs, 
                                                  overshoot_state = "mainstem, MCN to ICH or PRA",
                                                  origin = "Umatilla River",
                                                  rear_type = "hatchery")

WAWA_hatchery_home_state_probs <- one_state_movement_probs_MCH(niter = 1000, from_state = 3, states_dates = get_origin_states_dates(envir = MCH_envir, origin_select = "Walla Walla River", rear = "hatchery"),
                                                               origin_select = "Walla Walla River",
                                                               origin_param_map = origin_param_map,
                                                               temp_data = MCH_envir$data$temperature_data, 
                                                               spillwindow_data = MCH_envir$data$spill_window_data, 
                                                               winterspill_data = MCH_envir$data$winter_spill_days_data)

WAWA_hatchery_overshoot_PRA <- reformat_overshoot_prob(WAWA_hatchery_home_state_probs, 
                                                       overshoot_state = "mainstem, PRA to RIS",
                                                       origin = "Walla Walla River",
                                                       rear_type = "hatchery")

WAWA_hatchery_overshoot_ICH <- reformat_overshoot_prob(WAWA_hatchery_home_state_probs, 
                                                       overshoot_state = "mainstem, ICH to LGR",
                                                       origin = "Walla Walla River",
                                                       rear_type = "hatchery")

# rm(MCH_states_dates)

# Run for Snake River, wild populations
SRW_states_dates <- get_states_dates_direction(envir = SRW_envir)

TUC_wild_home_state_probs <- one_state_movement_probs_SRW(niter = 1000, from_state = 8, states_dates = get_origin_states_dates(envir = SRW_envir, origin_select = "Tucannon River", rear = "wild"),
                                                          origin_select = "Tucannon River",
                                                          origin_param_map = origin_param_map,
                                                          temp_data = SRW_envir$data$temperature_data, 
                                                          spillwindow_data = SRW_envir$data$spill_window_data, 
                                                          winterspill_data = SRW_envir$data$winter_spill_days_data)

TUC_wild_overshoot <- reformat_overshoot_prob(TUC_wild_home_state_probs, 
                                              overshoot_state = "mainstem, upstream of LGR",
                                              origin = "Tucannon River",
                                              rear_type = "wild")

# rm(SRW_states_dates)
# Run for Snake River, hatchery populations
SRH_states_dates <- get_states_dates_direction(envir = SRH_envir)

TUC_hatchery_home_state_probs <- one_state_movement_probs_SRH(niter = 1000, from_state = 8, states_dates = get_origin_states_dates(envir = SRH_envir, origin_select = "Tucannon River", rear = "hatchery"),
                                                              origin_select = "Tucannon River",
                                                              origin_param_map = origin_param_map,
                                                              temp_data = SRH_envir$data$temperature_data, 
                                                              spillwindow_data = SRH_envir$data$spill_window_data, 
                                                              winterspill_data = SRH_envir$data$winter_spill_days_data)

TUC_hatchery_overshoot <- reformat_overshoot_prob(TUC_hatchery_home_state_probs, 
                                                  overshoot_state = "mainstem, upstream of LGR",
                                                  origin = "Tucannon River",
                                                  rear_type = "hatchery")

# rm(SRH_states_dates)


#### Get final fates (homing) information ####
wen_ff_comp_median <- FF_comp_data$wen_ff_comp_median
wen_ff_comp_median$origin <- "Wenatchee River"
ent_ff_comp_median <- FF_comp_data$ent_ff_comp_median
ent_ff_comp_median$origin <- "Entiat River"
oka_ff_comp_median <- FF_comp_data$oka_ff_comp_median
oka_ff_comp_median$origin <- "Okanogan River"
met_ff_comp_median <- FF_comp_data$met_ff_comp_median
met_ff_comp_median$origin <- "Methow River"
des_ff_comp_median <- FF_comp_data$des_ff_comp_median
des_ff_comp_median$origin <- "Deschutes River"
jdr_ff_comp_median <- FF_comp_data$jdr_ff_comp_median
jdr_ff_comp_median$origin <- "John Day River"
fif_ff_comp_median <- FF_comp_data$fif_ff_comp_median
fif_ff_comp_median$origin <- "Fifteenmile Creek"
uma_ff_comp_median <- FF_comp_data$uma_ff_comp_median
uma_ff_comp_median$origin <- "Umatilla River"
yak_ff_comp_median <- FF_comp_data$yak_ff_comp_median
yak_ff_comp_median$origin <- "Yakima River"
wawa_ff_comp_median <- FF_comp_data$wawa_ff_comp_median
wawa_ff_comp_median$origin <- "Walla Walla River"
aso_ff_comp_median <- FF_comp_data$aso_ff_comp_median
aso_ff_comp_median$origin <- "Asotin Creek"
cle_ff_comp_median <- FF_comp_data$cle_ff_comp_median
cle_ff_comp_median$origin <- "Clearwater River"
sal_ff_comp_median <- FF_comp_data$sal_ff_comp_median
sal_ff_comp_median$origin <- "Salmon River"
gr_ff_comp_median <- FF_comp_data$gr_ff_comp_median
gr_ff_comp_median$origin <- "Grande Ronde River"
imn_ff_comp_median <- FF_comp_data$imn_ff_comp_median
imn_ff_comp_median$origin <- "Imnaha River"
tuc_ff_comp_median <- FF_comp_data$tuc_ff_comp_median
tuc_ff_comp_median$origin <- "Tucannon River"

wen_ff_comp_median  %>% 
  bind_rows(., ent_ff_comp_median)  %>% 
  bind_rows(., oka_ff_comp_median)   %>% 
  bind_rows(., met_ff_comp_median) %>% 
  bind_rows(., des_ff_comp_median) %>% 
  bind_rows(., jdr_ff_comp_median) %>% 
  bind_rows(., fif_ff_comp_median) %>% 
  bind_rows(., uma_ff_comp_median) %>% 
  bind_rows(., yak_ff_comp_median) %>% 
  bind_rows(., wawa_ff_comp_median) %>% 
  bind_rows(., aso_ff_comp_median) %>% 
  bind_rows(., cle_ff_comp_median) %>% 
  bind_rows(., sal_ff_comp_median) %>% 
  bind_rows(., gr_ff_comp_median) %>% 
  bind_rows(., imn_ff_comp_median) %>% 
  bind_rows(., tuc_ff_comp_median) -> FF_all

FF_all %>% 
  ungroup() %>% 
  mutate(homing_state = ifelse(state == origin, TRUE, FALSE)) %>% 
  filter(homing_state == TRUE) %>% 
  dplyr::select(-homing_state) %>% 
  dplyr::select(-state) %>% 
  dplyr::rename(rear = rear_type) %>% 
  dplyr::rename(homing_lower = `0.025`,
                homing_median = `0.5`,
                homing_upper = `0.975`) -> FF_homing



#### Combine outputs ####
WEN_wild_overshoot %>% 
  bind_rows(., ENT_wild_overshoot) %>% 
  bind_rows(., WEN_hatchery_overshoot) %>% 
  bind_rows(., DES_wild_overshoot) %>% 
  bind_rows(., FIF_wild_overshoot) %>% 
  bind_rows(., JDR_wild_overshoot) %>% 
  bind_rows(., UMA_wild_overshoot) %>% 
  bind_rows(., WAWA_wild_overshoot_PRA) %>% 
  bind_rows(., WAWA_wild_overshoot_ICH) %>% 
  bind_rows(., YAK_wild_overshoot_PRA) %>% 
  bind_rows(., YAK_wild_overshoot_ICH) %>% 
  bind_rows(., UMA_hatchery_overshoot) %>% 
  bind_rows(., WAWA_hatchery_overshoot_PRA) %>% 
  bind_rows(., WAWA_hatchery_overshoot_ICH) %>% 
  bind_rows(., TUC_wild_overshoot) %>% 
  bind_rows(., TUC_hatchery_overshoot) %>% 
  dplyr::rename(overshoot_lower = `0.025`,
                overshoot_median = `0.5`,
                overshoot_upper = `0.975`) -> overshoot_probs_df


# sum overshoot at PRA at ICH for WAWA/YAK origins
overshoot_probs_df %>% 
  filter(origin %in% c("Walla Walla River", "Yakima River")) %>% 
  group_by(origin, rear) %>% 
  summarise_at(c("overshoot_lower", "overshoot_median", "overshoot_upper"), sum) %>% 
  mutate(state = "mainstem, ICH to LGR OR mainstem, PRA to RIS") -> double_overshoot_origins

overshoot_probs_df %>% 
  bind_rows(., double_overshoot_origins) -> overshoot_probs_df

# Join with the homing data
FF_homing %>% 
  left_join(., overshoot_probs_df, by = c("origin", "rear")) %>% 
  relocate(origin, rear) -> FF_overshoot_comp

trib_abbrev <- data.frame(
  origin = origin_param_map$natal_origin,
  abbrev = c("DES", "JDR", "15M",
             "UMA", "YAK", "WAWA", "WEN", "ENT",
             "OKA", "MET", "TUC", "ASO",
             "CLE", "SAL", "GRRO", "IMN")
)

FF_overshoot_comp %>% 
  left_join(., trib_abbrev, by = "origin") -> FF_overshoot_comp

# add DPS info
trib_DPS <- data.frame(origin = origin_param_map$natal_origin,
                       DPS = c("Middle Columbia", "Middle Columbia", "Middle Columbia",
                               "Middle Columbia", "Middle Columbia", "Middle Columbia", 
                               "Upper Columbia", "Upper Columbia",
                               "Upper Columbia", "Upper Columbia",
                               "Snake River", "Snake River",
                               "Snake River", "Snake River",
                               "Snake River", "Snake River"))

FF_overshoot_comp %>% 
  left_join(., trib_DPS, by = "origin") -> FF_overshoot_comp

# Leave NAs; instead, let's separate these out
# Change NA to -0.1
# FF_overshoot_comp %>% 
#   mutate(overshoot_median = ifelse(is.na(overshoot_median), -0.1, overshoot_median)) -> FF_overshoot_comp

# Drop ones that we don't want to plot
FF_overshoot_comp %>% 
  # first, keep only the combined overshoot probs
  filter(!(origin %in% c("Walla Walla River", "Yakima River") & state %in% c("mainstem, ICH to LGR", "mainstem, PRA to RIS"))) %>% 
  # drop Snake River origins that don't have detection efficiency corrections
  filter(!(origin %in% c("Salmon River", "Clearwater River", "Grande Ronde River")))-> FF_overshoot_comp_filter

# make some manual edits to labels


FF_overshoot_comp_filter %>% 
  mutate(label_position = ifelse(abbrev == "UMA" & rear == "hatchery", "topleft",
                                 ifelse(abbrev == "YAK" & rear == "wild" |
                                          abbrev == "DES" & rear == "wild" |
                                          abbrev == "ENT" & rear == "wild" |
                                          abbrev == "JDR" & rear == "wild" |
                                          abbrev == "UMA" & rear == "wild" |
                                          abbrev == "TUC" & rear == "wild" |
                                          abbrev == "WAWA" & rear == "hatchery", "topright",
                                        ifelse(abbrev == "TUC" & rear == "hatchery" |
                                                 abbrev == "WEN" & rear == "hatchery", "bottomleft",
                                               ifelse(abbrev == "WEN" & rear == "wild" |
                                                        abbrev == "15M" & rear == "wild" |
                                                        abbrev == "WAWA" & rear == "wild", "bottomright", NA))))) %>% 
  # change the non-overshooting label positions to be just left or right
  mutate(label_position = ifelse(abbrev == "MET" & rear == "wild" |
                                   abbrev == "OKA" & rear == "hatchery", "left",
                                 ifelse(abbrev == "IMN" & rear == "wild" |
                                          abbrev == "ASO" & rear == "wild" |
                                          abbrev == "IMN" & rear == "hatchery" |
                                          abbrev == "MET" & rear == "hatchery", "right", label_position)))-> FF_overshoot_comp_filter

FF_overshoot_comp_filter$DPS <- factor(FF_overshoot_comp_filter$DPS, levels = c("Middle Columbia",
                                                                                "Upper Columbia",
                                                                                "Snake River"))

FF_overshoot_comp_filter %>% 
  mutate(rear = ifelse(rear == "wild", "natural", rear)) -> FF_overshoot_comp_filter

# export this as a CSV
write.csv(FF_overshoot_comp_filter, here::here("figures", "paper_figures", "overshoot_probability_table.csv"))


#### Split dataset into two: those that can overshoot and those that can't ####

FF_overshoot_comp_filter %>% 
  filter(is.na(overshoot_median)) %>% 
  mutate(overshoot_median = 0) -> FF_overshoot_comp_filter_no_overshoot

# manually jitter these so that we can see them better
FF_overshoot_comp_filter_no_overshoot %>% 
  mutate(overshoot_median = ifelse(abbrev == "IMN" & rear == "natural" |
                                     abbrev == "ASO" & rear == "natural" |
                                     abbrev == "MET" & rear == "hatchery", -0.02, 0.02)) -> FF_overshoot_comp_filter_no_overshoot

FF_overshoot_comp_filter_no_overshoot$overshoot_median <- seq(-0.06, 0.06, length.out = nrow(FF_overshoot_comp_filter_no_overshoot))


# different approach: just make a new column
# FF_overshoot_comp_filter_no_overshoot %>% 
#   arrange(desc(homing_median)) -> FF_overshoot_comp_filter_no_overshoot
# FF_overshoot_comp_filter_no_overshoot$overshoot_median <- seq(-0.1, 0.1, length.out = 6)

FF_overshoot_comp_filter %>% 
  filter(!(is.na(overshoot_median))) -> FF_overshoot_comp_filter_yes_overshoot

# add a new column: rear-DPS, which we can use to add fills
FF_overshoot_comp_filter_no_overshoot %>% 
  mutate(rear_DPS = paste(rear, DPS)) -> FF_overshoot_comp_filter_no_overshoot
FF_overshoot_comp_filter_yes_overshoot %>% 
  mutate(rear_DPS = paste(rear, DPS)) -> FF_overshoot_comp_filter_yes_overshoot

#### Generate plot v2: two panels, one for overshooting and one for non-overshooting
library(ggrepel)

# create the plot
rear_colors <- c(hatchery = "#ff7f00", natural = "#33a02c")
DPS_colors <- c("Middle Columbia" = "#66c2a5", "Upper Columbia" = "#8da0cb", "Snake River" = "#fc8d62")
rear_DPS_fills <- c("natural Upper Columbia" =  "#8da0cb",  "hatchery Upper Columbia" = "white",
                    "natural Middle Columbia" = "#66c2a5",  "hatchery Middle Columbia" = "white",
                    "natural Snake River" = "#fc8d62", "hatchery Snake River" = "white")
rear_shapes <- c(21, 24)

FF_non_overshoot_comp_errorbar_plot <- ggplot(FF_overshoot_comp_filter_no_overshoot, 
                                              aes(x = overshoot_median, y = homing_median, shape = rear, 
                                                  color = DPS, fill = rear_DPS, label = abbrev))  +
  geom_errorbar(aes(ymin = homing_lower, ymax = homing_upper), width = 0, show.legend=FALSE) +
  geom_errorbar(aes(xmin = overshoot_lower, xmax = overshoot_upper), width = 0, show.legend=FALSE) +
  geom_point(aes(size = rear)) +
  # geom_text_repel(max.overlaps = 100, show.legend=FALSE) +
  geom_text(data = subset(FF_overshoot_comp_filter_no_overshoot, label_position == "left"), aes(label = abbrev, x = overshoot_median-0.12, y = homing_median), 
            hjust = 0.5, show.legend=FALSE) +
  geom_text(data = subset(FF_overshoot_comp_filter_no_overshoot, label_position == "right"), aes(label = abbrev, x = overshoot_median+0.12, y = homing_median), 
            hjust = 0.5, show.legend=FALSE) +
  scale_color_manual(values = DPS_colors) +
  scale_fill_manual(values = rear_DPS_fills, guide = "none") +
  scale_size_manual(values = c("hatchery" = 3, "natural" = 2.5), guide = "none") +
  scale_shape_manual(values = rear_shapes) +
  scale_y_continuous(lim = c(0, 1.0)) +
  scale_x_continuous(lim = c(-0.25, 0.25),
                     breaks = c(0),
                     labels = c("NA")) +
  xlab("Probability of Overshoot") +
  ylab("Probability of Homing") +
  theme(legend.position = "none",
        # legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_text(size = 10),
        # axis.text.x = element_text(color = "white"),
        axis.ticks.x = element_line(color = "white"),
        axis.title.x = element_text(color = "white"),
        legend.key = element_blank(),
        plot.margin = unit(c(0.5, 0, 0.5, 0.5),"cm"))

FF_overshoot_comp_errorbar_plot <- ggplot(FF_overshoot_comp_filter_yes_overshoot, 
                                          aes(x = overshoot_median, y = homing_median, 
                                              shape = rear, color = DPS, fill = rear_DPS, 
                                              label = abbrev))  +
  geom_errorbar(aes(ymin = homing_lower, ymax = homing_upper), width = 0, show.legend=FALSE) +
  geom_errorbar(aes(xmin = overshoot_lower, xmax = overshoot_upper), width = 0, show.legend=FALSE) +
  geom_point(aes(size = rear)) +
  # geom_text_repel(max.overlaps = 100, show.legend=FALSE) +
  geom_text(data = subset(FF_overshoot_comp_filter_yes_overshoot, label_position == "topleft"), aes(label = abbrev, x = overshoot_median-0.01, y = homing_median+0.025), 
            hjust = 1, show.legend=FALSE) +
  geom_text(data = subset(FF_overshoot_comp_filter_yes_overshoot, label_position == "topright"), aes(label = abbrev, x = overshoot_median+0.01, y = homing_median+0.025), 
            hjust = 0, show.legend=FALSE) +
  geom_text(data = subset(FF_overshoot_comp_filter_yes_overshoot, label_position == "bottomleft"), aes(label = abbrev, x = overshoot_median-0.01, y = homing_median-0.025), 
            hjust = 1, show.legend=FALSE) +
  geom_text(data = subset(FF_overshoot_comp_filter_yes_overshoot, label_position == "bottomright"), aes(label = abbrev, x = overshoot_median+0.01, y = homing_median-0.025), 
            hjust = 0, show.legend=FALSE) +
  scale_color_manual(values = DPS_colors) +
  scale_shape_manual(values = rear_shapes) +
  scale_size_manual(values = c("hatchery" = 3, "natural" = 2.5), guide = "none") +
  scale_fill_manual(values = rear_DPS_fills, guide = "none") +
  scale_y_continuous(lim = c(0, 1.0)) +
  scale_x_continuous(lim = c(0, 1.0), breaks = seq(0,1,0.2),
                     labels = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  xlab("Probability of Overshoot") +
  ylab("Probability of Homing") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.83, 0.75),
        # legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length.y = unit(0, "cm"),
        axis.line.y.left = element_line(size = 2, colour = "gray80"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.key = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0),"cm")) +
  guides(shape=guide_legend(title="Rearing Type", override.aes = list(fill = c("white", "black"),
                                                                      size = c(3, 2.5))),
         color = guide_legend(override.aes = list(size = 3)))

fig4 <- ggarrange(FF_non_overshoot_comp_errorbar_plot, FF_overshoot_comp_errorbar_plot, ncol = 2, nrow = 1,
                  widths = c(1.1,4)) + bgcolor("white") + border(color = NA)


ggsave(here::here("figures", "paper_figures", "fig3_homing_v_overshoot_plot.png"), fig4, height = 6, width = 8)
