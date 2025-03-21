# 07-03_final_fates_simulation

# This script takes the output from the stan model runs 
# generates estimates of final fate distributions

# First, need to load in all of the model runs and all of the packages.
source("R/07_model_analysis/alternative_spill_treatment/07-01_load_stan_models_marchspill.R")

#### Final fates functions ####

# This function takes a model fit object, a number of simulated fish,
# the identities of those fish (origin, rear), the conditions throughout the
# basin (spill, temperature, year)
# and then simulates the final fates of those fish
# By default, the conditions (spill and temperature) will be taken from the data
# itself, but note that new values of these conditions could be simulated in order

# Note that another option here would be to use average (median) conditions experienced
# in each state - and that would remove the variability associated with taking random 
# conditions. But note that you'll have to decide what the data is that you're taking
# the median of - is it the whole DPS? Just the population? Currently, we have it set up
# to run by DPS, not population
# Another option would be to make the number of fish per simulation smaller and
# run it more times - that would give you a wider spread of covariate values
# and therefore reduce the chance that a handful of outlier covariate values
# lead to crazy distributions


# Final states simulation
# This function will have to be re-written for each DPS, because each DPS has different params
# currently we are not including a random effect of year, but we could easily
# amend the code to simulate final fates in a particular year
# origin1/origin2/origin3 are indicator variables, so set the origin you want to be 1
final_fates_simulation_UCW <- function(nsim,
                                       start_state = 2, states_dates,
                                       origin1 = 0, origin2 = 0, origin3 = 0,
                                       temp_data, spillwindow_data, winterspill_data,
                                       condition_jitter = FALSE){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(UCW_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      winterspill[i] <- winterspill_data[sample_year[i],i] 
    } else {
      # approach 2: Take the median conditions for each state across all years
      winterspill[i] <- median(winterspill_data[,i])
    }
    
    
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
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(UCW_states_dates, state == i  & direction == "upstream")$season)/length(subset(UCW_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(UCW_states_dates, state == i  & direction == "downstream")$season)/length(subset(UCW_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      # temp_upstream[i] <- temp_data[sample_date[i],i]
      # temp_downstream[i] <- temp_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(UCW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(UCW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(UCW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(UCW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
      }
      
    }
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      # spillwindow_upstream[i] <- spillwindow_data[sample_date[i],i]
      # spillwindow_downstream[i] <- spillwindow_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(UCW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(UCW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(UCW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(UCW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
      }
      
    }
  }
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCW[i,possible_movements[j],iter] +
                                                                    btemp1_array_UCW[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_UCW[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_UCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_UCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_UCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_UCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    borigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_UCW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_UCW[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCW[i,possible_movements,iter] +
                    btemp1_array_UCW[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_UCW[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_UCW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_UCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_UCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_UCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    borigin1_array_UCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCW[i,possible_movements,iter]*origin3))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCW[i,possible_movements[j],iter] +
                                                                    btemp0_array_UCW[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_UCW[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_UCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_UCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_UCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_UCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    borigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_UCW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_UCW[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCW[i,possible_movements,iter] +
                    btemp0_array_UCW[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_UCW[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_UCW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_UCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_UCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_UCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    borigin1_array_UCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCW[i,possible_movements,iter]*origin3))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCW[i,possible_movements[j],iter] +
                                                                      btemp1_array_UCW[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_UCW[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_UCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_UCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_UCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_UCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      borigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_UCW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_UCW[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCW[i,possible_movements,iter] +
                    btemp1_array_UCW[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_UCW[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_UCW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_UCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_UCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_UCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    borigin1_array_UCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCW[i,possible_movements,iter]*origin3))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCW[i,possible_movements[j],iter] +
                                                                      btemp0_array_UCW[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_UCW[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_UCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_UCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_UCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_UCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      borigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_UCW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_UCW[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCW[i,possible_movements,iter] +
                    btemp0_array_UCW[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_UCW[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_UCW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_UCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_UCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_UCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    borigin1_array_UCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCW[i,possible_movements,iter]*origin3))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_UCH <- function(nsim,
                                       start_state = 2, states_dates,
                                       origin1 = 0, origin2 = 0, origin3 = 0,
                                       temp_data, spillwindow_data, winterspill_data,
                                       condition_jitter = FALSE){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(UCH_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      winterspill[i] <- winterspill_data[sample_year[i],i] 
    } else {
      # approach 2: Take the median conditions for each state across all years
      winterspill[i] <- median(winterspill_data[,i])
    }
    
    
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
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(UCH_states_dates, state == i  & direction == "upstream")$season)/length(subset(UCH_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(UCH_states_dates, state == i  & direction == "downstream")$season)/length(subset(UCH_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      # temp_upstream[i] <- temp_data[sample_date[i],i]
      # temp_downstream[i] <- temp_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(UCH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(UCH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(UCH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(UCH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
      }
      
    }
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      # spillwindow_upstream[i] <- spillwindow_data[sample_date[i],i]
      # spillwindow_downstream[i] <- spillwindow_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(UCH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(UCH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(UCH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(UCH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
      }
      
    }
  }
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCH[i,possible_movements[j],iter] +
                                                                    btemp1_array_UCH[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_UCH[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_UCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_UCH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_UCH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_UCH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    borigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_UCH[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_UCH[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCH[i,possible_movements,iter] +
                    btemp1_array_UCH[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_UCH[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_UCH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_UCH[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_UCH[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_UCH[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    borigin1_array_UCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCH[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCH[i,possible_movements,iter]*origin3))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCH[i,possible_movements[j],iter] +
                                                                    btemp0_array_UCH[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_UCH[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_UCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_UCH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_UCH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_UCH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    borigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_UCH[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_UCH[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCH[i,possible_movements,iter] +
                    btemp0_array_UCH[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_UCH[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_UCH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_UCH[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_UCH[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_UCH[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    borigin1_array_UCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCH[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCH[i,possible_movements,iter]*origin3))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCH[i,possible_movements[j],iter] +
                                                                      btemp1_array_UCH[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_UCH[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_UCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_UCH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_UCH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_UCH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      borigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_UCH[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_UCH[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCH[i,possible_movements,iter] +
                    btemp1_array_UCH[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_UCH[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_UCH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_UCH[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_UCH[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_UCH[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    borigin1_array_UCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCH[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCH[i,possible_movements,iter]*origin3))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCH[i,possible_movements[j],iter] +
                                                                      btemp0_array_UCH[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_UCH[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_UCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_UCH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_UCH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_UCH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      borigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_UCH[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_UCH[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCH[i,possible_movements,iter] +
                    btemp0_array_UCH[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_UCH[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_UCH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_UCH[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_UCH[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_UCH[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    borigin1_array_UCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCH[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCH[i,possible_movements,iter]*origin3))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_MCW <- function(nsim,
                                       start_state = 2, states_dates,
                                       origin1 = 0, origin2 = 0, origin3 = 0,
                                       origin4 = 0, origin5 = 0, origin6 = 0,
                                       temp_data, spillwindow_data, winterspill_data,
                                       condition_jitter = FALSE){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(MCW_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      winterspill[i] <- winterspill_data[sample_year[i],i] 
    } else {
      # approach 2: Take the median conditions for each state across all years
      winterspill[i] <- median(winterspill_data[,i])
    }
    
    
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
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(MCW_states_dates, state == i  & direction == "upstream")$season)/length(subset(MCW_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(MCW_states_dates, state == i  & direction == "downstream")$season)/length(subset(MCW_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      # temp_upstream[i] <- temp_data[sample_date[i],i]
      # temp_downstream[i] <- temp_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(MCW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(MCW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(MCW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(MCW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
      }
      
    }
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      # spillwindow_upstream[i] <- spillwindow_data[sample_date[i],i]
      # spillwindow_downstream[i] <- spillwindow_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(MCW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(MCW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(MCW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(MCW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
      }
      
    }
  }
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCW[i,possible_movements[j],iter] +
                                                                    btemp1_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_MCW[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_MCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    btemp1xorigin6_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin6 +
                                                                    borigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_MCW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_MCW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_MCW[i,possible_movements,iter] +
                    btemp1_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_MCW[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_MCW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    btemp1xorigin6_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin6 +
                    borigin1_array_MCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_MCW[i,possible_movements,iter]*origin3 +
                    borigin4_array_MCW[i,possible_movements,iter]*origin4 +
                    borigin5_array_MCW[i,possible_movements,iter]*origin5 +
                    borigin6_array_MCW[i,possible_movements,iter]*origin6))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCW[i,possible_movements[j],iter] +
                                                                    btemp0_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_MCW[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_MCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    btemp0xorigin6_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin6 +
                                                                    borigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_MCW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_MCW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_MCW[i,possible_movements,iter] +
                    btemp0_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_MCW[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_MCW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    btemp0xorigin6_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin6 +
                    borigin1_array_MCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_MCW[i,possible_movements,iter]*origin3 +
                    borigin4_array_MCW[i,possible_movements,iter]*origin4 +
                    borigin5_array_MCW[i,possible_movements,iter]*origin5 +
                    borigin6_array_MCW[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCW[i,possible_movements[j],iter] +
                                                                      btemp1_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_MCW[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_MCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      btemp1xorigin6_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin6 +
                                                                      borigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_MCW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_MCW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_MCW[i,possible_movements,iter] +
                    btemp1_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_MCW[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_MCW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    btemp1xorigin6_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin6 +
                    borigin1_array_MCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_MCW[i,possible_movements,iter]*origin3 +
                    borigin4_array_MCW[i,possible_movements,iter]*origin4 +
                    borigin5_array_MCW[i,possible_movements,iter]*origin5 +
                    borigin6_array_MCW[i,possible_movements,iter]*origin6))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCW[i,possible_movements[j],iter] +
                                                                      btemp0_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_MCW[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_MCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      btemp0xorigin6_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin6 +
                                                                      borigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_MCW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_MCW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_MCW[i,possible_movements,iter] +
                    btemp0_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_MCW[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_MCW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    btemp0xorigin6_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin6 +
                    borigin1_array_MCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_MCW[i,possible_movements,iter]*origin3 +
                    borigin4_array_MCW[i,possible_movements,iter]*origin4 +
                    borigin5_array_MCW[i,possible_movements,iter]*origin5 +
                    borigin6_array_MCW[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_MCH <- function(nsim,
                                       start_state = 2, states_dates,
                                       origin1 = 0, origin2 = 0,
                                       temp_data, spillwindow_data, winterspill_data,
                                       condition_jitter = FALSE){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(MCH_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      winterspill[i] <- winterspill_data[sample_year[i],i] 
    } else {
      # approach 2: Take the median conditions for each state across all years
      winterspill[i] <- median(winterspill_data[,i])
    }
    
    
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
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(MCH_states_dates, state == i  & direction == "upstream")$season)/length(subset(MCH_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(MCH_states_dates, state == i  & direction == "downstream")$season)/length(subset(MCH_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      # temp_upstream[i] <- temp_data[sample_date[i],i]
      # temp_downstream[i] <- temp_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(MCH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(MCH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(MCH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(MCH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
      }
      
    }
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      # spillwindow_upstream[i] <- spillwindow_data[sample_date[i],i]
      # spillwindow_downstream[i] <- spillwindow_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(MCH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(MCH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(MCH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(MCH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
      }
      
    }
  }
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCH[i,possible_movements[j],iter] +
                                                                    btemp1_array_MCH[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_MCH[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_MCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_MCH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_MCH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    borigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_MCH[i,possible_movements[j],iter]*origin2)/
          sum(exp(b0_array_MCH[i,possible_movements,iter] +
                    btemp1_array_MCH[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_MCH[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_MCH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_MCH[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_MCH[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    borigin1_array_MCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCH[i,possible_movements,iter]*origin2))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCH[i,possible_movements[j],iter] +
                                                                    btemp0_array_MCH[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_MCH[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_MCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_MCH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_MCH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    borigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_MCH[i,possible_movements[j],iter]*origin2)/
          sum(exp(b0_array_MCH[i,possible_movements,iter] +
                    btemp0_array_MCH[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_MCH[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_MCH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_MCH[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_MCH[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    borigin1_array_MCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCH[i,possible_movements,iter]*origin2))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCH[i,possible_movements[j],iter] +
                                                                      btemp1_array_MCH[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_MCH[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_MCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_MCH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_MCH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      borigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_MCH[i,possible_movements[j],iter]*origin2)/
          sum(exp(b0_array_MCH[i,possible_movements,iter] +
                    btemp1_array_MCH[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_MCH[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_MCH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_MCH[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_MCH[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    borigin1_array_MCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCH[i,possible_movements,iter]*origin2))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCH[i,possible_movements[j],iter] +
                                                                      btemp0_array_MCH[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_MCH[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_MCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_MCH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_MCH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      borigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_MCH[i,possible_movements[j],iter]*origin2)/
          sum(exp(b0_array_MCH[i,possible_movements,iter] +
                    btemp0_array_MCH[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_MCH[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_MCH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_MCH[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_MCH[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    borigin1_array_MCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCH[i,possible_movements,iter]*origin2))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_SRW <- function(nsim,
                                       start_state = 2, states_dates,
                                       origin1 = 0, origin2 = 0, origin3 = 0,
                                       origin4 = 0, origin5 = 0, origin6 = 0,
                                       temp_data, spillwindow_data, winterspill_data,
                                       condition_jitter = FALSE){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(SRW_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      winterspill[i] <- winterspill_data[sample_year[i],i] 
    } else {
      # approach 2: Take the median conditions for each state across all years
      winterspill[i] <- median(winterspill_data[,i])
    }
    
    
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
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRW_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRW_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRW_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRW_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      # temp_upstream[i] <- temp_data[sample_date[i],i]
      # temp_downstream[i] <- temp_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
      }
      
    }
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      # spillwindow_upstream[i] <- spillwindow_data[sample_date[i],i]
      # spillwindow_downstream[i] <- spillwindow_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
      }
      
    }
  }
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRW[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    btemp1xorigin6_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin6 +
                                                                    borigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW[i,possible_movements,iter] +
                    btemp1_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRW[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin6 +
                    borigin1_array_SRW[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW[i,possible_movements,iter]*origin6))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRW[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    btemp0xorigin6_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin6 +
                                                                    borigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW[i,possible_movements,iter] +
                    btemp0_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRW[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin6 +
                    borigin1_array_SRW[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRW[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      btemp1xorigin6_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin6 +
                                                                      borigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW[i,possible_movements,iter] +
                    btemp1_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRW[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin6 +
                    borigin1_array_SRW[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW[i,possible_movements,iter]*origin6))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRW[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      btemp0xorigin6_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin6 +
                                                                      borigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW[i,possible_movements,iter] +
                    btemp0_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRW[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin6 +
                    borigin1_array_SRW[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}
final_fates_simulation_SRH <- function(nsim,
                                       start_state = 2, states_dates,
                                       origin1 = 0, origin2 = 0, origin3 = 0,
                                       origin4 = 0, origin5 = 0, 
                                       temp_data, spillwindow_data, winterspill_data,
                                       condition_jitter = FALSE){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(SRH_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      winterspill[i] <- winterspill_data[sample_year[i],i] 
    } else {
      # approach 2: Take the median conditions for each state across all years
      winterspill[i] <- median(winterspill_data[,i])
    }
    
    
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
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRH_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRH_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRH_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRH_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      # temp_upstream[i] <- temp_data[sample_date[i],i]
      # temp_downstream[i] <- temp_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
      }
      
    }
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      # spillwindow_upstream[i] <- spillwindow_data[sample_date[i],i]
      # spillwindow_downstream[i] <- spillwindow_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
      }
      
    }
  }
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRH[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    borigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH[i,possible_movements,iter] +
                    btemp1_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRH[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    borigin1_array_SRH[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH[i,possible_movements,iter]*origin5))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRH[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    borigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH[i,possible_movements,iter] +
                    btemp0_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRH[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    borigin1_array_SRH[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH[i,possible_movements,iter]*origin5))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))    
    
    for (j in 1:length(possible_movements)){
      
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRH[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      borigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH[i,possible_movements,iter] +
                    btemp1_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRH[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    borigin1_array_SRH[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH[i,possible_movements,iter]*origin5))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRH[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      borigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH[i,possible_movements,iter] +
                    btemp0_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRH[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    borigin1_array_SRH[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH[i,possible_movements,iter]*origin5))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}





#### Run final fates simulation - H v W comparisons ####


# In order to simulate covariate values, we are going to determine the dates
# where the fish were in each state, and then sample from those dates to 
# get a representative value for temperature/spill when fish are in those states
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


UCW_states_dates <- get_states_dates_direction(envir = UCW_envir)
UCH_states_dates <- get_states_dates_direction(envir = UCH_envir)
MCW_states_dates <- get_states_dates_direction(envir = MCW_envir)
MCH_states_dates <- get_states_dates_direction(envir = MCH_envir)
SRW_states_dates <- get_states_dates_direction(envir = SRW_envir)
SRH_states_dates <- get_states_dates_direction(envir = SRH_envir)

# Get states/dates by origin as well
# Also, track the seasonality of movement
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

compare_final_fate_rear_type_UC <- function(niter, nsim, condition_jitter,
                                            origin_select){
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- c(0,0,0)
    
    # index to the right origin param to turn it on
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_UCW(nsim = nsim,
                                             start_state = 2, states_dates = get_origin_states_dates(envir = UCW_envir, origin_select = origin_select, rear = "wild"),
                                             origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                             temp_data = UCW_envir$data$temperature_data, spillwindow_data = UCW_envir$data$spill_window_data, 
                                             winterspill_data = UCW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    ff_wild_quantiles -> ff_rear_quantiles
    
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    
    
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    hatchery_origin_params <- c(0,0,0)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_UCH(nsim = nsim,
                                                 start_state = 2, states_dates = get_origin_states_dates(envir = UCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                 origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], origin3 = hatchery_origin_params[3],
                                                 temp_data = UCH_envir$data$temperature_data, spillwindow_data = UCH_envir$data$spill_window_data, 
                                                 winterspill_data = UCH_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_hatchery_quantiles -> ff_rear_quantiles
  } 
  # else run both
  else {
    
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- c(0,0,0)
    hatchery_origin_params <- c(0,0,0)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_UCW(nsim = nsim,
                                             start_state = 2, states_dates = get_origin_states_dates(envir = UCW_envir, origin_select = origin_select, rear = "wild"),
                                             origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                             temp_data = UCW_envir$data$temperature_data, spillwindow_data = UCW_envir$data$spill_window_data, 
                                             winterspill_data = UCW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_UCH(nsim = nsim,
                                                 start_state = 2, states_dates = get_origin_states_dates(envir = UCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                 origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], origin3 = hatchery_origin_params[3],
                                                 temp_data = UCH_envir$data$temperature_data, spillwindow_data = UCH_envir$data$spill_window_data, 
                                                 winterspill_data = UCH_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_wild_quantiles %>% 
      bind_rows(ff_hatchery_quantiles) -> ff_rear_quantiles
  }
  
  
  return(ff_rear_quantiles)
  
}

compare_final_fate_rear_type_MC <- function(niter, nsim, condition_jitter,
                                            origin_select){
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    
    # index to the right origin param to turn it on
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_MCW(nsim = nsim,
                                             start_state = 2, states_dates = get_origin_states_dates(envir = MCW_envir, origin_select = origin_select, rear = "wild"),
                                             origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                             origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                             temp_data = MCW_envir$data$temperature_data, spillwindow_data = MCW_envir$data$spill_window_data, 
                                             winterspill_data = MCW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    ff_wild_quantiles -> ff_rear_quantiles
    
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    
    
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    hatchery_origin_params <- c(0,0)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_MCH(nsim = nsim,
                                                 start_state = 2, states_dates = get_origin_states_dates(envir = MCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                 origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                 temp_data = MCH_envir$data$temperature_data, spillwindow_data = MCH_envir$data$spill_window_data, 
                                                 winterspill_data = MCH_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_hatchery_quantiles -> ff_rear_quantiles
  } 
  # else run both
  else {
    
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    hatchery_origin_params <- c(0,0)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_MCW(nsim = nsim,
                                             start_state = 2, states_dates = get_origin_states_dates(envir = MCW_envir, origin_select = origin_select, rear = "wild"),
                                             origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                             origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                             temp_data = MCW_envir$data$temperature_data, spillwindow_data = MCW_envir$data$spill_window_data, 
                                             winterspill_data = MCW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_MCH(nsim = nsim,
                                                 start_state = 2, states_dates = get_origin_states_dates(envir = MCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                 origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2],
                                                 temp_data = MCH_envir$data$temperature_data, spillwindow_data = MCH_envir$data$spill_window_data, 
                                                 winterspill_data = MCH_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_wild_quantiles %>% 
      bind_rows(ff_hatchery_quantiles) -> ff_rear_quantiles
  }
  
  
  return(ff_rear_quantiles)
  
}

compare_final_fate_rear_type_SR <- function(niter, nsim, condition_jitter,
                                            origin_select){
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    
    # index to the right origin param to turn it on
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_SRW(nsim = nsim,
                                             start_state = 2, states_dates = get_origin_states_dates(envir = SRW_envir, origin_select = origin_select, rear = "wild"),
                                             origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                             origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                             temp_data = SRW_envir$data$temperature_data, spillwindow_data = SRW_envir$data$spill_window_data, 
                                             winterspill_data = SRW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    ff_wild_quantiles -> ff_rear_quantiles
    
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    
    
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    hatchery_origin_params <- rep(0,5)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_SRH(nsim = nsim,
                                                 start_state = 2, states_dates = get_origin_states_dates(envir = SRH_envir, origin_select = origin_select, rear = "hatchery"),
                                                 origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                 origin3 = hatchery_origin_params[3], origin4 = hatchery_origin_params[4], 
                                                 origin5 = hatchery_origin_params[5],
                                                 temp_data = SRH_envir$data$temperature_data, spillwindow_data = SRH_envir$data$spill_window_data, 
                                                 winterspill_data = SRH_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_hatchery_quantiles -> ff_rear_quantiles
  } 
  # else run both
  else {
    
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    hatchery_origin_params <- rep(0,5)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_SRW(nsim = nsim,
                                             start_state = 2, states_dates = get_origin_states_dates(envir = SRW_envir, origin_select = origin_select, rear = "wild"),
                                             origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                             origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                             temp_data = SRW_envir$data$temperature_data, spillwindow_data = SRW_envir$data$spill_window_data, 
                                             winterspill_data = SRW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_SRH(nsim = nsim,
                                                 start_state = 2, states_dates = get_origin_states_dates(envir = SRH_envir, origin_select = origin_select, rear = "hatchery"),
                                                 origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                 origin3 = hatchery_origin_params[3], origin4 = hatchery_origin_params[4], 
                                                 origin5 = hatchery_origin_params[5],
                                                 temp_data = SRH_envir$data$temperature_data, spillwindow_data = SRH_envir$data$spill_window_data, 
                                                 winterspill_data = SRH_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_wild_quantiles %>% 
      bind_rows(ff_hatchery_quantiles) -> ff_rear_quantiles
  }
  
  
  return(ff_rear_quantiles)
  
}


# ORDER THE STATES FOR PLOTTING
states_order_for_plot <- gsub(" Mouth| Upstream", "", model_states)
states_order_for_plot <- states_order_for_plot[!(duplicated(states_order_for_plot))]
# Make a couple of changes to make them be in the order from most downstream to most upstream
states_order_for_plot[10] <- "Fifteenmile Creek"
states_order_for_plot[11] <- "Deschutes River"
states_order_for_plot[12] <- "John Day River"
states_order_for_plot[14] <- "Walla Walla River"
states_order_for_plot[15] <- "Yakima River"
states_order_for_plot[18] <- "Methow River"
states_order_for_plot[19] <- "Okanogan River"

states_order_for_plot[15:28] <- states_order_for_plot[14:27]
states_order_for_plot[14] <- "BON to MCN other tributaries"
states_order_for_plot[22:28] <- states_order_for_plot[21:27]
states_order_for_plot[21] <- "Upstream WEL other tributaries"
states_order_for_plot[28] <- "loss"

states_order_for_plot[23] <- "Clearwater River"
states_order_for_plot[24] <- "Asotin Creek"
states_order_for_plot[25] <- "Grande Ronde River"
states_order_for_plot[26] <- "Salmon River"


plot_final_fate_rear_type <- function(ff_comp, natal_origin){
  rear_colors <- c(hatchery = "#ff7f00", natural = "#33a02c")
  
  ff_comp %>% 
    mutate(rear_type = ifelse(rear_type == "wild", "natural", rear_type)) -> ff_comp
  
  ff_comp$state <- fct_rev(factor(ff_comp$state, levels = states_order_for_plot))
  
  ff_comp_plot <- ggplot(ff_comp, aes(x = state, y = `0.5`, ymin = `0.025`, ymax = `0.975`, color = rear_type)) +
    geom_point(size = 3.5, shape = 18, position=position_dodge(width=0.5)) +
    geom_linerange(position=position_dodge(width=0.5)) +
    coord_flip() +
    ylab("Final Distribution Probability") +
    xlab("Model State") +
    # ggtitle(" ") +
    # Create a scale common to all
    scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_color_manual(values = rear_colors) +
    theme(panel.grid.major = element_line(color = "gray90"),
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
          plot.title = element_text(size = 12),
          # axis.text.y = element_text(color = rev(state_significance_colors)),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12)) +
    guides(color=guide_legend(title="Rearing Type")) +
    ggtitle(natal_origin)
  
  return(ff_comp_plot)
  
}


ff_iter <- 1000
ff_nsim <- 1000

# If you just want to regenerate plots:
# load(file = here::here("figures", "alternative_spill_treatment", "final_fates", "FF_comp_data.rda"))
# wen_ff_comp_median <- FF_comp_data$wen_ff_comp_median
# ent_ff_comp_median <- FF_comp_data$ent_ff_comp_median
# oka_ff_comp_median <- FF_comp_data$oka_ff_comp_median
# met_ff_comp_median <- FF_comp_data$met_ff_comp_median
# des_ff_comp_median <- FF_comp_data$des_ff_comp_median
# jdr_ff_comp_median <- FF_comp_data$jdr_ff_comp_median
# fif_ff_comp_median <- FF_comp_data$fif_ff_comp_median
# uma_ff_comp_median <- FF_comp_data$uma_ff_comp_median
# yak_ff_comp_median <- FF_comp_data$yak_ff_comp_median
# wawa_ff_comp_median <- FF_comp_data$wawa_ff_comp_median
# aso_ff_comp_median <- FF_comp_data$aso_ff_comp_median
# cle_ff_comp_median <- FF_comp_data$cle_ff_comp_median
# sal_ff_comp_median <- FF_comp_data$sal_ff_comp_median
# gr_ff_comp_median <- FF_comp_data$gr_ff_comp_median
# imn_ff_comp_median <- FF_comp_data$imn_ff_comp_median
# tuc_ff_comp_median <- FF_comp_data$tuc_ff_comp_median

## Upper Columbia

# Wenatchee comparison
# wen_ff_comp_jitter <- compare_final_fate_rear_type_UC(niter = ff_iter, nsim = ff_nsim, origin_select = "Wenatchee River", condition_jitter = TRUE)
# wen_ff_comp_jitter_plot <- plot_final_fate_rear_type(wen_ff_comp, natal_origin = "Wenatchee River")
# ggsave(here::here("figures", "alternative_spill_treatment", "final_fates", "wen_ff_comp_jitter_plot.png"), wen_ff_comp_jitter_plot, height = 8, width = 8)
wen_ff_comp_median <- compare_final_fate_rear_type_UC(niter = ff_iter, nsim = ff_nsim, origin_select = "Wenatchee River", condition_jitter = FALSE)
save(wen_ff_comp_median, file = here::here("figures", "alternative_spill_treatment", "final_fates", "simulation_runs", "WEN_comp_data.rda"))
wen_ff_comp_median_plot <- plot_final_fate_rear_type(wen_ff_comp_median, natal_origin = "Wenatchee River")
ggsave(here::here("figures", "alternative_spill_treatment", "final_fates", "wen_ff_comp_median_plot.png"), wen_ff_comp_median_plot, height = 8, width = 8)

# Entiat comparison
ent_ff_comp_median <- compare_final_fate_rear_type_UC(niter = ff_iter, nsim = ff_nsim, origin_select = "Entiat River", condition_jitter = FALSE)
save(ent_ff_comp_median, file = here::here("figures", "alternative_spill_treatment", "final_fates", "simulation_runs", "ENT_comp_data.rda"))
ent_ff_comp_median_plot <- plot_final_fate_rear_type(ent_ff_comp_median, natal_origin = "Entiat River")
ggsave(here::here("figures", "alternative_spill_treatment", "final_fates", "ent_ff_comp_median_plot.png"), ent_ff_comp_median_plot, height = 8, width = 8)

# Okanogan comparison
oka_ff_comp_median <- compare_final_fate_rear_type_UC(niter = ff_iter, nsim = ff_nsim, origin_select = "Okanogan River", condition_jitter = FALSE)
save(oka_ff_comp_median, file = here::here("figures", "alternative_spill_treatment", "final_fates", "simulation_runs", "OKA_comp_data.rda"))
oka_ff_comp_median_plot <- plot_final_fate_rear_type(oka_ff_comp_median, natal_origin = "Okanogan River")
ggsave(here::here("figures", "alternative_spill_treatment", "final_fates", "oka_ff_comp_median_plot.png"), oka_ff_comp_median_plot, height = 8, width = 8)

# Methow comparison
# this used to crash with jitter - let's see if it works now
met_ff_comp_median <- compare_final_fate_rear_type_UC(niter = ff_iter, nsim = ff_nsim, origin_select = "Methow River", condition_jitter = FALSE)
save(met_ff_comp_median, file = here::here("figures", "alternative_spill_treatment", "final_fates", "simulation_runs", "MET_comp_data.rda"))
met_ff_comp_median_plot <- plot_final_fate_rear_type(met_ff_comp_median, natal_origin = "Methow River")
ggsave(here::here("figures", "alternative_spill_treatment", "final_fates", "met_ff_comp_median_plot.png"), met_ff_comp_median_plot, height = 8, width = 8)



# Middle Columbia
# Deschutes comparison
des_ff_comp_median <- compare_final_fate_rear_type_MC(niter = ff_iter, nsim = ff_nsim, origin_select = "Deschutes River", condition_jitter = FALSE)
save(des_ff_comp_median, file = here::here("figures", "alternative_spill_treatment", "final_fates", "simulation_runs", "DES_comp_data.rda"))
des_ff_comp_median_plot <- plot_final_fate_rear_type(des_ff_comp_median, natal_origin = "Deschutes River")
ggsave(here::here("figures", "alternative_spill_treatment", "final_fates", "des_ff_comp_median_plot.png"), des_ff_comp_median_plot, height = 8, width = 8)

# John Day comparison
jdr_ff_comp_median <- compare_final_fate_rear_type_MC(niter = ff_iter, nsim = ff_nsim, origin_select = "John Day River", condition_jitter = FALSE)
save(jdr_ff_comp_median, file = here::here("figures", "alternative_spill_treatment", "final_fates", "simulation_runs", "JDR_comp_data.rda"))
jdr_ff_comp_median_plot <- plot_final_fate_rear_type(jdr_ff_comp_median, natal_origin = "John Day River")
ggsave(here::here("figures", "alternative_spill_treatment", "final_fates", "jdr_ff_comp_median_plot.png"), jdr_ff_comp_median_plot, height = 8, width = 8)

# Fifteenmile Creek comparison
fif_ff_comp_median <- compare_final_fate_rear_type_MC(niter = ff_iter, nsim = ff_nsim, origin_select = "Fifteenmile Creek", condition_jitter = FALSE)
save(fif_ff_comp_median, file = here::here("figures", "alternative_spill_treatment", "final_fates", "simulation_runs", "FIF_comp_data.rda"))
fif_ff_comp_median_plot <- plot_final_fate_rear_type(fif_ff_comp_median, natal_origin = "Fifteenmile Creek")
ggsave(here::here("figures", "alternative_spill_treatment", "final_fates", "fif_ff_comp_median_plot.png"), fif_ff_comp_median_plot, height = 8, width = 8)

# Umatilla comparison
uma_ff_comp_median <- compare_final_fate_rear_type_MC(niter = ff_iter, nsim = ff_nsim, origin_select = "Umatilla River", condition_jitter = FALSE)
save(uma_ff_comp_median, file = here::here("figures", "alternative_spill_treatment", "final_fates", "simulation_runs", "UMA_comp_data.rda"))
uma_ff_comp_median_plot <- plot_final_fate_rear_type(uma_ff_comp_median, natal_origin = "Umatilla River")
ggsave(here::here("figures", "alternative_spill_treatment", "final_fates", "uma_ff_comp_median_plot.png"), uma_ff_comp_median_plot, height = 8, width = 8)

# Yakima comparison
yak_ff_comp_median <- compare_final_fate_rear_type_MC(niter = ff_iter, nsim = ff_nsim, origin_select = "Yakima River", condition_jitter = FALSE)
save(yak_ff_comp_median, file = here::here("figures", "alternative_spill_treatment", "final_fates", "simulation_runs", "YAK_comp_data.rda"))
yak_ff_comp_median_plot <- plot_final_fate_rear_type(yak_ff_comp_median, natal_origin = "Yakima River")
ggsave(here::here("figures", "alternative_spill_treatment", "final_fates", "yak_ff_comp_median_plot.png"), yak_ff_comp_median_plot, height = 8, width = 8)

# Walla Walla comparison
wawa_ff_comp_median <- compare_final_fate_rear_type_MC(niter = ff_iter, nsim = ff_nsim, origin_select = "Walla Walla River", condition_jitter = FALSE)
save(wawa_ff_comp_median, file = here::here("figures", "alternative_spill_treatment", "final_fates", "simulation_runs", "WAWA_comp_data.rda"))
wawa_ff_comp_median_plot <- plot_final_fate_rear_type(wawa_ff_comp_median, natal_origin = "Walla Walla River")
ggsave(here::here("figures", "alternative_spill_treatment", "final_fates", "wawa_ff_comp_median_plot.png"), wawa_ff_comp_median_plot, height = 8, width = 8)


## Snake River
# Asotin Creek comparison
aso_ff_comp_median <- compare_final_fate_rear_type_SR(niter = ff_iter, nsim = ff_nsim, origin_select = "Asotin Creek", condition_jitter = FALSE)
save(aso_ff_comp_median, file = here::here("figures", "alternative_spill_treatment", "final_fates", "simulation_runs", "ASO_comp_data.rda"))
aso_ff_comp_median_plot <- plot_final_fate_rear_type(aso_ff_comp_median, natal_origin = "Asotin Creek")
ggsave(here::here("figures", "alternative_spill_treatment", "final_fates", "aso_ff_comp_median_plot.png"), aso_ff_comp_median_plot, height = 8, width = 8)

# Clearwater comparison
cle_ff_comp_median <- compare_final_fate_rear_type_SR(niter = ff_iter, nsim = ff_nsim, origin_select = "Clearwater River", condition_jitter = FALSE)
save(cle_ff_comp_median, file = here::here("figures", "alternative_spill_treatment", "final_fates", "simulation_runs", "CLE_comp_data.rda"))
cle_ff_comp_median_plot <- plot_final_fate_rear_type(cle_ff_comp_median, natal_origin = "Clearwater River")
ggsave(here::here("figures", "alternative_spill_treatment", "final_fates", "cle_ff_comp_median_plot.png"), cle_ff_comp_median_plot, height = 8, width = 8)

# Salmon comparison
sal_ff_comp_median <- compare_final_fate_rear_type_SR(niter = ff_iter, nsim = ff_nsim, origin_select = "Salmon River", condition_jitter = FALSE)
save(sal_ff_comp_median, file = here::here("figures", "alternative_spill_treatment", "final_fates", "simulation_runs", "SAL_comp_data.rda"))
sal_ff_comp_median_plot <- plot_final_fate_rear_type(sal_ff_comp_median, natal_origin = "Salmon River")
ggsave(here::here("figures", "alternative_spill_treatment", "final_fates", "sal_ff_comp_median_plot.png"), sal_ff_comp_median_plot, height = 8, width = 8)

# Grande Ronde comparison
gr_ff_comp_median <- compare_final_fate_rear_type_SR(niter = ff_iter, nsim = ff_nsim, origin_select = "Grande Ronde River", condition_jitter = FALSE)
save(gr_ff_comp_median, file = here::here("figures", "alternative_spill_treatment", "final_fates", "simulation_runs", "GR_comp_data.rda"))
gr_ff_comp_median_plot <- plot_final_fate_rear_type(gr_ff_comp_median, natal_origin = "Grande Ronde River")
ggsave(here::here("figures", "alternative_spill_treatment", "final_fates", "gr_ff_comp_median_plot.png"), gr_ff_comp_median_plot, height = 8, width = 8)

# Imnaha comparison
imn_ff_comp_median <- compare_final_fate_rear_type_SR(niter = ff_iter, nsim = ff_nsim, origin_select = "Imnaha River", condition_jitter = FALSE)
save(imn_ff_comp_median, file = here::here("figures", "alternative_spill_treatment", "final_fates", "simulation_runs", "IMN_comp_data.rda"))
imn_ff_comp_median_plot <- plot_final_fate_rear_type(imn_ff_comp_median, natal_origin = "Imnaha River")
ggsave(here::here("figures", "alternative_spill_treatment", "final_fates", "imn_ff_comp_median_plot.png"), imn_ff_comp_median_plot, height = 8, width = 8)

# Tucannon comparison
tuc_ff_comp_median <- compare_final_fate_rear_type_SR(niter = ff_iter, nsim = ff_nsim, origin_select = "Tucannon River", condition_jitter = FALSE)
save(tuc_ff_comp_median, file = here::here("figures", "alternative_spill_treatment", "final_fates", "simulation_runs", "TUC_comp_data.rda"))
tuc_ff_comp_median_plot <- plot_final_fate_rear_type(tuc_ff_comp_median, natal_origin = "Tucannon River")
ggsave(here::here("figures", "alternative_spill_treatment", "final_fates", "tuc_ff_comp_median_plot.png"), tuc_ff_comp_median_plot, height = 8, width = 8)



# SAVE ALL FINAL FATES OUTPUTS
FF_comp_data <- list(wen_ff_comp_median = wen_ff_comp_median,
                     ent_ff_comp_median = ent_ff_comp_median,
                     oka_ff_comp_median = oka_ff_comp_median,
                     met_ff_comp_median = met_ff_comp_median,
                     des_ff_comp_median = des_ff_comp_median,
                     jdr_ff_comp_median = jdr_ff_comp_median,
                     fif_ff_comp_median = fif_ff_comp_median,
                     uma_ff_comp_median = uma_ff_comp_median,
                     yak_ff_comp_median = yak_ff_comp_median,
                     wawa_ff_comp_median = wawa_ff_comp_median,
                     aso_ff_comp_median = aso_ff_comp_median,
                     cle_ff_comp_median = cle_ff_comp_median,
                     sal_ff_comp_median = sal_ff_comp_median,
                     gr_ff_comp_median = gr_ff_comp_median,
                     imn_ff_comp_median = imn_ff_comp_median,
                     tuc_ff_comp_median = tuc_ff_comp_median)

save(FF_comp_data, file = here::here("figures", "alternative_spill_treatment", "final_fates", "simulation_runs", "FF_comp_data.rda"))




#### Investigate seasonality of movement for temp0 vs. temp1 ####

# load temperature data
window_temp_summary <- read.csv(here::here("Data", "covariate_data", "model_inputs", "window_temps_summary.csv"), row.names = 1)

extract_seasonality_origin <- function(origin_select, rear_type,
                                       envir_select){
  
  states_dates <- get_origin_states_dates(envir = envir_select, origin_select = origin_select, rear = rear_type)
  
  season_upstream <- rep(0, 8)
  season_downstream <- rep(0, 8)
  upstream_counts <- rep(0, 8)
  downstream_counts <- rep(0, 8)
  
  for (i in 2:9){ #2:9 because 9 mainstem states, but we don't have conditions from mouth to BON
    
    season_upstream[i-1] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    season_downstream[i-1] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    
    upstream_counts[i-1] <- length(subset(states_dates, state == i  & direction == "upstream")$date)
    downstream_counts[i-1] <- length(subset(states_dates, state == i  & direction == "downstream")$date)
    
  }
  
  updown_seasons <- data.frame(state = model_states[2:9],
                               upstream = paste0(round(season_upstream,2), " (", upstream_counts, ")"),
                               downstream = paste0(round(season_downstream,2), " (", downstream_counts, ")"))
  
  return(updown_seasons)
}

extract_seasonality_population <- function(envir_select){
  season_data <- envir_select$data$transition_seasons_vector
  
  states_dates <- get_states_dates_direction(envir = envir_select)
  
  season_upstream <- rep(0, 8)
  season_downstream <- rep(0, 8)
  upstream_counts <- rep(0, 8)
  downstream_counts <- rep(0, 8)
  
  for (i in 2:9){ #2:9 because 9 mainstem states, but we don't have conditions from mouth to BON
    
    season_upstream[i-1] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    season_downstream[i-1] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    
    upstream_counts[i-1] <- length(subset(states_dates, state == i  & direction == "upstream")$date)
    downstream_counts[i-1] <- length(subset(states_dates, state == i  & direction == "downstream")$date)
    
  }
  
  updown_seasons <- data.frame(state = model_states[2:9],
                               upstream = paste0(round(season_upstream,2), " (", upstream_counts, ")"),
                               downstream = paste0(round(season_downstream,2), " (", downstream_counts, ")"))
  
  return(updown_seasons)
  
}


#### Export median temperature conditions on their own ####

# Median conditions will be different for every origin/rear type combination.
# We will export a different median condition for every combination of mainstem state and
# origin/rear type, as well as median conditions across the whole DPS/rear type combination.

# A function to extract median temperature upstream/downstream for each mainstem state
# for one combination of rear type, origin, and mainstem state

extract_median_temp_origin <- function(origin_select, rear_type,
                                       envir_select){
  
  temp_data <- envir_select$data$temperature_data
  
  states_dates <- get_origin_states_dates(envir = envir_select, origin_select = origin_select, rear = rear_type)
  
  temp_upstream <- rep(0, 8)
  temp_downstream <- rep(0, 8)
  upstream_counts <- rep(0, 8)
  downstream_counts <- rep(0, 8)
  
  for (i in 2:9){ #2:9 because 9 mainstem states, but we don't have conditions from mouth to BON
    temp_upstream[i-1] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream")$date,i])
    temp_downstream[i-1] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream")$date,i])
    
    upstream_counts[i-1] <- length(subset(states_dates, state == i  & direction == "upstream")$date)
    downstream_counts[i-1] <- length(subset(states_dates, state == i  & direction == "downstream")$date)
    
    # change temps to actual temps
    temp_upstream[i-1] <- as.numeric(temp_upstream[i-1]*window_temp_summary[(i-1)*2] + window_temp_summary[(i-1)*2-1])
    temp_downstream[i-1] <- as.numeric(temp_downstream[i-1]*window_temp_summary[(i-1)*2] + window_temp_summary[(i-1)*2-1])
    
    
  }
  
  updown_temps <- data.frame(state = model_states[2:9],
                             upstream = paste0(round(temp_upstream,2), " (", upstream_counts, ")"),
                             downstream = paste0(round(temp_downstream,2), " (", downstream_counts, ")"))
  
  return(updown_temps)
}

extract_median_temp_population <- function(envir_select){
  temp_data = envir_select$data$temperature_data
  
  states_dates <- get_states_dates_direction(envir = envir_select)
  
  temp_upstream <- rep(0, 8)
  temp_downstream <- rep(0, 8)
  upstream_counts <- rep(0, 8)
  downstream_counts <- rep(0, 8)
  
  for (i in 2:9){ #2:9 because 9 mainstem states, but we don't have conditions from mouth to BON
    temp_upstream[i-1] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream")$date,i])
    temp_downstream[i-1] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream")$date,i])
    
    upstream_counts[i-1] <- length(subset(states_dates, state == i  & direction == "upstream")$date)
    downstream_counts[i-1] <- length(subset(states_dates, state == i  & direction == "downstream")$date)
    
    # change temps to actual temps
    temp_upstream[i-1] <- as.numeric(temp_upstream[i-1]*window_temp_summary[(i-1)*2] + window_temp_summary[(i-1)*2-1])
    temp_downstream[i-1] <- as.numeric(temp_downstream[i-1]*window_temp_summary[(i-1)*2] + window_temp_summary[(i-1)*2-1])
    
  }
  
  updown_temps <- data.frame(state = model_states[2:9],
                             upstream = paste0(round(temp_upstream,2), " (", upstream_counts, ")"),
                             downstream = paste0(round(temp_downstream,2), " (", downstream_counts, ")"))
  
  return(updown_temps)
  
}

# Middle Columbia
DES_wild_updown_temps <- extract_median_temp_origin(origin_select = "Deschutes River", 
                                                    rear_type = "wild", 
                                                    envir_select = MCW_envir)

JDR_wild_updown_temps <- extract_median_temp_origin(origin_select = "John Day River", 
                                                    rear_type = "wild", 
                                                    envir_select = MCW_envir)

FIF_wild_updown_temps <- extract_median_temp_origin(origin_select = "Fifteenmile Creek", 
                                                    rear_type = "wild", 
                                                    envir_select = MCW_envir)

UMA_wild_updown_temps <- extract_median_temp_origin(origin_select = "Umatilla River", 
                                                    rear_type = "wild", 
                                                    envir_select = MCW_envir)

UMA_hatchery_updown_temps <- extract_median_temp_origin(origin_select = "Umatilla River", 
                                                        rear_type = "hatchery", 
                                                        envir_select = MCH_envir)

YAK_wild_updown_temps <- extract_median_temp_origin(origin_select = "Yakima River", 
                                                    rear_type = "wild", 
                                                    envir_select = MCW_envir)

WAWA_wild_updown_temps <- extract_median_temp_origin(origin_select = "Walla Walla River", 
                                                     rear_type = "wild", 
                                                     envir_select = MCW_envir)

WAWA_hatchery_updown_temps <- extract_median_temp_origin(origin_select = "Walla Walla River", 
                                                         rear_type = "hatchery", 
                                                         envir_select = MCH_envir)

# Upper Columbia
WEN_wild_updown_temps <- extract_median_temp_origin(origin_select = "Wenatchee River", 
                                                    rear_type = "wild", 
                                                    envir_select = UCW_envir)

WEN_hatchery_updown_temps <- extract_median_temp_origin(origin_select = "Wenatchee River", 
                                                        rear_type = "hatchery", 
                                                        envir_select = UCH_envir)

ENT_wild_updown_temps <- extract_median_temp_origin(origin_select = "Entiat River", 
                                                    rear_type = "wild", 
                                                    envir_select = UCW_envir)

OKA_hatchery_updown_temps <- extract_median_temp_origin(origin_select = "Okanogan River", 
                                                        rear_type = "hatchery", 
                                                        envir_select = UCH_envir)

MET_wild_updown_temps <- extract_median_temp_origin(origin_select = "Methow River", 
                                                    rear_type = "wild", 
                                                    envir_select = UCW_envir)

MET_hatchery_updown_temps <- extract_median_temp_origin(origin_select = "Methow River", 
                                                        rear_type = "hatchery", 
                                                        envir_select = UCH_envir)

# Snake River
ASO_wild_updown_temps <- extract_median_temp_origin(origin_select = "Asotin Creek", 
                                                    rear_type = "wild", 
                                                    envir_select = SRW_envir)

CLE_wild_updown_temps <- extract_median_temp_origin(origin_select = "Clearwater River", 
                                                    rear_type = "wild", 
                                                    envir_select = SRW_envir)

CLE_hatchery_updown_temps <- extract_median_temp_origin(origin_select = "Clearwater River", 
                                                        rear_type = "hatchery", 
                                                        envir_select = SRH_envir)

IMN_wild_updown_temps <- extract_median_temp_origin(origin_select = "Imnaha River", 
                                                    rear_type = "wild", 
                                                    envir_select = SRW_envir)

IMN_hatchery_updown_temps <- extract_median_temp_origin(origin_select = "Imnaha River", 
                                                        rear_type = "hatchery", 
                                                        envir_select = SRH_envir)

GR_wild_updown_temps <- extract_median_temp_origin(origin_select = "Grande Ronde River", 
                                                   rear_type = "wild", 
                                                   envir_select = SRW_envir)

GR_hatchery_updown_temps <- extract_median_temp_origin(origin_select = "Grande Ronde River", 
                                                       rear_type = "hatchery", 
                                                       envir_select = SRH_envir)

SAL_wild_updown_temps <- extract_median_temp_origin(origin_select = "Salmon River", 
                                                    rear_type = "wild", 
                                                    envir_select = SRW_envir)

SAL_hatchery_updown_temps <- extract_median_temp_origin(origin_select = "Salmon River", 
                                                        rear_type = "hatchery", 
                                                        envir_select = SRH_envir)

TUC_wild_updown_temps <- extract_median_temp_origin(origin_select = "Tucannon River", 
                                                    rear_type = "wild", 
                                                    envir_select = SRW_envir)

TUC_hatchery_updown_temps <- extract_median_temp_origin(origin_select = "Tucannon River", 
                                                        rear_type = "hatchery", 
                                                        envir_select = SRH_envir)


# Population wide
MCW_updown_temps <- extract_median_temp_population(envir_select = MCW_envir)
MCH_updown_temps <- extract_median_temp_population(envir_select = MCH_envir)
UCW_updown_temps <- extract_median_temp_population(envir_select = UCW_envir)
UCH_updown_temps <- extract_median_temp_population(envir_select = UCH_envir)
SRW_updown_temps <- extract_median_temp_population(envir_select = SRW_envir)
SRH_updown_temps <- extract_median_temp_population(envir_select = SRH_envir)

# export all of these median temps as .rda, for download to local

origin_rear_updown_temps <- list(DES_wild_updown_temps = DES_wild_updown_temps,
                                 JDR_wild_updown_temps = JDR_wild_updown_temps,
                                 FIF_wild_updown_temps = FIF_wild_updown_temps,
                                 UMA_wild_updown_temps = UMA_wild_updown_temps,
                                 UMA_hatchery_updown_temps = UMA_hatchery_updown_temps,
                                 YAK_wild_updown_temps = YAK_wild_updown_temps,
                                 WAWA_wild_updown_temps = WAWA_wild_updown_temps,
                                 WAWA_hatchery_updown_temps = WAWA_hatchery_updown_temps,
                                 WEN_wild_updown_temps = WEN_wild_updown_temps,
                                 WEN_hatchery_updown_temps = WEN_hatchery_updown_temps,
                                 ENT_wild_updown_temps = ENT_wild_updown_temps,
                                 OKA_hatchery_updown_temps = OKA_hatchery_updown_temps,
                                 MET_wild_updown_temps = MET_wild_updown_temps,
                                 MET_hatchery_updown_temps = MET_hatchery_updown_temps,
                                 ASO_wild_updown_temps = ASO_wild_updown_temps,
                                 CLE_wild_updown_temps = CLE_wild_updown_temps,
                                 CLE_hatchery_updown_temps = CLE_hatchery_updown_temps,
                                 IMN_wild_updown_temps = IMN_wild_updown_temps,
                                 IMN_hatchery_updown_temps = IMN_hatchery_updown_temps,
                                 GR_wild_updown_temps = GR_wild_updown_temps,
                                 GR_hatchery_updown_temps = GR_hatchery_updown_temps,
                                 SAL_wild_updown_temps = SAL_wild_updown_temps,
                                 SAL_hatchery_updown_temps = SAL_hatchery_updown_temps,
                                 TUC_wild_updown_temps = TUC_wild_updown_temps,
                                 TUC_hatchery_updown_temps = TUC_hatchery_updown_temps)

DPS_rear_updown_temps <- list(MCW_updown_temps = MCW_updown_temps,
                              MCH_updown_temps = MCH_updown_temps,
                              UCW_updown_temps = UCW_updown_temps,
                              UCH_updown_temps = UCH_updown_temps,
                              SRW_updown_temps = SRW_updown_temps,
                              SRH_updown_temps = SRH_updown_temps)

save(origin_rear_updown_temps, file = here::here("figures", "alternative_spill_treatment", "final_fates", "origin_rear_updown_temps.rda"))
save(DPS_rear_updown_temps, file = here::here("figures", "alternative_spill_treatment", "final_fates", "DPS_rear_updown_temps.rda"))

#### Export median spill window conditions on their own ####

# Median conditions will be different for every origin/rear type combination.
# We will export a different median condition for every combination of mainstem state and
# origin/rear type, as well as median conditions across the whole DPS/rear type combination.

# A function to extract median spill window upstream/downstream for each mainstem state
# for one combination of rear type, origin, and mainstem state

extract_median_spillwindow_origin <- function(origin_select, rear_type,
                                              envir_select){
  
  spillwindow_data <- envir_select$data$spill_window_data
  
  states_dates <- get_origin_states_dates(envir = envir_select, origin_select = origin_select, rear = rear_type)
  
  spillwindow_upstream <- rep(0, 8)
  spillwindow_downstream <- rep(0, 8)
  upstream_counts <- rep(0, 8)
  downstream_counts <- rep(0, 8)
  
  for (i in 2:9){ #2:9 because 9 mainstem states, but we don't have conditions from mouth to BON
    spillwindow_upstream[i-1] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream")$date,i])
    spillwindow_downstream[i-1] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream")$date,i])
    
    upstream_counts[i-1] <- length(subset(states_dates, state == i  & direction == "upstream")$date)
    downstream_counts[i-1] <- length(subset(states_dates, state == i  & direction == "downstream")$date)
    
    # change spillwindows to actual spillwindows
    spillwindow_upstream[i-1] <- spillwindow_upstream[i-1]*100
    spillwindow_downstream[i-1] <- spillwindow_downstream[i-1]*100
    
    
  }
  
  updown_spillwindows <- data.frame(state = model_states[2:9],
                                    upstream = paste0(round(spillwindow_upstream,2), " (", upstream_counts, ")"),
                                    downstream = paste0(round(spillwindow_downstream,2), " (", downstream_counts, ")"))
  
  return(updown_spillwindows)
}

extract_median_spillwindow_population <- function(envir_select){
  spillwindow_data = envir_select$data$spill_window_data
  
  states_dates <- get_states_dates_direction(envir = envir_select)
  
  spillwindow_upstream <- rep(0, 8)
  spillwindow_downstream <- rep(0, 8)
  upstream_counts <- rep(0, 8)
  downstream_counts <- rep(0, 8)
  
  for (i in 2:9){ #2:9 because 9 mainstem states, but we don't have conditions from mouth to BON
    spillwindow_upstream[i-1] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream")$date,i])
    spillwindow_downstream[i-1] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream")$date,i])
    
    upstream_counts[i-1] <- length(subset(states_dates, state == i  & direction == "upstream")$date)
    downstream_counts[i-1] <- length(subset(states_dates, state == i  & direction == "downstream")$date)
    
    # change spillwindows to actual spillwindows
    spillwindow_upstream[i-1] <- spillwindow_upstream[i-1]*100
    spillwindow_downstream[i-1] <- spillwindow_downstream[i-1]*100
    
  }
  
  updown_spillwindows <- data.frame(state = model_states[2:9],
                                    upstream = paste0(round(spillwindow_upstream,2), " (", upstream_counts, ")"),
                                    downstream = paste0(round(spillwindow_downstream,2), " (", downstream_counts, ")"))
  
  return(updown_spillwindows)
  
}

# Middle Columbia
DES_wild_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Deschutes River", 
                                                                  rear_type = "wild", 
                                                                  envir_select = MCW_envir)

JDR_wild_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "John Day River", 
                                                                  rear_type = "wild", 
                                                                  envir_select = MCW_envir)

FIF_wild_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Fifteenmile Creek", 
                                                                  rear_type = "wild", 
                                                                  envir_select = MCW_envir)

UMA_wild_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Umatilla River", 
                                                                  rear_type = "wild", 
                                                                  envir_select = MCW_envir)

UMA_hatchery_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Umatilla River", 
                                                                      rear_type = "hatchery", 
                                                                      envir_select = MCH_envir)

YAK_wild_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Yakima River", 
                                                                  rear_type = "wild", 
                                                                  envir_select = MCW_envir)

WAWA_wild_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Walla Walla River", 
                                                                   rear_type = "wild", 
                                                                   envir_select = MCW_envir)

WAWA_hatchery_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Walla Walla River", 
                                                                       rear_type = "hatchery", 
                                                                       envir_select = MCH_envir)

# Upper Columbia
WEN_wild_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Wenatchee River", 
                                                                  rear_type = "wild", 
                                                                  envir_select = UCW_envir)

WEN_hatchery_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Wenatchee River", 
                                                                      rear_type = "hatchery", 
                                                                      envir_select = UCH_envir)

ENT_wild_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Entiat River", 
                                                                  rear_type = "wild", 
                                                                  envir_select = UCW_envir)

OKA_hatchery_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Okanogan River", 
                                                                      rear_type = "hatchery", 
                                                                      envir_select = UCH_envir)

MET_wild_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Methow River", 
                                                                  rear_type = "wild", 
                                                                  envir_select = UCW_envir)

MET_hatchery_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Methow River", 
                                                                      rear_type = "hatchery", 
                                                                      envir_select = UCH_envir)

# Snake River
ASO_wild_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Asotin Creek", 
                                                                  rear_type = "wild", 
                                                                  envir_select = SRW_envir)

CLE_wild_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Clearwater River", 
                                                                  rear_type = "wild", 
                                                                  envir_select = SRW_envir)

CLE_hatchery_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Clearwater River", 
                                                                      rear_type = "hatchery", 
                                                                      envir_select = SRH_envir)

IMN_wild_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Imnaha River", 
                                                                  rear_type = "wild", 
                                                                  envir_select = SRW_envir)

IMN_hatchery_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Imnaha River", 
                                                                      rear_type = "hatchery", 
                                                                      envir_select = SRH_envir)

GR_wild_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Grande Ronde River", 
                                                                 rear_type = "wild", 
                                                                 envir_select = SRW_envir)

GR_hatchery_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Grande Ronde River", 
                                                                     rear_type = "hatchery", 
                                                                     envir_select = SRH_envir)

SAL_wild_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Salmon River", 
                                                                  rear_type = "wild", 
                                                                  envir_select = SRW_envir)

SAL_hatchery_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Salmon River", 
                                                                      rear_type = "hatchery", 
                                                                      envir_select = SRH_envir)

TUC_wild_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Tucannon River", 
                                                                  rear_type = "wild", 
                                                                  envir_select = SRW_envir)

TUC_hatchery_updown_spillwindows <- extract_median_spillwindow_origin(origin_select = "Tucannon River", 
                                                                      rear_type = "hatchery", 
                                                                      envir_select = SRH_envir)


# Population wide
MCW_updown_spillwindows <- extract_median_spillwindow_population(envir_select = MCW_envir)
MCH_updown_spillwindows <- extract_median_spillwindow_population(envir_select = MCH_envir)
UCW_updown_spillwindows <- extract_median_spillwindow_population(envir_select = UCW_envir)
UCH_updown_spillwindows <- extract_median_spillwindow_population(envir_select = UCH_envir)
SRW_updown_spillwindows <- extract_median_spillwindow_population(envir_select = SRW_envir)
SRH_updown_spillwindows <- extract_median_spillwindow_population(envir_select = SRH_envir)

# export all of these median spillwindows as .rda, for download to local

origin_rear_updown_spillwindows <- list(DES_wild_updown_spillwindows = DES_wild_updown_spillwindows,
                                        JDR_wild_updown_spillwindows = JDR_wild_updown_spillwindows,
                                        FIF_wild_updown_spillwindows = FIF_wild_updown_spillwindows,
                                        UMA_wild_updown_spillwindows = UMA_wild_updown_spillwindows,
                                        UMA_hatchery_updown_spillwindows = UMA_hatchery_updown_spillwindows,
                                        YAK_wild_updown_spillwindows = YAK_wild_updown_spillwindows,
                                        WAWA_wild_updown_spillwindows = WAWA_wild_updown_spillwindows,
                                        WAWA_hatchery_updown_spillwindows = WAWA_hatchery_updown_spillwindows,
                                        WEN_wild_updown_spillwindows = WEN_wild_updown_spillwindows,
                                        WEN_hatchery_updown_spillwindows = WEN_hatchery_updown_spillwindows,
                                        ENT_wild_updown_spillwindows = ENT_wild_updown_spillwindows,
                                        OKA_hatchery_updown_spillwindows = OKA_hatchery_updown_spillwindows,
                                        MET_wild_updown_spillwindows = MET_wild_updown_spillwindows,
                                        MET_hatchery_updown_spillwindows = MET_hatchery_updown_spillwindows,
                                        ASO_wild_updown_spillwindows = ASO_wild_updown_spillwindows,
                                        CLE_wild_updown_spillwindows = CLE_wild_updown_spillwindows,
                                        CLE_hatchery_updown_spillwindows = CLE_hatchery_updown_spillwindows,
                                        IMN_wild_updown_spillwindows = IMN_wild_updown_spillwindows,
                                        IMN_hatchery_updown_spillwindows = IMN_hatchery_updown_spillwindows,
                                        GR_wild_updown_spillwindows = GR_wild_updown_spillwindows,
                                        GR_hatchery_updown_spillwindows = GR_hatchery_updown_spillwindows,
                                        SAL_wild_updown_spillwindows = SAL_wild_updown_spillwindows,
                                        SAL_hatchery_updown_spillwindows = SAL_hatchery_updown_spillwindows,
                                        TUC_wild_updown_spillwindows = TUC_wild_updown_spillwindows,
                                        TUC_hatchery_updown_spillwindows = TUC_hatchery_updown_spillwindows)

DPS_rear_updown_spillwindows <- list(MCW_updown_spillwindows = MCW_updown_spillwindows,
                                     MCH_updown_spillwindows = MCH_updown_spillwindows,
                                     UCW_updown_spillwindows = UCW_updown_spillwindows,
                                     UCH_updown_spillwindows = UCH_updown_spillwindows,
                                     SRW_updown_spillwindows = SRW_updown_spillwindows,
                                     SRH_updown_spillwindows = SRH_updown_spillwindows)

save(origin_rear_updown_spillwindows, file = here::here("figures", "alternative_spill_treatment", "final_fates", "origin_rear_updown_spillwindows.rda"))
save(DPS_rear_updown_spillwindows, file = here::here("figures", "alternative_spill_treatment", "final_fates", "DPS_rear_updown_spillwindows.rda"))


#### Analyze spill vs. temperature data and correlations ####

get_states_dates_cov <- function(envir_select){
  temp_data = envir_select$data$temperature_data
  spillwindow_data = envir_select$data$spill_window_data
  
  state_index <- data.frame(dam = c("MOUTH", "BON", "MCN", "PRA", "RIS", "RRE", "WEL_quick", "ICH", "LGR_quick"),
                            state_numeric = seq(1,9))
  
  temp_data %>% 
    as.data.frame() %>% 
    rownames_to_column("date") %>% 
    mutate(date = as.numeric(date)) %>% 
    pivot_longer(., cols = -c(date), values_to = "temp", names_to = "dam") %>% 
    left_join(state_index, by = "dam") %>% 
    filter(!(is.na(state_numeric))) -> temp_long
  
  spillwindow_data %>% 
    as.data.frame() %>% 
    rownames_to_column("date") %>% 
    mutate(date = as.numeric(date)) %>% 
    pivot_longer(., cols = -c(date), values_to = "spillwindow", names_to = "dam") %>% 
    left_join(state_index, by = "dam") %>% 
    filter(!(is.na(state_numeric)))-> spillwindow_long
  
  states_dates <- get_states_dates_direction(envir = envir_select)
  
  # add the conditions for each date
  states_dates %>% 
    dplyr::rename(state_numeric = state) %>% 
    left_join(dplyr::select(temp_long, date, temp, state_numeric), by = c("state_numeric", "date")) %>% 
    left_join(dplyr::select(spillwindow_long, date, spillwindow, state_numeric), by = c("state_numeric", "date")) -> states_dates_cov
  
  
  # drop the non-mainstem states (keep only those with covariate data)
  # rename states
  model_states_df <- data.frame(state_numeric = 1:41, state = model_states)
  states_dates_cov %>% 
    filter(state_numeric %in% 2:9) %>% 
    left_join(model_states_df, by = "state_numeric") -> states_dates_cov
  
  return(states_dates_cov)
  
}

MCW_states_dates_cov <- get_states_dates_cov(envir_select = MCW_envir)

MCW_states_dates_cov_plot <- ggplot(MCW_states_dates_cov, aes(x = spillwindow, y = temp)) +
  geom_point() +
  facet_wrap(~state) +
  ggtitle("Middle Columbia, Wild")

ggsave(here::here("figures", "alternative_spill_treatment", "data_plots", "MCW_states_dates_cov_plot.png"), MCW_states_dates_cov_plot, height = 8, width = 8)

MCH_states_dates_cov <- get_states_dates_cov(envir_select = MCH_envir)

MCH_states_dates_cov_plot <- ggplot(MCH_states_dates_cov, aes(x = spillwindow, y = temp)) +
  geom_point() +
  facet_wrap(~state) +
  ggtitle("Middle Columbia, Hatchery")

ggsave(here::here("figures", "alternative_spill_treatment", "data_plots", "MCH_states_dates_cov_plot.png"), MCH_states_dates_cov_plot, height = 8, width = 8)

UCW_states_dates_cov <- get_states_dates_cov(envir_select = UCW_envir)

UCW_states_dates_cov_plot <- ggplot(UCW_states_dates_cov, aes(x = spillwindow, y = temp)) +
  geom_point() +
  facet_wrap(~state) +
  ggtitle("Upper Columbia, Wild")

ggsave(here::here("figures", "alternative_spill_treatment", "data_plots", "UCW_states_dates_cov_plot.png"), UCW_states_dates_cov_plot, height = 8, width = 8)

UCH_states_dates_cov <- get_states_dates_cov(envir_select = UCH_envir)

UCH_states_dates_cov_plot <- ggplot(UCH_states_dates_cov, aes(x = spillwindow, y = temp)) +
  geom_point() +
  facet_wrap(~state) +
  ggtitle("Upper Columbia, Hatchery")

ggsave(here::here("figures", "alternative_spill_treatment", "data_plots", "UCH_states_dates_cov_plot.png"), UCH_states_dates_cov_plot, height = 8, width = 8)

SRW_states_dates_cov <- get_states_dates_cov(envir_select = SRW_envir)

SRW_states_dates_cov_plot <- ggplot(SRW_states_dates_cov, aes(x = spillwindow, y = temp)) +
  geom_point() +
  facet_wrap(~state) +
  ggtitle("Snake River, Wild")

ggsave(here::here("figures", "alternative_spill_treatment", "data_plots", "SRW_states_dates_cov_plot.png"), SRW_states_dates_cov_plot, height = 8, width = 8)

SRH_states_dates_cov <- get_states_dates_cov(envir_select = SRH_envir)

SRH_states_dates_cov_plot <- ggplot(SRH_states_dates_cov, aes(x = spillwindow, y = temp)) +
  geom_point() +
  facet_wrap(~state) +
  ggtitle("Snake River, Hatchery")

ggsave(here::here("figures", "alternative_spill_treatment", "data_plots", "SRH_states_dates_cov_plot.png"), SRH_states_dates_cov_plot, height = 8, width = 8)


