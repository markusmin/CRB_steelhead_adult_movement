# 07-03_final_fates_simulation

# This script takes the output from the stan model runs 
# generates estimates of final fate distributions

# First, need to load in all of the model runs and all of the packages.
source("R/07_model_analysis/07-01_load_stan_models.R")
source("R/07_model_analysis/transport/07-01_load_stan_models_transport.R")

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

final_fates_simulation_individual_MCW <- function(nsim,
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
  
  # 2025-07-28: code below is all new, final fates but tracking individual fish
  # create a matrix to track fish
  # every fish starts in the start state
  indiv_fish_state_matrix <- matrix(data = start_state, nrow = nsim, ncol = 1)
  colnames(indiv_fish_state_matrix) <- "t1"
  
  # also track the direction of each fish
  # 1 = upstream, 2 = downstream, 0 = dead
  # every fish is going upstream to start
  indiv_fish_direction_matrix <- matrix(data = 1, nrow = nsim, ncol = 1)
  colnames(indiv_fish_direction_matrix) <- "t1"
  
  # use a while loop to keep running until all fish have entered loss state
  i <- 1
  condition <- "not finished"
  while (condition != "finished"){
    # if there are still fish moving, add a column to track their next movement
    next_time_step <- matrix(nrow = nsim, ncol = 1)
    colnames(next_time_step) <- paste0("t", i+1)
    
    indiv_fish_state_matrix <- cbind(indiv_fish_state_matrix, next_time_step)
    
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_moving_fish <- which(indiv_fish_direction_matrix[,i] == 1)
    if(sum(upstream_moving_fish)>0){
      for (j in upstream_moving_fish){
        indiv_fish_state_matrix[j,i+1] <- which(rmultinom(n = 1, size = 1, prob = upstream_move_prob_matrix[indiv_fish_state_matrix[j,i],]) == 1)
      }
    }
    
    
    downstream_moving_fish <- which(indiv_fish_direction_matrix[,i] == 2)
    if(sum(downstream_moving_fish)>0){
      for (j in downstream_moving_fish){
        indiv_fish_state_matrix[j,i+1] <- which(rmultinom(n = 1, size = 1, prob = downstream_move_prob_matrix[indiv_fish_state_matrix[j,i],]) == 1)
      }
    }
    
    # now, populate direction of movement for each fish
    indiv_fish_direction_matrix <- cbind(indiv_fish_direction_matrix, next_time_step)
    
    indiv_fish_direction_matrix[,i+1] <- ifelse(indiv_fish_state_matrix[,i+1] > indiv_fish_state_matrix[,i], 1, 2)
    
    # if it's in the loss state, direction = 0
    indiv_fish_direction_matrix[which(indiv_fish_state_matrix[,i+1] == 41),i+1] <- 0
    
    # if they're all in the loss state (direction = 0), then we're finish
    if (sum(indiv_fish_direction_matrix[,i+1], na.rm = T) == 0){
      condition <- "finished"
    }
    
    # increase counter
    i <- i + 1
  }
  
  
  return(indiv_fish_state_matrix)
  
}

final_fates_simulation_individual_MCH <- function(nsim,
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
  
  # 2025-07-28: code below is all new, final fates but tracking individual fish
  # create a matrix to track fish
  # every fish starts in the start state
  indiv_fish_state_matrix <- matrix(data = start_state, nrow = nsim, ncol = 1)
  colnames(indiv_fish_state_matrix) <- "t1"
  
  # also track the direction of each fish
  # 1 = upstream, 2 = downstream, 0 = dead
  # every fish is going upstream to start
  indiv_fish_direction_matrix <- matrix(data = 1, nrow = nsim, ncol = 1)
  colnames(indiv_fish_direction_matrix) <- "t1"
  
  # use a while loop to keep running until all fish have entered loss state
  i <- 1
  condition <- "not finished"
  while (condition != "finished"){
    # if there are still fish moving, add a column to track their next movement
    next_time_step <- matrix(nrow = nsim, ncol = 1)
    colnames(next_time_step) <- paste0("t", i+1)
    
    indiv_fish_state_matrix <- cbind(indiv_fish_state_matrix, next_time_step)
    
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_moving_fish <- which(indiv_fish_direction_matrix[,i] == 1)
    if(sum(upstream_moving_fish)>0){
      for (j in upstream_moving_fish){
        indiv_fish_state_matrix[j,i+1] <- which(rmultinom(n = 1, size = 1, prob = upstream_move_prob_matrix[indiv_fish_state_matrix[j,i],]) == 1)
      }
    }
    
    
    downstream_moving_fish <- which(indiv_fish_direction_matrix[,i] == 2)
    if(sum(downstream_moving_fish)>0){
      for (j in downstream_moving_fish){
        indiv_fish_state_matrix[j,i+1] <- which(rmultinom(n = 1, size = 1, prob = downstream_move_prob_matrix[indiv_fish_state_matrix[j,i],]) == 1)
      }
    }
    
    # now, populate direction of movement for each fish
    indiv_fish_direction_matrix <- cbind(indiv_fish_direction_matrix, next_time_step)
    
    indiv_fish_direction_matrix[,i+1] <- ifelse(indiv_fish_state_matrix[,i+1] > indiv_fish_state_matrix[,i], 1, 2)
    
    # if it's in the loss state, direction = 0
    indiv_fish_direction_matrix[which(indiv_fish_state_matrix[,i+1] == 41),i+1] <- 0
    
    # if they're all in the loss state (direction = 0), then we're finish
    if (sum(indiv_fish_direction_matrix[,i+1], na.rm = T) == 0){
      condition <- "finished"
    }
    
    # increase counter
    i <- i + 1
  }
  
  return(indiv_fish_state_matrix)
  
}


final_fates_simulation_individual_UCW <- function(nsim,
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
  
  # 2025-07-28: code below is all new, final fates but tracking individual fish
  # create a matrix to track fish
  # every fish starts in the start state
  indiv_fish_state_matrix <- matrix(data = start_state, nrow = nsim, ncol = 1)
  colnames(indiv_fish_state_matrix) <- "t1"
  
  # also track the direction of each fish
  # 1 = upstream, 2 = downstream, 0 = dead
  # every fish is going upstream to start
  indiv_fish_direction_matrix <- matrix(data = 1, nrow = nsim, ncol = 1)
  colnames(indiv_fish_direction_matrix) <- "t1"
  
  # use a while loop to keep running until all fish have entered loss state
  i <- 1
  condition <- "not finished"
  while (condition != "finished"){
    # if there are still fish moving, add a column to track their next movement
    next_time_step <- matrix(nrow = nsim, ncol = 1)
    colnames(next_time_step) <- paste0("t", i+1)
    
    indiv_fish_state_matrix <- cbind(indiv_fish_state_matrix, next_time_step)
    
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_moving_fish <- which(indiv_fish_direction_matrix[,i] == 1)
    if(sum(upstream_moving_fish)>0){
      for (j in upstream_moving_fish){
        indiv_fish_state_matrix[j,i+1] <- which(rmultinom(n = 1, size = 1, prob = upstream_move_prob_matrix[indiv_fish_state_matrix[j,i],]) == 1)
      }
    }
    
    
    downstream_moving_fish <- which(indiv_fish_direction_matrix[,i] == 2)
    if(sum(downstream_moving_fish)>0){
      for (j in downstream_moving_fish){
        indiv_fish_state_matrix[j,i+1] <- which(rmultinom(n = 1, size = 1, prob = downstream_move_prob_matrix[indiv_fish_state_matrix[j,i],]) == 1)
      }
    }
    
    # now, populate direction of movement for each fish
    indiv_fish_direction_matrix <- cbind(indiv_fish_direction_matrix, next_time_step)
    
    indiv_fish_direction_matrix[,i+1] <- ifelse(indiv_fish_state_matrix[,i+1] > indiv_fish_state_matrix[,i], 1, 2)
    
    # if it's in the loss state, direction = 0
    indiv_fish_direction_matrix[which(indiv_fish_state_matrix[,i+1] == 41),i+1] <- 0
    
    # if they're all in the loss state (direction = 0), then we're finish
    if (sum(indiv_fish_direction_matrix[,i+1], na.rm = T) == 0){
      condition <- "finished"
    }
    
    # increase counter
    i <- i + 1
  }
  
  return(indiv_fish_state_matrix)
  
}

final_fates_simulation_individual_UCH <- function(nsim,
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
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
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
  
  # 2025-07-28: code below is all new, final fates but tracking individual fish
  # create a matrix to track fish
  # every fish starts in the start state
  indiv_fish_state_matrix <- matrix(data = start_state, nrow = nsim, ncol = 1)
  colnames(indiv_fish_state_matrix) <- "t1"
  
  # also track the direction of each fish
  # 1 = upstream, 2 = downstream, 0 = dead
  # every fish is going upstream to start
  indiv_fish_direction_matrix <- matrix(data = 1, nrow = nsim, ncol = 1)
  colnames(indiv_fish_direction_matrix) <- "t1"
  
  # use a while loop to keep running until all fish have entered loss state
  i <- 1
  condition <- "not finished"
  while (condition != "finished"){
    # if there are still fish moving, add a column to track their next movement
    next_time_step <- matrix(nrow = nsim, ncol = 1)
    colnames(next_time_step) <- paste0("t", i+1)
    
    indiv_fish_state_matrix <- cbind(indiv_fish_state_matrix, next_time_step)
    
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_moving_fish <- which(indiv_fish_direction_matrix[,i] == 1)
    if(sum(upstream_moving_fish)>0){
      for (j in upstream_moving_fish){
        indiv_fish_state_matrix[j,i+1] <- which(rmultinom(n = 1, size = 1, prob = upstream_move_prob_matrix[indiv_fish_state_matrix[j,i],]) == 1)
      }
    }
    
    
    downstream_moving_fish <- which(indiv_fish_direction_matrix[,i] == 2)
    if(sum(downstream_moving_fish)>0){
      for (j in downstream_moving_fish){
        indiv_fish_state_matrix[j,i+1] <- which(rmultinom(n = 1, size = 1, prob = downstream_move_prob_matrix[indiv_fish_state_matrix[j,i],]) == 1)
      }
    }
    
    # now, populate direction of movement for each fish
    indiv_fish_direction_matrix <- cbind(indiv_fish_direction_matrix, next_time_step)
    
    indiv_fish_direction_matrix[,i+1] <- ifelse(indiv_fish_state_matrix[,i+1] > indiv_fish_state_matrix[,i], 1, 2)
    
    # if it's in the loss state, direction = 0
    indiv_fish_direction_matrix[which(indiv_fish_state_matrix[,i+1] == 41),i+1] <- 0
    
    # if they're all in the loss state (direction = 0), then we're finish
    if (sum(indiv_fish_direction_matrix[,i+1], na.rm = T) == 0){
      condition <- "finished"
    }
    
    # increase counter
    i <- i + 1
  }
  
  return(indiv_fish_state_matrix)
  
}


final_fates_simulation_individual_SRW_T <- function(nsim,
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
      sample_date[i] <- sample(subset(SRW_T_states_dates, state == i)$date, 1)
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
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRW_T_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRW_T_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRW_T_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRW_T_states_dates, state == i  & direction == "downstream")$season)
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
        if(nrow(subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
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
        if(nrow(subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
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
    possible_movements <- SRW_T_envir$data$movements[, "col"][SRW_T_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_T[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRW_T[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRW_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    btemp1xorigin6_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin6 +
                                                                    borigin1_array_SRW_T[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW_T[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW_T[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW_T[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW_T[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW_T[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_T[i,possible_movements,iter] +
                    btemp1_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRW_T[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRW_T[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin6 +
                    borigin1_array_SRW_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_T[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_T[i,possible_movements,iter]*origin6))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_T[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRW_T[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRW_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    btemp0xorigin6_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin6 +
                                                                    borigin1_array_SRW_T[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW_T[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW_T[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW_T[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW_T[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW_T[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_T[i,possible_movements,iter] +
                    btemp0_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRW_T[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRW_T[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin6 +
                    borigin1_array_SRW_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_T[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_T[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRW_T_envir$data$movements[, "col"][SRW_T_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_T[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRW_T[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRW_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      btemp1xorigin6_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin6 +
                                                                      borigin1_array_SRW_T[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW_T[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW_T[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW_T[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW_T[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW_T[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_T[i,possible_movements,iter] +
                    btemp1_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRW_T[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRW_T[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin6 +
                    borigin1_array_SRW_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_T[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_T[i,possible_movements,iter]*origin6))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_T[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRW_T[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRW_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      btemp0xorigin6_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin6 +
                                                                      borigin1_array_SRW_T[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW_T[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW_T[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW_T[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW_T[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW_T[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_T[i,possible_movements,iter] +
                    btemp0_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRW_T[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRW_T[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin6 +
                    borigin1_array_SRW_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_T[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_T[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  # 2025-07-28: code below is all new, final fates but tracking individual fish
  # create a matrix to track fish
  # every fish starts in the start state
  indiv_fish_state_matrix <- matrix(data = start_state, nrow = nsim, ncol = 1)
  colnames(indiv_fish_state_matrix) <- "t1"
  
  # also track the direction of each fish
  # 1 = upstream, 2 = downstream, 0 = dead
  # every fish is going upstream to start
  indiv_fish_direction_matrix <- matrix(data = 1, nrow = nsim, ncol = 1)
  colnames(indiv_fish_direction_matrix) <- "t1"
  
  # use a while loop to keep running until all fish have entered loss state
  i <- 1
  condition <- "not finished"
  while (condition != "finished"){
    # if there are still fish moving, add a column to track their next movement
    next_time_step <- matrix(nrow = nsim, ncol = 1)
    colnames(next_time_step) <- paste0("t", i+1)
    
    indiv_fish_state_matrix <- cbind(indiv_fish_state_matrix, next_time_step)
    
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_moving_fish <- which(indiv_fish_direction_matrix[,i] == 1)
    if(sum(upstream_moving_fish)>0){
      for (j in upstream_moving_fish){
        indiv_fish_state_matrix[j,i+1] <- which(rmultinom(n = 1, size = 1, prob = upstream_move_prob_matrix[indiv_fish_state_matrix[j,i],]) == 1)
      }
    }
    
    
    downstream_moving_fish <- which(indiv_fish_direction_matrix[,i] == 2)
    if(sum(downstream_moving_fish)>0){
      for (j in downstream_moving_fish){
        indiv_fish_state_matrix[j,i+1] <- which(rmultinom(n = 1, size = 1, prob = downstream_move_prob_matrix[indiv_fish_state_matrix[j,i],]) == 1)
      }
    }
    
    # now, populate direction of movement for each fish
    indiv_fish_direction_matrix <- cbind(indiv_fish_direction_matrix, next_time_step)
    
    indiv_fish_direction_matrix[,i+1] <- ifelse(indiv_fish_state_matrix[,i+1] > indiv_fish_state_matrix[,i], 1, 2)
    
    # if it's in the loss state, direction = 0
    indiv_fish_direction_matrix[which(indiv_fish_state_matrix[,i+1] == 41),i+1] <- 0
    
    # if they're all in the loss state (direction = 0), then we're finish
    if (sum(indiv_fish_direction_matrix[,i+1], na.rm = T) == 0){
      condition <- "finished"
    }

    # increase counter
    i <- i + 1
  }
  
  return(indiv_fish_state_matrix)
  
}

final_fates_simulation_individual_SRW_NT <- function(nsim,
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
      sample_date[i] <- sample(subset(SRW_NT_states_dates, state == i)$date, 1)
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
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRW_NT_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRW_NT_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRW_NT_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRW_NT_states_dates, state == i  & direction == "downstream")$season)
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
        if(nrow(subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
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
        if(nrow(subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
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
    possible_movements <- SRW_NT_envir$data$movements[, "col"][SRW_NT_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_NT[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRW_NT[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRW_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    btemp1xorigin6_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin6 +
                                                                    borigin1_array_SRW_NT[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW_NT[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW_NT[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW_NT[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW_NT[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW_NT[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_NT[i,possible_movements,iter] +
                    btemp1_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRW_NT[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRW_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin6 +
                    borigin1_array_SRW_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_NT[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_NT[i,possible_movements,iter]*origin6))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_NT[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRW_NT[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRW_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    btemp0xorigin6_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin6 +
                                                                    borigin1_array_SRW_NT[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW_NT[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW_NT[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW_NT[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW_NT[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW_NT[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_NT[i,possible_movements,iter] +
                    btemp0_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRW_NT[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRW_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin6 +
                    borigin1_array_SRW_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_NT[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_NT[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRW_NT_envir$data$movements[, "col"][SRW_NT_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_NT[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRW_NT[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRW_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      btemp1xorigin6_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin6 +
                                                                      borigin1_array_SRW_NT[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW_NT[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW_NT[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW_NT[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW_NT[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW_NT[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_NT[i,possible_movements,iter] +
                    btemp1_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRW_NT[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRW_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin6 +
                    borigin1_array_SRW_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_NT[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_NT[i,possible_movements,iter]*origin6))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_NT[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRW_NT[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRW_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      btemp0xorigin6_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin6 +
                                                                      borigin1_array_SRW_NT[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW_NT[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW_NT[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW_NT[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW_NT[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW_NT[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_NT[i,possible_movements,iter] +
                    btemp0_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRW_NT[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRW_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin6 +
                    borigin1_array_SRW_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_NT[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_NT[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  # 2025-07-28: code below is all new, final fates but tracking individual fish
  # create a matrix to track fish
  # every fish starts in the start state
  indiv_fish_state_matrix <- matrix(data = start_state, nrow = nsim, ncol = 1)
  colnames(indiv_fish_state_matrix) <- "t1"
  
  # also track the direction of each fish
  # 1 = upstream, 2 = downstream, 0 = dead
  # every fish is going upstream to start
  indiv_fish_direction_matrix <- matrix(data = 1, nrow = nsim, ncol = 1)
  colnames(indiv_fish_direction_matrix) <- "t1"
  
  # use a while loop to keep running until all fish have entered loss state
  i <- 1
  condition <- "not finished"
  while (condition != "finished"){
    # if there are still fish moving, add a column to track their next movement
    next_time_step <- matrix(nrow = nsim, ncol = 1)
    colnames(next_time_step) <- paste0("t", i+1)
    
    indiv_fish_state_matrix <- cbind(indiv_fish_state_matrix, next_time_step)
    
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_moving_fish <- which(indiv_fish_direction_matrix[,i] == 1)
    if(sum(upstream_moving_fish)>0){
      for (j in upstream_moving_fish){
        indiv_fish_state_matrix[j,i+1] <- which(rmultinom(n = 1, size = 1, prob = upstream_move_prob_matrix[indiv_fish_state_matrix[j,i],]) == 1)
      }
    }
    
    
    downstream_moving_fish <- which(indiv_fish_direction_matrix[,i] == 2)
    if(sum(downstream_moving_fish)>0){
      for (j in downstream_moving_fish){
        indiv_fish_state_matrix[j,i+1] <- which(rmultinom(n = 1, size = 1, prob = downstream_move_prob_matrix[indiv_fish_state_matrix[j,i],]) == 1)
      }
    }
    
    # now, populate direction of movement for each fish
    indiv_fish_direction_matrix <- cbind(indiv_fish_direction_matrix, next_time_step)
    
    indiv_fish_direction_matrix[,i+1] <- ifelse(indiv_fish_state_matrix[,i+1] > indiv_fish_state_matrix[,i], 1, 2)
    
    # if it's in the loss state, direction = 0
    indiv_fish_direction_matrix[which(indiv_fish_state_matrix[,i+1] == 41),i+1] <- 0
    
    # if they're all in the loss state (direction = 0), then we're finish
    if (sum(indiv_fish_direction_matrix[,i+1], na.rm = T) == 0){
      condition <- "finished"
    }
    
    # increase counter
    i <- i + 1
  }
  
  return(indiv_fish_state_matrix)
  
}

final_fates_simulation_individual_SRH_T <- function(nsim,
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
      sample_date[i] <- sample(subset(SRH_T_states_dates, state == i)$date, 1)
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
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRH_T_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRH_T_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRH_T_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRH_T_states_dates, state == i  & direction == "downstream")$season)
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
        if(nrow(subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
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
        if(nrow(subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
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
    possible_movements <- SRH_T_envir$data$movements[, "col"][SRH_T_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_T[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRH_T[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRH_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    borigin1_array_SRH_T[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH_T[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH_T[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH_T[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH_T[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_T[i,possible_movements,iter] +
                    btemp1_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRH_T[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRH_T[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    borigin1_array_SRH_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_T[i,possible_movements,iter]*origin5))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_T[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRH_T[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRH_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    borigin1_array_SRH_T[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH_T[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH_T[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH_T[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH_T[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_T[i,possible_movements,iter] +
                    btemp0_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRH_T[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRH_T[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    borigin1_array_SRH_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_T[i,possible_movements,iter]*origin5))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRH_T_envir$data$movements[, "col"][SRH_T_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_T[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRH_T[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRH_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      borigin1_array_SRH_T[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH_T[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH_T[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH_T[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH_T[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_T[i,possible_movements,iter] +
                    btemp1_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRH_T[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRH_T[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    borigin1_array_SRH_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_T[i,possible_movements,iter]*origin5))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_T[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRH_T[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRH_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      borigin1_array_SRH_T[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH_T[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH_T[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH_T[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH_T[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_T[i,possible_movements,iter] +
                    btemp0_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRH_T[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRH_T[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    borigin1_array_SRH_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_T[i,possible_movements,iter]*origin5))
      }
      
      
      
    }
    
  }
  
  # 2025-07-28: code below is all new, final fates but tracking individual fish
  # create a matrix to track fish
  # every fish starts in the start state
  indiv_fish_state_matrix <- matrix(data = start_state, nrow = nsim, ncol = 1)
  colnames(indiv_fish_state_matrix) <- "t1"
  
  # also track the direction of each fish
  # 1 = upstream, 2 = downstream, 0 = dead
  # every fish is going upstream to start
  indiv_fish_direction_matrix <- matrix(data = 1, nrow = nsim, ncol = 1)
  colnames(indiv_fish_direction_matrix) <- "t1"
  
  # use a while loop to keep running until all fish have entered loss state
  i <- 1
  condition <- "not finished"
  while (condition != "finished"){
    # if there are still fish moving, add a column to track their next movement
    next_time_step <- matrix(nrow = nsim, ncol = 1)
    colnames(next_time_step) <- paste0("t", i+1)
    
    indiv_fish_state_matrix <- cbind(indiv_fish_state_matrix, next_time_step)
    
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_moving_fish <- which(indiv_fish_direction_matrix[,i] == 1)
    if(sum(upstream_moving_fish)>0){
      for (j in upstream_moving_fish){
        indiv_fish_state_matrix[j,i+1] <- which(rmultinom(n = 1, size = 1, prob = upstream_move_prob_matrix[indiv_fish_state_matrix[j,i],]) == 1)
      }
    }
    
    
    downstream_moving_fish <- which(indiv_fish_direction_matrix[,i] == 2)
    if(sum(downstream_moving_fish)>0){
      for (j in downstream_moving_fish){
        indiv_fish_state_matrix[j,i+1] <- which(rmultinom(n = 1, size = 1, prob = downstream_move_prob_matrix[indiv_fish_state_matrix[j,i],]) == 1)
      }
    }
    
    # now, populate direction of movement for each fish
    indiv_fish_direction_matrix <- cbind(indiv_fish_direction_matrix, next_time_step)
    
    indiv_fish_direction_matrix[,i+1] <- ifelse(indiv_fish_state_matrix[,i+1] > indiv_fish_state_matrix[,i], 1, 2)
    
    # if it's in the loss state, direction = 0
    indiv_fish_direction_matrix[which(indiv_fish_state_matrix[,i+1] == 41),i+1] <- 0
    
    # if they're all in the loss state (direction = 0), then we're finish
    if (sum(indiv_fish_direction_matrix[,i+1], na.rm = T) == 0){
      condition <- "finished"
    }
    
    # increase counter
    i <- i + 1
  }
  
  return(indiv_fish_state_matrix)
  
}

final_fates_simulation_individual_SRH_NT <- function(nsim,
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
      sample_date[i] <- sample(subset(SRH_NT_states_dates, state == i)$date, 1)
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
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRH_NT_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRH_NT_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRH_NT_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRH_NT_states_dates, state == i  & direction == "downstream")$season)
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
        if(nrow(subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
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
        if(nrow(subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
      }
      # season 1
      if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
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
    possible_movements <- SRH_NT_envir$data$movements[, "col"][SRH_NT_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_NT[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRH_NT[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRH_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    borigin1_array_SRH_NT[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH_NT[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH_NT[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH_NT[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH_NT[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_NT[i,possible_movements,iter] +
                    btemp1_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRH_NT[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRH_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    borigin1_array_SRH_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_NT[i,possible_movements,iter]*origin5))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_NT[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRH_NT[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRH_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    borigin1_array_SRH_NT[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH_NT[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH_NT[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH_NT[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH_NT[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_NT[i,possible_movements,iter] +
                    btemp0_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRH_NT[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRH_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    borigin1_array_SRH_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_NT[i,possible_movements,iter]*origin5))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRH_NT_envir$data$movements[, "col"][SRH_NT_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_NT[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRH_NT[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRH_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      borigin1_array_SRH_NT[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH_NT[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH_NT[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH_NT[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH_NT[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_NT[i,possible_movements,iter] +
                    btemp1_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRH_NT[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRH_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    borigin1_array_SRH_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_NT[i,possible_movements,iter]*origin5))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_NT[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRH_NT[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRH_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      borigin1_array_SRH_NT[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH_NT[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH_NT[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH_NT[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH_NT[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_NT[i,possible_movements,iter] +
                    btemp0_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRH_NT[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRH_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    borigin1_array_SRH_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_NT[i,possible_movements,iter]*origin5))
      }
      
      
      
    }
    
  }
  
  # 2025-07-28: code below is all new, final fates but tracking individual fish
  # create a matrix to track fish
  # every fish starts in the start state
  indiv_fish_state_matrix <- matrix(data = start_state, nrow = nsim, ncol = 1)
  colnames(indiv_fish_state_matrix) <- "t1"
  
  # also track the direction of each fish
  # 1 = upstream, 2 = downstream, 0 = dead
  # every fish is going upstream to start
  indiv_fish_direction_matrix <- matrix(data = 1, nrow = nsim, ncol = 1)
  colnames(indiv_fish_direction_matrix) <- "t1"
  
  # use a while loop to keep running until all fish have entered loss state
  i <- 1
  condition <- "not finished"
  while (condition != "finished"){
    # if there are still fish moving, add a column to track their next movement
    next_time_step <- matrix(nrow = nsim, ncol = 1)
    colnames(next_time_step) <- paste0("t", i+1)
    
    indiv_fish_state_matrix <- cbind(indiv_fish_state_matrix, next_time_step)
    
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_moving_fish <- which(indiv_fish_direction_matrix[,i] == 1)
    if(sum(upstream_moving_fish)>0){
      for (j in upstream_moving_fish){
        indiv_fish_state_matrix[j,i+1] <- which(rmultinom(n = 1, size = 1, prob = upstream_move_prob_matrix[indiv_fish_state_matrix[j,i],]) == 1)
      }
    }
    
    
    downstream_moving_fish <- which(indiv_fish_direction_matrix[,i] == 2)
    if(sum(downstream_moving_fish)>0){
      for (j in downstream_moving_fish){
        indiv_fish_state_matrix[j,i+1] <- which(rmultinom(n = 1, size = 1, prob = downstream_move_prob_matrix[indiv_fish_state_matrix[j,i],]) == 1)
      }
    }
    
    # now, populate direction of movement for each fish
    indiv_fish_direction_matrix <- cbind(indiv_fish_direction_matrix, next_time_step)
    
    indiv_fish_direction_matrix[,i+1] <- ifelse(indiv_fish_state_matrix[,i+1] > indiv_fish_state_matrix[,i], 1, 2)
    
    # if it's in the loss state, direction = 0
    indiv_fish_direction_matrix[which(indiv_fish_state_matrix[,i+1] == 41),i+1] <- 0
    
    # if they're all in the loss state (direction = 0), then we're finish
    if (sum(indiv_fish_direction_matrix[,i+1], na.rm = T) == 0){
      condition <- "finished"
    }
    
    # increase counter
    i <- i + 1
  }
  
  return(indiv_fish_state_matrix)
  
}





#### Run final fates simulation - H v W, T vs NT comparisons ####


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
SRW_T_states_dates <- get_states_dates_direction(envir = SRW_T_envir)
SRW_NT_states_dates <- get_states_dates_direction(envir = SRW_NT_envir)
SRH_T_states_dates <- get_states_dates_direction(envir = SRH_T_envir)
SRH_NT_states_dates <- get_states_dates_direction(envir = SRH_NT_envir)

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

compare_sim_movement_histories_MC <- function(niter, nsim, condition_jitter,
                                              origin_select){
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    
    
    sim_history_wild <- data.frame()
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    
    # index to the right origin param to turn it on
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_individual_MCW(nsim = nsim,
                                                        start_state = 2, states_dates = get_origin_states_dates(envir = MCW_envir, origin_select = origin_select, rear = "wild"),
                                                        origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                        origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                        temp_data = MCW_envir$data$temperature_data, spillwindow_data = MCW_envir$data$spill_window_data, 
                                                        winterspill_data = MCW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      sim_history_wild %>% 
        bind_rows(as.data.frame(sim_wild)) -> sim_history_wild
    }
    
    sim_history_wild$iter <- rep(1:niter, each = nsim)
    
    return(list(sim_history_wild, NULL))
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    
    
    sim_history_hatchery <- data.frame()
    
    # set up vectors to control which origin parameters are turned on
    hatchery_origin_params <- rep(0,2)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_individual_MCH(nsim = nsim,
                                                            start_state = 2, states_dates = get_origin_states_dates(envir = MCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                            origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                            temp_data = MCH_envir$data$temperature_data, spillwindow_data = MCH_envir$data$spill_window_data, 
                                                            winterspill_data = MCH_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      sim_history_hatchery %>% 
        bind_rows(as.data.frame(sim_hatchery)) -> sim_history_hatchery
    }
    
    sim_history_hatchery$iter <- rep(1:niter, each = nsim)
    
    return(list(NULL, sim_history_hatchery))
    
  } 
  # else run both
  else {
    
    
    
    sim_history_hatchery <- data.frame()
    sim_history_wild <- data.frame()
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    hatchery_origin_params <- rep(0,2)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_individual_MCW(nsim = nsim,
                                                        start_state = 2, states_dates = get_origin_states_dates(envir = MCW_envir, origin_select = origin_select, rear = "wild"),
                                                        origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                        origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                        temp_data = MCW_envir$data$temperature_data, spillwindow_data = MCW_envir$data$spill_window_data, 
                                                        winterspill_data = MCW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      sim_history_wild %>% 
        bind_rows(as.data.frame(sim_wild)) -> sim_history_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_individual_MCH(nsim = nsim,
                                                            start_state = 2, states_dates = get_origin_states_dates(envir = MCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                            origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                            temp_data = MCH_envir$data$temperature_data, spillwindow_data = MCH_envir$data$spill_window_data, 
                                                            winterspill_data = MCH_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      sim_history_hatchery %>% 
        bind_rows(as.data.frame(sim_hatchery)) -> sim_history_hatchery
    }
    sim_history_wild$iter <- rep(1:niter, each = nsim)
    sim_history_hatchery$iter <- rep(1:niter, each = nsim)
    
    return(list(sim_history_wild, sim_history_hatchery))
    
  }
  
}

compare_sim_movement_histories_UC <- function(niter, nsim, condition_jitter,
                                              origin_select){
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    
    
    sim_history_wild <- data.frame()
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    
    # index to the right origin param to turn it on
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_individual_UCW(nsim = nsim,
                                                        start_state = 2, states_dates = get_origin_states_dates(envir = UCW_envir, origin_select = origin_select, rear = "wild"),
                                                        origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                        temp_data = UCW_envir$data$temperature_data, spillwindow_data = UCW_envir$data$spill_window_data, 
                                                        winterspill_data = UCW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      sim_history_wild %>% 
        bind_rows(as.data.frame(sim_wild)) -> sim_history_wild
    }
    
    sim_history_wild$iter <- rep(1:niter, each = nsim)
    
    return(list(sim_history_wild, NULL))
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    
    
    sim_history_hatchery <- data.frame()
    
    # set up vectors to control which origin parameters are turned on
    hatchery_origin_params <- rep(0,5)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_individual_UCH(nsim = nsim,
                                                            start_state = 2, states_dates = get_origin_states_dates(envir = UCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                            origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                            origin3 = hatchery_origin_params[3], 
                                                            temp_data = UCH_envir$data$temperature_data, spillwindow_data = UCH_envir$data$spill_window_data, 
                                                            winterspill_data = UCH_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      sim_history_hatchery %>% 
        bind_rows(as.data.frame(sim_hatchery)) -> sim_history_hatchery
    }
    
    sim_history_hatchery$iter <- rep(1:niter, each = nsim)
    
    return(list(NULL, sim_history_hatchery))
    
  } 
  # else run both
  else {
    
    
    
    sim_history_hatchery <- data.frame()
    sim_history_wild <- data.frame()
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    hatchery_origin_params <- rep(0,5)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_individual_UCW(nsim = nsim,
                                                        start_state = 2, states_dates = get_origin_states_dates(envir = UCW_envir, origin_select = origin_select, rear = "wild"),
                                                        origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                        temp_data = UCW_envir$data$temperature_data, spillwindow_data = UCW_envir$data$spill_window_data, 
                                                        winterspill_data = UCW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      sim_history_wild %>% 
        bind_rows(as.data.frame(sim_wild)) -> sim_history_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_individual_UCH(nsim = nsim,
                                                            start_state = 2, states_dates = get_origin_states_dates(envir = UCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                            origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                            origin3 = hatchery_origin_params[3], 
                                                            temp_data = UCH_envir$data$temperature_data, spillwindow_data = UCH_envir$data$spill_window_data, 
                                                            winterspill_data = UCH_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      sim_history_hatchery %>% 
        bind_rows(as.data.frame(sim_hatchery)) -> sim_history_hatchery
    }
    
    sim_history_wild$iter <- rep(1:niter, each = nsim)
    sim_history_hatchery$iter <- rep(1:niter, each = nsim)
    
    return(list(sim_history_wild, sim_history_hatchery))
    
  }
  
}

compare_sim_movement_histories_SR_T <- function(niter, nsim, condition_jitter,
                                            origin_select){
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    
    
    sim_history_wild <- data.frame()
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    
    # index to the right origin param to turn it on
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_individual_SRW_T(nsim = nsim,
                                             start_state = 2, states_dates = get_origin_states_dates(envir = SRW_T_envir, origin_select = origin_select, rear = "wild"),
                                             origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                             origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                             temp_data = SRW_T_envir$data$temperature_data, spillwindow_data = SRW_T_envir$data$spill_window_data, 
                                             winterspill_data = SRW_T_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      sim_history_wild %>% 
        bind_rows(as.data.frame(sim_wild)) -> sim_history_wild
    }
    
    sim_history_wild$iter <- rep(1:niter, each = nsim)
    
    return(list(sim_history_wild, NULL))
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    
    
    sim_history_hatchery <- data.frame()
    
    # set up vectors to control which origin parameters are turned on
    hatchery_origin_params <- rep(0,5)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_individual_SRH_T(nsim = nsim,
                                                 start_state = 2, states_dates = get_origin_states_dates(envir = SRH_T_envir, origin_select = origin_select, rear = "hatchery"),
                                                 origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                 origin3 = hatchery_origin_params[3], origin4 = hatchery_origin_params[4], 
                                                 origin5 = hatchery_origin_params[5],
                                                 temp_data = SRH_T_envir$data$temperature_data, spillwindow_data = SRH_T_envir$data$spill_window_data, 
                                                 winterspill_data = SRH_T_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      sim_history_hatchery %>% 
        bind_rows(as.data.frame(sim_hatchery)) -> sim_history_hatchery
    }
    
    sim_history_hatchery$iter <- rep(1:niter, each = nsim)
    
    return(list(NULL, sim_history_hatchery))
    
  } 
  # else run both
  else {
    
    
    
    sim_history_hatchery <- data.frame()
    sim_history_wild <- data.frame()
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    hatchery_origin_params <- rep(0,5)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_individual_SRW_T(nsim = nsim,
                                                          start_state = 2, states_dates = get_origin_states_dates(envir = SRW_T_envir, origin_select = origin_select, rear = "wild"),
                                                          origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                          origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                          temp_data = SRW_T_envir$data$temperature_data, spillwindow_data = SRW_T_envir$data$spill_window_data, 
                                                          winterspill_data = SRW_T_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      sim_history_wild %>% 
        bind_rows(as.data.frame(sim_wild)) -> sim_history_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_individual_SRH_T(nsim = nsim,
                                                              start_state = 2, states_dates = get_origin_states_dates(envir = SRH_T_envir, origin_select = origin_select, rear = "hatchery"),
                                                              origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                              origin3 = hatchery_origin_params[3], origin4 = hatchery_origin_params[4], 
                                                              origin5 = hatchery_origin_params[5],
                                                              temp_data = SRH_T_envir$data$temperature_data, spillwindow_data = SRH_T_envir$data$spill_window_data, 
                                                              winterspill_data = SRH_T_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      sim_history_hatchery %>% 
        bind_rows(as.data.frame(sim_hatchery)) -> sim_history_hatchery
    }
    
    sim_history_wild$iter <- rep(1:niter, each = nsim)
    sim_history_hatchery$iter <- rep(1:niter, each = nsim)
    
    return(list(sim_history_wild, sim_history_hatchery))

  }
  
}

compare_sim_movement_histories_SR_NT <- function(niter, nsim, condition_jitter,
                                                 origin_select){
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    
    
    sim_history_wild <- data.frame()
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    
    # index to the right origin param to turn it on
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_individual_SRW_NT(nsim = nsim,
                                                           start_state = 2, states_dates = get_origin_states_dates(envir = SRW_NT_envir, origin_select = origin_select, rear = "wild"),
                                                           origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                           origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                           temp_data = SRW_NT_envir$data$temperature_data, spillwindow_data = SRW_NT_envir$data$spill_window_data, 
                                                           winterspill_data = SRW_NT_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      sim_history_wild %>% 
        bind_rows(as.data.frame(sim_wild)) -> sim_history_wild
    }
    return(list(sim_history_wild, NULL))
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    
    
    sim_history_hatchery <- data.frame()
    
    # set up vectors to control which origin parameters are turned on
    hatchery_origin_params <- rep(0,5)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_individual_SRH_NT(nsim = nsim,
                                                               start_state = 2, states_dates = get_origin_states_dates(envir = SRH_NT_envir, origin_select = origin_select, rear = "hatchery"),
                                                               origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                               origin3 = hatchery_origin_params[3], origin4 = hatchery_origin_params[4], 
                                                               origin5 = hatchery_origin_params[5],
                                                               temp_data = SRH_NT_envir$data$temperature_data, spillwindow_data = SRH_NT_envir$data$spill_window_data, 
                                                               winterspill_data = SRH_NT_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      sim_history_hatchery %>% 
        bind_rows(as.data.frame(sim_hatchery)) -> sim_history_hatchery
    }
    
    return(list(NULL, sim_history_hatchery))
    
  } 
  # else run both
  else {
    
    
    
    sim_history_hatchery <- data.frame()
    sim_history_wild <- data.frame()
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    hatchery_origin_params <- rep(0,5)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_individual_SRW_NT(nsim = nsim,
                                                           start_state = 2, states_dates = get_origin_states_dates(envir = SRW_NT_envir, origin_select = origin_select, rear = "wild"),
                                                           origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                           origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                           temp_data = SRW_NT_envir$data$temperature_data, spillwindow_data = SRW_NT_envir$data$spill_window_data, 
                                                           winterspill_data = SRW_NT_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      sim_history_wild %>% 
        bind_rows(as.data.frame(sim_wild)) -> sim_history_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_individual_SRH_NT(nsim = nsim,
                                                               start_state = 2, states_dates = get_origin_states_dates(envir = SRH_NT_envir, origin_select = origin_select, rear = "hatchery"),
                                                               origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                               origin3 = hatchery_origin_params[3], origin4 = hatchery_origin_params[4], 
                                                               origin5 = hatchery_origin_params[5],
                                                               temp_data = SRH_NT_envir$data$temperature_data, spillwindow_data = SRH_NT_envir$data$spill_window_data, 
                                                               winterspill_data = SRH_NT_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      sim_history_hatchery %>% 
        bind_rows(as.data.frame(sim_hatchery)) -> sim_history_hatchery
    }
    
    
    return(list(sim_history_wild, sim_history_hatchery))
    
  }
  
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



ff_iter <- 1000
ff_nsim <- 1000


#### run individual tracking simulation ####

jdr_rear_sim_histories <- compare_sim_movement_histories_MC(niter = ff_iter, 
                                                                          nsim = ff_nsim, 
                                                                          origin_select = "John Day River",
                                                                          condition_jitter = FALSE)
save(jdr_rear_sim_histories, file = here::here("individual_simulation_runs", "jdr_rear_sim_histories.rda"))

uma_rear_sim_histories <- compare_sim_movement_histories_MC(niter = ff_iter, 
                                                                        nsim = ff_nsim, 
                                                                        origin_select = "Umatilla River",
                                                                        condition_jitter = FALSE)
save(uma_rear_sim_histories, file = here::here("individual_simulation_runs", "uma_rear_sim_histories.rda"))

wawa_rear_sim_histories <- compare_sim_movement_histories_MC(niter = ff_iter, 
                                                                        nsim = ff_nsim, 
                                                                        origin_select = "Walla Walla River",
                                                                        condition_jitter = FALSE)
save(wawa_rear_sim_histories, file = here::here("individual_simulation_runs", "wawa_rear_sim_histories.rda"))

yak_rear_sim_histories <- compare_sim_movement_histories_MC(niter = ff_iter, 
                                                                        nsim = ff_nsim, 
                                                                        origin_select = "Yakima River",
                                                                        condition_jitter = FALSE)
save(yak_rear_sim_histories, file = here::here("individual_simulation_runs", "yak_rear_sim_histories.rda"))

wen_rear_sim_histories <- compare_sim_movement_histories_UC(niter = ff_iter, 
                                                                        nsim = ff_nsim, 
                                                                        origin_select = "Wenatchee River",
                                                                        condition_jitter = FALSE)
save(wen_rear_sim_histories, file = here::here("individual_simulation_runs", "wen_rear_sim_histories.rda"))

ent_rear_sim_histories <- compare_sim_movement_histories_UC(niter = ff_iter, 
                                                                        nsim = ff_nsim, 
                                                                        origin_select = "Entiat River",
                                                                        condition_jitter = FALSE)
save(ent_rear_sim_histories, file = here::here("individual_simulation_runs", "ent_rear_sim_histories.rda"))

tuc_transported_rear_sim_histories <- compare_sim_movement_histories_SR_T(niter = ff_iter, 
                                    nsim = ff_nsim, 
                                    origin_select = "Tucannon River",
                                    condition_jitter = FALSE)
save(tuc_transported_rear_sim_histories, file = here::here("individual_simulation_runs", "tuc_transported_rear_sim_histories.rda"))


tuc_not_transported_rear_sim_histories <- compare_sim_movement_histories_SR_NT(niter = ff_iter, 
                                                                          nsim = ff_nsim, 
                                                                          origin_select = "Tucannon River",
                                                                          condition_jitter = FALSE)
save(tuc_not_transported_rear_sim_histories, file = here::here("individual_simulation_runs", "tuc_not_transported_rear_sim_histories.rda"))


# load files
load(file = here::here("individual_simulation_runs", "jdr_rear_sim_histories.rda"))
load(file = here::here("individual_simulation_runs", "uma_rear_sim_histories.rda"))
load(file = here::here("individual_simulation_runs", "wawa_rear_sim_histories.rda"))
load(file = here::here("individual_simulation_runs", "yak_rear_sim_histories.rda"))
load(file = here::here("individual_simulation_runs", "wen_rear_sim_histories.rda"))
load(file = here::here("individual_simulation_runs", "ent_rear_sim_histories.rda"))
load(file = here::here("individual_simulation_runs", "tuc_transported_rear_sim_histories.rda"))
load(file = here::here("individual_simulation_runs", "tuc_not_transported_rear_sim_histories.rda"))

# # temp fix
# jdr_rear_sim_histories[[1]]$iter <- rep(1:ff_iter, each = ff_nsim)
# uma_rear_sim_histories[[1]]$iter <- rep(1:ff_iter, each = ff_nsim)
# uma_rear_sim_histories[[2]]$iter <- rep(1:ff_iter, each = ff_nsim)
# wawa_rear_sim_histories[[1]]$iter <- rep(1:ff_iter, each = ff_nsim)
# wawa_rear_sim_histories[[2]]$iter <- rep(1:ff_iter, each = ff_nsim)
# yak_rear_sim_histories[[1]]$iter <- rep(1:ff_iter, each = ff_nsim)
wen_rear_sim_histories[[1]]$iter <- rep(1:ff_iter, each = ff_nsim)
wen_rear_sim_histories[[2]]$iter <- rep(1:ff_iter, each = ff_nsim)
ent_rear_sim_histories[[1]]$iter <- rep(1:ff_iter, each = ff_nsim)
tuc_transported_rear_sim_histories[[1]]$iter <- rep(1:ff_iter, each = ff_nsim)
tuc_not_transported_rear_sim_histories[[1]]$iter <- rep(1:ff_iter, each = ff_nsim)
tuc_transported_rear_sim_histories[[2]]$iter <- rep(1:ff_iter, each = ff_nsim)
tuc_not_transported_rear_sim_histories[[2]]$iter <- rep(1:ff_iter, each = ff_nsim)



#### probability of overshooting multiple dams ####
multiple_overshoot_summary <- function(rear_sim_history, natal_origin_abbrev){
  
  # if it ever happens that a movement history has over 100 movements, drop it on the 
  # grounds that that's a biologically unrealistic movement history
  if (ncol(rear_sim_history) > 100){
    rear_sim_history[-which(!(is.na(rear_sim_history[,101]))),] -> rear_sim_history
    rear_sim_history <- rear_sim_history[,c(1:100, ncol(rear_sim_history))]
  }
  
  
  
  
  # note which states are overshoot dams
  overshoot_states_list <- list("jdr" = model_states[3:9],
                                "uma" = model_states[3:9],
                                "yak" = model_states[4:9],
                                "wawa" = model_states[4:9],
                                "wawa" = model_states[4:9],
                                "wen" = model_states[6:7],
                                "ent" = model_states[7],
                                "tuc" = model_states[9])
  
  overshoot_states <- overshoot_states_list[[natal_origin_abbrev]]
  
  # reformat input histories
  as.data.frame(rear_sim_history) %>% 
    rownames_to_column("sim_fish_ID") %>% 
    mutate(sim_fish_ID = as.numeric(sim_fish_ID)) %>% 
    pivot_longer(cols = starts_with("t"), names_to = "t", values_to = "state_number") %>% 
    filter(!(is.na(state_number))) -> sim_history
  
  model_states_df <- data.frame(state = model_states, state_number = 1:length(model_states))
  
  sim_history %>% 
    left_join(model_states_df, by = "state_number") %>% 
    mutate(state = gsub(" Upstream", "", state)) %>% 
    mutate(state = gsub(" Mouth", "", state)) -> sim_history
  
  # count number of unique overshoot states (this is overshooting multiple dams)
  sim_history %>% 
    filter(state %in% overshoot_states) %>% 
    group_by(sim_fish_ID) %>% 
    summarise(n_dams_overshot = n_distinct(state)) -> n_overshoot_table
  
  all_fish <- data.frame(sim_fish_ID = as.numeric(rownames(rear_sim_history)), iter = rear_sim_history$iter)
  
  all_fish %>% 
    left_join(n_overshoot_table, by = "sim_fish_ID") %>% 
    mutate(n_dams_overshot = ifelse(is.na(n_dams_overshot), 0, n_dams_overshot))-> all_fish
  
  total_fish <- length(unique(all_fish$sim_fish_ID))
  
  all_fish %>% 
    mutate(dam_overshoots = as.character(ifelse(n_dams_overshot >= 2, "2+", n_dams_overshot))) %>% 
    group_by(iter) %>% 
    count(dam_overshoots) %>% 
    mutate(prop = n/ff_nsim) %>% 
    ungroup() %>% 
    group_by(dam_overshoots) %>% 
    summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) -> overshoot_table
    
  
  # as.data.frame(table(all_fish$n_dams_overshot)) %>% 
  #   dplyr::rename(n_overshoot_dams = Var1, nfish = Freq) %>% 
  #   mutate(prop = nfish/total_fish) -> overshoot_table
  
  
  return(overshoot_table)
}


jdr_wild_overshoot_summary <- multiple_overshoot_summary(rear_sim_history = jdr_rear_sim_histories[[1]], natal_origin_abbrev = "jdr")
uma_wild_overshoot_summary <- multiple_overshoot_summary(rear_sim_history = uma_rear_sim_histories[[1]], natal_origin_abbrev = "uma")
uma_hatchery_overshoot_summary <- multiple_overshoot_summary(rear_sim_history = uma_rear_sim_histories[[2]], natal_origin_abbrev = "uma")
wawa_wild_overshoot_summary <- multiple_overshoot_summary(rear_sim_history = wawa_rear_sim_histories[[1]], natal_origin_abbrev = "wawa")
wawa_hatchery_overshoot_summary <- multiple_overshoot_summary(rear_sim_history = wawa_rear_sim_histories[[2]], natal_origin_abbrev = "wawa")
yak_wild_overshoot_summary <- multiple_overshoot_summary(rear_sim_history = yak_rear_sim_histories[[1]], natal_origin_abbrev = "yak")
wen_wild_overshoot_summary <- multiple_overshoot_summary(rear_sim_history = wen_rear_sim_histories[[1]], natal_origin_abbrev = "wen")
wen_hatchery_overshoot_summary <- multiple_overshoot_summary(rear_sim_history = wen_rear_sim_histories[[2]], natal_origin_abbrev = "wen")
ent_wild_overshoot_summary <- multiple_overshoot_summary(rear_sim_history = ent_rear_sim_histories[[1]], natal_origin_abbrev = "ent")
tuc_wild_transported_overshoot_summary <- multiple_overshoot_summary(rear_sim_history = tuc_transported_rear_sim_histories[[1]], natal_origin_abbrev = "tuc")
tuc_wild_not_transported_overshoot_summary <- multiple_overshoot_summary(rear_sim_history = tuc_not_transported_rear_sim_histories[[1]], natal_origin_abbrev = "tuc")
tuc_hatchery_transported_overshoot_summary <- multiple_overshoot_summary(rear_sim_history = tuc_transported_rear_sim_histories[[2]], natal_origin_abbrev = "tuc")
tuc_hatchery_not_transported_overshoot_summary <- multiple_overshoot_summary(rear_sim_history = tuc_not_transported_rear_sim_histories[[2]], natal_origin_abbrev = "tuc")

# make a summary table - join these outputs

mutate(jdr_wild_overshoot_summary, origin_rear = "John Day River, Natural") %>% 
  bind_rows(., mutate(uma_wild_overshoot_summary, origin_rear = "Umatilla River, Natural")) %>% 
  bind_rows(., mutate(uma_hatchery_overshoot_summary, origin_rear = "Umatilla River, Hatchery")) %>% 
  bind_rows(., mutate(wawa_wild_overshoot_summary, origin_rear = "Walla Walla River, Natural")) %>% 
  bind_rows(., mutate(wawa_hatchery_overshoot_summary, origin_rear = "Walla Walla River, Hatchery")) %>% 
  bind_rows(., mutate(yak_wild_overshoot_summary, origin_rear = "Yakima River, Natural")) %>% 
  bind_rows(., mutate(wen_wild_overshoot_summary, origin_rear = "Wenatchee River, Natural")) %>% 
  bind_rows(., mutate(wen_hatchery_overshoot_summary, origin_rear = "Wenatchee River, Hatchery")) %>% 
  bind_rows(., mutate(ent_wild_overshoot_summary, origin_rear = "Entiat River, Natural")) %>% 
  bind_rows(., mutate(tuc_wild_transported_overshoot_summary, origin_rear = "Tucannon River, Natural (transported)")) %>% 
  bind_rows(., mutate(tuc_wild_not_transported_overshoot_summary, origin_rear = "Tucannon River, Natural (not transported)")) %>% 
  bind_rows(., mutate(tuc_hatchery_transported_overshoot_summary, origin_rear = "Tucannon River, Hatchery (transported)")) %>% 
  bind_rows(., mutate(tuc_hatchery_not_transported_overshoot_summary, origin_rear = "Tucannon River, Hatchery (not transported)")) %>% 
  ungroup() %>% 
  complete(dam_overshoots, origin_rear, fill = list(NA)) -> all_overshoot_summary
  
# reformat as a table
all_overshoot_summary %>% 
  pivot_wider(names_from = q, values_from = prob) %>% 
  mutate(prob_for_table = paste0(round(`0.5`, 2), " (", round(`0.025`, 2), " - ", round(`0.975`, 2), ")")) %>% 
  mutate(prob_for_table = ifelse(prob_for_table == "NA (NA - NA)", "-", prob_for_table)) %>% 
  dplyr::select(dam_overshoots, origin_rear, prob_for_table) %>% 
  pivot_wider(names_from = dam_overshoots, values_from = prob_for_table) -> all_overshoot_summary_table
  
write.csv(all_overshoot_summary_table, file = here::here("figures", "transport", "final_fates", "all_overshoot_summary_table.csv"))


#### probabilty of homing, conditional on multiple overshoots ####

overshoots_homing_summary <- function(rear_sim_history, natal_origin_abbrev, natal_origin){
  # if it ever happens that a movement history has over 100 movements, drop it on the 
  # grounds that that's a biologically unrealistic movement history
  if (ncol(rear_sim_history) > 100){
    rear_sim_history[-which(!(is.na(rear_sim_history[,101]))),] -> rear_sim_history
    rear_sim_history <- rear_sim_history[,c(1:100, ncol(rear_sim_history))]
  }
  
  # note which states are overshoot dams
  overshoot_states_list <- list("jdr" = model_states[3:9],
                                "uma" = model_states[3:9],
                                "yak" = model_states[4:9],
                                "wawa" = model_states[4:9],
                                "wawa" = model_states[4:9],
                                "wen" = model_states[6:7],
                                "ent" = model_states[7],
                                "tuc" = model_states[9])
  
  overshoot_states <- overshoot_states_list[[natal_origin_abbrev]]
  
  # reformat input histories
  as.data.frame(rear_sim_history) %>% 
    rownames_to_column("sim_fish_ID") %>% 
    mutate(sim_fish_ID = as.numeric(sim_fish_ID)) %>% 
    pivot_longer(cols = starts_with("t"), names_to = "t", values_to = "state_number") %>% 
    filter(!(is.na(state_number))) -> sim_history
  
  model_states_df <- data.frame(state = model_states, state_number = 1:length(model_states))
  
  sim_history %>% 
    left_join(model_states_df, by = "state_number") %>% 
    mutate(state = gsub(" Upstream", "", state)) %>% 
    mutate(state = gsub(" Mouth", "", state)) -> sim_history
  
  # count number of unique overshoot states (this is overshooting multiple dams)
  sim_history %>% 
    filter(state %in% overshoot_states) %>% 
    group_by(sim_fish_ID) %>% 
    summarise(n_dams_overshot = n_distinct(state)) -> n_overshoot_table
  
  # determine which fish successfully homed
  sim_history %>% 
    filter(state == natal_origin) %>% 
    filter(!(duplicated(sim_fish_ID))) %>% 
    mutate(home = 1) %>% 
    dplyr::select(sim_fish_ID, home) -> homing_fish
  
  
  all_fish <- data.frame(sim_fish_ID = as.numeric(rownames(rear_sim_history)), iter = rear_sim_history$iter)
  
  all_fish %>% 
    left_join(homing_fish, by = "sim_fish_ID") %>% 
    left_join(n_overshoot_table, by = "sim_fish_ID") %>% 
    mutate(n_dams_overshot = ifelse(is.na(n_dams_overshot), 0, n_dams_overshot),
           home = ifelse(is.na(home), 0, home))-> all_fish
  
  total_fish <- length(unique(all_fish$sim_fish_ID))
  
  all_fish %>% 
    mutate(dam_overshoots = as.character(ifelse(n_dams_overshot >= 2, "2+", n_dams_overshot))) %>% 
    dplyr::select(-n_dams_overshot) %>%
    #
    group_by(dam_overshoots, iter) %>%
    # complete(dam_overshoots, iter, fill = list(NA)) %>% 
    count(home) %>% 
    group_by(dam_overshoots, iter) %>% 
    mutate(n_per_iter_novershoots = sum(n)) %>% 
    mutate(prop = n/n_per_iter_novershoots) %>% 
    ungroup() %>% 
    group_by(dam_overshoots) %>% 
    filter(home == 1) %>% 
    summarise(prob_of_homing = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) -> homing_n_overshoots_table
  
  
  # as.data.frame(table(all_fish$n_dams_overshot)) %>% 
  #   dplyr::rename(n_overshoot_dams = Var1, nfish = Freq) %>% 
  #   mutate(prop = nfish/total_fish) -> overshoot_table
  
  
  return(homing_n_overshoots_table)
}


jdr_wild_overshoots_homing_summary <- overshoots_homing_summary(rear_sim_history = jdr_rear_sim_histories[[1]], natal_origin_abbrev = "jdr",
                                                        natal_origin = "John Day River")
uma_wild_overshoots_homing_summary <- overshoots_homing_summary(rear_sim_history = uma_rear_sim_histories[[1]], natal_origin_abbrev = "uma",
                                                        natal_origin = "Umatilla River")
uma_hatchery_overshoots_homing_summary <- overshoots_homing_summary(rear_sim_history = uma_rear_sim_histories[[2]], natal_origin_abbrev = "uma",
                                                            natal_origin = "Umatilla River")
wawa_wild_overshoots_homing_summary <- overshoots_homing_summary(rear_sim_history = wawa_rear_sim_histories[[1]], natal_origin_abbrev = "wawa",
                                                         natal_origin = "Walla Walla River")
wawa_hatchery_overshoots_homing_summary <- overshoots_homing_summary(rear_sim_history = wawa_rear_sim_histories[[2]], natal_origin_abbrev = "wawa",
                                                             natal_origin = "Walla Walla River")
yak_wild_overshoots_homing_summary <- overshoots_homing_summary(rear_sim_history = yak_rear_sim_histories[[1]], natal_origin_abbrev = "yak",
                                                        natal_origin = "Yakima River")
wen_wild_overshoots_homing_summary <- overshoots_homing_summary(rear_sim_history = wen_rear_sim_histories[[1]], natal_origin_abbrev = "wen",
                                                        natal_origin = "Wenatchee River")
wen_hatchery_overshoots_homing_summary <- overshoots_homing_summary(rear_sim_history = wen_rear_sim_histories[[2]], natal_origin_abbrev = "wen",
                                                            natal_origin = "Wenatchee River")
ent_wild_overshoots_homing_summary <- overshoots_homing_summary(rear_sim_history = ent_rear_sim_histories[[1]], natal_origin_abbrev = "ent",
                                                        natal_origin = "Entiat River")
tuc_wild_transported_overshoots_homing_summary <- overshoots_homing_summary(rear_sim_history = tuc_transported_rear_sim_histories[[1]], natal_origin_abbrev = "tuc",
                                                                    natal_origin = "Tucannon River")
tuc_wild_not_transported_overshoots_homing_summary <- overshoots_homing_summary(rear_sim_history = tuc_not_transported_rear_sim_histories[[1]], 
                                                                        natal_origin_abbrev = "tuc", natal_origin = "Tucannon River")
tuc_hatchery_transported_overshoots_homing_summary <- overshoots_homing_summary(rear_sim_history = tuc_transported_rear_sim_histories[[2]], 
                                                                        natal_origin_abbrev = "tuc", natal_origin = "Tucannon River")
tuc_hatchery_not_transported_overshoots_homing_summary <- overshoots_homing_summary(rear_sim_history = tuc_not_transported_rear_sim_histories[[2]], 
                                                                            natal_origin_abbrev = "tuc", natal_origin = "Tucannon River")

# make a summary table - join these outputs

mutate(jdr_wild_overshoots_homing_summary, origin_rear = "John Day River, Natural") %>% 
  bind_rows(., mutate(uma_wild_overshoots_homing_summary, origin_rear = "Umatilla River, Natural")) %>% 
  bind_rows(., mutate(uma_hatchery_overshoots_homing_summary, origin_rear = "Umatilla River, Hatchery")) %>% 
  bind_rows(., mutate(wawa_wild_overshoots_homing_summary, origin_rear = "Walla Walla River, Natural")) %>% 
  bind_rows(., mutate(wawa_hatchery_overshoots_homing_summary, origin_rear = "Walla Walla River, Hatchery")) %>% 
  bind_rows(., mutate(yak_wild_overshoots_homing_summary, origin_rear = "Yakima River, Natural")) %>% 
  bind_rows(., mutate(wen_wild_overshoots_homing_summary, origin_rear = "Wenatchee River, Natural")) %>% 
  bind_rows(., mutate(wen_hatchery_overshoots_homing_summary, origin_rear = "Wenatchee River, Hatchery")) %>% 
  bind_rows(., mutate(ent_wild_overshoots_homing_summary, origin_rear = "Entiat River, Natural")) %>% 
  bind_rows(., mutate(tuc_wild_transported_overshoots_homing_summary, origin_rear = "Tucannon River, Natural (transported)")) %>% 
  bind_rows(., mutate(tuc_wild_not_transported_overshoots_homing_summary, origin_rear = "Tucannon River, Natural (not transported)")) %>% 
  bind_rows(., mutate(tuc_hatchery_transported_overshoots_homing_summary, origin_rear = "Tucannon River, Hatchery (transported)")) %>% 
  bind_rows(., mutate(tuc_hatchery_not_transported_overshoots_homing_summary, origin_rear = "Tucannon River, Hatchery (not transported)")) %>% 
  ungroup() %>% 
  complete(dam_overshoots, origin_rear, fill = list(NA)) -> all_overshoots_homing_summary

# reformat as a table
all_overshoots_homing_summary %>% 
  pivot_wider(names_from = q, values_from = prob_of_homing) %>% 
  mutate(prob_for_table = paste0(round(`0.5`, 2), " (", round(`0.025`, 2), " - ", round(`0.975`, 2), ")")) %>% 
  mutate(prob_for_table = ifelse(prob_for_table == "NA (NA - NA)", "-", prob_for_table)) %>% 
  dplyr::select(dam_overshoots, origin_rear, prob_for_table) %>% 
  pivot_wider(names_from = dam_overshoots, values_from = prob_for_table) -> all_overshoots_homing_summary_table

write.csv(all_overshoots_homing_summary_table, file = here::here("figures", "transport", "final_fates", "all_overshoots_homing_summary_table.csv"))


#### probabilty of homing, conditional on multiple overshoots AND having made it to the right mainstem state ####

overshoots_homing_conditional_summary <- function(rear_sim_history, natal_origin_abbrev, natal_origin){
  # if it ever happens that a movement history has over 100 movements, drop it on the 
  # grounds that that's a biologically unrealistic movement history
  if (ncol(rear_sim_history) > 100){
    rear_sim_history[-which(!(is.na(rear_sim_history[,101]))),] -> rear_sim_history
    rear_sim_history <- rear_sim_history[,c(1:100, ncol(rear_sim_history))]
  }
  
  # note which states are overshoot dams
  overshoot_states_list <- list("jdr" = model_states[3:9],
                                "uma" = model_states[3:9],
                                "yak" = model_states[4:9],
                                "wawa" = model_states[4:9],
                                "wawa" = model_states[4:9],
                                "wen" = model_states[6:7],
                                "ent" = model_states[7],
                                "tuc" = model_states[9])
  
  mainstem_states_homing_list <- list("jdr" = model_states[2],
                                "uma" = model_states[2],
                                "yak" = model_states[3],
                                "wawa" = model_states[3],
                                "wawa" = model_states[3],
                                "wen" = model_states[5],
                                "ent" = model_states[6],
                                "tuc" = model_states[8])
  
  
  overshoot_states <- overshoot_states_list[[natal_origin_abbrev]]
  mainstem_state <- mainstem_states_homing_list[[natal_origin_abbrev]]
  
  # reformat input histories
  as.data.frame(rear_sim_history) %>% 
    rownames_to_column("sim_fish_ID") %>% 
    mutate(sim_fish_ID = as.numeric(sim_fish_ID)) %>% 
    pivot_longer(cols = starts_with("t"), names_to = "t", values_to = "state_number") %>% 
    filter(!(is.na(state_number))) -> sim_history
  
  model_states_df <- data.frame(state = model_states, state_number = 1:length(model_states))
  
  sim_history %>% 
    left_join(model_states_df, by = "state_number") %>% 
    mutate(state = gsub(" Upstream", "", state)) %>% 
    mutate(state = gsub(" Mouth", "", state)) -> sim_history
  
  # remove any fish that did not make it to the mainstem state
  mainstem_fish <- unique(subset(sim_history, state == mainstem_state)$sim_fish_ID)
  mainstem_fish_history <- subset(sim_history, state == mainstem_state)
  mainstem_fish_iter <- mainstem_fish_history[!(duplicated(mainstem_fish_history$sim_fish_ID)),]$iter
  
  
  sim_history <- subset(sim_history, sim_fish_ID %in% mainstem_fish)
  
  
  # count number of unique overshoot states (this is overshooting multiple dams)
  sim_history %>% 
    filter(state %in% overshoot_states) %>% 
    group_by(sim_fish_ID) %>% 
    summarise(n_dams_overshot = n_distinct(state)) -> n_overshoot_table
  
  # determine which fish successfully homed
  sim_history %>% 
    filter(state == natal_origin) %>% 
    filter(!(duplicated(sim_fish_ID))) %>% 
    mutate(home = 1) %>% 
    dplyr::select(sim_fish_ID, home) -> homing_fish
  
  all_fish <- data.frame(sim_fish_ID = mainstem_fish, iter = mainstem_fish_iter)
  
  all_fish %>% 
    left_join(homing_fish, by = "sim_fish_ID") %>% 
    left_join(n_overshoot_table, by = "sim_fish_ID") %>% 
    mutate(n_dams_overshot = ifelse(is.na(n_dams_overshot), 0, n_dams_overshot),
           home = ifelse(is.na(home), 0, home))-> all_fish
  
  total_fish <- length(unique(all_fish$sim_fish_ID))
  
  all_fish %>% 
    mutate(dam_overshoots = as.character(ifelse(n_dams_overshot >= 2, "2+", n_dams_overshot))) %>% 
    dplyr::select(-n_dams_overshot) %>%
    #
    group_by(dam_overshoots, iter) %>%
    # complete(dam_overshoots, iter, fill = list(NA)) %>% 
    count(home) %>% 
    group_by(dam_overshoots, iter) %>% 
    mutate(n_per_iter_novershoots = sum(n)) %>% 
    mutate(prop = n/n_per_iter_novershoots) %>% 
    ungroup() %>% 
    group_by(dam_overshoots) %>% 
    filter(home == 1) %>% 
    summarise(prob_of_homing = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) -> homing_n_overshoots_table
  
  
  # as.data.frame(table(all_fish$n_dams_overshot)) %>% 
  #   dplyr::rename(n_overshoot_dams = Var1, nfish = Freq) %>% 
  #   mutate(prop = nfish/total_fish) -> overshoot_table
  
  
  return(homing_n_overshoots_table)
}


jdr_wild_overshoots_homing_conditional_summary <- overshoots_homing_conditional_summary(rear_sim_history = jdr_rear_sim_histories[[1]], natal_origin_abbrev = "jdr",
                                                                natal_origin = "John Day River")
uma_wild_overshoots_homing_conditional_summary <- overshoots_homing_conditional_summary(rear_sim_history = uma_rear_sim_histories[[1]], natal_origin_abbrev = "uma",
                                                                natal_origin = "Umatilla River")
uma_hatchery_overshoots_homing_conditional_summary <- overshoots_homing_conditional_summary(rear_sim_history = uma_rear_sim_histories[[2]], natal_origin_abbrev = "uma",
                                                                    natal_origin = "Umatilla River")
wawa_wild_overshoots_homing_conditional_summary <- overshoots_homing_conditional_summary(rear_sim_history = wawa_rear_sim_histories[[1]], natal_origin_abbrev = "wawa",
                                                                 natal_origin = "Walla Walla River")
wawa_hatchery_overshoots_homing_conditional_summary <- overshoots_homing_conditional_summary(rear_sim_history = wawa_rear_sim_histories[[2]], natal_origin_abbrev = "wawa",
                                                                     natal_origin = "Walla Walla River")
yak_wild_overshoots_homing_conditional_summary <- overshoots_homing_conditional_summary(rear_sim_history = yak_rear_sim_histories[[1]], natal_origin_abbrev = "yak",
                                                                natal_origin = "Yakima River")
wen_wild_overshoots_homing_conditional_summary <- overshoots_homing_conditional_summary(rear_sim_history = wen_rear_sim_histories[[1]], natal_origin_abbrev = "wen",
                                                                natal_origin = "Wenatchee River")
wen_hatchery_overshoots_homing_conditional_summary <- overshoots_homing_conditional_summary(rear_sim_history = wen_rear_sim_histories[[2]], natal_origin_abbrev = "wen",
                                                                    natal_origin = "Wenatchee River")
ent_wild_overshoots_homing_conditional_summary <- overshoots_homing_conditional_summary(rear_sim_history = ent_rear_sim_histories[[1]], natal_origin_abbrev = "ent",
                                                                natal_origin = "Entiat River")
tuc_wild_transported_overshoots_homing_conditional_summary <- overshoots_homing_conditional_summary(rear_sim_history = tuc_transported_rear_sim_histories[[1]], natal_origin_abbrev = "tuc",
                                                                            natal_origin = "Tucannon River")
tuc_wild_not_transported_overshoots_homing_conditional_summary <- overshoots_homing_conditional_summary(rear_sim_history = tuc_not_transported_rear_sim_histories[[1]], 
                                                                                natal_origin_abbrev = "tuc", natal_origin = "Tucannon River")
tuc_hatchery_transported_overshoots_homing_conditional_summary <- overshoots_homing_conditional_summary(rear_sim_history = tuc_transported_rear_sim_histories[[2]], 
                                                                                natal_origin_abbrev = "tuc", natal_origin = "Tucannon River")
tuc_hatchery_not_transported_overshoots_homing_conditional_summary <- overshoots_homing_conditional_summary(rear_sim_history = tuc_not_transported_rear_sim_histories[[2]], 
                                                                                    natal_origin_abbrev = "tuc", natal_origin = "Tucannon River")

# make a summary table - join these outputs

mutate(jdr_wild_overshoots_homing_conditional_summary, origin_rear = "John Day River, Natural") %>% 
  bind_rows(., mutate(uma_wild_overshoots_homing_conditional_summary, origin_rear = "Umatilla River, Natural")) %>% 
  bind_rows(., mutate(uma_hatchery_overshoots_homing_conditional_summary, origin_rear = "Umatilla River, Hatchery")) %>% 
  bind_rows(., mutate(wawa_wild_overshoots_homing_conditional_summary, origin_rear = "Walla Walla River, Natural")) %>% 
  bind_rows(., mutate(wawa_hatchery_overshoots_homing_conditional_summary, origin_rear = "Walla Walla River, Hatchery")) %>% 
  bind_rows(., mutate(yak_wild_overshoots_homing_conditional_summary, origin_rear = "Yakima River, Natural")) %>% 
  bind_rows(., mutate(wen_wild_overshoots_homing_conditional_summary, origin_rear = "Wenatchee River, Natural")) %>% 
  bind_rows(., mutate(wen_hatchery_overshoots_homing_conditional_summary, origin_rear = "Wenatchee River, Hatchery")) %>% 
  bind_rows(., mutate(ent_wild_overshoots_homing_conditional_summary, origin_rear = "Entiat River, Natural")) %>% 
  bind_rows(., mutate(tuc_wild_transported_overshoots_homing_conditional_summary, origin_rear = "Tucannon River, Natural (transported)")) %>% 
  bind_rows(., mutate(tuc_wild_not_transported_overshoots_homing_conditional_summary, origin_rear = "Tucannon River, Natural (not transported)")) %>% 
  bind_rows(., mutate(tuc_hatchery_transported_overshoots_homing_conditional_summary, origin_rear = "Tucannon River, Hatchery (transported)")) %>% 
  bind_rows(., mutate(tuc_hatchery_not_transported_overshoots_homing_conditional_summary, origin_rear = "Tucannon River, Hatchery (not transported)")) %>% 
  ungroup() %>% 
  complete(dam_overshoots, origin_rear, fill = list(NA)) -> all_overshoots_homing_conditional_summary

# reformat as a table
all_overshoots_homing_conditional_summary %>% 
  pivot_wider(names_from = q, values_from = prob_of_homing) %>% 
  mutate(prob_for_table = paste0(round(`0.5`, 2), " (", round(`0.025`, 2), " - ", round(`0.975`, 2), ")")) %>% 
  mutate(prob_for_table = ifelse(prob_for_table == "NA (NA - NA)", "-", prob_for_table)) %>% 
  dplyr::select(dam_overshoots, origin_rear, prob_for_table) %>% 
  pivot_wider(names_from = dam_overshoots, values_from = prob_for_table) -> all_overshoots_homing_conditional_summary_table

write.csv(all_overshoots_homing_conditional_summary_table, file = here::here("figures", "transport", "final_fates", "all_overshoots_homing_conditional_summary_table.csv"))
