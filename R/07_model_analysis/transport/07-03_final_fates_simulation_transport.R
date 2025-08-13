# 07-03_final_fates_simulation

# This script takes the output from the stan model runs 
# generates estimates of final fate distributions

# First, need to load in all of the model runs and all of the packages.
source("R/07_model_analysis/alternative_spill_treatment/07-01_load_stan_models_transport.R")

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


final_fates_simulation_SRW_T <- function(nsim,
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

final_fates_simulation_SRW_NT <- function(nsim,
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

final_fates_simulation_SRH_T <- function(nsim,
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

final_fates_simulation_SRH_NT <- function(nsim,
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

compare_final_fate_rear_type_SR_T <- function(niter, nsim, condition_jitter,
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
      sim_wild <- final_fates_simulation_SRW_T(nsim = nsim,
                                             start_state = 2, states_dates = get_origin_states_dates(envir = SRW_T_envir, origin_select = origin_select, rear = "wild"),
                                             origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                             origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                             temp_data = SRW_T_envir$data$temperature_data, spillwindow_data = SRW_T_envir$data$spill_window_data, 
                                             winterspill_data = SRW_T_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
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
      mutate(rear_type = "wild",
             transport_type = "transported") -> ff_wild_quantiles
    
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
      sim_hatchery <- final_fates_simulation_SRH_T(nsim = nsim,
                                                 start_state = 2, states_dates = get_origin_states_dates(envir = SRH_T_envir, origin_select = origin_select, rear = "hatchery"),
                                                 origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                 origin3 = hatchery_origin_params[3], origin4 = hatchery_origin_params[4], 
                                                 origin5 = hatchery_origin_params[5],
                                                 temp_data = SRH_T_envir$data$temperature_data, spillwindow_data = SRH_T_envir$data$spill_window_data, 
                                                 winterspill_data = SRH_T_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
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
      mutate(rear_type = "hatchery",
             transport_type = "transported") -> ff_hatchery_quantiles
    
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
      sim_wild <- final_fates_simulation_SRW_T(nsim = nsim,
                                             start_state = 2, states_dates = get_origin_states_dates(envir = SRW_T_envir, origin_select = origin_select, rear = "wild"),
                                             origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                             origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                             temp_data = SRW_T_envir$data$temperature_data, spillwindow_data = SRW_T_envir$data$spill_window_data, 
                                             winterspill_data = SRW_T_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_SRH_T(nsim = nsim,
                                                 start_state = 2, states_dates = get_origin_states_dates(envir = SRH_T_envir, origin_select = origin_select, rear = "hatchery"),
                                                 origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                 origin3 = hatchery_origin_params[3], origin4 = hatchery_origin_params[4], 
                                                 origin5 = hatchery_origin_params[5],
                                                 temp_data = SRH_T_envir$data$temperature_data, spillwindow_data = SRH_T_envir$data$spill_window_data, 
                                                 winterspill_data = SRH_T_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
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
      mutate(rear_type = "wild",
             transport_type = "transported") -> ff_wild_quantiles
    
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
      mutate(rear_type = "hatchery",
             transport_type = "transported") -> ff_hatchery_quantiles
    
    ff_wild_quantiles %>% 
      bind_rows(ff_hatchery_quantiles) -> ff_rear_quantiles
  }
  
  
  return(ff_rear_quantiles)
  
}

compare_final_fate_rear_type_SR_NT <- function(niter, nsim, condition_jitter,
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
      sim_wild <- final_fates_simulation_SRW_NT(nsim = nsim,
                                                start_state = 2, states_dates = get_origin_states_dates(envir = SRW_NT_envir, origin_select = origin_select, rear = "wild"),
                                                origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                temp_data = SRW_NT_envir$data$temperature_data, spillwindow_data = SRW_NT_envir$data$spill_window_data, 
                                                winterspill_data = SRW_NT_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
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
      mutate(rear_type = "wild",
             transport_type = "not_transported") -> ff_wild_quantiles
    
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
      sim_hatchery <- final_fates_simulation_SRH_NT(nsim = nsim,
                                                    start_state = 2, states_dates = get_origin_states_dates(envir = SRH_NT_envir, origin_select = origin_select, rear = "hatchery"),
                                                    origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                    origin3 = hatchery_origin_params[3], origin4 = hatchery_origin_params[4], 
                                                    origin5 = hatchery_origin_params[5],
                                                    temp_data = SRH_NT_envir$data$temperature_data, spillwindow_data = SRH_NT_envir$data$spill_window_data, 
                                                    winterspill_data = SRH_NT_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
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
      mutate(rear_type = "hatchery",
             transport_type = "not_transported") -> ff_hatchery_quantiles
    
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
      sim_wild <- final_fates_simulation_SRW_NT(nsim = nsim,
                                                start_state = 2, states_dates = get_origin_states_dates(envir = SRW_NT_envir, origin_select = origin_select, rear = "wild"),
                                                origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                temp_data = SRW_NT_envir$data$temperature_data, spillwindow_data = SRW_NT_envir$data$spill_window_data, 
                                                winterspill_data = SRW_NT_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_SRH_NT(nsim = nsim,
                                                    start_state = 2, states_dates = get_origin_states_dates(envir = SRH_NT_envir, origin_select = origin_select, rear = "hatchery"),
                                                    origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                    origin3 = hatchery_origin_params[3], origin4 = hatchery_origin_params[4], 
                                                    origin5 = hatchery_origin_params[5],
                                                    temp_data = SRH_NT_envir$data$temperature_data, spillwindow_data = SRH_NT_envir$data$spill_window_data, 
                                                    winterspill_data = SRH_NT_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
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
      mutate(rear_type = "wild",
             transport_type = "not_transported") -> ff_wild_quantiles
    
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
      mutate(rear_type = "hatchery",
             transport_type = "not_transported") -> ff_hatchery_quantiles
    
    ff_wild_quantiles %>% 
      bind_rows(ff_hatchery_quantiles) -> ff_rear_quantiles
  }
  
  
  return(ff_rear_quantiles)
  
}

compare_final_fate_transport_type_SRH <- function(niter, nsim, condition_jitter,
                                                  origin_select){
  # If this origin doesn't have a hatchery population, skip it
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    # do nothing!
  }


  else {

    ff_hatchery_transport <- data.frame(state = model_states[1:(length(model_states)-1)])
    ff_hatchery_not_transport <- data.frame(state = model_states[1:(length(model_states)-1)])

    # set up vectors to control which origin parameters are turned on
    # wild_origin_params <- rep(0,6)
    hatchery_origin_params <- rep(0,5)

    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    # wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1

    for(i in 1:niter) {
      # Run final fates simulation for transported hatchery
      sim_hatchery_transport <- final_fates_simulation_SRH_T(nsim = nsim,
                                                             start_state = 2, states_dates = get_origin_states_dates(envir = SRH_T_envir, origin_select = origin_select, rear = "hatchery"),
                                                             origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2],
                                                             origin3 = hatchery_origin_params[3], origin4 = hatchery_origin_params[4],
                                                             origin5 = hatchery_origin_params[5],
                                                             temp_data = SRH_T_envir$data$temperature_data, spillwindow_data = SRH_T_envir$data$spill_window_data,
                                                             winterspill_data = SRH_T_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_hatchery_transport %>%
        bind_cols(sim_hatchery_transport[[2]]) -> ff_hatchery_transport

      # Run final fates simulation for not transported hatchery
      sim_hatchery_not_transport <- final_fates_simulation_SRH_NT(nsim = nsim,
                                                                  start_state = 2, states_dates = get_origin_states_dates(envir = SRH_NT_envir, origin_select = origin_select, rear = "hatchery"),
                                                                  origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2],
                                                                  origin3 = hatchery_origin_params[3], origin4 = hatchery_origin_params[4],
                                                                  origin5 = hatchery_origin_params[5],
                                                                  temp_data = SRH_NT_envir$data$temperature_data, spillwindow_data = SRH_NT_envir$data$spill_window_data,
                                                                  winterspill_data = SRH_NT_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_hatchery_not_transport %>%
        bind_cols(sim_hatchery_not_transport[[2]]) -> ff_hatchery_not_transport
    }


    # Reformat final fates simulation for hatchery, transported
    rownames(ff_hatchery_transport) <- NULL
    column_to_rownames(ff_hatchery_transport, "state") -> ff_hatchery_transport
    colnames(ff_hatchery_transport) <- paste0("iter", 1:niter)

    rownames_to_column(ff_hatchery_transport, "state") -> ff_hatchery_transport

    ff_hatchery_transport %>%
      mutate(state = gsub(" Upstream", "", state)) %>%
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery_transport

    ff_hatchery_transport %>%
      group_by(state) %>%
      summarise(across(where(is.numeric), sum)) -> ff_hatchery_transport

    ff_hatchery_transport %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_transport_long

    ff_hatchery_transport_long %>%
      mutate(prop = count/nsim) %>%
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>%
      pivot_wider(names_from = q, values_from = prob) %>%
      mutate(rear_type = "hatchery",
             transport_type = "transported") -> ff_hatchery_transport_quantiles

    # Reformat final fates simulation for hatchery, not transported
    rownames(ff_hatchery_not_transport) <- NULL
    column_to_rownames(ff_hatchery_not_transport, "state") -> ff_hatchery_not_transport
    colnames(ff_hatchery_not_transport) <- paste0("iter", 1:niter)

    rownames_to_column(ff_hatchery_not_transport, "state") -> ff_hatchery_not_transport

    ff_hatchery_not_transport %>%
      mutate(state = gsub(" Upstream", "", state)) %>%
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery_not_transport

    ff_hatchery_not_transport %>%
      group_by(state) %>%
      summarise(across(where(is.numeric), sum)) -> ff_hatchery_not_transport

    ff_hatchery_not_transport %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_not_transport_long

    ff_hatchery_not_transport_long %>%
      mutate(prop = count/nsim) %>%
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>%
      pivot_wider(names_from = q, values_from = prob) %>%
      mutate(rear_type = "hatchery",
             transport_type = "not_transported") -> ff_hatchery_not_transport_quantiles

    ff_hatchery_transport_quantiles %>%
      bind_rows(ff_hatchery_not_transport_quantiles) -> ff_hatchery_transport_type_quantiles

    return(ff_hatchery_transport_type_quantiles)
  }


}

compare_final_fate_transport_type_SRW <- function(niter, nsim, condition_jitter,
                                                  origin_select){
  # If this origin doesn't have a wild population, skip it
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    # do nothing!
  }


  else {

    ff_wild_transport <- data.frame(state = model_states[1:(length(model_states)-1)])
    ff_wild_not_transport <- data.frame(state = model_states[1:(length(model_states)-1)])

    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    # hatchery_origin_params <- rep(0,5)

    # index to the right origin param to turn it on
    # hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1

    for(i in 1:niter) {
      # Run final fates simulation for transported wild
      sim_wild_transport <- final_fates_simulation_SRW_T(nsim = nsim,
                                                         start_state = 2, states_dates = get_origin_states_dates(envir = SRW_T_envir, origin_select = origin_select, rear = "wild"),
                                                         origin1 = wild_origin_params[1], origin2 = wild_origin_params[2],
                                                         origin3 = wild_origin_params[3], origin4 = wild_origin_params[4],
                                                         origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                         temp_data = SRW_T_envir$data$temperature_data, spillwindow_data = SRW_T_envir$data$spill_window_data,
                                                         winterspill_data = SRW_T_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_wild_transport %>%
        bind_cols(sim_wild_transport[[2]]) -> ff_wild_transport

      # Run final fates simulation for not transported wild
      sim_wild_not_transport <- final_fates_simulation_SRW_NT(nsim = nsim,
                                                              start_state = 2, states_dates = get_origin_states_dates(envir = SRW_NT_envir, origin_select = origin_select, rear = "wild"),
                                                              origin1 = wild_origin_params[1], origin2 = wild_origin_params[2],
                                                              origin3 = wild_origin_params[3], origin4 = wild_origin_params[4],
                                                              origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                              temp_data = SRW_NT_envir$data$temperature_data, spillwindow_data = SRW_NT_envir$data$spill_window_data,
                                                              winterspill_data = SRW_NT_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_wild_not_transport %>%
        bind_cols(sim_wild_not_transport[[2]]) -> ff_wild_not_transport
    }


    # Reformat final fates simulation for wild, transported
    rownames(ff_wild_transport) <- NULL
    column_to_rownames(ff_wild_transport, "state") -> ff_wild_transport
    colnames(ff_wild_transport) <- paste0("iter", 1:niter)

    rownames_to_column(ff_wild_transport, "state") -> ff_wild_transport

    ff_wild_transport %>%
      mutate(state = gsub(" Upstream", "", state)) %>%
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild_transport

    ff_wild_transport %>%
      group_by(state) %>%
      summarise(across(where(is.numeric), sum)) -> ff_wild_transport

    ff_wild_transport %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_transport_long

    ff_wild_transport_long %>%
      mutate(prop = count/nsim) %>%
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>%
      pivot_wider(names_from = q, values_from = prob) %>%
      mutate(rear_type = "wild",
             transport_type = "transported") -> ff_wild_transport_quantiles

    # Reformat final fates simulation for wild, not transported
    rownames(ff_wild_not_transport) <- NULL
    column_to_rownames(ff_wild_not_transport, "state") -> ff_wild_not_transport
    colnames(ff_wild_not_transport) <- paste0("iter", 1:niter)

    rownames_to_column(ff_wild_not_transport, "state") -> ff_wild_not_transport

    ff_wild_not_transport %>%
      mutate(state = gsub(" Upstream", "", state)) %>%
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild_not_transport

    ff_wild_not_transport %>%
      group_by(state) %>%
      summarise(across(where(is.numeric), sum)) -> ff_wild_not_transport

    ff_wild_not_transport %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_not_transport_long

    ff_wild_not_transport_long %>%
      mutate(prop = count/nsim) %>%
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>%
      pivot_wider(names_from = q, values_from = prob) %>%
      mutate(rear_type = "wild",
             transport_type = "not_transported") -> ff_wild_not_transport_quantiles

    ff_wild_transport_quantiles %>%
      bind_rows(ff_wild_not_transport_quantiles) -> ff_wild_transport_type_quantiles

    return(ff_wild_transport_type_quantiles)
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


plot_final_fate_rear_type <- function(ff_comp, natal_origin, transport_type){
  rear_colors <- c(hatchery = "#ff7f00", natural = "#33a02c")
  
  ff_comp %>% 
    mutate(rear_type = ifelse(rear_type == "wild", "natural", rear_type)) -> ff_comp
  
  ff_comp$state <- fct_rev(factor(ff_comp$state, levels = states_order_for_plot))
  
  transport_for_plot <- ifelse(transport_type == "not_transported", "in-river", "transported")
  
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
    ggtitle(paste0(natal_origin, ", ", transport_for_plot))
  
  return(ff_comp_plot)
  
}

plot_final_fate_transport_type <- function(ff_comp_transported, ff_comp_not_transported,
                                           natal_origin, rear_type_select){
  transport_colors <- c("transported" = "#b15928", "not transported" = "#1f78b4")
  
  ff_comp_transported %>% 
    bind_rows(ff_comp_not_transported) -> ff_transport_comp
  
  ff_transport_comp %>% 
    mutate(transport_type = ifelse(transport_type == "not_transported", "not transported", transport_type)) -> ff_transport_comp
  
  ff_transport_comp$state <- fct_rev(factor(ff_transport_comp$state, levels = states_order_for_plot))
  
  ff_transport_comp_one_rear_type <- subset(ff_transport_comp, rear_type == rear_type_select)
  
  rear_for_plot <- ifelse(rear_type_select == "wild", "natural", "hatchery")

  
  ff_comp_plot <- ggplot(ff_transport_comp_one_rear_type, aes(x = state, y = `0.5`, ymin = `0.025`, ymax = `0.975`, color = transport_type)) +
    geom_point(size = 3.5, shape = 18, position=position_dodge(width=0.5)) +
    geom_linerange(position=position_dodge(width=0.5)) +
    coord_flip() +
    ylab("Final Distribution Probability") +
    xlab("Model State") +
    # ggtitle(" ") +
    # Create a scale common to all
    scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_color_manual(values = transport_colors, labels = c("in-river", "transported")) +
    theme(panel.grid.major = element_line(color = "gray90"),
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
          plot.title = element_text(size = 12),
          # axis.text.y = element_text(color = rev(state_significance_colors)),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12)) +
    guides(color=guide_legend(title="Juvenile Pathway")) +
    ggtitle(paste0(natal_origin, ", ", rear_for_plot))
  
  return(ff_comp_plot)
  
}



# new function version to plot rear type and transport type on the same plot

plot_final_fate_rear_and_transport_type <- function(ff_comp, natal_origin, transport_type){
  rear_colors <- c(hatchery = "#ff7f00", natural = "#33a02c")
  
  ff_comp %>% 
    mutate(rear_type = ifelse(rear_type == "wild", "natural", rear_type)) -> ff_comp
  
  ff_comp$state <- fct_rev(factor(ff_comp$state, levels = states_order_for_plot))
  
  transport_for_plot <- ifelse(transport_type == "not_transported", "not transported", "transported ")
  
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
    ggtitle(paste0(natal_origin, ", ", transport_for_plot))
  
  return(ff_comp_plot)
  
}

plot_final_fate_rear_and_transport_type <- function(ff_comp_transported, ff_comp_not_transported,
                                           natal_origin){
  # transport_colors <- c("transported" = "#b15928", "not transported" = "#1f78b4")
  transport_shapes <- c("transported" = 16, "not transported" = 17)
  transport_shape_Sizes <- c("transported" = 2, "not transported" = 3.5)
  rear_colors <- c(hatchery = "#ff7f00", natural = "#33a02c")
  
  ff_comp_transported %>% 
    bind_rows(ff_comp_not_transported) -> ff_transport_comp
  
  ff_transport_comp %>% 
    mutate(rear_type = ifelse(rear_type == "wild", "natural", rear_type)) -> ff_transport_comp
  
  ff_transport_comp %>% 
    mutate(transport_type = ifelse(transport_type == "not_transported", "not transported", transport_type)) -> ff_transport_comp
  
  # turn each combination of rear and transport into a factor level for plotting
  ff_transport_comp %>% 
    mutate(rear_transport_combo = paste0(rear_type, ", ", transport_type)) %>% 
    mutate(rear_transport_combo = factor(rear_transport_combo, levels = c("hatchery, transported",
                                                                          "hatchery, not transported",
                                                                          "natural, transported",
                                                                          "natural, not transported")))-> ff_transport_comp
  
  ff_transport_comp$state <- fct_rev(factor(ff_transport_comp$state, levels = states_order_for_plot))
  
  # change order of observations for plotting
  ff_transport_comp %>% 
    arrange(state, rear_transport_combo) %>% 
    ungroup() -> ff_transport_comp
  
  
  ff_comp_plot <- ggplot(ff_transport_comp, aes(x = state, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                color = rear_type, shape = transport_type, group = rear_transport_combo)) +
    # geom_point(aes(size = transport_type), position=position_dodge(width=0.5)) +
    geom_point(size = 2.5, position=position_dodge(width=0.75)) +
    geom_linerange(position=position_dodge(width=0.75)) +
    coord_flip() +
    ylab("Final Distribution Probability") +
    xlab("Model State") +
    # ggtitle(" ") +
    # Create a scale common to all
    scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_color_manual(values = rear_colors) +
    scale_shape_manual(values = transport_shapes, labels = c("in-river", "transported")) +
    theme(panel.grid.major = element_line(color = "gray90"),
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
          plot.title = element_text(size = 12),
          # axis.text.y = element_text(color = rev(state_significance_colors)),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12)) +
    guides(color=guide_legend(title="Rear Type", override.aes = list(shape = NA,
                                                                     linewidth = 3)),
           shape=guide_legend(title="Juvenile Pathway")) +
    ggtitle(natal_origin)
  
  return(ff_comp_plot)
  
}



ff_iter <- 1000
ff_nsim <- 1000

# If you just want to regenerate plots:
# load(file = here::here("figures", "transport", "final_fates", "FF_comp_data.rda"))
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
load(file = here::here("figures", "transport", "final_fates", "simulation_runs", "aso_not_transported_rear_comp_data.rda"))
load(file = here::here("figures", "transport", "final_fates", "simulation_runs", "aso_transported_rear_comp_data.rda"))
load(file = here::here("figures", "transport", "final_fates", "simulation_runs", "cle_not_transported_rear_comp_data.rda"))
load(file = here::here("figures", "transport", "final_fates", "simulation_runs", "cle_transported_rear_comp_data.rda"))
load(file = here::here("figures", "transport", "final_fates", "simulation_runs", "gr_not_transported_rear_comp_data.rda"))
load(file = here::here("figures", "transport", "final_fates", "simulation_runs", "gr_transported_rear_comp_data.rda"))
load(file = here::here("figures", "transport", "final_fates", "simulation_runs", "imn_not_transported_rear_comp_data.rda"))
load(file = here::here("figures", "transport", "final_fates", "simulation_runs", "imn_transported_rear_comp_data.rda"))
load(file = here::here("figures", "transport", "final_fates", "simulation_runs", "sal_not_transported_rear_comp_data.rda"))
load(file = here::here("figures", "transport", "final_fates", "simulation_runs", "sal_transported_rear_comp_data.rda"))
load(file = here::here("figures", "transport", "final_fates", "simulation_runs", "tuc_not_transported_rear_comp_data.rda"))
load(file = here::here("figures", "transport", "final_fates", "simulation_runs", "tuc_transported_rear_comp_data.rda"))

aso_transported_ff_rear_comp_median$transport_type <- "transported"
aso_not_transported_ff_rear_comp_median$transport_type <- "not_transported"


## Snake River
# Asotin Creek rear comparison - transported
aso_transported_ff_rear_comp_median <- compare_final_fate_rear_type_SR_T(niter = ff_iter, nsim = ff_nsim, origin_select = "Asotin Creek", condition_jitter = FALSE)
save(aso_transported_ff_rear_comp_median, file = here::here("figures", "transport", "final_fates", "simulation_runs", "aso_transported_rear_comp_data.rda"))
aso_transported_ff_rear_comp_median_plot <- plot_final_fate_rear_type(aso_transported_ff_rear_comp_median, natal_origin = "Asotin Creek", transport_type = "transported")
ggsave(here::here("figures", "transport", "final_fates", "aso_transported_ff_rear_comp_median_plot.png"), aso_transported_ff_rear_comp_median_plot, height = 8, width = 8)

# Asotin Creek rear comparison - not transported
aso_not_transported_ff_rear_comp_median <- compare_final_fate_rear_type_SR_NT(niter = ff_iter, nsim = ff_nsim, origin_select = "Asotin Creek", condition_jitter = FALSE)
save(aso_not_transported_ff_rear_comp_median, file = here::here("figures", "transport", "final_fates", "simulation_runs", "aso_not_transported_rear_comp_data.rda"))
aso_not_transported_ff_rear_comp_median_plot <- plot_final_fate_rear_type(aso_not_transported_ff_rear_comp_median, natal_origin = "Asotin Creek", transport_type = "not_transported")
ggsave(here::here("figures", "transport", "final_fates", "aso_not_transported_ff_rear_comp_median_plot.png"), aso_not_transported_ff_rear_comp_median_plot, height = 8, width = 8)

# Asotin Creek transport type comparison - hatchery
# no hatchery fish from asotin!
# aso_hatchery_ff_transport_comp_median <- compare_final_fate_transport_type_SRH(niter = ff_iter, nsim = ff_nsim, origin_select = "Asotin Creek", condition_jitter = FALSE)
# save(aso_hatchery_ff_transport_comp_median, file = here::here("figures", "transport", "final_fates", "simulation_runs", "aso_hatchery_transport_comp_data.rda"))
# aso_hatchery_ff_transport_comp_median_plot <- plot_final_fate_transport_type(aso_hatchery_ff_transport_comp_median, natal_origin = "Asotin Creek", rear_type = "hatchery")
# ggsave(here::here("figures", "transport", "final_fates", "aso_hatchery_ff_transport_comp_median_plot.png"), aso_hatchery_ff_transport_comp_median_plot, height = 8, width = 8)

# Asotin Creek transport type comparison - wild
aso_wild_ff_transport_comp_median_plot <- plot_final_fate_transport_type(ff_comp_transported = aso_transported_ff_rear_comp_median, 
                                                                         ff_comp_not_transported = aso_not_transported_ff_rear_comp_median, 
                                                                         natal_origin = "Asotin Creek", rear_type_select = "wild")
ggsave(here::here("figures", "transport", "final_fates", "aso_wild_ff_transport_comp_median_plot.png"), aso_wild_ff_transport_comp_median_plot, height = 8, width = 8)

# Asotin Creek rear and transport comparison plot
aso_ff_rear_transport_comp_plot <- plot_final_fate_rear_and_transport_type(ff_comp_transported = aso_transported_ff_rear_comp_median, 
                                                                           ff_comp_not_transported = aso_not_transported_ff_rear_comp_median,
                                                                           natal_origin = "Asotin Creek")

ggsave(here::here("figures", "transport", "final_fates", "aso_ff_rear_transport_comp_plot.png"), aso_ff_rear_transport_comp_plot, height = 8, width = 8)


# Clearwater River rear comparison - transported
cle_transported_ff_rear_comp_median <- compare_final_fate_rear_type_SR_T(niter = ff_iter, nsim = ff_nsim, origin_select = "Clearwater River", condition_jitter = FALSE)
save(cle_transported_ff_rear_comp_median, file = here::here("figures", "transport", "final_fates", "simulation_runs", "cle_transported_rear_comp_data.rda"))
cle_transported_ff_rear_comp_median_plot <- plot_final_fate_rear_type(cle_transported_ff_rear_comp_median, natal_origin = "Clearwater River", transport_type = "transported")
ggsave(here::here("figures", "transport", "final_fates", "cle_transported_ff_rear_comp_median_plot.png"), cle_transported_ff_rear_comp_median_plot, height = 8, width = 8)

# Clearwater River rear comparison - not transported
cle_not_transported_ff_rear_comp_median <- compare_final_fate_rear_type_SR_NT(niter = ff_iter, nsim = ff_nsim, origin_select = "Clearwater River", condition_jitter = FALSE)
save(cle_not_transported_ff_rear_comp_median, file = here::here("figures", "transport", "final_fates", "simulation_runs", "cle_not_transported_rear_comp_data.rda"))
cle_not_transported_ff_rear_comp_median_plot <- plot_final_fate_rear_type(cle_not_transported_ff_rear_comp_median, natal_origin = "Clearwater River", transport_type = "not_transported")
ggsave(here::here("figures", "transport", "final_fates", "cle_not_transported_ff_rear_comp_median_plot.png"), cle_not_transported_ff_rear_comp_median_plot, height = 8, width = 8)

# Clearwater River transport type comparison - hatchery
cle_hatchery_ff_transport_comp_median_plot <- plot_final_fate_transport_type(ff_comp_transported = cle_transported_ff_rear_comp_median, 
                                                                         ff_comp_not_transported = cle_not_transported_ff_rear_comp_median, 
                                                                         natal_origin = "Clearwater River", rear_type_select = "hatchery")
ggsave(here::here("figures", "transport", "final_fates", "cle_hatchery_ff_transport_comp_median_plot.png"), cle_hatchery_ff_transport_comp_median_plot, height = 8, width = 8)

# Clearwater River transport type comparison - wild
cle_wild_ff_transport_comp_median_plot <- plot_final_fate_transport_type(ff_comp_transported = cle_transported_ff_rear_comp_median, 
                                                                         ff_comp_not_transported = cle_not_transported_ff_rear_comp_median, 
                                                                         natal_origin = "Clearwater River", rear_type_select = "wild")
ggsave(here::here("figures", "transport", "final_fates", "cle_wild_ff_transport_comp_median_plot.png"), cle_wild_ff_transport_comp_median_plot, height = 8, width = 8)

# Clearwater River rear and transport comparison plot
cle_ff_rear_transport_comp_plot <- plot_final_fate_rear_and_transport_type(ff_comp_transported = cle_transported_ff_rear_comp_median, 
                                                                           ff_comp_not_transported = cle_not_transported_ff_rear_comp_median,
                                                                           natal_origin = "Clearwater River")

ggsave(here::here("figures", "transport", "final_fates", "cle_ff_rear_transport_comp_plot.png"), cle_ff_rear_transport_comp_plot, height = 8, width = 8)


# Salmon River rear comparison - transported
sal_transported_ff_rear_comp_median <- compare_final_fate_rear_type_SR_T(niter = ff_iter, nsim = ff_nsim, origin_select = "Salmon River", condition_jitter = FALSE)
save(sal_transported_ff_rear_comp_median, file = here::here("figures", "transport", "final_fates", "simulation_runs", "sal_transported_rear_comp_data.rda"))
sal_transported_ff_rear_comp_median_plot <- plot_final_fate_rear_type(sal_transported_ff_rear_comp_median, natal_origin = "Salmon River", transport_type = "transported")
ggsave(here::here("figures", "transport", "final_fates", "sal_transported_ff_rear_comp_median_plot.png"), sal_transported_ff_rear_comp_median_plot, height = 8, width = 8)

# Salmon River rear comparison - not transported
sal_not_transported_ff_rear_comp_median <- compare_final_fate_rear_type_SR_NT(niter = ff_iter, nsim = ff_nsim, origin_select = "Salmon River", condition_jitter = FALSE)
save(sal_not_transported_ff_rear_comp_median, file = here::here("figures", "transport", "final_fates", "simulation_runs", "sal_not_transported_rear_comp_data.rda"))
sal_not_transported_ff_rear_comp_median_plot <- plot_final_fate_rear_type(sal_not_transported_ff_rear_comp_median, natal_origin = "Salmon River", transport_type = "not_transported")
ggsave(here::here("figures", "transport", "final_fates", "sal_not_transported_ff_rear_comp_median_plot.png"), sal_not_transported_ff_rear_comp_median_plot, height = 8, width = 8)

# Salmon River transport type comparison - hatchery
sal_hatchery_ff_transport_comp_median_plot <- plot_final_fate_transport_type(ff_comp_transported = sal_transported_ff_rear_comp_median, 
                                                                             ff_comp_not_transported = sal_not_transported_ff_rear_comp_median, 
                                                                             natal_origin = "Salmon River", rear_type_select = "hatchery")
ggsave(here::here("figures", "transport", "final_fates", "sal_hatchery_ff_transport_comp_median_plot.png"), sal_hatchery_ff_transport_comp_median_plot, height = 8, width = 8)

# Salmon River transport type comparison - wild
sal_wild_ff_transport_comp_median_plot <- plot_final_fate_transport_type(ff_comp_transported = sal_transported_ff_rear_comp_median, 
                                                                         ff_comp_not_transported = sal_not_transported_ff_rear_comp_median, 
                                                                         natal_origin = "Salmon River", rear_type_select = "wild")
ggsave(here::here("figures", "transport", "final_fates", "sal_wild_ff_transport_comp_median_plot.png"), sal_wild_ff_transport_comp_median_plot, height = 8, width = 8)

# Salmon River rear and transport comparison plot
sal_ff_rear_transport_comp_plot <- plot_final_fate_rear_and_transport_type(ff_comp_transported = sal_transported_ff_rear_comp_median, 
                                                                           ff_comp_not_transported = sal_not_transported_ff_rear_comp_median,
                                                                           natal_origin = "Salmon River")

ggsave(here::here("figures", "transport", "final_fates", "sal_ff_rear_transport_comp_plot.png"), sal_ff_rear_transport_comp_plot, height = 8, width = 8)

# Grande Ronde River rear comparison - transported
gr_transported_ff_rear_comp_median <- compare_final_fate_rear_type_SR_T(niter = ff_iter, nsim = ff_nsim, origin_select = "Grande Ronde River", condition_jitter = FALSE)
save(gr_transported_ff_rear_comp_median, file = here::here("figures", "transport", "final_fates", "simulation_runs", "gr_transported_rear_comp_data.rda"))
gr_transported_ff_rear_comp_median_plot <- plot_final_fate_rear_type(gr_transported_ff_rear_comp_median, natal_origin = "Grande Ronde River", transport_type = "transported")
ggsave(here::here("figures", "transport", "final_fates", "gr_transported_ff_rear_comp_median_plot.png"), gr_transported_ff_rear_comp_median_plot, height = 8, width = 8)

# Grande Ronde River rear comparison - not transported
gr_not_transported_ff_rear_comp_median <- compare_final_fate_rear_type_SR_NT(niter = ff_iter, nsim = ff_nsim, origin_select = "Grande Ronde River", condition_jitter = FALSE)
save(gr_not_transported_ff_rear_comp_median, file = here::here("figures", "transport", "final_fates", "simulation_runs", "gr_not_transported_rear_comp_data.rda"))
gr_not_transported_ff_rear_comp_median_plot <- plot_final_fate_rear_type(gr_not_transported_ff_rear_comp_median, natal_origin = "Grande Ronde River", transport_type = "not_transported")
ggsave(here::here("figures", "transport", "final_fates", "gr_not_transported_ff_rear_comp_median_plot.png"), gr_not_transported_ff_rear_comp_median_plot, height = 8, width = 8)

# Grande Ronde River transport type comparison - hatchery
gr_hatchery_ff_transport_comp_median_plot <- plot_final_fate_transport_type(ff_comp_transported = gr_transported_ff_rear_comp_median, 
                                                                            ff_comp_not_transported = gr_not_transported_ff_rear_comp_median, 
                                                                            natal_origin = "Grande Ronde River", rear_type_select = "hatchery")
ggsave(here::here("figures", "transport", "final_fates", "gr_hatchery_ff_transport_comp_median_plot.png"), gr_hatchery_ff_transport_comp_median_plot, height = 8, width = 8)

# Grande Ronde River transport type comparison - wild
gr_wild_ff_transport_comp_median_plot <- plot_final_fate_transport_type(ff_comp_transported = gr_transported_ff_rear_comp_median, 
                                                                        ff_comp_not_transported = gr_not_transported_ff_rear_comp_median, 
                                                                        natal_origin = "Grande Ronde River", rear_type_select = "wild")
ggsave(here::here("figures", "transport", "final_fates", "gr_wild_ff_transport_comp_median_plot.png"), gr_wild_ff_transport_comp_median_plot, height = 8, width = 8)

# Grande Ronde River rear and transport comparison plot
gr_ff_rear_transport_comp_plot <- plot_final_fate_rear_and_transport_type(ff_comp_transported = gr_transported_ff_rear_comp_median, 
                                                                           ff_comp_not_transported = gr_not_transported_ff_rear_comp_median,
                                                                           natal_origin = "Grande Ronde River")

ggsave(here::here("figures", "transport", "final_fates", "gr_ff_rear_transport_comp_plot.png"), gr_ff_rear_transport_comp_plot, height = 8, width = 8)

# Imnaha River rear comparison - transported
imn_transported_ff_rear_comp_median <- compare_final_fate_rear_type_SR_T(niter = ff_iter, nsim = ff_nsim, origin_select = "Imnaha River", condition_jitter = FALSE)
save(imn_transported_ff_rear_comp_median, file = here::here("figures", "transport", "final_fates", "simulation_runs", "imn_transported_rear_comp_data.rda"))
imn_transported_ff_rear_comp_median_plot <- plot_final_fate_rear_type(imn_transported_ff_rear_comp_median, natal_origin = "Imnaha River", transport_type = "transported")
ggsave(here::here("figures", "transport", "final_fates", "imn_transported_ff_rear_comp_median_plot.png"), imn_transported_ff_rear_comp_median_plot, height = 8, width = 8)

# Imnaha River rear comparison - not transported
imn_not_transported_ff_rear_comp_median <- compare_final_fate_rear_type_SR_NT(niter = ff_iter, nsim = ff_nsim, origin_select = "Imnaha River", condition_jitter = FALSE)
save(imn_not_transported_ff_rear_comp_median, file = here::here("figures", "transport", "final_fates", "simulation_runs", "imn_not_transported_rear_comp_data.rda"))
imn_not_transported_ff_rear_comp_median_plot <- plot_final_fate_rear_type(imn_not_transported_ff_rear_comp_median, natal_origin = "Imnaha River", transport_type = "not_transported")
ggsave(here::here("figures", "transport", "final_fates", "imn_not_transported_ff_rear_comp_median_plot.png"), imn_not_transported_ff_rear_comp_median_plot, height = 8, width = 8)

# Imnaha River transport type comparison - hatchery
imn_hatchery_ff_transport_comp_median_plot <- plot_final_fate_transport_type(ff_comp_transported = imn_transported_ff_rear_comp_median, 
                                                                             ff_comp_not_transported = imn_not_transported_ff_rear_comp_median, 
                                                                             natal_origin = "Imnaha River", rear_type_select = "hatchery")
ggsave(here::here("figures", "transport", "final_fates", "imn_hatchery_ff_transport_comp_median_plot.png"), imn_hatchery_ff_transport_comp_median_plot, height = 8, width = 8)

# Imnaha River transport type comparison - wild
imn_wild_ff_transport_comp_median_plot <- plot_final_fate_transport_type(ff_comp_transported = imn_transported_ff_rear_comp_median, 
                                                                         ff_comp_not_transported = imn_not_transported_ff_rear_comp_median, 
                                                                         natal_origin = "Imnaha River", rear_type_select = "wild")
ggsave(here::here("figures", "transport", "final_fates", "imn_wild_ff_transport_comp_median_plot.png"), imn_wild_ff_transport_comp_median_plot, height = 8, width = 8)

# Imnaha River rear and transport comparison plot
imn_ff_rear_transport_comp_plot <- plot_final_fate_rear_and_transport_type(ff_comp_transported = imn_transported_ff_rear_comp_median, 
                                                                           ff_comp_not_transported = imn_not_transported_ff_rear_comp_median,
                                                                           natal_origin = "Imnaha River")

ggsave(here::here("figures", "transport", "final_fates", "imn_ff_rear_transport_comp_plot.png"), imn_ff_rear_transport_comp_plot, height = 8, width = 8)

# Tucannon River rear comparison - transported
tuc_transported_ff_rear_comp_median <- compare_final_fate_rear_type_SR_T(niter = ff_iter, nsim = ff_nsim, origin_select = "Tucannon River", condition_jitter = FALSE)
save(tuc_transported_ff_rear_comp_median, file = here::here("figures", "transport", "final_fates", "simulation_runs", "tuc_transported_rear_comp_data.rda"))
tuc_transported_ff_rear_comp_median_plot <- plot_final_fate_rear_type(tuc_transported_ff_rear_comp_median, natal_origin = "Tucannon River", transport_type = "transported")
ggsave(here::here("figures", "transport", "final_fates", "tuc_transported_ff_rear_comp_median_plot.png"), tuc_transported_ff_rear_comp_median_plot, height = 8, width = 8)

# Tucannon River rear comparison - not transported
tuc_not_transported_ff_rear_comp_median <- compare_final_fate_rear_type_SR_NT(niter = ff_iter, nsim = ff_nsim, origin_select = "Tucannon River", condition_jitter = FALSE)
save(tuc_not_transported_ff_rear_comp_median, file = here::here("figures", "transport", "final_fates", "simulation_runs", "tuc_not_transported_rear_comp_data.rda"))
tuc_not_transported_ff_rear_comp_median_plot <- plot_final_fate_rear_type(tuc_not_transported_ff_rear_comp_median, natal_origin = "Tucannon River", transport_type = "not_transported")
ggsave(here::here("figures", "transport", "final_fates", "tuc_not_transported_ff_rear_comp_median_plot.png"), tuc_not_transported_ff_rear_comp_median_plot, height = 8, width = 8)

# Tucannon River transport type comparison - hatchery
tuc_hatchery_ff_transport_comp_median_plot <- plot_final_fate_transport_type(ff_comp_transported = tuc_transported_ff_rear_comp_median, 
                                                                             ff_comp_not_transported = tuc_not_transported_ff_rear_comp_median, 
                                                                             natal_origin = "Tucannon River", rear_type_select = "hatchery")
ggsave(here::here("figures", "transport", "final_fates", "tuc_hatchery_ff_transport_comp_median_plot.png"), tuc_hatchery_ff_transport_comp_median_plot, height = 8, width = 8)

# Tucannon River transport type comparison - wild
tuc_wild_ff_transport_comp_median_plot <- plot_final_fate_transport_type(ff_comp_transported = tuc_transported_ff_rear_comp_median, 
                                                                         ff_comp_not_transported = tuc_not_transported_ff_rear_comp_median, 
                                                                         natal_origin = "Tucannon River", rear_type_select = "wild")
ggsave(here::here("figures", "transport", "final_fates", "tuc_wild_ff_transport_comp_median_plot.png"), tuc_wild_ff_transport_comp_median_plot, height = 8, width = 8)


# Tucannon River rear and transport comparison plot
tuc_ff_rear_transport_comp_plot <- plot_final_fate_rear_and_transport_type(ff_comp_transported = tuc_transported_ff_rear_comp_median, 
                                                                           ff_comp_not_transported = tuc_not_transported_ff_rear_comp_median,
                                                                           natal_origin = "Tucannon River")

ggsave(here::here("figures", "transport", "final_fates", "tuc_ff_rear_transport_comp_plot.png"), tuc_ff_rear_transport_comp_plot, height = 8, width = 8)

# SAVE ALL FINAL FATES OUTPUTS
# FF_comp_data <- list(wen_ff_comp_median = wen_ff_comp_median,
#                      ent_ff_comp_median = ent_ff_comp_median,
#                      oka_ff_comp_median = oka_ff_comp_median,
#                      met_ff_comp_median = met_ff_comp_median,
#                      des_ff_comp_median = des_ff_comp_median,
#                      jdr_ff_comp_median = jdr_ff_comp_median,
#                      fif_ff_comp_median = fif_ff_comp_median,
#                      uma_ff_comp_median = uma_ff_comp_median,
#                      yak_ff_comp_median = yak_ff_comp_median,
#                      wawa_ff_comp_median = wawa_ff_comp_median,
#                      aso_ff_comp_median = aso_ff_comp_median,
#                      cle_ff_comp_median = cle_ff_comp_median,
#                      sal_ff_comp_median = sal_ff_comp_median,
#                      gr_ff_comp_median = gr_ff_comp_median,
#                      imn_ff_comp_median = imn_ff_comp_median,
#                      tuc_ff_comp_median = tuc_ff_comp_median)

# save(FF_comp_data, file = here::here("figures", "transport", "final_fates", "simulation_runs", "FF_comp_data.rda"))
