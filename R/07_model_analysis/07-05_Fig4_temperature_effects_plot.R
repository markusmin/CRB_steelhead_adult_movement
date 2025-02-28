# 07-05_Fig4_temperature_effects_plot

# This script takes the output from the stan model runs and
# plots the effects of temperature

# First, need to load in all of the model runs and all of the packages.
source("R/07_model_analysis/07-01_load_stan_models.R")

# Load the temperature summary
window_temp_summary <- read.csv(here::here("Data", "covariate_data", "model_inputs", "window_temps_summary.csv"), row.names = 1)

#### Function that extracts the DE parameters
extract_DE_parameters <- function(fit, fit_summary){
  # extract b0 as an array
  parameters <- fit_summary$variable
  parameters[grepl("_alpha|_beta" , parameters)] -> param_subset
  
  # add in fifteenmile_beta and imnaha_beta
  DE_params <- c(param_subset[1:23], "fifteenmile_beta", "imnaha_beta", param_subset[24:length(param_subset)])
  
  
  # arrange all parameter values together in a df
  DE_param_matrix <- matrix(data = 0, nrow = length(DE_params),
                            length(as.matrix(fit[,,1])))
  
  rownames(DE_param_matrix) <- DE_params
  
  
  
  # extract the draws for each and store in a matrix
  for(i in 1:nrow(DE_param_matrix)){
    if (DE_params[i] %in% c("fifteenmile_beta", "imnaha_beta")){
      DE_param_matrix[i, ] <- 0
    } else {
      DE_param_matrix[i, ] <- as.matrix(fit[,,DE_params[i]])
    }
    
  }
  
  return(DE_param_matrix)
}

calculate_weighted_DE <- function(DE_param_matrix, tributary, tributary_design_matrices_array,
                                  covariate_experiences){
  
  ## calculate estimated detection efficiency by year
  # get the index for the tributary state (needs to be mouth since that's where we're correcting for DE)
  tributary_state <- intersect(grep(tributary, model_states, ignore.case = TRUE), grep("Mouth", model_states, ignore.case = TRUE))
  
  tributary_from_state_df <- data.frame(from = c(rep(2,4),3,3,5,6,7,7,8,9,9),
                                        to = grep("Mouth", model_states, ignore.case = FALSE))
  
  tributary_mainstem_state <- subset(tributary_from_state_df, to == tributary_state)$from
  
  # use the state indexing to get the right design matrix
  trib_design_matrix <- tributary_design_matrices_array[,,tributary_state]
  
  # create a matrix to store DE per year, per iter
  niter <- 4000
  DE_matrix <- matrix(nrow = nrow(trib_design_matrix), ncol = niter)
  
  # for each run year, get a confidence interval around detection efficiency by using the different draws
  for (i in 1:nrow(trib_design_matrix)){
    # if there is no intercept term, that means there is no DE correction for that year - so skip and leave as NA
    if(sum(trib_design_matrix[i,1:20]) == 0){
      
      
    } else {
      for (j in 1:niter){
        eta <- sum(trib_design_matrix[i,] * DE_param_matrix[,j])
        DE_matrix[i,j] <- exp(eta)/(1 + exp(eta))
      }
      
    }
    
    
  }
  
  colnames(DE_matrix) <- paste0("iter", 1:niter)
  
  DE_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("year") %>% 
    mutate(year = as.numeric(year)) %>% 
    pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "detection_probability") %>% 
    group_by(year) %>% 
    summarise(prob = quantile(detection_probability, c(0.025, 0.5, 0.975), na.rm = T), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, values_from = prob) %>% 
    mutate(rear = "wild") -> DE_matrix_long
  
  DE_matrix_long$year_actual <- 2004:2023
  
  ## calculate number of transitions into that tributary by year
  covariate_experiences %>% 
    mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 41)) -> covariate_experiences
  
  # get a table of counts by run year, into that tributary
  as.data.frame(table(subset(covariate_experiences, state == tributary_mainstem_state & next_state == tributary_state)$year)) %>% 
    dplyr::rename(index = Var1, count = Freq) %>% 
    mutate(index = as.numeric(as.character(index))) -> trib_entries_by_year
  
  # add in the actual year (not index year). 2004 = year 1
  year_indices <- data.frame(index = 1:21, year_actual = 2004:2024)
  trib_entries_by_year %>% 
    left_join(year_indices, by = "index") -> trib_entries_by_year
  
  DE_matrix_long %>% 
    left_join(dplyr::select(trib_entries_by_year, count, year_actual), by = "year_actual") -> DE_by_year
  
  # create a weighted average
  # make sure to drop any years where DE isn't estimated in the model
  DE_by_year %>% 
    filter(!(is.na(`0.5`))) %>% 
    group_by(rear) %>% 
    mutate(weight = count/sum(count, na.rm = T)) %>% 
    summarise(weighted_DE = sum(`0.5`*weight, na.rm = T)) -> weighted_DE
  
  return(list(DE_by_year = DE_by_year, weighted_average = weighted_DE$weighted_DE))
}

### extract all params ###
UCW_DE_param_matrix <- extract_DE_parameters(fit = UCW_fit, fit_summary = UCW_fit_summary)
UCH_DE_param_matrix <- extract_DE_parameters(fit = UCH_fit, fit_summary = UCH_fit_summary)
MCW_DE_param_matrix <- extract_DE_parameters(fit = MCW_fit, fit_summary = MCW_fit_summary)
MCH_DE_param_matrix <- extract_DE_parameters(fit = MCH_fit, fit_summary = MCH_fit_summary)
SRW_DE_param_matrix <- extract_DE_parameters(fit = SRW_fit, fit_summary = SRW_fit_summary)
SRH_DE_param_matrix <- extract_DE_parameters(fit = SRH_fit, fit_summary = SRH_fit_summary)

#### Temperature plot ####

# Here, we will probably want to look at only certain movements (?)
# For example, overshoot, holding in Deschutes, straying
# But we could plot all of them

# We have DPS-wide temperature effects and origin-specific ones, 
# and we have a winter/spring and summer/fall param (temp0 or temp1)

# to capture posterior correlations, we'll need to select that same iteration
# for all parameters that affect movement

# Add the DPS
DPS_map <- data.frame(
  natal_origin = c("Middle Columbia", "Upper Columbia", "Snake River"),
  hatchery = rep(999, 3),
  wild = rep(999, 3)
)
origin_param_map %>% 
  bind_rows(DPS_map) -> origin_param_map


extract_covariate_experiences <- function(envir, rear, origin_select = NULL){
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
  # WEL and LGR because we have slow v fast there
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
  
  
  # Now, keep only the origin selected, unless you want the whole DPS
  if (!(is.null(origin_select))){
    if(rear == "wild"){
      origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$wild
    } else {
      origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$hatchery
    }
    
    subset(pop_states_dates_covariates, origin == origin_numeric) -> origin_states_dates_covariates
  } else {
    # if we leave origin select as null, return the full DPS
    pop_states_dates_covariates -> origin_states_dates_covariates
  }
  
  
  return(origin_states_dates_covariates)
  
}



# Tell it the movements for which you want to estimate temperature effects
# movements are formatted as matrix, with column for from and column for to
estimate_temp_effect_UCW <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- rep(0,3)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = wild_origin_params[1]
  origin2 = wild_origin_params[2]
  origin3 = wild_origin_params[3]
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    
    # If it's a mainstem state, get the spill data. If it's not, ignore it.
    
    if (from %in% c(1:9)){
      # get median winter spill for this state
      UCW_states_dates <- data.frame(state = as.vector(UCW_envir$data$y),
                                     date = as.vector(UCW_envir$data$transition_dates))
      spillwindow_data <- UCW_envir$data$spill_window_data
      med_spillwindow <- median(spillwindow_data[subset(UCW_states_dates, state == from)$date,from])
      
      # get median spill window for this state
      winterspill_data <- UCW_envir$data$winter_spill_days_data
      med_winterspill <- median(winterspill_data[,from])      
      
    } else {
      med_spillwindow <- 0
      med_winterspill <- 0
    }
    
    
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j, iter, i] <- exp(b0_array_UCW[from,to,iter] +
                                                  # btemp0_array_UCW[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                  btemp1_array_UCW[from,to,iter]*temp_predict[j] + 
                                                  bspillwindow_array_UCW[from,to,iter]*med_spillwindow +
                                                  bwinterspill_array_UCW[from,to,iter]*med_winterspill +
                                                  # btemp0xorigin1_array_UCW[from,to,iter]*origin1 +
                                                  btemp1xorigin1_array_UCW[from,to,iter]*temp_predict[j]*origin1 +
                                                  # btemp0xorigin2_array_UCW[from,to,iter]*origin2 + 
                                                  btemp1xorigin2_array_UCW[from,to,iter]*temp_predict[j]*origin2 + 
                                                  # btemp0xorigin3_array_UCW[from,to,iter]*origin3 +
                                                  btemp1xorigin3_array_UCW[from,to,iter]*temp_predict[j]*origin3 +
                                                  borigin1_array_UCW[from,to,iter]*origin1 +
                                                  borigin2_array_UCW[from,to,iter]*origin2 +
                                                  borigin3_array_UCW[from,to,iter]*origin3)/
          sum(exp(b0_array_UCW[from,possible_movements,iter] +
                    # btemp0_array_UCW[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_UCW[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_UCW[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_UCW[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_UCW[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_UCW[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_UCW[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_UCW[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    # btemp0xorigin3_array_UCW[from,possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_UCW[from,possible_movements,iter]*temp_predict[j]*origin3 +
                    borigin1_array_UCW[from,possible_movements,iter]*origin1 +
                    borigin2_array_UCW[from,possible_movements,iter]*origin2 +
                    borigin3_array_UCW[from,possible_movements,iter]*origin3))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

estimate_temp_effect_UCH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,3)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  origin3 = hatchery_origin_params[3]
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # If it's a mainstem state, get the spill data. If it's not, ignore it.
    
    if (from %in% c(1:9)){
      # get median winter spill for this state
      UCH_states_dates <- data.frame(state = as.vector(UCH_envir$data$y),
                                     date = as.vector(UCH_envir$data$transition_dates))
      spillwindow_data <- UCH_envir$data$spill_window_data
      med_spillwindow <- median(spillwindow_data[subset(UCH_states_dates, state == from)$date,from])
      
      # get median spill window for this state
      winterspill_data <- UCH_envir$data$winter_spill_days_data
      med_winterspill <- median(winterspill_data[,from])      
      
    } else {
      med_spillwindow <- 0
      med_winterspill <- 0
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j, iter, i] <- exp(b0_array_UCH[from,to,iter] +
                                                  # btemp0_array_UCH[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                  btemp1_array_UCH[from,to,iter]*temp_predict[j] + 
                                                  bspillwindow_array_UCH[from,to,iter]*med_spillwindow +
                                                  bwinterspill_array_UCH[from,to,iter]*med_winterspill +
                                                  # btemp0xorigin1_array_UCH[from,to,iter]*origin1 +
                                                  btemp1xorigin1_array_UCH[from,to,iter]*temp_predict[j]*origin1 +
                                                  # btemp0xorigin2_array_UCH[from,to,iter]*origin2 + 
                                                  btemp1xorigin2_array_UCH[from,to,iter]*temp_predict[j]*origin2 + 
                                                  # btemp0xorigin3_array_UCH[from,to,iter]*origin3 +
                                                  btemp1xorigin3_array_UCH[from,to,iter]*temp_predict[j]*origin3 +
                                                  borigin1_array_UCH[from,to,iter]*origin1 +
                                                  borigin2_array_UCH[from,to,iter]*origin2 +
                                                  borigin3_array_UCH[from,to,iter]*origin3)/
          sum(exp(b0_array_UCH[from,possible_movements,iter] +
                    # btemp0_array_UCH[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_UCH[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_UCH[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_UCH[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_UCH[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_UCH[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_UCH[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_UCH[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    # btemp0xorigin3_array_UCH[from,possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_UCH[from,possible_movements,iter]*temp_predict[j]*origin3 +
                    borigin1_array_UCH[from,possible_movements,iter]*origin1 +
                    borigin2_array_UCH[from,possible_movements,iter]*origin2 +
                    borigin3_array_UCH[from,possible_movements,iter]*origin3))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

estimate_temp_effect_MCW <- function(origin_select, movements){
  
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
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # If it's a mainstem state, get the spill data. If it's not, ignore it.
    
    if (from %in% c(1:9)){
      # get median winter spill for this state
      MCW_states_dates <- data.frame(state = as.vector(MCW_envir$data$y),
                                     date = as.vector(MCW_envir$data$transition_dates))
      spillwindow_data <- MCW_envir$data$spill_window_data
      med_spillwindow <- median(spillwindow_data[subset(MCW_states_dates, state == from)$date,from])
      
      # get median spill window for this state
      winterspill_data <- MCW_envir$data$winter_spill_days_data
      med_winterspill <- median(winterspill_data[,from])      
      
    } else {
      med_spillwindow <- 0
      med_winterspill <- 0
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j, iter,i] <- exp(b0_array_MCW[from,to,iter] +
                                                 # btemp0_array_MCW[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                 btemp1_array_MCW[from,to,iter]*temp_predict[j] + 
                                                 bspillwindow_array_MCW[from,to,iter]*med_spillwindow +
                                                 bwinterspill_array_MCW[from,to,iter]*med_winterspill +
                                                 # btemp0xorigin1_array_MCW[from,to,iter]*origin1 +
                                                 btemp1xorigin1_array_MCW[from,to,iter]*temp_predict[j]*origin1 +
                                                 # btemp0xorigin2_array_MCW[from,to,iter]*origin2 + 
                                                 btemp1xorigin2_array_MCW[from,to,iter]*temp_predict[j]*origin2 + 
                                                 # btemp0xorigin3_array_MCW[from,to,iter]*origin3 +
                                                 btemp1xorigin3_array_MCW[from,to,iter]*temp_predict[j]*origin3 +
                                                 # btemp0xorigin4_array_MCW[from,to,iter]*origin4 +
                                                 btemp1xorigin4_array_MCW[from,to,iter]*temp_predict[j]*origin4 +
                                                 # btemp0xorigin5_array_MCW[from,to,iter]*origin5 +
                                                 btemp1xorigin5_array_MCW[from,to,iter]*temp_predict[j]*origin5 +
                                                 # btemp0xorigin6_array_MCW[from,to,iter]*origin6 +
                                                 btemp1xorigin6_array_MCW[from,to,iter]*temp_predict[j]*origin6 +
                                                 borigin1_array_MCW[from,to,iter]*origin1 +
                                                 borigin2_array_MCW[from,to,iter]*origin2 +
                                                 borigin3_array_MCW[from,to,iter]*origin3 +
                                                 borigin4_array_MCW[from,to,iter]*origin4 +
                                                 borigin5_array_MCW[from,to,iter]*origin5 +
                                                 borigin6_array_MCW[from,to,iter]*origin6)/
          sum(exp(b0_array_MCW[from,possible_movements,iter] +
                    # btemp0_array_MCW[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_MCW[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_MCW[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_MCW[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_MCW[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_MCW[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    # btemp0xorigin3_array_MCW[from,possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin3 +
                    # btemp0xorigin4_array_MCW[from,possible_movements,iter]*origin4 +
                    btemp1xorigin4_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin4 +
                    # btemp0xorigin5_array_MCW[from,possible_movements,iter]*origin5 +
                    btemp1xorigin5_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin5 +
                    # btemp0xorigin6_array_MCW[from,possible_movements,iter]*origin6 +
                    btemp1xorigin6_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin6 +
                    borigin1_array_MCW[from,possible_movements,iter]*origin1 +
                    borigin2_array_MCW[from,possible_movements,iter]*origin2 +
                    borigin3_array_MCW[from,possible_movements,iter]*origin3 +
                    borigin4_array_MCW[from,possible_movements,iter]*origin4 +
                    borigin5_array_MCW[from,possible_movements,iter]*origin5 +
                    borigin6_array_MCW[from,possible_movements,iter]*origin6))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

estimate_temp_effect_MCH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,2)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # If it's a mainstem state, get the spill data. If it's not, ignore it.
    
    if (from %in% c(1:9)){
      # get median winter spill for this state
      MCH_states_dates <- data.frame(state = as.vector(MCH_envir$data$y),
                                     date = as.vector(MCH_envir$data$transition_dates))
      spillwindow_data <- MCH_envir$data$spill_window_data
      med_spillwindow <- median(spillwindow_data[subset(MCH_states_dates, state == from)$date,from])
      
      # get median spill window for this state
      winterspill_data <- MCH_envir$data$winter_spill_days_data
      med_winterspill <- median(winterspill_data[,from])      
      
    } else {
      med_spillwindow <- 0
      med_winterspill <- 0
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j,iter,i] <- exp(b0_array_MCH[from,to,iter] +
                                                # btemp0_array_MCH[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                btemp1_array_MCH[from,to,iter]*temp_predict[j] + 
                                                bspillwindow_array_MCH[from,to,iter]*med_spillwindow +
                                                bwinterspill_array_MCH[from,to,iter]*med_winterspill +
                                                # btemp0xorigin1_array_MCH[from,to,iter]*origin1 +
                                                btemp1xorigin1_array_MCH[from,to,iter]*temp_predict[j]*origin1 +
                                                # btemp0xorigin2_array_MCH[from,to,iter]*origin2 + 
                                                btemp1xorigin2_array_MCH[from,to,iter]*temp_predict[j]*origin2 +
                                                borigin1_array_MCH[from,to,iter]*origin1 +
                                                borigin2_array_MCH[from,to,iter]*origin2)/
          sum(exp(b0_array_MCH[from,possible_movements,iter] +
                    # btemp0_array_MCH[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_MCH[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_MCH[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_MCH[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_MCH[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_MCH[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_MCH[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_MCH[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    borigin1_array_MCH[from,possible_movements,iter]*origin1 +
                    borigin2_array_MCH[from,possible_movements,iter]*origin2))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

estimate_temp_effect_SRW <- function(origin_select, movements){
  
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
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # If it's a mainstem state, get the spill data. If it's not, ignore it.
    
    if (from %in% c(1:9)){
      # get median winter spill for this state
      SRW_states_dates <- data.frame(state = as.vector(SRW_envir$data$y),
                                     date = as.vector(SRW_envir$data$transition_dates))
      spillwindow_data <- SRW_envir$data$spill_window_data
      med_spillwindow <- median(spillwindow_data[subset(SRW_states_dates, state == from)$date,from])
      
      # get median spill window for this state
      winterspill_data <- SRW_envir$data$winter_spill_days_data
      med_winterspill <- median(winterspill_data[,from])      
      
    } else {
      med_spillwindow <- 0
      med_winterspill <- 0
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j,iter,i] <- exp(b0_array_SRW[from,to,iter] +
                                                # btemp0_array_SRW[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                btemp1_array_SRW[from,to,iter]*temp_predict[j] + 
                                                bspillwindow_array_SRW[from,to,iter]*med_spillwindow +
                                                bwinterspill_array_SRW[from,to,iter]*med_winterspill +
                                                # btemp0xorigin1_array_SRW[from,to,iter]*origin1 +
                                                btemp1xorigin1_array_SRW[from,to,iter]*temp_predict[j]*origin1 +
                                                # btemp0xorigin2_array_SRW[from,to,iter]*origin2 + 
                                                btemp1xorigin2_array_SRW[from,to,iter]*temp_predict[j]*origin2 + 
                                                # btemp0xorigin3_array_SRW[from,to,iter]*origin3 +
                                                btemp1xorigin3_array_SRW[from,to,iter]*temp_predict[j]*origin3 +
                                                # btemp0xorigin4_array_SRW[from,to,iter]*origin4 +
                                                btemp1xorigin4_array_SRW[from,to,iter]*temp_predict[j]*origin4 +
                                                # btemp0xorigin5_array_SRW[from,to,iter]*origin5 +
                                                btemp1xorigin5_array_SRW[from,to,iter]*temp_predict[j]*origin5 +
                                                # btemp0xorigin6_array_SRW[from,to,iter]*origin6 +
                                                btemp1xorigin6_array_SRW[from,to,iter]*temp_predict[j]*origin6 +
                                                borigin1_array_SRW[from,to,iter]*origin1 +
                                                borigin2_array_SRW[from,to,iter]*origin2 +
                                                borigin3_array_SRW[from,to,iter]*origin3 +
                                                borigin4_array_SRW[from,to,iter]*origin4 +
                                                borigin5_array_SRW[from,to,iter]*origin5 +
                                                borigin6_array_SRW[from,to,iter]*origin6)/
          sum(exp(b0_array_SRW[from,possible_movements,iter] +
                    # btemp0_array_SRW[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_SRW[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_SRW[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_SRW[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_SRW[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_SRW[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    # btemp0xorigin3_array_SRW[from,possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin3 +
                    # btemp0xorigin4_array_SRW[from,possible_movements,iter]*origin4 +
                    btemp1xorigin4_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin4 +
                    # btemp0xorigin5_array_SRW[from,possible_movements,iter]*origin5 +
                    btemp1xorigin5_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin5 +
                    # btemp0xorigin6_array_SRW[from,possible_movements,iter]*origin6 +
                    btemp1xorigin6_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin6 +
                    borigin1_array_SRW[from,possible_movements,iter]*origin1 +
                    borigin2_array_SRW[from,possible_movements,iter]*origin2 +
                    borigin3_array_SRW[from,possible_movements,iter]*origin3 +
                    borigin4_array_SRW[from,possible_movements,iter]*origin4 +
                    borigin5_array_SRW[from,possible_movements,iter]*origin5 +
                    borigin6_array_SRW[from,possible_movements,iter]*origin6))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

estimate_temp_effect_SRH <- function(origin_select, movements){
  
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
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # If it's a mainstem state, get the spill data. If it's not, ignore it.
    
    if (from %in% c(1:9)){
      # get median winter spill for this state
      SRH_states_dates <- data.frame(state = as.vector(SRH_envir$data$y),
                                     date = as.vector(SRH_envir$data$transition_dates))
      spillwindow_data <- SRH_envir$data$spill_window_data
      med_spillwindow <- median(spillwindow_data[subset(SRH_states_dates, state == from)$date,from])
      
      # get median spill window for this state
      winterspill_data <- SRH_envir$data$winter_spill_days_data
      med_winterspill <- median(winterspill_data[,from])      
      
    } else {
      med_spillwindow <- 0
      med_winterspill <- 0
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j, iter,i] <- exp(b0_array_SRH[from,to,iter] +
                                                 # btemp0_array_SRH[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                 btemp1_array_SRH[from,to,iter]*temp_predict[j] + 
                                                 bspillwindow_array_SRH[from,to,iter]*med_spillwindow +
                                                 bwinterspill_array_SRH[from,to,iter]*med_winterspill +
                                                 # btemp0xorigin1_array_SRH[from,to,iter]*origin1 +
                                                 btemp1xorigin1_array_SRH[from,to,iter]*temp_predict[j]*origin1 +
                                                 # btemp0xorigin2_array_SRH[from,to,iter]*origin2 + 
                                                 btemp1xorigin2_array_SRH[from,to,iter]*temp_predict[j]*origin2 + 
                                                 # btemp0xorigin3_array_SRH[from,to,iter]*origin3 +
                                                 btemp1xorigin3_array_SRH[from,to,iter]*temp_predict[j]*origin3 +
                                                 # btemp0xorigin4_array_SRH[from,to,iter]*origin4 +
                                                 btemp1xorigin4_array_SRH[from,to,iter]*temp_predict[j]*origin4 +
                                                 # btemp0xorigin5_array_SRH[from,to,iter]*origin5 +
                                                 btemp1xorigin5_array_SRH[from,to,iter]*temp_predict[j]*origin5 +
                                                 borigin1_array_SRH[from,to,iter]*origin1 +
                                                 borigin2_array_SRH[from,to,iter]*origin2 +
                                                 borigin3_array_SRH[from,to,iter]*origin3 +
                                                 borigin4_array_SRH[from,to,iter]*origin4 +
                                                 borigin5_array_SRH[from,to,iter]*origin5)/
          sum(exp(b0_array_SRH[from,possible_movements,iter] +
                    # btemp0_array_SRH[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_SRH[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_SRH[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_SRH[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_SRH[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_SRH[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_SRH[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_SRH[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    # btemp0xorigin3_array_SRH[from,possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_SRH[from,possible_movements,iter]*temp_predict[j]*origin3 +
                    # btemp0xorigin4_array_SRH[from,possible_movements,iter]*origin4 +
                    btemp1xorigin4_array_SRH[from,possible_movements,iter]*temp_predict[j]*origin4 +
                    # btemp0xorigin5_array_SRH[from,possible_movements,iter]*origin5 +
                    btemp1xorigin5_array_SRH[from,possible_movements,iter]*temp_predict[j]*origin5 +
                    borigin1_array_SRH[from,possible_movements,iter]*origin1 +
                    borigin2_array_SRH[from,possible_movements,iter]*origin2 +
                    borigin3_array_SRH[from,possible_movements,iter]*origin3 +
                    borigin4_array_SRH[from,possible_movements,iter]*origin4 +
                    borigin5_array_SRH[from,possible_movements,iter]*origin5))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

#### Plot temperature effect on multiple movements on same plot ####

# For some origins, there are multiple interesting and deleterious movements - 
# for example, Middle Columbia fish could go to loss, Deschutes, or overshoot

# plotting function
plot_compare_rear_temp_effect_multiple_movements <- function(origin_select,
                                                             wild_move_prob_array = NULL, hatchery_move_prob_array = NULL,
                                                             wild_covariate_experiences = NULL, hatchery_covariate_experiences = NULL,
                                                             movements_evaluated,
                                                             plot_title = NULL, plot_legend = FALSE){
  
  # create df to index to right dam temp joining; here we need to change names of 
  # WEL and LGR because we have slow v fast there
  dam_index <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL_quick", "ICH", "LGR_quick"),
                          state = seq(2,9))
  
  niter <- 4000 # this is the number of draws we have
  
  # Initialize a df to store every movement, then loop through each movement and bind
  # rear_temp_move_prob_quantiles <- data.frame(temp_actual = NA, `0.025` = NA, `0.5` = NA, `0.975` = NA,
  #                                             rear = NA, from = NA, to = NA)
  rear_temp_move_prob_quantiles <- data.frame()
  for (i in 1:nrow(movements_evaluated)){
    # First, determine if this origin has both a hatchery and a wild population
    
    # If hatchery is NA, run wild only
    if (is.null(hatchery_covariate_experiences)){
      wild_temp_move_prob <- as.data.frame(wild_move_prob_array[,,i])
      
      colnames(wild_temp_move_prob) <- paste0("iter", 1:niter) 
      temp_predict <- seq(-2,2,length = 100)
      wild_temp_move_prob$temp <- temp_predict
      
      # Add a column with the actual temperatures
      if (movements_evaluated$from[i] %in% c(1:9)){
        wild_temp_move_prob %>% 
          mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_mean")] + 
                   window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_sd")]*temp) -> wild_temp_move_prob
      } else {
        wild_temp_move_prob %>% 
          mutate(temp_actual = 1) -> wild_temp_move_prob
      }
      
      
      
      
      wild_temp_move_prob %>% 
        pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
        group_by(temp_actual) %>% 
        dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
        pivot_wider(names_from = q, values_from = prob) %>% 
        mutate(rear = "wild") %>% 
        mutate(from = movements_evaluated$from[i],
               to = movements_evaluated$to[i]) -> wild_temp_move_prob_quantiles
      
      # get data organized for rug plot
      wild_covariate_experiences %>% 
        mutate(rear = "wild") -> wild_covariate_experiences
      
      wild_covariate_experiences %>% 
        # keep only the state that is our from state
        filter(state == movements_evaluated$from[i]) -> covariate_experiences
      
      rear_temp_move_prob_quantiles %>% 
        bind_rows(., wild_temp_move_prob_quantiles) -> rear_temp_move_prob_quantiles
      
    } 
    # If wild is NA, run hatchery only
    else if (is.null(wild_covariate_experiences)){
      hatchery_temp_move_prob <- as.data.frame(hatchery_move_prob_array[,,i])
      
      colnames(hatchery_temp_move_prob) <- paste0("iter", 1:niter) 
      temp_predict <- seq(-2,2,length = 100)
      hatchery_temp_move_prob$temp <- temp_predict
      
      # Add a column with the actual temperatures
      if (movements_evaluated$from[i] %in% c(1:9)){
        hatchery_temp_move_prob %>% 
          mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_mean")] + 
                   window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_sd")]*temp) -> hatchery_temp_move_prob
      } else {
        hatchery_temp_move_prob %>% 
          mutate(temp_actual = 1) -> hatchery_temp_move_prob
      }
      
      
      
      
      hatchery_temp_move_prob %>% 
        pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
        group_by(temp_actual) %>% 
        summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
        pivot_wider(names_from = q, values_from = prob) %>% 
        mutate(rear = "hatchery") %>% 
        mutate(from = movements_evaluated$from[i],
               to = movements_evaluated$to[i]) -> hatchery_temp_move_prob_quantiles
      
      # get data organized for rug plot
      hatchery_covariate_experiences %>% 
        mutate(rear = "hatchery") -> hatchery_covariate_experiences
      
      hatchery_covariate_experiences %>% 
        # keep only the state that is our from state
        filter(state == movements_evaluated$from[i]) -> covariate_experiences
      
      # combine wild and hatchery
      rear_temp_move_prob_quantiles %>% 
        bind_rows(., hatchery_temp_move_prob_quantiles) -> rear_temp_move_prob_quantiles
      
    } 
    # else run both
    else {
      wild_temp_move_prob <- as.data.frame(wild_move_prob_array[,,i])
      
      colnames(wild_temp_move_prob) <- paste0("iter", 1:niter) 
      temp_predict <- seq(-2,2,length = 100)
      wild_temp_move_prob$temp <- temp_predict
      
      # Add a column with the actual temperatures
      if (movements_evaluated$from[i] %in% c(1:9)){
        wild_temp_move_prob %>% 
          mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_mean")] + 
                   window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_sd")]*temp) -> wild_temp_move_prob
      } else {
        wild_temp_move_prob %>% 
          mutate(temp_actual = 1) -> wild_temp_move_prob
      }
      
      
      wild_temp_move_prob %>% 
        pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
        group_by(temp_actual) %>% 
        summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
        pivot_wider(names_from = q, values_from = prob) %>% 
        mutate(rear = "wild") %>% 
        mutate(from = movements_evaluated$from[i],
               to = movements_evaluated$to[i]) -> wild_temp_move_prob_quantiles
      
      # get data organized for rug plot
      wild_covariate_experiences %>% 
        mutate(rear = "wild") -> wild_covariate_experiences
      
      wild_covariate_experiences %>% 
        # keep only the state that is our from state
        filter(state == movements_evaluated$from[i]) -> wild_covariate_experiences
      
      hatchery_temp_move_prob <- as.data.frame(hatchery_move_prob_array[,,i])
      
      colnames(hatchery_temp_move_prob) <- paste0("iter", 1:niter) 
      hatchery_temp_move_prob$temp <- temp_predict
      
      # Add a column with the actual temperatures
      if (movements_evaluated$from[i] %in% c(1:9)){
        hatchery_temp_move_prob %>% 
          mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_mean")] + 
                   window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_sd")]*temp) -> hatchery_temp_move_prob
      } else {
        hatchery_temp_move_prob %>% 
          mutate(temp_actual = 1) -> hatchery_temp_move_prob
      }
      
      hatchery_temp_move_prob %>% 
        pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
        group_by(temp_actual) %>% 
        summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
        pivot_wider(names_from = q, values_from = prob) %>% 
        mutate(rear = "hatchery") %>% 
        mutate(from = movements_evaluated$from[i],
               to = movements_evaluated$to[i]) -> hatchery_temp_move_prob_quantiles
      
      # get data organized for rug plot
      hatchery_covariate_experiences %>% 
        mutate(rear = "hatchery") -> hatchery_covariate_experiences
      
      hatchery_covariate_experiences %>% 
        # keep only the state that is our from state
        filter(state == movements_evaluated$from[i]) -> hatchery_covariate_experiences
      
      # combine wild and hatchery
      rear_temp_move_prob_quantiles %>% 
        bind_rows(., wild_temp_move_prob_quantiles) %>% 
        bind_rows(., hatchery_temp_move_prob_quantiles) -> rear_temp_move_prob_quantiles
      
      wild_covariate_experiences %>% 
        bind_rows(., hatchery_covariate_experiences) -> covariate_experiences
      
    }
    
  }
  
  
  
  
  # convert temp to temp_actual in covariate experiences
  if (movements_evaluated$from[i] %in% c(1:9)){
    covariate_experiences %>% 
      mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_mean")] + 
               window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_sd")]*temperature) -> covariate_experiences
  } else {
    covariate_experiences %>% 
      mutate(temp_actual = 1) -> covariate_experiences
  }
  
  
  # Remove any instances for rug plot of fish that would have received temp0 covariate
  covariate_experiences %>% 
    dplyr::rename(date_numeric = date) %>%
    # keep only jan/feb/mar 
    mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
    mutate(month = month(date)) %>% 
    filter(!(month %in% c(1,2,3,4,5))) -> covariate_experiences
  
  # rear_temp_move_prob_quantiles$to <- as.character(rear_temp_move_prob_quantiles$to)
  
  
  movement_colors <- c("Overshoot" = "#ff7f00", "Overshoot - ICH" = "#ff7f00",
                       "Overshoot - PRA" = "#ff7f00",
                       "Deschutes River" = "#a6cee3",
                       "Home" = "#1f78b4", "Loss" = "#e31a1c")
  
  rear_temp_move_prob_quantiles %>% 
    left_join(., movements_evaluated, by = "to") -> rear_temp_move_prob_quantiles
  
  if (plot_legend == TRUE){
    
    combined_plot <- ggplot(rear_temp_move_prob_quantiles, aes(x = temp_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                               color = Movement, fill = Movement)) +
      geom_line(linewidth = 2.5) +
      geom_ribbon(alpha = 0.2, color = NA) +
      # geom_rug(data = covariate_experiences, aes(x = temp_actual), inherit.aes = FALSE,
      #          sides = "t", length = unit(0.3, "cm"), outside = TRUE) +
      scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
      # scale_x_continuous(lim = c(0,ceiling(max(covariate_experiences$temp_actual))), expand = c(0,0)) +
      # common x-axis scale across all populations
      scale_x_continuous(lim = c(0, 22.5), expand = c(0,0)) + 
      scale_color_manual(values = movement_colors) +
      scale_fill_manual(values =  movement_colors) +
      xlab(expression(~"Temperature" ~ ("°C"))) +
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
    temp_legend <- ggpubr::get_legend(combined_plot)
    # temp_legend$widths[[1]] <- unit(0, "cm")
    # temp_legend$widths[[5]] <- unit(0, "cm")
    # temp_legend$heights[[1]] <- unit(0, "cm")
    # temp_legend$heights[[5]] <- unit(0, "cm")
    temp_plot_legend_gg <- as_ggplot(temp_legend)
    
    # for testing
    # ggsave(here::here("figures", "paper_figures", "01_legend_test.png"), temp_plot_legend_gg, height = 4, width = 4)
    
  } else {
    # suppress common legend - for combined plot
    rear_temp_move_prob_plot <- ggplot(rear_temp_move_prob_quantiles, aes(x = temp_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                          color = Movement, fill = Movement)) +
      geom_line() +
      geom_ribbon(alpha = 0.2, color = NA) +
      # geom_rug(data = covariate_experiences, aes(x = temp_actual), inherit.aes = FALSE,
      #          sides = "t", length = unit(0.3, "cm"), outside = TRUE) +
      scale_y_continuous(lim = c(0,1)) +
      # common x-axis scale across all populations
      coord_cartesian(xlim = c(4, 22.5), expand = FALSE) +
      # scale_x_continuous(lim = c(floor(min(covariate_experiences$temp_actual)),ceiling(max(covariate_experiences$temp_actual))), expand = c(0,0)) +
      scale_color_manual(values = movement_colors) +
      scale_fill_manual(values =  movement_colors) +
      xlab(expression(~"Temperature" ~ ("°C"))) +
      ylab("Movement probability") +
      theme(legend.position = "none",
            panel.grid.major = element_line(color = "gray90"),
            panel.background = element_rect(fill = "white", color = "black"),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=0.4),
            # turn off the axis titles on each individual plot and just show one for whole plot
            axis.title = element_blank(),
            # these plot margins are to leave space for the population name on the big figure
            plot.margin = unit(c(0, 0.2, 0.2, 0.2),"cm"))
    
    density_plot <- ggplot(data = covariate_experiences, aes(temp_actual))+
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
      # coord_cartesian(xlim = c(floor(min(covariate_experiences$temp_actual)),ceiling(max(covariate_experiences$temp_actual)))) +
      coord_cartesian(xlim = c(4, 22.5), expand = FALSE) +
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
            # we actually don't need this anymore I think, because of the histogram
            plot.margin = unit(c(0.2, 0.2, 0, 0.2),"cm"))
    # testing theme
    # theme(
    #       # these plot margins are to leave space for the population name on the big figure
    #       plot.margin = unit(c(0.2, 0.2, 0, 0.2),"cm"))
    
    combined_plot <- ggarrange(density_plot, rear_temp_move_prob_plot, nrow = 2, ncol = 1,
                               heights = c(2,6))
    
    # for testing
    # ggsave(here::here("figures", "paper_figures", "01_test.png"), combined_plot, height = 6, width = 6)
    
    
  }
  
  
  
  return(combined_plot)
}


#### Prepare data to create figures for individual origins: ####
# John Day River (W), Umatilla River (H, W), Yakima River (W), Walla Walla River (H, W),
# Wenatchee River (H, W), Entiat River (W), Tucannon River (H, W), Imnaha River (H, W)

### John Day River ###
JDR_movements <- data.frame(from = c(2, 2, 2, 2), to = c(3, 10, 12, 41),
                            Movement = c("Overshoot", "Deschutes River",
                                         "Home", "Loss"))
JDR_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "John Day River", movements = JDR_movements)
JDR_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "John Day River")

### Umatilla River ###
UMA_movements <- data.frame(from = c(2,2,2,2), to = c(3, 10, 16, 41),
                            Movement = c("Overshoot", "Deschutes River",
                                         "Home", "Loss"))
UMA_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Umatilla River", movements = UMA_movements)
UMA_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Umatilla River", movements = UMA_movements)
UMA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Umatilla River")
UMA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Umatilla River")

### Yakima River ###
YAK_movements <- data.frame(from = c(3,3,3,3), to = c(8, 4, 18, 41),
                            Movement = c("Overshoot - ICH", "Overshoot - PRA",
                                         "Home", "Loss"))

YAK_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Yakima River", movements = YAK_movements)
YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Yakima River")

### Walla Walla River ###
WAWA_movements <- data.frame(from = c(3,3,3,3), to = c(8, 4, 20, 41),
                             Movement = c("Overshoot - ICH", "Overshoot - PRA",
                                          "Home", "Loss"))
WAWA_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Walla Walla River", movements = WAWA_movements)
WAWA_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Walla Walla River", movements = WAWA_movements)
WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Walla Walla River")
WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Walla Walla River")

### Wenatchee River ###
WEN_movements <- data.frame(from = c(5,5,5), to = c(6, 22, 41),
                            Movement = c("Overshoot",
                                         "Home", "Loss"))
WEN_wild_temp_move_prob_array <- estimate_temp_effect_UCW(origin_select = "Wenatchee River", movements = WEN_movements)
WEN_hatchery_temp_move_prob_array <- estimate_temp_effect_UCH(origin_select = "Wenatchee River", movements = WEN_movements)
WEN_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Wenatchee River")
WEN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Wenatchee River")

### Entiat River ###
ENT_movements <- data.frame(from = c(6,6,6), to = c(7, 24, 41),
                            Movement = c("Overshoot",
                                         "Home", "Loss"))
ENT_wild_temp_move_prob_array <- estimate_temp_effect_UCW(origin_select = "Entiat River", movements = ENT_movements)
ENT_hatchery_temp_move_prob_array <- estimate_temp_effect_UCH(origin_select = "Entiat River", movements = ENT_movements)
ENT_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Entiat River")
ENT_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Entiat River")

### Tucannon River ###
TUC_movements <- data.frame(from = c(8,8,8), to = c(9, 30, 41),
                            Movement = c("Overshoot", 
                                         "Home", "Loss"))
TUC_wild_temp_move_prob_array <- estimate_temp_effect_SRW(origin_select = "Tucannon River", movements = TUC_movements)
TUC_hatchery_temp_move_prob_array <- estimate_temp_effect_SRH(origin_select = "Tucannon River", movements = TUC_movements)
TUC_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Tucannon River")
TUC_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Tucannon River")

### Imnaha River ###
IMN_movements <- data.frame(from = c(9,9), to = c(37, 41),
                            Movement = c("Home", "Loss"))
IMN_wild_temp_move_prob_array <- estimate_temp_effect_SRW(origin_select = "Imnaha River", movements = IMN_movements)
IMN_hatchery_temp_move_prob_array <- estimate_temp_effect_SRH(origin_select = "Imnaha River", movements = IMN_movements)
IMN_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Imnaha River")
IMN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Imnaha River")




#### Create figures for individual origins: ####
# John Day River (W), Umatilla River (H, W), Yakima River (W), Walla Walla River (H, W),
# Wenatchee River (H, W), Entiat River (W), Tucannon River (H, W), Imnaha River (H, W)
### John Day River ###
JDR_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "John Day River",
                                                                                   wild_move_prob_array = JDR_wild_temp_move_prob_array,
                                                                                   
                                                                                   wild_covariate_experiences = JDR_wild_covariate_experiences,
                                                                                   
                                                                                   movements_evaluated = JDR_movements)

ggsave(here::here("figures", "temperature_effects", "JDR_wild_compare_movement_temp.png"), JDR_wild_compare_movement_temp, height = 8, width = 8)

### Umatilla River ###

# UMA wild plot
UMA_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Umatilla River",
                                                                                   wild_move_prob_array = UMA_wild_temp_move_prob_array,
                                                                                   wild_covariate_experiences = UMA_wild_covariate_experiences,
                                                                                   movements_evaluated = UMA_movements)

ggsave(here::here("figures", "temperature_effects", "UMA_wild_compare_movement_temp.png"), UMA_wild_compare_movement_temp, height = 8, width = 8)

# UMA hatchery plot
UMA_hatchery_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Umatilla River",
                                                                                       hatchery_move_prob_array = UMA_hatchery_temp_move_prob_array,
                                                                                       hatchery_covariate_experiences = UMA_hatchery_covariate_experiences,
                                                                                       movements_evaluated = UMA_movements)

ggsave(here::here("figures", "temperature_effects", "UMA_hatchery_compare_movement_temp.png"), UMA_hatchery_compare_movement_temp, height = 8, width = 8)

### Yakima River ###

YAK_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Yakima River",
                                                                                   wild_move_prob_array = YAK_wild_temp_move_prob_array,
                                                                                   
                                                                                   wild_covariate_experiences = YAK_wild_covariate_experiences,
                                                                                   
                                                                                   movements_evaluated = YAK_movements)

ggsave(here::here("figures", "temperature_effects", "YAK_wild_compare_movement_temp.png"), YAK_wild_compare_movement_temp, height = 8, width = 8)



### Walla Walla River ###

# WAWA wild plot
WAWA_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Walla Walla River",
                                                                                    wild_move_prob_array = WAWA_wild_temp_move_prob_array,
                                                                                    wild_covariate_experiences = WAWA_wild_covariate_experiences,
                                                                                    movements_evaluated = WAWA_movements)

ggsave(here::here("figures", "temperature_effects", "WAWA_wild_compare_movement_temp.png"), WAWA_wild_compare_movement_temp, height = 8, width = 8)

# WAWA hatchery plot
WAWA_hatchery_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Walla Walla River",
                                                                                        hatchery_move_prob_array = WAWA_hatchery_temp_move_prob_array,
                                                                                        hatchery_covariate_experiences = WAWA_hatchery_covariate_experiences,
                                                                                        movements_evaluated = WAWA_movements)

ggsave(here::here("figures", "temperature_effects", "WAWA_hatchery_compare_movement_temp.png"), WAWA_hatchery_compare_movement_temp, height = 8, width = 8)

### Wenatchee River ###

# WEN wild plot
WEN_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Wenatchee River",
                                                                                   wild_move_prob_array = WEN_wild_temp_move_prob_array,
                                                                                   wild_covariate_experiences = WEN_wild_covariate_experiences,
                                                                                   movements_evaluated = WEN_movements)

ggsave(here::here("figures", "temperature_effects", "WEN_wild_compare_movement_temp.png"), WEN_wild_compare_movement_temp, height = 8, width = 8)

# WEN hatchery plot
WEN_hatchery_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Wenatchee River",
                                                                                       hatchery_move_prob_array = WEN_hatchery_temp_move_prob_array,
                                                                                       hatchery_covariate_experiences = WEN_hatchery_covariate_experiences,
                                                                                       movements_evaluated = WEN_movements)

ggsave(here::here("figures", "temperature_effects", "WEN_hatchery_compare_movement_temp.png"), WEN_hatchery_compare_movement_temp, height = 8, width = 8)

### Entiat River ###

# ENT wild plot
ENT_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Entiat River",
                                                                                   wild_move_prob_array = ENT_wild_temp_move_prob_array,
                                                                                   wild_covariate_experiences = ENT_wild_covariate_experiences,
                                                                                   movements_evaluated = ENT_movements)

ggsave(here::here("figures", "temperature_effects", "ENT_wild_compare_movement_temp.png"), ENT_wild_compare_movement_temp, height = 8, width = 8)


### Tucannon River ###

# TUC wild plot
TUC_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Tucannon River",
                                                                                   wild_move_prob_array = TUC_wild_temp_move_prob_array,
                                                                                   wild_covariate_experiences = TUC_wild_covariate_experiences,
                                                                                   movements_evaluated = TUC_movements)

ggsave(here::here("figures", "temperature_effects", "TUC_wild_compare_movement_temp.png"), TUC_wild_compare_movement_temp, height = 8, width = 8)

# TUC hatchery plot
TUC_hatchery_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Tucannon River",
                                                                                       hatchery_move_prob_array = TUC_hatchery_temp_move_prob_array,
                                                                                       hatchery_covariate_experiences = TUC_hatchery_covariate_experiences,
                                                                                       movements_evaluated = TUC_movements)

ggsave(here::here("figures", "temperature_effects", "TUC_hatchery_compare_movement_temp.png"), TUC_hatchery_compare_movement_temp, height = 8, width = 8)

### Imnaha River ###

# IMN wild plot
IMN_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Imnaha River",
                                                                                   wild_move_prob_array = IMN_wild_temp_move_prob_array,
                                                                                   wild_covariate_experiences = IMN_wild_covariate_experiences,
                                                                                   movements_evaluated = IMN_movements)

ggsave(here::here("figures", "temperature_effects", "IMN_wild_compare_movement_temp.png"), IMN_wild_compare_movement_temp, height = 8, width = 8)

# IMN hatchery plot
IMN_hatchery_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Imnaha River",
                                                                                       hatchery_move_prob_array = IMN_hatchery_temp_move_prob_array,
                                                                                       hatchery_covariate_experiences = IMN_hatchery_covariate_experiences,
                                                                                       movements_evaluated = IMN_movements)

ggsave(here::here("figures", "temperature_effects", "IMN_hatchery_compare_movement_temp.png"), IMN_hatchery_compare_movement_temp, height = 8, width = 8)


# Create the legend figure by itself
temp_plot_for_legend <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "John Day River",
                                                                         wild_move_prob_array = JDR_wild_temp_move_prob_array,
                                                                         
                                                                         wild_covariate_experiences = JDR_wild_covariate_experiences,
                                                                         
                                                                         movements_evaluated = JDR_movements,
                                                                         plot_legend= TRUE)
temp_legend <- ggpubr::get_legend(temp_plot_for_legend)
temp_plot_legend_gg <- as_ggplot(temp_legend) + theme(panel.background = element_rect(fill = "white", color = "white"))


#### Generate Figure 4 using the figures above ####

combined_movement_temp_plot <- ggarrange(JDR_wild_compare_movement_temp, UMA_wild_compare_movement_temp, UMA_hatchery_compare_movement_temp,
                                         YAK_wild_compare_movement_temp, WAWA_wild_compare_movement_temp, WAWA_hatchery_compare_movement_temp,
                                         WEN_wild_compare_movement_temp, WEN_hatchery_compare_movement_temp, ENT_wild_compare_movement_temp,
                                         TUC_wild_compare_movement_temp, TUC_hatchery_compare_movement_temp,
                                         temp_plot_legend_gg, nrow = 3, ncol = 4,
                                         labels = c("(A) JDR, Natural", "(B) UMA, Natural", 
                                                    "(C) UMA, Hatchery", "(D) YAK, Natural", 
                                                    "(E) WAWA, Natural", "(F) WAWA, Hatchery", 
                                                    "(G) WEN, Natural", "(H) WEN, Hatchery", 
                                                    "(I) ENT, Natural",
                                                    "(J) TUC, Natural", "(K) TUC, Hatchery"),
                                         label.x = 0.05, label.y = 0.925, font.label = list(size = 14, face = "plain"),
                                         hjust = 0, vjust = 0)

# combined_movement_temp_plot <- annotate_figure(combined_movement_temp_plot,
#                                                bottom = textGrob(expression(~"Temperature" ~ ("°C")), gp = gpar(cex = 1.3)),
#                                                  left = textGrob("Movement probability", rot = 90, gp = gpar(cex = 1.3))) + bgcolor("white")

# Let's try this again, this time using a cowplot solution since ggpubr is
# struggling with the background color
combined_movement_temp_plot <- cowplot::ggdraw(annotate_figure(combined_movement_temp_plot,
                                                               bottom = textGrob(expression(~"Temperature" ~ ("°C")), gp = gpar(cex = 1.3)),
                                                               left = textGrob("Movement probability", rot = 90, gp = gpar(cex = 1.3)))) +
  theme(plot.background = element_rect(fill="white", color = NA))


ggsave(here::here("figures", "paper_figures", "Fig4_temperature_home_mainstem.png"), combined_movement_temp_plot, height = 12, width = 16)

#### Figure 5: Deschutes river use vs. temperature ####

### Step 1: Evaluate probability of movement into Deschutes for different origins ###

# establish df for movements into Deschutes only
des_movement <- data.frame(from = c(2), to = c(10),
                            Movement = c("Deschutes River"))

# estimate probability of movements into Deschutes across all origins

# Middle Columbia, Wild
JDR_wild_temp_DES_move_prob_array <- estimate_temp_effect_MCW(origin_select = "John Day River", movements = des_movement)
UMA_wild_temp_DES_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Umatilla River", movements = des_movement)
YAK_wild_temp_DES_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Yakima River", movements = des_movement)
WAWA_wild_temp_DES_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Walla Walla River", movements = des_movement)
FIF_wild_temp_DES_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Fifteenmile Creek", movements = des_movement)
DES_wild_temp_DES_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Deschutes River", movements = des_movement)

# Middle Columbia, Hatchery
UMA_hatchery_temp_DES_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Umatilla River", movements = des_movement)
WAWA_hatchery_temp_DES_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Walla Walla River", movements = des_movement)

# Upper Columbia, wild
# Note that since this movement is outside the DPS boundaries, it'll be shared across all origins
# so we can just choose one and it'll be the same for all
UC_wild_temp_DES_move_prob_array <- estimate_temp_effect_UCW(origin_select = "Wenatchee River", movements = des_movement)

# Upper Columbia, hatchery
# Note that since this movement is outside the DPS boundaries, it'll be shared across all origins
# so we can just choose one and it'll be the same for all
UC_hatchery_temp_DES_move_prob_array <- estimate_temp_effect_UCH(origin_select = "Wenatchee River", movements = des_movement)

# Snake River, wild
# Note that since this movement is outside the DPS boundaries, it'll be shared across all origins
# so we can just choose one and it'll be the same for all
SR_wild_temp_DES_move_prob_array <- estimate_temp_effect_SRW(origin_select = "Tucannon River", movements = des_movement)

# Snake River, hatchery
# Note that since this movement is outside the DPS boundaries, it'll be shared across all origins
# so we can just choose one and it'll be the same for all
SR_hatchery_temp_DES_move_prob_array <- estimate_temp_effect_SRH(origin_select = "Tucannon River", movements = des_movement)


### Step 2: Reformat predicted movement probabilities for plot ###
# little function to reformat movement prob arrays
reformat_DES_prob_array <- function(move_prob_array, origin_name, rear_type){
  movement_df <- as.data.frame(move_prob_array[,,1])
  colnames(movement_df) <- 1:ncol(movement_df)
  temp_predict <- seq(-2,2,length = 100)
  movement_df$temp <- temp_predict
  
  movement_df %>% 
    pivot_longer(!temp, names_to = "iter", values_to = "prob") %>% 
    mutate(origin = origin_name,
           rear = rear_type) -> movement_df
  
  dam_index <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL", "ICH", "LGR"),
                          state = seq(2,9))
  
  movement_df %>% 
    mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == 2)$dam, "_mean")] + 
           window_temp_summary[, paste0(subset(dam_index, state == 2)$dam, "_sd")]*temp) -> movement_df
  
  return(movement_df)
}

# reformat data for each origin
# Middle Columbia, Wild
JDR_natural_DES_move_long <- reformat_DES_prob_array(JDR_wild_temp_DES_move_prob_array, origin_name = "John Day River", rear_type = "natural")
UMA_natural_DES_move_long <- reformat_DES_prob_array(UMA_wild_temp_DES_move_prob_array, origin_name = "Umatilla River", rear_type = "natural")
YAK_natural_DES_move_long <- reformat_DES_prob_array(YAK_wild_temp_DES_move_prob_array, origin_name = "Yakima River", rear_type = "natural")
WAWA_natural_DES_move_long <- reformat_DES_prob_array(WAWA_wild_temp_DES_move_prob_array, origin_name = "Walla Walla River", rear_type = "natural")
FIF_natural_DES_move_long <- reformat_DES_prob_array(FIF_wild_temp_DES_move_prob_array, origin_name = "Fifteenmile Creek", rear_type = "natural")
DES_natural_DES_move_long <- reformat_DES_prob_array(DES_wild_temp_DES_move_prob_array, origin_name = "Deschutes River", rear_type = "natural")

# Middle Columbia, Hatchery
UMA_hatchery_DES_move_long <- reformat_DES_prob_array(UMA_hatchery_temp_DES_move_prob_array, origin_name = "Umatilla River", rear_type = "hatchery")
WAWA_hatchery_DES_move_long <- reformat_DES_prob_array(WAWA_hatchery_temp_DES_move_prob_array, origin_name = "Walla Walla River", rear_type = "hatchery")

# Upper Columbia, Wild
UC_wild_DES_move_long <- reformat_DES_prob_array(UC_wild_temp_DES_move_prob_array, origin_name = "Upper Columbia", rear_type = "natural")

# Upper Columbia, Hatchery
UC_hatchery_DES_move_long <- reformat_DES_prob_array(UC_hatchery_temp_DES_move_prob_array, origin_name = "Upper Columbia", rear_type = "hatchery")

# Snake River, Wild
SR_wild_DES_move_long <- reformat_DES_prob_array(SR_wild_temp_DES_move_prob_array, origin_name = "Snake River", rear_type = "natural")

# Snake River, Hatchery
SR_hatchery_DES_move_long <- reformat_DES_prob_array(SR_hatchery_temp_DES_move_prob_array, origin_name = "Snake River", rear_type = "hatchery")

JDR_natural_DES_move_long %>% 
  bind_rows(., UMA_natural_DES_move_long, YAK_natural_DES_move_long, WAWA_natural_DES_move_long, 
            FIF_natural_DES_move_long, DES_natural_DES_move_long,
            UMA_hatchery_DES_move_long, WAWA_hatchery_DES_move_long,
            UC_wild_DES_move_long, UC_hatchery_DES_move_long,
            SR_wild_DES_move_long, SR_hatchery_DES_move_long,) -> combined_DES_move_long

# summarise into median and 95% CI
combined_DES_move_long %>% 
  group_by(temp_actual, origin, rear) %>% 
  summarise(prob = quantile(prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
  pivot_wider(names_from = q, values_from = prob) %>% 
  mutate(origin_rear = paste0(origin, ", ", rear)) -> combined_DES_move_summary

# Drop Fifteenmile Creek and Deschutes
# Deschutes is the home tributary so it's not what we're interested in; Fifteenmile
# Creek is before Deschutes to it's of a different nature than the other tributaries
combined_DES_move_summary %>% 
  filter(!(origin %in% c("Fifteenmile Creek", "Deschutes River"))) -> combined_DES_move_summary

# Reorder these for plotting
combined_DES_move_summary$origin_rear <- factor(combined_DES_move_summary$origin_rear, levels = 
                                                  c("John Day River, natural",
                                                    "Umatilla River, hatchery",    
                                                    "Umatilla River, natural",
                                                    "Walla Walla River, hatchery",
                                                    "Walla Walla River, natural",  
                                                    "Yakima River, natural",
                                                    "Upper Columbia, hatchery",    
                                                    "Upper Columbia, natural",
                                                    "Snake River, hatchery",       
                                                    "Snake River, natural"      
                                                    ))

combined_DES_move_summary$origin <- factor(combined_DES_move_summary$origin, levels = 
                                                  c("John Day River",
                                                    "Umatilla River",    
                                                    "Walla Walla River",  
                                                    "Yakima River",
                                                    "Upper Columbia",
                                                    "Snake River"
                                                  ))


### Step 3: Generate and export the plot ###

DES_temp_move_prob_plot <- ggplot(combined_DES_move_summary, aes(x = temp_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                      color = origin_rear, fill = origin_rear)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, color = NA) +
  scale_y_continuous(lim = c(0,1)) +
  # common x-axis scale across all populations
  coord_cartesian(xlim = c(4, 22.5), expand = FALSE) +
  # scale_x_continuous(lim = c(floor(min(covariate_experiences$temp_actual)),ceiling(max(covariate_experiences$temp_actual))), expand = c(0,0)) +
  # scale_color_manual(values = movement_colors) +
  # scale_fill_manual(values =  movement_colors) +
  xlab(expression(~"Temperature" ~ ("°C"))) +
  ylab("Movement probability") +
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.4),
        # turn off the axis titles on each individual plot and just show one for whole plot
        axis.title = element_blank(),
        # these plot margins are to leave space for the population name on the big figure
        plot.margin = unit(c(0, 0.2, 0.2, 0.2),"cm"))


ggsave(here::here("figures", "paper_figures", "Fig5_deschutes_temp_movements.png"), DES_temp_move_prob_plot, height = 8, width = 10)

# Version 2 - without the credible intervals, to make it easier to interpret
# Set MC to one
# Set SR and UC to a different one
origin_rear_linetypes <- c("John Day River, natural" = 1,  
                           "Snake River, hatchery" = 2,       
                           "Snake River, natural" = 2,        
                           "Umatilla River, hatchery" = 1,    
                           "Umatilla River, natural" = 1,
                           "Upper Columbia, hatchery" = 2,    
                           "Upper Columbia, natural" = 2,     
                           "Walla Walla River, hatchery" = 1,
                           "Walla Walla River, natural" = 1,  
                           "Yakima River, natural" = 1)


DES_temp_move_prob_plot <- ggplot(combined_DES_move_summary, aes(x = temp_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                 color = origin, fill = origin, linetype = rear)) +
  geom_line() +
  # geom_ribbon(alpha = 0.2, color = NA) +
  scale_y_continuous(lim = c(0,0.55)) +
  scale_linetype_manual(values = c("natural" = 1, "hatchery" = 2), name = "Rearing Type") +
  # common x-axis scale across all populations
  coord_cartesian(xlim = c(4, 22.5), expand = FALSE) +
  # scale_x_continuous(lim = c(floor(min(covariate_experiences$temp_actual)),ceiling(max(covariate_experiences$temp_actual))), expand = c(0,0)) +
  # scale_color_manual(values = movement_colors) +
  scale_color_tableau(palette = "Tableau 10", name = "Natal Origin") +
  # scale_fill_manual(values =  movement_colors) +
  xlab(expression(~"Temperature" ~ ("°C"))) +
  ylab("Movement probability") +
  theme(panel.grid.major = element_blank(), #element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = "black"),
        legend.key = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.4),
        legend.position = c(0.15, 0.75),
        legend.key.height = unit(0.55, "cm"),
        legend.key.width = unit(0.55, "cm"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 10),
        # turn off the axis titles on each individual plot and just show one for whole plot
        # axis.title = element_blank(),
        # these plot margins are to leave space for the population name on the big figure
        plot.margin = unit(c(0, 0.2, 0.2, 0.2),"cm"))


ggsave(here::here("figures", "paper_figures", "Fig5_deschutes_temp_movements.png"), DES_temp_move_prob_plot, height = 8, width = 10)

