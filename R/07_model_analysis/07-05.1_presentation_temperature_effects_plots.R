# 07-05.1 Temperature effects figure, for presentation

# Start by loading packages and model files
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


### extract all params for Middle Columbia hatchery and wild ###
MCW_DE_param_matrix <- extract_DE_parameters(fit = MCW_fit, fit_summary = MCW_fit_summary)
MCH_DE_param_matrix <- extract_DE_parameters(fit = MCH_fit, fit_summary = MCH_fit_summary)

#### Plot temperature effect on multiple movements on same plot ####

# For some origins, there are multiple interesting and deleterious movements - 
# for example, Middle Columbia fish could go to loss, Deschutes, or overshoot

# plotting function
plot_compare_rear_temp_effect_multiple_movements_pres <- function(origin_select,
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
            axis.text = element_text(size = 15),
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
    # ggsave(here::here("stan_actual", "output", "paper_figures", "01_legend_test.png"), temp_plot_legend_gg, height = 4, width = 4)
    
  } else {
    # suppress common legend - for combined plot
    rear_temp_move_prob_plot <- ggplot(rear_temp_move_prob_quantiles, aes(x = temp_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                          color = Movement, fill = Movement)) +
      geom_line(linewidth = 3) +
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
            axis.text = element_text(size = 18),
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
    # ggsave(here::here("stan_actual", "output", "paper_figures", "01_test.png"), combined_plot, height = 6, width = 6)
    
    
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

#### Create figures for individual origins: ####
# John Day River (W), Umatilla River (H, W), Yakima River (W), Walla Walla River (H, W),
# Wenatchee River (H, W), Entiat River (W), Tucannon River (H, W), Imnaha River (H, W)
### John Day River ###
JDR_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements_pres(origin_select = "John Day River",
                                                                                   wild_move_prob_array = JDR_wild_temp_move_prob_array,
                                                                                   
                                                                                   wild_covariate_experiences = JDR_wild_covariate_experiences,
                                                                                   
                                                                                   movements_evaluated = JDR_movements)

ggsave(here::here("figures", "presentation", "JDR_wild_compare_movement_temp.png"), JDR_wild_compare_movement_temp, height = 8, width = 8)

### Umatilla River ###

# UMA wild plot
UMA_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements_pres(origin_select = "Umatilla River",
                                                                                   wild_move_prob_array = UMA_wild_temp_move_prob_array,
                                                                                   wild_covariate_experiences = UMA_wild_covariate_experiences,
                                                                                   movements_evaluated = UMA_movements)

ggsave(here::here("figures", "presentation", "UMA_wild_compare_movement_temp.png"), UMA_wild_compare_movement_temp, height = 8, width = 8)

# UMA hatchery plot
UMA_hatchery_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements_pres(origin_select = "Umatilla River",
                                                                                       hatchery_move_prob_array = UMA_hatchery_temp_move_prob_array,
                                                                                       hatchery_covariate_experiences = UMA_hatchery_covariate_experiences,
                                                                                       movements_evaluated = UMA_movements)

ggsave(here::here("figures", "presentation", "UMA_hatchery_compare_movement_temp.png"), UMA_hatchery_compare_movement_temp, height = 8, width = 8)

### Yakima River ###

YAK_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements_pres(origin_select = "Yakima River",
                                                                                   wild_move_prob_array = YAK_wild_temp_move_prob_array,
                                                                                   
                                                                                   wild_covariate_experiences = YAK_wild_covariate_experiences,
                                                                                   
                                                                                   movements_evaluated = YAK_movements)

ggsave(here::here("figures", "presentation", "YAK_wild_compare_movement_temp.png"), YAK_wild_compare_movement_temp, height = 8, width = 8)



### Walla Walla River ###

# WAWA wild plot
WAWA_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements_pres(origin_select = "Walla Walla River",
                                                                                    wild_move_prob_array = WAWA_wild_temp_move_prob_array,
                                                                                    wild_covariate_experiences = WAWA_wild_covariate_experiences,
                                                                                    movements_evaluated = WAWA_movements)

ggsave(here::here("figures", "presentation", "WAWA_wild_compare_movement_temp.png"), WAWA_wild_compare_movement_temp, height = 8, width = 8)

# WAWA hatchery plot
WAWA_hatchery_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements_pres(origin_select = "Walla Walla River",
                                                                                        hatchery_move_prob_array = WAWA_hatchery_temp_move_prob_array,
                                                                                        hatchery_covariate_experiences = WAWA_hatchery_covariate_experiences,
                                                                                        movements_evaluated = WAWA_movements)

ggsave(here::here("figures", "presentation", "WAWA_hatchery_compare_movement_temp.png"), WAWA_hatchery_compare_movement_temp, height = 8, width = 8)

#### Combine these all into one plot ####

# Create the legend figure by itself
temp_plot_for_legend <- plot_compare_rear_temp_effect_multiple_movements_pres(origin_select = "John Day River",
                                                                         wild_move_prob_array = JDR_wild_temp_move_prob_array,
                                                                         
                                                                         wild_covariate_experiences = JDR_wild_covariate_experiences,
                                                                         
                                                                         movements_evaluated = JDR_movements,
                                                                         plot_legend= TRUE)
temp_legend <- ggpubr::get_legend(temp_plot_for_legend)
temp_plot_legend_gg <- as_ggplot(temp_legend) + theme(panel.background = element_rect(fill = "white", color = "white"))

combined_movement_temp_plot_pres <- ggarrange(JDR_wild_compare_movement_temp, UMA_wild_compare_movement_temp, UMA_hatchery_compare_movement_temp,
                                         YAK_wild_compare_movement_temp, WAWA_wild_compare_movement_temp, WAWA_hatchery_compare_movement_temp,
                                         nrow = 2, ncol = 3,
                                         labels = c("(A) JDR, Natural", "(B) UMA, Natural", 
                                                    "(C) UMA, Hatchery", "(D) YAK, Natural", 
                                                    "(E) WAWA, Natural", "(F) WAWA, Hatchery"),
                                         label.x = 0.05, label.y = 0.915, font.label = list(size = 24, face = "plain"),
                                         hjust = 0, vjust = 0)

# combined_movement_temp_plot <- annotate_figure(combined_movement_temp_plot,
#                                                bottom = textGrob(expression(~"Temperature" ~ ("°C")), gp = gpar(cex = 1.3)),
#                                                  left = textGrob("Movement probability", rot = 90, gp = gpar(cex = 1.3))) + bgcolor("white")

# Let's try this again, this time using a cowplot solution since ggpubr is
# struggling with the background color
# add a buffer for the y-axis
combined_movement_temp_plot_pres <- cowplot::ggdraw(annotate_figure(combined_movement_temp_plot_pres,
                                                                    left = textGrob(" ", rot = 90, gp = gpar(cex = 0.4)))) +
  theme(plot.background = element_rect(fill="white", color = NA))

combined_movement_temp_plot_pres <- cowplot::ggdraw(annotate_figure(combined_movement_temp_plot_pres,
                                                               bottom = textGrob(expression(~"Temperature" ~ ("°C")), gp = gpar(cex = 2.6)),
                                                               left = textGrob("Movement probability", rot = 90, gp = gpar(cex = 2.6)))) +
  theme(plot.background = element_rect(fill="white", color = NA))


ggsave(here::here("figures", "presentation", "combined_movement_temp_plot_pres.png"), combined_movement_temp_plot_pres, height = 12, width = 16)

