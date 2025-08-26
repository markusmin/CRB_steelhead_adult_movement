// Middle Columbia Wild Steelhead model

// 2025-02-14: Modifications to add run years, drop Hood River as its own state

functions{
    real partial_sum_lpmf( // int[] slice_n_fish, // I don't think that we need this, given that we are using start and end to index the detection histories
                        array[,] int slice_y,
                        int start, int end,
                        // I think we need to re-declare things for our function that are also already declared in our data?
                        int n_ind,
                        int n_states,
                        int n_temp_days,
                        int max_visits,
                        int nmovements,
                        array[,] int parameter_indices_matrix,
                        
                        // declare an array with the date of transition for fish
                        // this will have the same dimensions/structure as states_mat
                        // rows = individual fish, columns = transitions
                        // dates will not be formatted as actual dates, but instead as integers
                        // because this is necessary to use this to index covariates
                        array[,] int transition_dates,
                        
                        // transition seasons (to index temperature parameters)
                        array[] int transition_seasons_vector,
                        
                        // declare temperature data
                        array[,] real temperature_data,
                        
                        // categorical covariates (in this case, just origin parameters)
                        // array[n_ind, 7] int cat_X_mat,
                        array[,] int cat_X_mat,
                        array[,] int temp_X_mat,
                        
                        array[,] int year_X_mat,
                        
                        // array[n_ind, max_visits-1] int states_mat,
                        array[,] int states_mat,
                        // array[n_ind,max_visits] int y, // I don't think that we need this for the same reason as above
                        // array[,] int y,
                        array[] int n_obs,
                        
                        // spill data
                        array[,] real spill_window_data,
                        array[,] real winter_spill_days_data,
                        
                        // spill indexing
                        array[] int winter_post_overshoot_vector,
                        
                        // For detection efficiency calculation purposes: We need two of each of these one for transitions
                        // where we can calculate detection efficiency, one for where we can't. Most parameters
                        // will be shared between the two, but transitions into tributaries won't.
                        // DE = detection efficiency, NDE = no detection efficiency
                        // array[,] real b0_matrix_DE,
                        array[] real b0_vector_DE,
                        // array[,] real borigin1_matrix_DE,
                        array[] real borigin1_vector_DE,
                        array[] real borigin2_vector_DE,
                        array[] real borigin3_vector_DE,
                        array[] real borigin4_vector_DE,
                        array[] real borigin5_vector_DE,
                        array[] real borigin6_vector_DE,
                        
                        // array[,] real b0_matrix_NDE,
                        array[] real b0_vector_NDE,
                        array[] real borigin1_vector_NDE,
                        array[] real borigin2_vector_NDE,
                        array[] real borigin3_vector_NDE,
                        array[] real borigin4_vector_NDE,
                        array[] real borigin5_vector_NDE,
                        array[] real borigin6_vector_NDE,
                        
                        array[] real btemp0_vector_DE, 
                        array[] real btemp0xorigin1_vector_DE, 
                        array[] real btemp0xorigin2_vector_DE, 
                        array[] real btemp0xorigin3_vector_DE, 
                        array[] real btemp0xorigin4_vector_DE, 
                        array[] real btemp0xorigin5_vector_DE, 
                        array[] real btemp0xorigin6_vector_DE, 
                        
                        array[] real btemp0_vector_NDE,
                        array[] real btemp0xorigin1_vector_NDE, 
                        array[] real btemp0xorigin2_vector_NDE, 
                        array[] real btemp0xorigin3_vector_NDE, 
                        array[] real btemp0xorigin4_vector_NDE, 
                        array[] real btemp0xorigin5_vector_NDE, 
                        array[] real btemp0xorigin6_vector_NDE, 
                        
                        array[] real btemp1_vector_DE, 
                        array[] real btemp1xorigin1_vector_DE, 
                        array[] real btemp1xorigin2_vector_DE, 
                        array[] real btemp1xorigin3_vector_DE, 
                        array[] real btemp1xorigin4_vector_DE, 
                        array[] real btemp1xorigin5_vector_DE, 
                        array[] real btemp1xorigin6_vector_DE, 
                        
                        array[] real btemp1_vector_NDE, 
                        array[] real btemp1xorigin1_vector_NDE, 
                        array[] real btemp1xorigin2_vector_NDE, 
                        array[] real btemp1xorigin3_vector_NDE, 
                        array[] real btemp1xorigin4_vector_NDE, 
                        array[] real btemp1xorigin5_vector_NDE, 
                        array[] real btemp1xorigin6_vector_NDE, 
                        
                        // spill parameters
                        array[] real bspillwindow_vector, 
                        array[] real bwinterspill_vector,
                        
                        // matrices to store sigma_year_matrices (DE and NDE)
                        array [] real sigma_yearxorigin1_vector_DE, 
                        array [] real sigma_yearxorigin1_vector_NDE,
                        array [] real sigma_yearxorigin2_vector_DE, 
                        array [] real sigma_yearxorigin2_vector_NDE,
                        array [] real sigma_yearxorigin3_vector_DE, 
                        array [] real sigma_yearxorigin3_vector_NDE,
                        array [] real sigma_yearxorigin4_vector_DE, 
                        array [] real sigma_yearxorigin4_vector_NDE,
                        array [] real sigma_yearxorigin5_vector_DE, 
                        array [] real sigma_yearxorigin5_vector_NDE,
                        array [] real sigma_yearxorigin6_vector_DE, 
                        array [] real sigma_yearxorigin6_vector_NDE,
                        
                        
                        // arrays to store byear parameters
                        // first, create a temporary df, with each of the individual byear parameters (drawn from normal, using sigma_year)
                        // create an empty array with 3 dimensions: matrix for parameters, slices for years
                        // matt trick, transform these outside
                        // array [,,] real byearxorigin1_raw_parameters_array_DE, 
                        // array [,,] real byearxorigin1_raw_parameters_array_NDE,
                        // array [,,] real byearxorigin2_raw_parameters_array_DE, 
                        // array [,,] real byearxorigin2_raw_parameters_array_NDE,
                        // array [,,] real byearxorigin3_raw_parameters_array_DE, 
                        // array [,,] real byearxorigin3_raw_parameters_array_NDE,
                        
                        // these are also arguments, since they're transformed parameters
                        // array [,,] real byearxorigin1_actual_parameters_array_DE, 
                        // array [,,] real byearxorigin1_actual_parameters_array_NDE,
                        // array [,,] real byearxorigin2_actual_parameters_array_DE,
                        // array [,,] real byearxorigin2_actual_parameters_array_NDE,
                        // array [,,] real byearxorigin3_actual_parameters_array_DE,
                        // array [,,] real byearxorigin3_actual_parameters_array_NDE,
                        
                        // reparameterization: Turn the array into a matrix. We can have one dimension be years,
                        // and the other dimension be the parameters, with the same indexing as they are in vectors;
                        // now the parameter indexing is the row of the matrix, rather than the spot in the vector
                        array [,] real byearxorigin1_raw_parameters_matrix_DE, 
                        array [,] real byearxorigin1_raw_parameters_matrix_NDE,
                        array [,] real byearxorigin2_raw_parameters_matrix_DE, 
                        array [,] real byearxorigin2_raw_parameters_matrix_NDE,
                        array [,] real byearxorigin3_raw_parameters_matrix_DE, 
                        array [,] real byearxorigin3_raw_parameters_matrix_NDE,
                        array [,] real byearxorigin4_raw_parameters_matrix_DE, 
                        array [,] real byearxorigin4_raw_parameters_matrix_NDE,
                        array [,] real byearxorigin5_raw_parameters_matrix_DE, 
                        array [,] real byearxorigin5_raw_parameters_matrix_NDE,
                        array [,] real byearxorigin6_raw_parameters_matrix_DE, 
                        array [,] real byearxorigin6_raw_parameters_matrix_NDE,
                        
                        // these are also arguments, since they're transformed parameters
                        array [,] real byearxorigin1_actual_parameters_matrix_DE, 
                        array [,] real byearxorigin1_actual_parameters_matrix_NDE,
                        array [,] real byearxorigin2_actual_parameters_matrix_DE,
                        array [,] real byearxorigin2_actual_parameters_matrix_NDE,
                        array [,] real byearxorigin3_actual_parameters_matrix_DE,
                        array [,] real byearxorigin3_actual_parameters_matrix_NDE,
                        array [,] real byearxorigin4_actual_parameters_matrix_DE,
                        array [,] real byearxorigin4_actual_parameters_matrix_NDE,
                        array [,] real byearxorigin5_actual_parameters_matrix_DE,
                        array [,] real byearxorigin5_actual_parameters_matrix_NDE,
                        array [,] real byearxorigin6_actual_parameters_matrix_DE,
                        array [,] real byearxorigin6_actual_parameters_matrix_NDE,
                        
                        
                        // below is new data for detection efficiency
                        // declare an array that you can put detection probabilities in, 
                        // indexed by rows = run years, columns = tributaries.
                        // We will calculate estimated detection probability
                        // This matrix will contain the parameters (both alpha/intercept from era
                        // and beta, for slope of discharge relationship)
                        // array[,] real det_eff_param_matrix;
                        // array[,] real discharge_matrix;
                        
                        // different approach - include it directly in the model
                        // need to declare an array with dimensions (number of tributaries for which we can estimate detection efficiency) x 
                        // (number of run years) x (number of detection efficiency parameters)
                        // this is in effect an array of design matrices
                        array[,,] real tributary_design_matrices_array,
                        
                        // declare a vector to store the run year that each transition occurs in
                        array[] int transition_run_years, // or does this need to be array[] int?
                        
                        // declare how many run years there are
                        int nyears,
                        
                        // declare the array that says in which run years in which states we calculate DE
                        array[,,] int run_year_DE_array,
                        
                        // declare the vector that contains the parameters for detection efficiency
                        vector det_eff_param_vector
                        ) { // I don't think we need this either? Since we're just indexing it again with start and end
                          
  // First, declare the total lp (log probability) variable that we will be returning
  real total_lp = 0;
  
  
  // Now, loop through the detection matrices for each individual IN THIS SLICE OF THE DATA
  // for (i in 1:slice_n_fish){
    
    // start - end apparently doesn't work because reduce_sum() resets the indices for each slice (confusing) - according to this post: https://discourse.mc-stan.org/t/parallelization-in-cmdstanr-using-reduce-sum-on-a-multinomial-likelihood/24607/7
  // or maybe not, this post seems to contradict this: https://discourse.mc-stan.org/t/help-with-multi-threading-a-simple-ordinal-probit-model-using-reduce-sum/15353
  // i is the number of fish
  for (i in start:end){
    // for (i in 1:end-start+1){
    // Loop through each of the observations, stopping at the loss column
    
    // print("fish # ", i);
    
    // Here, create a vector to store the lp at each observation for a fish
    // vector[n_obs[i]] lp_fish;
    // Let's initialize this instead as a real value starting at zero
    real lp_fish = 0;
    // j is the index of the observation (i.e., each individual state transition)
    for (j in 1:n_obs[i]){
      // for (j in 1:n_obs[i - start + 1]){
          // print("i = ",i);
            // print("j = ",j);

        // vector for logits
        vector[41] logits;
        
        // Get the index of the current state
        int current;
        current = states_mat[i,j];
        // current = states_mat[i - start + 1,j];
        
        // Get the current date
        int date;
        date = transition_dates[i,j];
        
        // Get the temperature in the current state at the current time, based on the state x date
        // temperature only for eight mainstem sites
        real temp;
        if (current >= 2 && current <= 9){
            temp = temperature_data[date, current];
            
        } else {
          temp = 0; 
        }
        
        // Get the spill (window) in the current state at the current time, based on the state x date
        // spill only for eight mainstem sites
        real spill_window;
        if (current >= 2 && current <= 9){
            spill_window = spill_window_data[date, current];
            
        } else {
          spill_window = 0; 
        }
        
        // Get the spill (days in each month)
        // use the vector that we create for whether or not fish could possibly 
        // have experienced these conditions
        array[nyears] real winter_spill_days;
        if (current >= 2 && current <= 9){
            winter_spill_days = winter_spill_days_data[, current];
            
        } else {
          winter_spill_days = rep_array(0, nyears); 
        }
        
        
        // here: create a vector with the length of the states, where we will store whether or not they are affected by DE
        vector[40] DE_correction;
        
            // Get the current run year
            
            // indexing is different if i == 1
              int run_year_index;      
              int run_year_actual;
              // index for transition season as well
              int season_actual;
              
              // index for monthly spill effects
              int winter_spill_index;
              
              
            if (i == 1) {
              // now, sum the vector to get a single integer value for indexing the run year vector
              run_year_index = j;
              // get the actual run year
              run_year_actual = transition_run_years[run_year_index];
              // index for transition season as well
              season_actual = transition_seasons_vector[run_year_index];
              
              // index for monthly spill effects as well
              winter_spill_index = winter_post_overshoot_vector[run_year_index];
              
            } else {
                // indexing tweaks
                // first: declare a vector that stores all of the indices that we need to sum to get the right run year
                array[i] int run_year_indices_vector;
                // then populate: first elements are all the number of transitions from all previous fish
                run_year_indices_vector[1:(i-1)] = n_obs[1:(i-1)];
                // last element is the transition number of the current fish  
                run_year_indices_vector[i] = j;
                // now, sum the vector to get a single integer value for indexing the run year vector
                run_year_index = sum(run_year_indices_vector);
                // get the actual run year
                run_year_actual = transition_run_years[run_year_index];
                
                // index for transition season as well
              season_actual = transition_seasons_vector[run_year_index];
              
              // index for monthly spill effects as well
              winter_spill_index = winter_post_overshoot_vector[run_year_index];

            }
          

        // Populate each of the first 40 (non-loss)
        for (k in 1:40){
          
          // here - add an if else statement. For each transition (combination of current state, and to state),
          // use the run year of the fish to decide whether to pull from the DE or NDE matrix.
          // we need an array that we can index - current (from, rows) x k (to, columns) x run years (slices) - and the values
          // indicate which of the two matrices to use
          
          // what we'll do:
          // create a vector (or array?). Every time we need to select from the DE matrix, note
          // which transition this is. Then, once we're done, we can then index to those and do
          // the correction, where we multiply by a p, and then subtract that from the loss category
          
          // If we are in a run year in a state transition where we are estimating DE, use that matrix;
          // for indexing by run year - some up all n_obs prior to this fish, plus whatever index we're on for the current fish
          // need to make an exception for the first fish, which uses only the index j

            // If we're in a state/run year combo that needs DE correction, correct for it
                    // if (run_year_DE_array[current,k,transition_run_years[sum(n_obs[1:i-1], j)]] == 1){
                      // need to declare another intermediate value: an integer that takes either 0 or 1, for whether it's DE or NDE
                      int DE_index;
                      // now, determine if it's 0 or 1
                      DE_index = run_year_DE_array[current,k,run_year_actual];
                      // if (run_year_DE_array[current,k,transition_run_years[run_year_indices_vector]] == 1){
                  if (DE_index == 1){
                    int param_index;
                    param_index = parameter_indices_matrix[current, k];
                    logits[k] = b0_vector_DE[param_index] + 
                    
                    // temp effects
                    temp_X_mat[i,1] * ((1 - season_actual) * btemp0_vector_DE[param_index] + season_actual * btemp1_vector_DE[param_index]) * temp +
                    temp_X_mat[i,2] * ((1 - season_actual) * btemp0xorigin1_vector_DE[param_index] + season_actual * btemp1xorigin1_vector_DE[param_index]) * temp +
                    temp_X_mat[i,3] * ((1 - season_actual) * btemp0xorigin2_vector_DE[param_index] + season_actual * btemp1xorigin2_vector_DE[param_index]) * temp +
                    temp_X_mat[i,4] * ((1 - season_actual) * btemp0xorigin3_vector_DE[param_index] + season_actual * btemp1xorigin3_vector_DE[param_index]) * temp +
                    temp_X_mat[i,5] * ((1 - season_actual) * btemp0xorigin4_vector_DE[param_index] + season_actual * btemp1xorigin4_vector_DE[param_index]) * temp +
                    temp_X_mat[i,6] * ((1 - season_actual) * btemp0xorigin5_vector_DE[param_index] + season_actual * btemp1xorigin5_vector_DE[param_index]) * temp +
                    temp_X_mat[i,7] * ((1 - season_actual) * btemp0xorigin6_vector_DE[param_index] + season_actual * btemp1xorigin6_vector_DE[param_index]) * temp +
                    
                    // spill effects
                    // this next line is necessary to make sure that post-overshoot fish
                    // don't also get the spill_window covariate
                    // 2024-02-05
                    max(1 - (winter_spill_index), 0) *
                    bspillwindow_vector[param_index] * spill_window +
                    winter_spill_index * bwinterspill_vector[param_index] * winter_spill_days[run_year_actual] +
                    
                    // origin effects
                    cat_X_mat[i,1] * borigin1_vector_DE[param_index] +
                    cat_X_mat[i,2] * borigin2_vector_DE[param_index] +
                    cat_X_mat[i,3] * borigin3_vector_DE[param_index] +
                    cat_X_mat[i,4] * borigin4_vector_DE[param_index] +
                    cat_X_mat[i,5] * borigin5_vector_DE[param_index] +
                    cat_X_mat[i,6] * borigin6_vector_DE[param_index] +
                    
                    // add year as a random effect
                    year_X_mat[i,1] * byearxorigin1_actual_parameters_matrix_DE[param_index, run_year_actual] +
                    year_X_mat[i,2] * byearxorigin2_actual_parameters_matrix_DE[param_index, run_year_actual] +
                    year_X_mat[i,3] * byearxorigin3_actual_parameters_matrix_DE[param_index, run_year_actual] +
                    year_X_mat[i,4] * byearxorigin4_actual_parameters_matrix_DE[param_index, run_year_actual] +
                    year_X_mat[i,5] * byearxorigin5_actual_parameters_matrix_DE[param_index, run_year_actual] +
                    year_X_mat[i,6] * byearxorigin6_actual_parameters_matrix_DE[param_index, run_year_actual];
                    
                    // store in the vector that it's DE
                    DE_correction[k] = 1;
                    
                    // Else, use the NDE matrix, and don't correct for detection efficiency
                  } else {
                    int param_index;
                    param_index = parameter_indices_matrix[current, k];
                    logits[k] = b0_vector_NDE[param_index] + 
                    
                    // temp effects
                    temp_X_mat[i,1] * ((1 - season_actual) * btemp0_vector_NDE[param_index] + season_actual * btemp1_vector_NDE[param_index]) * temp +
                    temp_X_mat[i,2] * ((1 - season_actual) * btemp0xorigin1_vector_NDE[param_index] + season_actual * btemp1xorigin1_vector_NDE[param_index]) * temp +
                    temp_X_mat[i,3] * ((1 - season_actual) * btemp0xorigin2_vector_NDE[param_index] + season_actual * btemp1xorigin2_vector_NDE[param_index]) * temp +
                    temp_X_mat[i,4] * ((1 - season_actual) * btemp0xorigin3_vector_NDE[param_index] + season_actual * btemp1xorigin3_vector_NDE[param_index]) * temp +
                    temp_X_mat[i,5] * ((1 - season_actual) * btemp0xorigin4_vector_NDE[param_index] + season_actual * btemp1xorigin4_vector_NDE[param_index]) * temp +
                    temp_X_mat[i,6] * ((1 - season_actual) * btemp0xorigin5_vector_NDE[param_index] + season_actual * btemp1xorigin5_vector_NDE[param_index]) * temp +
                    temp_X_mat[i,7] * ((1 - season_actual) * btemp0xorigin6_vector_NDE[param_index] + season_actual * btemp1xorigin6_vector_NDE[param_index]) * temp +
                    
                    // spill effects
                    // this next line is necessary to make sure that post-overshoot fish
                    // don't also get the spill_window covariate
                    // 2024-02-05
                    max(1 - (winter_spill_index), 0) *
                    bspillwindow_vector[param_index] * spill_window +
                    winter_spill_index * bwinterspill_vector[param_index] * winter_spill_days[run_year_actual] +
                    
                    // origin effects
                    cat_X_mat[i,1] * borigin1_vector_NDE[param_index] +
                    cat_X_mat[i,2] * borigin2_vector_NDE[param_index] +
                    cat_X_mat[i,3] * borigin3_vector_NDE[param_index] +
                    cat_X_mat[i,4] * borigin4_vector_NDE[param_index] +
                    cat_X_mat[i,5] * borigin5_vector_NDE[param_index] +
                    cat_X_mat[i,6] * borigin6_vector_NDE[param_index] +
                    
                    // add year as a random effect
                    year_X_mat[i,1] * byearxorigin1_actual_parameters_matrix_NDE[param_index, run_year_actual] +
                    year_X_mat[i,2] * byearxorigin2_actual_parameters_matrix_NDE[param_index, run_year_actual] +
                    year_X_mat[i,3] * byearxorigin3_actual_parameters_matrix_NDE[param_index, run_year_actual] +
                    year_X_mat[i,4] * byearxorigin4_actual_parameters_matrix_NDE[param_index, run_year_actual] +
                    year_X_mat[i,5] * byearxorigin5_actual_parameters_matrix_NDE[param_index, run_year_actual] +
                    year_X_mat[i,6] * byearxorigin6_actual_parameters_matrix_NDE[param_index, run_year_actual];
                    
                    
                    
                    // otherwise, store in the vector that it's not DE
                    DE_correction[k] = 0;
                    }
                    
                    
          
          
        }
          // loss param
          logits[41] = 0;
                    
          // declare vector to store true movement probabilities
          vector[41] p_vec_actual;
          
          // declare vector to store observed movement probabilities, which are a product of the true probabilities, mulitplied by detection efficiency
          vector[41] p_vec_observed;
          
          // print("i!=1; logits = ",logits);
          
          
          // proportion vector, uncorrected for detection efficiency
          p_vec_actual = softmax(logits);
          
          // now, loop through transitions again, and calculate detection efficiency for each where DE_correction == 1
          // declare one vector for linear predictors and one for actual detection efficiency
          
          vector[40] det_eff_eta;
          vector[40] det_eff;
        
        for (k in 1:40){
          if (DE_correction[k] == 1) {

            // to calculate detection efficiency by indexing:
            // The tributary design matrices array will have the same number of slices as states (40, for 41 - loss).
            // Only 14 of these slices (the 14 tributaries that have DE calculations) will have non-zero values, but this will make the indexing simpler.
            // we will index slices using k.
            // Rows will be run years, indexed using transition_run_years[sum(n_obs[1:i-1], j)], same as above
            // the columns will then be the appropriate row of the design matrix, for that state and tributary.
            // That will then be multiplied by the full, 34 length parameter vector. This will be written out in
            // the transformed parameters section, and will have 20 alpha terms and 14 beta terms.

            // vector[34] trib_design_matrix_result;
            // trib_design_matrix_result = to_row_vector(tributary_design_matrices_array[run_year_actual,,k]) * det_eff_param_vector;
            // det_eff_eta[k] = sum(trib_design_matrix_result);
            det_eff_eta[k] = to_row_vector(tributary_design_matrices_array[run_year_actual,,k]) * to_vector(det_eff_param_vector);
            det_eff[k] = exp(det_eff_eta[k])/(1 + exp(det_eff_eta[k]));


          } else {
            // If we don't have to calculate a detection efficiency, don't do anything

          }


        }
          
          
          // now loop through transitions and modify p_vec to account for this
        
          
          // // Create a vector to store all of the loss corrections for detection efficiency
          vector[40] loss_term_DE_corrections;
          
       
        for (k in 1:40){
          if (DE_correction[k] == 1) {
            p_vec_observed[k] = p_vec_actual[k] * det_eff[k];
            // p_vec_observed[k] = p_vec_actual[k]; // again this is just for testing - remove det eff correction for now

            // each time you modify a term, modify the loss term by the same amount
            loss_term_DE_corrections[k] = p_vec_actual[k] * (1 - det_eff[k]);

          } else {
            p_vec_observed[k] = p_vec_actual[k];

            // Just put a zero there, otherwise it'll be NAs
            loss_term_DE_corrections[k] = 0;

          }


        }
        
        // Once you've looped through all states, correct loss term
        p_vec_observed[41] = p_vec_actual[41] + sum(loss_term_DE_corrections);
        
                   // print("i = ",i);
           // print("j = ",j);
           // print("current state: ", current);
           // print("logits = ", logits);
           // print("p_vec_actual: ", p_vec_actual);
           // print("p_vec_observed: ", p_vec_observed);
           // print("current state: ", current);
        // print("next state: ", slice_y[i - start + 1,j+1], "; prob of next state = ", p_vec_observed[slice_y[i - start + 1,j+1]]);

        
        lp_fish += categorical_lpmf(slice_y[i - start + 1,j+1] | p_vec_observed);
      
    }
    
    // Now, increment the log density
    total_lp += lp_fish;
    
    
} // end looping through fish in this slice

return total_lp;

} // end partial_sum_lpmf function

} // end functions block

data {
  int max_visits;
  int n_ind; // number of individuals (this is nfish)
  array[n_ind,max_visits] int y; // array with dimensions nfish (number of fish), 41 (maximum number of site visits)
  array[n_ind] int n_obs; // n_obs is a vector that contains the number of site visits for each individual - but needs to be declared as array to take integer values
	int n_states; // the number of states in our model
  vector[n_states] possible_movements; // a vector containing the number of possible transitions out of each state
  array[n_ind, max_visits] int transition_dates; // a matrix of the dates (as numeric string with 2005-06-01 as 1) that each transition took place
  int n_temp_days; // the number of days of temperature data that we have
  array[n_temp_days,11] real temperature_data; // a matrix of the temperature data for each of eight states
  array[n_ind, max_visits-1] int states_mat; // a matrix (array to take integer values) with rows = each fish and columns = the site visits
  // array[54, 2] int movements; // a matrix that contains the possible transitions out of each state
  int nmovements; // an integer value containing the number of possible transitions (same as rows in movements data)
  array[n_states,n_states] int parameter_indices_matrix; // this is the matrix that maps the movements to their spots in a vector containing the parameters
  array[nmovements, 2] int movements; // a matrix that contains the possible transitions out of each state
  // array[758,2] int not_movements; // a matrix containing all non-allowed state transitions
  int n_notmovements; // an integer value containing the number of non- possible transitions (same as rows in not_movements data)
  array[n_notmovements,2] int not_movements; // a matrix containing all non-allowed state transitions
  // array[n_ind, 48] int dates; // a matrix containing dates (as integers) where rows = number of fish and columns = site visits
  matrix[n_states, n_states] possible_states; // the transition matrix with 1s for allowable transitions and 0s for non-allowable
  array[n_ind, 6] int cat_X_mat; // a matrix with rows = individual fish and columns = various categorical covariates (origin)
  array[n_ind, 7] int temp_X_mat; // a matrix with rows = individual fish and columns = design matrix for temperature
  array[n_ind, 6] int year_X_mat; // a matrix with rows = individual fish and columns = origins for year effects
  
  int nyears;
  
  // the spill data - 5 dfs
  array[n_temp_days,11] real spill_window_data; // the window-based spill estimates for en-route fallback
  array[nyears,9] real winter_spill_days_data; // the days of spill in the months of january, february, and march
  
  
  matrix[33,2] det_eff_param_posteriors; // declare the matrix that contains the posteriors from the det eff script:
  
  // array that contains design matrices for each tributary
  array[20,33,40] int tributary_design_matrices_array;
  
  // vector that contains the run years in which each individual transition occurred
  int ntransitions;
  array[ntransitions] int transition_run_years;
  
  // a vector of the seasons that each transition took place (0 = winter/spring, before June 1; 1 = summer/fall, June 1 of later)
  array[ntransitions] int transition_seasons_vector; 
  
  // vectors that contain whether or not fish were in states to experience post-overshoot spill conditions
  array[ntransitions] int winter_post_overshoot_vector;
  
  
  // array that contains which run years in which transitions need DE correction
  array[n_states,n_states,20] int run_year_DE_array;
  
  // Declare data for parallelization (reduce sum function)
  int<lower=0> N;
  // array[N,max_visits] int slice_y; # doesn't need to be declared as data!
  // array[N,max_visits] int n_fish;
  int<lower=1> grainsize;
}

transformed data {
real sigma_yearxorigin1_vector_3_2 = 0;
real sigma_yearxorigin1_vector_4_3 = 0;
real sigma_yearxorigin3_vector_3_2 = 0;
real sigma_yearxorigin3_vector_4_3 = 0;
real sigma_yearxorigin2_vector_3_2 = 0;
real sigma_yearxorigin2_vector_4_3 = 0;
real sigma_yearxorigin4_vector_3_2 = 0;
real sigma_yearxorigin4_vector_4_3 = 0;
real sigma_yearxorigin6_vector_4_3 = 0;
real sigma_yearxorigin5_vector_4_3 = 0;
}

// The parameters accepted by the model. Our model has only one matrix of parameters
parameters {
  // matrix[10,10] b0_matrix; // parameters we are monitoring - in this case it's the intercept matrix
  // Let's instead only select certain parameters, then put them in a matrix later
real b0_matrix_4_5;
real b0_matrix_5_4;
real b0_matrix_5_6;
real b0_matrix_5_22_DE;
real b0_matrix_5_22_NDE;
real b0_matrix_6_5;
real b0_matrix_6_7;
real b0_matrix_7_6;
real b0_matrix_7_26_DE;
real b0_matrix_7_26_NDE;
real b0_matrix_7_28_DE;
real b0_matrix_7_28_NDE;
real b0_matrix_8_9;
real b0_matrix_8_30_DE;
real b0_matrix_8_30_NDE;
real b0_matrix_9_8;
real b0_matrix_9_32_DE;
real b0_matrix_9_32_NDE;
real b0_matrix_9_34;
real b0_matrix_9_35;
real b0_matrix_9_36;
real b0_matrix_9_37_DE;
real b0_matrix_9_37_NDE;
real b0_matrix_22_5;
real b0_matrix_26_7;
real b0_matrix_28_7;
real b0_matrix_30_8;
real b0_matrix_32_9;
real b0_matrix_34_9;
real b0_matrix_35_9;
real b0_matrix_36_9;
real b0_matrix_37_9;

// spill parameters
real bspillwindow_matrix_2_1;
real bspillwindow_matrix_3_2;
real bspillwindow_matrix_4_3;
real bspillwindow_matrix_5_4;
real bspillwindow_matrix_6_5;
real bspillwindow_matrix_7_6;
real bspillwindow_matrix_8_3;
real bspillwindow_matrix_9_8;

real bwinterspill_matrix_3_2;
real bwinterspill_matrix_4_3;
real bwinterspill_matrix_5_4;
real bwinterspill_matrix_6_5;
real bwinterspill_matrix_7_6;
real bwinterspill_matrix_8_3;
real bwinterspill_matrix_9_8;

// temperature parameters
real btemp0_matrix_4_5;
real btemp1_matrix_4_5;
real btemp0_matrix_5_6;
real btemp1_matrix_5_6;
real btemp0_matrix_5_22_DE;
real btemp0_matrix_5_22_NDE;
real btemp1_matrix_5_22_DE;
real btemp1_matrix_5_22_NDE;
real btemp0_matrix_6_7;
real btemp1_matrix_6_7;
real btemp0_matrix_7_26_DE;
real btemp0_matrix_7_26_NDE;
real btemp1_matrix_7_26_DE;
real btemp1_matrix_7_26_NDE;
real btemp0_matrix_7_28_DE;
real btemp0_matrix_7_28_NDE;
real btemp1_matrix_7_28_DE;
real btemp1_matrix_7_28_NDE;
real btemp0_matrix_8_9;
real btemp1_matrix_8_9;
real btemp0_matrix_8_30_DE;
real btemp0_matrix_8_30_NDE;
real btemp1_matrix_8_30_DE;
real btemp1_matrix_8_30_NDE;
real btemp0_matrix_9_32_DE;
real btemp0_matrix_9_32_NDE;
real btemp1_matrix_9_32_DE;
real btemp1_matrix_9_32_NDE;
real btemp0_matrix_9_34;
real btemp1_matrix_9_34;
real btemp0_matrix_9_35;
real btemp1_matrix_9_35;
real btemp0_matrix_9_36;
real btemp1_matrix_9_36;
real btemp0_matrix_9_37_DE;
real btemp0_matrix_9_37_NDE;
real btemp1_matrix_9_37_DE;
real btemp1_matrix_9_37_NDE;

real btemp0xorigin1_matrix_1_2;
real btemp1xorigin1_matrix_1_2;
real btemp0xorigin1_matrix_2_3;
real btemp1xorigin1_matrix_2_3;
real btemp0xorigin1_matrix_2_10_DE;
real btemp0xorigin1_matrix_2_10_NDE;
real btemp1xorigin1_matrix_2_10_DE;
real btemp1xorigin1_matrix_2_10_NDE;
real btemp0xorigin1_matrix_2_12_DE;
real btemp0xorigin1_matrix_2_12_NDE;
real btemp1xorigin1_matrix_2_12_DE;
real btemp1xorigin1_matrix_2_12_NDE;
real btemp0xorigin1_matrix_2_14_DE;
real btemp0xorigin1_matrix_2_14_NDE;
real btemp1xorigin1_matrix_2_14_DE;
real btemp1xorigin1_matrix_2_14_NDE;
real btemp0xorigin1_matrix_2_16_DE;
real btemp0xorigin1_matrix_2_16_NDE;
real btemp1xorigin1_matrix_2_16_DE;
real btemp1xorigin1_matrix_2_16_NDE;
real btemp0xorigin1_matrix_2_39;
real btemp1xorigin1_matrix_2_39;
real btemp0xorigin1_matrix_3_4;
real btemp1xorigin1_matrix_3_4;
real btemp0xorigin1_matrix_3_8;
real btemp1xorigin1_matrix_3_8;
real btemp0xorigin1_matrix_3_18_DE;
real btemp0xorigin1_matrix_3_18_NDE;
real btemp1xorigin1_matrix_3_18_DE;
real btemp1xorigin1_matrix_3_18_NDE;
real btemp0xorigin1_matrix_3_20_DE;
real btemp0xorigin1_matrix_3_20_NDE;
real btemp1xorigin1_matrix_3_20_DE;
real btemp1xorigin1_matrix_3_20_NDE;

real btemp0xorigin2_matrix_1_2;
real btemp1xorigin2_matrix_1_2;
real btemp0xorigin2_matrix_2_3;
real btemp1xorigin2_matrix_2_3;
real btemp0xorigin2_matrix_2_10_DE;
real btemp0xorigin2_matrix_2_10_NDE;
real btemp1xorigin2_matrix_2_10_DE;
real btemp1xorigin2_matrix_2_10_NDE;
real btemp0xorigin2_matrix_2_12_DE;
real btemp0xorigin2_matrix_2_12_NDE;
real btemp1xorigin2_matrix_2_12_DE;
real btemp1xorigin2_matrix_2_12_NDE;
real btemp0xorigin2_matrix_2_14_DE;
real btemp0xorigin2_matrix_2_14_NDE;
real btemp1xorigin2_matrix_2_14_DE;
real btemp1xorigin2_matrix_2_14_NDE;
real btemp0xorigin2_matrix_2_16_DE;
real btemp0xorigin2_matrix_2_16_NDE;
real btemp1xorigin2_matrix_2_16_DE;
real btemp1xorigin2_matrix_2_16_NDE;
real btemp0xorigin2_matrix_2_39;
real btemp1xorigin2_matrix_2_39;
real btemp0xorigin2_matrix_3_4;
real btemp1xorigin2_matrix_3_4;
real btemp0xorigin2_matrix_3_8;
real btemp1xorigin2_matrix_3_8;
real btemp0xorigin2_matrix_3_18_DE;
real btemp0xorigin2_matrix_3_18_NDE;
real btemp1xorigin2_matrix_3_18_DE;
real btemp1xorigin2_matrix_3_18_NDE;
real btemp0xorigin2_matrix_3_20_DE;
real btemp0xorigin2_matrix_3_20_NDE;
real btemp1xorigin2_matrix_3_20_DE;
real btemp1xorigin2_matrix_3_20_NDE;

real btemp0xorigin3_matrix_1_2;
real btemp1xorigin3_matrix_1_2;
real btemp0xorigin3_matrix_2_3;
real btemp1xorigin3_matrix_2_3;
real btemp0xorigin3_matrix_2_10_DE;
real btemp0xorigin3_matrix_2_10_NDE;
real btemp1xorigin3_matrix_2_10_DE;
real btemp1xorigin3_matrix_2_10_NDE;
real btemp0xorigin3_matrix_2_12_DE;
real btemp0xorigin3_matrix_2_12_NDE;
real btemp1xorigin3_matrix_2_12_DE;
real btemp1xorigin3_matrix_2_12_NDE;
real btemp0xorigin3_matrix_2_14_DE;
real btemp0xorigin3_matrix_2_14_NDE;
real btemp1xorigin3_matrix_2_14_DE;
real btemp1xorigin3_matrix_2_14_NDE;
real btemp0xorigin3_matrix_2_16_DE;
real btemp0xorigin3_matrix_2_16_NDE;
real btemp1xorigin3_matrix_2_16_DE;
real btemp1xorigin3_matrix_2_16_NDE;
real btemp0xorigin3_matrix_2_39;
real btemp1xorigin3_matrix_2_39;
real btemp0xorigin3_matrix_3_4;
real btemp1xorigin3_matrix_3_4;
real btemp0xorigin3_matrix_3_8;
real btemp1xorigin3_matrix_3_8;
real btemp0xorigin3_matrix_3_18_DE;
real btemp0xorigin3_matrix_3_18_NDE;
real btemp1xorigin3_matrix_3_18_DE;
real btemp1xorigin3_matrix_3_18_NDE;
real btemp0xorigin3_matrix_3_20_DE;
real btemp0xorigin3_matrix_3_20_NDE;
real btemp1xorigin3_matrix_3_20_DE;
real btemp1xorigin3_matrix_3_20_NDE;

real btemp0xorigin4_matrix_1_2;
real btemp1xorigin4_matrix_1_2;
real btemp0xorigin4_matrix_2_3;
real btemp1xorigin4_matrix_2_3;
real btemp0xorigin4_matrix_2_10_DE;
real btemp0xorigin4_matrix_2_10_NDE;
real btemp1xorigin4_matrix_2_10_DE;
real btemp1xorigin4_matrix_2_10_NDE;
real btemp0xorigin4_matrix_2_12_DE;
real btemp0xorigin4_matrix_2_12_NDE;
real btemp1xorigin4_matrix_2_12_DE;
real btemp1xorigin4_matrix_2_12_NDE;
real btemp0xorigin4_matrix_2_14_DE;
real btemp0xorigin4_matrix_2_14_NDE;
real btemp1xorigin4_matrix_2_14_DE;
real btemp1xorigin4_matrix_2_14_NDE;
real btemp0xorigin4_matrix_2_16_DE;
real btemp0xorigin4_matrix_2_16_NDE;
real btemp1xorigin4_matrix_2_16_DE;
real btemp1xorigin4_matrix_2_16_NDE;
real btemp0xorigin4_matrix_2_39;
real btemp1xorigin4_matrix_2_39;
real btemp0xorigin4_matrix_3_4;
real btemp1xorigin4_matrix_3_4;
real btemp0xorigin4_matrix_3_8;
real btemp1xorigin4_matrix_3_8;
real btemp0xorigin4_matrix_3_18_DE;
real btemp0xorigin4_matrix_3_18_NDE;
real btemp1xorigin4_matrix_3_18_DE;
real btemp1xorigin4_matrix_3_18_NDE;
real btemp0xorigin4_matrix_3_20_DE;
real btemp0xorigin4_matrix_3_20_NDE;
real btemp1xorigin4_matrix_3_20_DE;
real btemp1xorigin4_matrix_3_20_NDE;

real btemp0xorigin5_matrix_1_2;
real btemp1xorigin5_matrix_1_2;
real btemp0xorigin5_matrix_2_3;
real btemp1xorigin5_matrix_2_3;
real btemp0xorigin5_matrix_2_10_DE;
real btemp0xorigin5_matrix_2_10_NDE;
real btemp1xorigin5_matrix_2_10_DE;
real btemp1xorigin5_matrix_2_10_NDE;
real btemp0xorigin5_matrix_2_12_DE;
real btemp0xorigin5_matrix_2_12_NDE;
real btemp1xorigin5_matrix_2_12_DE;
real btemp1xorigin5_matrix_2_12_NDE;
real btemp0xorigin5_matrix_2_14_DE;
real btemp0xorigin5_matrix_2_14_NDE;
real btemp1xorigin5_matrix_2_14_DE;
real btemp1xorigin5_matrix_2_14_NDE;
real btemp0xorigin5_matrix_2_16_DE;
real btemp0xorigin5_matrix_2_16_NDE;
real btemp1xorigin5_matrix_2_16_DE;
real btemp1xorigin5_matrix_2_16_NDE;
real btemp0xorigin5_matrix_2_39;
real btemp1xorigin5_matrix_2_39;
real btemp0xorigin5_matrix_3_4;
real btemp1xorigin5_matrix_3_4;
real btemp0xorigin5_matrix_3_8;
real btemp1xorigin5_matrix_3_8;
real btemp0xorigin5_matrix_3_18_DE;
real btemp0xorigin5_matrix_3_18_NDE;
real btemp1xorigin5_matrix_3_18_DE;
real btemp1xorigin5_matrix_3_18_NDE;
real btemp0xorigin5_matrix_3_20_DE;
real btemp0xorigin5_matrix_3_20_NDE;
real btemp1xorigin5_matrix_3_20_DE;
real btemp1xorigin5_matrix_3_20_NDE;

real btemp0xorigin6_matrix_1_2;
real btemp1xorigin6_matrix_1_2;
real btemp0xorigin6_matrix_2_3;
real btemp1xorigin6_matrix_2_3;
real btemp0xorigin6_matrix_2_10_DE;
real btemp0xorigin6_matrix_2_10_NDE;
real btemp1xorigin6_matrix_2_10_DE;
real btemp1xorigin6_matrix_2_10_NDE;
real btemp0xorigin6_matrix_2_12_DE;
real btemp0xorigin6_matrix_2_12_NDE;
real btemp1xorigin6_matrix_2_12_DE;
real btemp1xorigin6_matrix_2_12_NDE;
real btemp0xorigin6_matrix_2_14_DE;
real btemp0xorigin6_matrix_2_14_NDE;
real btemp1xorigin6_matrix_2_14_DE;
real btemp1xorigin6_matrix_2_14_NDE;
real btemp0xorigin6_matrix_2_16_DE;
real btemp0xorigin6_matrix_2_16_NDE;
real btemp1xorigin6_matrix_2_16_DE;
real btemp1xorigin6_matrix_2_16_NDE;
real btemp0xorigin6_matrix_2_39;
real btemp1xorigin6_matrix_2_39;
real btemp0xorigin6_matrix_3_4;
real btemp1xorigin6_matrix_3_4;
real btemp0xorigin6_matrix_3_8;
real btemp1xorigin6_matrix_3_8;
real btemp0xorigin6_matrix_3_18_DE;
real btemp0xorigin6_matrix_3_18_NDE;
real btemp1xorigin6_matrix_3_18_DE;
real btemp1xorigin6_matrix_3_18_NDE;
real btemp0xorigin6_matrix_3_20_DE;
real btemp0xorigin6_matrix_3_20_NDE;
real btemp1xorigin6_matrix_3_20_DE;
real btemp1xorigin6_matrix_3_20_NDE;

// origins parameters

real borigin1_matrix_1_2;
real borigin1_matrix_2_1;
real borigin1_matrix_2_3;
real borigin1_matrix_2_10_DE;
real borigin1_matrix_2_10_NDE;
real borigin1_matrix_2_12_DE;
real borigin1_matrix_2_12_NDE;
real borigin1_matrix_2_14_DE;
real borigin1_matrix_2_14_NDE;
real borigin1_matrix_2_16_DE;
real borigin1_matrix_2_16_NDE;
real borigin1_matrix_2_39;
real borigin1_matrix_3_2;
real borigin1_matrix_3_4;
real borigin1_matrix_3_8;
real borigin1_matrix_3_18_DE;
real borigin1_matrix_3_18_NDE;
real borigin1_matrix_3_20_DE;
real borigin1_matrix_3_20_NDE;
real borigin1_matrix_4_3;
real borigin1_matrix_8_3;
real borigin1_matrix_10_2;
real borigin1_matrix_12_2;
real borigin1_matrix_14_2;
real borigin1_matrix_16_2;
real borigin1_matrix_18_3;
real borigin1_matrix_20_3;
real borigin1_matrix_39_2;

real borigin2_matrix_1_2;
real borigin2_matrix_2_1;
real borigin2_matrix_2_3;
real borigin2_matrix_2_10_DE;
real borigin2_matrix_2_10_NDE;
real borigin2_matrix_2_12_DE;
real borigin2_matrix_2_12_NDE;
real borigin2_matrix_2_14_DE;
real borigin2_matrix_2_14_NDE;
real borigin2_matrix_2_16_DE;
real borigin2_matrix_2_16_NDE;
real borigin2_matrix_2_39;
real borigin2_matrix_3_2;
real borigin2_matrix_3_4;
real borigin2_matrix_3_8;
real borigin2_matrix_3_18_DE;
real borigin2_matrix_3_18_NDE;
real borigin2_matrix_3_20_DE;
real borigin2_matrix_3_20_NDE;
real borigin2_matrix_4_3;
real borigin2_matrix_8_3;
real borigin2_matrix_10_2;
real borigin2_matrix_12_2;
real borigin2_matrix_14_2;
real borigin2_matrix_16_2;
real borigin2_matrix_18_3;
real borigin2_matrix_20_3;
real borigin2_matrix_39_2;

real borigin3_matrix_1_2;
real borigin3_matrix_2_1;
real borigin3_matrix_2_3;
real borigin3_matrix_2_10_DE;
real borigin3_matrix_2_10_NDE;
real borigin3_matrix_2_12_DE;
real borigin3_matrix_2_12_NDE;
real borigin3_matrix_2_14_DE;
real borigin3_matrix_2_14_NDE;
real borigin3_matrix_2_16_DE;
real borigin3_matrix_2_16_NDE;
real borigin3_matrix_2_39;
real borigin3_matrix_3_2;
real borigin3_matrix_3_4;
real borigin3_matrix_3_8;
real borigin3_matrix_3_18_DE;
real borigin3_matrix_3_18_NDE;
real borigin3_matrix_3_20_DE;
real borigin3_matrix_3_20_NDE;
real borigin3_matrix_4_3;
real borigin3_matrix_8_3;
real borigin3_matrix_10_2;
real borigin3_matrix_12_2;
real borigin3_matrix_14_2;
real borigin3_matrix_16_2;
real borigin3_matrix_18_3;
real borigin3_matrix_20_3;
real borigin3_matrix_39_2;

real borigin4_matrix_1_2;
real borigin4_matrix_2_1;
real borigin4_matrix_2_3;
real borigin4_matrix_2_10_DE;
real borigin4_matrix_2_10_NDE;
real borigin4_matrix_2_12_DE;
real borigin4_matrix_2_12_NDE;
real borigin4_matrix_2_14_DE;
real borigin4_matrix_2_14_NDE;
real borigin4_matrix_2_16_DE;
real borigin4_matrix_2_16_NDE;
real borigin4_matrix_2_39;
real borigin4_matrix_3_2;
real borigin4_matrix_3_4;
real borigin4_matrix_3_8;
real borigin4_matrix_3_18_DE;
real borigin4_matrix_3_18_NDE;
real borigin4_matrix_3_20_DE;
real borigin4_matrix_3_20_NDE;
real borigin4_matrix_4_3;
real borigin4_matrix_8_3;
real borigin4_matrix_10_2;
real borigin4_matrix_12_2;
real borigin4_matrix_14_2;
real borigin4_matrix_16_2;
real borigin4_matrix_18_3;
real borigin4_matrix_20_3;
real borigin4_matrix_39_2;

real borigin5_matrix_1_2;
real borigin5_matrix_2_1;
real borigin5_matrix_2_3;
real borigin5_matrix_2_10_DE;
real borigin5_matrix_2_10_NDE;
real borigin5_matrix_2_12_DE;
real borigin5_matrix_2_12_NDE;
real borigin5_matrix_2_14_DE;
real borigin5_matrix_2_14_NDE;
real borigin5_matrix_2_16_DE;
real borigin5_matrix_2_16_NDE;
real borigin5_matrix_2_39;
real borigin5_matrix_3_2;
real borigin5_matrix_3_4;
real borigin5_matrix_3_8;
real borigin5_matrix_3_18_DE;
real borigin5_matrix_3_18_NDE;
real borigin5_matrix_3_20_DE;
real borigin5_matrix_3_20_NDE;
real borigin5_matrix_4_3;
real borigin5_matrix_8_3;
real borigin5_matrix_10_2;
real borigin5_matrix_12_2;
real borigin5_matrix_14_2;
real borigin5_matrix_16_2;
real borigin5_matrix_18_3;
real borigin5_matrix_20_3;
real borigin5_matrix_39_2;

real borigin6_matrix_1_2;
real borigin6_matrix_2_1;
real borigin6_matrix_2_3;
real borigin6_matrix_2_10_DE;
real borigin6_matrix_2_10_NDE;
real borigin6_matrix_2_12_DE;
real borigin6_matrix_2_12_NDE;
real borigin6_matrix_2_14_DE;
real borigin6_matrix_2_14_NDE;
real borigin6_matrix_2_16_DE;
real borigin6_matrix_2_16_NDE;
real borigin6_matrix_2_39;
real borigin6_matrix_3_2;
real borigin6_matrix_3_4;
real borigin6_matrix_3_8;
real borigin6_matrix_3_18_DE;
real borigin6_matrix_3_18_NDE;
real borigin6_matrix_3_20_DE;
real borigin6_matrix_3_20_NDE;
real borigin6_matrix_4_3;
real borigin6_matrix_8_3;
real borigin6_matrix_10_2;
real borigin6_matrix_12_2;
real borigin6_matrix_14_2;
real borigin6_matrix_16_2;
real borigin6_matrix_18_3;
real borigin6_matrix_20_3;
real borigin6_matrix_39_2;

// random effect of year parameters
// vectors to store byear parameters for each movement
vector[nyears] byearxorigin1_raw_vector_1_2;
vector[nyears] byearxorigin1_raw_vector_2_1;
vector[nyears] byearxorigin1_raw_vector_2_3;
vector[nyears] byearxorigin1_raw_vector_2_10_DE;
vector[nyears] byearxorigin1_raw_vector_2_10_NDE;
vector[nyears] byearxorigin1_raw_vector_2_12_DE;
vector[nyears] byearxorigin1_raw_vector_2_12_NDE;
vector[nyears] byearxorigin1_raw_vector_2_14_DE;
vector[nyears] byearxorigin1_raw_vector_2_14_NDE;
vector[nyears] byearxorigin1_raw_vector_2_16_DE;
vector[nyears] byearxorigin1_raw_vector_2_16_NDE;
vector[nyears] byearxorigin1_raw_vector_2_39;
vector[nyears] byearxorigin1_raw_vector_3_2;
vector[nyears] byearxorigin1_raw_vector_3_4;
vector[nyears] byearxorigin1_raw_vector_3_8;
vector[nyears] byearxorigin1_raw_vector_3_18_DE;
vector[nyears] byearxorigin1_raw_vector_3_18_NDE;
vector[nyears] byearxorigin1_raw_vector_3_20_DE;
vector[nyears] byearxorigin1_raw_vector_3_20_NDE;
vector[nyears] byearxorigin1_raw_vector_4_3;
vector[nyears] byearxorigin1_raw_vector_8_3;

vector[nyears] byearxorigin2_raw_vector_1_2;
vector[nyears] byearxorigin2_raw_vector_2_1;
vector[nyears] byearxorigin2_raw_vector_2_3;
vector[nyears] byearxorigin2_raw_vector_2_10_DE;
vector[nyears] byearxorigin2_raw_vector_2_10_NDE;
vector[nyears] byearxorigin2_raw_vector_2_12_DE;
vector[nyears] byearxorigin2_raw_vector_2_12_NDE;
vector[nyears] byearxorigin2_raw_vector_2_14_DE;
vector[nyears] byearxorigin2_raw_vector_2_14_NDE;
vector[nyears] byearxorigin2_raw_vector_2_16_DE;
vector[nyears] byearxorigin2_raw_vector_2_16_NDE;
vector[nyears] byearxorigin2_raw_vector_2_39;
vector[nyears] byearxorigin2_raw_vector_3_2;
vector[nyears] byearxorigin2_raw_vector_3_4;
vector[nyears] byearxorigin2_raw_vector_3_8;
vector[nyears] byearxorigin2_raw_vector_3_18_DE;
vector[nyears] byearxorigin2_raw_vector_3_18_NDE;
vector[nyears] byearxorigin2_raw_vector_3_20_DE;
vector[nyears] byearxorigin2_raw_vector_3_20_NDE;
vector[nyears] byearxorigin2_raw_vector_4_3;
vector[nyears] byearxorigin2_raw_vector_8_3;

vector[nyears] byearxorigin3_raw_vector_1_2;
vector[nyears] byearxorigin3_raw_vector_2_1;
vector[nyears] byearxorigin3_raw_vector_2_3;
vector[nyears] byearxorigin3_raw_vector_2_10_DE;
vector[nyears] byearxorigin3_raw_vector_2_10_NDE;
vector[nyears] byearxorigin3_raw_vector_2_12_DE;
vector[nyears] byearxorigin3_raw_vector_2_12_NDE;
vector[nyears] byearxorigin3_raw_vector_2_14_DE;
vector[nyears] byearxorigin3_raw_vector_2_14_NDE;
vector[nyears] byearxorigin3_raw_vector_2_16_DE;
vector[nyears] byearxorigin3_raw_vector_2_16_NDE;
vector[nyears] byearxorigin3_raw_vector_2_39;
vector[nyears] byearxorigin3_raw_vector_3_2;
vector[nyears] byearxorigin3_raw_vector_3_4;
vector[nyears] byearxorigin3_raw_vector_3_8;
vector[nyears] byearxorigin3_raw_vector_3_18_DE;
vector[nyears] byearxorigin3_raw_vector_3_18_NDE;
vector[nyears] byearxorigin3_raw_vector_3_20_DE;
vector[nyears] byearxorigin3_raw_vector_3_20_NDE;
vector[nyears] byearxorigin3_raw_vector_4_3;
vector[nyears] byearxorigin3_raw_vector_8_3;

vector[nyears] byearxorigin4_raw_vector_1_2;
vector[nyears] byearxorigin4_raw_vector_2_1;
vector[nyears] byearxorigin4_raw_vector_2_3;
vector[nyears] byearxorigin4_raw_vector_2_10_DE;
vector[nyears] byearxorigin4_raw_vector_2_10_NDE;
vector[nyears] byearxorigin4_raw_vector_2_12_DE;
vector[nyears] byearxorigin4_raw_vector_2_12_NDE;
vector[nyears] byearxorigin4_raw_vector_2_14_DE;
vector[nyears] byearxorigin4_raw_vector_2_14_NDE;
vector[nyears] byearxorigin4_raw_vector_2_16_DE;
vector[nyears] byearxorigin4_raw_vector_2_16_NDE;
vector[nyears] byearxorigin4_raw_vector_2_39;
vector[nyears] byearxorigin4_raw_vector_3_2;
vector[nyears] byearxorigin4_raw_vector_3_4;
vector[nyears] byearxorigin4_raw_vector_3_8;
vector[nyears] byearxorigin4_raw_vector_3_18_DE;
vector[nyears] byearxorigin4_raw_vector_3_18_NDE;
vector[nyears] byearxorigin4_raw_vector_3_20_DE;
vector[nyears] byearxorigin4_raw_vector_3_20_NDE;
vector[nyears] byearxorigin4_raw_vector_4_3;
vector[nyears] byearxorigin4_raw_vector_8_3;

vector[nyears] byearxorigin5_raw_vector_1_2;
vector[nyears] byearxorigin5_raw_vector_2_1;
vector[nyears] byearxorigin5_raw_vector_2_3;
vector[nyears] byearxorigin5_raw_vector_2_10_DE;
vector[nyears] byearxorigin5_raw_vector_2_10_NDE;
vector[nyears] byearxorigin5_raw_vector_2_12_DE;
vector[nyears] byearxorigin5_raw_vector_2_12_NDE;
vector[nyears] byearxorigin5_raw_vector_2_14_DE;
vector[nyears] byearxorigin5_raw_vector_2_14_NDE;
vector[nyears] byearxorigin5_raw_vector_2_16_DE;
vector[nyears] byearxorigin5_raw_vector_2_16_NDE;
vector[nyears] byearxorigin5_raw_vector_2_39;
vector[nyears] byearxorigin5_raw_vector_3_2;
vector[nyears] byearxorigin5_raw_vector_3_4;
vector[nyears] byearxorigin5_raw_vector_3_8;
vector[nyears] byearxorigin5_raw_vector_3_18_DE;
vector[nyears] byearxorigin5_raw_vector_3_18_NDE;
vector[nyears] byearxorigin5_raw_vector_3_20_DE;
vector[nyears] byearxorigin5_raw_vector_3_20_NDE;
vector[nyears] byearxorigin5_raw_vector_4_3;
vector[nyears] byearxorigin5_raw_vector_8_3;

vector[nyears] byearxorigin6_raw_vector_1_2;
vector[nyears] byearxorigin6_raw_vector_2_1;
vector[nyears] byearxorigin6_raw_vector_2_3;
vector[nyears] byearxorigin6_raw_vector_2_10_DE;
vector[nyears] byearxorigin6_raw_vector_2_10_NDE;
vector[nyears] byearxorigin6_raw_vector_2_12_DE;
vector[nyears] byearxorigin6_raw_vector_2_12_NDE;
vector[nyears] byearxorigin6_raw_vector_2_14_DE;
vector[nyears] byearxorigin6_raw_vector_2_14_NDE;
vector[nyears] byearxorigin6_raw_vector_2_16_DE;
vector[nyears] byearxorigin6_raw_vector_2_16_NDE;
vector[nyears] byearxorigin6_raw_vector_2_39;
vector[nyears] byearxorigin6_raw_vector_3_2;
vector[nyears] byearxorigin6_raw_vector_3_4;
vector[nyears] byearxorigin6_raw_vector_3_8;
vector[nyears] byearxorigin6_raw_vector_3_18_DE;
vector[nyears] byearxorigin6_raw_vector_3_18_NDE;
vector[nyears] byearxorigin6_raw_vector_3_20_DE;
vector[nyears] byearxorigin6_raw_vector_3_20_NDE;
vector[nyears] byearxorigin6_raw_vector_4_3;
vector[nyears] byearxorigin6_raw_vector_8_3;

real<lower=0> sigma_yearxorigin1_vector_1_2;
real<lower=0> sigma_yearxorigin1_vector_2_1;
real<lower=0> sigma_yearxorigin1_vector_2_3;
real<lower=0> sigma_yearxorigin1_vector_2_10_DE;
real<lower=0> sigma_yearxorigin1_vector_2_10_NDE;
real<lower=0> sigma_yearxorigin1_vector_2_12_DE;
real<lower=0> sigma_yearxorigin1_vector_2_12_NDE;
real<lower=0> sigma_yearxorigin1_vector_2_14_DE;
real<lower=0> sigma_yearxorigin1_vector_2_14_NDE;
real<lower=0> sigma_yearxorigin1_vector_2_16_DE;
real<lower=0> sigma_yearxorigin1_vector_2_16_NDE;
real<lower=0> sigma_yearxorigin1_vector_2_39;
// real<lower=0> sigma_yearxorigin1_vector_3_2_DE;
real<lower=0> sigma_yearxorigin1_vector_3_4;
real<lower=0> sigma_yearxorigin1_vector_3_8;
real<lower=0> sigma_yearxorigin1_vector_3_18_DE;
real<lower=0> sigma_yearxorigin1_vector_3_18_NDE;
real<lower=0> sigma_yearxorigin1_vector_3_20_DE;
real<lower=0> sigma_yearxorigin1_vector_3_20_NDE;
// real<lower=0> sigma_yearxorigin1_vector_4_3_DE;
real<lower=0> sigma_yearxorigin1_vector_8_3;

real<lower=0> sigma_yearxorigin2_vector_1_2;
real<lower=0> sigma_yearxorigin2_vector_2_1;
real<lower=0> sigma_yearxorigin2_vector_2_3;
real<lower=0> sigma_yearxorigin2_vector_2_10_DE;
real<lower=0> sigma_yearxorigin2_vector_2_10_NDE;
real<lower=0> sigma_yearxorigin2_vector_2_12_DE;
real<lower=0> sigma_yearxorigin2_vector_2_12_NDE;
real<lower=0> sigma_yearxorigin2_vector_2_14_DE;
real<lower=0> sigma_yearxorigin2_vector_2_14_NDE;
real<lower=0> sigma_yearxorigin2_vector_2_16_DE;
real<lower=0> sigma_yearxorigin2_vector_2_16_NDE;
real<lower=0> sigma_yearxorigin2_vector_2_39;
// real<lower=0> sigma_yearxorigin2_vector_3_2_DE;
real<lower=0> sigma_yearxorigin2_vector_3_4;
real<lower=0> sigma_yearxorigin2_vector_3_8;
real<lower=0> sigma_yearxorigin2_vector_3_18_DE;
real<lower=0> sigma_yearxorigin2_vector_3_18_NDE;
real<lower=0> sigma_yearxorigin2_vector_3_20_DE;
real<lower=0> sigma_yearxorigin2_vector_3_20_NDE;
// real<lower=0> sigma_yearxorigin2_vector_4_3_DE;
real<lower=0> sigma_yearxorigin2_vector_8_3;

real<lower=0> sigma_yearxorigin3_vector_1_2;
real<lower=0> sigma_yearxorigin3_vector_2_1;
real<lower=0> sigma_yearxorigin3_vector_2_3;
real<lower=0> sigma_yearxorigin3_vector_2_10_DE;
real<lower=0> sigma_yearxorigin3_vector_2_10_NDE;
real<lower=0> sigma_yearxorigin3_vector_2_12_DE;
real<lower=0> sigma_yearxorigin3_vector_2_12_NDE;
real<lower=0> sigma_yearxorigin3_vector_2_14_DE;
real<lower=0> sigma_yearxorigin3_vector_2_14_NDE;
real<lower=0> sigma_yearxorigin3_vector_2_16_DE;
real<lower=0> sigma_yearxorigin3_vector_2_16_NDE;
real<lower=0> sigma_yearxorigin3_vector_2_39;
// real<lower=0> sigma_yearxorigin3_vector_3_2_DE;
real<lower=0> sigma_yearxorigin3_vector_3_4;
real<lower=0> sigma_yearxorigin3_vector_3_8;
real<lower=0> sigma_yearxorigin3_vector_3_18_DE;
real<lower=0> sigma_yearxorigin3_vector_3_18_NDE;
real<lower=0> sigma_yearxorigin3_vector_3_20_DE;
real<lower=0> sigma_yearxorigin3_vector_3_20_NDE;
// real<lower=0> sigma_yearxorigin3_vector_4_3_DE;
real<lower=0> sigma_yearxorigin3_vector_8_3;

real<lower=0> sigma_yearxorigin4_vector_1_2;
real<lower=0> sigma_yearxorigin4_vector_2_1;
real<lower=0> sigma_yearxorigin4_vector_2_3;
real<lower=0> sigma_yearxorigin4_vector_2_10_DE;
real<lower=0> sigma_yearxorigin4_vector_2_10_NDE;
real<lower=0> sigma_yearxorigin4_vector_2_12_DE;
real<lower=0> sigma_yearxorigin4_vector_2_12_NDE;
real<lower=0> sigma_yearxorigin4_vector_2_14_DE;
real<lower=0> sigma_yearxorigin4_vector_2_14_NDE;
real<lower=0> sigma_yearxorigin4_vector_2_16_DE;
real<lower=0> sigma_yearxorigin4_vector_2_16_NDE;
real<lower=0> sigma_yearxorigin4_vector_2_39;
// real<lower=0> sigma_yearxorigin4_vector_3_2_DE;
real<lower=0> sigma_yearxorigin4_vector_3_4;
real<lower=0> sigma_yearxorigin4_vector_3_8;
real<lower=0> sigma_yearxorigin4_vector_3_18_DE;
real<lower=0> sigma_yearxorigin4_vector_3_18_NDE;
real<lower=0> sigma_yearxorigin4_vector_3_20_DE;
real<lower=0> sigma_yearxorigin4_vector_3_20_NDE;
// real<lower=0> sigma_yearxorigin4_vector_4_3_DE;
real<lower=0> sigma_yearxorigin4_vector_8_3;

real<lower=0> sigma_yearxorigin5_vector_1_2;
real<lower=0> sigma_yearxorigin5_vector_2_1;
real<lower=0> sigma_yearxorigin5_vector_2_3;
real<lower=0> sigma_yearxorigin5_vector_2_10_DE;
real<lower=0> sigma_yearxorigin5_vector_2_10_NDE;
real<lower=0> sigma_yearxorigin5_vector_2_12_DE;
real<lower=0> sigma_yearxorigin5_vector_2_12_NDE;
real<lower=0> sigma_yearxorigin5_vector_2_14_DE;
real<lower=0> sigma_yearxorigin5_vector_2_14_NDE;
real<lower=0> sigma_yearxorigin5_vector_2_16_DE;
real<lower=0> sigma_yearxorigin5_vector_2_16_NDE;
real<lower=0> sigma_yearxorigin5_vector_2_39;
real<lower=0> sigma_yearxorigin5_vector_3_2;
real<lower=0> sigma_yearxorigin5_vector_3_4;
real<lower=0> sigma_yearxorigin5_vector_3_8;
real<lower=0> sigma_yearxorigin5_vector_3_18_DE;
real<lower=0> sigma_yearxorigin5_vector_3_18_NDE;
real<lower=0> sigma_yearxorigin5_vector_3_20_DE;
real<lower=0> sigma_yearxorigin5_vector_3_20_NDE;
// real<lower=0> sigma_yearxorigin5_vector_4_3_DE;
real<lower=0> sigma_yearxorigin5_vector_8_3;

real<lower=0> sigma_yearxorigin6_vector_1_2;
real<lower=0> sigma_yearxorigin6_vector_2_1;
real<lower=0> sigma_yearxorigin6_vector_2_3;
real<lower=0> sigma_yearxorigin6_vector_2_10_DE;
real<lower=0> sigma_yearxorigin6_vector_2_10_NDE;
real<lower=0> sigma_yearxorigin6_vector_2_12_DE;
real<lower=0> sigma_yearxorigin6_vector_2_12_NDE;
real<lower=0> sigma_yearxorigin6_vector_2_14_DE;
real<lower=0> sigma_yearxorigin6_vector_2_14_NDE;
real<lower=0> sigma_yearxorigin6_vector_2_16_DE;
real<lower=0> sigma_yearxorigin6_vector_2_16_NDE;
real<lower=0> sigma_yearxorigin6_vector_2_39;
real<lower=0> sigma_yearxorigin6_vector_3_2;
real<lower=0> sigma_yearxorigin6_vector_3_4;
real<lower=0> sigma_yearxorigin6_vector_3_8;
real<lower=0> sigma_yearxorigin6_vector_3_18_DE;
real<lower=0> sigma_yearxorigin6_vector_3_18_NDE;
real<lower=0> sigma_yearxorigin6_vector_3_20_DE;
real<lower=0> sigma_yearxorigin6_vector_3_20_NDE;
// real<lower=0> sigma_yearxorigin6_vector_4_3_DE;
real<lower=0> sigma_yearxorigin6_vector_8_3;


// here, write out all of the parameters for detection efficiency
// twenty terms for intercepts for different eras (configurations of antennas) in the different tributaries
real asotin_alpha1;
real asotin_alpha2;
real deschutes_alpha1;
real entiat_alpha1;
real fifteenmile_alpha1;
real imnaha_alpha1;
real john_day_alpha1;
real methow_alpha1;
real methow_alpha2;
real okanogan_alpha1;
real tucannon_alpha1;
real tucannon_alpha2;
real umatilla_alpha1;
real umatilla_alpha2;
real walla_walla_alpha1;
real walla_walla_alpha2;
real walla_walla_alpha3;
real walla_walla_alpha4;
real wenatchee_alpha1;
real yakima_alpha1;

// 14 terms for discharge relationship, one for each tributary
real asotin_beta;
real deschutes_beta;
real entiat_beta;
// real fifteenmile_beta;
// real imnaha_beta;
real john_day_beta;
real methow_beta;
real okanogan_beta;
real tucannon_beta;
real umatilla_beta;
real walla_walla_beta;
real wenatchee_beta;
real yakima_beta;

  
}

transformed parameters {
  // We now need two matrices always - one for DE params, and one for NDE params
  
  
  // Declare a matrix to store b0 params - intercept terms
  // Every allowable movement will be replaced by a parameter value. The -100000
  // (effectively a 0 in log space) forces impossible movements to be assigned zero probability.
  // Almost all will be overwritten, except in the case of some movements that are
  // allowable in NDE but not in DE (i.e., certain tributary movements)
  array[nmovements + 1] real b0_vector_DE;
  b0_vector_DE = rep_array(-100000, nmovements + 1);
  array[nmovements + 1] real b0_vector_NDE;
  b0_vector_NDE = rep_array(-100000, nmovements + 1);
  
  // declare matrices to store spill parameters
  
  // spill windows
  array[nmovements + 1] real bspillwindow_vector;
  bspillwindow_vector = rep_array(0, nmovements + 1);
  
  // days of spill across winter months
  array[nmovements + 1] real bwinterspill_vector;
  bwinterspill_vector = rep_array(0, nmovements + 1);
  
  
  // Declare vectors to store sigma_year params - these scale the year effects
  // These don't need the + 1 for length, because they aren't being directly used in
  // the calculation of movement probabilities
  array[nmovements] real sigma_yearxorigin1_vector_DE;
  sigma_yearxorigin1_vector_DE = rep_array(0,nmovements);
  array[nmovements] real sigma_yearxorigin1_vector_NDE;
  sigma_yearxorigin1_vector_NDE = rep_array(0,nmovements);

  array[nmovements] real sigma_yearxorigin2_vector_DE;
  sigma_yearxorigin2_vector_DE = rep_array(0,nmovements);
  array[nmovements] real sigma_yearxorigin2_vector_NDE;
  sigma_yearxorigin2_vector_NDE = rep_array(0,nmovements);
  
  array[nmovements] real sigma_yearxorigin3_vector_DE;
  sigma_yearxorigin3_vector_DE = rep_array(0,nmovements);
  array[nmovements] real sigma_yearxorigin3_vector_NDE;
  sigma_yearxorigin3_vector_NDE = rep_array(0,nmovements);
  
  array[nmovements] real sigma_yearxorigin4_vector_DE;
  sigma_yearxorigin4_vector_DE = rep_array(0,nmovements);
  array[nmovements] real sigma_yearxorigin4_vector_NDE;
  sigma_yearxorigin4_vector_NDE = rep_array(0,nmovements);
  
  array[nmovements] real sigma_yearxorigin5_vector_DE;
  sigma_yearxorigin5_vector_DE = rep_array(0,nmovements);
  array[nmovements] real sigma_yearxorigin5_vector_NDE;
  sigma_yearxorigin5_vector_NDE = rep_array(0,nmovements);
  
  array[nmovements] real sigma_yearxorigin6_vector_DE;
  sigma_yearxorigin6_vector_DE = rep_array(0,nmovements);
  array[nmovements] real sigma_yearxorigin6_vector_NDE;
  sigma_yearxorigin6_vector_NDE = rep_array(0,nmovements);

  // Declare matrices to contain the year effects (each year is a column) for each movement (rows)
  // First declare containers for the actual parameters (this is raw multiplied by scaling (sigma))
  // These year parameters are actually being accessed, so they need the +1
  array [nmovements + 1,nyears] real byearxorigin1_actual_parameters_matrix_DE;
  byearxorigin1_actual_parameters_matrix_DE = rep_array(0, nmovements + 1, nyears);
  array [nmovements + 1,nyears] real byearxorigin1_actual_parameters_matrix_NDE;
  byearxorigin1_actual_parameters_matrix_NDE = rep_array(0, nmovements + 1, nyears);
  
  array [nmovements + 1,nyears] real byearxorigin2_actual_parameters_matrix_DE;
  byearxorigin2_actual_parameters_matrix_DE = rep_array(0, nmovements + 1, nyears);
  array [nmovements + 1,nyears] real byearxorigin2_actual_parameters_matrix_NDE;
  byearxorigin2_actual_parameters_matrix_NDE = rep_array(0, nmovements + 1, nyears);
  
  array [nmovements + 1,nyears] real byearxorigin3_actual_parameters_matrix_DE;
  byearxorigin3_actual_parameters_matrix_DE = rep_array(0, nmovements + 1, nyears);
  array [nmovements + 1,nyears] real byearxorigin3_actual_parameters_matrix_NDE;
  byearxorigin3_actual_parameters_matrix_NDE = rep_array(0, nmovements + 1, nyears);
  
  array [nmovements + 1,nyears] real byearxorigin4_actual_parameters_matrix_DE;
  byearxorigin4_actual_parameters_matrix_DE = rep_array(0, nmovements + 1, nyears);
  array [nmovements + 1,nyears] real byearxorigin4_actual_parameters_matrix_NDE;
  byearxorigin4_actual_parameters_matrix_NDE = rep_array(0, nmovements + 1, nyears);
  
  array [nmovements + 1,nyears] real byearxorigin5_actual_parameters_matrix_DE;
  byearxorigin5_actual_parameters_matrix_DE = rep_array(0, nmovements + 1, nyears);
  array [nmovements + 1,nyears] real byearxorigin5_actual_parameters_matrix_NDE;
  byearxorigin5_actual_parameters_matrix_NDE = rep_array(0, nmovements + 1, nyears);
  
  array [nmovements + 1,nyears] real byearxorigin6_actual_parameters_matrix_DE;
  byearxorigin6_actual_parameters_matrix_DE = rep_array(0, nmovements + 1, nyears);
  array [nmovements + 1,nyears] real byearxorigin6_actual_parameters_matrix_NDE;
  byearxorigin6_actual_parameters_matrix_NDE = rep_array(0, nmovements + 1, nyears);
  
  // Then declare matrices to store raw year effects
  // don't need +1, they're not being accessed via parameter_indices_matrix
  array [nmovements,nyears] real byearxorigin1_raw_parameters_matrix_DE;
  byearxorigin1_raw_parameters_matrix_DE = rep_array(0, nmovements, nyears);
  array [nmovements,nyears] real byearxorigin1_raw_parameters_matrix_NDE;
  byearxorigin1_raw_parameters_matrix_NDE = rep_array(0, nmovements, nyears);

  array [nmovements,nyears] real byearxorigin2_raw_parameters_matrix_DE;
  byearxorigin2_raw_parameters_matrix_DE = rep_array(0, nmovements, nyears);
  array [nmovements,nyears] real byearxorigin2_raw_parameters_matrix_NDE;
  byearxorigin2_raw_parameters_matrix_NDE = rep_array(0, nmovements, nyears);
  
  array [nmovements,nyears] real byearxorigin3_raw_parameters_matrix_DE;
  byearxorigin3_raw_parameters_matrix_DE = rep_array(0, nmovements, nyears);
  array [nmovements,nyears] real byearxorigin3_raw_parameters_matrix_NDE;
  byearxorigin3_raw_parameters_matrix_NDE = rep_array(0, nmovements, nyears);
  
  array [nmovements,nyears] real byearxorigin4_raw_parameters_matrix_DE;
  byearxorigin4_raw_parameters_matrix_DE = rep_array(0, nmovements, nyears);
  array [nmovements,nyears] real byearxorigin4_raw_parameters_matrix_NDE;
  byearxorigin4_raw_parameters_matrix_NDE = rep_array(0, nmovements, nyears);
  
  array [nmovements,nyears] real byearxorigin5_raw_parameters_matrix_DE;
  byearxorigin5_raw_parameters_matrix_DE = rep_array(0, nmovements, nyears);
  array [nmovements,nyears] real byearxorigin5_raw_parameters_matrix_NDE;
  byearxorigin5_raw_parameters_matrix_NDE = rep_array(0, nmovements, nyears);
  
  array [nmovements,nyears] real byearxorigin6_raw_parameters_matrix_DE;
  byearxorigin6_raw_parameters_matrix_DE = rep_array(0, nmovements, nyears);
  array [nmovements,nyears] real byearxorigin6_raw_parameters_matrix_NDE;
  byearxorigin6_raw_parameters_matrix_NDE = rep_array(0, nmovements, nyears);
  
    // Declare a vector to store btemp0 params
  array[nmovements + 1] real btemp0_vector_DE;
  btemp0_vector_DE = rep_array(0, nmovements + 1);
  
  array[nmovements + 1] real btemp0_vector_NDE;
  btemp0_vector_NDE = rep_array(0, nmovements + 1);
  
    // Declare a vector to store btemp1 params
  array[nmovements + 1] real btemp1_vector_DE;
  btemp1_vector_DE = rep_array(0, nmovements + 1);
  
  array[nmovements + 1] real btemp1_vector_NDE;
  btemp1_vector_NDE = rep_array(0, nmovements + 1);
  
  // Declare matrices to store origin x temperature interactions (DE and NDE)
  // Origin 1 x temp0
  array[nmovements + 1] real btemp0xorigin1_vector_DE;
  btemp0xorigin1_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real btemp0xorigin1_vector_NDE;
  btemp0xorigin1_vector_NDE = rep_array(0, nmovements + 1);
  // Origin 2 x temp0
  array[nmovements + 1] real btemp0xorigin2_vector_DE;
  btemp0xorigin2_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real btemp0xorigin2_vector_NDE;
  btemp0xorigin2_vector_NDE = rep_array(0, nmovements + 1);
  // Origin 3 x temp0
  array[nmovements + 1] real btemp0xorigin3_vector_DE;
  btemp0xorigin3_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real btemp0xorigin3_vector_NDE;
  btemp0xorigin3_vector_NDE = rep_array(0, nmovements + 1);
  // Origin 4 x temp0
  array[nmovements + 1] real btemp0xorigin4_vector_DE;
  btemp0xorigin4_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real btemp0xorigin4_vector_NDE;
  btemp0xorigin4_vector_NDE = rep_array(0, nmovements + 1);
  // Origin 5 x temp0
  array[nmovements + 1] real btemp0xorigin5_vector_DE;
  btemp0xorigin5_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real btemp0xorigin5_vector_NDE;
  btemp0xorigin5_vector_NDE = rep_array(0, nmovements + 1);
  // Origin 6 x temp0
  array[nmovements + 1] real btemp0xorigin6_vector_DE;
  btemp0xorigin6_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real btemp0xorigin6_vector_NDE;
  btemp0xorigin6_vector_NDE = rep_array(0, nmovements + 1);
  
  
  // Origin 1 x temp1
  array[nmovements + 1] real btemp1xorigin1_vector_DE;
  btemp1xorigin1_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real btemp1xorigin1_vector_NDE;
  btemp1xorigin1_vector_NDE = rep_array(0, nmovements + 1);
  // Origin 2 x temp1
  array[nmovements + 1] real btemp1xorigin2_vector_DE;
  btemp1xorigin2_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real btemp1xorigin2_vector_NDE;
  btemp1xorigin2_vector_NDE = rep_array(0, nmovements + 1);
  // Origin 3 x temp1
  array[nmovements + 1] real btemp1xorigin3_vector_DE;
  btemp1xorigin3_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real btemp1xorigin3_vector_NDE;
  btemp1xorigin3_vector_NDE = rep_array(0, nmovements + 1);
  // Origin 4 x temp1
  array[nmovements + 1] real btemp1xorigin4_vector_DE;
  btemp1xorigin4_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real btemp1xorigin4_vector_NDE;
  btemp1xorigin4_vector_NDE = rep_array(0, nmovements + 1);
  // Origin 5 x temp1
  array[nmovements + 1] real btemp1xorigin5_vector_DE;
  btemp1xorigin5_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real btemp1xorigin5_vector_NDE;
  btemp1xorigin5_vector_NDE = rep_array(0, nmovements + 1);
  // Origin 6 x temp1
  array[nmovements + 1] real btemp1xorigin6_vector_DE;
  btemp1xorigin6_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real btemp1xorigin6_vector_NDE;
  btemp1xorigin6_vector_NDE = rep_array(0, nmovements + 1);  
  
  // Declare a vector to store borigin1 params
  array[nmovements + 1] real borigin1_vector_DE;
  borigin1_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real borigin1_vector_NDE;
  borigin1_vector_NDE = rep_array(0, nmovements + 1);
  // Declare a vector to store borigin2 params
  array[nmovements + 1] real borigin2_vector_DE;
  borigin2_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real borigin2_vector_NDE;
  borigin2_vector_NDE = rep_array(0, nmovements + 1);
  // Declare a vector to store borigin3 params
  array[nmovements + 1] real borigin3_vector_DE;
  borigin3_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real borigin3_vector_NDE;
  borigin3_vector_NDE = rep_array(0, nmovements + 1);
  // Declare a vector to store borigin4 params
  array[nmovements + 1] real borigin4_vector_DE;
  borigin4_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real borigin4_vector_NDE;
  borigin4_vector_NDE = rep_array(0, nmovements + 1);
  // Declare a vector to store borigin5 params
  array[nmovements + 1] real borigin5_vector_DE;
  borigin5_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real borigin5_vector_NDE;
  borigin5_vector_NDE = rep_array(0, nmovements + 1);
  // Declare a vector to store borigin6 params
  array[nmovements + 1] real borigin6_vector_DE;
  borigin6_vector_DE = rep_array(0, nmovements + 1);
  array[nmovements + 1] real borigin6_vector_NDE;
  borigin6_vector_NDE = rep_array(0, nmovements + 1);
  
  // Populate each of those containers (vectors and arrays) with the actual
  // parameters, using indexing

// Intercept terms
b0_vector_DE[18] = b0_matrix_4_5;
b0_vector_NDE[18] = b0_matrix_4_5;
b0_vector_DE[19] = b0_matrix_5_4;
b0_vector_NDE[19] = b0_matrix_5_4;
b0_vector_DE[20] = b0_matrix_5_6;
b0_vector_NDE[20] = b0_matrix_5_6;
b0_vector_DE[21] = b0_matrix_5_22_DE;
b0_vector_NDE[21] = b0_matrix_5_22_NDE;
b0_vector_DE[22] = b0_matrix_6_5;
b0_vector_NDE[22] = b0_matrix_6_5;
b0_vector_DE[23] = b0_matrix_6_7;
b0_vector_NDE[23] = b0_matrix_6_7;
b0_vector_DE[24] = b0_matrix_7_6;
b0_vector_NDE[24] = b0_matrix_7_6;
b0_vector_DE[25] = b0_matrix_7_26_DE;
b0_vector_NDE[25] = b0_matrix_7_26_NDE;
b0_vector_DE[26] = b0_matrix_7_28_DE;
b0_vector_NDE[26] = b0_matrix_7_28_NDE;
b0_vector_DE[28] = b0_matrix_8_9;
b0_vector_NDE[28] = b0_matrix_8_9;
b0_vector_DE[29] = b0_matrix_8_30_DE;
b0_vector_NDE[29] = b0_matrix_8_30_NDE;
b0_vector_DE[30] = b0_matrix_9_8;
b0_vector_NDE[30] = b0_matrix_9_8;
b0_vector_DE[31] = b0_matrix_9_32_DE;
b0_vector_NDE[31] = b0_matrix_9_32_NDE;
b0_vector_DE[32] = b0_matrix_9_34;
b0_vector_NDE[32] = b0_matrix_9_34;
b0_vector_DE[33] = b0_matrix_9_35;
b0_vector_NDE[33] = b0_matrix_9_35;
b0_vector_DE[34] = b0_matrix_9_36;
b0_vector_NDE[34] = b0_matrix_9_36;
b0_vector_DE[35] = b0_matrix_9_37_DE;
b0_vector_NDE[35] = b0_matrix_9_37_NDE;
b0_vector_DE[45] = b0_matrix_22_5;
b0_vector_NDE[45] = b0_matrix_22_5;
b0_vector_DE[46] = b0_matrix_26_7;
b0_vector_NDE[46] = b0_matrix_26_7;
b0_vector_DE[47] = b0_matrix_28_7;
b0_vector_NDE[47] = b0_matrix_28_7;
b0_vector_DE[48] = b0_matrix_30_8;
b0_vector_NDE[48] = b0_matrix_30_8;
b0_vector_DE[49] = b0_matrix_32_9;
b0_vector_NDE[49] = b0_matrix_32_9;
b0_vector_DE[50] = b0_matrix_34_9;
b0_vector_NDE[50] = b0_matrix_34_9;
b0_vector_DE[51] = b0_matrix_35_9;
b0_vector_NDE[51] = b0_matrix_35_9;
b0_vector_DE[52] = b0_matrix_36_9;
b0_vector_NDE[52] = b0_matrix_36_9;
b0_vector_DE[53] = b0_matrix_37_9;
b0_vector_NDE[53] = b0_matrix_37_9;
b0_vector_DE[1] = 0;
b0_vector_NDE[1] = 0;
b0_vector_DE[2] = 0;
b0_vector_NDE[2] = 0;
b0_vector_DE[3] = 0;
b0_vector_NDE[3] = 0;
b0_vector_DE[4] = 0;
b0_vector_NDE[4] = 0;
b0_vector_DE[5] = 0;
b0_vector_NDE[5] = 0;
b0_vector_DE[6] = 0;
b0_vector_NDE[6] = 0;
b0_vector_DE[7] = 0;
b0_vector_NDE[7] = 0;
b0_vector_DE[8] = 0;
b0_vector_NDE[8] = 0;
b0_vector_DE[9] = 0;
b0_vector_NDE[9] = 0;
b0_vector_DE[10] = 0;
b0_vector_NDE[10] = 0;
b0_vector_DE[11] = 0;
b0_vector_NDE[11] = 0;
b0_vector_DE[12] = 0;
b0_vector_NDE[12] = 0;
b0_vector_DE[13] = 0;
b0_vector_NDE[13] = 0;
b0_vector_DE[14] = 0;
b0_vector_NDE[14] = 0;
b0_vector_DE[15] = 0;
b0_vector_NDE[15] = 0;
b0_vector_DE[16] = 0;
b0_vector_NDE[16] = 0;
b0_vector_DE[17] = 0;
b0_vector_NDE[17] = 0;
b0_vector_DE[27] = 0;
b0_vector_NDE[27] = 0;
b0_vector_DE[36] = 0;
b0_vector_NDE[36] = 0;
b0_vector_DE[37] = 0;
b0_vector_NDE[37] = 0;
b0_vector_DE[38] = 0;
b0_vector_NDE[38] = 0;
b0_vector_DE[39] = 0;
b0_vector_NDE[39] = 0;
b0_vector_DE[40] = 0;
b0_vector_NDE[40] = 0;
b0_vector_DE[41] = 0;
b0_vector_NDE[41] = 0;
b0_vector_DE[42] = 0;
b0_vector_NDE[42] = 0;
b0_vector_DE[43] = 0;
b0_vector_NDE[43] = 0;
b0_vector_DE[44] = 0;
b0_vector_NDE[44] = 0;
b0_vector_DE[54] = 0;
b0_vector_NDE[54] = 0;
  
// Spill window
bspillwindow_vector[2] = bspillwindow_matrix_2_1;
bspillwindow_vector[12] = bspillwindow_matrix_3_2;
bspillwindow_vector[17] = bspillwindow_matrix_4_3;
bspillwindow_vector[19] = bspillwindow_matrix_5_4;
bspillwindow_vector[22] = bspillwindow_matrix_6_5;
bspillwindow_vector[24] = bspillwindow_matrix_7_6;
bspillwindow_vector[27] = bspillwindow_matrix_8_3;
bspillwindow_vector[30] = bspillwindow_matrix_9_8;

// Winter spill days
bwinterspill_vector[12] = bwinterspill_matrix_3_2;
bwinterspill_vector[17] = bwinterspill_matrix_4_3;
bwinterspill_vector[19] = bwinterspill_matrix_5_4;
bwinterspill_vector[22] = bwinterspill_matrix_6_5;
bwinterspill_vector[24] = bwinterspill_matrix_7_6;
bwinterspill_vector[27] = bwinterspill_matrix_8_3;
bwinterspill_vector[30] = bwinterspill_matrix_9_8;
  
// Year random effects - three blocks

// 1. store the byear_raw_vectors in the array
byearxorigin1_raw_parameters_matrix_DE[1, ] = to_array_1d(byearxorigin1_raw_vector_1_2);
byearxorigin1_raw_parameters_matrix_NDE[1, ] = to_array_1d(byearxorigin1_raw_vector_1_2);
byearxorigin1_raw_parameters_matrix_DE[2, ] = to_array_1d(byearxorigin1_raw_vector_2_1);
byearxorigin1_raw_parameters_matrix_NDE[2, ] = to_array_1d(byearxorigin1_raw_vector_2_1);
byearxorigin1_raw_parameters_matrix_DE[3, ] = to_array_1d(byearxorigin1_raw_vector_2_3);
byearxorigin1_raw_parameters_matrix_NDE[3, ] = to_array_1d(byearxorigin1_raw_vector_2_3);
byearxorigin1_raw_parameters_matrix_DE[4, ] = to_array_1d(byearxorigin1_raw_vector_2_10_DE);
byearxorigin1_raw_parameters_matrix_NDE[4, ] = to_array_1d(byearxorigin1_raw_vector_2_10_NDE);
byearxorigin1_raw_parameters_matrix_NDE[5, ] = to_array_1d(byearxorigin1_raw_vector_2_10_NDE);
byearxorigin1_raw_parameters_matrix_DE[6, ] = to_array_1d(byearxorigin1_raw_vector_2_12_DE);
byearxorigin1_raw_parameters_matrix_NDE[6, ] = to_array_1d(byearxorigin1_raw_vector_2_12_NDE);
byearxorigin1_raw_parameters_matrix_DE[7, ] = to_array_1d(byearxorigin1_raw_vector_2_14_DE);
byearxorigin1_raw_parameters_matrix_NDE[7, ] = to_array_1d(byearxorigin1_raw_vector_2_14_NDE);
byearxorigin1_raw_parameters_matrix_NDE[8, ] = to_array_1d(byearxorigin1_raw_vector_2_14_NDE);
byearxorigin1_raw_parameters_matrix_DE[9, ] = to_array_1d(byearxorigin1_raw_vector_2_16_DE);
byearxorigin1_raw_parameters_matrix_NDE[9, ] = to_array_1d(byearxorigin1_raw_vector_2_16_NDE);
byearxorigin1_raw_parameters_matrix_NDE[10, ] = to_array_1d(byearxorigin1_raw_vector_2_16_NDE);
byearxorigin1_raw_parameters_matrix_DE[11, ] = to_array_1d(byearxorigin1_raw_vector_2_39);
byearxorigin1_raw_parameters_matrix_NDE[11, ] = to_array_1d(byearxorigin1_raw_vector_2_39);
byearxorigin1_raw_parameters_matrix_DE[12, ] = to_array_1d(byearxorigin1_raw_vector_3_2);
byearxorigin1_raw_parameters_matrix_NDE[12, ] = to_array_1d(byearxorigin1_raw_vector_3_2);
byearxorigin1_raw_parameters_matrix_DE[13, ] = to_array_1d(byearxorigin1_raw_vector_3_4);
byearxorigin1_raw_parameters_matrix_NDE[13, ] = to_array_1d(byearxorigin1_raw_vector_3_4);
byearxorigin1_raw_parameters_matrix_DE[14, ] = to_array_1d(byearxorigin1_raw_vector_3_8);
byearxorigin1_raw_parameters_matrix_NDE[14, ] = to_array_1d(byearxorigin1_raw_vector_3_8);
byearxorigin1_raw_parameters_matrix_DE[15, ] = to_array_1d(byearxorigin1_raw_vector_3_18_DE);
byearxorigin1_raw_parameters_matrix_NDE[15, ] = to_array_1d(byearxorigin1_raw_vector_3_18_NDE);
byearxorigin1_raw_parameters_matrix_DE[16, ] = to_array_1d(byearxorigin1_raw_vector_3_20_DE);
byearxorigin1_raw_parameters_matrix_NDE[16, ] = to_array_1d(byearxorigin1_raw_vector_3_20_NDE);
byearxorigin1_raw_parameters_matrix_DE[17, ] = to_array_1d(byearxorigin1_raw_vector_4_3);
byearxorigin1_raw_parameters_matrix_NDE[17, ] = to_array_1d(byearxorigin1_raw_vector_4_3);
byearxorigin1_raw_parameters_matrix_DE[27, ] = to_array_1d(byearxorigin1_raw_vector_8_3);
byearxorigin1_raw_parameters_matrix_NDE[27, ] = to_array_1d(byearxorigin1_raw_vector_8_3);

byearxorigin2_raw_parameters_matrix_DE[1, ] = to_array_1d(byearxorigin2_raw_vector_1_2);
byearxorigin2_raw_parameters_matrix_NDE[1, ] = to_array_1d(byearxorigin2_raw_vector_1_2);
byearxorigin2_raw_parameters_matrix_DE[2, ] = to_array_1d(byearxorigin2_raw_vector_2_1);
byearxorigin2_raw_parameters_matrix_NDE[2, ] = to_array_1d(byearxorigin2_raw_vector_2_1);
byearxorigin2_raw_parameters_matrix_DE[3, ] = to_array_1d(byearxorigin2_raw_vector_2_3);
byearxorigin2_raw_parameters_matrix_NDE[3, ] = to_array_1d(byearxorigin2_raw_vector_2_3);
byearxorigin2_raw_parameters_matrix_DE[4, ] = to_array_1d(byearxorigin2_raw_vector_2_10_DE);
byearxorigin2_raw_parameters_matrix_NDE[4, ] = to_array_1d(byearxorigin2_raw_vector_2_10_NDE);
byearxorigin2_raw_parameters_matrix_NDE[5, ] = to_array_1d(byearxorigin2_raw_vector_2_10_NDE);
byearxorigin2_raw_parameters_matrix_DE[6, ] = to_array_1d(byearxorigin2_raw_vector_2_12_DE);
byearxorigin2_raw_parameters_matrix_NDE[6, ] = to_array_1d(byearxorigin2_raw_vector_2_12_NDE);
byearxorigin2_raw_parameters_matrix_DE[7, ] = to_array_1d(byearxorigin2_raw_vector_2_14_DE);
byearxorigin2_raw_parameters_matrix_NDE[7, ] = to_array_1d(byearxorigin2_raw_vector_2_14_NDE);
byearxorigin2_raw_parameters_matrix_NDE[8, ] = to_array_1d(byearxorigin2_raw_vector_2_14_NDE);
byearxorigin2_raw_parameters_matrix_DE[9, ] = to_array_1d(byearxorigin2_raw_vector_2_16_DE);
byearxorigin2_raw_parameters_matrix_NDE[9, ] = to_array_1d(byearxorigin2_raw_vector_2_16_NDE);
byearxorigin2_raw_parameters_matrix_NDE[10, ] = to_array_1d(byearxorigin2_raw_vector_2_16_NDE);
byearxorigin2_raw_parameters_matrix_DE[11, ] = to_array_1d(byearxorigin2_raw_vector_2_39);
byearxorigin2_raw_parameters_matrix_NDE[11, ] = to_array_1d(byearxorigin2_raw_vector_2_39);
byearxorigin2_raw_parameters_matrix_DE[12, ] = to_array_1d(byearxorigin2_raw_vector_3_2);
byearxorigin2_raw_parameters_matrix_NDE[12, ] = to_array_1d(byearxorigin2_raw_vector_3_2);
byearxorigin2_raw_parameters_matrix_DE[13, ] = to_array_1d(byearxorigin2_raw_vector_3_4);
byearxorigin2_raw_parameters_matrix_NDE[13, ] = to_array_1d(byearxorigin2_raw_vector_3_4);
byearxorigin2_raw_parameters_matrix_DE[14, ] = to_array_1d(byearxorigin2_raw_vector_3_8);
byearxorigin2_raw_parameters_matrix_NDE[14, ] = to_array_1d(byearxorigin2_raw_vector_3_8);
byearxorigin2_raw_parameters_matrix_DE[15, ] = to_array_1d(byearxorigin2_raw_vector_3_18_DE);
byearxorigin2_raw_parameters_matrix_NDE[15, ] = to_array_1d(byearxorigin2_raw_vector_3_18_NDE);
byearxorigin2_raw_parameters_matrix_DE[16, ] = to_array_1d(byearxorigin2_raw_vector_3_20_DE);
byearxorigin2_raw_parameters_matrix_NDE[16, ] = to_array_1d(byearxorigin2_raw_vector_3_20_NDE);
byearxorigin2_raw_parameters_matrix_DE[17, ] = to_array_1d(byearxorigin2_raw_vector_4_3);
byearxorigin2_raw_parameters_matrix_NDE[17, ] = to_array_1d(byearxorigin2_raw_vector_4_3);
byearxorigin2_raw_parameters_matrix_DE[27, ] = to_array_1d(byearxorigin2_raw_vector_8_3);
byearxorigin2_raw_parameters_matrix_NDE[27, ] = to_array_1d(byearxorigin2_raw_vector_8_3);

byearxorigin3_raw_parameters_matrix_DE[1, ] = to_array_1d(byearxorigin3_raw_vector_1_2);
byearxorigin3_raw_parameters_matrix_NDE[1, ] = to_array_1d(byearxorigin3_raw_vector_1_2);
byearxorigin3_raw_parameters_matrix_DE[2, ] = to_array_1d(byearxorigin3_raw_vector_2_1);
byearxorigin3_raw_parameters_matrix_NDE[2, ] = to_array_1d(byearxorigin3_raw_vector_2_1);
byearxorigin3_raw_parameters_matrix_DE[3, ] = to_array_1d(byearxorigin3_raw_vector_2_3);
byearxorigin3_raw_parameters_matrix_NDE[3, ] = to_array_1d(byearxorigin3_raw_vector_2_3);
byearxorigin3_raw_parameters_matrix_DE[4, ] = to_array_1d(byearxorigin3_raw_vector_2_10_DE);
byearxorigin3_raw_parameters_matrix_NDE[4, ] = to_array_1d(byearxorigin3_raw_vector_2_10_NDE);
byearxorigin3_raw_parameters_matrix_NDE[5, ] = to_array_1d(byearxorigin3_raw_vector_2_10_NDE);
byearxorigin3_raw_parameters_matrix_DE[6, ] = to_array_1d(byearxorigin3_raw_vector_2_12_DE);
byearxorigin3_raw_parameters_matrix_NDE[6, ] = to_array_1d(byearxorigin3_raw_vector_2_12_NDE);
byearxorigin3_raw_parameters_matrix_DE[7, ] = to_array_1d(byearxorigin3_raw_vector_2_14_DE);
byearxorigin3_raw_parameters_matrix_NDE[7, ] = to_array_1d(byearxorigin3_raw_vector_2_14_NDE);
byearxorigin3_raw_parameters_matrix_NDE[8, ] = to_array_1d(byearxorigin3_raw_vector_2_14_NDE);
byearxorigin3_raw_parameters_matrix_DE[9, ] = to_array_1d(byearxorigin3_raw_vector_2_16_DE);
byearxorigin3_raw_parameters_matrix_NDE[9, ] = to_array_1d(byearxorigin3_raw_vector_2_16_NDE);
byearxorigin3_raw_parameters_matrix_NDE[10, ] = to_array_1d(byearxorigin3_raw_vector_2_16_NDE);
byearxorigin3_raw_parameters_matrix_DE[11, ] = to_array_1d(byearxorigin3_raw_vector_2_39);
byearxorigin3_raw_parameters_matrix_NDE[11, ] = to_array_1d(byearxorigin3_raw_vector_2_39);
byearxorigin3_raw_parameters_matrix_DE[12, ] = to_array_1d(byearxorigin3_raw_vector_3_2);
byearxorigin3_raw_parameters_matrix_NDE[12, ] = to_array_1d(byearxorigin3_raw_vector_3_2);
byearxorigin3_raw_parameters_matrix_DE[13, ] = to_array_1d(byearxorigin3_raw_vector_3_4);
byearxorigin3_raw_parameters_matrix_NDE[13, ] = to_array_1d(byearxorigin3_raw_vector_3_4);
byearxorigin3_raw_parameters_matrix_DE[14, ] = to_array_1d(byearxorigin3_raw_vector_3_8);
byearxorigin3_raw_parameters_matrix_NDE[14, ] = to_array_1d(byearxorigin3_raw_vector_3_8);
byearxorigin3_raw_parameters_matrix_DE[15, ] = to_array_1d(byearxorigin3_raw_vector_3_18_DE);
byearxorigin3_raw_parameters_matrix_NDE[15, ] = to_array_1d(byearxorigin3_raw_vector_3_18_NDE);
byearxorigin3_raw_parameters_matrix_DE[16, ] = to_array_1d(byearxorigin3_raw_vector_3_20_DE);
byearxorigin3_raw_parameters_matrix_NDE[16, ] = to_array_1d(byearxorigin3_raw_vector_3_20_NDE);
byearxorigin3_raw_parameters_matrix_DE[17, ] = to_array_1d(byearxorigin3_raw_vector_4_3);
byearxorigin3_raw_parameters_matrix_NDE[17, ] = to_array_1d(byearxorigin3_raw_vector_4_3);
byearxorigin3_raw_parameters_matrix_DE[27, ] = to_array_1d(byearxorigin3_raw_vector_8_3);
byearxorigin3_raw_parameters_matrix_NDE[27, ] = to_array_1d(byearxorigin3_raw_vector_8_3);

byearxorigin4_raw_parameters_matrix_DE[1, ] = to_array_1d(byearxorigin4_raw_vector_1_2);
byearxorigin4_raw_parameters_matrix_NDE[1, ] = to_array_1d(byearxorigin4_raw_vector_1_2);
byearxorigin4_raw_parameters_matrix_DE[2, ] = to_array_1d(byearxorigin4_raw_vector_2_1);
byearxorigin4_raw_parameters_matrix_NDE[2, ] = to_array_1d(byearxorigin4_raw_vector_2_1);
byearxorigin4_raw_parameters_matrix_DE[3, ] = to_array_1d(byearxorigin4_raw_vector_2_3);
byearxorigin4_raw_parameters_matrix_NDE[3, ] = to_array_1d(byearxorigin4_raw_vector_2_3);
byearxorigin4_raw_parameters_matrix_DE[4, ] = to_array_1d(byearxorigin4_raw_vector_2_10_DE);
byearxorigin4_raw_parameters_matrix_NDE[4, ] = to_array_1d(byearxorigin4_raw_vector_2_10_NDE);
byearxorigin4_raw_parameters_matrix_NDE[5, ] = to_array_1d(byearxorigin4_raw_vector_2_10_NDE);
byearxorigin4_raw_parameters_matrix_DE[6, ] = to_array_1d(byearxorigin4_raw_vector_2_12_DE);
byearxorigin4_raw_parameters_matrix_NDE[6, ] = to_array_1d(byearxorigin4_raw_vector_2_12_NDE);
byearxorigin4_raw_parameters_matrix_DE[7, ] = to_array_1d(byearxorigin4_raw_vector_2_14_DE);
byearxorigin4_raw_parameters_matrix_NDE[7, ] = to_array_1d(byearxorigin4_raw_vector_2_14_NDE);
byearxorigin4_raw_parameters_matrix_NDE[8, ] = to_array_1d(byearxorigin4_raw_vector_2_14_NDE);
byearxorigin4_raw_parameters_matrix_DE[9, ] = to_array_1d(byearxorigin4_raw_vector_2_16_DE);
byearxorigin4_raw_parameters_matrix_NDE[9, ] = to_array_1d(byearxorigin4_raw_vector_2_16_NDE);
byearxorigin4_raw_parameters_matrix_NDE[10, ] = to_array_1d(byearxorigin4_raw_vector_2_16_NDE);
byearxorigin4_raw_parameters_matrix_DE[11, ] = to_array_1d(byearxorigin4_raw_vector_2_39);
byearxorigin4_raw_parameters_matrix_NDE[11, ] = to_array_1d(byearxorigin4_raw_vector_2_39);
byearxorigin4_raw_parameters_matrix_DE[12, ] = to_array_1d(byearxorigin4_raw_vector_3_2);
byearxorigin4_raw_parameters_matrix_NDE[12, ] = to_array_1d(byearxorigin4_raw_vector_3_2);
byearxorigin4_raw_parameters_matrix_DE[13, ] = to_array_1d(byearxorigin4_raw_vector_3_4);
byearxorigin4_raw_parameters_matrix_NDE[13, ] = to_array_1d(byearxorigin4_raw_vector_3_4);
byearxorigin4_raw_parameters_matrix_DE[14, ] = to_array_1d(byearxorigin4_raw_vector_3_8);
byearxorigin4_raw_parameters_matrix_NDE[14, ] = to_array_1d(byearxorigin4_raw_vector_3_8);
byearxorigin4_raw_parameters_matrix_DE[15, ] = to_array_1d(byearxorigin4_raw_vector_3_18_DE);
byearxorigin4_raw_parameters_matrix_NDE[15, ] = to_array_1d(byearxorigin4_raw_vector_3_18_NDE);
byearxorigin4_raw_parameters_matrix_DE[16, ] = to_array_1d(byearxorigin4_raw_vector_3_20_DE);
byearxorigin4_raw_parameters_matrix_NDE[16, ] = to_array_1d(byearxorigin4_raw_vector_3_20_NDE);
byearxorigin4_raw_parameters_matrix_DE[17, ] = to_array_1d(byearxorigin4_raw_vector_4_3);
byearxorigin4_raw_parameters_matrix_NDE[17, ] = to_array_1d(byearxorigin4_raw_vector_4_3);
byearxorigin4_raw_parameters_matrix_DE[27, ] = to_array_1d(byearxorigin4_raw_vector_8_3);
byearxorigin4_raw_parameters_matrix_NDE[27, ] = to_array_1d(byearxorigin4_raw_vector_8_3);

byearxorigin5_raw_parameters_matrix_DE[1, ] = to_array_1d(byearxorigin5_raw_vector_1_2);
byearxorigin5_raw_parameters_matrix_NDE[1, ] = to_array_1d(byearxorigin5_raw_vector_1_2);
byearxorigin5_raw_parameters_matrix_DE[2, ] = to_array_1d(byearxorigin5_raw_vector_2_1);
byearxorigin5_raw_parameters_matrix_NDE[2, ] = to_array_1d(byearxorigin5_raw_vector_2_1);
byearxorigin5_raw_parameters_matrix_DE[3, ] = to_array_1d(byearxorigin5_raw_vector_2_3);
byearxorigin5_raw_parameters_matrix_NDE[3, ] = to_array_1d(byearxorigin5_raw_vector_2_3);
byearxorigin5_raw_parameters_matrix_DE[4, ] = to_array_1d(byearxorigin5_raw_vector_2_10_DE);
byearxorigin5_raw_parameters_matrix_NDE[4, ] = to_array_1d(byearxorigin5_raw_vector_2_10_NDE);
byearxorigin5_raw_parameters_matrix_NDE[5, ] = to_array_1d(byearxorigin5_raw_vector_2_10_NDE);
byearxorigin5_raw_parameters_matrix_DE[6, ] = to_array_1d(byearxorigin5_raw_vector_2_12_DE);
byearxorigin5_raw_parameters_matrix_NDE[6, ] = to_array_1d(byearxorigin5_raw_vector_2_12_NDE);
byearxorigin5_raw_parameters_matrix_DE[7, ] = to_array_1d(byearxorigin5_raw_vector_2_14_DE);
byearxorigin5_raw_parameters_matrix_NDE[7, ] = to_array_1d(byearxorigin5_raw_vector_2_14_NDE);
byearxorigin5_raw_parameters_matrix_NDE[8, ] = to_array_1d(byearxorigin5_raw_vector_2_14_NDE);
byearxorigin5_raw_parameters_matrix_DE[9, ] = to_array_1d(byearxorigin5_raw_vector_2_16_DE);
byearxorigin5_raw_parameters_matrix_NDE[9, ] = to_array_1d(byearxorigin5_raw_vector_2_16_NDE);
byearxorigin5_raw_parameters_matrix_NDE[10, ] = to_array_1d(byearxorigin5_raw_vector_2_16_NDE);
byearxorigin5_raw_parameters_matrix_DE[11, ] = to_array_1d(byearxorigin5_raw_vector_2_39);
byearxorigin5_raw_parameters_matrix_NDE[11, ] = to_array_1d(byearxorigin5_raw_vector_2_39);
byearxorigin5_raw_parameters_matrix_DE[12, ] = to_array_1d(byearxorigin5_raw_vector_3_2);
byearxorigin5_raw_parameters_matrix_NDE[12, ] = to_array_1d(byearxorigin5_raw_vector_3_2);
byearxorigin5_raw_parameters_matrix_DE[13, ] = to_array_1d(byearxorigin5_raw_vector_3_4);
byearxorigin5_raw_parameters_matrix_NDE[13, ] = to_array_1d(byearxorigin5_raw_vector_3_4);
byearxorigin5_raw_parameters_matrix_DE[14, ] = to_array_1d(byearxorigin5_raw_vector_3_8);
byearxorigin5_raw_parameters_matrix_NDE[14, ] = to_array_1d(byearxorigin5_raw_vector_3_8);
byearxorigin5_raw_parameters_matrix_DE[15, ] = to_array_1d(byearxorigin5_raw_vector_3_18_DE);
byearxorigin5_raw_parameters_matrix_NDE[15, ] = to_array_1d(byearxorigin5_raw_vector_3_18_NDE);
byearxorigin5_raw_parameters_matrix_DE[16, ] = to_array_1d(byearxorigin5_raw_vector_3_20_DE);
byearxorigin5_raw_parameters_matrix_NDE[16, ] = to_array_1d(byearxorigin5_raw_vector_3_20_NDE);
byearxorigin5_raw_parameters_matrix_DE[17, ] = to_array_1d(byearxorigin5_raw_vector_4_3);
byearxorigin5_raw_parameters_matrix_NDE[17, ] = to_array_1d(byearxorigin5_raw_vector_4_3);
byearxorigin5_raw_parameters_matrix_DE[27, ] = to_array_1d(byearxorigin5_raw_vector_8_3);
byearxorigin5_raw_parameters_matrix_NDE[27, ] = to_array_1d(byearxorigin5_raw_vector_8_3);

byearxorigin6_raw_parameters_matrix_DE[1, ] = to_array_1d(byearxorigin6_raw_vector_1_2);
byearxorigin6_raw_parameters_matrix_NDE[1, ] = to_array_1d(byearxorigin6_raw_vector_1_2);
byearxorigin6_raw_parameters_matrix_DE[2, ] = to_array_1d(byearxorigin6_raw_vector_2_1);
byearxorigin6_raw_parameters_matrix_NDE[2, ] = to_array_1d(byearxorigin6_raw_vector_2_1);
byearxorigin6_raw_parameters_matrix_DE[3, ] = to_array_1d(byearxorigin6_raw_vector_2_3);
byearxorigin6_raw_parameters_matrix_NDE[3, ] = to_array_1d(byearxorigin6_raw_vector_2_3);
byearxorigin6_raw_parameters_matrix_DE[4, ] = to_array_1d(byearxorigin6_raw_vector_2_10_DE);
byearxorigin6_raw_parameters_matrix_NDE[4, ] = to_array_1d(byearxorigin6_raw_vector_2_10_NDE);
byearxorigin6_raw_parameters_matrix_NDE[5, ] = to_array_1d(byearxorigin6_raw_vector_2_10_NDE);
byearxorigin6_raw_parameters_matrix_DE[6, ] = to_array_1d(byearxorigin6_raw_vector_2_12_DE);
byearxorigin6_raw_parameters_matrix_NDE[6, ] = to_array_1d(byearxorigin6_raw_vector_2_12_NDE);
byearxorigin6_raw_parameters_matrix_DE[7, ] = to_array_1d(byearxorigin6_raw_vector_2_14_DE);
byearxorigin6_raw_parameters_matrix_NDE[7, ] = to_array_1d(byearxorigin6_raw_vector_2_14_NDE);
byearxorigin6_raw_parameters_matrix_NDE[8, ] = to_array_1d(byearxorigin6_raw_vector_2_14_NDE);
byearxorigin6_raw_parameters_matrix_DE[9, ] = to_array_1d(byearxorigin6_raw_vector_2_16_DE);
byearxorigin6_raw_parameters_matrix_NDE[9, ] = to_array_1d(byearxorigin6_raw_vector_2_16_NDE);
byearxorigin6_raw_parameters_matrix_NDE[10, ] = to_array_1d(byearxorigin6_raw_vector_2_16_NDE);
byearxorigin6_raw_parameters_matrix_DE[11, ] = to_array_1d(byearxorigin6_raw_vector_2_39);
byearxorigin6_raw_parameters_matrix_NDE[11, ] = to_array_1d(byearxorigin6_raw_vector_2_39);
byearxorigin6_raw_parameters_matrix_DE[12, ] = to_array_1d(byearxorigin6_raw_vector_3_2);
byearxorigin6_raw_parameters_matrix_NDE[12, ] = to_array_1d(byearxorigin6_raw_vector_3_2);
byearxorigin6_raw_parameters_matrix_DE[13, ] = to_array_1d(byearxorigin6_raw_vector_3_4);
byearxorigin6_raw_parameters_matrix_NDE[13, ] = to_array_1d(byearxorigin6_raw_vector_3_4);
byearxorigin6_raw_parameters_matrix_DE[14, ] = to_array_1d(byearxorigin6_raw_vector_3_8);
byearxorigin6_raw_parameters_matrix_NDE[14, ] = to_array_1d(byearxorigin6_raw_vector_3_8);
byearxorigin6_raw_parameters_matrix_DE[15, ] = to_array_1d(byearxorigin6_raw_vector_3_18_DE);
byearxorigin6_raw_parameters_matrix_NDE[15, ] = to_array_1d(byearxorigin6_raw_vector_3_18_NDE);
byearxorigin6_raw_parameters_matrix_DE[16, ] = to_array_1d(byearxorigin6_raw_vector_3_20_DE);
byearxorigin6_raw_parameters_matrix_NDE[16, ] = to_array_1d(byearxorigin6_raw_vector_3_20_NDE);
byearxorigin6_raw_parameters_matrix_DE[17, ] = to_array_1d(byearxorigin6_raw_vector_4_3);
byearxorigin6_raw_parameters_matrix_NDE[17, ] = to_array_1d(byearxorigin6_raw_vector_4_3);
byearxorigin6_raw_parameters_matrix_DE[27, ] = to_array_1d(byearxorigin6_raw_vector_8_3);
byearxorigin6_raw_parameters_matrix_NDE[27, ] = to_array_1d(byearxorigin6_raw_vector_8_3);
  
// 2. now, take the byear_raw_vector parameters and multiply them by the sigma (scale) to store them)
// take the raw matrices and transform them to the actual matrices

byearxorigin1_actual_parameters_matrix_DE[1, ] = to_array_1d(byearxorigin1_raw_vector_1_2 * sigma_yearxorigin1_vector_1_2);
byearxorigin1_actual_parameters_matrix_NDE[1, ] = to_array_1d(byearxorigin1_raw_vector_1_2 * sigma_yearxorigin1_vector_1_2);
byearxorigin1_actual_parameters_matrix_DE[2, ] = to_array_1d(byearxorigin1_raw_vector_2_1 * sigma_yearxorigin1_vector_2_1);
byearxorigin1_actual_parameters_matrix_NDE[2, ] = to_array_1d(byearxorigin1_raw_vector_2_1 * sigma_yearxorigin1_vector_2_1);
byearxorigin1_actual_parameters_matrix_DE[3, ] = to_array_1d(byearxorigin1_raw_vector_2_3 * sigma_yearxorigin1_vector_2_3);
byearxorigin1_actual_parameters_matrix_NDE[3, ] = to_array_1d(byearxorigin1_raw_vector_2_3 * sigma_yearxorigin1_vector_2_3);
byearxorigin1_actual_parameters_matrix_DE[4, ] = to_array_1d(byearxorigin1_raw_vector_2_10_DE * sigma_yearxorigin1_vector_2_10_DE);
byearxorigin1_actual_parameters_matrix_NDE[4, ] = to_array_1d(byearxorigin1_raw_vector_2_10_NDE * sigma_yearxorigin1_vector_2_10_NDE);
byearxorigin1_actual_parameters_matrix_NDE[5, ] = to_array_1d(byearxorigin1_raw_vector_2_10_NDE * sigma_yearxorigin1_vector_2_10_NDE);
byearxorigin1_actual_parameters_matrix_DE[6, ] = to_array_1d(byearxorigin1_raw_vector_2_12_DE * sigma_yearxorigin1_vector_2_12_DE);
byearxorigin1_actual_parameters_matrix_NDE[6, ] = to_array_1d(byearxorigin1_raw_vector_2_12_NDE * sigma_yearxorigin1_vector_2_12_NDE);
byearxorigin1_actual_parameters_matrix_DE[7, ] = to_array_1d(byearxorigin1_raw_vector_2_14_DE * sigma_yearxorigin1_vector_2_14_DE);
byearxorigin1_actual_parameters_matrix_NDE[7, ] = to_array_1d(byearxorigin1_raw_vector_2_14_NDE * sigma_yearxorigin1_vector_2_14_NDE);
byearxorigin1_actual_parameters_matrix_NDE[8, ] = to_array_1d(byearxorigin1_raw_vector_2_14_NDE * sigma_yearxorigin1_vector_2_14_NDE);
byearxorigin1_actual_parameters_matrix_DE[9, ] = to_array_1d(byearxorigin1_raw_vector_2_16_DE * sigma_yearxorigin1_vector_2_16_DE);
byearxorigin1_actual_parameters_matrix_NDE[9, ] = to_array_1d(byearxorigin1_raw_vector_2_16_NDE * sigma_yearxorigin1_vector_2_16_NDE);
byearxorigin1_actual_parameters_matrix_NDE[10, ] = to_array_1d(byearxorigin1_raw_vector_2_16_NDE * sigma_yearxorigin1_vector_2_16_NDE);
byearxorigin1_actual_parameters_matrix_DE[11, ] = to_array_1d(byearxorigin1_raw_vector_2_39 * sigma_yearxorigin1_vector_2_39);
byearxorigin1_actual_parameters_matrix_NDE[11, ] = to_array_1d(byearxorigin1_raw_vector_2_39 * sigma_yearxorigin1_vector_2_39);
byearxorigin1_actual_parameters_matrix_DE[12, ] = to_array_1d(byearxorigin1_raw_vector_3_2 * sigma_yearxorigin1_vector_3_2);
byearxorigin1_actual_parameters_matrix_NDE[12, ] = to_array_1d(byearxorigin1_raw_vector_3_2 * sigma_yearxorigin1_vector_3_2);
byearxorigin1_actual_parameters_matrix_DE[13, ] = to_array_1d(byearxorigin1_raw_vector_3_4 * sigma_yearxorigin1_vector_3_4);
byearxorigin1_actual_parameters_matrix_NDE[13, ] = to_array_1d(byearxorigin1_raw_vector_3_4 * sigma_yearxorigin1_vector_3_4);
byearxorigin1_actual_parameters_matrix_DE[14, ] = to_array_1d(byearxorigin1_raw_vector_3_8 * sigma_yearxorigin1_vector_3_8);
byearxorigin1_actual_parameters_matrix_NDE[14, ] = to_array_1d(byearxorigin1_raw_vector_3_8 * sigma_yearxorigin1_vector_3_8);
byearxorigin1_actual_parameters_matrix_DE[15, ] = to_array_1d(byearxorigin1_raw_vector_3_18_DE * sigma_yearxorigin1_vector_3_18_DE);
byearxorigin1_actual_parameters_matrix_NDE[15, ] = to_array_1d(byearxorigin1_raw_vector_3_18_NDE * sigma_yearxorigin1_vector_3_18_NDE);
byearxorigin1_actual_parameters_matrix_DE[16, ] = to_array_1d(byearxorigin1_raw_vector_3_20_DE * sigma_yearxorigin1_vector_3_20_DE);
byearxorigin1_actual_parameters_matrix_NDE[16, ] = to_array_1d(byearxorigin1_raw_vector_3_20_NDE * sigma_yearxorigin1_vector_3_20_NDE);
byearxorigin1_actual_parameters_matrix_DE[17, ] = to_array_1d(byearxorigin1_raw_vector_4_3 * sigma_yearxorigin1_vector_4_3);
byearxorigin1_actual_parameters_matrix_NDE[17, ] = to_array_1d(byearxorigin1_raw_vector_4_3 * sigma_yearxorigin1_vector_4_3);
byearxorigin1_actual_parameters_matrix_DE[27, ] = to_array_1d(byearxorigin1_raw_vector_8_3 * sigma_yearxorigin1_vector_8_3);
byearxorigin1_actual_parameters_matrix_NDE[27, ] = to_array_1d(byearxorigin1_raw_vector_8_3 * sigma_yearxorigin1_vector_8_3);

byearxorigin2_actual_parameters_matrix_DE[1, ] = to_array_1d(byearxorigin2_raw_vector_1_2 * sigma_yearxorigin2_vector_1_2);
byearxorigin2_actual_parameters_matrix_NDE[1, ] = to_array_1d(byearxorigin2_raw_vector_1_2 * sigma_yearxorigin2_vector_1_2);
byearxorigin2_actual_parameters_matrix_DE[2, ] = to_array_1d(byearxorigin2_raw_vector_2_1 * sigma_yearxorigin2_vector_2_1);
byearxorigin2_actual_parameters_matrix_NDE[2, ] = to_array_1d(byearxorigin2_raw_vector_2_1 * sigma_yearxorigin2_vector_2_1);
byearxorigin2_actual_parameters_matrix_DE[3, ] = to_array_1d(byearxorigin2_raw_vector_2_3 * sigma_yearxorigin2_vector_2_3);
byearxorigin2_actual_parameters_matrix_NDE[3, ] = to_array_1d(byearxorigin2_raw_vector_2_3 * sigma_yearxorigin2_vector_2_3);
byearxorigin2_actual_parameters_matrix_DE[4, ] = to_array_1d(byearxorigin2_raw_vector_2_10_DE * sigma_yearxorigin2_vector_2_10_DE);
byearxorigin2_actual_parameters_matrix_NDE[4, ] = to_array_1d(byearxorigin2_raw_vector_2_10_NDE * sigma_yearxorigin2_vector_2_10_NDE);
byearxorigin2_actual_parameters_matrix_NDE[5, ] = to_array_1d(byearxorigin2_raw_vector_2_10_NDE * sigma_yearxorigin2_vector_2_10_NDE);
byearxorigin2_actual_parameters_matrix_DE[6, ] = to_array_1d(byearxorigin2_raw_vector_2_12_DE * sigma_yearxorigin2_vector_2_12_DE);
byearxorigin2_actual_parameters_matrix_NDE[6, ] = to_array_1d(byearxorigin2_raw_vector_2_12_NDE * sigma_yearxorigin2_vector_2_12_NDE);
byearxorigin2_actual_parameters_matrix_DE[7, ] = to_array_1d(byearxorigin2_raw_vector_2_14_DE * sigma_yearxorigin2_vector_2_14_DE);
byearxorigin2_actual_parameters_matrix_NDE[7, ] = to_array_1d(byearxorigin2_raw_vector_2_14_NDE * sigma_yearxorigin2_vector_2_14_NDE);
byearxorigin2_actual_parameters_matrix_NDE[8, ] = to_array_1d(byearxorigin2_raw_vector_2_14_NDE * sigma_yearxorigin2_vector_2_14_NDE);
byearxorigin2_actual_parameters_matrix_DE[9, ] = to_array_1d(byearxorigin2_raw_vector_2_16_DE * sigma_yearxorigin2_vector_2_16_DE);
byearxorigin2_actual_parameters_matrix_NDE[9, ] = to_array_1d(byearxorigin2_raw_vector_2_16_NDE * sigma_yearxorigin2_vector_2_16_NDE);
byearxorigin2_actual_parameters_matrix_NDE[10, ] = to_array_1d(byearxorigin2_raw_vector_2_16_NDE * sigma_yearxorigin2_vector_2_16_NDE);
byearxorigin2_actual_parameters_matrix_DE[11, ] = to_array_1d(byearxorigin2_raw_vector_2_39 * sigma_yearxorigin2_vector_2_39);
byearxorigin2_actual_parameters_matrix_NDE[11, ] = to_array_1d(byearxorigin2_raw_vector_2_39 * sigma_yearxorigin2_vector_2_39);
byearxorigin2_actual_parameters_matrix_DE[12, ] = to_array_1d(byearxorigin2_raw_vector_3_2 * sigma_yearxorigin2_vector_3_2);
byearxorigin2_actual_parameters_matrix_NDE[12, ] = to_array_1d(byearxorigin2_raw_vector_3_2 * sigma_yearxorigin2_vector_3_2);
byearxorigin2_actual_parameters_matrix_DE[13, ] = to_array_1d(byearxorigin2_raw_vector_3_4 * sigma_yearxorigin2_vector_3_4);
byearxorigin2_actual_parameters_matrix_NDE[13, ] = to_array_1d(byearxorigin2_raw_vector_3_4 * sigma_yearxorigin2_vector_3_4);
byearxorigin2_actual_parameters_matrix_DE[14, ] = to_array_1d(byearxorigin2_raw_vector_3_8 * sigma_yearxorigin2_vector_3_8);
byearxorigin2_actual_parameters_matrix_NDE[14, ] = to_array_1d(byearxorigin2_raw_vector_3_8 * sigma_yearxorigin2_vector_3_8);
byearxorigin2_actual_parameters_matrix_DE[15, ] = to_array_1d(byearxorigin2_raw_vector_3_18_DE * sigma_yearxorigin2_vector_3_18_DE);
byearxorigin2_actual_parameters_matrix_NDE[15, ] = to_array_1d(byearxorigin2_raw_vector_3_18_NDE * sigma_yearxorigin2_vector_3_18_NDE);
byearxorigin2_actual_parameters_matrix_DE[16, ] = to_array_1d(byearxorigin2_raw_vector_3_20_DE * sigma_yearxorigin2_vector_3_20_DE);
byearxorigin2_actual_parameters_matrix_NDE[16, ] = to_array_1d(byearxorigin2_raw_vector_3_20_NDE * sigma_yearxorigin2_vector_3_20_NDE);
byearxorigin2_actual_parameters_matrix_DE[17, ] = to_array_1d(byearxorigin2_raw_vector_4_3 * sigma_yearxorigin2_vector_4_3);
byearxorigin2_actual_parameters_matrix_NDE[17, ] = to_array_1d(byearxorigin2_raw_vector_4_3 * sigma_yearxorigin2_vector_4_3);
byearxorigin2_actual_parameters_matrix_DE[27, ] = to_array_1d(byearxorigin2_raw_vector_8_3 * sigma_yearxorigin2_vector_8_3);
byearxorigin2_actual_parameters_matrix_NDE[27, ] = to_array_1d(byearxorigin2_raw_vector_8_3 * sigma_yearxorigin2_vector_8_3);

byearxorigin3_actual_parameters_matrix_DE[1, ] = to_array_1d(byearxorigin3_raw_vector_1_2 * sigma_yearxorigin3_vector_1_2);
byearxorigin3_actual_parameters_matrix_NDE[1, ] = to_array_1d(byearxorigin3_raw_vector_1_2 * sigma_yearxorigin3_vector_1_2);
byearxorigin3_actual_parameters_matrix_DE[2, ] = to_array_1d(byearxorigin3_raw_vector_2_1 * sigma_yearxorigin3_vector_2_1);
byearxorigin3_actual_parameters_matrix_NDE[2, ] = to_array_1d(byearxorigin3_raw_vector_2_1 * sigma_yearxorigin3_vector_2_1);
byearxorigin3_actual_parameters_matrix_DE[3, ] = to_array_1d(byearxorigin3_raw_vector_2_3 * sigma_yearxorigin3_vector_2_3);
byearxorigin3_actual_parameters_matrix_NDE[3, ] = to_array_1d(byearxorigin3_raw_vector_2_3 * sigma_yearxorigin3_vector_2_3);
byearxorigin3_actual_parameters_matrix_DE[4, ] = to_array_1d(byearxorigin3_raw_vector_2_10_DE * sigma_yearxorigin3_vector_2_10_DE);
byearxorigin3_actual_parameters_matrix_NDE[4, ] = to_array_1d(byearxorigin3_raw_vector_2_10_NDE * sigma_yearxorigin3_vector_2_10_NDE);
byearxorigin3_actual_parameters_matrix_NDE[5, ] = to_array_1d(byearxorigin3_raw_vector_2_10_NDE * sigma_yearxorigin3_vector_2_10_NDE);
byearxorigin3_actual_parameters_matrix_DE[6, ] = to_array_1d(byearxorigin3_raw_vector_2_12_DE * sigma_yearxorigin3_vector_2_12_DE);
byearxorigin3_actual_parameters_matrix_NDE[6, ] = to_array_1d(byearxorigin3_raw_vector_2_12_NDE * sigma_yearxorigin3_vector_2_12_NDE);
byearxorigin3_actual_parameters_matrix_DE[7, ] = to_array_1d(byearxorigin3_raw_vector_2_14_DE * sigma_yearxorigin3_vector_2_14_DE);
byearxorigin3_actual_parameters_matrix_NDE[7, ] = to_array_1d(byearxorigin3_raw_vector_2_14_NDE * sigma_yearxorigin3_vector_2_14_NDE);
byearxorigin3_actual_parameters_matrix_NDE[8, ] = to_array_1d(byearxorigin3_raw_vector_2_14_NDE * sigma_yearxorigin3_vector_2_14_NDE);
byearxorigin3_actual_parameters_matrix_DE[9, ] = to_array_1d(byearxorigin3_raw_vector_2_16_DE * sigma_yearxorigin3_vector_2_16_DE);
byearxorigin3_actual_parameters_matrix_NDE[9, ] = to_array_1d(byearxorigin3_raw_vector_2_16_NDE * sigma_yearxorigin3_vector_2_16_NDE);
byearxorigin3_actual_parameters_matrix_NDE[10, ] = to_array_1d(byearxorigin3_raw_vector_2_16_NDE * sigma_yearxorigin3_vector_2_16_NDE);
byearxorigin3_actual_parameters_matrix_DE[11, ] = to_array_1d(byearxorigin3_raw_vector_2_39 * sigma_yearxorigin3_vector_2_39);
byearxorigin3_actual_parameters_matrix_NDE[11, ] = to_array_1d(byearxorigin3_raw_vector_2_39 * sigma_yearxorigin3_vector_2_39);
byearxorigin3_actual_parameters_matrix_DE[12, ] = to_array_1d(byearxorigin3_raw_vector_3_2 * sigma_yearxorigin3_vector_3_2);
byearxorigin3_actual_parameters_matrix_NDE[12, ] = to_array_1d(byearxorigin3_raw_vector_3_2 * sigma_yearxorigin3_vector_3_2);
byearxorigin3_actual_parameters_matrix_DE[13, ] = to_array_1d(byearxorigin3_raw_vector_3_4 * sigma_yearxorigin3_vector_3_4);
byearxorigin3_actual_parameters_matrix_NDE[13, ] = to_array_1d(byearxorigin3_raw_vector_3_4 * sigma_yearxorigin3_vector_3_4);
byearxorigin3_actual_parameters_matrix_DE[14, ] = to_array_1d(byearxorigin3_raw_vector_3_8 * sigma_yearxorigin3_vector_3_8);
byearxorigin3_actual_parameters_matrix_NDE[14, ] = to_array_1d(byearxorigin3_raw_vector_3_8 * sigma_yearxorigin3_vector_3_8);
byearxorigin3_actual_parameters_matrix_DE[15, ] = to_array_1d(byearxorigin3_raw_vector_3_18_DE * sigma_yearxorigin3_vector_3_18_DE);
byearxorigin3_actual_parameters_matrix_NDE[15, ] = to_array_1d(byearxorigin3_raw_vector_3_18_NDE * sigma_yearxorigin3_vector_3_18_NDE);
byearxorigin3_actual_parameters_matrix_DE[16, ] = to_array_1d(byearxorigin3_raw_vector_3_20_DE * sigma_yearxorigin3_vector_3_20_DE);
byearxorigin3_actual_parameters_matrix_NDE[16, ] = to_array_1d(byearxorigin3_raw_vector_3_20_NDE * sigma_yearxorigin3_vector_3_20_NDE);
byearxorigin3_actual_parameters_matrix_DE[17, ] = to_array_1d(byearxorigin3_raw_vector_4_3 * sigma_yearxorigin3_vector_4_3);
byearxorigin3_actual_parameters_matrix_NDE[17, ] = to_array_1d(byearxorigin3_raw_vector_4_3 * sigma_yearxorigin3_vector_4_3);
byearxorigin3_actual_parameters_matrix_DE[27, ] = to_array_1d(byearxorigin3_raw_vector_8_3 * sigma_yearxorigin3_vector_8_3);
byearxorigin3_actual_parameters_matrix_NDE[27, ] = to_array_1d(byearxorigin3_raw_vector_8_3 * sigma_yearxorigin3_vector_8_3);

byearxorigin4_actual_parameters_matrix_DE[1, ] = to_array_1d(byearxorigin4_raw_vector_1_2 * sigma_yearxorigin4_vector_1_2);
byearxorigin4_actual_parameters_matrix_NDE[1, ] = to_array_1d(byearxorigin4_raw_vector_1_2 * sigma_yearxorigin4_vector_1_2);
byearxorigin4_actual_parameters_matrix_DE[2, ] = to_array_1d(byearxorigin4_raw_vector_2_1 * sigma_yearxorigin4_vector_2_1);
byearxorigin4_actual_parameters_matrix_NDE[2, ] = to_array_1d(byearxorigin4_raw_vector_2_1 * sigma_yearxorigin4_vector_2_1);
byearxorigin4_actual_parameters_matrix_DE[3, ] = to_array_1d(byearxorigin4_raw_vector_2_3 * sigma_yearxorigin4_vector_2_3);
byearxorigin4_actual_parameters_matrix_NDE[3, ] = to_array_1d(byearxorigin4_raw_vector_2_3 * sigma_yearxorigin4_vector_2_3);
byearxorigin4_actual_parameters_matrix_DE[4, ] = to_array_1d(byearxorigin4_raw_vector_2_10_DE * sigma_yearxorigin4_vector_2_10_DE);
byearxorigin4_actual_parameters_matrix_NDE[4, ] = to_array_1d(byearxorigin4_raw_vector_2_10_NDE * sigma_yearxorigin4_vector_2_10_NDE);
byearxorigin4_actual_parameters_matrix_NDE[5, ] = to_array_1d(byearxorigin4_raw_vector_2_10_NDE * sigma_yearxorigin4_vector_2_10_NDE);
byearxorigin4_actual_parameters_matrix_DE[6, ] = to_array_1d(byearxorigin4_raw_vector_2_12_DE * sigma_yearxorigin4_vector_2_12_DE);
byearxorigin4_actual_parameters_matrix_NDE[6, ] = to_array_1d(byearxorigin4_raw_vector_2_12_NDE * sigma_yearxorigin4_vector_2_12_NDE);
byearxorigin4_actual_parameters_matrix_DE[7, ] = to_array_1d(byearxorigin4_raw_vector_2_14_DE * sigma_yearxorigin4_vector_2_14_DE);
byearxorigin4_actual_parameters_matrix_NDE[7, ] = to_array_1d(byearxorigin4_raw_vector_2_14_NDE * sigma_yearxorigin4_vector_2_14_NDE);
byearxorigin4_actual_parameters_matrix_NDE[8, ] = to_array_1d(byearxorigin4_raw_vector_2_14_NDE * sigma_yearxorigin4_vector_2_14_NDE);
byearxorigin4_actual_parameters_matrix_DE[9, ] = to_array_1d(byearxorigin4_raw_vector_2_16_DE * sigma_yearxorigin4_vector_2_16_DE);
byearxorigin4_actual_parameters_matrix_NDE[9, ] = to_array_1d(byearxorigin4_raw_vector_2_16_NDE * sigma_yearxorigin4_vector_2_16_NDE);
byearxorigin4_actual_parameters_matrix_NDE[10, ] = to_array_1d(byearxorigin4_raw_vector_2_16_NDE * sigma_yearxorigin4_vector_2_16_NDE);
byearxorigin4_actual_parameters_matrix_DE[11, ] = to_array_1d(byearxorigin4_raw_vector_2_39 * sigma_yearxorigin4_vector_2_39);
byearxorigin4_actual_parameters_matrix_NDE[11, ] = to_array_1d(byearxorigin4_raw_vector_2_39 * sigma_yearxorigin4_vector_2_39);
byearxorigin4_actual_parameters_matrix_DE[12, ] = to_array_1d(byearxorigin4_raw_vector_3_2 * sigma_yearxorigin4_vector_3_2);
byearxorigin4_actual_parameters_matrix_NDE[12, ] = to_array_1d(byearxorigin4_raw_vector_3_2 * sigma_yearxorigin4_vector_3_2);
byearxorigin4_actual_parameters_matrix_DE[13, ] = to_array_1d(byearxorigin4_raw_vector_3_4 * sigma_yearxorigin4_vector_3_4);
byearxorigin4_actual_parameters_matrix_NDE[13, ] = to_array_1d(byearxorigin4_raw_vector_3_4 * sigma_yearxorigin4_vector_3_4);
byearxorigin4_actual_parameters_matrix_DE[14, ] = to_array_1d(byearxorigin4_raw_vector_3_8 * sigma_yearxorigin4_vector_3_8);
byearxorigin4_actual_parameters_matrix_NDE[14, ] = to_array_1d(byearxorigin4_raw_vector_3_8 * sigma_yearxorigin4_vector_3_8);
byearxorigin4_actual_parameters_matrix_DE[15, ] = to_array_1d(byearxorigin4_raw_vector_3_18_DE * sigma_yearxorigin4_vector_3_18_DE);
byearxorigin4_actual_parameters_matrix_NDE[15, ] = to_array_1d(byearxorigin4_raw_vector_3_18_NDE * sigma_yearxorigin4_vector_3_18_NDE);
byearxorigin4_actual_parameters_matrix_DE[16, ] = to_array_1d(byearxorigin4_raw_vector_3_20_DE * sigma_yearxorigin4_vector_3_20_DE);
byearxorigin4_actual_parameters_matrix_NDE[16, ] = to_array_1d(byearxorigin4_raw_vector_3_20_NDE * sigma_yearxorigin4_vector_3_20_NDE);
byearxorigin4_actual_parameters_matrix_DE[17, ] = to_array_1d(byearxorigin4_raw_vector_4_3 * sigma_yearxorigin4_vector_4_3);
byearxorigin4_actual_parameters_matrix_NDE[17, ] = to_array_1d(byearxorigin4_raw_vector_4_3 * sigma_yearxorigin4_vector_4_3);
byearxorigin4_actual_parameters_matrix_DE[27, ] = to_array_1d(byearxorigin4_raw_vector_8_3 * sigma_yearxorigin4_vector_8_3);
byearxorigin4_actual_parameters_matrix_NDE[27, ] = to_array_1d(byearxorigin4_raw_vector_8_3 * sigma_yearxorigin4_vector_8_3);

byearxorigin5_actual_parameters_matrix_DE[1, ] = to_array_1d(byearxorigin5_raw_vector_1_2 * sigma_yearxorigin5_vector_1_2);
byearxorigin5_actual_parameters_matrix_NDE[1, ] = to_array_1d(byearxorigin5_raw_vector_1_2 * sigma_yearxorigin5_vector_1_2);
byearxorigin5_actual_parameters_matrix_DE[2, ] = to_array_1d(byearxorigin5_raw_vector_2_1 * sigma_yearxorigin5_vector_2_1);
byearxorigin5_actual_parameters_matrix_NDE[2, ] = to_array_1d(byearxorigin5_raw_vector_2_1 * sigma_yearxorigin5_vector_2_1);
byearxorigin5_actual_parameters_matrix_DE[3, ] = to_array_1d(byearxorigin5_raw_vector_2_3 * sigma_yearxorigin5_vector_2_3);
byearxorigin5_actual_parameters_matrix_NDE[3, ] = to_array_1d(byearxorigin5_raw_vector_2_3 * sigma_yearxorigin5_vector_2_3);
byearxorigin5_actual_parameters_matrix_DE[4, ] = to_array_1d(byearxorigin5_raw_vector_2_10_DE * sigma_yearxorigin5_vector_2_10_DE);
byearxorigin5_actual_parameters_matrix_NDE[4, ] = to_array_1d(byearxorigin5_raw_vector_2_10_NDE * sigma_yearxorigin5_vector_2_10_NDE);
byearxorigin5_actual_parameters_matrix_NDE[5, ] = to_array_1d(byearxorigin5_raw_vector_2_10_NDE * sigma_yearxorigin5_vector_2_10_NDE);
byearxorigin5_actual_parameters_matrix_DE[6, ] = to_array_1d(byearxorigin5_raw_vector_2_12_DE * sigma_yearxorigin5_vector_2_12_DE);
byearxorigin5_actual_parameters_matrix_NDE[6, ] = to_array_1d(byearxorigin5_raw_vector_2_12_NDE * sigma_yearxorigin5_vector_2_12_NDE);
byearxorigin5_actual_parameters_matrix_DE[7, ] = to_array_1d(byearxorigin5_raw_vector_2_14_DE * sigma_yearxorigin5_vector_2_14_DE);
byearxorigin5_actual_parameters_matrix_NDE[7, ] = to_array_1d(byearxorigin5_raw_vector_2_14_NDE * sigma_yearxorigin5_vector_2_14_NDE);
byearxorigin5_actual_parameters_matrix_NDE[8, ] = to_array_1d(byearxorigin5_raw_vector_2_14_NDE * sigma_yearxorigin5_vector_2_14_NDE);
byearxorigin5_actual_parameters_matrix_DE[9, ] = to_array_1d(byearxorigin5_raw_vector_2_16_DE * sigma_yearxorigin5_vector_2_16_DE);
byearxorigin5_actual_parameters_matrix_NDE[9, ] = to_array_1d(byearxorigin5_raw_vector_2_16_NDE * sigma_yearxorigin5_vector_2_16_NDE);
byearxorigin5_actual_parameters_matrix_NDE[10, ] = to_array_1d(byearxorigin5_raw_vector_2_16_NDE * sigma_yearxorigin5_vector_2_16_NDE);
byearxorigin5_actual_parameters_matrix_DE[11, ] = to_array_1d(byearxorigin5_raw_vector_2_39 * sigma_yearxorigin5_vector_2_39);
byearxorigin5_actual_parameters_matrix_NDE[11, ] = to_array_1d(byearxorigin5_raw_vector_2_39 * sigma_yearxorigin5_vector_2_39);
byearxorigin5_actual_parameters_matrix_DE[12, ] = to_array_1d(byearxorigin5_raw_vector_3_2 * sigma_yearxorigin5_vector_3_2);
byearxorigin5_actual_parameters_matrix_NDE[12, ] = to_array_1d(byearxorigin5_raw_vector_3_2 * sigma_yearxorigin5_vector_3_2);
byearxorigin5_actual_parameters_matrix_DE[13, ] = to_array_1d(byearxorigin5_raw_vector_3_4 * sigma_yearxorigin5_vector_3_4);
byearxorigin5_actual_parameters_matrix_NDE[13, ] = to_array_1d(byearxorigin5_raw_vector_3_4 * sigma_yearxorigin5_vector_3_4);
byearxorigin5_actual_parameters_matrix_DE[14, ] = to_array_1d(byearxorigin5_raw_vector_3_8 * sigma_yearxorigin5_vector_3_8);
byearxorigin5_actual_parameters_matrix_NDE[14, ] = to_array_1d(byearxorigin5_raw_vector_3_8 * sigma_yearxorigin5_vector_3_8);
byearxorigin5_actual_parameters_matrix_DE[15, ] = to_array_1d(byearxorigin5_raw_vector_3_18_DE * sigma_yearxorigin5_vector_3_18_DE);
byearxorigin5_actual_parameters_matrix_NDE[15, ] = to_array_1d(byearxorigin5_raw_vector_3_18_NDE * sigma_yearxorigin5_vector_3_18_NDE);
byearxorigin5_actual_parameters_matrix_DE[16, ] = to_array_1d(byearxorigin5_raw_vector_3_20_DE * sigma_yearxorigin5_vector_3_20_DE);
byearxorigin5_actual_parameters_matrix_NDE[16, ] = to_array_1d(byearxorigin5_raw_vector_3_20_NDE * sigma_yearxorigin5_vector_3_20_NDE);
byearxorigin5_actual_parameters_matrix_DE[17, ] = to_array_1d(byearxorigin5_raw_vector_4_3 * sigma_yearxorigin5_vector_4_3);
byearxorigin5_actual_parameters_matrix_NDE[17, ] = to_array_1d(byearxorigin5_raw_vector_4_3 * sigma_yearxorigin5_vector_4_3);
byearxorigin5_actual_parameters_matrix_DE[27, ] = to_array_1d(byearxorigin5_raw_vector_8_3 * sigma_yearxorigin5_vector_8_3);
byearxorigin5_actual_parameters_matrix_NDE[27, ] = to_array_1d(byearxorigin5_raw_vector_8_3 * sigma_yearxorigin5_vector_8_3);

byearxorigin6_actual_parameters_matrix_DE[1, ] = to_array_1d(byearxorigin6_raw_vector_1_2 * sigma_yearxorigin6_vector_1_2);
byearxorigin6_actual_parameters_matrix_NDE[1, ] = to_array_1d(byearxorigin6_raw_vector_1_2 * sigma_yearxorigin6_vector_1_2);
byearxorigin6_actual_parameters_matrix_DE[2, ] = to_array_1d(byearxorigin6_raw_vector_2_1 * sigma_yearxorigin6_vector_2_1);
byearxorigin6_actual_parameters_matrix_NDE[2, ] = to_array_1d(byearxorigin6_raw_vector_2_1 * sigma_yearxorigin6_vector_2_1);
byearxorigin6_actual_parameters_matrix_DE[3, ] = to_array_1d(byearxorigin6_raw_vector_2_3 * sigma_yearxorigin6_vector_2_3);
byearxorigin6_actual_parameters_matrix_NDE[3, ] = to_array_1d(byearxorigin6_raw_vector_2_3 * sigma_yearxorigin6_vector_2_3);
byearxorigin6_actual_parameters_matrix_DE[4, ] = to_array_1d(byearxorigin6_raw_vector_2_10_DE * sigma_yearxorigin6_vector_2_10_DE);
byearxorigin6_actual_parameters_matrix_NDE[4, ] = to_array_1d(byearxorigin6_raw_vector_2_10_NDE * sigma_yearxorigin6_vector_2_10_NDE);
byearxorigin6_actual_parameters_matrix_NDE[5, ] = to_array_1d(byearxorigin6_raw_vector_2_10_NDE * sigma_yearxorigin6_vector_2_10_NDE);
byearxorigin6_actual_parameters_matrix_DE[6, ] = to_array_1d(byearxorigin6_raw_vector_2_12_DE * sigma_yearxorigin6_vector_2_12_DE);
byearxorigin6_actual_parameters_matrix_NDE[6, ] = to_array_1d(byearxorigin6_raw_vector_2_12_NDE * sigma_yearxorigin6_vector_2_12_NDE);
byearxorigin6_actual_parameters_matrix_DE[7, ] = to_array_1d(byearxorigin6_raw_vector_2_14_DE * sigma_yearxorigin6_vector_2_14_DE);
byearxorigin6_actual_parameters_matrix_NDE[7, ] = to_array_1d(byearxorigin6_raw_vector_2_14_NDE * sigma_yearxorigin6_vector_2_14_NDE);
byearxorigin6_actual_parameters_matrix_NDE[8, ] = to_array_1d(byearxorigin6_raw_vector_2_14_NDE * sigma_yearxorigin6_vector_2_14_NDE);
byearxorigin6_actual_parameters_matrix_DE[9, ] = to_array_1d(byearxorigin6_raw_vector_2_16_DE * sigma_yearxorigin6_vector_2_16_DE);
byearxorigin6_actual_parameters_matrix_NDE[9, ] = to_array_1d(byearxorigin6_raw_vector_2_16_NDE * sigma_yearxorigin6_vector_2_16_NDE);
byearxorigin6_actual_parameters_matrix_NDE[10, ] = to_array_1d(byearxorigin6_raw_vector_2_16_NDE * sigma_yearxorigin6_vector_2_16_NDE);
byearxorigin6_actual_parameters_matrix_DE[11, ] = to_array_1d(byearxorigin6_raw_vector_2_39 * sigma_yearxorigin6_vector_2_39);
byearxorigin6_actual_parameters_matrix_NDE[11, ] = to_array_1d(byearxorigin6_raw_vector_2_39 * sigma_yearxorigin6_vector_2_39);
byearxorigin6_actual_parameters_matrix_DE[12, ] = to_array_1d(byearxorigin6_raw_vector_3_2 * sigma_yearxorigin6_vector_3_2);
byearxorigin6_actual_parameters_matrix_NDE[12, ] = to_array_1d(byearxorigin6_raw_vector_3_2 * sigma_yearxorigin6_vector_3_2);
byearxorigin6_actual_parameters_matrix_DE[13, ] = to_array_1d(byearxorigin6_raw_vector_3_4 * sigma_yearxorigin6_vector_3_4);
byearxorigin6_actual_parameters_matrix_NDE[13, ] = to_array_1d(byearxorigin6_raw_vector_3_4 * sigma_yearxorigin6_vector_3_4);
byearxorigin6_actual_parameters_matrix_DE[14, ] = to_array_1d(byearxorigin6_raw_vector_3_8 * sigma_yearxorigin6_vector_3_8);
byearxorigin6_actual_parameters_matrix_NDE[14, ] = to_array_1d(byearxorigin6_raw_vector_3_8 * sigma_yearxorigin6_vector_3_8);
byearxorigin6_actual_parameters_matrix_DE[15, ] = to_array_1d(byearxorigin6_raw_vector_3_18_DE * sigma_yearxorigin6_vector_3_18_DE);
byearxorigin6_actual_parameters_matrix_NDE[15, ] = to_array_1d(byearxorigin6_raw_vector_3_18_NDE * sigma_yearxorigin6_vector_3_18_NDE);
byearxorigin6_actual_parameters_matrix_DE[16, ] = to_array_1d(byearxorigin6_raw_vector_3_20_DE * sigma_yearxorigin6_vector_3_20_DE);
byearxorigin6_actual_parameters_matrix_NDE[16, ] = to_array_1d(byearxorigin6_raw_vector_3_20_NDE * sigma_yearxorigin6_vector_3_20_NDE);
byearxorigin6_actual_parameters_matrix_DE[17, ] = to_array_1d(byearxorigin6_raw_vector_4_3 * sigma_yearxorigin6_vector_4_3);
byearxorigin6_actual_parameters_matrix_NDE[17, ] = to_array_1d(byearxorigin6_raw_vector_4_3 * sigma_yearxorigin6_vector_4_3);
byearxorigin6_actual_parameters_matrix_DE[27, ] = to_array_1d(byearxorigin6_raw_vector_8_3 * sigma_yearxorigin6_vector_8_3);
byearxorigin6_actual_parameters_matrix_NDE[27, ] = to_array_1d(byearxorigin6_raw_vector_8_3 * sigma_yearxorigin6_vector_8_3);

// 3. Store the sigma_year parameters in the vector that contains all of them
// I think this is necessary so that sigma_year_matrix_DE[,] can be passed
// as an argument to reduce_sum, and therefore those parameters can change
sigma_yearxorigin1_vector_DE[1] = sigma_yearxorigin1_vector_1_2;
sigma_yearxorigin1_vector_NDE[1] = sigma_yearxorigin1_vector_1_2;
sigma_yearxorigin1_vector_DE[2] = sigma_yearxorigin1_vector_2_1;
sigma_yearxorigin1_vector_NDE[2] = sigma_yearxorigin1_vector_2_1;
sigma_yearxorigin1_vector_DE[3] = sigma_yearxorigin1_vector_2_3;
sigma_yearxorigin1_vector_NDE[3] = sigma_yearxorigin1_vector_2_3;
sigma_yearxorigin1_vector_DE[4] = sigma_yearxorigin1_vector_2_10_DE;
sigma_yearxorigin1_vector_NDE[4] = sigma_yearxorigin1_vector_2_10_NDE;
sigma_yearxorigin1_vector_NDE[5] = sigma_yearxorigin1_vector_2_10_NDE;
sigma_yearxorigin1_vector_DE[6] = sigma_yearxorigin1_vector_2_12_DE;
sigma_yearxorigin1_vector_NDE[6] = sigma_yearxorigin1_vector_2_12_NDE;
sigma_yearxorigin1_vector_DE[7] = sigma_yearxorigin1_vector_2_14_DE;
sigma_yearxorigin1_vector_NDE[7] = sigma_yearxorigin1_vector_2_14_NDE;
sigma_yearxorigin1_vector_NDE[8] = sigma_yearxorigin1_vector_2_14_NDE;
sigma_yearxorigin1_vector_DE[9] = sigma_yearxorigin1_vector_2_16_DE;
sigma_yearxorigin1_vector_NDE[9] = sigma_yearxorigin1_vector_2_16_NDE;
sigma_yearxorigin1_vector_NDE[10] = sigma_yearxorigin1_vector_2_16_NDE;
sigma_yearxorigin1_vector_DE[11] = sigma_yearxorigin1_vector_2_39;
sigma_yearxorigin1_vector_NDE[11] = sigma_yearxorigin1_vector_2_39;
sigma_yearxorigin1_vector_DE[12] = sigma_yearxorigin1_vector_3_2;
sigma_yearxorigin1_vector_NDE[12] = sigma_yearxorigin1_vector_3_2;
sigma_yearxorigin1_vector_DE[13] = sigma_yearxorigin1_vector_3_4;
sigma_yearxorigin1_vector_NDE[13] = sigma_yearxorigin1_vector_3_4;
sigma_yearxorigin1_vector_DE[14] = sigma_yearxorigin1_vector_3_8;
sigma_yearxorigin1_vector_NDE[14] = sigma_yearxorigin1_vector_3_8;
sigma_yearxorigin1_vector_DE[15] = sigma_yearxorigin1_vector_3_18_DE;
sigma_yearxorigin1_vector_NDE[15] = sigma_yearxorigin1_vector_3_18_NDE;
sigma_yearxorigin1_vector_DE[16] = sigma_yearxorigin1_vector_3_20_DE;
sigma_yearxorigin1_vector_NDE[16] = sigma_yearxorigin1_vector_3_20_NDE;
sigma_yearxorigin1_vector_DE[17] = sigma_yearxorigin1_vector_4_3;
sigma_yearxorigin1_vector_NDE[17] = sigma_yearxorigin1_vector_4_3;
sigma_yearxorigin1_vector_DE[27] = sigma_yearxorigin1_vector_8_3;
sigma_yearxorigin1_vector_NDE[27] = sigma_yearxorigin1_vector_8_3;

sigma_yearxorigin2_vector_DE[1] = sigma_yearxorigin2_vector_1_2;
sigma_yearxorigin2_vector_NDE[1] = sigma_yearxorigin2_vector_1_2;
sigma_yearxorigin2_vector_DE[2] = sigma_yearxorigin2_vector_2_1;
sigma_yearxorigin2_vector_NDE[2] = sigma_yearxorigin2_vector_2_1;
sigma_yearxorigin2_vector_DE[3] = sigma_yearxorigin2_vector_2_3;
sigma_yearxorigin2_vector_NDE[3] = sigma_yearxorigin2_vector_2_3;
sigma_yearxorigin2_vector_DE[4] = sigma_yearxorigin2_vector_2_10_DE;
sigma_yearxorigin2_vector_NDE[4] = sigma_yearxorigin2_vector_2_10_NDE;
sigma_yearxorigin2_vector_NDE[5] = sigma_yearxorigin2_vector_2_10_NDE;
sigma_yearxorigin2_vector_DE[6] = sigma_yearxorigin2_vector_2_12_DE;
sigma_yearxorigin2_vector_NDE[6] = sigma_yearxorigin2_vector_2_12_NDE;
sigma_yearxorigin2_vector_DE[7] = sigma_yearxorigin2_vector_2_14_DE;
sigma_yearxorigin2_vector_NDE[7] = sigma_yearxorigin2_vector_2_14_NDE;
sigma_yearxorigin2_vector_NDE[8] = sigma_yearxorigin2_vector_2_14_NDE;
sigma_yearxorigin2_vector_DE[9] = sigma_yearxorigin2_vector_2_16_DE;
sigma_yearxorigin2_vector_NDE[9] = sigma_yearxorigin2_vector_2_16_NDE;
sigma_yearxorigin2_vector_NDE[10] = sigma_yearxorigin2_vector_2_16_NDE;
sigma_yearxorigin2_vector_DE[11] = sigma_yearxorigin2_vector_2_39;
sigma_yearxorigin2_vector_NDE[11] = sigma_yearxorigin2_vector_2_39;
sigma_yearxorigin2_vector_DE[12] = sigma_yearxorigin2_vector_3_2;
sigma_yearxorigin2_vector_NDE[12] = sigma_yearxorigin2_vector_3_2;
sigma_yearxorigin2_vector_DE[13] = sigma_yearxorigin2_vector_3_4;
sigma_yearxorigin2_vector_NDE[13] = sigma_yearxorigin2_vector_3_4;
sigma_yearxorigin2_vector_DE[14] = sigma_yearxorigin2_vector_3_8;
sigma_yearxorigin2_vector_NDE[14] = sigma_yearxorigin2_vector_3_8;
sigma_yearxorigin2_vector_DE[15] = sigma_yearxorigin2_vector_3_18_DE;
sigma_yearxorigin2_vector_NDE[15] = sigma_yearxorigin2_vector_3_18_NDE;
sigma_yearxorigin2_vector_DE[16] = sigma_yearxorigin2_vector_3_20_DE;
sigma_yearxorigin2_vector_NDE[16] = sigma_yearxorigin2_vector_3_20_NDE;
sigma_yearxorigin2_vector_DE[17] = sigma_yearxorigin2_vector_4_3;
sigma_yearxorigin2_vector_NDE[17] = sigma_yearxorigin2_vector_4_3;
sigma_yearxorigin2_vector_DE[27] = sigma_yearxorigin2_vector_8_3;
sigma_yearxorigin2_vector_NDE[27] = sigma_yearxorigin2_vector_8_3;

sigma_yearxorigin3_vector_DE[1] = sigma_yearxorigin3_vector_1_2;
sigma_yearxorigin3_vector_NDE[1] = sigma_yearxorigin3_vector_1_2;
sigma_yearxorigin3_vector_DE[2] = sigma_yearxorigin3_vector_2_1;
sigma_yearxorigin3_vector_NDE[2] = sigma_yearxorigin3_vector_2_1;
sigma_yearxorigin3_vector_DE[3] = sigma_yearxorigin3_vector_2_3;
sigma_yearxorigin3_vector_NDE[3] = sigma_yearxorigin3_vector_2_3;
sigma_yearxorigin3_vector_DE[4] = sigma_yearxorigin3_vector_2_10_DE;
sigma_yearxorigin3_vector_NDE[4] = sigma_yearxorigin3_vector_2_10_NDE;
sigma_yearxorigin3_vector_NDE[5] = sigma_yearxorigin3_vector_2_10_NDE;
sigma_yearxorigin3_vector_DE[6] = sigma_yearxorigin3_vector_2_12_DE;
sigma_yearxorigin3_vector_NDE[6] = sigma_yearxorigin3_vector_2_12_NDE;
sigma_yearxorigin3_vector_DE[7] = sigma_yearxorigin3_vector_2_14_DE;
sigma_yearxorigin3_vector_NDE[7] = sigma_yearxorigin3_vector_2_14_NDE;
sigma_yearxorigin3_vector_NDE[8] = sigma_yearxorigin3_vector_2_14_NDE;
sigma_yearxorigin3_vector_DE[9] = sigma_yearxorigin3_vector_2_16_DE;
sigma_yearxorigin3_vector_NDE[9] = sigma_yearxorigin3_vector_2_16_NDE;
sigma_yearxorigin3_vector_NDE[10] = sigma_yearxorigin3_vector_2_16_NDE;
sigma_yearxorigin3_vector_DE[11] = sigma_yearxorigin3_vector_2_39;
sigma_yearxorigin3_vector_NDE[11] = sigma_yearxorigin3_vector_2_39;
sigma_yearxorigin3_vector_DE[12] = sigma_yearxorigin3_vector_3_2;
sigma_yearxorigin3_vector_NDE[12] = sigma_yearxorigin3_vector_3_2;
sigma_yearxorigin3_vector_DE[13] = sigma_yearxorigin3_vector_3_4;
sigma_yearxorigin3_vector_NDE[13] = sigma_yearxorigin3_vector_3_4;
sigma_yearxorigin3_vector_DE[14] = sigma_yearxorigin3_vector_3_8;
sigma_yearxorigin3_vector_NDE[14] = sigma_yearxorigin3_vector_3_8;
sigma_yearxorigin3_vector_DE[15] = sigma_yearxorigin3_vector_3_18_DE;
sigma_yearxorigin3_vector_NDE[15] = sigma_yearxorigin3_vector_3_18_NDE;
sigma_yearxorigin3_vector_DE[16] = sigma_yearxorigin3_vector_3_20_DE;
sigma_yearxorigin3_vector_NDE[16] = sigma_yearxorigin3_vector_3_20_NDE;
sigma_yearxorigin3_vector_DE[17] = sigma_yearxorigin3_vector_4_3;
sigma_yearxorigin3_vector_NDE[17] = sigma_yearxorigin3_vector_4_3;
sigma_yearxorigin3_vector_DE[27] = sigma_yearxorigin3_vector_8_3;
sigma_yearxorigin3_vector_NDE[27] = sigma_yearxorigin3_vector_8_3;

sigma_yearxorigin4_vector_DE[1] = sigma_yearxorigin4_vector_1_2;
sigma_yearxorigin4_vector_NDE[1] = sigma_yearxorigin4_vector_1_2;
sigma_yearxorigin4_vector_DE[2] = sigma_yearxorigin4_vector_2_1;
sigma_yearxorigin4_vector_NDE[2] = sigma_yearxorigin4_vector_2_1;
sigma_yearxorigin4_vector_DE[3] = sigma_yearxorigin4_vector_2_3;
sigma_yearxorigin4_vector_NDE[3] = sigma_yearxorigin4_vector_2_3;
sigma_yearxorigin4_vector_DE[4] = sigma_yearxorigin4_vector_2_10_DE;
sigma_yearxorigin4_vector_NDE[4] = sigma_yearxorigin4_vector_2_10_NDE;
sigma_yearxorigin4_vector_NDE[5] = sigma_yearxorigin4_vector_2_10_NDE;
sigma_yearxorigin4_vector_DE[6] = sigma_yearxorigin4_vector_2_12_DE;
sigma_yearxorigin4_vector_NDE[6] = sigma_yearxorigin4_vector_2_12_NDE;
sigma_yearxorigin4_vector_DE[7] = sigma_yearxorigin4_vector_2_14_DE;
sigma_yearxorigin4_vector_NDE[7] = sigma_yearxorigin4_vector_2_14_NDE;
sigma_yearxorigin4_vector_NDE[8] = sigma_yearxorigin4_vector_2_14_NDE;
sigma_yearxorigin4_vector_DE[9] = sigma_yearxorigin4_vector_2_16_DE;
sigma_yearxorigin4_vector_NDE[9] = sigma_yearxorigin4_vector_2_16_NDE;
sigma_yearxorigin4_vector_NDE[10] = sigma_yearxorigin4_vector_2_16_NDE;
sigma_yearxorigin4_vector_DE[11] = sigma_yearxorigin4_vector_2_39;
sigma_yearxorigin4_vector_NDE[11] = sigma_yearxorigin4_vector_2_39;
sigma_yearxorigin4_vector_DE[12] = sigma_yearxorigin4_vector_3_2;
sigma_yearxorigin4_vector_NDE[12] = sigma_yearxorigin4_vector_3_2;
sigma_yearxorigin4_vector_DE[13] = sigma_yearxorigin4_vector_3_4;
sigma_yearxorigin4_vector_NDE[13] = sigma_yearxorigin4_vector_3_4;
sigma_yearxorigin4_vector_DE[14] = sigma_yearxorigin4_vector_3_8;
sigma_yearxorigin4_vector_NDE[14] = sigma_yearxorigin4_vector_3_8;
sigma_yearxorigin4_vector_DE[15] = sigma_yearxorigin4_vector_3_18_DE;
sigma_yearxorigin4_vector_NDE[15] = sigma_yearxorigin4_vector_3_18_NDE;
sigma_yearxorigin4_vector_DE[16] = sigma_yearxorigin4_vector_3_20_DE;
sigma_yearxorigin4_vector_NDE[16] = sigma_yearxorigin4_vector_3_20_NDE;
sigma_yearxorigin4_vector_DE[17] = sigma_yearxorigin4_vector_4_3;
sigma_yearxorigin4_vector_NDE[17] = sigma_yearxorigin4_vector_4_3;
sigma_yearxorigin4_vector_DE[27] = sigma_yearxorigin4_vector_8_3;
sigma_yearxorigin4_vector_NDE[27] = sigma_yearxorigin4_vector_8_3;

sigma_yearxorigin5_vector_DE[1] = sigma_yearxorigin5_vector_1_2;
sigma_yearxorigin5_vector_NDE[1] = sigma_yearxorigin5_vector_1_2;
sigma_yearxorigin5_vector_DE[2] = sigma_yearxorigin5_vector_2_1;
sigma_yearxorigin5_vector_NDE[2] = sigma_yearxorigin5_vector_2_1;
sigma_yearxorigin5_vector_DE[3] = sigma_yearxorigin5_vector_2_3;
sigma_yearxorigin5_vector_NDE[3] = sigma_yearxorigin5_vector_2_3;
sigma_yearxorigin5_vector_DE[4] = sigma_yearxorigin5_vector_2_10_DE;
sigma_yearxorigin5_vector_NDE[4] = sigma_yearxorigin5_vector_2_10_NDE;
sigma_yearxorigin5_vector_NDE[5] = sigma_yearxorigin5_vector_2_10_NDE;
sigma_yearxorigin5_vector_DE[6] = sigma_yearxorigin5_vector_2_12_DE;
sigma_yearxorigin5_vector_NDE[6] = sigma_yearxorigin5_vector_2_12_NDE;
sigma_yearxorigin5_vector_DE[7] = sigma_yearxorigin5_vector_2_14_DE;
sigma_yearxorigin5_vector_NDE[7] = sigma_yearxorigin5_vector_2_14_NDE;
sigma_yearxorigin5_vector_NDE[8] = sigma_yearxorigin5_vector_2_14_NDE;
sigma_yearxorigin5_vector_DE[9] = sigma_yearxorigin5_vector_2_16_DE;
sigma_yearxorigin5_vector_NDE[9] = sigma_yearxorigin5_vector_2_16_NDE;
sigma_yearxorigin5_vector_NDE[10] = sigma_yearxorigin5_vector_2_16_NDE;
sigma_yearxorigin5_vector_DE[11] = sigma_yearxorigin5_vector_2_39;
sigma_yearxorigin5_vector_NDE[11] = sigma_yearxorigin5_vector_2_39;
sigma_yearxorigin5_vector_DE[12] = sigma_yearxorigin5_vector_3_2;
sigma_yearxorigin5_vector_NDE[12] = sigma_yearxorigin5_vector_3_2;
sigma_yearxorigin5_vector_DE[13] = sigma_yearxorigin5_vector_3_4;
sigma_yearxorigin5_vector_NDE[13] = sigma_yearxorigin5_vector_3_4;
sigma_yearxorigin5_vector_DE[14] = sigma_yearxorigin5_vector_3_8;
sigma_yearxorigin5_vector_NDE[14] = sigma_yearxorigin5_vector_3_8;
sigma_yearxorigin5_vector_DE[15] = sigma_yearxorigin5_vector_3_18_DE;
sigma_yearxorigin5_vector_NDE[15] = sigma_yearxorigin5_vector_3_18_NDE;
sigma_yearxorigin5_vector_DE[16] = sigma_yearxorigin5_vector_3_20_DE;
sigma_yearxorigin5_vector_NDE[16] = sigma_yearxorigin5_vector_3_20_NDE;
sigma_yearxorigin5_vector_DE[17] = sigma_yearxorigin5_vector_4_3;
sigma_yearxorigin5_vector_NDE[17] = sigma_yearxorigin5_vector_4_3;
sigma_yearxorigin5_vector_DE[27] = sigma_yearxorigin5_vector_8_3;
sigma_yearxorigin5_vector_NDE[27] = sigma_yearxorigin5_vector_8_3;

sigma_yearxorigin6_vector_DE[1] = sigma_yearxorigin6_vector_1_2;
sigma_yearxorigin6_vector_NDE[1] = sigma_yearxorigin6_vector_1_2;
sigma_yearxorigin6_vector_DE[2] = sigma_yearxorigin6_vector_2_1;
sigma_yearxorigin6_vector_NDE[2] = sigma_yearxorigin6_vector_2_1;
sigma_yearxorigin6_vector_DE[3] = sigma_yearxorigin6_vector_2_3;
sigma_yearxorigin6_vector_NDE[3] = sigma_yearxorigin6_vector_2_3;
sigma_yearxorigin6_vector_DE[4] = sigma_yearxorigin6_vector_2_10_DE;
sigma_yearxorigin6_vector_NDE[4] = sigma_yearxorigin6_vector_2_10_NDE;
sigma_yearxorigin6_vector_NDE[5] = sigma_yearxorigin6_vector_2_10_NDE;
sigma_yearxorigin6_vector_DE[6] = sigma_yearxorigin6_vector_2_12_DE;
sigma_yearxorigin6_vector_NDE[6] = sigma_yearxorigin6_vector_2_12_NDE;
sigma_yearxorigin6_vector_DE[7] = sigma_yearxorigin6_vector_2_14_DE;
sigma_yearxorigin6_vector_NDE[7] = sigma_yearxorigin6_vector_2_14_NDE;
sigma_yearxorigin6_vector_NDE[8] = sigma_yearxorigin6_vector_2_14_NDE;
sigma_yearxorigin6_vector_DE[9] = sigma_yearxorigin6_vector_2_16_DE;
sigma_yearxorigin6_vector_NDE[9] = sigma_yearxorigin6_vector_2_16_NDE;
sigma_yearxorigin6_vector_NDE[10] = sigma_yearxorigin6_vector_2_16_NDE;
sigma_yearxorigin6_vector_DE[11] = sigma_yearxorigin6_vector_2_39;
sigma_yearxorigin6_vector_NDE[11] = sigma_yearxorigin6_vector_2_39;
sigma_yearxorigin6_vector_DE[12] = sigma_yearxorigin6_vector_3_2;
sigma_yearxorigin6_vector_NDE[12] = sigma_yearxorigin6_vector_3_2;
sigma_yearxorigin6_vector_DE[13] = sigma_yearxorigin6_vector_3_4;
sigma_yearxorigin6_vector_NDE[13] = sigma_yearxorigin6_vector_3_4;
sigma_yearxorigin6_vector_DE[14] = sigma_yearxorigin6_vector_3_8;
sigma_yearxorigin6_vector_NDE[14] = sigma_yearxorigin6_vector_3_8;
sigma_yearxorigin6_vector_DE[15] = sigma_yearxorigin6_vector_3_18_DE;
sigma_yearxorigin6_vector_NDE[15] = sigma_yearxorigin6_vector_3_18_NDE;
sigma_yearxorigin6_vector_DE[16] = sigma_yearxorigin6_vector_3_20_DE;
sigma_yearxorigin6_vector_NDE[16] = sigma_yearxorigin6_vector_3_20_NDE;
sigma_yearxorigin6_vector_DE[17] = sigma_yearxorigin6_vector_4_3;
sigma_yearxorigin6_vector_NDE[17] = sigma_yearxorigin6_vector_4_3;
sigma_yearxorigin6_vector_DE[27] = sigma_yearxorigin6_vector_8_3;
sigma_yearxorigin6_vector_NDE[27] = sigma_yearxorigin6_vector_8_3;

// Temperature parameters
btemp0_vector_DE[18] = btemp0_matrix_4_5;
btemp0_vector_NDE[18] = btemp0_matrix_4_5;
btemp1_vector_DE[18] = btemp1_matrix_4_5;
btemp1_vector_NDE[18] = btemp1_matrix_4_5;
btemp0_vector_DE[20] = btemp0_matrix_5_6;
btemp0_vector_NDE[20] = btemp0_matrix_5_6;
btemp1_vector_DE[20] = btemp1_matrix_5_6;
btemp1_vector_NDE[20] = btemp1_matrix_5_6;
btemp0_vector_DE[21] = btemp0_matrix_5_22_DE;
btemp0_vector_NDE[21] = btemp0_matrix_5_22_NDE;
btemp1_vector_DE[21] = btemp1_matrix_5_22_DE;
btemp1_vector_NDE[21] = btemp1_matrix_5_22_NDE;
btemp0_vector_DE[23] = btemp0_matrix_6_7;
btemp0_vector_NDE[23] = btemp0_matrix_6_7;
btemp1_vector_DE[23] = btemp1_matrix_6_7;
btemp1_vector_NDE[23] = btemp1_matrix_6_7;
btemp0_vector_DE[25] = btemp0_matrix_7_26_DE;
btemp0_vector_NDE[25] = btemp0_matrix_7_26_NDE;
btemp1_vector_DE[25] = btemp1_matrix_7_26_DE;
btemp1_vector_NDE[25] = btemp1_matrix_7_26_NDE;
btemp0_vector_DE[26] = btemp0_matrix_7_28_DE;
btemp0_vector_NDE[26] = btemp0_matrix_7_28_NDE;
btemp1_vector_DE[26] = btemp1_matrix_7_28_DE;
btemp1_vector_NDE[26] = btemp1_matrix_7_28_NDE;
btemp0_vector_DE[28] = btemp0_matrix_8_9;
btemp0_vector_NDE[28] = btemp0_matrix_8_9;
btemp1_vector_DE[28] = btemp1_matrix_8_9;
btemp1_vector_NDE[28] = btemp1_matrix_8_9;
btemp0_vector_DE[29] = btemp0_matrix_8_30_DE;
btemp0_vector_NDE[29] = btemp0_matrix_8_30_NDE;
btemp1_vector_DE[29] = btemp1_matrix_8_30_DE;
btemp1_vector_NDE[29] = btemp1_matrix_8_30_NDE;
btemp0_vector_DE[31] = btemp0_matrix_9_32_DE;
btemp0_vector_NDE[31] = btemp0_matrix_9_32_NDE;
btemp1_vector_DE[31] = btemp1_matrix_9_32_DE;
btemp1_vector_NDE[31] = btemp1_matrix_9_32_NDE;
btemp0_vector_DE[32] = btemp0_matrix_9_34;
btemp0_vector_NDE[32] = btemp0_matrix_9_34;
btemp1_vector_DE[32] = btemp1_matrix_9_34;
btemp1_vector_NDE[32] = btemp1_matrix_9_34;
btemp0_vector_DE[33] = btemp0_matrix_9_35;
btemp0_vector_NDE[33] = btemp0_matrix_9_35;
btemp1_vector_DE[33] = btemp1_matrix_9_35;
btemp1_vector_NDE[33] = btemp1_matrix_9_35;
btemp0_vector_DE[34] = btemp0_matrix_9_36;
btemp0_vector_NDE[34] = btemp0_matrix_9_36;
btemp1_vector_DE[34] = btemp1_matrix_9_36;
btemp1_vector_NDE[34] = btemp1_matrix_9_36;
btemp0_vector_DE[35] = btemp0_matrix_9_37_DE;
btemp0_vector_NDE[35] = btemp0_matrix_9_37_NDE;
btemp1_vector_DE[35] = btemp1_matrix_9_37_DE;
btemp1_vector_NDE[35] = btemp1_matrix_9_37_NDE;

btemp0xorigin1_vector_DE[1] = btemp0xorigin1_matrix_1_2;
btemp0xorigin1_vector_NDE[1] = btemp0xorigin1_matrix_1_2;
btemp1xorigin1_vector_DE[1] = btemp1xorigin1_matrix_1_2;
btemp1xorigin1_vector_NDE[1] = btemp1xorigin1_matrix_1_2;
btemp0xorigin1_vector_DE[3] = btemp0xorigin1_matrix_2_3;
btemp0xorigin1_vector_NDE[3] = btemp0xorigin1_matrix_2_3;
btemp1xorigin1_vector_DE[3] = btemp1xorigin1_matrix_2_3;
btemp1xorigin1_vector_NDE[3] = btemp1xorigin1_matrix_2_3;
btemp0xorigin1_vector_DE[4] = btemp0xorigin1_matrix_2_10_DE;
btemp0xorigin1_vector_NDE[4] = btemp0xorigin1_matrix_2_10_NDE;
btemp1xorigin1_vector_DE[4] = btemp1xorigin1_matrix_2_10_DE;
btemp1xorigin1_vector_NDE[4] = btemp1xorigin1_matrix_2_10_NDE;
btemp0xorigin1_vector_NDE[5] = btemp0xorigin1_matrix_2_10_NDE;
btemp1xorigin1_vector_NDE[5] = btemp1xorigin1_matrix_2_10_NDE;
btemp0xorigin1_vector_DE[6] = btemp0xorigin1_matrix_2_12_DE;
btemp0xorigin1_vector_NDE[6] = btemp0xorigin1_matrix_2_12_NDE;
btemp1xorigin1_vector_DE[6] = btemp1xorigin1_matrix_2_12_DE;
btemp1xorigin1_vector_NDE[6] = btemp1xorigin1_matrix_2_12_NDE;
btemp0xorigin1_vector_DE[7] = btemp0xorigin1_matrix_2_14_DE;
btemp0xorigin1_vector_NDE[7] = btemp0xorigin1_matrix_2_14_NDE;
btemp1xorigin1_vector_DE[7] = btemp1xorigin1_matrix_2_14_DE;
btemp1xorigin1_vector_NDE[7] = btemp1xorigin1_matrix_2_14_NDE;
btemp0xorigin1_vector_NDE[8] = btemp0xorigin1_matrix_2_14_NDE;
btemp1xorigin1_vector_NDE[8] = btemp1xorigin1_matrix_2_14_NDE;
btemp0xorigin1_vector_DE[9] = btemp0xorigin1_matrix_2_16_DE;
btemp0xorigin1_vector_NDE[9] = btemp0xorigin1_matrix_2_16_NDE;
btemp1xorigin1_vector_DE[9] = btemp1xorigin1_matrix_2_16_DE;
btemp1xorigin1_vector_NDE[9] = btemp1xorigin1_matrix_2_16_NDE;
btemp0xorigin1_vector_NDE[10] = btemp0xorigin1_matrix_2_16_NDE;
btemp1xorigin1_vector_NDE[10] = btemp1xorigin1_matrix_2_16_NDE;
btemp0xorigin1_vector_DE[11] = btemp0xorigin1_matrix_2_39;
btemp0xorigin1_vector_NDE[11] = btemp0xorigin1_matrix_2_39;
btemp1xorigin1_vector_DE[11] = btemp1xorigin1_matrix_2_39;
btemp1xorigin1_vector_NDE[11] = btemp1xorigin1_matrix_2_39;
btemp0xorigin1_vector_DE[13] = btemp0xorigin1_matrix_3_4;
btemp0xorigin1_vector_NDE[13] = btemp0xorigin1_matrix_3_4;
btemp1xorigin1_vector_DE[13] = btemp1xorigin1_matrix_3_4;
btemp1xorigin1_vector_NDE[13] = btemp1xorigin1_matrix_3_4;
btemp0xorigin1_vector_DE[14] = btemp0xorigin1_matrix_3_8;
btemp0xorigin1_vector_NDE[14] = btemp0xorigin1_matrix_3_8;
btemp1xorigin1_vector_DE[14] = btemp1xorigin1_matrix_3_8;
btemp1xorigin1_vector_NDE[14] = btemp1xorigin1_matrix_3_8;
btemp0xorigin1_vector_DE[15] = btemp0xorigin1_matrix_3_18_DE;
btemp0xorigin1_vector_NDE[15] = btemp0xorigin1_matrix_3_18_NDE;
btemp1xorigin1_vector_DE[15] = btemp1xorigin1_matrix_3_18_DE;
btemp1xorigin1_vector_NDE[15] = btemp1xorigin1_matrix_3_18_NDE;
btemp0xorigin1_vector_DE[16] = btemp0xorigin1_matrix_3_20_DE;
btemp0xorigin1_vector_NDE[16] = btemp0xorigin1_matrix_3_20_NDE;
btemp1xorigin1_vector_DE[16] = btemp1xorigin1_matrix_3_20_DE;
btemp1xorigin1_vector_NDE[16] = btemp1xorigin1_matrix_3_20_NDE;

btemp0xorigin2_vector_DE[1] = btemp0xorigin2_matrix_1_2;
btemp0xorigin2_vector_NDE[1] = btemp0xorigin2_matrix_1_2;
btemp1xorigin2_vector_DE[1] = btemp1xorigin2_matrix_1_2;
btemp1xorigin2_vector_NDE[1] = btemp1xorigin2_matrix_1_2;
btemp0xorigin2_vector_DE[3] = btemp0xorigin2_matrix_2_3;
btemp0xorigin2_vector_NDE[3] = btemp0xorigin2_matrix_2_3;
btemp1xorigin2_vector_DE[3] = btemp1xorigin2_matrix_2_3;
btemp1xorigin2_vector_NDE[3] = btemp1xorigin2_matrix_2_3;
btemp0xorigin2_vector_DE[4] = btemp0xorigin2_matrix_2_10_DE;
btemp0xorigin2_vector_NDE[4] = btemp0xorigin2_matrix_2_10_NDE;
btemp1xorigin2_vector_DE[4] = btemp1xorigin2_matrix_2_10_DE;
btemp1xorigin2_vector_NDE[4] = btemp1xorigin2_matrix_2_10_NDE;
btemp0xorigin2_vector_NDE[5] = btemp0xorigin2_matrix_2_10_NDE;
btemp1xorigin2_vector_NDE[5] = btemp1xorigin2_matrix_2_10_NDE;
btemp0xorigin2_vector_DE[6] = btemp0xorigin2_matrix_2_12_DE;
btemp0xorigin2_vector_NDE[6] = btemp0xorigin2_matrix_2_12_NDE;
btemp1xorigin2_vector_DE[6] = btemp1xorigin2_matrix_2_12_DE;
btemp1xorigin2_vector_NDE[6] = btemp1xorigin2_matrix_2_12_NDE;
btemp0xorigin2_vector_DE[7] = btemp0xorigin2_matrix_2_14_DE;
btemp0xorigin2_vector_NDE[7] = btemp0xorigin2_matrix_2_14_NDE;
btemp1xorigin2_vector_DE[7] = btemp1xorigin2_matrix_2_14_DE;
btemp1xorigin2_vector_NDE[7] = btemp1xorigin2_matrix_2_14_NDE;
btemp0xorigin2_vector_NDE[8] = btemp0xorigin2_matrix_2_14_NDE;
btemp1xorigin2_vector_NDE[8] = btemp1xorigin2_matrix_2_14_NDE;
btemp0xorigin2_vector_DE[9] = btemp0xorigin2_matrix_2_16_DE;
btemp0xorigin2_vector_NDE[9] = btemp0xorigin2_matrix_2_16_NDE;
btemp1xorigin2_vector_DE[9] = btemp1xorigin2_matrix_2_16_DE;
btemp1xorigin2_vector_NDE[9] = btemp1xorigin2_matrix_2_16_NDE;
btemp0xorigin2_vector_NDE[10] = btemp0xorigin2_matrix_2_16_NDE;
btemp1xorigin2_vector_NDE[10] = btemp1xorigin2_matrix_2_16_NDE;
btemp0xorigin2_vector_DE[11] = btemp0xorigin2_matrix_2_39;
btemp0xorigin2_vector_NDE[11] = btemp0xorigin2_matrix_2_39;
btemp1xorigin2_vector_DE[11] = btemp1xorigin2_matrix_2_39;
btemp1xorigin2_vector_NDE[11] = btemp1xorigin2_matrix_2_39;
btemp0xorigin2_vector_DE[13] = btemp0xorigin2_matrix_3_4;
btemp0xorigin2_vector_NDE[13] = btemp0xorigin2_matrix_3_4;
btemp1xorigin2_vector_DE[13] = btemp1xorigin2_matrix_3_4;
btemp1xorigin2_vector_NDE[13] = btemp1xorigin2_matrix_3_4;
btemp0xorigin2_vector_DE[14] = btemp0xorigin2_matrix_3_8;
btemp0xorigin2_vector_NDE[14] = btemp0xorigin2_matrix_3_8;
btemp1xorigin2_vector_DE[14] = btemp1xorigin2_matrix_3_8;
btemp1xorigin2_vector_NDE[14] = btemp1xorigin2_matrix_3_8;
btemp0xorigin2_vector_DE[15] = btemp0xorigin2_matrix_3_18_DE;
btemp0xorigin2_vector_NDE[15] = btemp0xorigin2_matrix_3_18_NDE;
btemp1xorigin2_vector_DE[15] = btemp1xorigin2_matrix_3_18_DE;
btemp1xorigin2_vector_NDE[15] = btemp1xorigin2_matrix_3_18_NDE;
btemp0xorigin2_vector_DE[16] = btemp0xorigin2_matrix_3_20_DE;
btemp0xorigin2_vector_NDE[16] = btemp0xorigin2_matrix_3_20_NDE;
btemp1xorigin2_vector_DE[16] = btemp1xorigin2_matrix_3_20_DE;
btemp1xorigin2_vector_NDE[16] = btemp1xorigin2_matrix_3_20_NDE;

btemp0xorigin3_vector_DE[1] = btemp0xorigin3_matrix_1_2;
btemp0xorigin3_vector_NDE[1] = btemp0xorigin3_matrix_1_2;
btemp1xorigin3_vector_DE[1] = btemp1xorigin3_matrix_1_2;
btemp1xorigin3_vector_NDE[1] = btemp1xorigin3_matrix_1_2;
btemp0xorigin3_vector_DE[3] = btemp0xorigin3_matrix_2_3;
btemp0xorigin3_vector_NDE[3] = btemp0xorigin3_matrix_2_3;
btemp1xorigin3_vector_DE[3] = btemp1xorigin3_matrix_2_3;
btemp1xorigin3_vector_NDE[3] = btemp1xorigin3_matrix_2_3;
btemp0xorigin3_vector_DE[4] = btemp0xorigin3_matrix_2_10_DE;
btemp0xorigin3_vector_NDE[4] = btemp0xorigin3_matrix_2_10_NDE;
btemp1xorigin3_vector_DE[4] = btemp1xorigin3_matrix_2_10_DE;
btemp1xorigin3_vector_NDE[4] = btemp1xorigin3_matrix_2_10_NDE;
btemp0xorigin3_vector_NDE[5] = btemp0xorigin3_matrix_2_10_NDE;
btemp1xorigin3_vector_NDE[5] = btemp1xorigin3_matrix_2_10_NDE;
btemp0xorigin3_vector_DE[6] = btemp0xorigin3_matrix_2_12_DE;
btemp0xorigin3_vector_NDE[6] = btemp0xorigin3_matrix_2_12_NDE;
btemp1xorigin3_vector_DE[6] = btemp1xorigin3_matrix_2_12_DE;
btemp1xorigin3_vector_NDE[6] = btemp1xorigin3_matrix_2_12_NDE;
btemp0xorigin3_vector_DE[7] = btemp0xorigin3_matrix_2_14_DE;
btemp0xorigin3_vector_NDE[7] = btemp0xorigin3_matrix_2_14_NDE;
btemp1xorigin3_vector_DE[7] = btemp1xorigin3_matrix_2_14_DE;
btemp1xorigin3_vector_NDE[7] = btemp1xorigin3_matrix_2_14_NDE;
btemp0xorigin3_vector_NDE[8] = btemp0xorigin3_matrix_2_14_NDE;
btemp1xorigin3_vector_NDE[8] = btemp1xorigin3_matrix_2_14_NDE;
btemp0xorigin3_vector_DE[9] = btemp0xorigin3_matrix_2_16_DE;
btemp0xorigin3_vector_NDE[9] = btemp0xorigin3_matrix_2_16_NDE;
btemp1xorigin3_vector_DE[9] = btemp1xorigin3_matrix_2_16_DE;
btemp1xorigin3_vector_NDE[9] = btemp1xorigin3_matrix_2_16_NDE;
btemp0xorigin3_vector_NDE[10] = btemp0xorigin3_matrix_2_16_NDE;
btemp1xorigin3_vector_NDE[10] = btemp1xorigin3_matrix_2_16_NDE;
btemp0xorigin3_vector_DE[11] = btemp0xorigin3_matrix_2_39;
btemp0xorigin3_vector_NDE[11] = btemp0xorigin3_matrix_2_39;
btemp1xorigin3_vector_DE[11] = btemp1xorigin3_matrix_2_39;
btemp1xorigin3_vector_NDE[11] = btemp1xorigin3_matrix_2_39;
btemp0xorigin3_vector_DE[13] = btemp0xorigin3_matrix_3_4;
btemp0xorigin3_vector_NDE[13] = btemp0xorigin3_matrix_3_4;
btemp1xorigin3_vector_DE[13] = btemp1xorigin3_matrix_3_4;
btemp1xorigin3_vector_NDE[13] = btemp1xorigin3_matrix_3_4;
btemp0xorigin3_vector_DE[14] = btemp0xorigin3_matrix_3_8;
btemp0xorigin3_vector_NDE[14] = btemp0xorigin3_matrix_3_8;
btemp1xorigin3_vector_DE[14] = btemp1xorigin3_matrix_3_8;
btemp1xorigin3_vector_NDE[14] = btemp1xorigin3_matrix_3_8;
btemp0xorigin3_vector_DE[15] = btemp0xorigin3_matrix_3_18_DE;
btemp0xorigin3_vector_NDE[15] = btemp0xorigin3_matrix_3_18_NDE;
btemp1xorigin3_vector_DE[15] = btemp1xorigin3_matrix_3_18_DE;
btemp1xorigin3_vector_NDE[15] = btemp1xorigin3_matrix_3_18_NDE;
btemp0xorigin3_vector_DE[16] = btemp0xorigin3_matrix_3_20_DE;
btemp0xorigin3_vector_NDE[16] = btemp0xorigin3_matrix_3_20_NDE;
btemp1xorigin3_vector_DE[16] = btemp1xorigin3_matrix_3_20_DE;
btemp1xorigin3_vector_NDE[16] = btemp1xorigin3_matrix_3_20_NDE;

btemp0xorigin4_vector_DE[1] = btemp0xorigin4_matrix_1_2;
btemp0xorigin4_vector_NDE[1] = btemp0xorigin4_matrix_1_2;
btemp1xorigin4_vector_DE[1] = btemp1xorigin4_matrix_1_2;
btemp1xorigin4_vector_NDE[1] = btemp1xorigin4_matrix_1_2;
btemp0xorigin4_vector_DE[3] = btemp0xorigin4_matrix_2_3;
btemp0xorigin4_vector_NDE[3] = btemp0xorigin4_matrix_2_3;
btemp1xorigin4_vector_DE[3] = btemp1xorigin4_matrix_2_3;
btemp1xorigin4_vector_NDE[3] = btemp1xorigin4_matrix_2_3;
btemp0xorigin4_vector_DE[4] = btemp0xorigin4_matrix_2_10_DE;
btemp0xorigin4_vector_NDE[4] = btemp0xorigin4_matrix_2_10_NDE;
btemp1xorigin4_vector_DE[4] = btemp1xorigin4_matrix_2_10_DE;
btemp1xorigin4_vector_NDE[4] = btemp1xorigin4_matrix_2_10_NDE;
btemp0xorigin4_vector_NDE[5] = btemp0xorigin4_matrix_2_10_NDE;
btemp1xorigin4_vector_NDE[5] = btemp1xorigin4_matrix_2_10_NDE;
btemp0xorigin4_vector_DE[6] = btemp0xorigin4_matrix_2_12_DE;
btemp0xorigin4_vector_NDE[6] = btemp0xorigin4_matrix_2_12_NDE;
btemp1xorigin4_vector_DE[6] = btemp1xorigin4_matrix_2_12_DE;
btemp1xorigin4_vector_NDE[6] = btemp1xorigin4_matrix_2_12_NDE;
btemp0xorigin4_vector_DE[7] = btemp0xorigin4_matrix_2_14_DE;
btemp0xorigin4_vector_NDE[7] = btemp0xorigin4_matrix_2_14_NDE;
btemp1xorigin4_vector_DE[7] = btemp1xorigin4_matrix_2_14_DE;
btemp1xorigin4_vector_NDE[7] = btemp1xorigin4_matrix_2_14_NDE;
btemp0xorigin4_vector_NDE[8] = btemp0xorigin4_matrix_2_14_NDE;
btemp1xorigin4_vector_NDE[8] = btemp1xorigin4_matrix_2_14_NDE;
btemp0xorigin4_vector_DE[9] = btemp0xorigin4_matrix_2_16_DE;
btemp0xorigin4_vector_NDE[9] = btemp0xorigin4_matrix_2_16_NDE;
btemp1xorigin4_vector_DE[9] = btemp1xorigin4_matrix_2_16_DE;
btemp1xorigin4_vector_NDE[9] = btemp1xorigin4_matrix_2_16_NDE;
btemp0xorigin4_vector_NDE[10] = btemp0xorigin4_matrix_2_16_NDE;
btemp1xorigin4_vector_NDE[10] = btemp1xorigin4_matrix_2_16_NDE;
btemp0xorigin4_vector_DE[11] = btemp0xorigin4_matrix_2_39;
btemp0xorigin4_vector_NDE[11] = btemp0xorigin4_matrix_2_39;
btemp1xorigin4_vector_DE[11] = btemp1xorigin4_matrix_2_39;
btemp1xorigin4_vector_NDE[11] = btemp1xorigin4_matrix_2_39;
btemp0xorigin4_vector_DE[13] = btemp0xorigin4_matrix_3_4;
btemp0xorigin4_vector_NDE[13] = btemp0xorigin4_matrix_3_4;
btemp1xorigin4_vector_DE[13] = btemp1xorigin4_matrix_3_4;
btemp1xorigin4_vector_NDE[13] = btemp1xorigin4_matrix_3_4;
btemp0xorigin4_vector_DE[14] = btemp0xorigin4_matrix_3_8;
btemp0xorigin4_vector_NDE[14] = btemp0xorigin4_matrix_3_8;
btemp1xorigin4_vector_DE[14] = btemp1xorigin4_matrix_3_8;
btemp1xorigin4_vector_NDE[14] = btemp1xorigin4_matrix_3_8;
btemp0xorigin4_vector_DE[15] = btemp0xorigin4_matrix_3_18_DE;
btemp0xorigin4_vector_NDE[15] = btemp0xorigin4_matrix_3_18_NDE;
btemp1xorigin4_vector_DE[15] = btemp1xorigin4_matrix_3_18_DE;
btemp1xorigin4_vector_NDE[15] = btemp1xorigin4_matrix_3_18_NDE;
btemp0xorigin4_vector_DE[16] = btemp0xorigin4_matrix_3_20_DE;
btemp0xorigin4_vector_NDE[16] = btemp0xorigin4_matrix_3_20_NDE;
btemp1xorigin4_vector_DE[16] = btemp1xorigin4_matrix_3_20_DE;
btemp1xorigin4_vector_NDE[16] = btemp1xorigin4_matrix_3_20_NDE;

btemp0xorigin5_vector_DE[1] = btemp0xorigin5_matrix_1_2;
btemp0xorigin5_vector_NDE[1] = btemp0xorigin5_matrix_1_2;
btemp1xorigin5_vector_DE[1] = btemp1xorigin5_matrix_1_2;
btemp1xorigin5_vector_NDE[1] = btemp1xorigin5_matrix_1_2;
btemp0xorigin5_vector_DE[3] = btemp0xorigin5_matrix_2_3;
btemp0xorigin5_vector_NDE[3] = btemp0xorigin5_matrix_2_3;
btemp1xorigin5_vector_DE[3] = btemp1xorigin5_matrix_2_3;
btemp1xorigin5_vector_NDE[3] = btemp1xorigin5_matrix_2_3;
btemp0xorigin5_vector_DE[4] = btemp0xorigin5_matrix_2_10_DE;
btemp0xorigin5_vector_NDE[4] = btemp0xorigin5_matrix_2_10_NDE;
btemp1xorigin5_vector_DE[4] = btemp1xorigin5_matrix_2_10_DE;
btemp1xorigin5_vector_NDE[4] = btemp1xorigin5_matrix_2_10_NDE;
btemp0xorigin5_vector_NDE[5] = btemp0xorigin5_matrix_2_10_NDE;
btemp1xorigin5_vector_NDE[5] = btemp1xorigin5_matrix_2_10_NDE;
btemp0xorigin5_vector_DE[6] = btemp0xorigin5_matrix_2_12_DE;
btemp0xorigin5_vector_NDE[6] = btemp0xorigin5_matrix_2_12_NDE;
btemp1xorigin5_vector_DE[6] = btemp1xorigin5_matrix_2_12_DE;
btemp1xorigin5_vector_NDE[6] = btemp1xorigin5_matrix_2_12_NDE;
btemp0xorigin5_vector_DE[7] = btemp0xorigin5_matrix_2_14_DE;
btemp0xorigin5_vector_NDE[7] = btemp0xorigin5_matrix_2_14_NDE;
btemp1xorigin5_vector_DE[7] = btemp1xorigin5_matrix_2_14_DE;
btemp1xorigin5_vector_NDE[7] = btemp1xorigin5_matrix_2_14_NDE;
btemp0xorigin5_vector_NDE[8] = btemp0xorigin5_matrix_2_14_NDE;
btemp1xorigin5_vector_NDE[8] = btemp1xorigin5_matrix_2_14_NDE;
btemp0xorigin5_vector_DE[9] = btemp0xorigin5_matrix_2_16_DE;
btemp0xorigin5_vector_NDE[9] = btemp0xorigin5_matrix_2_16_NDE;
btemp1xorigin5_vector_DE[9] = btemp1xorigin5_matrix_2_16_DE;
btemp1xorigin5_vector_NDE[9] = btemp1xorigin5_matrix_2_16_NDE;
btemp0xorigin5_vector_NDE[10] = btemp0xorigin5_matrix_2_16_NDE;
btemp1xorigin5_vector_NDE[10] = btemp1xorigin5_matrix_2_16_NDE;
btemp0xorigin5_vector_DE[11] = btemp0xorigin5_matrix_2_39;
btemp0xorigin5_vector_NDE[11] = btemp0xorigin5_matrix_2_39;
btemp1xorigin5_vector_DE[11] = btemp1xorigin5_matrix_2_39;
btemp1xorigin5_vector_NDE[11] = btemp1xorigin5_matrix_2_39;
btemp0xorigin5_vector_DE[13] = btemp0xorigin5_matrix_3_4;
btemp0xorigin5_vector_NDE[13] = btemp0xorigin5_matrix_3_4;
btemp1xorigin5_vector_DE[13] = btemp1xorigin5_matrix_3_4;
btemp1xorigin5_vector_NDE[13] = btemp1xorigin5_matrix_3_4;
btemp0xorigin5_vector_DE[14] = btemp0xorigin5_matrix_3_8;
btemp0xorigin5_vector_NDE[14] = btemp0xorigin5_matrix_3_8;
btemp1xorigin5_vector_DE[14] = btemp1xorigin5_matrix_3_8;
btemp1xorigin5_vector_NDE[14] = btemp1xorigin5_matrix_3_8;
btemp0xorigin5_vector_DE[15] = btemp0xorigin5_matrix_3_18_DE;
btemp0xorigin5_vector_NDE[15] = btemp0xorigin5_matrix_3_18_NDE;
btemp1xorigin5_vector_DE[15] = btemp1xorigin5_matrix_3_18_DE;
btemp1xorigin5_vector_NDE[15] = btemp1xorigin5_matrix_3_18_NDE;
btemp0xorigin5_vector_DE[16] = btemp0xorigin5_matrix_3_20_DE;
btemp0xorigin5_vector_NDE[16] = btemp0xorigin5_matrix_3_20_NDE;
btemp1xorigin5_vector_DE[16] = btemp1xorigin5_matrix_3_20_DE;
btemp1xorigin5_vector_NDE[16] = btemp1xorigin5_matrix_3_20_NDE;

btemp0xorigin6_vector_DE[1] = btemp0xorigin6_matrix_1_2;
btemp0xorigin6_vector_NDE[1] = btemp0xorigin6_matrix_1_2;
btemp1xorigin6_vector_DE[1] = btemp1xorigin6_matrix_1_2;
btemp1xorigin6_vector_NDE[1] = btemp1xorigin6_matrix_1_2;
btemp0xorigin6_vector_DE[3] = btemp0xorigin6_matrix_2_3;
btemp0xorigin6_vector_NDE[3] = btemp0xorigin6_matrix_2_3;
btemp1xorigin6_vector_DE[3] = btemp1xorigin6_matrix_2_3;
btemp1xorigin6_vector_NDE[3] = btemp1xorigin6_matrix_2_3;
btemp0xorigin6_vector_DE[4] = btemp0xorigin6_matrix_2_10_DE;
btemp0xorigin6_vector_NDE[4] = btemp0xorigin6_matrix_2_10_NDE;
btemp1xorigin6_vector_DE[4] = btemp1xorigin6_matrix_2_10_DE;
btemp1xorigin6_vector_NDE[4] = btemp1xorigin6_matrix_2_10_NDE;
btemp0xorigin6_vector_NDE[5] = btemp0xorigin6_matrix_2_10_NDE;
btemp1xorigin6_vector_NDE[5] = btemp1xorigin6_matrix_2_10_NDE;
btemp0xorigin6_vector_DE[6] = btemp0xorigin6_matrix_2_12_DE;
btemp0xorigin6_vector_NDE[6] = btemp0xorigin6_matrix_2_12_NDE;
btemp1xorigin6_vector_DE[6] = btemp1xorigin6_matrix_2_12_DE;
btemp1xorigin6_vector_NDE[6] = btemp1xorigin6_matrix_2_12_NDE;
btemp0xorigin6_vector_DE[7] = btemp0xorigin6_matrix_2_14_DE;
btemp0xorigin6_vector_NDE[7] = btemp0xorigin6_matrix_2_14_NDE;
btemp1xorigin6_vector_DE[7] = btemp1xorigin6_matrix_2_14_DE;
btemp1xorigin6_vector_NDE[7] = btemp1xorigin6_matrix_2_14_NDE;
btemp0xorigin6_vector_NDE[8] = btemp0xorigin6_matrix_2_14_NDE;
btemp1xorigin6_vector_NDE[8] = btemp1xorigin6_matrix_2_14_NDE;
btemp0xorigin6_vector_DE[9] = btemp0xorigin6_matrix_2_16_DE;
btemp0xorigin6_vector_NDE[9] = btemp0xorigin6_matrix_2_16_NDE;
btemp1xorigin6_vector_DE[9] = btemp1xorigin6_matrix_2_16_DE;
btemp1xorigin6_vector_NDE[9] = btemp1xorigin6_matrix_2_16_NDE;
btemp0xorigin6_vector_NDE[10] = btemp0xorigin6_matrix_2_16_NDE;
btemp1xorigin6_vector_NDE[10] = btemp1xorigin6_matrix_2_16_NDE;
btemp0xorigin6_vector_DE[11] = btemp0xorigin6_matrix_2_39;
btemp0xorigin6_vector_NDE[11] = btemp0xorigin6_matrix_2_39;
btemp1xorigin6_vector_DE[11] = btemp1xorigin6_matrix_2_39;
btemp1xorigin6_vector_NDE[11] = btemp1xorigin6_matrix_2_39;
btemp0xorigin6_vector_DE[13] = btemp0xorigin6_matrix_3_4;
btemp0xorigin6_vector_NDE[13] = btemp0xorigin6_matrix_3_4;
btemp1xorigin6_vector_DE[13] = btemp1xorigin6_matrix_3_4;
btemp1xorigin6_vector_NDE[13] = btemp1xorigin6_matrix_3_4;
btemp0xorigin6_vector_DE[14] = btemp0xorigin6_matrix_3_8;
btemp0xorigin6_vector_NDE[14] = btemp0xorigin6_matrix_3_8;
btemp1xorigin6_vector_DE[14] = btemp1xorigin6_matrix_3_8;
btemp1xorigin6_vector_NDE[14] = btemp1xorigin6_matrix_3_8;
btemp0xorigin6_vector_DE[15] = btemp0xorigin6_matrix_3_18_DE;
btemp0xorigin6_vector_NDE[15] = btemp0xorigin6_matrix_3_18_NDE;
btemp1xorigin6_vector_DE[15] = btemp1xorigin6_matrix_3_18_DE;
btemp1xorigin6_vector_NDE[15] = btemp1xorigin6_matrix_3_18_NDE;
btemp0xorigin6_vector_DE[16] = btemp0xorigin6_matrix_3_20_DE;
btemp0xorigin6_vector_NDE[16] = btemp0xorigin6_matrix_3_20_NDE;
btemp1xorigin6_vector_DE[16] = btemp1xorigin6_matrix_3_20_DE;
btemp1xorigin6_vector_NDE[16] = btemp1xorigin6_matrix_3_20_NDE;

// origin effects
borigin1_vector_DE[1] = borigin1_matrix_1_2;
borigin1_vector_NDE[1] = borigin1_matrix_1_2;
borigin1_vector_DE[2] = borigin1_matrix_2_1;
borigin1_vector_NDE[2] = borigin1_matrix_2_1;
borigin1_vector_DE[3] = borigin1_matrix_2_3;
borigin1_vector_NDE[3] = borigin1_matrix_2_3;
borigin1_vector_DE[4] = borigin1_matrix_2_10_DE;
borigin1_vector_NDE[4] = borigin1_matrix_2_10_NDE;
borigin1_vector_NDE[5] = borigin1_matrix_2_10_NDE;
borigin1_vector_DE[6] = borigin1_matrix_2_12_DE;
borigin1_vector_NDE[6] = borigin1_matrix_2_12_NDE;
borigin1_vector_DE[7] = borigin1_matrix_2_14_DE;
borigin1_vector_NDE[7] = borigin1_matrix_2_14_NDE;
borigin1_vector_NDE[8] = borigin1_matrix_2_14_NDE;
borigin1_vector_DE[9] = borigin1_matrix_2_16_DE;
borigin1_vector_NDE[9] = borigin1_matrix_2_16_NDE;
borigin1_vector_NDE[10] = borigin1_matrix_2_16_NDE;
borigin1_vector_DE[11] = borigin1_matrix_2_39;
borigin1_vector_NDE[11] = borigin1_matrix_2_39;
borigin1_vector_DE[12] = borigin1_matrix_3_2;
borigin1_vector_NDE[12] = borigin1_matrix_3_2;
borigin1_vector_DE[13] = borigin1_matrix_3_4;
borigin1_vector_NDE[13] = borigin1_matrix_3_4;
borigin1_vector_DE[14] = borigin1_matrix_3_8;
borigin1_vector_NDE[14] = borigin1_matrix_3_8;
borigin1_vector_DE[15] = borigin1_matrix_3_18_DE;
borigin1_vector_NDE[15] = borigin1_matrix_3_18_NDE;
borigin1_vector_DE[16] = borigin1_matrix_3_20_DE;
borigin1_vector_NDE[16] = borigin1_matrix_3_20_NDE;
borigin1_vector_DE[17] = borigin1_matrix_4_3;
borigin1_vector_NDE[17] = borigin1_matrix_4_3;
borigin1_vector_DE[27] = borigin1_matrix_8_3;
borigin1_vector_NDE[27] = borigin1_matrix_8_3;
borigin1_vector_DE[36] = borigin1_matrix_10_2;
borigin1_vector_NDE[36] = borigin1_matrix_10_2;
borigin1_vector_NDE[37] = borigin1_matrix_10_2;
borigin1_vector_DE[38] = borigin1_matrix_12_2;
borigin1_vector_NDE[38] = borigin1_matrix_12_2;
borigin1_vector_DE[39] = borigin1_matrix_14_2;
borigin1_vector_NDE[39] = borigin1_matrix_14_2;
borigin1_vector_NDE[40] = borigin1_matrix_14_2;
borigin1_vector_DE[41] = borigin1_matrix_16_2;
borigin1_vector_NDE[41] = borigin1_matrix_16_2;
borigin1_vector_NDE[42] = borigin1_matrix_16_2;
borigin1_vector_DE[43] = borigin1_matrix_18_3;
borigin1_vector_NDE[43] = borigin1_matrix_18_3;
borigin1_vector_DE[44] = borigin1_matrix_20_3;
borigin1_vector_NDE[44] = borigin1_matrix_20_3;
borigin1_vector_DE[54] = borigin1_matrix_39_2;
borigin1_vector_NDE[54] = borigin1_matrix_39_2;

borigin2_vector_DE[1] = borigin2_matrix_1_2;
borigin2_vector_NDE[1] = borigin2_matrix_1_2;
borigin2_vector_DE[2] = borigin2_matrix_2_1;
borigin2_vector_NDE[2] = borigin2_matrix_2_1;
borigin2_vector_DE[3] = borigin2_matrix_2_3;
borigin2_vector_NDE[3] = borigin2_matrix_2_3;
borigin2_vector_DE[4] = borigin2_matrix_2_10_DE;
borigin2_vector_NDE[4] = borigin2_matrix_2_10_NDE;
borigin2_vector_NDE[5] = borigin2_matrix_2_10_NDE;
borigin2_vector_DE[6] = borigin2_matrix_2_12_DE;
borigin2_vector_NDE[6] = borigin2_matrix_2_12_NDE;
borigin2_vector_DE[7] = borigin2_matrix_2_14_DE;
borigin2_vector_NDE[7] = borigin2_matrix_2_14_NDE;
borigin2_vector_NDE[8] = borigin2_matrix_2_14_NDE;
borigin2_vector_DE[9] = borigin2_matrix_2_16_DE;
borigin2_vector_NDE[9] = borigin2_matrix_2_16_NDE;
borigin2_vector_NDE[10] = borigin2_matrix_2_16_NDE;
borigin2_vector_DE[11] = borigin2_matrix_2_39;
borigin2_vector_NDE[11] = borigin2_matrix_2_39;
borigin2_vector_DE[12] = borigin2_matrix_3_2;
borigin2_vector_NDE[12] = borigin2_matrix_3_2;
borigin2_vector_DE[13] = borigin2_matrix_3_4;
borigin2_vector_NDE[13] = borigin2_matrix_3_4;
borigin2_vector_DE[14] = borigin2_matrix_3_8;
borigin2_vector_NDE[14] = borigin2_matrix_3_8;
borigin2_vector_DE[15] = borigin2_matrix_3_18_DE;
borigin2_vector_NDE[15] = borigin2_matrix_3_18_NDE;
borigin2_vector_DE[16] = borigin2_matrix_3_20_DE;
borigin2_vector_NDE[16] = borigin2_matrix_3_20_NDE;
borigin2_vector_DE[17] = borigin2_matrix_4_3;
borigin2_vector_NDE[17] = borigin2_matrix_4_3;
borigin2_vector_DE[27] = borigin2_matrix_8_3;
borigin2_vector_NDE[27] = borigin2_matrix_8_3;
borigin2_vector_DE[36] = borigin2_matrix_10_2;
borigin2_vector_NDE[36] = borigin2_matrix_10_2;
borigin2_vector_NDE[37] = borigin2_matrix_10_2;
borigin2_vector_DE[38] = borigin2_matrix_12_2;
borigin2_vector_NDE[38] = borigin2_matrix_12_2;
borigin2_vector_DE[39] = borigin2_matrix_14_2;
borigin2_vector_NDE[39] = borigin2_matrix_14_2;
borigin2_vector_NDE[40] = borigin2_matrix_14_2;
borigin2_vector_DE[41] = borigin2_matrix_16_2;
borigin2_vector_NDE[41] = borigin2_matrix_16_2;
borigin2_vector_NDE[42] = borigin2_matrix_16_2;
borigin2_vector_DE[43] = borigin2_matrix_18_3;
borigin2_vector_NDE[43] = borigin2_matrix_18_3;
borigin2_vector_DE[44] = borigin2_matrix_20_3;
borigin2_vector_NDE[44] = borigin2_matrix_20_3;
borigin2_vector_DE[54] = borigin2_matrix_39_2;
borigin2_vector_NDE[54] = borigin2_matrix_39_2;

borigin3_vector_DE[1] = borigin3_matrix_1_2;
borigin3_vector_NDE[1] = borigin3_matrix_1_2;
borigin3_vector_DE[2] = borigin3_matrix_2_1;
borigin3_vector_NDE[2] = borigin3_matrix_2_1;
borigin3_vector_DE[3] = borigin3_matrix_2_3;
borigin3_vector_NDE[3] = borigin3_matrix_2_3;
borigin3_vector_DE[4] = borigin3_matrix_2_10_DE;
borigin3_vector_NDE[4] = borigin3_matrix_2_10_NDE;
borigin3_vector_NDE[5] = borigin3_matrix_2_10_NDE;
borigin3_vector_DE[6] = borigin3_matrix_2_12_DE;
borigin3_vector_NDE[6] = borigin3_matrix_2_12_NDE;
borigin3_vector_DE[7] = borigin3_matrix_2_14_DE;
borigin3_vector_NDE[7] = borigin3_matrix_2_14_NDE;
borigin3_vector_NDE[8] = borigin3_matrix_2_14_NDE;
borigin3_vector_DE[9] = borigin3_matrix_2_16_DE;
borigin3_vector_NDE[9] = borigin3_matrix_2_16_NDE;
borigin3_vector_NDE[10] = borigin3_matrix_2_16_NDE;
borigin3_vector_DE[11] = borigin3_matrix_2_39;
borigin3_vector_NDE[11] = borigin3_matrix_2_39;
borigin3_vector_DE[12] = borigin3_matrix_3_2;
borigin3_vector_NDE[12] = borigin3_matrix_3_2;
borigin3_vector_DE[13] = borigin3_matrix_3_4;
borigin3_vector_NDE[13] = borigin3_matrix_3_4;
borigin3_vector_DE[14] = borigin3_matrix_3_8;
borigin3_vector_NDE[14] = borigin3_matrix_3_8;
borigin3_vector_DE[15] = borigin3_matrix_3_18_DE;
borigin3_vector_NDE[15] = borigin3_matrix_3_18_NDE;
borigin3_vector_DE[16] = borigin3_matrix_3_20_DE;
borigin3_vector_NDE[16] = borigin3_matrix_3_20_NDE;
borigin3_vector_DE[17] = borigin3_matrix_4_3;
borigin3_vector_NDE[17] = borigin3_matrix_4_3;
borigin3_vector_DE[27] = borigin3_matrix_8_3;
borigin3_vector_NDE[27] = borigin3_matrix_8_3;
borigin3_vector_DE[36] = borigin3_matrix_10_2;
borigin3_vector_NDE[36] = borigin3_matrix_10_2;
borigin3_vector_NDE[37] = borigin3_matrix_10_2;
borigin3_vector_DE[38] = borigin3_matrix_12_2;
borigin3_vector_NDE[38] = borigin3_matrix_12_2;
borigin3_vector_DE[39] = borigin3_matrix_14_2;
borigin3_vector_NDE[39] = borigin3_matrix_14_2;
borigin3_vector_NDE[40] = borigin3_matrix_14_2;
borigin3_vector_DE[41] = borigin3_matrix_16_2;
borigin3_vector_NDE[41] = borigin3_matrix_16_2;
borigin3_vector_NDE[42] = borigin3_matrix_16_2;
borigin3_vector_DE[43] = borigin3_matrix_18_3;
borigin3_vector_NDE[43] = borigin3_matrix_18_3;
borigin3_vector_DE[44] = borigin3_matrix_20_3;
borigin3_vector_NDE[44] = borigin3_matrix_20_3;
borigin3_vector_DE[54] = borigin3_matrix_39_2;
borigin3_vector_NDE[54] = borigin3_matrix_39_2;

borigin4_vector_DE[1] = borigin4_matrix_1_2;
borigin4_vector_NDE[1] = borigin4_matrix_1_2;
borigin4_vector_DE[2] = borigin4_matrix_2_1;
borigin4_vector_NDE[2] = borigin4_matrix_2_1;
borigin4_vector_DE[3] = borigin4_matrix_2_3;
borigin4_vector_NDE[3] = borigin4_matrix_2_3;
borigin4_vector_DE[4] = borigin4_matrix_2_10_DE;
borigin4_vector_NDE[4] = borigin4_matrix_2_10_NDE;
borigin4_vector_NDE[5] = borigin4_matrix_2_10_NDE;
borigin4_vector_DE[6] = borigin4_matrix_2_12_DE;
borigin4_vector_NDE[6] = borigin4_matrix_2_12_NDE;
borigin4_vector_DE[7] = borigin4_matrix_2_14_DE;
borigin4_vector_NDE[7] = borigin4_matrix_2_14_NDE;
borigin4_vector_NDE[8] = borigin4_matrix_2_14_NDE;
borigin4_vector_DE[9] = borigin4_matrix_2_16_DE;
borigin4_vector_NDE[9] = borigin4_matrix_2_16_NDE;
borigin4_vector_NDE[10] = borigin4_matrix_2_16_NDE;
borigin4_vector_DE[11] = borigin4_matrix_2_39;
borigin4_vector_NDE[11] = borigin4_matrix_2_39;
borigin4_vector_DE[12] = borigin4_matrix_3_2;
borigin4_vector_NDE[12] = borigin4_matrix_3_2;
borigin4_vector_DE[13] = borigin4_matrix_3_4;
borigin4_vector_NDE[13] = borigin4_matrix_3_4;
borigin4_vector_DE[14] = borigin4_matrix_3_8;
borigin4_vector_NDE[14] = borigin4_matrix_3_8;
borigin4_vector_DE[15] = borigin4_matrix_3_18_DE;
borigin4_vector_NDE[15] = borigin4_matrix_3_18_NDE;
borigin4_vector_DE[16] = borigin4_matrix_3_20_DE;
borigin4_vector_NDE[16] = borigin4_matrix_3_20_NDE;
borigin4_vector_DE[17] = borigin4_matrix_4_3;
borigin4_vector_NDE[17] = borigin4_matrix_4_3;
borigin4_vector_DE[27] = borigin4_matrix_8_3;
borigin4_vector_NDE[27] = borigin4_matrix_8_3;
borigin4_vector_DE[36] = borigin4_matrix_10_2;
borigin4_vector_NDE[36] = borigin4_matrix_10_2;
borigin4_vector_NDE[37] = borigin4_matrix_10_2;
borigin4_vector_DE[38] = borigin4_matrix_12_2;
borigin4_vector_NDE[38] = borigin4_matrix_12_2;
borigin4_vector_DE[39] = borigin4_matrix_14_2;
borigin4_vector_NDE[39] = borigin4_matrix_14_2;
borigin4_vector_NDE[40] = borigin4_matrix_14_2;
borigin4_vector_DE[41] = borigin4_matrix_16_2;
borigin4_vector_NDE[41] = borigin4_matrix_16_2;
borigin4_vector_NDE[42] = borigin4_matrix_16_2;
borigin4_vector_DE[43] = borigin4_matrix_18_3;
borigin4_vector_NDE[43] = borigin4_matrix_18_3;
borigin4_vector_DE[44] = borigin4_matrix_20_3;
borigin4_vector_NDE[44] = borigin4_matrix_20_3;
borigin4_vector_DE[54] = borigin4_matrix_39_2;
borigin4_vector_NDE[54] = borigin4_matrix_39_2;

borigin5_vector_DE[1] = borigin5_matrix_1_2;
borigin5_vector_NDE[1] = borigin5_matrix_1_2;
borigin5_vector_DE[2] = borigin5_matrix_2_1;
borigin5_vector_NDE[2] = borigin5_matrix_2_1;
borigin5_vector_DE[3] = borigin5_matrix_2_3;
borigin5_vector_NDE[3] = borigin5_matrix_2_3;
borigin5_vector_DE[4] = borigin5_matrix_2_10_DE;
borigin5_vector_NDE[4] = borigin5_matrix_2_10_NDE;
borigin5_vector_NDE[5] = borigin5_matrix_2_10_NDE;
borigin5_vector_DE[6] = borigin5_matrix_2_12_DE;
borigin5_vector_NDE[6] = borigin5_matrix_2_12_NDE;
borigin5_vector_DE[7] = borigin5_matrix_2_14_DE;
borigin5_vector_NDE[7] = borigin5_matrix_2_14_NDE;
borigin5_vector_NDE[8] = borigin5_matrix_2_14_NDE;
borigin5_vector_DE[9] = borigin5_matrix_2_16_DE;
borigin5_vector_NDE[9] = borigin5_matrix_2_16_NDE;
borigin5_vector_NDE[10] = borigin5_matrix_2_16_NDE;
borigin5_vector_DE[11] = borigin5_matrix_2_39;
borigin5_vector_NDE[11] = borigin5_matrix_2_39;
borigin5_vector_DE[12] = borigin5_matrix_3_2;
borigin5_vector_NDE[12] = borigin5_matrix_3_2;
borigin5_vector_DE[13] = borigin5_matrix_3_4;
borigin5_vector_NDE[13] = borigin5_matrix_3_4;
borigin5_vector_DE[14] = borigin5_matrix_3_8;
borigin5_vector_NDE[14] = borigin5_matrix_3_8;
borigin5_vector_DE[15] = borigin5_matrix_3_18_DE;
borigin5_vector_NDE[15] = borigin5_matrix_3_18_NDE;
borigin5_vector_DE[16] = borigin5_matrix_3_20_DE;
borigin5_vector_NDE[16] = borigin5_matrix_3_20_NDE;
borigin5_vector_DE[17] = borigin5_matrix_4_3;
borigin5_vector_NDE[17] = borigin5_matrix_4_3;
borigin5_vector_DE[27] = borigin5_matrix_8_3;
borigin5_vector_NDE[27] = borigin5_matrix_8_3;
borigin5_vector_DE[36] = borigin5_matrix_10_2;
borigin5_vector_NDE[36] = borigin5_matrix_10_2;
borigin5_vector_NDE[37] = borigin5_matrix_10_2;
borigin5_vector_DE[38] = borigin5_matrix_12_2;
borigin5_vector_NDE[38] = borigin5_matrix_12_2;
borigin5_vector_DE[39] = borigin5_matrix_14_2;
borigin5_vector_NDE[39] = borigin5_matrix_14_2;
borigin5_vector_NDE[40] = borigin5_matrix_14_2;
borigin5_vector_DE[41] = borigin5_matrix_16_2;
borigin5_vector_NDE[41] = borigin5_matrix_16_2;
borigin5_vector_NDE[42] = borigin5_matrix_16_2;
borigin5_vector_DE[43] = borigin5_matrix_18_3;
borigin5_vector_NDE[43] = borigin5_matrix_18_3;
borigin5_vector_DE[44] = borigin5_matrix_20_3;
borigin5_vector_NDE[44] = borigin5_matrix_20_3;
borigin5_vector_DE[54] = borigin5_matrix_39_2;
borigin5_vector_NDE[54] = borigin5_matrix_39_2;

borigin6_vector_DE[1] = borigin6_matrix_1_2;
borigin6_vector_NDE[1] = borigin6_matrix_1_2;
borigin6_vector_DE[2] = borigin6_matrix_2_1;
borigin6_vector_NDE[2] = borigin6_matrix_2_1;
borigin6_vector_DE[3] = borigin6_matrix_2_3;
borigin6_vector_NDE[3] = borigin6_matrix_2_3;
borigin6_vector_DE[4] = borigin6_matrix_2_10_DE;
borigin6_vector_NDE[4] = borigin6_matrix_2_10_NDE;
borigin6_vector_NDE[5] = borigin6_matrix_2_10_NDE;
borigin6_vector_DE[6] = borigin6_matrix_2_12_DE;
borigin6_vector_NDE[6] = borigin6_matrix_2_12_NDE;
borigin6_vector_DE[7] = borigin6_matrix_2_14_DE;
borigin6_vector_NDE[7] = borigin6_matrix_2_14_NDE;
borigin6_vector_NDE[8] = borigin6_matrix_2_14_NDE;
borigin6_vector_DE[9] = borigin6_matrix_2_16_DE;
borigin6_vector_NDE[9] = borigin6_matrix_2_16_NDE;
borigin6_vector_NDE[10] = borigin6_matrix_2_16_NDE;
borigin6_vector_DE[11] = borigin6_matrix_2_39;
borigin6_vector_NDE[11] = borigin6_matrix_2_39;
borigin6_vector_DE[12] = borigin6_matrix_3_2;
borigin6_vector_NDE[12] = borigin6_matrix_3_2;
borigin6_vector_DE[13] = borigin6_matrix_3_4;
borigin6_vector_NDE[13] = borigin6_matrix_3_4;
borigin6_vector_DE[14] = borigin6_matrix_3_8;
borigin6_vector_NDE[14] = borigin6_matrix_3_8;
borigin6_vector_DE[15] = borigin6_matrix_3_18_DE;
borigin6_vector_NDE[15] = borigin6_matrix_3_18_NDE;
borigin6_vector_DE[16] = borigin6_matrix_3_20_DE;
borigin6_vector_NDE[16] = borigin6_matrix_3_20_NDE;
borigin6_vector_DE[17] = borigin6_matrix_4_3;
borigin6_vector_NDE[17] = borigin6_matrix_4_3;
borigin6_vector_DE[27] = borigin6_matrix_8_3;
borigin6_vector_NDE[27] = borigin6_matrix_8_3;
borigin6_vector_DE[36] = borigin6_matrix_10_2;
borigin6_vector_NDE[36] = borigin6_matrix_10_2;
borigin6_vector_NDE[37] = borigin6_matrix_10_2;
borigin6_vector_DE[38] = borigin6_matrix_12_2;
borigin6_vector_NDE[38] = borigin6_matrix_12_2;
borigin6_vector_DE[39] = borigin6_matrix_14_2;
borigin6_vector_NDE[39] = borigin6_matrix_14_2;
borigin6_vector_NDE[40] = borigin6_matrix_14_2;
borigin6_vector_DE[41] = borigin6_matrix_16_2;
borigin6_vector_NDE[41] = borigin6_matrix_16_2;
borigin6_vector_NDE[42] = borigin6_matrix_16_2;
borigin6_vector_DE[43] = borigin6_matrix_18_3;
borigin6_vector_NDE[43] = borigin6_matrix_18_3;
borigin6_vector_DE[44] = borigin6_matrix_20_3;
borigin6_vector_NDE[44] = borigin6_matrix_20_3;
borigin6_vector_DE[54] = borigin6_matrix_39_2;
borigin6_vector_NDE[54] = borigin6_matrix_39_2;

// detection efficiency - create a vector that stores all parameters
vector[33] det_eff_param_vector; // this is length 33 because that's how many det eff params we have

// populate the vector
det_eff_param_vector[1] = asotin_alpha1;
det_eff_param_vector[2] = asotin_alpha2;
det_eff_param_vector[3] = deschutes_alpha1;
det_eff_param_vector[4] = entiat_alpha1;
det_eff_param_vector[5] = fifteenmile_alpha1;
det_eff_param_vector[6] = imnaha_alpha1;
det_eff_param_vector[7] = john_day_alpha1;
det_eff_param_vector[8] = methow_alpha1;
det_eff_param_vector[9] = methow_alpha2;
det_eff_param_vector[10] = okanogan_alpha1;
det_eff_param_vector[11] = tucannon_alpha1;
det_eff_param_vector[12] = tucannon_alpha2;
det_eff_param_vector[13] = umatilla_alpha1;
det_eff_param_vector[14] = umatilla_alpha2;
det_eff_param_vector[15] = walla_walla_alpha1;
det_eff_param_vector[16] = walla_walla_alpha2;
det_eff_param_vector[17] = walla_walla_alpha3;
det_eff_param_vector[18] = walla_walla_alpha4;
det_eff_param_vector[19] = wenatchee_alpha1;
det_eff_param_vector[20] = yakima_alpha1;

// 14 terms for discharge relationship, one for each tributary
det_eff_param_vector[21] = asotin_beta;
det_eff_param_vector[22] = deschutes_beta;
det_eff_param_vector[23] = entiat_beta;
// det_eff_param_vector[24] = fifteenmile_beta;
// // fix these to zero
det_eff_param_vector[24] = 0;
// det_eff_param_vector[25] = imnaha_beta;
det_eff_param_vector[25] = 0;
det_eff_param_vector[26] = john_day_beta;
det_eff_param_vector[27] = methow_beta;
det_eff_param_vector[28] = okanogan_beta;
det_eff_param_vector[29] = tucannon_beta;
det_eff_param_vector[30] = umatilla_beta;
det_eff_param_vector[31] = walla_walla_beta;
det_eff_param_vector[32] = wenatchee_beta;
det_eff_param_vector[33] = yakima_beta;

// Calculate each individual byear parameter per movement as a derived parameter
// this is because in actuality, each byear parameter for a movement is 
// calculated as a sigma

}


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

  // Set the priors for each of the non-zero elements of the b0 matrix
b0_matrix_4_5 ~ normal(0,5);
b0_matrix_5_4 ~ normal(0,5);
b0_matrix_5_6 ~ normal(0,5);
b0_matrix_5_22_DE ~ normal(0,5);
b0_matrix_5_22_NDE ~ normal(0,5);
b0_matrix_6_5 ~ normal(0,5);
b0_matrix_6_7 ~ normal(0,5);
b0_matrix_7_6 ~ normal(0,5);
b0_matrix_7_26_DE ~ normal(0,5);
b0_matrix_7_26_NDE ~ normal(0,5);
b0_matrix_7_28_DE ~ normal(0,5);
b0_matrix_7_28_NDE ~ normal(0,5);
b0_matrix_8_9 ~ normal(0,5);
b0_matrix_8_30_DE ~ normal(0,5);
b0_matrix_8_30_NDE ~ normal(0,5);
b0_matrix_9_8 ~ normal(0,5);
b0_matrix_9_32_DE ~ normal(0,5);
b0_matrix_9_32_NDE ~ normal(0,5);
b0_matrix_9_34 ~ normal(0,5);
b0_matrix_9_35 ~ normal(0,5);
b0_matrix_9_36 ~ normal(0,5);
b0_matrix_9_37_DE ~ normal(0,5);
b0_matrix_9_37_NDE ~ normal(0,5);
b0_matrix_22_5 ~ normal(0,5);
b0_matrix_26_7 ~ normal(0,5);
b0_matrix_28_7 ~ normal(0,5);
b0_matrix_30_8 ~ normal(0,5);
b0_matrix_32_9 ~ normal(0,5);
b0_matrix_34_9 ~ normal(0,5);
b0_matrix_35_9 ~ normal(0,5);
b0_matrix_36_9 ~ normal(0,5);
b0_matrix_37_9 ~ normal(0,5);

// spill priors
bspillwindow_matrix_2_1 ~ normal(0,5);
bspillwindow_matrix_3_2 ~ normal(0,5);
bspillwindow_matrix_4_3 ~ normal(0,5);
bspillwindow_matrix_5_4 ~ normal(0,5);
bspillwindow_matrix_6_5 ~ normal(0,5);
bspillwindow_matrix_7_6 ~ normal(0,5);
bspillwindow_matrix_8_3 ~ normal(0,5);
bspillwindow_matrix_9_8 ~ normal(0,5);

bwinterspill_matrix_3_2 ~ normal(0,5);
bwinterspill_matrix_4_3 ~ normal(0,5);
bwinterspill_matrix_5_4 ~ normal(0,5);
bwinterspill_matrix_6_5 ~ normal(0,5);
bwinterspill_matrix_7_6 ~ normal(0,5);
bwinterspill_matrix_8_3 ~ normal(0,5);
bwinterspill_matrix_9_8 ~ normal(0,5);

// temp priors
btemp0_matrix_4_5 ~ normal(0,5);
btemp1_matrix_4_5 ~ normal(0,5);
btemp0_matrix_5_6 ~ normal(0,5);
btemp1_matrix_5_6 ~ normal(0,5);
btemp0_matrix_5_22_DE ~ normal(0,5);
btemp0_matrix_5_22_NDE ~ normal(0,5);
btemp1_matrix_5_22_DE ~ normal(0,5);
btemp1_matrix_5_22_NDE ~ normal(0,5);
btemp0_matrix_6_7 ~ normal(0,5);
btemp1_matrix_6_7 ~ normal(0,5);
btemp0_matrix_7_26_DE ~ normal(0,5);
btemp0_matrix_7_26_NDE ~ normal(0,5);
btemp1_matrix_7_26_DE ~ normal(0,5);
btemp1_matrix_7_26_NDE ~ normal(0,5);
btemp0_matrix_7_28_DE ~ normal(0,5);
btemp0_matrix_7_28_NDE ~ normal(0,5);
btemp1_matrix_7_28_DE ~ normal(0,5);
btemp1_matrix_7_28_NDE ~ normal(0,5);
btemp0_matrix_8_9 ~ normal(0,5);
btemp1_matrix_8_9 ~ normal(0,5);
btemp0_matrix_8_30_DE ~ normal(0,5);
btemp0_matrix_8_30_NDE ~ normal(0,5);
btemp1_matrix_8_30_DE ~ normal(0,5);
btemp1_matrix_8_30_NDE ~ normal(0,5);
btemp0_matrix_9_32_DE ~ normal(0,5);
btemp0_matrix_9_32_NDE ~ normal(0,5);
btemp1_matrix_9_32_DE ~ normal(0,5);
btemp1_matrix_9_32_NDE ~ normal(0,5);
btemp0_matrix_9_34 ~ normal(0,5);
btemp1_matrix_9_34 ~ normal(0,5);
btemp0_matrix_9_35 ~ normal(0,5);
btemp1_matrix_9_35 ~ normal(0,5);
btemp0_matrix_9_36 ~ normal(0,5);
btemp1_matrix_9_36 ~ normal(0,5);
btemp0_matrix_9_37_DE ~ normal(0,5);
btemp0_matrix_9_37_NDE ~ normal(0,5);
btemp1_matrix_9_37_DE ~ normal(0,5);
btemp1_matrix_9_37_NDE ~ normal(0,5);

btemp0xorigin1_matrix_1_2 ~ normal(0,5);
btemp1xorigin1_matrix_1_2 ~ normal(0,5);
btemp0xorigin1_matrix_2_3 ~ normal(0,5);
btemp1xorigin1_matrix_2_3 ~ normal(0,5);
btemp0xorigin1_matrix_2_10_DE ~ normal(0,5);
btemp0xorigin1_matrix_2_10_NDE ~ normal(0,5);
btemp1xorigin1_matrix_2_10_DE ~ normal(0,5);
btemp1xorigin1_matrix_2_10_NDE ~ normal(0,5);
btemp0xorigin1_matrix_2_12_DE ~ normal(0,5);
btemp0xorigin1_matrix_2_12_NDE ~ normal(0,5);
btemp1xorigin1_matrix_2_12_DE ~ normal(0,5);
btemp1xorigin1_matrix_2_12_NDE ~ normal(0,5);
btemp0xorigin1_matrix_2_14_DE ~ normal(0,5);
btemp0xorigin1_matrix_2_14_NDE ~ normal(0,5);
btemp1xorigin1_matrix_2_14_DE ~ normal(0,5);
btemp1xorigin1_matrix_2_14_NDE ~ normal(0,5);
btemp0xorigin1_matrix_2_16_DE ~ normal(0,5);
btemp0xorigin1_matrix_2_16_NDE ~ normal(0,5);
btemp1xorigin1_matrix_2_16_DE ~ normal(0,5);
btemp1xorigin1_matrix_2_16_NDE ~ normal(0,5);
btemp0xorigin1_matrix_2_39 ~ normal(0,5);
btemp1xorigin1_matrix_2_39 ~ normal(0,5);
btemp0xorigin1_matrix_3_4 ~ normal(0,5);
btemp1xorigin1_matrix_3_4 ~ normal(0,5);
btemp0xorigin1_matrix_3_8 ~ normal(0,5);
btemp1xorigin1_matrix_3_8 ~ normal(0,5);
btemp0xorigin1_matrix_3_18_DE ~ normal(0,5);
btemp0xorigin1_matrix_3_18_NDE ~ normal(0,5);
btemp1xorigin1_matrix_3_18_DE ~ normal(0,5);
btemp1xorigin1_matrix_3_18_NDE ~ normal(0,5);
btemp0xorigin1_matrix_3_20_DE ~ normal(0,5);
btemp0xorigin1_matrix_3_20_NDE ~ normal(0,5);
btemp1xorigin1_matrix_3_20_DE ~ normal(0,5);
btemp1xorigin1_matrix_3_20_NDE ~ normal(0,5);

btemp0xorigin2_matrix_1_2 ~ normal(0,5);
btemp1xorigin2_matrix_1_2 ~ normal(0,5);
btemp0xorigin2_matrix_2_3 ~ normal(0,5);
btemp1xorigin2_matrix_2_3 ~ normal(0,5);
btemp0xorigin2_matrix_2_10_DE ~ normal(0,5);
btemp0xorigin2_matrix_2_10_NDE ~ normal(0,5);
btemp1xorigin2_matrix_2_10_DE ~ normal(0,5);
btemp1xorigin2_matrix_2_10_NDE ~ normal(0,5);
btemp0xorigin2_matrix_2_12_DE ~ normal(0,5);
btemp0xorigin2_matrix_2_12_NDE ~ normal(0,5);
btemp1xorigin2_matrix_2_12_DE ~ normal(0,5);
btemp1xorigin2_matrix_2_12_NDE ~ normal(0,5);
btemp0xorigin2_matrix_2_14_DE ~ normal(0,5);
btemp0xorigin2_matrix_2_14_NDE ~ normal(0,5);
btemp1xorigin2_matrix_2_14_DE ~ normal(0,5);
btemp1xorigin2_matrix_2_14_NDE ~ normal(0,5);
btemp0xorigin2_matrix_2_16_DE ~ normal(0,5);
btemp0xorigin2_matrix_2_16_NDE ~ normal(0,5);
btemp1xorigin2_matrix_2_16_DE ~ normal(0,5);
btemp1xorigin2_matrix_2_16_NDE ~ normal(0,5);
btemp0xorigin2_matrix_2_39 ~ normal(0,5);
btemp1xorigin2_matrix_2_39 ~ normal(0,5);
btemp0xorigin2_matrix_3_4 ~ normal(0,5);
btemp1xorigin2_matrix_3_4 ~ normal(0,5);
btemp0xorigin2_matrix_3_8 ~ normal(0,5);
btemp1xorigin2_matrix_3_8 ~ normal(0,5);
btemp0xorigin2_matrix_3_18_DE ~ normal(0,5);
btemp0xorigin2_matrix_3_18_NDE ~ normal(0,5);
btemp1xorigin2_matrix_3_18_DE ~ normal(0,5);
btemp1xorigin2_matrix_3_18_NDE ~ normal(0,5);
btemp0xorigin2_matrix_3_20_DE ~ normal(0,5);
btemp0xorigin2_matrix_3_20_NDE ~ normal(0,5);
btemp1xorigin2_matrix_3_20_DE ~ normal(0,5);
btemp1xorigin2_matrix_3_20_NDE ~ normal(0,5);

btemp0xorigin3_matrix_1_2 ~ normal(0,5);
btemp1xorigin3_matrix_1_2 ~ normal(0,5);
btemp0xorigin3_matrix_2_3 ~ normal(0,5);
btemp1xorigin3_matrix_2_3 ~ normal(0,5);
btemp0xorigin3_matrix_2_10_DE ~ normal(0,5);
btemp0xorigin3_matrix_2_10_NDE ~ normal(0,5);
btemp1xorigin3_matrix_2_10_DE ~ normal(0,5);
btemp1xorigin3_matrix_2_10_NDE ~ normal(0,5);
btemp0xorigin3_matrix_2_12_DE ~ normal(0,5);
btemp0xorigin3_matrix_2_12_NDE ~ normal(0,5);
btemp1xorigin3_matrix_2_12_DE ~ normal(0,5);
btemp1xorigin3_matrix_2_12_NDE ~ normal(0,5);
btemp0xorigin3_matrix_2_14_DE ~ normal(0,5);
btemp0xorigin3_matrix_2_14_NDE ~ normal(0,5);
btemp1xorigin3_matrix_2_14_DE ~ normal(0,5);
btemp1xorigin3_matrix_2_14_NDE ~ normal(0,5);
btemp0xorigin3_matrix_2_16_DE ~ normal(0,5);
btemp0xorigin3_matrix_2_16_NDE ~ normal(0,5);
btemp1xorigin3_matrix_2_16_DE ~ normal(0,5);
btemp1xorigin3_matrix_2_16_NDE ~ normal(0,5);
btemp0xorigin3_matrix_2_39 ~ normal(0,5);
btemp1xorigin3_matrix_2_39 ~ normal(0,5);
btemp0xorigin3_matrix_3_4 ~ normal(0,5);
btemp1xorigin3_matrix_3_4 ~ normal(0,5);
btemp0xorigin3_matrix_3_8 ~ normal(0,5);
btemp1xorigin3_matrix_3_8 ~ normal(0,5);
btemp0xorigin3_matrix_3_18_DE ~ normal(0,5);
btemp0xorigin3_matrix_3_18_NDE ~ normal(0,5);
btemp1xorigin3_matrix_3_18_DE ~ normal(0,5);
btemp1xorigin3_matrix_3_18_NDE ~ normal(0,5);
btemp0xorigin3_matrix_3_20_DE ~ normal(0,5);
btemp0xorigin3_matrix_3_20_NDE ~ normal(0,5);
btemp1xorigin3_matrix_3_20_DE ~ normal(0,5);
btemp1xorigin3_matrix_3_20_NDE ~ normal(0,5);

btemp0xorigin4_matrix_1_2 ~ normal(0,5);
btemp1xorigin4_matrix_1_2 ~ normal(0,5);
btemp0xorigin4_matrix_2_3 ~ normal(0,5);
btemp1xorigin4_matrix_2_3 ~ normal(0,5);
btemp0xorigin4_matrix_2_10_DE ~ normal(0,5);
btemp0xorigin4_matrix_2_10_NDE ~ normal(0,5);
btemp1xorigin4_matrix_2_10_DE ~ normal(0,5);
btemp1xorigin4_matrix_2_10_NDE ~ normal(0,5);
btemp0xorigin4_matrix_2_12_DE ~ normal(0,5);
btemp0xorigin4_matrix_2_12_NDE ~ normal(0,5);
btemp1xorigin4_matrix_2_12_DE ~ normal(0,5);
btemp1xorigin4_matrix_2_12_NDE ~ normal(0,5);
btemp0xorigin4_matrix_2_14_DE ~ normal(0,5);
btemp0xorigin4_matrix_2_14_NDE ~ normal(0,5);
btemp1xorigin4_matrix_2_14_DE ~ normal(0,5);
btemp1xorigin4_matrix_2_14_NDE ~ normal(0,5);
btemp0xorigin4_matrix_2_16_DE ~ normal(0,5);
btemp0xorigin4_matrix_2_16_NDE ~ normal(0,5);
btemp1xorigin4_matrix_2_16_DE ~ normal(0,5);
btemp1xorigin4_matrix_2_16_NDE ~ normal(0,5);
btemp0xorigin4_matrix_2_39 ~ normal(0,5);
btemp1xorigin4_matrix_2_39 ~ normal(0,5);
btemp0xorigin4_matrix_3_4 ~ normal(0,5);
btemp1xorigin4_matrix_3_4 ~ normal(0,5);
btemp0xorigin4_matrix_3_8 ~ normal(0,5);
btemp1xorigin4_matrix_3_8 ~ normal(0,5);
btemp0xorigin4_matrix_3_18_DE ~ normal(0,5);
btemp0xorigin4_matrix_3_18_NDE ~ normal(0,5);
btemp1xorigin4_matrix_3_18_DE ~ normal(0,5);
btemp1xorigin4_matrix_3_18_NDE ~ normal(0,5);
btemp0xorigin4_matrix_3_20_DE ~ normal(0,5);
btemp0xorigin4_matrix_3_20_NDE ~ normal(0,5);
btemp1xorigin4_matrix_3_20_DE ~ normal(0,5);
btemp1xorigin4_matrix_3_20_NDE ~ normal(0,5);

btemp0xorigin5_matrix_1_2 ~ normal(0,5);
btemp1xorigin5_matrix_1_2 ~ normal(0,5);
btemp0xorigin5_matrix_2_3 ~ normal(0,5);
btemp1xorigin5_matrix_2_3 ~ normal(0,5);
btemp0xorigin5_matrix_2_10_DE ~ normal(0,5);
btemp0xorigin5_matrix_2_10_NDE ~ normal(0,5);
btemp1xorigin5_matrix_2_10_DE ~ normal(0,5);
btemp1xorigin5_matrix_2_10_NDE ~ normal(0,5);
btemp0xorigin5_matrix_2_12_DE ~ normal(0,5);
btemp0xorigin5_matrix_2_12_NDE ~ normal(0,5);
btemp1xorigin5_matrix_2_12_DE ~ normal(0,5);
btemp1xorigin5_matrix_2_12_NDE ~ normal(0,5);
btemp0xorigin5_matrix_2_14_DE ~ normal(0,5);
btemp0xorigin5_matrix_2_14_NDE ~ normal(0,5);
btemp1xorigin5_matrix_2_14_DE ~ normal(0,5);
btemp1xorigin5_matrix_2_14_NDE ~ normal(0,5);
btemp0xorigin5_matrix_2_16_DE ~ normal(0,5);
btemp0xorigin5_matrix_2_16_NDE ~ normal(0,5);
btemp1xorigin5_matrix_2_16_DE ~ normal(0,5);
btemp1xorigin5_matrix_2_16_NDE ~ normal(0,5);
btemp0xorigin5_matrix_2_39 ~ normal(0,5);
btemp1xorigin5_matrix_2_39 ~ normal(0,5);
btemp0xorigin5_matrix_3_4 ~ normal(0,5);
btemp1xorigin5_matrix_3_4 ~ normal(0,5);
btemp0xorigin5_matrix_3_8 ~ normal(0,5);
btemp1xorigin5_matrix_3_8 ~ normal(0,5);
btemp0xorigin5_matrix_3_18_DE ~ normal(0,5);
btemp0xorigin5_matrix_3_18_NDE ~ normal(0,5);
btemp1xorigin5_matrix_3_18_DE ~ normal(0,5);
btemp1xorigin5_matrix_3_18_NDE ~ normal(0,5);
btemp0xorigin5_matrix_3_20_DE ~ normal(0,5);
btemp0xorigin5_matrix_3_20_NDE ~ normal(0,5);
btemp1xorigin5_matrix_3_20_DE ~ normal(0,5);
btemp1xorigin5_matrix_3_20_NDE ~ normal(0,5);

btemp0xorigin6_matrix_1_2 ~ normal(0,5);
btemp1xorigin6_matrix_1_2 ~ normal(0,5);
btemp0xorigin6_matrix_2_3 ~ normal(0,5);
btemp1xorigin6_matrix_2_3 ~ normal(0,5);
btemp0xorigin6_matrix_2_10_DE ~ normal(0,5);
btemp0xorigin6_matrix_2_10_NDE ~ normal(0,5);
btemp1xorigin6_matrix_2_10_DE ~ normal(0,5);
btemp1xorigin6_matrix_2_10_NDE ~ normal(0,5);
btemp0xorigin6_matrix_2_12_DE ~ normal(0,5);
btemp0xorigin6_matrix_2_12_NDE ~ normal(0,5);
btemp1xorigin6_matrix_2_12_DE ~ normal(0,5);
btemp1xorigin6_matrix_2_12_NDE ~ normal(0,5);
btemp0xorigin6_matrix_2_14_DE ~ normal(0,5);
btemp0xorigin6_matrix_2_14_NDE ~ normal(0,5);
btemp1xorigin6_matrix_2_14_DE ~ normal(0,5);
btemp1xorigin6_matrix_2_14_NDE ~ normal(0,5);
btemp0xorigin6_matrix_2_16_DE ~ normal(0,5);
btemp0xorigin6_matrix_2_16_NDE ~ normal(0,5);
btemp1xorigin6_matrix_2_16_DE ~ normal(0,5);
btemp1xorigin6_matrix_2_16_NDE ~ normal(0,5);
btemp0xorigin6_matrix_2_39 ~ normal(0,5);
btemp1xorigin6_matrix_2_39 ~ normal(0,5);
btemp0xorigin6_matrix_3_4 ~ normal(0,5);
btemp1xorigin6_matrix_3_4 ~ normal(0,5);
btemp0xorigin6_matrix_3_8 ~ normal(0,5);
btemp1xorigin6_matrix_3_8 ~ normal(0,5);
btemp0xorigin6_matrix_3_18_DE ~ normal(0,5);
btemp0xorigin6_matrix_3_18_NDE ~ normal(0,5);
btemp1xorigin6_matrix_3_18_DE ~ normal(0,5);
btemp1xorigin6_matrix_3_18_NDE ~ normal(0,5);
btemp0xorigin6_matrix_3_20_DE ~ normal(0,5);
btemp0xorigin6_matrix_3_20_NDE ~ normal(0,5);
btemp1xorigin6_matrix_3_20_DE ~ normal(0,5);
btemp1xorigin6_matrix_3_20_NDE ~ normal(0,5);

// priors for origin parameters
borigin1_matrix_1_2 ~ normal(0,5);
borigin1_matrix_2_1 ~ normal(0,5);
borigin1_matrix_2_3 ~ normal(0,5);
borigin1_matrix_2_10_DE ~ normal(0,5);
borigin1_matrix_2_10_NDE ~ normal(0,5);
borigin1_matrix_2_12_DE ~ normal(0,5);
borigin1_matrix_2_12_NDE ~ normal(0,5);
borigin1_matrix_2_14_DE ~ normal(0,5);
borigin1_matrix_2_14_NDE ~ normal(0,5);
borigin1_matrix_2_16_DE ~ normal(0,5);
borigin1_matrix_2_16_NDE ~ normal(0,5);
borigin1_matrix_2_39 ~ normal(0,5);
borigin1_matrix_3_2 ~ normal(0,5);
borigin1_matrix_3_4 ~ normal(0,5);
borigin1_matrix_3_8 ~ normal(0,5);
borigin1_matrix_3_18_DE ~ normal(0,5);
borigin1_matrix_3_18_NDE ~ normal(0,5);
borigin1_matrix_3_20_DE ~ normal(0,5);
borigin1_matrix_3_20_NDE ~ normal(0,5);
borigin1_matrix_4_3 ~ normal(0,5);
borigin1_matrix_8_3 ~ normal(0,5);
borigin1_matrix_10_2 ~ normal(0,5);
borigin1_matrix_12_2 ~ normal(0,5);
borigin1_matrix_14_2 ~ normal(0,5);
borigin1_matrix_16_2 ~ normal(0,5);
borigin1_matrix_18_3 ~ normal(0,5);
borigin1_matrix_20_3 ~ normal(0,5);
borigin1_matrix_39_2 ~ normal(0,5);

borigin2_matrix_1_2 ~ normal(0,5);
borigin2_matrix_2_1 ~ normal(0,5);
borigin2_matrix_2_3 ~ normal(0,5);
borigin2_matrix_2_10_DE ~ normal(0,5);
borigin2_matrix_2_10_NDE ~ normal(0,5);
borigin2_matrix_2_12_DE ~ normal(0,5);
borigin2_matrix_2_12_NDE ~ normal(0,5);
borigin2_matrix_2_14_DE ~ normal(0,5);
borigin2_matrix_2_14_NDE ~ normal(0,5);
borigin2_matrix_2_16_DE ~ normal(0,5);
borigin2_matrix_2_16_NDE ~ normal(0,5);
borigin2_matrix_2_39 ~ normal(0,5);
borigin2_matrix_3_2 ~ normal(0,5);
borigin2_matrix_3_4 ~ normal(0,5);
borigin2_matrix_3_8 ~ normal(0,5);
borigin2_matrix_3_18_DE ~ normal(0,5);
borigin2_matrix_3_18_NDE ~ normal(0,5);
borigin2_matrix_3_20_DE ~ normal(0,5);
borigin2_matrix_3_20_NDE ~ normal(0,5);
borigin2_matrix_4_3 ~ normal(0,5);
borigin2_matrix_8_3 ~ normal(0,5);
borigin2_matrix_10_2 ~ normal(0,5);
borigin2_matrix_12_2 ~ normal(0,5);
borigin2_matrix_14_2 ~ normal(0,5);
borigin2_matrix_16_2 ~ normal(0,5);
borigin2_matrix_18_3 ~ normal(0,5);
borigin2_matrix_20_3 ~ normal(0,5);
borigin2_matrix_39_2 ~ normal(0,5);

borigin3_matrix_1_2 ~ normal(0,5);
borigin3_matrix_2_1 ~ normal(0,5);
borigin3_matrix_2_3 ~ normal(0,5);
borigin3_matrix_2_10_DE ~ normal(0,5);
borigin3_matrix_2_10_NDE ~ normal(0,5);
borigin3_matrix_2_12_DE ~ normal(0,5);
borigin3_matrix_2_12_NDE ~ normal(0,5);
borigin3_matrix_2_14_DE ~ normal(0,5);
borigin3_matrix_2_14_NDE ~ normal(0,5);
borigin3_matrix_2_16_DE ~ normal(0,5);
borigin3_matrix_2_16_NDE ~ normal(0,5);
borigin3_matrix_2_39 ~ normal(0,5);
borigin3_matrix_3_2 ~ normal(0,5);
borigin3_matrix_3_4 ~ normal(0,5);
borigin3_matrix_3_8 ~ normal(0,5);
borigin3_matrix_3_18_DE ~ normal(0,5);
borigin3_matrix_3_18_NDE ~ normal(0,5);
borigin3_matrix_3_20_DE ~ normal(0,5);
borigin3_matrix_3_20_NDE ~ normal(0,5);
borigin3_matrix_4_3 ~ normal(0,5);
borigin3_matrix_8_3 ~ normal(0,5);
borigin3_matrix_10_2 ~ normal(0,5);
borigin3_matrix_12_2 ~ normal(0,5);
borigin3_matrix_14_2 ~ normal(0,5);
borigin3_matrix_16_2 ~ normal(0,5);
borigin3_matrix_18_3 ~ normal(0,5);
borigin3_matrix_20_3 ~ normal(0,5);
borigin3_matrix_39_2 ~ normal(0,5);

borigin4_matrix_1_2 ~ normal(0,5);
borigin4_matrix_2_1 ~ normal(0,5);
borigin4_matrix_2_3 ~ normal(0,5);
borigin4_matrix_2_10_DE ~ normal(0,5);
borigin4_matrix_2_10_NDE ~ normal(0,5);
borigin4_matrix_2_12_DE ~ normal(0,5);
borigin4_matrix_2_12_NDE ~ normal(0,5);
borigin4_matrix_2_14_DE ~ normal(0,5);
borigin4_matrix_2_14_NDE ~ normal(0,5);
borigin4_matrix_2_16_DE ~ normal(0,5);
borigin4_matrix_2_16_NDE ~ normal(0,5);
borigin4_matrix_2_39 ~ normal(0,5);
borigin4_matrix_3_2 ~ normal(0,5);
borigin4_matrix_3_4 ~ normal(0,5);
borigin4_matrix_3_8 ~ normal(0,5);
borigin4_matrix_3_18_DE ~ normal(0,5);
borigin4_matrix_3_18_NDE ~ normal(0,5);
borigin4_matrix_3_20_DE ~ normal(0,5);
borigin4_matrix_3_20_NDE ~ normal(0,5);
borigin4_matrix_4_3 ~ normal(0,5);
borigin4_matrix_8_3 ~ normal(0,5);
borigin4_matrix_10_2 ~ normal(0,5);
borigin4_matrix_12_2 ~ normal(0,5);
borigin4_matrix_14_2 ~ normal(0,5);
borigin4_matrix_16_2 ~ normal(0,5);
borigin4_matrix_18_3 ~ normal(0,5);
borigin4_matrix_20_3 ~ normal(0,5);
borigin4_matrix_39_2 ~ normal(0,5);

borigin5_matrix_1_2 ~ normal(0,5);
borigin5_matrix_2_1 ~ normal(0,5);
borigin5_matrix_2_3 ~ normal(0,5);
borigin5_matrix_2_10_DE ~ normal(0,5);
borigin5_matrix_2_10_NDE ~ normal(0,5);
borigin5_matrix_2_12_DE ~ normal(0,5);
borigin5_matrix_2_12_NDE ~ normal(0,5);
borigin5_matrix_2_14_DE ~ normal(0,5);
borigin5_matrix_2_14_NDE ~ normal(0,5);
borigin5_matrix_2_16_DE ~ normal(0,5);
borigin5_matrix_2_16_NDE ~ normal(0,5);
borigin5_matrix_2_39 ~ normal(0,5);
borigin5_matrix_3_2 ~ normal(0,5);
borigin5_matrix_3_4 ~ normal(0,5);
borigin5_matrix_3_8 ~ normal(0,5);
borigin5_matrix_3_18_DE ~ normal(0,5);
borigin5_matrix_3_18_NDE ~ normal(0,5);
borigin5_matrix_3_20_DE ~ normal(0,5);
borigin5_matrix_3_20_NDE ~ normal(0,5);
borigin5_matrix_4_3 ~ normal(0,5);
borigin5_matrix_8_3 ~ normal(0,5);
borigin5_matrix_10_2 ~ normal(0,5);
borigin5_matrix_12_2 ~ normal(0,5);
borigin5_matrix_14_2 ~ normal(0,5);
borigin5_matrix_16_2 ~ normal(0,5);
borigin5_matrix_18_3 ~ normal(0,5);
borigin5_matrix_20_3 ~ normal(0,5);
borigin5_matrix_39_2 ~ normal(0,5);

borigin6_matrix_1_2 ~ normal(0,5);
borigin6_matrix_2_1 ~ normal(0,5);
borigin6_matrix_2_3 ~ normal(0,5);
borigin6_matrix_2_10_DE ~ normal(0,5);
borigin6_matrix_2_10_NDE ~ normal(0,5);
borigin6_matrix_2_12_DE ~ normal(0,5);
borigin6_matrix_2_12_NDE ~ normal(0,5);
borigin6_matrix_2_14_DE ~ normal(0,5);
borigin6_matrix_2_14_NDE ~ normal(0,5);
borigin6_matrix_2_16_DE ~ normal(0,5);
borigin6_matrix_2_16_NDE ~ normal(0,5);
borigin6_matrix_2_39 ~ normal(0,5);
borigin6_matrix_3_2 ~ normal(0,5);
borigin6_matrix_3_4 ~ normal(0,5);
borigin6_matrix_3_8 ~ normal(0,5);
borigin6_matrix_3_18_DE ~ normal(0,5);
borigin6_matrix_3_18_NDE ~ normal(0,5);
borigin6_matrix_3_20_DE ~ normal(0,5);
borigin6_matrix_3_20_NDE ~ normal(0,5);
borigin6_matrix_4_3 ~ normal(0,5);
borigin6_matrix_8_3 ~ normal(0,5);
borigin6_matrix_10_2 ~ normal(0,5);
borigin6_matrix_12_2 ~ normal(0,5);
borigin6_matrix_14_2 ~ normal(0,5);
borigin6_matrix_16_2 ~ normal(0,5);
borigin6_matrix_18_3 ~ normal(0,5);
borigin6_matrix_20_3 ~ normal(0,5);
borigin6_matrix_39_2 ~ normal(0,5);

// priors for year effects

// priors for raw vectors
// prior on the byear_raw parameters arrays
// byear_raw_parameters_array_DE ~ normal(0,1);
// byear_raw_parameters_array_NDE ~ normal(0,1);

// these are now all vectors that get a prior, which is normal(0,1) (and are tranformed using the sigma parameters)
byearxorigin1_raw_vector_1_2 ~ normal(0,1);
byearxorigin1_raw_vector_2_1 ~ normal(0,1);
byearxorigin1_raw_vector_2_3 ~ normal(0,1);
byearxorigin1_raw_vector_2_10_DE ~ normal(0,1);
byearxorigin1_raw_vector_2_10_NDE ~ normal(0,1);
byearxorigin1_raw_vector_2_12_DE ~ normal(0,1);
byearxorigin1_raw_vector_2_12_NDE ~ normal(0,1);
byearxorigin1_raw_vector_2_14_DE ~ normal(0,1);
byearxorigin1_raw_vector_2_14_NDE ~ normal(0,1);
byearxorigin1_raw_vector_2_16_DE ~ normal(0,1);
byearxorigin1_raw_vector_2_16_NDE ~ normal(0,1);
byearxorigin1_raw_vector_2_39 ~ normal(0,1);
byearxorigin1_raw_vector_3_2 ~ normal(0,1);
byearxorigin1_raw_vector_3_4 ~ normal(0,1);
byearxorigin1_raw_vector_3_8 ~ normal(0,1);
byearxorigin1_raw_vector_3_18_DE ~ normal(0,1);
byearxorigin1_raw_vector_3_18_NDE ~ normal(0,1);
byearxorigin1_raw_vector_3_20_DE ~ normal(0,1);
byearxorigin1_raw_vector_3_20_NDE ~ normal(0,1);
byearxorigin1_raw_vector_4_3 ~ normal(0,1);
byearxorigin1_raw_vector_8_3 ~ normal(0,1);

byearxorigin2_raw_vector_1_2 ~ normal(0,1);
byearxorigin2_raw_vector_2_1 ~ normal(0,1);
byearxorigin2_raw_vector_2_3 ~ normal(0,1);
byearxorigin2_raw_vector_2_10_DE ~ normal(0,1);
byearxorigin2_raw_vector_2_10_NDE ~ normal(0,1);
byearxorigin2_raw_vector_2_12_DE ~ normal(0,1);
byearxorigin2_raw_vector_2_12_NDE ~ normal(0,1);
byearxorigin2_raw_vector_2_14_DE ~ normal(0,1);
byearxorigin2_raw_vector_2_14_NDE ~ normal(0,1);
byearxorigin2_raw_vector_2_16_DE ~ normal(0,1);
byearxorigin2_raw_vector_2_16_NDE ~ normal(0,1);
byearxorigin2_raw_vector_2_39 ~ normal(0,1);
byearxorigin2_raw_vector_3_2 ~ normal(0,1);
byearxorigin2_raw_vector_3_4 ~ normal(0,1);
byearxorigin2_raw_vector_3_8 ~ normal(0,1);
byearxorigin2_raw_vector_3_18_DE ~ normal(0,1);
byearxorigin2_raw_vector_3_18_NDE ~ normal(0,1);
byearxorigin2_raw_vector_3_20_DE ~ normal(0,1);
byearxorigin2_raw_vector_3_20_NDE ~ normal(0,1);
byearxorigin2_raw_vector_4_3 ~ normal(0,1);
byearxorigin2_raw_vector_8_3 ~ normal(0,1);

byearxorigin3_raw_vector_1_2 ~ normal(0,1);
byearxorigin3_raw_vector_2_1 ~ normal(0,1);
byearxorigin3_raw_vector_2_3 ~ normal(0,1);
byearxorigin3_raw_vector_2_10_DE ~ normal(0,1);
byearxorigin3_raw_vector_2_10_NDE ~ normal(0,1);
byearxorigin3_raw_vector_2_12_DE ~ normal(0,1);
byearxorigin3_raw_vector_2_12_NDE ~ normal(0,1);
byearxorigin3_raw_vector_2_14_DE ~ normal(0,1);
byearxorigin3_raw_vector_2_14_NDE ~ normal(0,1);
byearxorigin3_raw_vector_2_16_DE ~ normal(0,1);
byearxorigin3_raw_vector_2_16_NDE ~ normal(0,1);
byearxorigin3_raw_vector_2_39 ~ normal(0,1);
byearxorigin3_raw_vector_3_2 ~ normal(0,1);
byearxorigin3_raw_vector_3_4 ~ normal(0,1);
byearxorigin3_raw_vector_3_8 ~ normal(0,1);
byearxorigin3_raw_vector_3_18_DE ~ normal(0,1);
byearxorigin3_raw_vector_3_18_NDE ~ normal(0,1);
byearxorigin3_raw_vector_3_20_DE ~ normal(0,1);
byearxorigin3_raw_vector_3_20_NDE ~ normal(0,1);
byearxorigin3_raw_vector_4_3 ~ normal(0,1);
byearxorigin3_raw_vector_8_3 ~ normal(0,1);

byearxorigin4_raw_vector_1_2 ~ normal(0,1);
byearxorigin4_raw_vector_2_1 ~ normal(0,1);
byearxorigin4_raw_vector_2_3 ~ normal(0,1);
byearxorigin4_raw_vector_2_10_DE ~ normal(0,1);
byearxorigin4_raw_vector_2_10_NDE ~ normal(0,1);
byearxorigin4_raw_vector_2_12_DE ~ normal(0,1);
byearxorigin4_raw_vector_2_12_NDE ~ normal(0,1);
byearxorigin4_raw_vector_2_14_DE ~ normal(0,1);
byearxorigin4_raw_vector_2_14_NDE ~ normal(0,1);
byearxorigin4_raw_vector_2_16_DE ~ normal(0,1);
byearxorigin4_raw_vector_2_16_NDE ~ normal(0,1);
byearxorigin4_raw_vector_2_39 ~ normal(0,1);
byearxorigin4_raw_vector_3_2 ~ normal(0,1);
byearxorigin4_raw_vector_3_4 ~ normal(0,1);
byearxorigin4_raw_vector_3_8 ~ normal(0,1);
byearxorigin4_raw_vector_3_18_DE ~ normal(0,1);
byearxorigin4_raw_vector_3_18_NDE ~ normal(0,1);
byearxorigin4_raw_vector_3_20_DE ~ normal(0,1);
byearxorigin4_raw_vector_3_20_NDE ~ normal(0,1);
byearxorigin4_raw_vector_4_3 ~ normal(0,1);
byearxorigin4_raw_vector_8_3 ~ normal(0,1);

byearxorigin5_raw_vector_1_2 ~ normal(0,1);
byearxorigin5_raw_vector_2_1 ~ normal(0,1);
byearxorigin5_raw_vector_2_3 ~ normal(0,1);
byearxorigin5_raw_vector_2_10_DE ~ normal(0,1);
byearxorigin5_raw_vector_2_10_NDE ~ normal(0,1);
byearxorigin5_raw_vector_2_12_DE ~ normal(0,1);
byearxorigin5_raw_vector_2_12_NDE ~ normal(0,1);
byearxorigin5_raw_vector_2_14_DE ~ normal(0,1);
byearxorigin5_raw_vector_2_14_NDE ~ normal(0,1);
byearxorigin5_raw_vector_2_16_DE ~ normal(0,1);
byearxorigin5_raw_vector_2_16_NDE ~ normal(0,1);
byearxorigin5_raw_vector_2_39 ~ normal(0,1);
byearxorigin5_raw_vector_3_2 ~ normal(0,1);
byearxorigin5_raw_vector_3_4 ~ normal(0,1);
byearxorigin5_raw_vector_3_8 ~ normal(0,1);
byearxorigin5_raw_vector_3_18_DE ~ normal(0,1);
byearxorigin5_raw_vector_3_18_NDE ~ normal(0,1);
byearxorigin5_raw_vector_3_20_DE ~ normal(0,1);
byearxorigin5_raw_vector_3_20_NDE ~ normal(0,1);
byearxorigin5_raw_vector_4_3 ~ normal(0,1);
byearxorigin5_raw_vector_8_3 ~ normal(0,1);

byearxorigin6_raw_vector_1_2 ~ normal(0,1);
byearxorigin6_raw_vector_2_1 ~ normal(0,1);
byearxorigin6_raw_vector_2_3 ~ normal(0,1);
byearxorigin6_raw_vector_2_10_DE ~ normal(0,1);
byearxorigin6_raw_vector_2_10_NDE ~ normal(0,1);
byearxorigin6_raw_vector_2_12_DE ~ normal(0,1);
byearxorigin6_raw_vector_2_12_NDE ~ normal(0,1);
byearxorigin6_raw_vector_2_14_DE ~ normal(0,1);
byearxorigin6_raw_vector_2_14_NDE ~ normal(0,1);
byearxorigin6_raw_vector_2_16_DE ~ normal(0,1);
byearxorigin6_raw_vector_2_16_NDE ~ normal(0,1);
byearxorigin6_raw_vector_2_39 ~ normal(0,1);
byearxorigin6_raw_vector_3_2 ~ normal(0,1);
byearxorigin6_raw_vector_3_4 ~ normal(0,1);
byearxorigin6_raw_vector_3_8 ~ normal(0,1);
byearxorigin6_raw_vector_3_18_DE ~ normal(0,1);
byearxorigin6_raw_vector_3_18_NDE ~ normal(0,1);
byearxorigin6_raw_vector_3_20_DE ~ normal(0,1);
byearxorigin6_raw_vector_3_20_NDE ~ normal(0,1);
byearxorigin6_raw_vector_4_3 ~ normal(0,1);
byearxorigin6_raw_vector_8_3 ~ normal(0,1);

// priors for sigma (scaling factor)
sigma_yearxorigin1_vector_1_2 ~ cauchy(0,1);
sigma_yearxorigin1_vector_2_1 ~ cauchy(0,1);
sigma_yearxorigin1_vector_2_3 ~ cauchy(0,1);
sigma_yearxorigin1_vector_2_10_DE ~ cauchy(0,1);
sigma_yearxorigin1_vector_2_10_NDE ~ cauchy(0,1);
sigma_yearxorigin1_vector_2_12_DE ~ cauchy(0,1);
sigma_yearxorigin1_vector_2_12_NDE ~ cauchy(0,1);
sigma_yearxorigin1_vector_2_14_DE ~ cauchy(0,1);
sigma_yearxorigin1_vector_2_14_NDE ~ cauchy(0,1);
sigma_yearxorigin1_vector_2_16_DE ~ cauchy(0,1);
sigma_yearxorigin1_vector_2_16_NDE ~ cauchy(0,1);
sigma_yearxorigin1_vector_2_39 ~ cauchy(0,1);
// sigma_yearxorigin1_vector_3_2 ~ cauchy(0,1);
sigma_yearxorigin1_vector_3_4 ~ cauchy(0,1);
sigma_yearxorigin1_vector_3_8 ~ cauchy(0,1);
sigma_yearxorigin1_vector_3_18_DE ~ cauchy(0,1);
sigma_yearxorigin1_vector_3_18_NDE ~ cauchy(0,1);
sigma_yearxorigin1_vector_3_20_DE ~ cauchy(0,1);
sigma_yearxorigin1_vector_3_20_NDE ~ cauchy(0,1);
// sigma_yearxorigin1_vector_4_3 ~ cauchy(0,1);
sigma_yearxorigin1_vector_8_3 ~ cauchy(0,1);

sigma_yearxorigin2_vector_1_2 ~ cauchy(0,1);
sigma_yearxorigin2_vector_2_1 ~ cauchy(0,1);
sigma_yearxorigin2_vector_2_3 ~ cauchy(0,1);
sigma_yearxorigin2_vector_2_10_DE ~ cauchy(0,1);
sigma_yearxorigin2_vector_2_10_NDE ~ cauchy(0,1);
sigma_yearxorigin2_vector_2_12_DE ~ cauchy(0,1);
sigma_yearxorigin2_vector_2_12_NDE ~ cauchy(0,1);
sigma_yearxorigin2_vector_2_14_DE ~ cauchy(0,1);
sigma_yearxorigin2_vector_2_14_NDE ~ cauchy(0,1);
sigma_yearxorigin2_vector_2_16_DE ~ cauchy(0,1);
sigma_yearxorigin2_vector_2_16_NDE ~ cauchy(0,1);
sigma_yearxorigin2_vector_2_39 ~ cauchy(0,1);
// sigma_yearxorigin2_vector_3_2 ~ cauchy(0,1);
sigma_yearxorigin2_vector_3_4 ~ cauchy(0,1);
sigma_yearxorigin2_vector_3_8 ~ cauchy(0,1);
sigma_yearxorigin2_vector_3_18_DE ~ cauchy(0,1);
sigma_yearxorigin2_vector_3_18_NDE ~ cauchy(0,1);
sigma_yearxorigin2_vector_3_20_DE ~ cauchy(0,1);
sigma_yearxorigin2_vector_3_20_NDE ~ cauchy(0,1);
// sigma_yearxorigin2_vector_4_3 ~ cauchy(0,1);
sigma_yearxorigin2_vector_8_3 ~ cauchy(0,1);

sigma_yearxorigin3_vector_1_2 ~ cauchy(0,1);
sigma_yearxorigin3_vector_2_1 ~ cauchy(0,1);
sigma_yearxorigin3_vector_2_3 ~ cauchy(0,1);
sigma_yearxorigin3_vector_2_10_DE ~ cauchy(0,1);
sigma_yearxorigin3_vector_2_10_NDE ~ cauchy(0,1);
sigma_yearxorigin3_vector_2_12_DE ~ cauchy(0,1);
sigma_yearxorigin3_vector_2_12_NDE ~ cauchy(0,1);
sigma_yearxorigin3_vector_2_14_DE ~ cauchy(0,1);
sigma_yearxorigin3_vector_2_14_NDE ~ cauchy(0,1);
sigma_yearxorigin3_vector_2_16_DE ~ cauchy(0,1);
sigma_yearxorigin3_vector_2_16_NDE ~ cauchy(0,1);
sigma_yearxorigin3_vector_2_39 ~ cauchy(0,1);
// sigma_yearxorigin3_vector_3_2 ~ cauchy(0,1);
sigma_yearxorigin3_vector_3_4 ~ cauchy(0,1);
sigma_yearxorigin3_vector_3_8 ~ cauchy(0,1);
sigma_yearxorigin3_vector_3_18_DE ~ cauchy(0,1);
sigma_yearxorigin3_vector_3_18_NDE ~ cauchy(0,1);
sigma_yearxorigin3_vector_3_20_DE ~ cauchy(0,1);
sigma_yearxorigin3_vector_3_20_NDE ~ cauchy(0,1);
// sigma_yearxorigin3_vector_4_3 ~ cauchy(0,1);
sigma_yearxorigin3_vector_8_3 ~ cauchy(0,1);

sigma_yearxorigin4_vector_1_2 ~ cauchy(0,1);
sigma_yearxorigin4_vector_2_1 ~ cauchy(0,1);
sigma_yearxorigin4_vector_2_3 ~ cauchy(0,1);
sigma_yearxorigin4_vector_2_10_DE ~ cauchy(0,1);
sigma_yearxorigin4_vector_2_10_NDE ~ cauchy(0,1);
sigma_yearxorigin4_vector_2_12_DE ~ cauchy(0,1);
sigma_yearxorigin4_vector_2_12_NDE ~ cauchy(0,1);
sigma_yearxorigin4_vector_2_14_DE ~ cauchy(0,1);
sigma_yearxorigin4_vector_2_14_NDE ~ cauchy(0,1);
sigma_yearxorigin4_vector_2_16_DE ~ cauchy(0,1);
sigma_yearxorigin4_vector_2_16_NDE ~ cauchy(0,1);
sigma_yearxorigin4_vector_2_39 ~ cauchy(0,1);
// sigma_yearxorigin4_vector_3_2 ~ cauchy(0,1);
sigma_yearxorigin4_vector_3_4 ~ cauchy(0,1);
sigma_yearxorigin4_vector_3_8 ~ cauchy(0,1);
sigma_yearxorigin4_vector_3_18_DE ~ cauchy(0,1);
sigma_yearxorigin4_vector_3_18_NDE ~ cauchy(0,1);
sigma_yearxorigin4_vector_3_20_DE ~ cauchy(0,1);
sigma_yearxorigin4_vector_3_20_NDE ~ cauchy(0,1);
// sigma_yearxorigin4_vector_4_3 ~ cauchy(0,1);
sigma_yearxorigin4_vector_8_3 ~ cauchy(0,1);

sigma_yearxorigin5_vector_1_2 ~ cauchy(0,1);
sigma_yearxorigin5_vector_2_1 ~ cauchy(0,1);
sigma_yearxorigin5_vector_2_3 ~ cauchy(0,1);
sigma_yearxorigin5_vector_2_10_DE ~ cauchy(0,1);
sigma_yearxorigin5_vector_2_10_NDE ~ cauchy(0,1);
sigma_yearxorigin5_vector_2_12_DE ~ cauchy(0,1);
sigma_yearxorigin5_vector_2_12_NDE ~ cauchy(0,1);
sigma_yearxorigin5_vector_2_14_DE ~ cauchy(0,1);
sigma_yearxorigin5_vector_2_14_NDE ~ cauchy(0,1);
sigma_yearxorigin5_vector_2_16_DE ~ cauchy(0,1);
sigma_yearxorigin5_vector_2_16_NDE ~ cauchy(0,1);
sigma_yearxorigin5_vector_2_39 ~ cauchy(0,1);
sigma_yearxorigin5_vector_3_2 ~ cauchy(0,1);
sigma_yearxorigin5_vector_3_4 ~ cauchy(0,1);
sigma_yearxorigin5_vector_3_8 ~ cauchy(0,1);
sigma_yearxorigin5_vector_3_18_DE ~ cauchy(0,1);
sigma_yearxorigin5_vector_3_18_NDE ~ cauchy(0,1);
sigma_yearxorigin5_vector_3_20_DE ~ cauchy(0,1);
sigma_yearxorigin5_vector_3_20_NDE ~ cauchy(0,1);
// sigma_yearxorigin5_vector_4_3 ~ cauchy(0,1);
sigma_yearxorigin5_vector_8_3 ~ cauchy(0,1);

sigma_yearxorigin6_vector_1_2 ~ cauchy(0,1);
sigma_yearxorigin6_vector_2_1 ~ cauchy(0,1);
sigma_yearxorigin6_vector_2_3 ~ cauchy(0,1);
sigma_yearxorigin6_vector_2_10_DE ~ cauchy(0,1);
sigma_yearxorigin6_vector_2_10_NDE ~ cauchy(0,1);
sigma_yearxorigin6_vector_2_12_DE ~ cauchy(0,1);
sigma_yearxorigin6_vector_2_12_NDE ~ cauchy(0,1);
sigma_yearxorigin6_vector_2_14_DE ~ cauchy(0,1);
sigma_yearxorigin6_vector_2_14_NDE ~ cauchy(0,1);
sigma_yearxorigin6_vector_2_16_DE ~ cauchy(0,1);
sigma_yearxorigin6_vector_2_16_NDE ~ cauchy(0,1);
sigma_yearxorigin6_vector_2_39 ~ cauchy(0,1);
sigma_yearxorigin6_vector_3_2 ~ cauchy(0,1);
sigma_yearxorigin6_vector_3_4 ~ cauchy(0,1);
sigma_yearxorigin6_vector_3_8 ~ cauchy(0,1);
sigma_yearxorigin6_vector_3_18_DE ~ cauchy(0,1);
sigma_yearxorigin6_vector_3_18_NDE ~ cauchy(0,1);
sigma_yearxorigin6_vector_3_20_DE ~ cauchy(0,1);
sigma_yearxorigin6_vector_3_20_NDE ~ cauchy(0,1);
// sigma_yearxorigin6_vector_4_3 ~ cauchy(0,1);
sigma_yearxorigin6_vector_8_3 ~ cauchy(0,1);







// Prior on detection efficiency parameters - from the other stan script for detection efficiency
// here, write out all of the parameters for detection efficiency
// twenty terms for intercepts for different eras (configurations of antennas) in the different tributaries
// the outputs from that model will be the first column [,1] containing central tendency (mean)
// and the second column [,2] containing the standard deviation
asotin_alpha1 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
asotin_alpha2 ~ normal(det_eff_param_posteriors[2,1], det_eff_param_posteriors[2,2]);
deschutes_alpha1 ~ normal(det_eff_param_posteriors[3,1], det_eff_param_posteriors[3,2]);
entiat_alpha1 ~ normal(det_eff_param_posteriors[4,1], det_eff_param_posteriors[4,2]);
fifteenmile_alpha1 ~ normal(det_eff_param_posteriors[5,1], det_eff_param_posteriors[5,2]);
imnaha_alpha1 ~ normal(det_eff_param_posteriors[6,1], det_eff_param_posteriors[6,2]);
john_day_alpha1 ~ normal(det_eff_param_posteriors[7,1], det_eff_param_posteriors[7,2]);
methow_alpha1 ~ normal(det_eff_param_posteriors[8,1], det_eff_param_posteriors[8,2]);
methow_alpha2 ~ normal(det_eff_param_posteriors[9,1], det_eff_param_posteriors[9,2]);
okanogan_alpha1 ~ normal(det_eff_param_posteriors[10,1], det_eff_param_posteriors[10,2]);
tucannon_alpha1 ~ normal(det_eff_param_posteriors[11,1], det_eff_param_posteriors[11,2]);
tucannon_alpha2 ~ normal(det_eff_param_posteriors[12,1], det_eff_param_posteriors[12,2]);
umatilla_alpha1 ~ normal(det_eff_param_posteriors[13,1], det_eff_param_posteriors[13,2]);
umatilla_alpha2 ~ normal(det_eff_param_posteriors[14,1], det_eff_param_posteriors[14,2]);
walla_walla_alpha1 ~ normal(det_eff_param_posteriors[15,1], det_eff_param_posteriors[15,2]);
walla_walla_alpha2 ~ normal(det_eff_param_posteriors[16,1], det_eff_param_posteriors[16,2]);
walla_walla_alpha3 ~ normal(det_eff_param_posteriors[17,1], det_eff_param_posteriors[17,2]);
walla_walla_alpha4 ~ normal(det_eff_param_posteriors[18,1], det_eff_param_posteriors[18,2]);
wenatchee_alpha1 ~ normal(det_eff_param_posteriors[19,1], det_eff_param_posteriors[19,2]);
yakima_alpha1 ~ normal(det_eff_param_posteriors[20,1], det_eff_param_posteriors[20,2]);

// 14 terms for discharge relationship, one for each tributary
asotin_beta ~ normal(det_eff_param_posteriors[21,1], det_eff_param_posteriors[21,2]);
deschutes_beta ~ normal(det_eff_param_posteriors[22,1], det_eff_param_posteriors[22,2]);
entiat_beta ~ normal(det_eff_param_posteriors[23,1], det_eff_param_posteriors[23,2]);
// // fifteenmile_beta ~ normal(det_eff_param_posteriors[25,1], det_eff_param_posteriors[25,2]);
// // note that fifteenmile and imnaha don't have any discharge data, so they just get intercepts and the beta term is fixed to zero
// // fifteenmile_beta ~ normal(0,0.000001);
// // imnaha_beta ~ normal(det_eff_param_posteriors[27,1], det_eff_param_posteriors[27,2]);
// // imnaha_beta ~ normal(0,0.000001);
john_day_beta ~ normal(det_eff_param_posteriors[26,1], det_eff_param_posteriors[26,2]);
methow_beta ~ normal(det_eff_param_posteriors[27,1], det_eff_param_posteriors[27,2]);
okanogan_beta ~ normal(det_eff_param_posteriors[28,1], det_eff_param_posteriors[28,2]);
tucannon_beta ~ normal(det_eff_param_posteriors[29,1], det_eff_param_posteriors[29,2]);
umatilla_beta ~ normal(det_eff_param_posteriors[30,1], det_eff_param_posteriors[30,2]);
walla_walla_beta ~ normal(det_eff_param_posteriors[31,1], det_eff_param_posteriors[31,2]);
wenatchee_beta ~ normal(det_eff_param_posteriors[32,1], det_eff_param_posteriors[32,2]);
yakima_beta ~ normal(det_eff_param_posteriors[33,1], det_eff_param_posteriors[33,2]);


// PARALLELIZATION EDITS 
// What we will do is modify the incremental log density statement, to bring target into the first loop by individual.
// This should allow us to increment log density across fish, rather than observations within a fish, allowing us to 
// break up the dataset by fish and therefore run different chunks of the dataset in parallel.

  target += reduce_sum(partial_sum_lupmf, y, grainsize, n_ind, n_states, n_temp_days, max_visits, nmovements, parameter_indices_matrix, transition_dates, // arguments 1-7
  transition_seasons_vector, temperature_data, cat_X_mat, temp_X_mat, year_X_mat, states_mat, n_obs, // arguments 8-13
  spill_window_data, winter_spill_days_data, // arguments 14-18
  winter_post_overshoot_vector, // arguments 19-22
  b0_vector_DE, borigin1_vector_DE, borigin2_vector_DE, borigin3_vector_DE, borigin4_vector_DE, borigin5_vector_DE, borigin6_vector_DE,  // arguments 14-16
  b0_vector_NDE, borigin1_vector_NDE, borigin2_vector_NDE, borigin3_vector_NDE, borigin4_vector_NDE, borigin5_vector_NDE, borigin6_vector_NDE, // arguments 17-19
  btemp0_vector_DE, btemp0xorigin1_vector_DE, btemp0xorigin2_vector_DE, 
  btemp0xorigin3_vector_DE, btemp0xorigin4_vector_DE, btemp0xorigin5_vector_DE, btemp0xorigin6_vector_DE,// arguments 20-23
  btemp0_vector_NDE, btemp0xorigin1_vector_NDE, btemp0xorigin2_vector_NDE, 
  btemp0xorigin3_vector_NDE, btemp0xorigin4_vector_NDE, btemp0xorigin5_vector_NDE, btemp0xorigin6_vector_NDE, // arguments 24-27
  btemp1_vector_DE, btemp1xorigin1_vector_DE, btemp1xorigin2_vector_DE,
  btemp1xorigin3_vector_DE, btemp1xorigin4_vector_DE,
  btemp1xorigin5_vector_DE, btemp1xorigin6_vector_DE,// arguments 28-31
  btemp1_vector_NDE, btemp1xorigin1_vector_NDE, btemp1xorigin2_vector_NDE, 
  btemp1xorigin3_vector_NDE, btemp1xorigin4_vector_NDE,
  btemp1xorigin5_vector_NDE, btemp1xorigin6_vector_NDE,// arguments 32-35
  
  bspillwindow_vector, bwinterspill_vector,
  
  sigma_yearxorigin1_vector_DE, sigma_yearxorigin1_vector_NDE,
  sigma_yearxorigin2_vector_DE, sigma_yearxorigin2_vector_NDE,
  sigma_yearxorigin3_vector_DE, sigma_yearxorigin3_vector_NDE,
  sigma_yearxorigin4_vector_DE, sigma_yearxorigin4_vector_NDE,
  sigma_yearxorigin5_vector_DE, sigma_yearxorigin5_vector_NDE,
  sigma_yearxorigin6_vector_DE, sigma_yearxorigin6_vector_NDE,// arguments 36-38
  byearxorigin1_raw_parameters_matrix_DE, byearxorigin1_raw_parameters_matrix_NDE,
  byearxorigin2_raw_parameters_matrix_DE, byearxorigin2_raw_parameters_matrix_NDE,
  byearxorigin3_raw_parameters_matrix_DE, byearxorigin3_raw_parameters_matrix_NDE,
  byearxorigin4_raw_parameters_matrix_DE, byearxorigin4_raw_parameters_matrix_NDE,
  byearxorigin5_raw_parameters_matrix_DE, byearxorigin5_raw_parameters_matrix_NDE,
  byearxorigin6_raw_parameters_matrix_DE, byearxorigin6_raw_parameters_matrix_NDE,// arguments 39-40
  byearxorigin1_actual_parameters_matrix_DE, byearxorigin1_actual_parameters_matrix_NDE,
  byearxorigin2_actual_parameters_matrix_DE, byearxorigin2_actual_parameters_matrix_NDE,
  byearxorigin3_actual_parameters_matrix_DE, byearxorigin3_actual_parameters_matrix_NDE,
  byearxorigin4_actual_parameters_matrix_DE, byearxorigin4_actual_parameters_matrix_NDE,
  byearxorigin5_actual_parameters_matrix_DE, byearxorigin5_actual_parameters_matrix_NDE,
  byearxorigin6_actual_parameters_matrix_DE, byearxorigin6_actual_parameters_matrix_NDE,// arguments 39-40
  tributary_design_matrices_array, transition_run_years, 
  nyears,
  run_year_DE_array, det_eff_param_vector); // arguments 43-47

}

