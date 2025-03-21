# 07-07-02_jdr_ff_cov.R

# This script runs final fates + covariates simulations.

#### Load libraries, state information ####
library(cmdstanr)
library(rstan)
library(bayesplot)
library(posterior)
library(dplyr)
library(tidyr)
library(ggpubr)
library(stringr)
library(rlang)
library(tibble)
library(forcats)
library(lubridate)

# set the working directory so that here() will cooperate
setwd("/gscratch/scrubbed/mmin/")

source("R/07_model_analysis/alternative_spill_treatment/07-07_final_fates_scenarios_marchspill/07-07-01-ff_cov_functions_marchspill.R")


#### Winter spill days comparison ####


# compare at the hottest, coldest, and average conditions
fix_run_years <- rep_years$fix_run_year[c(1,3,5)]

### John Day River ###

# John Day River Steelhead: compare how homing changes based on winter spill days at McNary Dam
# compare homing at 0, 10, 20, 30 March spill days
MCN_winter_spill_values <- c(0, 0.1, 0.2, 0.3)

JDR_MCN_winterspill_homing <- data.frame()

for (i in 1:length(MCN_winter_spill_values)){
  for (j in 1:length(fix_run_years)){
    JDR_FF <- compare_final_fate_fixcov_rear_type_MC(niter = ff_iter, nsim = ff_nsim,
                                                      origin_select = "John Day River",
                                                      fix_state = 8, fix_temp_season0_value = NA,
                                                      fix_temp_season1_value = NA, 
                                                      fix_spillwindow_season0_value = NA,
                                                      fix_spillwindow_season1_value = NA, 
                                                      fix_winterspill_value = MCN_winter_spill_values[i],
                                                      fix_run_year = fix_run_years[j])
    JDR_FF$MCN_winterspill <-  MCN_winter_spill_values[i]
    JDR_FF$fix_run_year <-  fix_run_years[j]
    
    JDR_homing <- subset(JDR_FF, state == "John Day River")
    
    JDR_MCN_winterspill_homing %>% 
      bind_rows(., JDR_homing) -> JDR_MCN_winterspill_homing
    
  }
}

JDR_MCN_winterspill_homing %>% 
  mutate(MCN_winterspill_actual = MCN_winterspill*100) -> JDR_MCN_winterspill_homing

JDR_MCN_winterspill_homing %>% 
  left_join(rep_years, by = "fix_run_year") -> JDR_MCN_winterspill_homing

# save data
save(JDR_MCN_winterspill_homing, file = "figures/alternative_spill_treatment/final_fates_scenarios/simulation_runs/JDR_MCN_marchspill_homing.rda")
