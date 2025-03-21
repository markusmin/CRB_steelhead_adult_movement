# 07-07-04_wawa_ff_cov.R

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

### Walla Walla River ###

# Walla Walla River Steelhead: compare how homing changes based on winter spill days at Ice Harbor Dam
# compare homing at 0, 10, 20, 30 March spill days
ICH_winter_spill_values <- c(0, 0.1, 0.2, 0.3)

# compare at the hottest, coldest, and average conditions
fix_run_years <- rep_years$fix_run_year[c(1,3,5)]


WAWA_ICH_winterspill_homing <- data.frame()

for (i in 1:length(ICH_winter_spill_values)){
  for (j in 1:length(fix_run_years)){
    WAWA_FF <- compare_final_fate_fixcov_rear_type_MC(niter = ff_iter, nsim = ff_nsim,
                                                      origin_select = "Walla Walla River",
                                                      fix_state = 8, fix_temp_season0_value = NA,
                                                      fix_temp_season1_value = NA, 
                                                      fix_spillwindow_season0_value = NA,
                                                      fix_spillwindow_season1_value = NA, 
                                                      fix_winterspill_value = ICH_winter_spill_values[i],
                                                      fix_run_year = fix_run_years[j])
    WAWA_FF$ICH_winterspill <-  ICH_winter_spill_values[i]
    WAWA_FF$fix_run_year <-  fix_run_years[j]
    
    WAWA_homing <- subset(WAWA_FF, state == "Walla Walla River")
    
    WAWA_ICH_winterspill_homing %>% 
      bind_rows(., WAWA_homing) -> WAWA_ICH_winterspill_homing
    
  }
}

WAWA_ICH_winterspill_homing %>% 
  mutate(ICH_winterspill_actual = ICH_winterspill*100) -> WAWA_ICH_winterspill_homing

WAWA_ICH_winterspill_homing %>% 
  left_join(rep_years, by = "fix_run_year") -> WAWA_ICH_winterspill_homing

# save data
save(WAWA_ICH_winterspill_homing, file = "figures/alternative_spill_treatment/final_fates_scenarios/simulation_runs/WAWA_ICH_marchspill_homing.rda")