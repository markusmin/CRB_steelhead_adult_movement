# 07-07-07_ent_ff_cov.R

# This script takesruns final fates + covariates simulations.

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

### Entiat River ###

# Entiat River Steelhead: compare how homing changes based on winter spill days at Wells Dam
# compare homing at 0, 10, 20, 30 March spill days
WEL_winter_spill_values <- c(0, 0.1, 0.2, 0.3)

# cool year: 08/09
# average year: 19/20
# warm year: 17/18

# coldest year on record: 07/08
# hottest year on record: 15/16
# compare at the hottest, coldest, and average conditions
fix_run_years <- rep_years$fix_run_year[c(1,3,5)]

ENT_WEL_winterspill_homing <- data.frame()

for (i in 1:length(WEL_winter_spill_values)){
  for (j in 1:length(fix_run_years)){
    ENT_FF <- compare_final_fate_fixcov_rear_type_UC(niter = ff_iter, nsim = ff_nsim,
                                                     origin_select = "Entiat River",
                                                     fix_state = 7, fix_temp_season0_value = NA,
                                                     fix_temp_season1_value = NA, 
                                                     fix_spillwindow_season0_value = NA,
                                                     fix_spillwindow_season1_value = NA, 
                                                     fix_winterspill_value = WEL_winter_spill_values[i],
                                                     fix_run_year = fix_run_years[j])
    ENT_FF$WEL_winterspill <-  WEL_winter_spill_values[i]
    ENT_FF$fix_run_year <-  fix_run_years[j]
    
    ENT_homing <- subset(ENT_FF, state == "Entiat River")
    
    ENT_WEL_winterspill_homing %>% 
      bind_rows(., ENT_homing) -> ENT_WEL_winterspill_homing
    
  }
}

ENT_WEL_winterspill_homing %>% 
  mutate(WEL_winterspill_actual = WEL_winterspill*100) -> ENT_WEL_winterspill_homing

ENT_WEL_winterspill_homing %>% 
  left_join(rep_years, by = "fix_run_year") -> ENT_WEL_winterspill_homing


# save output
save(ENT_WEL_winterspill_homing, file = "figures/alternative_spill_treatment/final_fates_scenarios/simulation_runs/ENT_WEL_marchspill_homing.rda")