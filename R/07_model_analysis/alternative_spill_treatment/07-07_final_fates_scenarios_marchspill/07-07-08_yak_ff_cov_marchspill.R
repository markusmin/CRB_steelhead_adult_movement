# 07-07-08_yak_ff_cov.R

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



### Yakima River ###

# Yakima River Steelhead: compare how homing changes based on winter spill days at Priest Rapids Dam
# compare homing at 0, 10, 20, 30 March spill days
PRA_winter_spill_values <- c(0, 0.1, 0.2, 0.3)

# cool year: 08/09
# average year: 19/20
# warm year: 17/18

# coldest year on record: 07/08
# hottest year on record: 15/16
# compare at the hottest, coldest, and average conditions
fix_run_years <- rep_years$fix_run_year[c(1,3,5)]

YAK_PRA_winterspill_homing <- data.frame()

for (i in 1:length(PRA_winter_spill_values)){
  for (j in 1:length(fix_run_years)){
    YAK_FF <- compare_final_fate_fixcov_rear_type_MC(niter = ff_iter, nsim = ff_nsim,
                                                     origin_select = "Yakima River",
                                                     fix_state = 4, fix_temp_season0_value = NA,
                                                     fix_temp_season1_value = NA, 
                                                     fix_spillwindow_season0_value = NA,
                                                     fix_spillwindow_season1_value = NA, 
                                                     fix_winterspill_value = PRA_winter_spill_values[i],
                                                     fix_run_year = fix_run_years[j])
    YAK_FF$PRA_winterspill <-  PRA_winter_spill_values[i]
    YAK_FF$fix_run_year <-  fix_run_years[j]
    
    YAK_homing <- subset(YAK_FF, state == "Yakima River")
    
    YAK_PRA_winterspill_homing %>% 
      bind_rows(., YAK_homing) -> YAK_PRA_winterspill_homing
    
  }
}

YAK_PRA_winterspill_homing %>% 
  mutate(PRA_winterspill_actual = PRA_winterspill*100) -> YAK_PRA_winterspill_homing

YAK_PRA_winterspill_homing %>% 
  left_join(rep_years, by = "fix_run_year") -> YAK_PRA_winterspill_homing


# save output
save(YAK_PRA_winterspill_homing, file = "figures/alternative_spill_treatment/final_fates_scenarios/simulation_runs/YAK_PRA_marchspill_homing.rda")