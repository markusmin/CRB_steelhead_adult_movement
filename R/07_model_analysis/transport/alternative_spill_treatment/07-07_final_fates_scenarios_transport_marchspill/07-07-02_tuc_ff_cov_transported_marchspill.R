# 07-07-06_tuc_ff_cov.R

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

source("R/07_model_analysis/transport/alternative_spill_treatment/07-07_final_fates_scenarios_transport_marchspill/07-07-01-ff_cov_functions_transport_marchspill.R")


#### March spill days comparison ####

### Tucannon River ###

# Tucannon River Steelhead: compare how homing changes based on winter spill days at Lower Granite Dam
# compare homing at 0, 10, 20, 30 March spill days
LGR_winter_spill_values <- c(0, 0.1, 0.2, 0.3)

# cool year: 08/09
# average year: 19/20
# warm year: 17/18

# coldest year on record: 07/08
# hottest year on record: 15/16
# compare at the hottest, coldest, and average conditions
fix_run_years <- rep_years$fix_run_year[c(1,3,5)]

TUC_LGR_winterspill_homing <- data.frame()

for (i in 1:length(LGR_winter_spill_values)){
  for (j in 1:length(fix_run_years)){
    TUC_FF <- compare_final_fate_fixcov_rear_type_SR_T(niter = ff_iter, nsim = ff_nsim,
                                                     origin_select = "Tucannon River",
                                                     fix_state = 9, fix_temp_season0_value = NA,
                                                     fix_temp_season1_value = NA, 
                                                     fix_spillwindow_season0_value = NA,
                                                     fix_spillwindow_season1_value = NA, 
                                                     fix_winterspill_value = LGR_winter_spill_values[i],
                                                     fix_run_year = fix_run_years[j])
    TUC_FF$LGR_winterspill <-  LGR_winter_spill_values[i]
    TUC_FF$fix_run_year <-  fix_run_years[j]
    
    TUC_homing <- subset(TUC_FF, state == "Tucannon River")
    
    TUC_LGR_winterspill_homing %>% 
      bind_rows(., TUC_homing) -> TUC_LGR_winterspill_homing
    
  }
}

TUC_LGR_winterspill_homing %>% 
  mutate(LGR_winterspill_actual = LGR_winterspill*100) -> TUC_LGR_winterspill_homing

TUC_LGR_winterspill_homing %>% 
  left_join(rep_years, by = "fix_run_year") -> TUC_LGR_winterspill_homing


# save output
save(TUC_LGR_winterspill_homing, file = "figures/transport/alternative_spill_treatment/final_fates_scenarios/simulation_runs/TUC_LGR_transported_homing.rda")