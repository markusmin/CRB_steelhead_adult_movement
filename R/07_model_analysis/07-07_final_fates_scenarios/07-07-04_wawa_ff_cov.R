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

source("R/07_model_analysis/07-07_final_fates_scenarios/07-07-01-ff_cov_functions.R")


#### Winter spill days comparison ####

### Walla Walla River ###

# Walla Walla River Steelhead: compare how homing changes based on winter spill days at Ice Harbor Dam
# compare homing at 0, 30, 60, 90 spill days
ICH_winter_spill_values <- c(0, 0.30, 0.60, 0.90)

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

save(WAWA_ICH_winterspill_homing, file = "figures/final_fates_scenarios/simulation_runs/WAWA_ICH_winterspill_homing.rda")

rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")
rear_shapes <- c(17, 19)
WAWA_ICH_winterspill_homing$conditions <- factor(WAWA_ICH_winterspill_homing$conditions, levels = c("coldest", "average", "warmest"))
condition_colors <- c("coldest" = "#2c7bb6", "average" = "#ffffbf", "warmest" = "#d7191c")

WAWA_ICH_winterspill_homing_plot <- ggplot(WAWA_ICH_winterspill_homing, aes(x = ICH_winterspill_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, color = conditions,
                                                                            shape = rear_type)) +
  geom_point(size = 3.5, position=position_dodge(width=3)) +
  geom_linerange(position=position_dodge(width=3)) +
  ylab("Homing Probability") +
  xlab("Days of Winter Spill at Ice Harbor Dam") +
  # ggtitle(" ") +
  # Create a scale common to all
  scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_shape_manual(values = rear_shapes) +
  scale_color_manual(values = condition_colors) +
  theme(plot.title = element_text(size = 12),
        # axis.text.y = element_text(color = rev(state_significance_colors)),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12)) +
  ggtitle("Homing by Walla Walla River Steelhead under different winter spill conditions at Ice Harbor Dam")

ggsave("figures/final_fates_scenarios/WAWA_ICH_winterspill_homing_plot_temps.png", WAWA_ICH_winterspill_homing_plot, height = 8, width = 8)