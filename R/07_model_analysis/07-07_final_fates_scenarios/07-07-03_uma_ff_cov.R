# 07-07-03_uma_ff_cov.R

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
library(here)

source(here::here("analysis", "analysis", "14.0-ff-cov-functions.R"))


#### Winter spill days comparison ####


# compare at the hottest, coldest, and average conditions
fix_run_years <- rep_years$fix_run_year[c(1,3,5)]

### Umatilla River ###

# Umatilla River Steelhead: compare how homing changes based on winter spill days at McNary Dam
# compare homing at 0, 30, 60, 90 winter spill days
MCN_winter_spill_values <- c(0, 0.30, 0.60, 0.90)

UMA_MCN_winterspill_homing <- data.frame()

for (i in 1:length(MCN_winter_spill_values)){
  for (j in 1:length(fix_run_years)){
    UMA_FF <- compare_final_fate_fixcov_rear_type_MC(niter = ff_iter, nsim = ff_nsim,
                                                      origin_select = "Umatilla River",
                                                      fix_state = 3, fix_temp_season0_value = NA,
                                                      fix_temp_season1_value = NA, 
                                                      fix_spillwindow_season0_value = NA,
                                                      fix_spillwindow_season1_value = NA, 
                                                      fix_winterspill_value = MCN_winter_spill_values[i],
                                                      fix_run_year = fix_run_years[j])
    UMA_FF$MCN_winterspill <-  MCN_winter_spill_values[i]
    UMA_FF$fix_run_year <-  fix_run_years[j]
    
    UMA_homing <- subset(UMA_FF, state == "Umatilla River")
    
    UMA_MCN_winterspill_homing %>% 
      bind_rows(., UMA_homing) -> UMA_MCN_winterspill_homing
    
  }
}

UMA_MCN_winterspill_homing %>% 
  mutate(MCN_winterspill_actual = MCN_winterspill*100) -> UMA_MCN_winterspill_homing

UMA_MCN_winterspill_homing %>% 
  left_join(rep_years, by = "fix_run_year") -> UMA_MCN_winterspill_homing

# save data
save(UMA_MCN_winterspill_homing, file = here::here("stan_actual", "output", "final_fates_covariates", "UMA_MCN_winterspill_homing.rda"))

rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")
rear_shapes <- c(17, 19)
condition_colors <- c("coldest" = "#2c7bb6", "average" = "#ffffbf", "warmest" = "#d7191c")
UMA_MCN_winterspill_homing_plot <- ggplot(UMA_MCN_winterspill_homing, aes(x = MCN_winterspill_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, color = conditions,
                                                                          shape = rear_type)) +
  geom_point(size = 3.5, position=position_dodge(width=3)) +
  geom_linerange(position=position_dodge(width=3)) +
  ylab("Homing Probability") +
  xlab("Days of Winter Spill at McNary Dam") +
  # ggtitle(" ") +
  # Create a scale common to all
  scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_shape_manual(values = rear_shapes) +
  scale_color_manual(values = condition_colors) +
  theme(plot.title = element_text(size = 12),
        # axis.text.y = element_text(color = rev(state_significance_colors)),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12)) +
  ggtitle("Homing by Umatilla River Steelhead under different winter spill conditions at McNary Dam")

ggsave(here::here("stan_actual", "output", "final_fates_covariates", "UMA_MCN_winterspill_homing_plot_temps.png"), UMA_MCN_winterspill_homing_plot, height = 8, width = 8)