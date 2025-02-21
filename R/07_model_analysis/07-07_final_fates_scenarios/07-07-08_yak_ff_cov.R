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
library(here)

source(here::here("analysis", "analysis", "14.0-ff-cov-functions.R"))


#### Winter spill days comparison ####



### Yakima River ###

# Yakima River Steelhead: compare how homing changes based on winter spill days at Priest Rapids Dam
# compare homing at 0, 30, 60, 90 spill days
PRA_winter_spill_values <- c(0, 0.30, 0.60, 0.90)

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
save(YAK_PRA_winterspill_homing, file = here::here("stan_actual", "output", "final_fates_covariates", "YAK_PRA_winterspill_homing.rda"))

rear_shapes <- c(17, 19)
rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")
condition_colors <- c("coldest" = "#2c7bb6", "average" = "#ffffbf", "warmest" = "#d7191c")

YAK_PRA_winterspill_homing_plot <- ggplot(YAK_PRA_winterspill_homing, aes(x = PRA_winterspill_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, color = conditions,
                                                                          shape = rear_type)) +
  geom_point(size = 3.5, position=position_dodge(width=3)) +
  geom_linerange(position=position_dodge(width=3)) +
  ylab("Homing Probability") +
  xlab("Days of Winter Spill at Priest Rapids Dam") +
  # ggtitle(" ") +
  # Create a scale common to all
  scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_shape_manual(values = rear_shapes) +
  scale_color_manual(values = condition_colors) +
  theme(plot.title = element_text(size = 12),
        # axis.text.y = element_text(color = rev(state_significance_colors)),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12)) +
  ggtitle("Homing by Yakima River Steelhead under different winter spill conditions at Priest Rapids Dam")

ggsave(here::here("stan_actual", "output", "final_fates_covariates", "YAK_PRA_winterspill_homing_plot_temps.png"), YAK_PRA_winterspill_homing_plot, height = 8, width = 8)