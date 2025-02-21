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
library(here)

source(here::here("analysis", "analysis", "14.0-ff-cov-functions.R"))


#### Winter spill days comparison ####

### Tucannon River ###

# Tucannon River Steelhead: compare how homing changes based on winter spill days at Lower Granite Dam
# compare homing at 0, 30, 60, 90 winter spill days
LGR_winter_spill_values <- c(0, 0.30, 0.60, 0.90)

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
    TUC_FF <- compare_final_fate_fixcov_rear_type_SR(niter = ff_iter, nsim = ff_nsim,
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
save(TUC_LGR_winterspill_homing, file = here::here("stan_actual", "output", "final_fates_covariates", "TUC_LGR_winterspill_homing.rda"))

rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")
rear_shapes <- c(17, 19)
condition_colors <- c("coldest" = "#2c7bb6", "average" = "#ffffbf", "warmest" = "#d7191c")

TUC_LGR_winterspill_homing_plot <- ggplot(TUC_LGR_winterspill_homing, aes(x = LGR_winterspill_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, color = conditions,
                                                                          shape = rear_type)) +
  geom_point(size = 3.5, position=position_dodge(width=3)) +
  geom_linerange(position=position_dodge(width=3)) +
  ylab("Homing Probability") +
  xlab("Days of Winter Spill at Lower Granite Dam") +
  # ggtitle(" ") +
  # Create a scale common to all
  scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_shape_manual(values = rear_shapes) +
  scale_color_manual(values = condition_colors) +
  theme(plot.title = element_text(size = 12),
        # axis.text.y = element_text(color = rev(state_significance_colors)),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12)) +
  ggtitle("Homing by Tucannon River Steelhead under different winter spill conditions at Lower Granite Dam")

ggsave(here::here("stan_actual", "output", "final_fates_covariates", "TUC_LGR_winterspill_homing_plot_temps.png"), TUC_LGR_winterspill_homing_plot, height = 8, width = 8)