# 07-01_load_stan_models - for transport models, and alternative (March) spill runs

# This script loads and reformats all of our model runs, so that we can run subsequent analyses on them.

# if this is being run on hyak, need to change wd:
setwd("/gscratch/scrubbed/mmin/")

#### Load libraries, state information ####
library(cmdstanr)
library(rstan)
library(bayesplot)
library(posterior)
library(dplyr)
library(tidyr)
library(here)
library(ggpubr)
library(stringr)
library(rlang)
library(tibble)
library(forcats)
library(lubridate)
library(ggthemes)
library(patchwork)
library(cowplot)
library(grid)
library(ggExtra)

# get the model states into a df, to help with interpretation
model_states = c(
  # Mainstem states (9)
  "mainstem, mouth to BON",
  "mainstem, BON to MCN",
  "mainstem, MCN to ICH or PRA",
  "mainstem, PRA to RIS",
  "mainstem, RIS to RRE",
  "mainstem, RRE to WEL",
  "mainstem, upstream of WEL",
  "mainstem, ICH to LGR",
  "mainstem, upstream of LGR",
  
  # Tributary states ()
  # With detection efficiencies in the model, we now have more tributary states,
  # since we have an upstream and a river mouth state
  
  # "Deschutes River", 
  "Deschutes River Mouth", "Deschutes River Upstream",
  # "John Day River", 
  "John Day River Mouth", "John Day River Upstream",
  # "Hood River",
  # "Hood River Mouth", "Hood River Upstream",
  # "Fifteenmile Creek", 
  "Fifteenmile Creek Mouth", "Fifteenmile Creek Upstream",
  # "Umatilla River",
  "Umatilla River Mouth", "Umatilla River Upstream",
  # "Yakima River",
  "Yakima River Mouth", "Yakima River Upstream",
  # "Walla Walla River",
  "Walla Walla River Mouth", "Walla Walla River Upstream",
  # "Wenatchee River", 
  "Wenatchee River Mouth", "Wenatchee River Upstream",
  # "Entiat River", 
  "Entiat River Mouth", "Entiat River Upstream",
  # "Okanogan River", 
  "Okanogan River Mouth", "Okanogan River Upstream",
  # "Methow River", 
  "Methow River Mouth", "Methow River Upstream",
  # "Tucannon River",
  "Tucannon River Mouth", "Tucannon River Upstream",
  # "Asotin Creek", 
  "Asotin Creek Mouth", "Asotin Creek Upstream",
  "Clearwater River",
  "Salmon River",
  "Grande Ronde River",
  # "Imnaha River",
  "Imnaha River Mouth", "Imnaha River Upstream",
  "BON to MCN other tributaries",
  "Upstream WEL other tributaries",
  
  # Loss
  "loss"
)


##### Load the model runs #####

# Load the model data associated with each run (necessary to load covariates)
# Store these each in an environment, because most things share names
SRW_T_marchspill_envir <- new.env()
SRW_NT_marchspill_envir <- new.env()
SRH_T_marchspill_envir <- new.env()
SRH_NT_marchspill_envir <- new.env()
load(here::here("Stan", "alternative_spill_treatment", "snake_river_wild_transported", "SRW_T_marchspill_model_data.rda"),
     envir = SRW_T_marchspill_envir)
load(here::here("Stan", "alternative_spill_treatment", "snake_river_wild_not_transported", "SRW_NT_marchspill_model_data.rda"),
     envir = SRW_NT_marchspill_envir)
load(here::here("Stan", "alternative_spill_treatment", "snake_river_hatchery_transported", "SRH_T_marchspill_model_data.rda"),
     envir = SRH_T_marchspill_envir)
load(here::here("Stan", "alternative_spill_treatment", "snake_river_hatchery_not_transported", "SRH_NT_marchspill_model_data.rda"),
     envir = SRH_NT_marchspill_envir)


# Function to bind four chains together
bind4chains <- function(chain1, chain2, chain3, chain4){
  bound_draws <- bind_draws(chain1$draws(),
                            chain2$draws(),
                            chain3$draws(),
                            chain4$draws(), along = "chain")
  
  return(bound_draws)
}


## Snake River, Wild transported
SRW_T_marchspill_chain1 <- readRDS(here::here("Stan",  "alternative_spill_treatment", "snake_river_wild_transported", "SRW_T_marchspill_chain1_fit.rds"))
SRW_T_marchspill_chain2 <- readRDS(here::here("Stan",  "alternative_spill_treatment", "snake_river_wild_transported", "SRW_T_marchspill_chain2_fit.rds"))
SRW_T_marchspill_chain3 <- readRDS(here::here("Stan",  "alternative_spill_treatment", "snake_river_wild_transported", "SRW_T_marchspill_chain3_fit.rds"))
SRW_T_marchspill_chain4 <- readRDS(here::here("Stan",  "alternative_spill_treatment", "snake_river_wild_transported", "SRW_T_marchspill_chain4_fit.rds"))

# bind chains together
SRW_T_marchspill_fit_raw <- bind4chains(SRW_T_marchspill_chain1, SRW_T_marchspill_chain2, SRW_T_marchspill_chain3, SRW_T_marchspill_chain4)
# thin2
thin_draws(SRW_T_marchspill_fit_raw, thin = 2) -> SRW_T_marchspill_fit
# summarise
SRW_T_marchspill_fit_summary <- summarise_draws(SRW_T_marchspill_fit)

## Snake River, Wild not transported
SRW_NT_marchspill_chain1 <- readRDS(here::here("Stan",  "alternative_spill_treatment", "snake_river_wild_not_transported", "SRW_NT_marchspill_chain1_fit.rds"))
SRW_NT_marchspill_chain2 <- readRDS(here::here("Stan",  "alternative_spill_treatment", "snake_river_wild_not_transported", "SRW_NT_marchspill_chain2_fit.rds"))
SRW_NT_marchspill_chain3 <- readRDS(here::here("Stan",  "alternative_spill_treatment", "snake_river_wild_not_transported", "SRW_NT_marchspill_chain3_fit.rds"))
SRW_NT_marchspill_chain4 <- readRDS(here::here("Stan",  "alternative_spill_treatment", "snake_river_wild_not_transported", "SRW_NT_marchspill_chain4_fit.rds"))

# bind chains together
SRW_NT_marchspill_fit_raw <- bind4chains(SRW_NT_marchspill_chain1, SRW_NT_marchspill_chain2, SRW_NT_marchspill_chain3, SRW_NT_marchspill_chain4)
# thin2
thin_draws(SRW_NT_marchspill_fit_raw, thin = 2) -> SRW_NT_marchspill_fit
# summarise
SRW_NT_marchspill_fit_summary <- summarise_draws(SRW_NT_marchspill_fit)

## Snake River, Hatchery transported
SRH_T_marchspill_chain1 <- readRDS(here::here("Stan",  "alternative_spill_treatment", "snake_river_hatchery_transported", "SRH_T_marchspill_chain1_fit.rds"))
SRH_T_marchspill_chain2 <- readRDS(here::here("Stan",  "alternative_spill_treatment", "snake_river_hatchery_transported", "SRH_T_marchspill_chain2_fit.rds"))
SRH_T_marchspill_chain3 <- readRDS(here::here("Stan",  "alternative_spill_treatment", "snake_river_hatchery_transported", "SRH_T_marchspill_chain3_fit.rds"))
SRH_T_marchspill_chain4 <- readRDS(here::here("Stan",  "alternative_spill_treatment", "snake_river_hatchery_transported", "SRH_T_marchspill_chain4_fit.rds"))

# bind chains together
SRH_T_marchspill_fit_raw <- bind4chains(SRH_T_marchspill_chain1, SRH_T_marchspill_chain2, SRH_T_marchspill_chain3, SRH_T_marchspill_chain4)
# thin2
thin_draws(SRH_T_marchspill_fit_raw, thin = 2) -> SRH_T_marchspill_fit
# summarise
SRH_T_marchspill_fit_summary <- summarise_draws(SRH_T_marchspill_fit)

## Snake River, Hatchery not transported
SRH_NT_marchspill_chain1 <- readRDS(here::here("Stan",  "alternative_spill_treatment", "snake_river_hatchery_not_transported", "SRH_NT_marchspill_chain1_fit.rds"))
SRH_NT_marchspill_chain2 <- readRDS(here::here("Stan",  "alternative_spill_treatment", "snake_river_hatchery_not_transported", "SRH_NT_marchspill_chain2_fit.rds"))
SRH_NT_marchspill_chain3 <- readRDS(here::here("Stan",  "alternative_spill_treatment", "snake_river_hatchery_not_transported", "SRH_NT_marchspill_chain3_fit.rds"))
SRH_NT_marchspill_chain4 <- readRDS(here::here("Stan",  "alternative_spill_treatment", "snake_river_hatchery_not_transported", "SRH_NT_marchspill_chain4_fit.rds"))

# bind chains together
SRH_NT_marchspill_fit_raw <- bind4chains(SRH_NT_marchspill_chain1, SRH_NT_marchspill_chain2, SRH_NT_marchspill_chain3, SRH_NT_marchspill_chain4)
# thin2
thin_draws(SRH_NT_marchspill_fit_raw, thin = 2) -> SRH_NT_marchspill_fit
# summarise
SRH_NT_marchspill_fit_summary <- summarise_draws(SRH_NT_marchspill_fit)

#### Extract all parameter values from the model fit objects ####

# Here, you need to make sure to check the origin params, because those change by DPS

# function to take a parameter type (fixed effect) and store all of them in an array
make_parameter_draws_array <- function(parameter_prefix, fit, fit_summary){
  # extract b0 as an array
  parameters <- fit_summary$variable
  parameters[grepl(paste0(parameter_prefix, "_"), parameters)] -> param_subset
  param_subset[!(grepl("vector", param_subset))] -> param_subset
  # drop the NDE parameters
  param_subset <- param_subset[!(grepl("_NDE", param_subset))]
  
  
  param_subset_from = as.numeric(sub("[^_]*_[^_]*_", "", str_extract(param_subset, "[^_]*_[^_]*_[^_]*")))
  param_subset_to = as.numeric(str_extract(sub("[^_]*_[^_]*_[^_]*_", "", param_subset), "\\d+"))
  
  param_subset_indices <- data.frame(parameter = param_subset, from = param_subset_from, to = param_subset_to)
  
  
  # arrange all parameter values into an array
  param_array <- array(data = 0, dim = c(length(model_states), length(model_states),
                                         length(as.matrix(fit[,,1]))))
  
  # 0 is meaningful for loss (this is what is used in the stan code)
  # 0s for all other movements that are not overwritten will not be used
  
  
  for(i in 1:nrow(param_subset_indices)){
    param_array[param_subset_indices[i, "from"], param_subset_indices[i, "to"], ] <- as.matrix(fit[,,param_subset_indices[i, "parameter"]])
  }
  
  return(param_array)
}

### SRW_T ###
SRW_T_marchspill_parameters <- SRW_T_marchspill_fit_summary$variable
SRW_T_marchspill_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", SRW_T_marchspill_parameters)] -> SRW_T_marchspill_fixed_effects
SRW_T_marchspill_fixed_effects[!(grepl("vector", SRW_T_marchspill_fixed_effects))] -> SRW_T_marchspill_fixed_effects

b0_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "b0", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
btemp0_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
btemp1_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
bspillwindow_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
bwinterspill_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
btemp0xorigin1_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
btemp1xorigin1_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
btemp0xorigin2_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
btemp1xorigin2_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
btemp0xorigin3_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
btemp1xorigin3_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
btemp0xorigin4_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin4", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
btemp1xorigin4_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin4", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
btemp0xorigin5_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin5", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
btemp1xorigin5_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin5", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
btemp0xorigin6_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin6", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
btemp1xorigin6_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin6", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
borigin1_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
borigin2_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
borigin3_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
borigin4_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "borigin4", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
borigin5_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "borigin5", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)
borigin6_array_SRW_T <- make_parameter_draws_array(parameter_prefix = "borigin6", fit = SRW_T_marchspill_fit, fit_summary = SRW_T_marchspill_fit_summary)

### SRW_NT ###
SRW_NT_marchspill_parameters <- SRW_NT_marchspill_fit_summary$variable
SRW_NT_marchspill_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", SRW_NT_marchspill_parameters)] -> SRW_NT_marchspill_fixed_effects
SRW_NT_marchspill_fixed_effects[!(grepl("vector", SRW_NT_marchspill_fixed_effects))] -> SRW_NT_marchspill_fixed_effects

b0_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "b0", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
btemp0_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
btemp1_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
bspillwindow_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
bwinterspill_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
btemp0xorigin1_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
btemp1xorigin1_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
btemp0xorigin2_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
btemp1xorigin2_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
btemp0xorigin3_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
btemp1xorigin3_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
btemp0xorigin4_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin4", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
btemp1xorigin4_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin4", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
btemp0xorigin5_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin5", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
btemp1xorigin5_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin5", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
btemp0xorigin6_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin6", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
btemp1xorigin6_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin6", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
borigin1_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
borigin2_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
borigin3_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
borigin4_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "borigin4", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
borigin5_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "borigin5", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)
borigin6_array_SRW_NT <- make_parameter_draws_array(parameter_prefix = "borigin6", fit = SRW_NT_marchspill_fit, fit_summary = SRW_NT_marchspill_fit_summary)

### SRH_T ###
SRH_T_marchspill_parameters <- SRH_T_marchspill_fit_summary$variable
SRH_T_marchspill_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", SRH_T_marchspill_parameters)] -> SRH_T_marchspill_fixed_effects
SRH_T_marchspill_fixed_effects[!(grepl("vector", SRH_T_marchspill_fixed_effects))] -> SRH_T_marchspill_fixed_effects

b0_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "b0", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
btemp0_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
btemp1_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
bspillwindow_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
bwinterspill_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
btemp0xorigin1_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
btemp1xorigin1_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
btemp0xorigin2_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
btemp1xorigin2_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
btemp0xorigin3_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
btemp1xorigin3_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
btemp0xorigin4_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin4", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
btemp1xorigin4_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin4", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
btemp0xorigin5_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin5", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
btemp1xorigin5_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin5", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
borigin1_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
borigin2_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
borigin3_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
borigin4_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "borigin4", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)
borigin5_array_SRH_T <- make_parameter_draws_array(parameter_prefix = "borigin5", fit = SRH_T_marchspill_fit, fit_summary = SRH_T_marchspill_fit_summary)

### SRH_NT ###
SRH_NT_marchspill_parameters <- SRH_NT_marchspill_fit_summary$variable
SRH_NT_marchspill_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", SRH_NT_marchspill_parameters)] -> SRH_NT_marchspill_fixed_effects
SRH_NT_marchspill_fixed_effects[!(grepl("vector", SRH_NT_marchspill_fixed_effects))] -> SRH_NT_marchspill_fixed_effects

b0_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "b0", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
btemp0_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
btemp1_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
bspillwindow_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
bwinterspill_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
btemp0xorigin1_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
btemp1xorigin1_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
btemp0xorigin2_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
btemp1xorigin2_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
btemp0xorigin3_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
btemp1xorigin3_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
btemp0xorigin4_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin4", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
btemp1xorigin4_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin4", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
btemp0xorigin5_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin5", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
btemp1xorigin5_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin5", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
borigin1_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
borigin2_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
borigin3_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
borigin4_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "borigin4", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)
borigin5_array_SRH_NT <- make_parameter_draws_array(parameter_prefix = "borigin5", fit = SRH_NT_marchspill_fit, fit_summary = SRH_NT_marchspill_fit_summary)

# Make the origin_param_map, since a lot of other scripts rely on this
# create a list that maps origin numbers (params) to what they actually are
natal_origins <- gsub(" Mouth| Upstream", "", model_states)
natal_origins <- natal_origins[!(duplicated(natal_origins))]
natal_origins <- natal_origins[!(grepl("mainstem", natal_origins))]
natal_origins <- natal_origins[!(grepl("other tributaries", natal_origins))]
natal_origins <- natal_origins[!(natal_origins == "loss")]

# Use the parameter map to index the right effects
origin_param_map <- data.frame(
  natal_origin = natal_origins,
  hatchery = c(NA, NA, NA, 1, NA, 2, # MC
               1,NA,2,3, # UC
               5,NA,1,4,2,3), # SR,
  wild = c(1,3,2,4,6,5, # MC
           1,2,NA,3, # UC
           6,1,2,5,3,4)) # SR

# Add the DPS
DPS_map <- data.frame(
  natal_origin = c("Middle Columbia", "Upper Columbia", "Snake River"),
  hatchery = rep(999, 3),
  wild = rep(999, 3)
)
origin_param_map %>% 
  bind_rows(DPS_map) -> origin_param_map
