# 07-07-01-ff_cov_functions.R

# This script contains all of the functions for running the simulations for final fates and covariates.
# This script will be sourced by the other scripts that run individual populations through these simulations.

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
  "Hood River Mouth", "Hood River Upstream",
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

# Get info about state names and numbers
from_state_number_names <- data.frame(from = seq(1,43,1), from_name = model_states)
to_state_number_names <- data.frame(to = seq(1,43,1), to_name = model_states)


# get the info on transitions
UCW_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v3", "upper_columbia_wild","UCW_transition_counts.csv"))
UCW_movements <- paste0("_", UCW_transition_counts$from, "_", UCW_transition_counts$to)
UCW_movements <- UCW_movements[!(grepl("NA", UCW_movements))]

UCH_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v3", "upper_columbia_hatchery","UCH_transition_counts.csv"))
UCH_movements <- paste0("_", UCH_transition_counts$from, "_", UCH_transition_counts$to)
UCH_movements <- UCH_movements[!(grepl("NA", UCH_movements))]

MCW_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v3", "middle_columbia_wild","MCW_transition_counts.csv"))
MCW_movements <- paste0("_", MCW_transition_counts$from, "_", MCW_transition_counts$to)
MCW_movements <- MCW_movements[!(grepl("NA", MCW_movements))]

MCH_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v3", "middle_columbia_hatchery","MCH_transition_counts.csv"))
MCH_movements <- paste0("_", MCH_transition_counts$from, "_", MCH_transition_counts$to)
MCH_movements <- MCH_movements[!(grepl("NA", MCH_movements))]

SRW_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v3", "snake_river_wild","SRW_transition_counts.csv"))
SRW_movements <- paste0("_", SRW_transition_counts$from, "_", SRW_transition_counts$to)
SRW_movements <- SRW_movements[!(grepl("NA", SRW_movements))]

SRH_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v3", "snake_river_hatchery","SRH_transition_counts.csv"))
SRH_movements <- paste0("_", SRH_transition_counts$from, "_", SRH_transition_counts$to)
SRH_movements <- SRH_movements[!(grepl("NA", SRH_movements))]


##### Load the model runs #####

# Load the model data associated with each run (necessary to load covariates)
# Store these each in an environment, because most things share names
UCW_envir <- new.env()
UCH_envir <- new.env()
MCW_envir <- new.env()
MCH_envir <- new.env()
SRW_envir <- new.env()
SRH_envir <- new.env()
load(here::here("stan_actual", "reparameterization_v3", "upper_columbia_wild", "model_data.rda"),
     envir = UCW_envir)
load(here::here("stan_actual", "reparameterization_v3", "upper_columbia_hatchery", "model_data.rda"),
     envir = UCH_envir)
load(here::here("stan_actual", "reparameterization_v3", "middle_columbia_wild", "model_data.rda"),
     envir = MCW_envir)
load(here::here("stan_actual", "reparameterization_v3", "middle_columbia_hatchery", "model_data.rda"),
     envir = MCH_envir)
load(here::here("stan_actual", "reparameterization_v3", "snake_river_wild", "model_data.rda"),
     envir = SRW_envir)
load(here::here("stan_actual", "reparameterization_v3", "snake_river_hatchery", "model_data.rda"),
     envir = SRH_envir)




# Function to bind four chains together
bind4chains <- function(chain1, chain2, chain3, chain4){
  bound_draws <- bind_draws(chain1$draws(),
                            chain2$draws(),
                            chain3$draws(),
                            chain4$draws(), along = "chain")
  
  return(bound_draws)
}

## Upper Columbia, Wild
UCW_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v3", "upper_columbia_wild", "chain1_UCW_reparam_v3_fit.rds"))
UCW_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v3", "upper_columbia_wild", "chain2_UCW_reparam_v3_fit.rds"))
UCW_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v3", "upper_columbia_wild", "chain3_UCW_reparam_v3_fit.rds"))
UCW_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v3", "upper_columbia_wild", "chain4_UCW_reparam_v3_fit.rds"))

# bind chains together
UCW_fit_raw <- bind4chains(UCW_chain1, UCW_chain2, UCW_chain3, UCW_chain4)
# thin2
thin_draws(UCW_fit_raw, thin = 2) -> UCW_fit
# summarise
UCW_fit_summary <- summarise_draws(UCW_fit)

## Upper Columbia, Hatchery
UCH_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v3", "upper_columbia_hatchery", "chain1_UCH_reparam_v3_fit.rds"))
UCH_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v3", "upper_columbia_hatchery", "chain2_UCH_reparam_v3_fit.rds"))
UCH_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v3", "upper_columbia_hatchery", "chain3_UCH_reparam_v3_fit.rds"))
UCH_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v3", "upper_columbia_hatchery", "chain4_UCH_reparam_v3_fit.rds"))

# bind chains together
UCH_fit_raw <- bind4chains(UCH_chain1, UCH_chain2, UCH_chain3, UCH_chain4)
# thin2
thin_draws(UCH_fit_raw, thin = 2) -> UCH_fit
# summarise
UCH_fit_summary <- summarise_draws(UCH_fit)

## Middle Columbia, Wild
MCW_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v3", "middle_columbia_wild", "chain1_MCW_reparam_v3_fit.rds"))
MCW_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v3", "middle_columbia_wild", "chain2_MCW_reparam_v3_fit.rds"))
MCW_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v3", "middle_columbia_wild", "chain3_MCW_reparam_v3_fit.rds"))
MCW_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v3", "middle_columbia_wild", "chain4_MCW_reparam_v3_fit.rds"))

# bind chains together
MCW_fit_raw <- bind4chains(MCW_chain1, MCW_chain2, MCW_chain3, MCW_chain4)
# thin2
thin_draws(MCW_fit_raw, thin = 2) -> MCW_fit
# summarise
MCW_fit_summary <- summarise_draws(MCW_fit)

## Middle Columbia, Hatchery
MCH_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v3", "middle_columbia_hatchery", "chain1_MCH_reparam_v3_fit.rds"))
MCH_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v3", "middle_columbia_hatchery", "chain2_MCH_reparam_v3_fit.rds"))
MCH_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v3", "middle_columbia_hatchery", "chain3_MCH_reparam_v3_fit.rds"))
MCH_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v3", "middle_columbia_hatchery", "chain4_MCH_reparam_v3_fit.rds"))

# bind chains together
MCH_fit_raw <- bind4chains(MCH_chain1, MCH_chain2, MCH_chain3, MCH_chain4)
# thin2
thin_draws(MCH_fit_raw, thin = 2) -> MCH_fit
# summarise
MCH_fit_summary <- summarise_draws(MCH_fit)

## Snake River, Wild
SRW_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v3", "snake_river_wild", "chain1_SRW_reparam_v3_fit.rds"))
SRW_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v3", "snake_river_wild", "chain2_SRW_reparam_v3_fit.rds"))
SRW_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v3", "snake_river_wild", "chain3_SRW_reparam_v3_fit.rds"))
SRW_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v3", "snake_river_wild", "chain4_SRW_reparam_v3_fit.rds"))

# bind chains together
SRW_fit_raw <- bind4chains(SRW_chain1, SRW_chain2, SRW_chain3, SRW_chain4)
# thin2
thin_draws(SRW_fit_raw, thin = 2) -> SRW_fit
# summarise
SRW_fit_summary <- summarise_draws(SRW_fit)

## Snake River, Hatchery
SRH_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v3", "snake_river_hatchery", "chain1_SRH_reparam_v3_fit.rds"))
SRH_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v3", "snake_river_hatchery", "chain2_SRH_reparam_v3_fit.rds"))
SRH_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v3", "snake_river_hatchery", "chain3_SRH_reparam_v3_fit.rds"))
SRH_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v3", "snake_river_hatchery", "chain4_SRH_reparam_v3_fit.rds"))

# bind chains together
SRH_fit_raw <- bind4chains(SRH_chain1, SRH_chain2, SRH_chain3, SRH_chain4)
# thin2
thin_draws(SRH_fit_raw, thin = 2) -> SRH_fit
# summarise
SRH_fit_summary <- summarise_draws(SRH_fit)

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

### UCW ###
UCW_parameters <- UCW_fit_summary$variable
UCW_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", UCW_parameters)] -> UCW_fixed_effects
UCW_fixed_effects[!(grepl("vector", UCW_fixed_effects))] -> UCW_fixed_effects

b0_array_UCW <- make_parameter_draws_array(parameter_prefix = "b0", fit = UCW_fit, fit_summary = UCW_fit_summary)
btemp0_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = UCW_fit, fit_summary = UCW_fit_summary)
btemp1_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = UCW_fit, fit_summary = UCW_fit_summary)
bspillwindow_array_UCW <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = UCW_fit, fit_summary = UCW_fit_summary)
bwinterspill_array_UCW <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = UCW_fit, fit_summary = UCW_fit_summary)
btemp0xorigin1_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = UCW_fit, fit_summary = UCW_fit_summary)
btemp1xorigin1_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = UCW_fit, fit_summary = UCW_fit_summary)
btemp0xorigin2_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = UCW_fit, fit_summary = UCW_fit_summary)
btemp1xorigin2_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = UCW_fit, fit_summary = UCW_fit_summary)
btemp0xorigin3_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = UCW_fit, fit_summary = UCW_fit_summary)
btemp1xorigin3_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = UCW_fit, fit_summary = UCW_fit_summary)
borigin1_array_UCW <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = UCW_fit, fit_summary = UCW_fit_summary)
borigin2_array_UCW <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = UCW_fit, fit_summary = UCW_fit_summary)
borigin3_array_UCW <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = UCW_fit, fit_summary = UCW_fit_summary)

### UCH ###
UCH_parameters <- UCH_fit_summary$variable
UCH_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", UCH_parameters)] -> UCH_fixed_effects
UCH_fixed_effects[!(grepl("vector", UCH_fixed_effects))] -> UCH_fixed_effects

b0_array_UCH <- make_parameter_draws_array(parameter_prefix = "b0", fit = UCH_fit, fit_summary = UCH_fit_summary)
btemp0_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = UCH_fit, fit_summary = UCH_fit_summary)
btemp1_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = UCH_fit, fit_summary = UCH_fit_summary)
bspillwindow_array_UCH <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = UCH_fit, fit_summary = UCH_fit_summary)
bwinterspill_array_UCH <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = UCH_fit, fit_summary = UCH_fit_summary)
btemp0xorigin1_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = UCH_fit, fit_summary = UCH_fit_summary)
btemp1xorigin1_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = UCH_fit, fit_summary = UCH_fit_summary)
btemp0xorigin2_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = UCH_fit, fit_summary = UCH_fit_summary)
btemp1xorigin2_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = UCH_fit, fit_summary = UCH_fit_summary)
btemp0xorigin3_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = UCH_fit, fit_summary = UCH_fit_summary)
btemp1xorigin3_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = UCH_fit, fit_summary = UCH_fit_summary)
borigin1_array_UCH <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = UCH_fit, fit_summary = UCH_fit_summary)
borigin2_array_UCH <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = UCH_fit, fit_summary = UCH_fit_summary)
borigin3_array_UCH <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = UCH_fit, fit_summary = UCH_fit_summary)

### MCW ###
MCW_parameters <- MCW_fit_summary$variable
MCW_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", MCW_parameters)] -> MCW_fixed_effects
MCW_fixed_effects[!(grepl("vector", MCW_fixed_effects))] -> MCW_fixed_effects

b0_array_MCW <- make_parameter_draws_array(parameter_prefix = "b0", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp0_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp1_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = MCW_fit, fit_summary = MCW_fit_summary)
bspillwindow_array_MCW <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = MCW_fit, fit_summary = MCW_fit_summary)
bwinterspill_array_MCW <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp0xorigin1_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp1xorigin1_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp0xorigin2_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp1xorigin2_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp0xorigin3_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp1xorigin3_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp0xorigin4_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin4", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp1xorigin4_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin4", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp0xorigin5_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin5", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp1xorigin5_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin5", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp0xorigin6_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin6", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp1xorigin6_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin6", fit = MCW_fit, fit_summary = MCW_fit_summary)
borigin1_array_MCW <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = MCW_fit, fit_summary = MCW_fit_summary)
borigin2_array_MCW <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = MCW_fit, fit_summary = MCW_fit_summary)
borigin3_array_MCW <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = MCW_fit, fit_summary = MCW_fit_summary)
borigin4_array_MCW <- make_parameter_draws_array(parameter_prefix = "borigin4", fit = MCW_fit, fit_summary = MCW_fit_summary)
borigin5_array_MCW <- make_parameter_draws_array(parameter_prefix = "borigin5", fit = MCW_fit, fit_summary = MCW_fit_summary)
borigin6_array_MCW <- make_parameter_draws_array(parameter_prefix = "borigin6", fit = MCW_fit, fit_summary = MCW_fit_summary)

### MCH ###
MCH_parameters <- MCH_fit_summary$variable
MCH_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", MCH_parameters)] -> MCH_fixed_effects
MCH_fixed_effects[!(grepl("vector", MCH_fixed_effects))] -> MCH_fixed_effects

b0_array_MCH <- make_parameter_draws_array(parameter_prefix = "b0", fit = MCH_fit, fit_summary = MCH_fit_summary)
btemp0_array_MCH <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = MCH_fit, fit_summary = MCH_fit_summary)
btemp1_array_MCH <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = MCH_fit, fit_summary = MCH_fit_summary)
bspillwindow_array_MCH <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = MCH_fit, fit_summary = MCH_fit_summary)
bwinterspill_array_MCH <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = MCH_fit, fit_summary = MCH_fit_summary)
btemp0xorigin1_array_MCH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = MCH_fit, fit_summary = MCH_fit_summary)
btemp1xorigin1_array_MCH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = MCH_fit, fit_summary = MCH_fit_summary)
btemp0xorigin2_array_MCH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = MCH_fit, fit_summary = MCH_fit_summary)
btemp1xorigin2_array_MCH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = MCH_fit, fit_summary = MCH_fit_summary)
borigin1_array_MCH <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = MCH_fit, fit_summary = MCH_fit_summary)
borigin2_array_MCH <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = MCH_fit, fit_summary = MCH_fit_summary)

### SRW ###
SRW_parameters <- SRW_fit_summary$variable
SRW_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", SRW_parameters)] -> SRW_fixed_effects
SRW_fixed_effects[!(grepl("vector", SRW_fixed_effects))] -> SRW_fixed_effects

b0_array_SRW <- make_parameter_draws_array(parameter_prefix = "b0", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp0_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp1_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = SRW_fit, fit_summary = SRW_fit_summary)
bspillwindow_array_SRW <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = SRW_fit, fit_summary = SRW_fit_summary)
bwinterspill_array_SRW <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp0xorigin1_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp1xorigin1_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp0xorigin2_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp1xorigin2_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp0xorigin3_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp1xorigin3_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp0xorigin4_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin4", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp1xorigin4_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin4", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp0xorigin5_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin5", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp1xorigin5_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin5", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp0xorigin6_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin6", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp1xorigin6_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin6", fit = SRW_fit, fit_summary = SRW_fit_summary)
borigin1_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = SRW_fit, fit_summary = SRW_fit_summary)
borigin2_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = SRW_fit, fit_summary = SRW_fit_summary)
borigin3_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = SRW_fit, fit_summary = SRW_fit_summary)
borigin4_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin4", fit = SRW_fit, fit_summary = SRW_fit_summary)
borigin5_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin5", fit = SRW_fit, fit_summary = SRW_fit_summary)
borigin6_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin6", fit = SRW_fit, fit_summary = SRW_fit_summary)

### SRH ###
SRH_parameters <- SRH_fit_summary$variable
SRH_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", SRH_parameters)] -> SRH_fixed_effects
SRH_fixed_effects[!(grepl("vector", SRH_fixed_effects))] -> SRH_fixed_effects

b0_array_SRH <- make_parameter_draws_array(parameter_prefix = "b0", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp0_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp1_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = SRH_fit, fit_summary = SRH_fit_summary)
bspillwindow_array_SRH <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = SRH_fit, fit_summary = SRH_fit_summary)
bwinterspill_array_SRH <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp0xorigin1_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp1xorigin1_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp0xorigin2_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp1xorigin2_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp0xorigin3_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp1xorigin3_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp0xorigin4_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin4", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp1xorigin4_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin4", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp0xorigin5_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin5", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp1xorigin5_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin5", fit = SRH_fit, fit_summary = SRH_fit_summary)
borigin1_array_SRH <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = SRH_fit, fit_summary = SRH_fit_summary)
borigin2_array_SRH <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = SRH_fit, fit_summary = SRH_fit_summary)
borigin3_array_SRH <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = SRH_fit, fit_summary = SRH_fit_summary)
borigin4_array_SRH <- make_parameter_draws_array(parameter_prefix = "borigin4", fit = SRH_fit, fit_summary = SRH_fit_summary)
borigin5_array_SRH <- make_parameter_draws_array(parameter_prefix = "borigin5", fit = SRH_fit, fit_summary = SRH_fit_summary)

#### Select representative warm vs. cold years ####

run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
run_year_start <- seq(ymd("2004-06-01"), ymd("2022-06-01"), by = "years")
run_year_end <- seq(ymd("2005-05-31"), ymd("2023-05-31"), by = "years")
run_year_numeric = seq(4, 22, 1)

run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)

temp_ts <- as.data.frame(UCW_envir$data$temperature_data)
dates <- seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days")
years <- year(dates)
temp_ts$date <- dates
temp_ts$year <- years

# Get the median for each run year for both winter/spring (temp0) and summer/fall (temp1)

temp_ts %>% 
  rowwise() %>%
  dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date & run_year_end >= date)$run_year) %>% 
  dplyr::mutate(season = ifelse(yday(date)<152, 0, 1))  %>% 
  # drop 22/23 because it's not a complete year (and it isn't modeled)
  filter(!(run_year == "22/23")) %>% 
  subset(season == 0) %>% 
  group_by(run_year) %>% 
  summarise_all(median) %>% 
  arrange(BON) -> season0_annual_temp_medians

temp_ts %>% 
  rowwise() %>%
  dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date & run_year_end >= date)$run_year) %>% 
  dplyr::mutate(season = ifelse(yday(date)<152, 0, 1))  %>% 
  # drop 22/23 because it's not a complete year (and it isn't modeled)
  filter(!(run_year == "22/23")) %>% 
  subset(season == 1) %>% 
  group_by(run_year) %>% 
  summarise_all(median) %>% 
  arrange(BON) -> season1_annual_temp_medians

spill_ts <- as.data.frame(UCW_envir$data$spill_window_data)
dates <- seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days")
years <- year(dates)
spill_ts$date <- dates
spill_ts$year <- years

spill_ts %>% 
  rowwise() %>%
  dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date & run_year_end >= date)$run_year) %>% 
  dplyr::mutate(season = ifelse(yday(date)<152, 0, 1))  %>% 
  # drop 22/23 because it's not a complete year (and it isn't modeled)
  filter(!(run_year == "22/23")) %>% 
  subset(season == 0) %>% 
  group_by(run_year) %>% 
  summarise_all(median) %>% 
  arrange(BON) -> season0_annual_spill_medians

spill_ts %>% 
  rowwise() %>%
  dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date & run_year_end >= date)$run_year) %>% 
  dplyr::mutate(season = ifelse(yday(date)<152, 0, 1))  %>% 
  # drop 22/23 because it's not a complete year (and it isn't modeled)
  filter(!(run_year == "22/23")) %>% 
  subset(season == 1) %>% 
  group_by(run_year) %>% 
  summarise_all(median) %>% 
  arrange(BON) -> season1_annual_spill_medians

ggplot(season0_annual_temp_medians, aes(x = year, y = BON)) +
  geom_point(color = "blue") +
  geom_point(data = season1_annual_temp_medians, aes(x = year, y = BON), color = "red")

ggplot(season0_annual_spill_medians, aes(x = year, y = BON)) +
  geom_point(color = "blue") +
  geom_point(data = season1_annual_spill_medians, aes(x = year, y = BON), color = "red")

plot(x = season1_annual_temp_medians$BON, y = season1_annual_spill_medians$BON)
# spill volume and temperature are clearly correlated: lower temperature = lower spill
# So you have to plot spill and temp together

# We'll judge warm v cool based on summer/fall temps, since that's when the majority of movements are happening
# pick out median, 25%, 75%
summary(1:17)
season1_annual_temp_medians[c(1,5,9,13,17),]

# coldest year on record: 11/12
# cool year (25%): 09/10
# average year (50%): 05/06
# warm year (75%): 18/19
# hottest year on record: 15/16

rep_years <- data.frame(fix_run_year = season1_annual_temp_medians[c(1,5,9,13,17),]$run_year,
                        conditions = c("coldest", "cool", "average", "warm", "warmest"))





#### Final fates functions under different covariate values ####

# This function takes a model fit object, a number of simulated fish,
# the identities of those fish (origin, rear), the conditions throughout the
# basin (spill, temperature, year)
# and then simulates the final fates of those fish
# By default, the conditions (spill and temperature) will be taken from the data
# itself, but note that new values of these conditions could be simulated in order

# Note that another option here would be to use average (median) conditions experienced
# in each state - and that would remove the variability associated with taking random 
# conditions. But note that you'll have to decide what the data is that you're taking
# the median of - is it the whole DPS? Just the population? Currently, we have it set up
# to run by DPS, not population
# Another option would be to make the number of fish per simulation smaller and
# run it more times - that would give you a wider spread of covariate values
# and therefore reduce the chance that a handful of outlier covariate values
# lead to crazy distributions


# Final states simulation under different covariate values

# this function will be very similar to the other final fates simulations, but
# with additional arguments that allow you to fix conditions in certain reaches




# This function will have to be re-written for each DPS, because each DPS has different params
# currently we are not including a random effect of year, but we could easily
# amend the code to simulate final fates in a particular year
# origin1/origin2/origin3 are indicator variables, so set the origin you want to be 1

final_fates_simulation_fixcov_UCW <- function(nsim,
                                              start_state = 2, states_dates,
                                              origin1 = 0, origin2 = 0, origin3 = 0,
                                              temp_data, spillwindow_data, winterspill_data,
                                              fix_state, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                              fix_spillwindow_season0_value = NA, fix_spillwindow_season1_value = NA,
                                              fix_winterspill_value = NA,
                                              fix_run_year = NA){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(UCW_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state across all years
    winterspill[i] <- median(winterspill_data[,i])
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(UCW_states_dates, state == i  & direction == "upstream")$season)/length(subset(UCW_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(UCW_states_dates, state == i  & direction == "downstream")$season)/length(subset(UCW_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    # get the date numeric for each day/month of that run year
    run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
    run_year_start <- seq(ymd("2004-06-01"), ymd("2022-06-01"), by = "years")
    run_year_end <- seq(ymd("2005-05-31"), ymd("2023-05-31"), by = "years")
    run_year_numeric = seq(4, 22, 1)
    
    run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)
    
    
    
    # rowwise() %>%
    # dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_actual & run_year_end >= date_actual)$run_year) %>% 
    # filter(run_year == fix_run_year)-> UCW_states_dates_run_year
    
    # subset states_dates for origin and population to pull out just month day, and then
    # use those months and days to index temps from that year
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    dates_actual %>% 
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) -> dates_actual
    
    
    subset(run_year_df, run_year == fix_run_year) -> fix_run_year_df
    subset(dates_actual, date_actual <= fix_run_year_df$run_year_end & date_actual >= fix_run_year_df$run_year_start) %>% 
      dplyr::rename(date_run_year = date) -> dates_fix_run_year
    
    # Deal with leap years by using 2/28 twice
    if(!("02-29" %in% dates_fix_run_year$month_day)){
      subset(dates_fix_run_year, month_day == "02-28") %>% 
        mutate(month_day = "02-29") -> leap_year
      
      dates_fix_run_year %>% 
        bind_rows(., leap_year) -> dates_fix_run_year
      
    }
    
    
    
    
    
    UCW_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> UCW_states_month_day
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> states_month_day
    
    
    temp_upstream_season0 <- rep(0, length(model_states))
    temp_upstream_season1 <- rep(0, length(model_states))
    temp_downstream_season0 <- rep(0, length(model_states))
    temp_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(UCW_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(UCW_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(UCW_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(UCW_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(UCW_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(UCW_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(UCW_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(UCW_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    spillwindow_upstream_season0 <- rep(0, length(model_states))
    spillwindow_upstream_season1 <- rep(0, length(model_states))
    spillwindow_downstream_season0 <- rep(0, length(model_states))
    spillwindow_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(UCW_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(UCW_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(UCW_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(UCW_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(UCW_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(UCW_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(UCW_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(UCW_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  
  # Fix covariate values to the ones that you are setting
  if(!(is.na(fix_temp_season0_value))){
    temp_upstream_season0[fix_state] <- fix_temp_season0_value
    temp_downstream_season0[fix_state] <- fix_temp_season0_value
  }
  
  if(!(is.na(fix_temp_season1_value))){
    temp_upstream_season1[fix_state] <- fix_temp_season1_value
    temp_downstream_season1[fix_state] <- fix_temp_season1_value
  }
  
  if(!(is.na(fix_spillwindow_season0_value))){
    spillwindow_upstream_season0[fix_state] <- fix_spillwindow_season0_value
    spillwindow_downstream_season0[fix_state] <- fix_spillwindow_season0_value
  }
  
  if(!(is.na(fix_spillwindow_season1_value))){
    spillwindow_upstream_season1[fix_state] <- fix_spillwindow_season1_value
    spillwindow_downstream_season1[fix_state] <- fix_spillwindow_season1_value
  }
  
  if(!(is.na(fix_winterspill_value))){
    winterspill[fix_state] <- fix_winterspill_value
  }
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[j], season_upstream[j]))
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCW[i,possible_movements[j],iter] +
                                                                    btemp1_array_UCW[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_UCW[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_UCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_UCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_UCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_UCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    borigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_UCW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_UCW[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCW[i,possible_movements,iter] +
                    btemp1_array_UCW[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_UCW[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_UCW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_UCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_UCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_UCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    borigin1_array_UCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCW[i,possible_movements,iter]*origin3))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCW[i,possible_movements[j],iter] +
                                                                    btemp0_array_UCW[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_UCW[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_UCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_UCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_UCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_UCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    borigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_UCW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_UCW[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCW[i,possible_movements,iter] +
                    btemp0_array_UCW[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_UCW[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_UCW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_UCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_UCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_UCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    borigin1_array_UCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCW[i,possible_movements,iter]*origin3))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[j], season_downstream[j]))
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCW[i,possible_movements[j],iter] +
                                                                      btemp1_array_UCW[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_UCW[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_UCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_UCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_UCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_UCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      borigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_UCW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_UCW[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCW[i,possible_movements,iter] +
                    btemp1_array_UCW[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_UCW[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_UCW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_UCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_UCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_UCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    borigin1_array_UCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCW[i,possible_movements,iter]*origin3))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCW[i,possible_movements[j],iter] +
                                                                      btemp0_array_UCW[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_UCW[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_UCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_UCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_UCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_UCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      borigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_UCW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_UCW[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCW[i,possible_movements,iter] +
                    btemp0_array_UCW[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_UCW[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_UCW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_UCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_UCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_UCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    borigin1_array_UCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCW[i,possible_movements,iter]*origin3))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_fixcov_UCH <- function(nsim,
                                              start_state = 2, states_dates,
                                              origin1 = 0, origin2 = 0, origin3 = 0,
                                              temp_data, spillwindow_data, winterspill_data,
                                              fix_state, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                              fix_spillwindow_season0_value = NA, fix_spillwindow_season1_value = NA,
                                              fix_winterspill_value = NA,
                                              fix_run_year = NA){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(UCH_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state across all years
    winterspill[i] <- median(winterspill_data[,i])
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(UCH_states_dates, state == i  & direction == "upstream")$season)/length(subset(UCH_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(UCH_states_dates, state == i  & direction == "downstream")$season)/length(subset(UCH_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    # get the date numeric for each day/month of that run year
    run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
    run_year_start <- seq(ymd("2004-06-01"), ymd("2022-06-01"), by = "years")
    run_year_end <- seq(ymd("2005-05-31"), ymd("2023-05-31"), by = "years")
    run_year_numeric = seq(4, 22, 1)
    
    run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)
    
    
    
    # rowwise() %>%
    # dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_actual & run_year_end >= date_actual)$run_year) %>% 
    # filter(run_year == fix_run_year)-> UCH_states_dates_run_year
    
    # subset states_dates for origin and population to pull out just month day, and then
    # use those months and days to index temps from that year
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    dates_actual %>% 
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) -> dates_actual
    
    
    subset(run_year_df, run_year == fix_run_year) -> fix_run_year_df
    subset(dates_actual, date_actual <= fix_run_year_df$run_year_end & date_actual >= fix_run_year_df$run_year_start) %>% 
      dplyr::rename(date_run_year = date) -> dates_fix_run_year
    
    # Deal with leap years by using 2/28 twice
    if(!("02-29" %in% dates_fix_run_year$month_day)){
      subset(dates_fix_run_year, month_day == "02-28") %>% 
        mutate(month_day = "02-29") -> leap_year
      
      dates_fix_run_year %>% 
        bind_rows(., leap_year) -> dates_fix_run_year
      
    }
    
    
    
    
    
    UCH_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> UCH_states_month_day
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> states_month_day
    
    
    temp_upstream_season0 <- rep(0, length(model_states))
    temp_upstream_season1 <- rep(0, length(model_states))
    temp_downstream_season0 <- rep(0, length(model_states))
    temp_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(UCH_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(UCH_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(UCH_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(UCH_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(UCH_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(UCH_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(UCH_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(UCH_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    spillwindow_upstream_season0 <- rep(0, length(model_states))
    spillwindow_upstream_season1 <- rep(0, length(model_states))
    spillwindow_downstream_season0 <- rep(0, length(model_states))
    spillwindow_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(UCH_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(UCH_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(UCH_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(UCH_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(UCH_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(UCH_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(UCH_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(UCH_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  
  # Fix covariate values to the ones that you are setting
  if(!(is.na(fix_temp_season0_value))){
    temp_upstream_season0[fix_state] <- fix_temp_season0_value
    temp_downstream_season0[fix_state] <- fix_temp_season0_value
  }
  
  if(!(is.na(fix_temp_season1_value))){
    temp_upstream_season1[fix_state] <- fix_temp_season1_value
    temp_downstream_season1[fix_state] <- fix_temp_season1_value
  }
  
  if(!(is.na(fix_spillwindow_season0_value))){
    spillwindow_upstream_season0[fix_state] <- fix_spillwindow_season0_value
    spillwindow_downstream_season0[fix_state] <- fix_spillwindow_season0_value
  }
  
  if(!(is.na(fix_spillwindow_season1_value))){
    spillwindow_upstream_season1[fix_state] <- fix_spillwindow_season1_value
    spillwindow_downstream_season1[fix_state] <- fix_spillwindow_season1_value
  }
  
  if(!(is.na(fix_winterspill_value))){
    winterspill[fix_state] <- fix_winterspill_value
  }
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[j], season_upstream[j]))
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCH[i,possible_movements[j],iter] +
                                                                    btemp1_array_UCH[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_UCH[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_UCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_UCH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_UCH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_UCH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    borigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_UCH[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_UCH[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCH[i,possible_movements,iter] +
                    btemp1_array_UCH[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_UCH[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_UCH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_UCH[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_UCH[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_UCH[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    borigin1_array_UCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCH[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCH[i,possible_movements,iter]*origin3))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCH[i,possible_movements[j],iter] +
                                                                    btemp0_array_UCH[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_UCH[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_UCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_UCH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_UCH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_UCH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    borigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_UCH[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_UCH[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCH[i,possible_movements,iter] +
                    btemp0_array_UCH[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_UCH[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_UCH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_UCH[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_UCH[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_UCH[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    borigin1_array_UCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCH[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCH[i,possible_movements,iter]*origin3))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[j], season_downstream[j]))
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCH[i,possible_movements[j],iter] +
                                                                      btemp1_array_UCH[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_UCH[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_UCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_UCH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_UCH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_UCH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      borigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_UCH[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_UCH[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCH[i,possible_movements,iter] +
                    btemp1_array_UCH[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_UCH[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_UCH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_UCH[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_UCH[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_UCH[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    borigin1_array_UCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCH[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCH[i,possible_movements,iter]*origin3))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCH[i,possible_movements[j],iter] +
                                                                      btemp0_array_UCH[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_UCH[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_UCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_UCH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_UCH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_UCH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      borigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_UCH[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_UCH[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCH[i,possible_movements,iter] +
                    btemp0_array_UCH[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_UCH[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_UCH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_UCH[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_UCH[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_UCH[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    borigin1_array_UCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCH[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCH[i,possible_movements,iter]*origin3))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_fixcov_MCW <- function(nsim,
                                              start_state = 2, states_dates,
                                              origin1 = 0, origin2 = 0, origin3 = 0,
                                              origin4 = 0, origin5 = 0, origin6 = 0,
                                              temp_data, spillwindow_data, winterspill_data,
                                              fix_state, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                              fix_spillwindow_season0_value = NA, fix_spillwindow_season1_value = NA,
                                              fix_winterspill_value = NA,
                                              fix_run_year = NA){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(MCW_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state across all years
    winterspill[i] <- median(winterspill_data[,i])
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(MCW_states_dates, state == i  & direction == "upstream")$season)/length(subset(MCW_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(MCW_states_dates, state == i  & direction == "downstream")$season)/length(subset(MCW_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    # get the date numeric for each day/month of that run year
    run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
    run_year_start <- seq(ymd("2004-06-01"), ymd("2022-06-01"), by = "years")
    run_year_end <- seq(ymd("2005-05-31"), ymd("2023-05-31"), by = "years")
    run_year_numeric = seq(4, 22, 1)
    
    run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)
    
    
    
    # rowwise() %>%
    # dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_actual & run_year_end >= date_actual)$run_year) %>% 
    # filter(run_year == fix_run_year)-> MCW_states_dates_run_year
    
    # subset states_dates for origin and population to pull out just month day, and then
    # use those months and days to index temps from that year
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    dates_actual %>% 
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) -> dates_actual
    
    
    subset(run_year_df, run_year == fix_run_year) -> fix_run_year_df
    subset(dates_actual, date_actual <= fix_run_year_df$run_year_end & date_actual >= fix_run_year_df$run_year_start) %>% 
      dplyr::rename(date_run_year = date) -> dates_fix_run_year
    
    # Deal with leap years by using 2/28 twice
    if(!("02-29" %in% dates_fix_run_year$month_day)){
      subset(dates_fix_run_year, month_day == "02-28") %>% 
        mutate(month_day = "02-29") -> leap_year
      
      dates_fix_run_year %>% 
        bind_rows(., leap_year) -> dates_fix_run_year
      
    }
    
    
    
    
    
    MCW_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> MCW_states_month_day
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> states_month_day
    
    
    temp_upstream_season0 <- rep(0, length(model_states))
    temp_upstream_season1 <- rep(0, length(model_states))
    temp_downstream_season0 <- rep(0, length(model_states))
    temp_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(MCW_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(MCW_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(MCW_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(MCW_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    spillwindow_upstream_season0 <- rep(0, length(model_states))
    spillwindow_upstream_season1 <- rep(0, length(model_states))
    spillwindow_downstream_season0 <- rep(0, length(model_states))
    spillwindow_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(MCW_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(MCW_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(MCW_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(MCW_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  
  # Fix covariate values to the ones that you are setting
  if(!(is.na(fix_temp_season0_value))){
    temp_upstream_season0[fix_state] <- fix_temp_season0_value
    temp_downstream_season0[fix_state] <- fix_temp_season0_value
  }
  
  if(!(is.na(fix_temp_season1_value))){
    temp_upstream_season1[fix_state] <- fix_temp_season1_value
    temp_downstream_season1[fix_state] <- fix_temp_season1_value
  }
  
  if(!(is.na(fix_spillwindow_season0_value))){
    spillwindow_upstream_season0[fix_state] <- fix_spillwindow_season0_value
    spillwindow_downstream_season0[fix_state] <- fix_spillwindow_season0_value
  }
  
  if(!(is.na(fix_spillwindow_season1_value))){
    spillwindow_upstream_season1[fix_state] <- fix_spillwindow_season1_value
    spillwindow_downstream_season1[fix_state] <- fix_spillwindow_season1_value
  }
  
  if(!(is.na(fix_winterspill_value))){
    winterspill[fix_state] <- fix_winterspill_value
  }
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[j], season_upstream[j]))
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCW[i,possible_movements[j],iter] +
                                                                    btemp1_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_MCW[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_MCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    btemp1xorigin6_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin6 +
                                                                    borigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_MCW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_MCW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_MCW[i,possible_movements,iter] +
                    btemp1_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_MCW[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_MCW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    btemp1xorigin6_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin6 +
                    borigin1_array_MCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_MCW[i,possible_movements,iter]*origin3 +
                    borigin4_array_MCW[i,possible_movements,iter]*origin4 +
                    borigin5_array_MCW[i,possible_movements,iter]*origin5 +
                    borigin6_array_MCW[i,possible_movements,iter]*origin6))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCW[i,possible_movements[j],iter] +
                                                                    btemp0_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_MCW[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_MCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    btemp0xorigin6_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin6 +
                                                                    borigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_MCW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_MCW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_MCW[i,possible_movements,iter] +
                    btemp0_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_MCW[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_MCW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    btemp0xorigin6_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin6 +
                    borigin1_array_MCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_MCW[i,possible_movements,iter]*origin3 +
                    borigin4_array_MCW[i,possible_movements,iter]*origin4 +
                    borigin5_array_MCW[i,possible_movements,iter]*origin5 +
                    borigin6_array_MCW[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[j], season_downstream[j]))
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCW[i,possible_movements[j],iter] +
                                                                      btemp1_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_MCW[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_MCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      btemp1xorigin6_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin6 +
                                                                      borigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_MCW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_MCW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_MCW[i,possible_movements,iter] +
                    btemp1_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_MCW[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_MCW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    btemp1xorigin6_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin6 +
                    borigin1_array_MCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_MCW[i,possible_movements,iter]*origin3 +
                    borigin4_array_MCW[i,possible_movements,iter]*origin4 +
                    borigin5_array_MCW[i,possible_movements,iter]*origin5 +
                    borigin6_array_MCW[i,possible_movements,iter]*origin6))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCW[i,possible_movements[j],iter] +
                                                                      btemp0_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_MCW[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_MCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      btemp0xorigin6_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin6 +
                                                                      borigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_MCW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_MCW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_MCW[i,possible_movements,iter] +
                    btemp0_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_MCW[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_MCW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    btemp0xorigin6_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin6 +
                    borigin1_array_MCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_MCW[i,possible_movements,iter]*origin3 +
                    borigin4_array_MCW[i,possible_movements,iter]*origin4 +
                    borigin5_array_MCW[i,possible_movements,iter]*origin5 +
                    borigin6_array_MCW[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_fixcov_MCH <- function(nsim,
                                              start_state = 2, states_dates,
                                              origin1 = 0, origin2 = 0, 
                                              temp_data, spillwindow_data, winterspill_data,
                                              fix_state, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                              fix_spillwindow_season0_value = NA, fix_spillwindow_season1_value = NA,
                                              fix_winterspill_value = NA,
                                              fix_run_year = NA){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(MCH_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # approach 2: Take the median conditions for each state across all years
    winterspill[i] <- median(winterspill_data[,i])
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(MCH_states_dates, state == i  & direction == "upstream")$season)/length(subset(MCH_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(MCH_states_dates, state == i  & direction == "downstream")$season)/length(subset(MCH_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    # get the date numeric for each day/month of that run year
    run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
    run_year_start <- seq(ymd("2004-06-01"), ymd("2022-06-01"), by = "years")
    run_year_end <- seq(ymd("2005-05-31"), ymd("2023-05-31"), by = "years")
    run_year_numeric = seq(4, 22, 1)
    
    run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)
    
    
    
    # rowwise() %>%
    # dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_actual & run_year_end >= date_actual)$run_year) %>% 
    # filter(run_year == fix_run_year)-> MCH_states_dates_run_year
    
    # subset states_dates for origin and population to pull out just month day, and then
    # use those months and days to index temps from that year
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    dates_actual %>% 
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) -> dates_actual
    
    
    subset(run_year_df, run_year == fix_run_year) -> fix_run_year_df
    subset(dates_actual, date_actual <= fix_run_year_df$run_year_end & date_actual >= fix_run_year_df$run_year_start) %>% 
      dplyr::rename(date_run_year = date) -> dates_fix_run_year
    
    # Deal with leap years by using 2/28 twice
    if(!("02-29" %in% dates_fix_run_year$month_day)){
      subset(dates_fix_run_year, month_day == "02-28") %>% 
        mutate(month_day = "02-29") -> leap_year
      
      dates_fix_run_year %>% 
        bind_rows(., leap_year) -> dates_fix_run_year
      
    }
    
    
    
    
    
    MCH_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> MCH_states_month_day
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> states_month_day
    
    
    temp_upstream_season0 <- rep(0, length(model_states))
    temp_upstream_season1 <- rep(0, length(model_states))
    temp_downstream_season0 <- rep(0, length(model_states))
    temp_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(MCH_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(MCH_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(MCH_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(MCH_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    spillwindow_upstream_season0 <- rep(0, length(model_states))
    spillwindow_upstream_season1 <- rep(0, length(model_states))
    spillwindow_downstream_season0 <- rep(0, length(model_states))
    spillwindow_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(MCH_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(MCH_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(MCH_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(MCH_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  
  # Fix covariate values to the ones that you are setting
  if(!(is.na(fix_temp_season0_value))){
    temp_upstream_season0[fix_state] <- fix_temp_season0_value
    temp_downstream_season0[fix_state] <- fix_temp_season0_value
  }
  
  if(!(is.na(fix_temp_season1_value))){
    temp_upstream_season1[fix_state] <- fix_temp_season1_value
    temp_downstream_season1[fix_state] <- fix_temp_season1_value
  }
  
  if(!(is.na(fix_spillwindow_season0_value))){
    spillwindow_upstream_season0[fix_state] <- fix_spillwindow_season0_value
    spillwindow_downstream_season0[fix_state] <- fix_spillwindow_season0_value
  }
  
  if(!(is.na(fix_spillwindow_season1_value))){
    spillwindow_upstream_season1[fix_state] <- fix_spillwindow_season1_value
    spillwindow_downstream_season1[fix_state] <- fix_spillwindow_season1_value
  }
  
  if(!(is.na(fix_winterspill_value))){
    winterspill[fix_state] <- fix_winterspill_value
  }
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[j], season_upstream[j]))
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCH[i,possible_movements[j],iter] +
                                                                    btemp1_array_MCH[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_MCH[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_MCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_MCH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_MCH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    borigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_MCH[i,possible_movements[j],iter]*origin2)/
          sum(exp(b0_array_MCH[i,possible_movements,iter] +
                    btemp1_array_MCH[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_MCH[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_MCH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_MCH[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_MCH[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    borigin1_array_MCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCH[i,possible_movements,iter]*origin2))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCH[i,possible_movements[j],iter] +
                                                                    btemp0_array_MCH[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_MCH[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_MCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_MCH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_MCH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    borigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_MCH[i,possible_movements[j],iter]*origin2)/
          sum(exp(b0_array_MCH[i,possible_movements,iter] +
                    btemp0_array_MCH[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_MCH[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_MCH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_MCH[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_MCH[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    borigin1_array_MCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCH[i,possible_movements,iter]*origin2))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[j], season_downstream[j]))
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCH[i,possible_movements[j],iter] +
                                                                      btemp1_array_MCH[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_MCH[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_MCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_MCH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_MCH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      borigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_MCH[i,possible_movements[j],iter]*origin2)/
          sum(exp(b0_array_MCH[i,possible_movements,iter] +
                    btemp1_array_MCH[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_MCH[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_MCH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_MCH[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_MCH[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    borigin1_array_MCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCH[i,possible_movements,iter]*origin2))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCH[i,possible_movements[j],iter] +
                                                                      btemp0_array_MCH[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_MCH[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_MCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_MCH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_MCH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      borigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_MCH[i,possible_movements[j],iter]*origin2)/
          sum(exp(b0_array_MCH[i,possible_movements,iter] +
                    btemp0_array_MCH[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_MCH[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_MCH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_MCH[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_MCH[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    borigin1_array_MCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCH[i,possible_movements,iter]*origin2))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_fixcov_SRW <- function(nsim,
                                              start_state = 2, states_dates,
                                              origin1 = 0, origin2 = 0, origin3 = 0,
                                              origin4 = 0, origin5 = 0, origin6 = 0,
                                              temp_data, spillwindow_data, winterspill_data,
                                              fix_state, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                              fix_spillwindow_season0_value = NA, fix_spillwindow_season1_value = NA,
                                              fix_winterspill_value = NA,
                                              fix_run_year = NA){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(SRW_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state across all years
    winterspill[i] <- median(winterspill_data[,i])
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRW_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRW_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRW_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRW_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    # get the date numeric for each day/month of that run year
    run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
    run_year_start <- seq(ymd("2004-06-01"), ymd("2022-06-01"), by = "years")
    run_year_end <- seq(ymd("2005-05-31"), ymd("2023-05-31"), by = "years")
    run_year_numeric = seq(4, 22, 1)
    
    run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)
    
    
    
    # rowwise() %>%
    # dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_actual & run_year_end >= date_actual)$run_year) %>% 
    # filter(run_year == fix_run_year)-> SRW_states_dates_run_year
    
    # subset states_dates for origin and population to pull out just month day, and then
    # use those months and days to index temps from that year
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    dates_actual %>% 
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) -> dates_actual
    
    
    subset(run_year_df, run_year == fix_run_year) -> fix_run_year_df
    subset(dates_actual, date_actual <= fix_run_year_df$run_year_end & date_actual >= fix_run_year_df$run_year_start) %>% 
      dplyr::rename(date_run_year = date) -> dates_fix_run_year
    
    # Deal with leap years by using 2/28 twice
    if(!("02-29" %in% dates_fix_run_year$month_day)){
      subset(dates_fix_run_year, month_day == "02-28") %>% 
        mutate(month_day = "02-29") -> leap_year
      
      dates_fix_run_year %>% 
        bind_rows(., leap_year) -> dates_fix_run_year
      
    }
    
    
    
    
    
    SRW_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> SRW_states_month_day
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> states_month_day
    
    
    temp_upstream_season0 <- rep(0, length(model_states))
    temp_upstream_season1 <- rep(0, length(model_states))
    temp_downstream_season0 <- rep(0, length(model_states))
    temp_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRW_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRW_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRW_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRW_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRW_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    spillwindow_upstream_season0 <- rep(0, length(model_states))
    spillwindow_upstream_season1 <- rep(0, length(model_states))
    spillwindow_downstream_season0 <- rep(0, length(model_states))
    spillwindow_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRW_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRW_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRW_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRW_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRW_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  
  # Fix covariate values to the ones that you are setting
  if(!(is.na(fix_temp_season0_value))){
    temp_upstream_season0[fix_state] <- fix_temp_season0_value
    temp_downstream_season0[fix_state] <- fix_temp_season0_value
  }
  
  if(!(is.na(fix_temp_season1_value))){
    temp_upstream_season1[fix_state] <- fix_temp_season1_value
    temp_downstream_season1[fix_state] <- fix_temp_season1_value
  }
  
  if(!(is.na(fix_spillwindow_season0_value))){
    spillwindow_upstream_season0[fix_state] <- fix_spillwindow_season0_value
    spillwindow_downstream_season0[fix_state] <- fix_spillwindow_season0_value
  }
  
  if(!(is.na(fix_spillwindow_season1_value))){
    spillwindow_upstream_season1[fix_state] <- fix_spillwindow_season1_value
    spillwindow_downstream_season1[fix_state] <- fix_spillwindow_season1_value
  }
  
  if(!(is.na(fix_winterspill_value))){
    winterspill[fix_state] <- fix_winterspill_value
  }
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[j], season_upstream[j]))
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRW[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    btemp1xorigin6_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin6 +
                                                                    borigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW[i,possible_movements,iter] +
                    btemp1_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRW[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin6 +
                    borigin1_array_SRW[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW[i,possible_movements,iter]*origin6))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRW[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    btemp0xorigin6_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin6 +
                                                                    borigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW[i,possible_movements,iter] +
                    btemp0_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRW[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin6 +
                    borigin1_array_SRW[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[j], season_downstream[j]))
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRW[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      btemp1xorigin6_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin6 +
                                                                      borigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW[i,possible_movements,iter] +
                    btemp1_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRW[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin6 +
                    borigin1_array_SRW[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW[i,possible_movements,iter]*origin6))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRW[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      btemp0xorigin6_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin6 +
                                                                      borigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW[i,possible_movements,iter] +
                    btemp0_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRW[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin6 +
                    borigin1_array_SRW[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_fixcov_SRH <- function(nsim,
                                              start_state = 2, states_dates,
                                              origin1 = 0, origin2 = 0, origin3 = 0,
                                              origin4 = 0, origin5 = 0,
                                              temp_data, spillwindow_data, winterspill_data,
                                              fix_state, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                              fix_spillwindow_season0_value = NA, fix_spillwindow_season1_value = NA,
                                              fix_winterspill_value = NA,
                                              fix_run_year = NA){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(SRH_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state across all years
    winterspill[i] <- median(winterspill_data[,i])
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRH_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRH_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRH_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRH_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    # get the date numeric for each day/month of that run year
    run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
    run_year_start <- seq(ymd("2004-06-01"), ymd("2022-06-01"), by = "years")
    run_year_end <- seq(ymd("2005-05-31"), ymd("2023-05-31"), by = "years")
    run_year_numeric = seq(4, 22, 1)
    
    run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)
    
    
    
    # rowwise() %>%
    # dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_actual & run_year_end >= date_actual)$run_year) %>% 
    # filter(run_year == fix_run_year)-> SRH_states_dates_run_year
    
    # subset states_dates for origin and population to pull out just month day, and then
    # use those months and days to index temps from that year
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    dates_actual %>% 
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) -> dates_actual
    
    
    subset(run_year_df, run_year == fix_run_year) -> fix_run_year_df
    subset(dates_actual, date_actual <= fix_run_year_df$run_year_end & date_actual >= fix_run_year_df$run_year_start) %>% 
      dplyr::rename(date_run_year = date) -> dates_fix_run_year
    
    # Deal with leap years by using 2/28 twice
    if(!("02-29" %in% dates_fix_run_year$month_day)){
      subset(dates_fix_run_year, month_day == "02-28") %>% 
        mutate(month_day = "02-29") -> leap_year
      
      dates_fix_run_year %>% 
        bind_rows(., leap_year) -> dates_fix_run_year
      
    }
    
    
    
    
    
    SRH_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> SRH_states_month_day
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> states_month_day
    
    
    temp_upstream_season0 <- rep(0, length(model_states))
    temp_upstream_season1 <- rep(0, length(model_states))
    temp_downstream_season0 <- rep(0, length(model_states))
    temp_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRH_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRH_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRH_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRH_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRH_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    spillwindow_upstream_season0 <- rep(0, length(model_states))
    spillwindow_upstream_season1 <- rep(0, length(model_states))
    spillwindow_downstream_season0 <- rep(0, length(model_states))
    spillwindow_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRH_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRH_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRH_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRH_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRH_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  
  # Fix covariate values to the ones that you are setting
  if(!(is.na(fix_temp_season0_value))){
    temp_upstream_season0[fix_state] <- fix_temp_season0_value
    temp_downstream_season0[fix_state] <- fix_temp_season0_value
  }
  
  if(!(is.na(fix_temp_season1_value))){
    temp_upstream_season1[fix_state] <- fix_temp_season1_value
    temp_downstream_season1[fix_state] <- fix_temp_season1_value
  }
  
  if(!(is.na(fix_spillwindow_season0_value))){
    spillwindow_upstream_season0[fix_state] <- fix_spillwindow_season0_value
    spillwindow_downstream_season0[fix_state] <- fix_spillwindow_season0_value
  }
  
  if(!(is.na(fix_spillwindow_season1_value))){
    spillwindow_upstream_season1[fix_state] <- fix_spillwindow_season1_value
    spillwindow_downstream_season1[fix_state] <- fix_spillwindow_season1_value
  }
  
  if(!(is.na(fix_winterspill_value))){
    winterspill[fix_state] <- fix_winterspill_value
  }
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[j], season_upstream[j]))
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRH[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    borigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH[i,possible_movements,iter] +
                    btemp1_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRH[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    borigin1_array_SRH[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH[i,possible_movements,iter]*origin5))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRH[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    borigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH[i,possible_movements,iter] +
                    btemp0_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRH[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    borigin1_array_SRH[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH[i,possible_movements,iter]*origin5))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[j], season_downstream[j]))
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRH[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      borigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH[i,possible_movements,iter] +
                    btemp1_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRH[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    borigin1_array_SRH[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH[i,possible_movements,iter]*origin5))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRH[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      borigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH[i,possible_movements,iter] +
                    btemp0_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRH[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    borigin1_array_SRH[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH[i,possible_movements,iter]*origin5))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}


#### Run final fates simulation - H v W comparisons ####


# In order to simulate covariate values, we are going to determine the dates
# where the fish were in each state, and then sample from those dates to 
# get a representative value for temperature/spill when fish are in those states
get_states_dates_direction <- function(envir){
  transitions <- envir$data$y
  transition_dates <- envir$data$transition_dates
  # convert the seasons vector to the same shape as the state transitions
  transition_seasons <- envir$data$transition_seasons_vector
  transition_seasons_matrix <- matrix(nrow = nrow(transitions), ncol = ncol(transitions))
  
  movements_counter <- 0
  for (i in 1:nrow(transition_seasons_matrix)){
    # first, count how many observed transitions for that fish
    ntransitions <- length(transitions[i,][which(!(transitions[i,] %in% c(0,43)))])
    # then, take that many of the transitions seasons vector to populate
    transition_seasons_matrix[i,1:ntransitions] <- transition_seasons[(movements_counter+1):(movements_counter+ntransitions)]
    
    # increase the counter
    movements_counter <- movements_counter + ntransitions
  }
  
  state_history <- as.vector(t(transitions))
  transition_date_history <- as.vector(t(transition_dates))
  transition_seasons_history <- as.vector(t(transition_seasons_matrix))
  
  # drop the non-data
  nostate <- which(state_history == 0)
  state_history <- state_history[-nostate]
  transition_date_history <- transition_date_history[-nostate]
  transition_seasons_history <- transition_seasons_history[-nostate]
  
  transitions_df <- data.frame(state = rep(NA, length(state_history)-1),
                               previous_state = rep(NA, length(state_history)-1),
                               date = rep(NA, length(state_history)-1),
                               season = rep(NA, length(state_history)-1))
  
  for (i in 1:nrow(transitions_df)){
    if (i == 1){ # first observation has to come from the ether
      transitions_df$state[i] <- state_history[i]
      transitions_df$previous_state[i] <- 0
      transitions_df$date[i] <- transition_date_history[i]
      transitions_df$season[i] <- transition_seasons_history[i]
    } else {
      transitions_df$state[i] <- state_history[i]
      transitions_df$previous_state[i] <- ifelse(state_history[i-1]==43, 0, state_history[i-1]) # all first detections come from state 0 (fish purgatory)
      transitions_df$date[i] <- transition_date_history[i]
      transitions_df$season[i] <- transition_seasons_history[i]
    }
  }
  
  # drop the non-data
  transitions_df %>% 
    filter(!(state %in% c(0,43))) -> transitions_df
  
  # determine if the fish is going upstream or downstream
  transitions_df %>% 
    mutate(direction = ifelse(state > previous_state, "upstream", "downstream")) -> transitions_df
  
  return(transitions_df)
}


UCW_states_dates <- get_states_dates_direction(envir = UCW_envir)
UCH_states_dates <- get_states_dates_direction(envir = UCH_envir)
MCW_states_dates <- get_states_dates_direction(envir = MCW_envir)
MCH_states_dates <- get_states_dates_direction(envir = MCH_envir)
SRW_states_dates <- get_states_dates_direction(envir = SRW_envir)
SRH_states_dates <- get_states_dates_direction(envir = SRH_envir)

# UCW_states_dates <- data.frame(state = as.vector(UCW_envir$data$y),
#                                date = as.vector(UCW_envir$data$transition_dates))
# UCH_states_dates <- data.frame(state = as.vector(UCH_envir$data$y),
#                                date = as.vector(UCH_envir$data$transition_dates))
# MCW_states_dates <- data.frame(state = as.vector(MCW_envir$data$y),
#                                date = as.vector(MCW_envir$data$transition_dates))
# MCH_states_dates <- data.frame(state = as.vector(MCH_envir$data$y),
#                                date = as.vector(MCH_envir$data$transition_dates))
# SRW_states_dates <- data.frame(state = as.vector(SRW_envir$data$y),
#                                date = as.vector(SRW_envir$data$transition_dates))
# SRH_states_dates <- data.frame(state = as.vector(SRH_envir$data$y),
#                                date = as.vector(SRH_envir$data$transition_dates))

# Get states/dates by origin as well
get_origin_states_dates <- function(envir, origin_select, rear){
  
  transitions <- envir$data$y
  transition_dates <- envir$data$transition_dates
  # convert the seasons vector to the same shape as the state transitions
  transition_seasons <- envir$data$transition_seasons_vector
  transition_seasons_matrix <- matrix(nrow = nrow(transitions), ncol = ncol(transitions))
  
  movements_counter <- 0
  for (i in 1:nrow(transition_seasons_matrix)){
    # first, count how many observed transitions for that fish
    ntransitions <- length(transitions[i,][which(!(transitions[i,] %in% c(0,43)))])
    # then, take that many of the transitions seasons vector to populate
    transition_seasons_matrix[i,1:ntransitions] <- transition_seasons[(movements_counter+1):(movements_counter+ntransitions)]
    
    # increase the counter
    movements_counter <- movements_counter + ntransitions
  }
  
  
  state_history <- as.vector(t(transitions))
  transition_date_history <- as.vector(t(transition_dates))
  transition_seasons_history <- as.vector(t(transition_seasons_matrix))
  
  # get origin info
  origin_vector <- vector(length = nrow(envir$data$cat_X_mat))
  for(i in 1:nrow(envir$data$cat_X_mat)){
    origin_vector[i] <- which(envir$data$cat_X_mat[i,]==1)
  }
  
  origin_vector_matched <- rep(origin_vector, each = ncol(envir$data$y))
  
  # Now, keep only the origin selected
  if(rear == "wild"){
    origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$wild
  } else {
    origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$hatchery
  }
  
  origin_to_keep <- which(origin_vector_matched == origin_numeric)
  
  state_history <- state_history[origin_to_keep]
  transition_date_history <- transition_date_history[origin_to_keep]
  transition_seasons_history <- transition_seasons_history[origin_to_keep]
  
  # drop the non-data
  nostate <- which(state_history == 0)
  state_history <- state_history[-nostate]
  transition_date_history <- transition_date_history[-nostate]
  transition_seasons_history <- transition_seasons_history[-nostate]
  
  transitions_df <- data.frame(state = rep(NA, length(state_history)-1),
                               previous_state = rep(NA, length(state_history)-1),
                               date = rep(NA, length(state_history)-1),
                               season = rep(NA, length(state_history)-1))
  
  for (i in 1:nrow(transitions_df)){
    if (i == 1){ # first observation has to come from the ether
      transitions_df$state[i] <- state_history[i]
      transitions_df$previous_state[i] <- 0
      transitions_df$date[i] <- transition_date_history[i]
      transitions_df$season[i] <- transition_seasons_history[i]
    } else {
      transitions_df$state[i] <- state_history[i]
      transitions_df$previous_state[i] <- ifelse(state_history[i-1]==43, 0, state_history[i-1]) # all first detections come from state 0 (fish purgatory)
      transitions_df$date[i] <- transition_date_history[i]
      transitions_df$season[i] <- transition_seasons_history[i]
    }
  }
  
  # drop the non-data
  transitions_df %>% 
    filter(!(state %in% c(0,43))) -> transitions_df
  
  # determine if the fish is going upstream or downstream
  transitions_df %>% 
    mutate(direction = ifelse(state > previous_state, "upstream", "downstream")) -> transitions_df
  
  return(transitions_df)
  
}

# create a list that maps origin numbers (params) to what they actually are
natal_origins <- gsub(" Mouth| Upstream", "", model_states)
natal_origins <- natal_origins[!(duplicated(natal_origins))]
natal_origins <- natal_origins[!(grepl("mainstem", natal_origins))]
natal_origins <- natal_origins[!(grepl("other tributaries", natal_origins))]
natal_origins <- natal_origins[!(natal_origins == "loss")]

# Use the parameter map to index the right effects
origin_param_map <- data.frame(
  natal_origin = natal_origins,
  hatchery = c(NA, NA, NA, NA, 1, NA, 2, # MC
               1,NA,2,3, # UC
               5,NA,1,4,2,3), # SR,
  wild = c(1,3,NA,2,4,6,5, # MC
           1,2,NA,3, # UC
           6,1,2,5,3,4)) # SR


# wen_wild_states_dates <- get_origin_states_dates(envir = UCW_envir, origin_select = "Wenatchee River", rear = "wild")


compare_final_fate_fixcov_rear_type_UC <- function(niter, nsim,
                                                   origin_select,
                                                   fix_state, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                                   fix_spillwindow_season0_value = NA, fix_spillwindow_season1_value = NA,
                                                   fix_winterspill_value = NA,
                                                   fix_run_year = NA){
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,3)
    
    # index to the right origin param to turn it on
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_fixcov_UCW(nsim = nsim,
                                                    start_state = 2, states_dates = get_origin_states_dates(envir = UCW_envir, origin_select = origin_select, rear = "wild"),
                                                    origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                    temp_data = UCW_envir$data$temperature_data, spillwindow_data = UCW_envir$data$spill_window_data, 
                                                    winterspill_data = UCW_envir$data$winter_spill_days_data,
                                                    fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                    fix_temp_season1_value = fix_temp_season1_value, 
                                                    fix_spillwindow_season0_value = fix_spillwindow_season0_value, fix_spillwindow_season1_value = fix_spillwindow_season1_value,
                                                    fix_winterspill_value = fix_winterspill_value,
                                                    fix_run_year = fix_run_year)
      
      
      
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    ff_wild_quantiles -> ff_rear_quantiles
    
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    
    
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    hatchery_origin_params <- rep(0,3)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_fixcov_UCH(nsim = nsim,
                                                        start_state = 2, states_dates = get_origin_states_dates(envir = UCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                        origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], origin3 = hatchery_origin_params[3], 
                                                        temp_data = UCH_envir$data$temperature_data, spillwindow_data = UCH_envir$data$spill_window_data, 
                                                        winterspill_data = UCH_envir$data$winter_spill_days_data,
                                                        fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                        fix_temp_season1_value = fix_temp_season1_value, 
                                                        fix_spillwindow_season0_value = fix_spillwindow_season0_value, fix_spillwindow_season1_value = fix_spillwindow_season1_value,
                                                        fix_winterspill_value = fix_winterspill_value,
                                                        fix_run_year = fix_run_year)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_hatchery_quantiles -> ff_rear_quantiles
  } 
  # else run both
  else {
    
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,3)
    hatchery_origin_params <- rep(0,3)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_fixcov_UCW(nsim = nsim,
                                                    start_state = 2, states_dates = get_origin_states_dates(envir = UCW_envir, origin_select = origin_select, rear = "wild"),
                                                    origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                    temp_data = UCW_envir$data$temperature_data, spillwindow_data = UCW_envir$data$spill_window_data, 
                                                    winterspill_data = UCW_envir$data$winter_spill_days_data,
                                                    fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                    fix_temp_season1_value = fix_temp_season1_value, 
                                                    fix_spillwindow_season0_value = fix_spillwindow_season0_value, fix_spillwindow_season1_value = fix_spillwindow_season1_value,
                                                    fix_winterspill_value = fix_winterspill_value,
                                                    fix_run_year = fix_run_year)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_fixcov_UCH(nsim = nsim,
                                                        start_state = 2, states_dates = get_origin_states_dates(envir = UCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                        origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], origin3 = hatchery_origin_params[3], 
                                                        temp_data = UCH_envir$data$temperature_data, spillwindow_data = UCH_envir$data$spill_window_data, 
                                                        winterspill_data = UCH_envir$data$winter_spill_days_data,
                                                        fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                        fix_temp_season1_value = fix_temp_season1_value, 
                                                        fix_spillwindow_season0_value = fix_spillwindow_season0_value, fix_spillwindow_season1_value = fix_spillwindow_season1_value,
                                                        fix_winterspill_value = fix_winterspill_value,
                                                        fix_run_year = fix_run_year)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_wild_quantiles %>% 
      bind_rows(ff_hatchery_quantiles) -> ff_rear_quantiles
  }
  
  
  return(ff_rear_quantiles)
  
}

compare_final_fate_fixcov_rear_type_MC <- function(niter, nsim,
                                                   origin_select,
                                                   fix_state, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                                   fix_spillwindow_season0_value = NA, fix_spillwindow_season1_value = NA,
                                                   fix_winterspill_value = NA,
                                                   fix_run_year = NA){
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    
    # index to the right origin param to turn it on
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_fixcov_MCW(nsim = nsim,
                                                    start_state = 2, states_dates = get_origin_states_dates(envir = MCW_envir, origin_select = origin_select, rear = "wild"),
                                                    origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                    origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                    temp_data = MCW_envir$data$temperature_data, spillwindow_data = MCW_envir$data$spill_window_data, 
                                                    winterspill_data = MCW_envir$data$winter_spill_days_data,
                                                    fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                    fix_temp_season1_value = fix_temp_season1_value, 
                                                    fix_spillwindow_season0_value = fix_spillwindow_season0_value, fix_spillwindow_season1_value = fix_spillwindow_season1_value,
                                                    fix_winterspill_value = fix_winterspill_value,
                                                    fix_run_year = fix_run_year)
      
      
      
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    ff_wild_quantiles -> ff_rear_quantiles
    
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    
    
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    hatchery_origin_params <- c(0,0)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_fixcov_MCH(nsim = nsim,
                                                        start_state = 2, states_dates = get_origin_states_dates(envir = MCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                        origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                        temp_data = MCH_envir$data$temperature_data, spillwindow_data = MCH_envir$data$spill_window_data, 
                                                        winterspill_data = MCH_envir$data$winter_spill_days_data,
                                                        fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                        fix_temp_season1_value = fix_temp_season1_value, 
                                                        fix_spillwindow_season0_value = fix_spillwindow_season0_value, fix_spillwindow_season1_value = fix_spillwindow_season1_value,
                                                        fix_winterspill_value = fix_winterspill_value,
                                                        fix_run_year = fix_run_year)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_hatchery_quantiles -> ff_rear_quantiles
  } 
  # else run both
  else {
    
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    hatchery_origin_params <- c(0,0)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_fixcov_MCW(nsim = nsim,
                                                    start_state = 2, states_dates = get_origin_states_dates(envir = MCW_envir, origin_select = origin_select, rear = "wild"),
                                                    origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                    origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                    temp_data = MCW_envir$data$temperature_data, spillwindow_data = MCW_envir$data$spill_window_data, 
                                                    winterspill_data = MCW_envir$data$winter_spill_days_data,
                                                    fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                    fix_temp_season1_value = fix_temp_season1_value, 
                                                    fix_spillwindow_season0_value = fix_spillwindow_season0_value, fix_spillwindow_season1_value = fix_spillwindow_season1_value,
                                                    fix_winterspill_value = fix_winterspill_value,
                                                    fix_run_year = fix_run_year)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_fixcov_MCH(nsim = nsim,
                                                        start_state = 2, states_dates = get_origin_states_dates(envir = MCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                        origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                        temp_data = MCH_envir$data$temperature_data, spillwindow_data = MCH_envir$data$spill_window_data, 
                                                        winterspill_data = MCH_envir$data$winter_spill_days_data,
                                                        fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                        fix_temp_season1_value = fix_temp_season1_value, 
                                                        fix_spillwindow_season0_value = fix_spillwindow_season0_value, fix_spillwindow_season1_value = fix_spillwindow_season1_value,
                                                        fix_winterspill_value = fix_winterspill_value,
                                                        fix_run_year = fix_run_year)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_wild_quantiles %>% 
      bind_rows(ff_hatchery_quantiles) -> ff_rear_quantiles
  }
  
  
  return(ff_rear_quantiles)
  
}

compare_final_fate_fixcov_rear_type_SR <- function(niter, nsim,
                                                   origin_select,
                                                   fix_state, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                                   fix_spillwindow_season0_value = NA, fix_spillwindow_season1_value = NA,
                                                   fix_winterspill_value = NA,
                                                   fix_run_year = NA){
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    
    # index to the right origin param to turn it on
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_fixcov_SRW(nsim = nsim,
                                                    start_state = 2, states_dates = get_origin_states_dates(envir = SRW_envir, origin_select = origin_select, rear = "wild"),
                                                    origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                    origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                    temp_data = SRW_envir$data$temperature_data, spillwindow_data = SRW_envir$data$spill_window_data, 
                                                    winterspill_data = SRW_envir$data$winter_spill_days_data,
                                                    fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                    fix_temp_season1_value = fix_temp_season1_value, 
                                                    fix_spillwindow_season0_value = fix_spillwindow_season0_value, fix_spillwindow_season1_value = fix_spillwindow_season1_value,
                                                    fix_winterspill_value = fix_winterspill_value,
                                                    fix_run_year = fix_run_year)
      
      
      
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    ff_wild_quantiles -> ff_rear_quantiles
    
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    
    
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    hatchery_origin_params <- rep(0,5)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_fixcov_SRH(nsim = nsim,
                                                        start_state = 2, states_dates = get_origin_states_dates(envir = SRH_envir, origin_select = origin_select, rear = "hatchery"),
                                                        origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], origin3 = hatchery_origin_params[3], 
                                                        origin4 = hatchery_origin_params[4], origin5 = hatchery_origin_params[5], 
                                                        temp_data = SRH_envir$data$temperature_data, spillwindow_data = SRH_envir$data$spill_window_data, 
                                                        winterspill_data = SRH_envir$data$winter_spill_days_data,
                                                        fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                        fix_temp_season1_value = fix_temp_season1_value, 
                                                        fix_spillwindow_season0_value = fix_spillwindow_season0_value, fix_spillwindow_season1_value = fix_spillwindow_season1_value,
                                                        fix_winterspill_value = fix_winterspill_value,
                                                        fix_run_year = fix_run_year)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_hatchery_quantiles -> ff_rear_quantiles
  } 
  # else run both
  else {
    
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    hatchery_origin_params <- rep(0,5)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_fixcov_SRW(nsim = nsim,
                                                    start_state = 2, states_dates = get_origin_states_dates(envir = SRW_envir, origin_select = origin_select, rear = "wild"),
                                                    origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                    origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                    temp_data = SRW_envir$data$temperature_data, spillwindow_data = SRW_envir$data$spill_window_data, 
                                                    winterspill_data = SRW_envir$data$winter_spill_days_data,
                                                    fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                    fix_temp_season1_value = fix_temp_season1_value, 
                                                    fix_spillwindow_season0_value = fix_spillwindow_season0_value, fix_spillwindow_season1_value = fix_spillwindow_season1_value,
                                                    fix_winterspill_value = fix_winterspill_value,
                                                    fix_run_year = fix_run_year)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_fixcov_SRH(nsim = nsim,
                                                        start_state = 2, states_dates = get_origin_states_dates(envir = SRH_envir, origin_select = origin_select, rear = "hatchery"),
                                                        origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], origin3 = hatchery_origin_params[3], 
                                                        origin4 = hatchery_origin_params[4], origin5 = hatchery_origin_params[5], 
                                                        temp_data = SRH_envir$data$temperature_data, spillwindow_data = SRH_envir$data$spill_window_data, 
                                                        winterspill_data = SRH_envir$data$winter_spill_days_data,
                                                        fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                        fix_temp_season1_value = fix_temp_season1_value, 
                                                        fix_spillwindow_season0_value = fix_spillwindow_season0_value, fix_spillwindow_season1_value = fix_spillwindow_season1_value,
                                                        fix_winterspill_value = fix_winterspill_value,
                                                        fix_run_year = fix_run_year)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_wild_quantiles %>% 
      bind_rows(ff_hatchery_quantiles) -> ff_rear_quantiles
  }
  
  
  return(ff_rear_quantiles)
  
}

compare_final_fate_fixcov_times_rear_type_UC <- function(niter, nsim,
                                                         origin_select,
                                                         fix_state, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                                         fix_spillwindow_season0_value = NA, fix_spillwindow_season1_value = NA,
                                                         fix_spillwindow_states = rep(F, 9),
                                                         fix_spillwindow_values = rep(NA, 9),
                                                         fix_spillwindow_start_days = rep(NA, 9),
                                                         fix_spillwindow_end_days = rep(NA, 9),
                                                         fix_winterspill_value = NA,
                                                         fix_run_year = NA){
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    
    # index to the right origin param to turn it on
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_fixcov_times_UCW(nsim = nsim,
                                                          start_state = 2, states_dates = get_origin_states_dates(envir = UCW_envir, origin_select = origin_select, rear = "wild"),
                                                          origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                          temp_data = UCW_envir$data$temperature_data, spillwindow_data = UCW_envir$data$spill_window_data, 
                                                          winterspill_data = UCW_envir$data$winter_spill_days_data,
                                                          fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                          fix_temp_season1_value = fix_temp_season1_value, 
                                                          fix_spillwindow_states = fix_spillwindow_states,
                                                          fix_spillwindow_values = fix_spillwindow_values,
                                                          fix_spillwindow_start_days = fix_spillwindow_start_days,
                                                          fix_spillwindow_end_days = fix_spillwindow_end_days,
                                                          fix_winterspill_value = fix_winterspill_value,
                                                          fix_run_year = fix_run_year)
      
      
      
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    ff_wild_quantiles -> ff_rear_quantiles
    
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    
    
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    hatchery_origin_params <- c(0,0)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_fixcov_times_UCH(nsim = nsim,
                                                              start_state = 2, states_dates = get_origin_states_dates(envir = UCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                              origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], origin3 = hatchery_origin_params[3], 
                                                              temp_data = UCH_envir$data$temperature_data, spillwindow_data = UCH_envir$data$spill_window_data, 
                                                              winterspill_data = UCH_envir$data$winter_spill_days_data,
                                                              fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                              fix_temp_season1_value = fix_temp_season1_value, 
                                                              fix_spillwindow_states = fix_spillwindow_states,
                                                              fix_spillwindow_values = fix_spillwindow_values,
                                                              fix_spillwindow_start_days = fix_spillwindow_start_days,
                                                              fix_spillwindow_end_days = fix_spillwindow_end_days,
                                                              fix_winterspill_value = fix_winterspill_value,
                                                              fix_run_year = fix_run_year)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_hatchery_quantiles -> ff_rear_quantiles
  } 
  # else run both
  else {
    
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    hatchery_origin_params <- c(0,0)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_fixcov_times_UCW(nsim = nsim,
                                                          start_state = 2, states_dates = get_origin_states_dates(envir = UCW_envir, origin_select = origin_select, rear = "wild"),
                                                          origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                          temp_data = UCW_envir$data$temperature_data, spillwindow_data = UCW_envir$data$spill_window_data, 
                                                          winterspill_data = UCW_envir$data$winter_spill_days_data,
                                                          fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                          fix_temp_season1_value = fix_temp_season1_value, 
                                                          fix_spillwindow_states = fix_spillwindow_states,
                                                          fix_spillwindow_values = fix_spillwindow_values,
                                                          fix_spillwindow_start_days = fix_spillwindow_start_days,
                                                          fix_spillwindow_end_days = fix_spillwindow_end_days,
                                                          fix_winterspill_value = fix_winterspill_value,
                                                          fix_run_year = fix_run_year)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_fixcov_times_UCH(nsim = nsim,
                                                              start_state = 2, states_dates = get_origin_states_dates(envir = UCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                              origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], origin3 = hatchery_origin_params[3], 
                                                              temp_data = UCH_envir$data$temperature_data, spillwindow_data = UCH_envir$data$spill_window_data, 
                                                              winterspill_data = UCH_envir$data$winter_spill_days_data,
                                                              fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                              fix_temp_season1_value = fix_temp_season1_value, 
                                                              fix_spillwindow_states = fix_spillwindow_states,
                                                              fix_spillwindow_values = fix_spillwindow_values,
                                                              fix_spillwindow_start_days = fix_spillwindow_start_days,
                                                              fix_spillwindow_end_days = fix_spillwindow_end_days,
                                                              fix_winterspill_value = fix_winterspill_value,
                                                              fix_run_year = fix_run_year)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_wild_quantiles %>% 
      bind_rows(ff_hatchery_quantiles) -> ff_rear_quantiles
  }
  
  
  return(ff_rear_quantiles)
  
}

compare_final_fate_fixcov_times_rear_type_MC <- function(niter, nsim,
                                                         origin_select,
                                                         fix_state, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                                         fix_spillwindow_season0_value = NA, fix_spillwindow_season1_value = NA,
                                                         fix_spillwindow_states = rep(F, 9),
                                                         fix_spillwindow_values = rep(NA, 9),
                                                         fix_spillwindow_start_days = rep(NA, 9),
                                                         fix_spillwindow_end_days = rep(NA, 9),
                                                         fix_winterspill_value = NA,
                                                         fix_run_year = NA){
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    
    # index to the right origin param to turn it on
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_fixcov_times_MCW(nsim = nsim,
                                                          start_state = 2, states_dates = get_origin_states_dates(envir = MCW_envir, origin_select = origin_select, rear = "wild"),
                                                          origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                          origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                          temp_data = MCW_envir$data$temperature_data, spillwindow_data = MCW_envir$data$spill_window_data, 
                                                          winterspill_data = MCW_envir$data$winter_spill_days_data,
                                                          fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                          fix_temp_season1_value = fix_temp_season1_value, 
                                                          fix_spillwindow_states = fix_spillwindow_states,
                                                          fix_spillwindow_values = fix_spillwindow_values,
                                                          fix_spillwindow_start_days = fix_spillwindow_start_days,
                                                          fix_spillwindow_end_days = fix_spillwindow_end_days,
                                                          fix_winterspill_value = fix_winterspill_value,
                                                          fix_run_year = fix_run_year)
      
      
      
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    ff_wild_quantiles -> ff_rear_quantiles
    
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    
    
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    hatchery_origin_params <- c(0,0)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_fixcov_times_MCH(nsim = nsim,
                                                              start_state = 2, states_dates = get_origin_states_dates(envir = MCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                              origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                              temp_data = MCH_envir$data$temperature_data, spillwindow_data = MCH_envir$data$spill_window_data, 
                                                              winterspill_data = MCH_envir$data$winter_spill_days_data,
                                                              fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                              fix_temp_season1_value = fix_temp_season1_value, 
                                                              fix_spillwindow_states = fix_spillwindow_states,
                                                              fix_spillwindow_values = fix_spillwindow_values,
                                                              fix_spillwindow_start_days = fix_spillwindow_start_days,
                                                              fix_spillwindow_end_days = fix_spillwindow_end_days,
                                                              fix_winterspill_value = fix_winterspill_value,
                                                              fix_run_year = fix_run_year)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_hatchery_quantiles -> ff_rear_quantiles
  } 
  # else run both
  else {
    
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    hatchery_origin_params <- c(0,0)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_fixcov_times_MCW(nsim = nsim,
                                                          start_state = 2, states_dates = get_origin_states_dates(envir = MCW_envir, origin_select = origin_select, rear = "wild"),
                                                          origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                          origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                          temp_data = MCW_envir$data$temperature_data, spillwindow_data = MCW_envir$data$spill_window_data, 
                                                          winterspill_data = MCW_envir$data$winter_spill_days_data,
                                                          fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                          fix_temp_season1_value = fix_temp_season1_value, 
                                                          fix_spillwindow_states = fix_spillwindow_states,
                                                          fix_spillwindow_values = fix_spillwindow_values,
                                                          fix_spillwindow_start_days = fix_spillwindow_start_days,
                                                          fix_spillwindow_end_days = fix_spillwindow_end_days,
                                                          fix_winterspill_value = fix_winterspill_value,
                                                          fix_run_year = fix_run_year)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_fixcov_times_MCH(nsim = nsim,
                                                              start_state = 2, states_dates = get_origin_states_dates(envir = MCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                              origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                              temp_data = MCH_envir$data$temperature_data, spillwindow_data = MCH_envir$data$spill_window_data, 
                                                              winterspill_data = MCH_envir$data$winter_spill_days_data,
                                                              fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                              fix_temp_season1_value = fix_temp_season1_value, 
                                                              fix_spillwindow_states = fix_spillwindow_states,
                                                              fix_spillwindow_values = fix_spillwindow_values,
                                                              fix_spillwindow_start_days = fix_spillwindow_start_days,
                                                              fix_spillwindow_end_days = fix_spillwindow_end_days,
                                                              fix_winterspill_value = fix_winterspill_value,
                                                              fix_run_year = fix_run_year)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_wild_quantiles %>% 
      bind_rows(ff_hatchery_quantiles) -> ff_rear_quantiles
  }
  
  
  return(ff_rear_quantiles)
  
}

compare_final_fate_fixcov_times_rear_type_SR <- function(niter, nsim,
                                                         origin_select,
                                                         fix_state, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                                         fix_spillwindow_season0_value = NA, fix_spillwindow_season1_value = NA,
                                                         fix_spillwindow_states = rep(F, 9),
                                                         fix_spillwindow_values = rep(NA, 9),
                                                         fix_spillwindow_start_days = rep(NA, 9),
                                                         fix_spillwindow_end_days = rep(NA, 9),
                                                         fix_winterspill_value = NA,
                                                         fix_run_year = NA){
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    
    # index to the right origin param to turn it on
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_fixcov_times_SRW(nsim = nsim,
                                                          start_state = 2, states_dates = get_origin_states_dates(envir = SRW_envir, origin_select = origin_select, rear = "wild"),
                                                          origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                          origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                          temp_data = SRW_envir$data$temperature_data, spillwindow_data = SRW_envir$data$spill_window_data, 
                                                          winterspill_data = SRW_envir$data$winter_spill_days_data,
                                                          fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                          fix_temp_season1_value = fix_temp_season1_value, 
                                                          fix_spillwindow_states = fix_spillwindow_states,
                                                          fix_spillwindow_values = fix_spillwindow_values,
                                                          fix_spillwindow_start_days = fix_spillwindow_start_days,
                                                          fix_spillwindow_end_days = fix_spillwindow_end_days,
                                                          fix_winterspill_value = fix_winterspill_value,
                                                          fix_run_year = fix_run_year)
      
      
      
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    ff_wild_quantiles -> ff_rear_quantiles
    
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    
    
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    hatchery_origin_params <- c(0,0)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_fixcov_times_SRH(nsim = nsim,
                                                              start_state = 2, states_dates = get_origin_states_dates(envir = SRH_envir, origin_select = origin_select, rear = "hatchery"),
                                                              origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], origin3 = hatchery_origin_params[3], 
                                                              origin4 = hatchery_origin_params[4], origin5 = hatchery_origin_params[5], 
                                                              temp_data = SRH_envir$data$temperature_data, spillwindow_data = SRH_envir$data$spill_window_data, 
                                                              winterspill_data = SRH_envir$data$winter_spill_days_data,
                                                              fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                              fix_temp_season1_value = fix_temp_season1_value, 
                                                              fix_spillwindow_states = fix_spillwindow_states,
                                                              fix_spillwindow_values = fix_spillwindow_values,
                                                              fix_spillwindow_start_days = fix_spillwindow_start_days,
                                                              fix_spillwindow_end_days = fix_spillwindow_end_days,
                                                              fix_winterspill_value = fix_winterspill_value,
                                                              fix_run_year = fix_run_year)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_hatchery_quantiles -> ff_rear_quantiles
  } 
  # else run both
  else {
    
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    hatchery_origin_params <- c(0,0)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_fixcov_times_SRW(nsim = nsim,
                                                          start_state = 2, states_dates = get_origin_states_dates(envir = SRW_envir, origin_select = origin_select, rear = "wild"),
                                                          origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                          origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                          temp_data = SRW_envir$data$temperature_data, spillwindow_data = SRW_envir$data$spill_window_data, 
                                                          winterspill_data = SRW_envir$data$winter_spill_days_data,
                                                          fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                          fix_temp_season1_value = fix_temp_season1_value, 
                                                          fix_spillwindow_states = fix_spillwindow_states,
                                                          fix_spillwindow_values = fix_spillwindow_values,
                                                          fix_spillwindow_start_days = fix_spillwindow_start_days,
                                                          fix_spillwindow_end_days = fix_spillwindow_end_days,
                                                          fix_winterspill_value = fix_winterspill_value,
                                                          fix_run_year = fix_run_year)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_fixcov_times_SRH(nsim = nsim,
                                                              start_state = 2, states_dates = get_origin_states_dates(envir = SRH_envir, origin_select = origin_select, rear = "hatchery"),
                                                              origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], origin3 = hatchery_origin_params[3], 
                                                              origin4 = hatchery_origin_params[4], origin5 = hatchery_origin_params[5], 
                                                              temp_data = SRH_envir$data$temperature_data, spillwindow_data = SRH_envir$data$spill_window_data, 
                                                              winterspill_data = SRH_envir$data$winter_spill_days_data,
                                                              fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                              fix_temp_season1_value = fix_temp_season1_value, 
                                                              fix_spillwindow_states = fix_spillwindow_states,
                                                              fix_spillwindow_values = fix_spillwindow_values,
                                                              fix_spillwindow_start_days = fix_spillwindow_start_days,
                                                              fix_spillwindow_end_days = fix_spillwindow_end_days,
                                                              fix_winterspill_value = fix_winterspill_value,
                                                              fix_run_year = fix_run_year)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_wild_quantiles %>% 
      bind_rows(ff_hatchery_quantiles) -> ff_rear_quantiles
  }
  
  
  return(ff_rear_quantiles)
  
}




#### Set spill window only for specific times ####

# this new function allows you to toggle conditions for each of the mainstem states (1:9),
# with conditions set in each state for specific times of the year

# The dates in the fix_spillwindow_start_days and fix_spillwindow_end_days arguments
# are formatted as julian day; if you want to set it to a specific day of year,
# just use yday("2024-06-01") for example, to extract the julian day for June 1
final_fates_simulation_fixcov_times_UCW <- function(nsim,
                                                    start_state = 2, states_dates,
                                                    origin1 = 0, origin2 = 0, origin3 = 0,
                                                    temp_data, spillwindow_data, winterspill_data,
                                                    fix_state = NA, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                                    fix_spillwindow_states = rep(F, 9),
                                                    fix_spillwindow_values = rep(NA, 9),
                                                    fix_spillwindow_start_days = rep(NA, 9),
                                                    fix_spillwindow_end_days = rep(NA, 9),
                                                    fix_winterspill_value = NA,
                                                    fix_run_year = NA){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(UCW_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state across all years
    winterspill[i] <- median(winterspill_data[,i])
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(UCW_states_dates, state == i  & direction == "upstream")$season)/length(subset(UCW_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(UCW_states_dates, state == i  & direction == "downstream")$season)/length(subset(UCW_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(UCW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    # get the date numeric for each day/month of that run year
    run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
    run_year_start <- seq(ymd("2004-06-01"), ymd("2022-06-01"), by = "years")
    run_year_end <- seq(ymd("2005-05-31"), ymd("2023-05-31"), by = "years")
    run_year_numeric = seq(4, 22, 1)
    
    run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)
    
    
    
    # rowwise() %>%
    # dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_actual & run_year_end >= date_actual)$run_year) %>% 
    # filter(run_year == fix_run_year)-> UCW_states_dates_run_year
    
    # subset states_dates for origin and population to pull out just month day, and then
    # use those months and days to index temps from that year
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    dates_actual %>% 
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) -> dates_actual
    
    
    subset(run_year_df, run_year == fix_run_year) -> fix_run_year_df
    subset(dates_actual, date_actual <= fix_run_year_df$run_year_end & date_actual >= fix_run_year_df$run_year_start) %>% 
      dplyr::rename(date_run_year = date) -> dates_fix_run_year
    
    # Deal with leap years by using 2/28 twice
    if(!("02-29" %in% dates_fix_run_year$month_day)){
      subset(dates_fix_run_year, month_day == "02-28") %>% 
        mutate(month_day = "02-29") -> leap_year
      
      dates_fix_run_year %>% 
        bind_rows(., leap_year) -> dates_fix_run_year
      
    }
    
    
    
    
    
    UCW_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> UCW_states_month_day
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> states_month_day
    
    
    temp_upstream_season0 <- rep(0, length(model_states))
    temp_upstream_season1 <- rep(0, length(model_states))
    temp_downstream_season0 <- rep(0, length(model_states))
    temp_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(UCW_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(UCW_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(UCW_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(UCW_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(UCW_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(UCW_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(UCW_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(UCW_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    spillwindow_upstream_season0 <- rep(0, length(model_states))
    spillwindow_upstream_season1 <- rep(0, length(model_states))
    spillwindow_downstream_season0 <- rep(0, length(model_states))
    spillwindow_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(UCW_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(UCW_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(UCW_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(UCW_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(UCW_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(UCW_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(UCW_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(UCW_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  
  # Fix covariate values to the ones that you are setting
  if(!(is.na(fix_temp_season0_value))){
    temp_upstream_season0[fix_state] <- fix_temp_season0_value
    temp_downstream_season0[fix_state] <- fix_temp_season0_value
  }
  
  if(!(is.na(fix_temp_season1_value))){
    temp_upstream_season1[fix_state] <- fix_temp_season1_value
    temp_downstream_season1[fix_state] <- fix_temp_season1_value
  }
  
  if(!(is.na(fix_winterspill_value))){
    winterspill[fix_state] <- fix_winterspill_value
  }
  
  # Fix covariate values for specific time periods
  
  if(any(!(is.na(fix_spillwindow_states))) & any(!(is.na(fix_spillwindow_values))) & 
     any(!(is.na(fix_spillwindow_start_days))) & any(!(is.na(fix_spillwindow_end_days)))){
    
    # calculate additional values:
    # get the median value for season0 and season1, excluding the time period in which
    # you are fixing spill
    # set season0_fix and season1_fix to the value that you are fixing it to
    # then, use the simulation to select between season0_fix and season0/season1_fix and season1,
    # depending on the relative frequency of movement in the different time periods
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    
    # now loop through the mainstem states, and determine if they dates fall within
    # the fix range for each state
    to_fix <- data.frame(state = 1:9, state_to_fix = fix_spillwindow_states,
                         fix_start_date = fix_spillwindow_start_days,
                         fix_end_date = fix_spillwindow_end_days)
    
    UCW_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(fix = ifelse(yday(date_actual) > fix_start_date &
                            yday(date_actual) < fix_end_date, 1, 0)) %>% 
      mutate(fix = ifelse(is.na(fix), 0, fix)) -> UCW_states_dates_fix
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(fix = ifelse(yday(date_actual) > fix_start_date &
                            yday(date_actual) < fix_end_date, 1, 0)) %>% 
      mutate(fix = ifelse(is.na(fix), 0, fix)) -> states_dates_fix
    
    # within season upstream vs. season downstream, calculate proportion of time to fix
    season0_upstream_fixprop <- rep(0, length(model_states))
    season0_downstream_fixprop <- rep(0, length(model_states))
    season1_upstream_fixprop <- rep(0, length(model_states))
    season1_downstream_fixprop <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # proportion of season0, upstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(UCW_states_dates_fix, state == i & direction == "upstream" & season == 0)) == 0){
          season0_upstream_fixprop[i] <- 0
        } else {
          season0_upstream_fixprop[i] <- sum(subset(UCW_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
            length(subset(UCW_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
        }
      } else {
        season0_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
      }
      # proportion of season1, upstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(UCW_states_dates_fix, state == i & direction == "upstream" & season == 1)) == 0){
          season1_upstream_fixprop[i] <- 0
        } else {
          season1_upstream_fixprop[i] <- sum(subset(UCW_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
            length(subset(UCW_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)
        }
      } else {
        season1_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1 )$fix)
      }
      # proportion of season0, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(UCW_states_dates_fix, state == i & direction == "downstream" & season == 0)) == 0){
          season0_downstream_fixprop[i] <- 0
        } else {
          season0_downstream_fixprop[i] <- sum(subset(UCW_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
            length(subset(UCW_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
        }
      } else {
        season0_downstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
      }
      # proportion of season1, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(UCW_states_dates_fix, state == i & direction == "downstream" & season == 1)) == 0){
          season1_downstream_fixprop[i] <- 0
        } else {
          season1_downstream_fixprop[i] <- sum(subset(UCW_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)/
            length(subset(UCW_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)
        }
      } else {
        season1_downstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1 )$fix)
      }
      
    }
    
    # Now, get the conditions outside of the fixed windows, for any state that you're fixing
    spillwindow_upstream_season0_notfixed <- rep(0, length(model_states))
    spillwindow_upstream_season1_notfixed <- rep(0, length(model_states))
    spillwindow_downstream_season0_notfixed <- rep(0, length(model_states))
    spillwindow_downstream_season1_notfixed <- rep(0, length(model_states))
    
    # Remove any states_month_day that are within those windows
    states_month_day %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(drop_date = ifelse(yday(date_actual.x) > fix_start_date &
                                  yday(date_actual.x) < fix_end_date, 1, 0)) %>% 
      mutate(drop_date = ifelse(is.na(drop_date), 0, drop_date)) %>% 
      filter(drop_date == 0) -> states_month_day_trimmed
    
    UCW_states_month_day %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(drop_date = ifelse(yday(date_actual.x) > fix_start_date &
                                  yday(date_actual.x) < fix_end_date, 1, 0)) %>% 
      mutate(drop_date = ifelse(is.na(drop_date), 0, drop_date)) %>% 
      filter(drop_date == 0) -> UCW_states_month_day_trimmed
    
    
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(UCW_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(UCW_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(UCW_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(UCW_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(UCW_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(UCW_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(UCW_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season1_notfixed[i] <- median(spillwindow_data[subset(UCW_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  # toggle what the spill value should be based on 1) how the function is run and 
  # 2) if spill is fixed, fixed spill vs. not fixed spill based on resampling fixprop
  
  # We are going to change the already created spillwindow_upstream_season1 etc.
  # vectors that contain the spillwindow values for
  # upstream/downstream, season0/season1, and each state
  
  
  # first - toggle based on whether or not we're choosing to fix spill here
  if(any(!(is.na(fix_spillwindow_states))) & any(!(is.na(fix_spillwindow_values))) & 
     any(!(is.na(fix_spillwindow_start_days))) & any(!(is.na(fix_spillwindow_end_days)))){
    
    # loop through this vector to determine if we're changing spill values
    for (k in 1:length(fix_spillwindow_states)){
      if (fix_spillwindow_states[k] == T){
        # resample to determine if we're going to use the fixed spill values
        fixspill_upstream_season0_sim <- sample(c(0,1), 1, prob = c(1-season0_upstream_fixprop[k], season0_upstream_fixprop[k]))
        fixspill_upstream_season1_sim <- sample(c(0,1), 1, prob = c(1-season1_upstream_fixprop[k], season1_upstream_fixprop[k]))
        fixspill_downstream_season0_sim <- sample(c(0,1), 1, prob = c(1-season0_downstream_fixprop[k], season0_downstream_fixprop[k]))
        fixspill_downstream_season1_sim <- sample(c(0,1), 1, prob = c(1-season1_downstream_fixprop[k], season1_downstream_fixprop[k]))
        
        # if we resampled that we are going to change to the fixed spill values, change them
        # if we resampled that we are not going to change to the fixed spill values because
        # we're not in the right time of year, then change the spill values to those
        # for the remainder of the year
        if(fixspill_upstream_season0_sim == 1){
          spillwindow_upstream_season0[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_upstream_season0[k] <- spillwindow_upstream_season0_notfixed[k]
        }
        if(fixspill_upstream_season1_sim == 1){
          spillwindow_upstream_season1[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_upstream_season1[k] <- spillwindow_upstream_season1_notfixed[k]
        }
        if(fixspill_downstream_season0_sim == 1){
          spillwindow_downstream_season0[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_downstream_season0[k] <- spillwindow_downstream_season0_notfixed[k]
        }
        if(fixspill_downstream_season1_sim == 1){
          spillwindow_downstream_season1[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_downstream_season1[k] <- spillwindow_downstream_season1_notfixed[k]
        }
        
      } else {
        # do nothing!
      }
      
    }
    
    
    # If we're not choosing to fix spill, we don't need to do anything - just leave spillwindow vectors as is  
  } else {
    # do nothing!
  }
  
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      
      if (upstream_season_sim == 1){
        
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCW[i,possible_movements[j],iter] +
                                                                    btemp1_array_UCW[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_UCW[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_UCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_UCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_UCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_UCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    borigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_UCW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_UCW[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCW[i,possible_movements,iter] +
                    btemp1_array_UCW[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_UCW[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_UCW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_UCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_UCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_UCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    borigin1_array_UCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCW[i,possible_movements,iter]*origin3))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCW[i,possible_movements[j],iter] +
                                                                    btemp0_array_UCW[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_UCW[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_UCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_UCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_UCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_UCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    borigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_UCW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_UCW[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCW[i,possible_movements,iter] +
                    btemp0_array_UCW[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_UCW[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_UCW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_UCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_UCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_UCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    borigin1_array_UCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCW[i,possible_movements,iter]*origin3))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCW[i,possible_movements[j],iter] +
                                                                      btemp1_array_UCW[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_UCW[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_UCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_UCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_UCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_UCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      borigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_UCW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_UCW[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCW[i,possible_movements,iter] +
                    btemp1_array_UCW[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_UCW[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_UCW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_UCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_UCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_UCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    borigin1_array_UCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCW[i,possible_movements,iter]*origin3))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCW[i,possible_movements[j],iter] +
                                                                      btemp0_array_UCW[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_UCW[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_UCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_UCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_UCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_UCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      borigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_UCW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_UCW[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCW[i,possible_movements,iter] +
                    btemp0_array_UCW[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_UCW[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_UCW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_UCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_UCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_UCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    borigin1_array_UCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCW[i,possible_movements,iter]*origin3))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_fixcov_times_UCH <- function(nsim,
                                                    start_state = 2, states_dates,
                                                    origin1 = 0, origin2 = 0, origin3 = 0,
                                                    temp_data, spillwindow_data, winterspill_data,
                                                    fix_state = NA, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                                    fix_spillwindow_states = rep(F, 9),
                                                    fix_spillwindow_values = rep(NA, 9),
                                                    fix_spillwindow_start_days = rep(NA, 9),
                                                    fix_spillwindow_end_days = rep(NA, 9),
                                                    fix_winterspill_value = NA,
                                                    fix_run_year = NA){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(UCH_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state across all years
    winterspill[i] <- median(winterspill_data[,i])
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(UCH_states_dates, state == i  & direction == "upstream")$season)/length(subset(UCH_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(UCH_states_dates, state == i  & direction == "downstream")$season)/length(subset(UCH_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(UCH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    # get the date numeric for each day/month of that run year
    run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
    run_year_start <- seq(ymd("2004-06-01"), ymd("2022-06-01"), by = "years")
    run_year_end <- seq(ymd("2005-05-31"), ymd("2023-05-31"), by = "years")
    run_year_numeric = seq(4, 22, 1)
    
    run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)
    
    
    
    # rowwise() %>%
    # dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_actual & run_year_end >= date_actual)$run_year) %>% 
    # filter(run_year == fix_run_year)-> UCH_states_dates_run_year
    
    # subset states_dates for origin and population to pull out just month day, and then
    # use those months and days to index temps from that year
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    dates_actual %>% 
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) -> dates_actual
    
    
    subset(run_year_df, run_year == fix_run_year) -> fix_run_year_df
    subset(dates_actual, date_actual <= fix_run_year_df$run_year_end & date_actual >= fix_run_year_df$run_year_start) %>% 
      dplyr::rename(date_run_year = date) -> dates_fix_run_year
    
    # Deal with leap years by using 2/28 twice
    if(!("02-29" %in% dates_fix_run_year$month_day)){
      subset(dates_fix_run_year, month_day == "02-28") %>% 
        mutate(month_day = "02-29") -> leap_year
      
      dates_fix_run_year %>% 
        bind_rows(., leap_year) -> dates_fix_run_year
      
    }
    
    
    
    
    
    UCH_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> UCH_states_month_day
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> states_month_day
    
    
    temp_upstream_season0 <- rep(0, length(model_states))
    temp_upstream_season1 <- rep(0, length(model_states))
    temp_downstream_season0 <- rep(0, length(model_states))
    temp_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(UCH_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(UCH_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(UCH_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(UCH_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(UCH_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(UCH_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(UCH_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(UCH_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    spillwindow_upstream_season0 <- rep(0, length(model_states))
    spillwindow_upstream_season1 <- rep(0, length(model_states))
    spillwindow_downstream_season0 <- rep(0, length(model_states))
    spillwindow_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(UCH_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(UCH_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(UCH_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(UCH_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(UCH_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(UCH_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(UCH_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(UCH_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  
  # Fix covariate values to the ones that you are setting
  if(!(is.na(fix_temp_season0_value))){
    temp_upstream_season0[fix_state] <- fix_temp_season0_value
    temp_downstream_season0[fix_state] <- fix_temp_season0_value
  }
  
  if(!(is.na(fix_temp_season1_value))){
    temp_upstream_season1[fix_state] <- fix_temp_season1_value
    temp_downstream_season1[fix_state] <- fix_temp_season1_value
  }
  
  if(!(is.na(fix_winterspill_value))){
    winterspill[fix_state] <- fix_winterspill_value
  }
  
  # Fix covariate values for specific time periods
  
  if(any(!(is.na(fix_spillwindow_states))) & any(!(is.na(fix_spillwindow_values))) & 
     any(!(is.na(fix_spillwindow_start_days))) & any(!(is.na(fix_spillwindow_end_days)))){
    
    # calculate additional values:
    # get the median value for season0 and season1, excluding the time period in which
    # you are fixing spill
    # set season0_fix and season1_fix to the value that you are fixing it to
    # then, use the simulation to select between season0_fix and season0/season1_fix and season1,
    # depending on the relative frequency of movement in the different time periods
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    
    # now loop through the mainstem states, and determine if they dates fall within
    # the fix range for each state
    to_fix <- data.frame(state = 1:9, state_to_fix = fix_spillwindow_states,
                         fix_start_date = fix_spillwindow_start_days,
                         fix_end_date = fix_spillwindow_end_days)
    
    UCH_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(fix = ifelse(yday(date_actual) > fix_start_date &
                            yday(date_actual) < fix_end_date, 1, 0)) %>% 
      mutate(fix = ifelse(is.na(fix), 0, fix)) -> UCH_states_dates_fix
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(fix = ifelse(yday(date_actual) > fix_start_date &
                            yday(date_actual) < fix_end_date, 1, 0)) %>% 
      mutate(fix = ifelse(is.na(fix), 0, fix)) -> states_dates_fix
    
    # within season upstream vs. season downstream, calculate proportion of time to fix
    season0_upstream_fixprop <- rep(0, length(model_states))
    season0_downstream_fixprop <- rep(0, length(model_states))
    season1_upstream_fixprop <- rep(0, length(model_states))
    season1_downstream_fixprop <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # proportion of season0, upstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(UCH_states_dates_fix, state == i & direction == "upstream" & season == 0)) == 0){
          season0_upstream_fixprop[i] <- 0
        } else {
          season0_upstream_fixprop[i] <- sum(subset(UCH_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
            length(subset(UCH_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
        }
      } else {
        season0_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
      }
      # proportion of season1, upstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(UCH_states_dates_fix, state == i & direction == "upstream" & season == 1)) == 0){
          season1_upstream_fixprop[i] <- 0
        } else {
          season1_upstream_fixprop[i] <- sum(subset(UCH_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
            length(subset(UCH_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)
        }
      } else {
        season1_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1 )$fix)
      }
      # proportion of season0, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(UCH_states_dates_fix, state == i & direction == "downstream" & season == 0)) == 0){
          season0_downstream_fixprop[i] <- 0
        } else {
          season0_downstream_fixprop[i] <- sum(subset(UCH_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
            length(subset(UCH_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
        }
      } else {
        season0_downstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
      }
      # proportion of season1, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(UCH_states_dates_fix, state == i & direction == "downstream" & season == 1)) == 0){
          season1_downstream_fixprop[i] <- 0
        } else {
          season1_downstream_fixprop[i] <- sum(subset(UCH_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)/
            length(subset(UCH_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)
        }
      } else {
        season1_downstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1 )$fix)
      }
      
    }
    
    # Now, get the conditions outside of the fixed windows, for any state that you're fixing
    spillwindow_upstream_season0_notfixed <- rep(0, length(model_states))
    spillwindow_upstream_season1_notfixed <- rep(0, length(model_states))
    spillwindow_downstream_season0_notfixed <- rep(0, length(model_states))
    spillwindow_downstream_season1_notfixed <- rep(0, length(model_states))
    
    # Remove any states_month_day that are within those windows
    states_month_day %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(drop_date = ifelse(yday(date_actual.x) > fix_start_date &
                                  yday(date_actual.x) < fix_end_date, 1, 0)) %>% 
      mutate(drop_date = ifelse(is.na(drop_date), 0, drop_date)) %>% 
      filter(drop_date == 0) -> states_month_day_trimmed
    
    UCH_states_month_day %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(drop_date = ifelse(yday(date_actual.x) > fix_start_date &
                                  yday(date_actual.x) < fix_end_date, 1, 0)) %>% 
      mutate(drop_date = ifelse(is.na(drop_date), 0, drop_date)) %>% 
      filter(drop_date == 0) -> UCH_states_month_day_trimmed
    
    
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(UCH_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(UCH_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(UCH_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(UCH_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(UCH_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(UCH_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(UCH_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season1_notfixed[i] <- median(spillwindow_data[subset(UCH_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  # toggle what the spill value should be based on 1) how the function is run and 
  # 2) if spill is fixed, fixed spill vs. not fixed spill based on resampling fixprop
  
  # We are going to change the already created spillwindow_upstream_season1 etc.
  # vectors that contain the spillwindow values for
  # upstream/downstream, season0/season1, and each state
  
  
  # first - toggle based on whether or not we're choosing to fix spill here
  if(any(!(is.na(fix_spillwindow_states))) & any(!(is.na(fix_spillwindow_values))) & 
     any(!(is.na(fix_spillwindow_start_days))) & any(!(is.na(fix_spillwindow_end_days)))){
    
    # loop through this vector to determine if we're changing spill values
    for (k in 1:length(fix_spillwindow_states)){
      if (fix_spillwindow_states[k] == T){
        # resample to determine if we're going to use the fixed spill values
        fixspill_upstream_season0_sim <- sample(c(0,1), 1, prob = c(1-season0_upstream_fixprop[k], season0_upstream_fixprop[k]))
        fixspill_upstream_season1_sim <- sample(c(0,1), 1, prob = c(1-season1_upstream_fixprop[k], season1_upstream_fixprop[k]))
        fixspill_downstream_season0_sim <- sample(c(0,1), 1, prob = c(1-season0_downstream_fixprop[k], season0_downstream_fixprop[k]))
        fixspill_downstream_season1_sim <- sample(c(0,1), 1, prob = c(1-season1_downstream_fixprop[k], season1_downstream_fixprop[k]))
        
        # if we resampled that we are going to change to the fixed spill values, change them
        # if we resampled that we are not going to change to the fixed spill values because
        # we're not in the right time of year, then change the spill values to those
        # for the remainder of the year
        if(fixspill_upstream_season0_sim == 1){
          spillwindow_upstream_season0[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_upstream_season0[k] <- spillwindow_upstream_season0_notfixed[k]
        }
        if(fixspill_upstream_season1_sim == 1){
          spillwindow_upstream_season1[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_upstream_season1[k] <- spillwindow_upstream_season1_notfixed[k]
        }
        if(fixspill_downstream_season0_sim == 1){
          spillwindow_downstream_season0[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_downstream_season0[k] <- spillwindow_downstream_season0_notfixed[k]
        }
        if(fixspill_downstream_season1_sim == 1){
          spillwindow_downstream_season1[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_downstream_season1[k] <- spillwindow_downstream_season1_notfixed[k]
        }
        
      } else {
        # do nothing!
      }
      
    }
    
    
    # If we're not choosing to fix spill, we don't need to do anything - just leave spillwindow vectors as is  
  } else {
    # do nothing!
  }
  
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      
      if (upstream_season_sim == 1){
        
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCH[i,possible_movements[j],iter] +
                                                                    btemp1_array_UCH[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_UCH[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_UCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_UCH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_UCH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_UCH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    borigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_UCH[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_UCH[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCH[i,possible_movements,iter] +
                    btemp1_array_UCH[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_UCH[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_UCH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_UCH[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_UCH[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_UCH[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    borigin1_array_UCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCH[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCH[i,possible_movements,iter]*origin3))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCH[i,possible_movements[j],iter] +
                                                                    btemp0_array_UCH[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_UCH[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_UCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_UCH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_UCH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_UCH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    borigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_UCH[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_UCH[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCH[i,possible_movements,iter] +
                    btemp0_array_UCH[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_UCH[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_UCH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_UCH[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_UCH[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_UCH[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    borigin1_array_UCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCH[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCH[i,possible_movements,iter]*origin3))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCH[i,possible_movements[j],iter] +
                                                                      btemp1_array_UCH[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_UCH[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_UCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_UCH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_UCH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_UCH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      borigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_UCH[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_UCH[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCH[i,possible_movements,iter] +
                    btemp1_array_UCH[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_UCH[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_UCH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_UCH[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_UCH[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_UCH[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    borigin1_array_UCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCH[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCH[i,possible_movements,iter]*origin3))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCH[i,possible_movements[j],iter] +
                                                                      btemp0_array_UCH[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_UCH[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_UCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_UCH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_UCH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_UCH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      borigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_UCH[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_UCH[i,possible_movements[j],iter]*origin3)/
          sum(exp(b0_array_UCH[i,possible_movements,iter] +
                    btemp0_array_UCH[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_UCH[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_UCH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_UCH[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_UCH[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_UCH[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    borigin1_array_UCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_UCH[i,possible_movements,iter]*origin2 +
                    borigin3_array_UCH[i,possible_movements,iter]*origin3))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_fixcov_times_MCW <- function(nsim,
                                                    start_state = 2, states_dates,
                                                    origin1 = 0, origin2 = 0, origin3 = 0,
                                                    origin4 = 0, origin5 = 0, origin6 = 0,
                                                    temp_data, spillwindow_data, winterspill_data,
                                                    fix_state = NA, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                                    fix_spillwindow_states = rep(F, 9),
                                                    fix_spillwindow_values = rep(NA, 9),
                                                    fix_spillwindow_start_days = rep(NA, 9),
                                                    fix_spillwindow_end_days = rep(NA, 9),
                                                    fix_winterspill_value = NA,
                                                    fix_run_year = NA){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(MCW_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state across all years
    winterspill[i] <- median(winterspill_data[,i])
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(MCW_states_dates, state == i  & direction == "upstream")$season)/length(subset(MCW_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(MCW_states_dates, state == i  & direction == "downstream")$season)/length(subset(MCW_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    # get the date numeric for each day/month of that run year
    run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
    run_year_start <- seq(ymd("2004-06-01"), ymd("2022-06-01"), by = "years")
    run_year_end <- seq(ymd("2005-05-31"), ymd("2023-05-31"), by = "years")
    run_year_numeric = seq(4, 22, 1)
    
    run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)
    
    
    
    # rowwise() %>%
    # dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_actual & run_year_end >= date_actual)$run_year) %>% 
    # filter(run_year == fix_run_year)-> MCW_states_dates_run_year
    
    # subset states_dates for origin and population to pull out just month day, and then
    # use those months and days to index temps from that year
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    dates_actual %>% 
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) -> dates_actual
    
    
    subset(run_year_df, run_year == fix_run_year) -> fix_run_year_df
    subset(dates_actual, date_actual <= fix_run_year_df$run_year_end & date_actual >= fix_run_year_df$run_year_start) %>% 
      dplyr::rename(date_run_year = date) -> dates_fix_run_year
    
    # Deal with leap years by using 2/28 twice
    if(!("02-29" %in% dates_fix_run_year$month_day)){
      subset(dates_fix_run_year, month_day == "02-28") %>% 
        mutate(month_day = "02-29") -> leap_year
      
      dates_fix_run_year %>% 
        bind_rows(., leap_year) -> dates_fix_run_year
      
    }
    
    
    
    
    
    MCW_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> MCW_states_month_day
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> states_month_day
    
    
    temp_upstream_season0 <- rep(0, length(model_states))
    temp_upstream_season1 <- rep(0, length(model_states))
    temp_downstream_season0 <- rep(0, length(model_states))
    temp_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(MCW_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(MCW_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(MCW_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(MCW_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    spillwindow_upstream_season0 <- rep(0, length(model_states))
    spillwindow_upstream_season1 <- rep(0, length(model_states))
    spillwindow_downstream_season0 <- rep(0, length(model_states))
    spillwindow_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(MCW_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(MCW_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(MCW_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(MCW_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  
  # Fix covariate values to the ones that you are setting
  if(!(is.na(fix_temp_season0_value))){
    temp_upstream_season0[fix_state] <- fix_temp_season0_value
    temp_downstream_season0[fix_state] <- fix_temp_season0_value
  }
  
  if(!(is.na(fix_temp_season1_value))){
    temp_upstream_season1[fix_state] <- fix_temp_season1_value
    temp_downstream_season1[fix_state] <- fix_temp_season1_value
  }
  
  if(!(is.na(fix_winterspill_value))){
    winterspill[fix_state] <- fix_winterspill_value
  }
  
  # Fix covariate values for specific time periods
  
  if(any(!(is.na(fix_spillwindow_states))) & any(!(is.na(fix_spillwindow_values))) & 
     any(!(is.na(fix_spillwindow_start_days))) & any(!(is.na(fix_spillwindow_end_days)))){
    
    # calculate additional values:
    # get the median value for season0 and season1, excluding the time period in which
    # you are fixing spill
    # set season0_fix and season1_fix to the value that you are fixing it to
    # then, use the simulation to select between season0_fix and season0/season1_fix and season1,
    # depending on the relative frequency of movement in the different time periods
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    
    # now loop through the mainstem states, and determine if they dates fall within
    # the fix range for each state
    to_fix <- data.frame(state = 1:9, state_to_fix = fix_spillwindow_states,
                         fix_start_date = fix_spillwindow_start_days,
                         fix_end_date = fix_spillwindow_end_days)
    
    MCW_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(fix = ifelse(yday(date_actual) > fix_start_date &
                            yday(date_actual) < fix_end_date, 1, 0)) %>% 
      mutate(fix = ifelse(is.na(fix), 0, fix)) -> MCW_states_dates_fix
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(fix = ifelse(yday(date_actual) > fix_start_date &
                            yday(date_actual) < fix_end_date, 1, 0)) %>% 
      mutate(fix = ifelse(is.na(fix), 0, fix)) -> states_dates_fix
    
    # within season upstream vs. season downstream, calculate proportion of time to fix
    season0_upstream_fixprop <- rep(0, length(model_states))
    season0_downstream_fixprop <- rep(0, length(model_states))
    season1_upstream_fixprop <- rep(0, length(model_states))
    season1_downstream_fixprop <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # proportion of season0, upstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(MCW_states_dates_fix, state == i & direction == "upstream" & season == 0)) == 0){
          season0_upstream_fixprop[i] <- 0
        } else {
          season0_upstream_fixprop[i] <- sum(subset(MCW_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
            length(subset(MCW_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
        }
      } else {
        season0_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
      }
      # proportion of season1, upstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(MCW_states_dates_fix, state == i & direction == "upstream" & season == 1)) == 0){
          season1_upstream_fixprop[i] <- 0
        } else {
          season1_upstream_fixprop[i] <- sum(subset(MCW_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
            length(subset(MCW_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)
        }
      } else {
        season1_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1 )$fix)
      }
      # proportion of season0, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(MCW_states_dates_fix, state == i & direction == "downstream" & season == 0)) == 0){
          season0_downstream_fixprop[i] <- 0
        } else {
          season0_downstream_fixprop[i] <- sum(subset(MCW_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
            length(subset(MCW_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
        }
      } else {
        season0_downstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
      }
      # proportion of season1, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(MCW_states_dates_fix, state == i & direction == "downstream" & season == 1)) == 0){
          season1_downstream_fixprop[i] <- 0
        } else {
          season1_downstream_fixprop[i] <- sum(subset(MCW_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)/
            length(subset(MCW_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)
        }
      } else {
        season1_downstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1 )$fix)
      }
      
    }
    
    # Now, get the conditions outside of the fixed windows, for any state that you're fixing
    spillwindow_upstream_season0_notfixed <- rep(0, length(model_states))
    spillwindow_upstream_season1_notfixed <- rep(0, length(model_states))
    spillwindow_downstream_season0_notfixed <- rep(0, length(model_states))
    spillwindow_downstream_season1_notfixed <- rep(0, length(model_states))
    
    # Remove any states_month_day that are within those windows
    states_month_day %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(drop_date = ifelse(yday(date_actual.x) > fix_start_date &
                                  yday(date_actual.x) < fix_end_date, 1, 0)) %>% 
      mutate(drop_date = ifelse(is.na(drop_date), 0, drop_date)) %>% 
      filter(drop_date == 0) -> states_month_day_trimmed
    
    MCW_states_month_day %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(drop_date = ifelse(yday(date_actual.x) > fix_start_date &
                                  yday(date_actual.x) < fix_end_date, 1, 0)) %>% 
      mutate(drop_date = ifelse(is.na(drop_date), 0, drop_date)) %>% 
      filter(drop_date == 0) -> MCW_states_month_day_trimmed
    
    
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(MCW_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(MCW_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(MCW_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(MCW_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(MCW_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(MCW_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(MCW_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season1_notfixed[i] <- median(spillwindow_data[subset(MCW_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  # toggle what the spill value should be based on 1) how the function is run and 
  # 2) if spill is fixed, fixed spill vs. not fixed spill based on resampling fixprop
  
  # We are going to change the already created spillwindow_upstream_season1 etc.
  # vectors that contain the spillwindow values for
  # upstream/downstream, season0/season1, and each state
  
  
  # first - toggle based on whether or not we're choosing to fix spill here
  if(any(!(is.na(fix_spillwindow_states))) & any(!(is.na(fix_spillwindow_values))) & 
     any(!(is.na(fix_spillwindow_start_days))) & any(!(is.na(fix_spillwindow_end_days)))){
    
    # loop through this vector to determine if we're changing spill values
    for (k in 1:length(fix_spillwindow_states)){
      if (fix_spillwindow_states[k] == T){
        # resample to determine if we're going to use the fixed spill values
        fixspill_upstream_season0_sim <- sample(c(0,1), 1, prob = c(1-season0_upstream_fixprop[k], season0_upstream_fixprop[k]))
        fixspill_upstream_season1_sim <- sample(c(0,1), 1, prob = c(1-season1_upstream_fixprop[k], season1_upstream_fixprop[k]))
        fixspill_downstream_season0_sim <- sample(c(0,1), 1, prob = c(1-season0_downstream_fixprop[k], season0_downstream_fixprop[k]))
        fixspill_downstream_season1_sim <- sample(c(0,1), 1, prob = c(1-season1_downstream_fixprop[k], season1_downstream_fixprop[k]))
        
        # if we resampled that we are going to change to the fixed spill values, change them
        # if we resampled that we are not going to change to the fixed spill values because
        # we're not in the right time of year, then change the spill values to those
        # for the remainder of the year
        if(fixspill_upstream_season0_sim == 1){
          spillwindow_upstream_season0[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_upstream_season0[k] <- spillwindow_upstream_season0_notfixed[k]
        }
        if(fixspill_upstream_season1_sim == 1){
          spillwindow_upstream_season1[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_upstream_season1[k] <- spillwindow_upstream_season1_notfixed[k]
        }
        if(fixspill_downstream_season0_sim == 1){
          spillwindow_downstream_season0[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_downstream_season0[k] <- spillwindow_downstream_season0_notfixed[k]
        }
        if(fixspill_downstream_season1_sim == 1){
          spillwindow_downstream_season1[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_downstream_season1[k] <- spillwindow_downstream_season1_notfixed[k]
        }
        
      } else {
        # do nothing!
      }
      
    }
    
    
    # If we're not choosing to fix spill, we don't need to do anything - just leave spillwindow vectors as is  
  } else {
    # do nothing!
  }
  
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      
      if (upstream_season_sim == 1){
        
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCW[i,possible_movements[j],iter] +
                                                                    btemp1_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_MCW[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_MCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    btemp1xorigin6_array_MCW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin6 +
                                                                    borigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_MCW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_MCW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_MCW[i,possible_movements,iter] +
                    btemp1_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_MCW[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_MCW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    btemp1xorigin6_array_MCW[i,possible_movements,iter]*temp_upstream_season1[i]*origin6 +
                    borigin1_array_MCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_MCW[i,possible_movements,iter]*origin3 +
                    borigin4_array_MCW[i,possible_movements,iter]*origin4 +
                    borigin5_array_MCW[i,possible_movements,iter]*origin5 +
                    borigin6_array_MCW[i,possible_movements,iter]*origin6))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCW[i,possible_movements[j],iter] +
                                                                    btemp0_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_MCW[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_MCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    btemp0xorigin6_array_MCW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin6 +
                                                                    borigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_MCW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_MCW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_MCW[i,possible_movements,iter] +
                    btemp0_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_MCW[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_MCW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    btemp0xorigin6_array_MCW[i,possible_movements,iter]*temp_upstream_season0[i]*origin6 +
                    borigin1_array_MCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_MCW[i,possible_movements,iter]*origin3 +
                    borigin4_array_MCW[i,possible_movements,iter]*origin4 +
                    borigin5_array_MCW[i,possible_movements,iter]*origin5 +
                    borigin6_array_MCW[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCW[i,possible_movements[j],iter] +
                                                                      btemp1_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_MCW[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_MCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      btemp1xorigin6_array_MCW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin6 +
                                                                      borigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_MCW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_MCW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_MCW[i,possible_movements,iter] +
                    btemp1_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_MCW[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_MCW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    btemp1xorigin6_array_MCW[i,possible_movements,iter]*temp_downstream_season1[i]*origin6 +
                    borigin1_array_MCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_MCW[i,possible_movements,iter]*origin3 +
                    borigin4_array_MCW[i,possible_movements,iter]*origin4 +
                    borigin5_array_MCW[i,possible_movements,iter]*origin5 +
                    borigin6_array_MCW[i,possible_movements,iter]*origin6))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCW[i,possible_movements[j],iter] +
                                                                      btemp0_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_MCW[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_MCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      btemp0xorigin6_array_MCW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin6 +
                                                                      borigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_MCW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_MCW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_MCW[i,possible_movements,iter] +
                    btemp0_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_MCW[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_MCW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    btemp0xorigin6_array_MCW[i,possible_movements,iter]*temp_downstream_season0[i]*origin6 +
                    borigin1_array_MCW[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCW[i,possible_movements,iter]*origin2 +
                    borigin3_array_MCW[i,possible_movements,iter]*origin3 +
                    borigin4_array_MCW[i,possible_movements,iter]*origin4 +
                    borigin5_array_MCW[i,possible_movements,iter]*origin5 +
                    borigin6_array_MCW[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_fixcov_times_MCH <- function(nsim,
                                                    start_state = 2, states_dates,
                                                    origin1 = 0, origin2 = 0, 
                                                    temp_data, spillwindow_data, winterspill_data,
                                                    fix_state = NA, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                                    fix_spillwindow_states = rep(F, 9),
                                                    fix_spillwindow_values = rep(NA, 9),
                                                    fix_spillwindow_start_days = rep(NA, 9),
                                                    fix_spillwindow_end_days = rep(NA, 9),
                                                    fix_winterspill_value = NA,
                                                    fix_run_year = NA){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(MCH_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state across all years
    winterspill[i] <- median(winterspill_data[,i])
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(MCH_states_dates, state == i  & direction == "upstream")$season)/length(subset(MCH_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(MCH_states_dates, state == i  & direction == "downstream")$season)/length(subset(MCH_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    # get the date numeric for each day/month of that run year
    run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
    run_year_start <- seq(ymd("2004-06-01"), ymd("2022-06-01"), by = "years")
    run_year_end <- seq(ymd("2005-05-31"), ymd("2023-05-31"), by = "years")
    run_year_numeric = seq(4, 22, 1)
    
    run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)
    
    
    
    # rowwise() %>%
    # dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_actual & run_year_end >= date_actual)$run_year) %>% 
    # filter(run_year == fix_run_year)-> MCH_states_dates_run_year
    
    # subset states_dates for origin and population to pull out just month day, and then
    # use those months and days to index temps from that year
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    dates_actual %>% 
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) -> dates_actual
    
    
    subset(run_year_df, run_year == fix_run_year) -> fix_run_year_df
    subset(dates_actual, date_actual <= fix_run_year_df$run_year_end & date_actual >= fix_run_year_df$run_year_start) %>% 
      dplyr::rename(date_run_year = date) -> dates_fix_run_year
    
    # Deal with leap years by using 2/28 twice
    if(!("02-29" %in% dates_fix_run_year$month_day)){
      subset(dates_fix_run_year, month_day == "02-28") %>% 
        mutate(month_day = "02-29") -> leap_year
      
      dates_fix_run_year %>% 
        bind_rows(., leap_year) -> dates_fix_run_year
      
    }
    
    
    
    
    
    MCH_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> MCH_states_month_day
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> states_month_day
    
    
    temp_upstream_season0 <- rep(0, length(model_states))
    temp_upstream_season1 <- rep(0, length(model_states))
    temp_downstream_season0 <- rep(0, length(model_states))
    temp_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(MCH_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(MCH_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(MCH_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(MCH_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    spillwindow_upstream_season0 <- rep(0, length(model_states))
    spillwindow_upstream_season1 <- rep(0, length(model_states))
    spillwindow_downstream_season0 <- rep(0, length(model_states))
    spillwindow_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(MCH_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(MCH_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(MCH_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(MCH_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  
  # Fix covariate values to the ones that you are setting
  if(!(is.na(fix_temp_season0_value))){
    temp_upstream_season0[fix_state] <- fix_temp_season0_value
    temp_downstream_season0[fix_state] <- fix_temp_season0_value
  }
  
  if(!(is.na(fix_temp_season1_value))){
    temp_upstream_season1[fix_state] <- fix_temp_season1_value
    temp_downstream_season1[fix_state] <- fix_temp_season1_value
  }
  
  if(!(is.na(fix_winterspill_value))){
    winterspill[fix_state] <- fix_winterspill_value
  }
  
  # Fix covariate values for specific time periods
  
  if(any(!(is.na(fix_spillwindow_states))) & any(!(is.na(fix_spillwindow_values))) & 
     any(!(is.na(fix_spillwindow_start_days))) & any(!(is.na(fix_spillwindow_end_days)))){
    
    # calculate additional values:
    # get the median value for season0 and season1, excluding the time period in which
    # you are fixing spill
    # set season0_fix and season1_fix to the value that you are fixing it to
    # then, use the simulation to select between season0_fix and season0/season1_fix and season1,
    # depending on the relative frequency of movement in the different time periods
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    
    # now loop through the mainstem states, and determine if they dates fall within
    # the fix range for each state
    to_fix <- data.frame(state = 1:9, state_to_fix = fix_spillwindow_states,
                         fix_start_date = fix_spillwindow_start_days,
                         fix_end_date = fix_spillwindow_end_days)
    
    MCH_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(fix = ifelse(yday(date_actual) > fix_start_date &
                            yday(date_actual) < fix_end_date, 1, 0)) %>% 
      mutate(fix = ifelse(is.na(fix), 0, fix)) -> MCH_states_dates_fix
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(fix = ifelse(yday(date_actual) > fix_start_date &
                            yday(date_actual) < fix_end_date, 1, 0)) %>% 
      mutate(fix = ifelse(is.na(fix), 0, fix)) -> states_dates_fix
    
    # within season upstream vs. season downstream, calculate proportion of time to fix
    season0_upstream_fixprop <- rep(0, length(model_states))
    season0_downstream_fixprop <- rep(0, length(model_states))
    season1_upstream_fixprop <- rep(0, length(model_states))
    season1_downstream_fixprop <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # proportion of season0, upstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(MCH_states_dates_fix, state == i & direction == "upstream" & season == 0)) == 0){
          season0_upstream_fixprop[i] <- 0
        } else {
          season0_upstream_fixprop[i] <- sum(subset(MCH_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
            length(subset(MCH_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
        }
      } else {
        season0_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
      }
      # proportion of season1, upstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(MCH_states_dates_fix, state == i & direction == "upstream" & season == 1)) == 0){
          season1_upstream_fixprop[i] <- 0
        } else {
          season1_upstream_fixprop[i] <- sum(subset(MCH_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
            length(subset(MCH_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)
        }
      } else {
        season1_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1 )$fix)
      }
      # proportion of season0, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(MCH_states_dates_fix, state == i & direction == "downstream" & season == 0)) == 0){
          season0_downstream_fixprop[i] <- 0
        } else {
          season0_downstream_fixprop[i] <- sum(subset(MCH_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
            length(subset(MCH_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
        }
      } else {
        season0_downstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
      }
      # proportion of season1, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(MCH_states_dates_fix, state == i & direction == "downstream" & season == 1)) == 0){
          season1_downstream_fixprop[i] <- 0
        } else {
          season1_downstream_fixprop[i] <- sum(subset(MCH_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)/
            length(subset(MCH_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)
        }
      } else {
        season1_downstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1 )$fix)
      }
      
    }
    
    # Now, get the conditions outside of the fixed windows, for any state that you're fixing
    spillwindow_upstream_season0_notfixed <- rep(0, length(model_states))
    spillwindow_upstream_season1_notfixed <- rep(0, length(model_states))
    spillwindow_downstream_season0_notfixed <- rep(0, length(model_states))
    spillwindow_downstream_season1_notfixed <- rep(0, length(model_states))
    
    # Remove any states_month_day that are within those windows
    states_month_day %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(drop_date = ifelse(yday(date_actual.x) > fix_start_date &
                                  yday(date_actual.x) < fix_end_date, 1, 0)) %>% 
      mutate(drop_date = ifelse(is.na(drop_date), 0, drop_date)) %>% 
      filter(drop_date == 0) -> states_month_day_trimmed
    
    MCH_states_month_day %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(drop_date = ifelse(yday(date_actual.x) > fix_start_date &
                                  yday(date_actual.x) < fix_end_date, 1, 0)) %>% 
      mutate(drop_date = ifelse(is.na(drop_date), 0, drop_date)) %>% 
      filter(drop_date == 0) -> MCH_states_month_day_trimmed
    
    
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(MCH_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(MCH_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(MCH_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(MCH_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(MCH_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(MCH_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(MCH_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season1_notfixed[i] <- median(spillwindow_data[subset(MCH_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  # toggle what the spill value should be based on 1) how the function is run and 
  # 2) if spill is fixed, fixed spill vs. not fixed spill based on resampling fixprop
  
  # We are going to change the already created spillwindow_upstream_season1 etc.
  # vectors that contain the spillwindow values for
  # upstream/downstream, season0/season1, and each state
  
  
  # first - toggle based on whether or not we're choosing to fix spill here
  if(any(!(is.na(fix_spillwindow_states))) & any(!(is.na(fix_spillwindow_values))) & 
     any(!(is.na(fix_spillwindow_start_days))) & any(!(is.na(fix_spillwindow_end_days)))){
    
    # loop through this vector to determine if we're changing spill values
    for (k in 1:length(fix_spillwindow_states)){
      if (fix_spillwindow_states[k] == T){
        # resample to determine if we're going to use the fixed spill values
        fixspill_upstream_season0_sim <- sample(c(0,1), 1, prob = c(1-season0_upstream_fixprop[k], season0_upstream_fixprop[k]))
        fixspill_upstream_season1_sim <- sample(c(0,1), 1, prob = c(1-season1_upstream_fixprop[k], season1_upstream_fixprop[k]))
        fixspill_downstream_season0_sim <- sample(c(0,1), 1, prob = c(1-season0_downstream_fixprop[k], season0_downstream_fixprop[k]))
        fixspill_downstream_season1_sim <- sample(c(0,1), 1, prob = c(1-season1_downstream_fixprop[k], season1_downstream_fixprop[k]))
        
        # if we resampled that we are going to change to the fixed spill values, change them
        # if we resampled that we are not going to change to the fixed spill values because
        # we're not in the right time of year, then change the spill values to those
        # for the remainder of the year
        if(fixspill_upstream_season0_sim == 1){
          spillwindow_upstream_season0[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_upstream_season0[k] <- spillwindow_upstream_season0_notfixed[k]
        }
        if(fixspill_upstream_season1_sim == 1){
          spillwindow_upstream_season1[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_upstream_season1[k] <- spillwindow_upstream_season1_notfixed[k]
        }
        if(fixspill_downstream_season0_sim == 1){
          spillwindow_downstream_season0[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_downstream_season0[k] <- spillwindow_downstream_season0_notfixed[k]
        }
        if(fixspill_downstream_season1_sim == 1){
          spillwindow_downstream_season1[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_downstream_season1[k] <- spillwindow_downstream_season1_notfixed[k]
        }
        
      } else {
        # do nothing!
      }
      
    }
    
    
    # If we're not choosing to fix spill, we don't need to do anything - just leave spillwindow vectors as is  
  } else {
    # do nothing!
  }
  
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      
      if (upstream_season_sim == 1){
        
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCH[i,possible_movements[j],iter] +
                                                                    btemp1_array_MCH[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_MCH[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_MCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_MCH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_MCH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 +
                                                                    borigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_MCH[i,possible_movements[j],iter]*origin2)/
          sum(exp(b0_array_MCH[i,possible_movements,iter] +
                    btemp1_array_MCH[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_MCH[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_MCH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_MCH[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_MCH[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    borigin1_array_MCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCH[i,possible_movements,iter]*origin2))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCH[i,possible_movements[j],iter] +
                                                                    btemp0_array_MCH[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_MCH[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_MCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_MCH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_MCH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    borigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_MCH[i,possible_movements[j],iter]*origin2)/
          sum(exp(b0_array_MCH[i,possible_movements,iter] +
                    btemp0_array_MCH[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_MCH[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_MCH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_MCH[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_MCH[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    borigin1_array_MCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCH[i,possible_movements,iter]*origin2))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCH[i,possible_movements[j],iter] +
                                                                      btemp1_array_MCH[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_MCH[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_MCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_MCH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_MCH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      borigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_MCH[i,possible_movements[j],iter]*origin2)/
          sum(exp(b0_array_MCH[i,possible_movements,iter] +
                    btemp1_array_MCH[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_MCH[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_MCH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_MCH[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_MCH[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    borigin1_array_MCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCH[i,possible_movements,iter]*origin2))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCH[i,possible_movements[j],iter] +
                                                                      btemp0_array_MCH[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_MCH[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_MCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_MCH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_MCH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      borigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_MCH[i,possible_movements[j],iter]*origin2)/
          sum(exp(b0_array_MCH[i,possible_movements,iter] +
                    btemp0_array_MCH[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_MCH[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_MCH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_MCH[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_MCH[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    borigin1_array_MCH[i,possible_movements,iter]*origin1 +
                    borigin2_array_MCH[i,possible_movements,iter]*origin2))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_fixcov_times_SRW <- function(nsim,
                                                    start_state = 2, states_dates,
                                                    origin1 = 0, origin2 = 0, origin3 = 0,
                                                    origin4 = 0, origin5 = 0, origin6 = 0,
                                                    temp_data, spillwindow_data, winterspill_data,
                                                    fix_state = NA, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                                    fix_spillwindow_states = rep(F, 9),
                                                    fix_spillwindow_values = rep(NA, 9),
                                                    fix_spillwindow_start_days = rep(NA, 9),
                                                    fix_spillwindow_end_days = rep(NA, 9),
                                                    fix_winterspill_value = NA,
                                                    fix_run_year = NA){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(SRW_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state across all years
    winterspill[i] <- median(winterspill_data[,i])
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRW_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRW_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRW_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRW_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRW_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    # get the date numeric for each day/month of that run year
    run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
    run_year_start <- seq(ymd("2004-06-01"), ymd("2022-06-01"), by = "years")
    run_year_end <- seq(ymd("2005-05-31"), ymd("2023-05-31"), by = "years")
    run_year_numeric = seq(4, 22, 1)
    
    run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)
    
    
    
    # rowwise() %>%
    # dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_actual & run_year_end >= date_actual)$run_year) %>% 
    # filter(run_year == fix_run_year)-> SRW_states_dates_run_year
    
    # subset states_dates for origin and population to pull out just month day, and then
    # use those months and days to index temps from that year
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    dates_actual %>% 
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) -> dates_actual
    
    
    subset(run_year_df, run_year == fix_run_year) -> fix_run_year_df
    subset(dates_actual, date_actual <= fix_run_year_df$run_year_end & date_actual >= fix_run_year_df$run_year_start) %>% 
      dplyr::rename(date_run_year = date) -> dates_fix_run_year
    
    # Deal with leap years by using 2/28 twice
    if(!("02-29" %in% dates_fix_run_year$month_day)){
      subset(dates_fix_run_year, month_day == "02-28") %>% 
        mutate(month_day = "02-29") -> leap_year
      
      dates_fix_run_year %>% 
        bind_rows(., leap_year) -> dates_fix_run_year
      
    }
    
    
    
    
    
    SRW_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> SRW_states_month_day
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> states_month_day
    
    
    temp_upstream_season0 <- rep(0, length(model_states))
    temp_upstream_season1 <- rep(0, length(model_states))
    temp_downstream_season0 <- rep(0, length(model_states))
    temp_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRW_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRW_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRW_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRW_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRW_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    spillwindow_upstream_season0 <- rep(0, length(model_states))
    spillwindow_upstream_season1 <- rep(0, length(model_states))
    spillwindow_downstream_season0 <- rep(0, length(model_states))
    spillwindow_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRW_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRW_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRW_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRW_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRW_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  
  # Fix covariate values to the ones that you are setting
  if(!(is.na(fix_temp_season0_value))){
    temp_upstream_season0[fix_state] <- fix_temp_season0_value
    temp_downstream_season0[fix_state] <- fix_temp_season0_value
  }
  
  if(!(is.na(fix_temp_season1_value))){
    temp_upstream_season1[fix_state] <- fix_temp_season1_value
    temp_downstream_season1[fix_state] <- fix_temp_season1_value
  }
  
  if(!(is.na(fix_winterspill_value))){
    winterspill[fix_state] <- fix_winterspill_value
  }
  
  # Fix covariate values for specific time periods
  
  if(any(!(is.na(fix_spillwindow_states))) & any(!(is.na(fix_spillwindow_values))) & 
     any(!(is.na(fix_spillwindow_start_days))) & any(!(is.na(fix_spillwindow_end_days)))){
    
    # calculate additional values:
    # get the median value for season0 and season1, excluding the time period in which
    # you are fixing spill
    # set season0_fix and season1_fix to the value that you are fixing it to
    # then, use the simulation to select between season0_fix and season0/season1_fix and season1,
    # depending on the relative frequency of movement in the different time periods
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    
    # now loop through the mainstem states, and determine if they dates fall within
    # the fix range for each state
    to_fix <- data.frame(state = 1:9, state_to_fix = fix_spillwindow_states,
                         fix_start_date = fix_spillwindow_start_days,
                         fix_end_date = fix_spillwindow_end_days)
    
    SRW_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(fix = ifelse(yday(date_actual) > fix_start_date &
                            yday(date_actual) < fix_end_date, 1, 0)) %>% 
      mutate(fix = ifelse(is.na(fix), 0, fix)) -> SRW_states_dates_fix
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(fix = ifelse(yday(date_actual) > fix_start_date &
                            yday(date_actual) < fix_end_date, 1, 0)) %>% 
      mutate(fix = ifelse(is.na(fix), 0, fix)) -> states_dates_fix
    
    # within season upstream vs. season downstream, calculate proportion of time to fix
    season0_upstream_fixprop <- rep(0, length(model_states))
    season0_downstream_fixprop <- rep(0, length(model_states))
    season1_upstream_fixprop <- rep(0, length(model_states))
    season1_downstream_fixprop <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # proportion of season0, upstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRW_states_dates_fix, state == i & direction == "upstream" & season == 0)) == 0){
          season0_upstream_fixprop[i] <- 0
        } else {
          season0_upstream_fixprop[i] <- sum(subset(SRW_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
            length(subset(SRW_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
        }
      } else {
        season0_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
      }
      # proportion of season1, upstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_states_dates_fix, state == i & direction == "upstream" & season == 1)) == 0){
          season1_upstream_fixprop[i] <- 0
        } else {
          season1_upstream_fixprop[i] <- sum(subset(SRW_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
            length(subset(SRW_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)
        }
      } else {
        season1_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1 )$fix)
      }
      # proportion of season0, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_states_dates_fix, state == i & direction == "downstream" & season == 0)) == 0){
          season0_downstream_fixprop[i] <- 0
        } else {
          season0_downstream_fixprop[i] <- sum(subset(SRW_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
            length(subset(SRW_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
        }
      } else {
        season0_downstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
      }
      # proportion of season1, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_states_dates_fix, state == i & direction == "downstream" & season == 1)) == 0){
          season1_downstream_fixprop[i] <- 0
        } else {
          season1_downstream_fixprop[i] <- sum(subset(SRW_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)/
            length(subset(SRW_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)
        }
      } else {
        season1_downstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1 )$fix)
      }
      
    }
    
    # Now, get the conditions outside of the fixed windows, for any state that you're fixing
    spillwindow_upstream_season0_notfixed <- rep(0, length(model_states))
    spillwindow_upstream_season1_notfixed <- rep(0, length(model_states))
    spillwindow_downstream_season0_notfixed <- rep(0, length(model_states))
    spillwindow_downstream_season1_notfixed <- rep(0, length(model_states))
    
    # Remove any states_month_day that are within those windows
    states_month_day %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(drop_date = ifelse(yday(date_actual.x) > fix_start_date &
                                  yday(date_actual.x) < fix_end_date, 1, 0)) %>% 
      mutate(drop_date = ifelse(is.na(drop_date), 0, drop_date)) %>% 
      filter(drop_date == 0) -> states_month_day_trimmed
    
    SRW_states_month_day %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(drop_date = ifelse(yday(date_actual.x) > fix_start_date &
                                  yday(date_actual.x) < fix_end_date, 1, 0)) %>% 
      mutate(drop_date = ifelse(is.na(drop_date), 0, drop_date)) %>% 
      filter(drop_date == 0) -> SRW_states_month_day_trimmed
    
    
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRW_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(SRW_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(SRW_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(SRW_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season1_notfixed[i] <- median(spillwindow_data[subset(SRW_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  # toggle what the spill value should be based on 1) how the function is run and 
  # 2) if spill is fixed, fixed spill vs. not fixed spill based on resampling fixprop
  
  # We are going to change the already created spillwindow_upstream_season1 etc.
  # vectors that contain the spillwindow values for
  # upstream/downstream, season0/season1, and each state
  
  
  # first - toggle based on whether or not we're choosing to fix spill here
  if(any(!(is.na(fix_spillwindow_states))) & any(!(is.na(fix_spillwindow_values))) & 
     any(!(is.na(fix_spillwindow_start_days))) & any(!(is.na(fix_spillwindow_end_days)))){
    
    # loop through this vector to determine if we're changing spill values
    for (k in 1:length(fix_spillwindow_states)){
      if (fix_spillwindow_states[k] == T){
        # resample to determine if we're going to use the fixed spill values
        fixspill_upstream_season0_sim <- sample(c(0,1), 1, prob = c(1-season0_upstream_fixprop[k], season0_upstream_fixprop[k]))
        fixspill_upstream_season1_sim <- sample(c(0,1), 1, prob = c(1-season1_upstream_fixprop[k], season1_upstream_fixprop[k]))
        fixspill_downstream_season0_sim <- sample(c(0,1), 1, prob = c(1-season0_downstream_fixprop[k], season0_downstream_fixprop[k]))
        fixspill_downstream_season1_sim <- sample(c(0,1), 1, prob = c(1-season1_downstream_fixprop[k], season1_downstream_fixprop[k]))
        
        # if we resampled that we are going to change to the fixed spill values, change them
        # if we resampled that we are not going to change to the fixed spill values because
        # we're not in the right time of year, then change the spill values to those
        # for the remainder of the year
        if(fixspill_upstream_season0_sim == 1){
          spillwindow_upstream_season0[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_upstream_season0[k] <- spillwindow_upstream_season0_notfixed[k]
        }
        if(fixspill_upstream_season1_sim == 1){
          spillwindow_upstream_season1[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_upstream_season1[k] <- spillwindow_upstream_season1_notfixed[k]
        }
        if(fixspill_downstream_season0_sim == 1){
          spillwindow_downstream_season0[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_downstream_season0[k] <- spillwindow_downstream_season0_notfixed[k]
        }
        if(fixspill_downstream_season1_sim == 1){
          spillwindow_downstream_season1[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_downstream_season1[k] <- spillwindow_downstream_season1_notfixed[k]
        }
        
      } else {
        # do nothing!
      }
      
    }
    
    
    # If we're not choosing to fix spill, we don't need to do anything - just leave spillwindow vectors as is  
  } else {
    # do nothing!
  }
  
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      
      if (upstream_season_sim == 1){
        
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRW[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    btemp1xorigin6_array_SRW[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin6 +
                                                                    borigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW[i,possible_movements,iter] +
                    btemp1_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRW[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW[i,possible_movements,iter]*temp_upstream_season1[i]*origin6 +
                    borigin1_array_SRW[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW[i,possible_movements,iter]*origin6))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRW[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    btemp0xorigin6_array_SRW[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin6 +
                                                                    borigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW[i,possible_movements,iter] +
                    btemp0_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRW[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW[i,possible_movements,iter]*temp_upstream_season0[i]*origin6 +
                    borigin1_array_SRW[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRW[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      btemp1xorigin6_array_SRW[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin6 +
                                                                      borigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW[i,possible_movements,iter] +
                    btemp1_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRW[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRW[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW[i,possible_movements,iter]*temp_downstream_season1[i]*origin6 +
                    borigin1_array_SRW[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW[i,possible_movements,iter]*origin6))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRW[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRW[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      btemp0xorigin6_array_SRW[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin6 +
                                                                      borigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW[i,possible_movements,iter] +
                    btemp0_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRW[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRW[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW[i,possible_movements,iter]*temp_downstream_season0[i]*origin6 +
                    borigin1_array_SRW[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_fixcov_times_SRH <- function(nsim,
                                                    start_state = 2, states_dates,
                                                    origin1 = 0, origin2 = 0, origin3 = 0,
                                                    origin4 = 0, origin5 = 0,
                                                    temp_data, spillwindow_data, winterspill_data,
                                                    fix_state = NA, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                                    fix_spillwindow_states = rep(F, 9),
                                                    fix_spillwindow_values = rep(NA, 9),
                                                    fix_spillwindow_start_days = rep(NA, 9),
                                                    fix_spillwindow_end_days = rep(NA, 9),
                                                    fix_winterspill_value = NA,
                                                    fix_run_year = NA){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(SRH_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state across all years
    winterspill[i] <- median(winterspill_data[,i])
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  # 2024-07-30 edit: We're going to include temp0 (not just temp1!) by selecting a season
  # of movement based on the proportion of fish that moved from that state in that season
  # 2024-08-02 edit: We're also going to change the temp values to not just be temp_upstream and temp_downstream,
  # but temp_upstream_season0, temp_upstream_season1, temp_downstream_season0, temp_downstream_season1
  # makes more sense for it to be divided this way, since one median temp is going to reflect primarily season1
  
  season_upstream <- rep(0, length(model_states))
  season_downstream <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRH_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRH_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRH_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRH_states_dates, state == i  & direction == "downstream")$season)
      }
    } else {
      season_downstream[i] <- sum(subset(states_dates, state == i  & direction == "downstream")$season)/length(subset(states_dates, state == i  & direction == "downstream")$season)
    }
    
  }
  
  
  temp_upstream_season0 <- rep(0, length(model_states))
  temp_upstream_season1 <- rep(0, length(model_states))
  temp_downstream_season0 <- rep(0, length(model_states))
  temp_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      temp_downstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  
  spillwindow_upstream_season0 <- rep(0, length(model_states))
  spillwindow_upstream_season1 <- rep(0, length(model_states))
  spillwindow_downstream_season0 <- rep(0, length(model_states))
  spillwindow_downstream_season1 <- rep(0, length(model_states))
  
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    # approach 2: Take the median conditions for each state, by direction and season
    # season 0
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 0)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRH_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    # get the date numeric for each day/month of that run year
    run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
    run_year_start <- seq(ymd("2004-06-01"), ymd("2022-06-01"), by = "years")
    run_year_end <- seq(ymd("2005-05-31"), ymd("2023-05-31"), by = "years")
    run_year_numeric = seq(4, 22, 1)
    
    run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)
    
    
    
    # rowwise() %>%
    # dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_actual & run_year_end >= date_actual)$run_year) %>% 
    # filter(run_year == fix_run_year)-> SRH_states_dates_run_year
    
    # subset states_dates for origin and population to pull out just month day, and then
    # use those months and days to index temps from that year
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    dates_actual %>% 
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) -> dates_actual
    
    
    subset(run_year_df, run_year == fix_run_year) -> fix_run_year_df
    subset(dates_actual, date_actual <= fix_run_year_df$run_year_end & date_actual >= fix_run_year_df$run_year_start) %>% 
      dplyr::rename(date_run_year = date) -> dates_fix_run_year
    
    # Deal with leap years by using 2/28 twice
    if(!("02-29" %in% dates_fix_run_year$month_day)){
      subset(dates_fix_run_year, month_day == "02-28") %>% 
        mutate(month_day = "02-29") -> leap_year
      
      dates_fix_run_year %>% 
        bind_rows(., leap_year) -> dates_fix_run_year
      
    }
    
    
    
    
    
    SRH_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> SRH_states_month_day
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> states_month_day
    
    
    temp_upstream_season0 <- rep(0, length(model_states))
    temp_upstream_season1 <- rep(0, length(model_states))
    temp_downstream_season0 <- rep(0, length(model_states))
    temp_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRH_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRH_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRH_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRH_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRH_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    spillwindow_upstream_season0 <- rep(0, length(model_states))
    spillwindow_upstream_season1 <- rep(0, length(model_states))
    spillwindow_downstream_season0 <- rep(0, length(model_states))
    spillwindow_downstream_season1 <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRH_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRH_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRH_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRH_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRH_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  
  # Fix covariate values to the ones that you are setting
  if(!(is.na(fix_temp_season0_value))){
    temp_upstream_season0[fix_state] <- fix_temp_season0_value
    temp_downstream_season0[fix_state] <- fix_temp_season0_value
  }
  
  if(!(is.na(fix_temp_season1_value))){
    temp_upstream_season1[fix_state] <- fix_temp_season1_value
    temp_downstream_season1[fix_state] <- fix_temp_season1_value
  }
  
  if(!(is.na(fix_winterspill_value))){
    winterspill[fix_state] <- fix_winterspill_value
  }
  
  # Fix covariate values for specific time periods
  
  if(any(!(is.na(fix_spillwindow_states))) & any(!(is.na(fix_spillwindow_values))) & 
     any(!(is.na(fix_spillwindow_start_days))) & any(!(is.na(fix_spillwindow_end_days)))){
    
    # calculate additional values:
    # get the median value for season0 and season1, excluding the time period in which
    # you are fixing spill
    # set season0_fix and season1_fix to the value that you are fixing it to
    # then, use the simulation to select between season0_fix and season0/season1_fix and season1,
    # depending on the relative frequency of movement in the different time periods
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    
    # now loop through the mainstem states, and determine if they dates fall within
    # the fix range for each state
    to_fix <- data.frame(state = 1:9, state_to_fix = fix_spillwindow_states,
                         fix_start_date = fix_spillwindow_start_days,
                         fix_end_date = fix_spillwindow_end_days)
    
    SRH_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(fix = ifelse(yday(date_actual) > fix_start_date &
                            yday(date_actual) < fix_end_date, 1, 0)) %>% 
      mutate(fix = ifelse(is.na(fix), 0, fix)) -> SRH_states_dates_fix
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(fix = ifelse(yday(date_actual) > fix_start_date &
                            yday(date_actual) < fix_end_date, 1, 0)) %>% 
      mutate(fix = ifelse(is.na(fix), 0, fix)) -> states_dates_fix
    
    # within season upstream vs. season downstream, calculate proportion of time to fix
    season0_upstream_fixprop <- rep(0, length(model_states))
    season0_downstream_fixprop <- rep(0, length(model_states))
    season1_upstream_fixprop <- rep(0, length(model_states))
    season1_downstream_fixprop <- rep(0, length(model_states))
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # proportion of season0, upstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRH_states_dates_fix, state == i & direction == "upstream" & season == 0)) == 0){
          season0_upstream_fixprop[i] <- 0
        } else {
          season0_upstream_fixprop[i] <- sum(subset(SRH_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
            length(subset(SRH_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
        }
      } else {
        season0_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
      }
      # proportion of season1, upstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_states_dates_fix, state == i & direction == "upstream" & season == 1)) == 0){
          season1_upstream_fixprop[i] <- 0
        } else {
          season1_upstream_fixprop[i] <- sum(subset(SRH_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
            length(subset(SRH_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)
        }
      } else {
        season1_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1 )$fix)
      }
      # proportion of season0, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_states_dates_fix, state == i & direction == "downstream" & season == 0)) == 0){
          season0_downstream_fixprop[i] <- 0
        } else {
          season0_downstream_fixprop[i] <- sum(subset(SRH_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
            length(subset(SRH_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
        }
      } else {
        season0_downstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
      }
      # proportion of season1, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_states_dates_fix, state == i & direction == "downstream" & season == 1)) == 0){
          season1_downstream_fixprop[i] <- 0
        } else {
          season1_downstream_fixprop[i] <- sum(subset(SRH_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)/
            length(subset(SRH_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)
        }
      } else {
        season1_downstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1 )$fix)
      }
      
    }
    
    # Now, get the conditions outside of the fixed windows, for any state that you're fixing
    spillwindow_upstream_season0_notfixed <- rep(0, length(model_states))
    spillwindow_upstream_season1_notfixed <- rep(0, length(model_states))
    spillwindow_downstream_season0_notfixed <- rep(0, length(model_states))
    spillwindow_downstream_season1_notfixed <- rep(0, length(model_states))
    
    # Remove any states_month_day that are within those windows
    states_month_day %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(drop_date = ifelse(yday(date_actual.x) > fix_start_date &
                                  yday(date_actual.x) < fix_end_date, 1, 0)) %>% 
      mutate(drop_date = ifelse(is.na(drop_date), 0, drop_date)) %>% 
      filter(drop_date == 0) -> states_month_day_trimmed
    
    SRH_states_month_day %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(drop_date = ifelse(yday(date_actual.x) > fix_start_date &
                                  yday(date_actual.x) < fix_end_date, 1, 0)) %>% 
      mutate(drop_date = ifelse(is.na(drop_date), 0, drop_date)) %>% 
      filter(drop_date == 0) -> SRH_states_month_day_trimmed
    
    
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRH_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(SRH_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(SRH_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(SRH_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season1_notfixed[i] <- median(spillwindow_data[subset(SRH_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season1_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 1)$date_run_year,i])
      }
      
    }
    
    
  }
  
  # toggle what the spill value should be based on 1) how the function is run and 
  # 2) if spill is fixed, fixed spill vs. not fixed spill based on resampling fixprop
  
  # We are going to change the already created spillwindow_upstream_season1 etc.
  # vectors that contain the spillwindow values for
  # upstream/downstream, season0/season1, and each state
  
  
  # first - toggle based on whether or not we're choosing to fix spill here
  if(any(!(is.na(fix_spillwindow_states))) & any(!(is.na(fix_spillwindow_values))) & 
     any(!(is.na(fix_spillwindow_start_days))) & any(!(is.na(fix_spillwindow_end_days)))){
    
    # loop through this vector to determine if we're changing spill values
    for (k in 1:length(fix_spillwindow_states)){
      if (fix_spillwindow_states[k] == T){
        # resample to determine if we're going to use the fixed spill values
        fixspill_upstream_season0_sim <- sample(c(0,1), 1, prob = c(1-season0_upstream_fixprop[k], season0_upstream_fixprop[k]))
        fixspill_upstream_season1_sim <- sample(c(0,1), 1, prob = c(1-season1_upstream_fixprop[k], season1_upstream_fixprop[k]))
        fixspill_downstream_season0_sim <- sample(c(0,1), 1, prob = c(1-season0_downstream_fixprop[k], season0_downstream_fixprop[k]))
        fixspill_downstream_season1_sim <- sample(c(0,1), 1, prob = c(1-season1_downstream_fixprop[k], season1_downstream_fixprop[k]))
        
        # if we resampled that we are going to change to the fixed spill values, change them
        # if we resampled that we are not going to change to the fixed spill values because
        # we're not in the right time of year, then change the spill values to those
        # for the remainder of the year
        if(fixspill_upstream_season0_sim == 1){
          spillwindow_upstream_season0[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_upstream_season0[k] <- spillwindow_upstream_season0_notfixed[k]
        }
        if(fixspill_upstream_season1_sim == 1){
          spillwindow_upstream_season1[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_upstream_season1[k] <- spillwindow_upstream_season1_notfixed[k]
        }
        if(fixspill_downstream_season0_sim == 1){
          spillwindow_downstream_season0[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_downstream_season0[k] <- spillwindow_downstream_season0_notfixed[k]
        }
        if(fixspill_downstream_season1_sim == 1){
          spillwindow_downstream_season1[k] <- fix_spillwindow_values[k]
        } else {
          spillwindow_downstream_season1[k] <- spillwindow_downstream_season1_notfixed[k]
        }
        
      } else {
        # do nothing!
      }
      
    }
    
    
    # If we're not choosing to fix spill, we don't need to do anything - just leave spillwindow vectors as is  
  } else {
    # do nothing!
  }
  
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      
      if (upstream_season_sim == 1){
        
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRH[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRH[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    borigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH[i,possible_movements,iter] +
                    btemp1_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRH[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    borigin1_array_SRH[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH[i,possible_movements,iter]*origin5))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRH[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRH[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    borigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH[i,possible_movements,iter] +
                    btemp0_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRH[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    borigin1_array_SRH[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH[i,possible_movements,iter]*origin5))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRH[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRH[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      borigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH[i,possible_movements,iter] +
                    btemp1_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRH[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRH[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    borigin1_array_SRH[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH[i,possible_movements,iter]*origin5))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRH[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRH[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRH[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      borigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH[i,possible_movements,iter] +
                    btemp0_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRH[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRH[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    borigin1_array_SRH[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH[i,possible_movements,iter]*origin5))
      }
      
      
      
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}



#### Plotting function ####
# ORDER THE STATES FOR PLOTTING
states_order_for_plot <- gsub(" Mouth| Upstream", "", model_states)
states_order_for_plot <- states_order_for_plot[!(duplicated(states_order_for_plot))]
# Make a couple of changes to make them be in the order from most downstream to most upstream
states_order_for_plot[10] <- "Hood River"
states_order_for_plot[11] <- "Fifteenmile Creek"
states_order_for_plot[12] <- "Deschutes River"
states_order_for_plot[13] <- "John Day River"
states_order_for_plot[15] <- "Walla Walla River"
states_order_for_plot[16] <- "Yakima River"
states_order_for_plot[19] <- "Methow River"
states_order_for_plot[20] <- "Okanogan River"

states_order_for_plot[16:29] <- states_order_for_plot[15:28]
states_order_for_plot[15] <- "BON to MCN other tributaries"
states_order_for_plot[23:29] <- states_order_for_plot[22:28]
states_order_for_plot[22] <- "Upstream WEL other tributaries"
states_order_for_plot[29] <- "loss"

states_order_for_plot[24] <- "Clearwater River"
states_order_for_plot[25] <- "Asotin Creek"
states_order_for_plot[26] <- "Grande Ronde River"
states_order_for_plot[27] <- "Salmon River"


plot_final_fate_rear_type <- function(ff_comp, natal_origin){
  rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")
  
  ff_comp$state <- fct_rev(factor(ff_comp$state, levels = states_order_for_plot))
  
  ff_comp_plot <- ggplot(ff_comp, aes(x = state, y = `0.5`, ymin = `0.025`, ymax = `0.975`, color = rear_type)) +
    geom_point(size = 3.5, shape = 18, position=position_dodge(width=0.5)) +
    geom_linerange(position=position_dodge(width=0.5)) +
    coord_flip() +
    ylab("Final Distribution Probability") +
    xlab("Model State") +
    # ggtitle(" ") +
    # Create a scale common to all
    scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_color_manual(values = rear_colors) +
    theme(plot.title = element_text(size = 12),
          # axis.text.y = element_text(color = rev(state_significance_colors)),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12)) +
    ggtitle(natal_origin)
  
  return(ff_comp_plot)
  
}

#### Set conditions for simulations ####
ff_iter <- 1000
ff_nsim <- 1000