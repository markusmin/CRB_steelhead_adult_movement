# 06 stan data prep: Upper Columbia Hatchery

# This script doesn't do any model fitting - it just creates the input files for stan

library(cmdstanr)
library(posterior)
library(tidyverse)
library(lubridate)
# library(here)

#### SECTION 1: MODEL SETUP ####
# Create a transition matrix of 1s and 0s for movements from (rows) to (columns)

# Create a starting matrix of sites, with the first occasion populated (BON to MCN)
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

nstates <- length(model_states)
# 41 states

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

# v3: Create a manual way to determine which parameters to fix
subset(origin_param_map, natal_origin %in% c("Wenatchee River", "Entiat River", 
                                             "Okanogan River", "Methow River")) -> UC_origin_param_map

UC_origin_param_map$home_mainstem_state <- c(5,6,7,7)

mainstem_states <- 1:9
DPS_mainstem_states <- 4:7

# find where the parameter is declared
for (i in 1:nrow(subset(UC_origin_param_map, !(is.na(hatchery))))){
  for (j in 1:length(DPS_mainstem_states)){
    if(DPS_mainstem_states[j] > subset(UC_origin_param_map, !(is.na(hatchery)))$home_mainstem_state[i]){
      print(paste0("real<lower=0> sigma_yearxorigin", subset(UC_origin_param_map, !(is.na(hatchery)))$hatchery[i], "_vector_", DPS_mainstem_states[j], "_", DPS_mainstem_states[j]-1))
    }
    
  }
}

# find where the prior is set
for (i in 1:nrow(subset(UC_origin_param_map, !(is.na(hatchery))))){
  for (j in 1:length(DPS_mainstem_states)){
    if(DPS_mainstem_states[j] > subset(UC_origin_param_map, !(is.na(hatchery)))$home_mainstem_state[i]){
      print(paste0("sigma_yearxorigin", subset(UC_origin_param_map, !(is.na(hatchery)))$hatchery[i], "_vector_", DPS_mainstem_states[j], "_", DPS_mainstem_states[j]-1, " ~ cauchy(0,1)"))
    }
    
  }
}

# declare the parameters to be zero (for transformed data section)
for (i in 1:nrow(subset(UC_origin_param_map, !(is.na(hatchery))))){
  for (j in 1:length(DPS_mainstem_states)){
    if(DPS_mainstem_states[j] > subset(UC_origin_param_map, !(is.na(hatchery)))$home_mainstem_state[i]){
      print(paste0("real sigma_yearxorigin", subset(UC_origin_param_map, !(is.na(hatchery)))$hatchery[i], "_vector_", DPS_mainstem_states[j], "_", DPS_mainstem_states[j]-1, " = 0;"))
    }
    
  }
}

# Create a rows = from, columns = to matrix for movement probabilities


transition_matrix <- matrix(0, nrow = nstates, ncol = nstates)
rownames(transition_matrix) <- model_states
colnames(transition_matrix) <- model_states

# Populate every possible option with a 1
# 1: mainstem, mouth to BON
transition_matrix["mainstem, mouth to BON", "loss"] <- 1
transition_matrix["mainstem, mouth to BON", "mainstem, BON to MCN"] <- 1


# 2: mainstem, BON to MCN
transition_matrix["mainstem, BON to MCN", "loss"] <- 1
transition_matrix["mainstem, BON to MCN", "mainstem, mouth to BON"] <- 1
transition_matrix["mainstem, BON to MCN", "mainstem, MCN to ICH or PRA"] <- 1
# transition_matrix["mainstem, BON to MCN", "Deschutes River"] <- 1
transition_matrix["mainstem, BON to MCN", "Deschutes River Mouth"] <- 1
transition_matrix["mainstem, BON to MCN", "Deschutes River Upstream"] <- 1
# transition_matrix["mainstem, BON to MCN", "John Day River"] <- 1
transition_matrix["mainstem, BON to MCN", "John Day River Mouth"] <- 1
transition_matrix["mainstem, BON to MCN", "John Day River Upstream"] <- 1
# transition_matrix["mainstem, BON to MCN", "Hood River"] <- 1
# transition_matrix["mainstem, BON to MCN", "Hood River Mouth"] <- 1
# transition_matrix["mainstem, BON to MCN", "Hood River Upstream"] <- 1
# transition_matrix["mainstem, BON to MCN", "Fifteenmile Creek"] <- 1
transition_matrix["mainstem, BON to MCN", "Fifteenmile Creek Mouth"] <- 1
transition_matrix["mainstem, BON to MCN", "Fifteenmile Creek Upstream"] <- 1
# transition_matrix["mainstem, BON to MCN", "Umatilla River"] <- 1
transition_matrix["mainstem, BON to MCN", "Umatilla River Mouth"] <- 1
transition_matrix["mainstem, BON to MCN", "Umatilla River Upstream"] <- 1
transition_matrix["mainstem, BON to MCN", "BON to MCN other tributaries"] <- 1


# 3: mainstem, MCN to ICH or PRA
transition_matrix["mainstem, MCN to ICH or PRA", "loss"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, BON to MCN"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, ICH to LGR"] <- 1
# transition_matrix["mainstem, MCN to ICH or PRA", "Yakima River"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "Yakima River Mouth"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "Yakima River Upstream"] <- 1
# transition_matrix["mainstem, MCN to ICH or PRA", "Walla Walla River"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "Walla Walla River Mouth"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "Walla Walla River Upstream"] <- 1


# 4: mainstem, PRA to RIS
transition_matrix["mainstem, PRA to RIS", "loss"] <- 1
transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"] <- 1
transition_matrix["mainstem, PRA to RIS", "mainstem, RIS to RRE"] <- 1


# 5: mainstem, RIS to RRE
transition_matrix["mainstem, RIS to RRE", "loss"] <- 1
transition_matrix["mainstem, RIS to RRE", "mainstem, PRA to RIS"] <- 1
transition_matrix["mainstem, RIS to RRE", "mainstem, RRE to WEL"] <- 1
# transition_matrix["mainstem, RIS to RRE", "Wenatchee River"] <- 1
transition_matrix["mainstem, RIS to RRE", "Wenatchee River Mouth"] <- 1
transition_matrix["mainstem, RIS to RRE", "Wenatchee River Upstream"] <- 1


# 6: mainstem, RRE to WEL
transition_matrix["mainstem, RRE to WEL", "loss"] <- 1
transition_matrix["mainstem, RRE to WEL", "mainstem, RIS to RRE"] <- 1
transition_matrix["mainstem, RRE to WEL", "mainstem, upstream of WEL"] <- 1
# transition_matrix["mainstem, RRE to WEL", "Entiat River"] <- 1
transition_matrix["mainstem, RRE to WEL", "Entiat River Mouth"] <- 1
transition_matrix["mainstem, RRE to WEL", "Entiat River Upstream"] <- 1


# 7: mainstem, upstream of WEL
transition_matrix["mainstem, upstream of WEL", "loss"] <- 1
transition_matrix["mainstem, upstream of WEL", "mainstem, RRE to WEL"] <- 1
# transition_matrix["mainstem, upstream of WEL", "Okanogan River"] <- 1
transition_matrix["mainstem, upstream of WEL", "Okanogan River Mouth"] <- 1
transition_matrix["mainstem, upstream of WEL", "Okanogan River Upstream"] <- 1
# transition_matrix["mainstem, upstream of WEL", "Methow River"] <- 1
transition_matrix["mainstem, upstream of WEL", "Methow River Mouth"] <- 1
transition_matrix["mainstem, upstream of WEL", "Methow River Upstream"] <- 1
transition_matrix["mainstem, upstream of WEL", "Upstream WEL other tributaries"] <- 1


# 8: mainstem, ICH to LGR
transition_matrix["mainstem, ICH to LGR", "loss"] <- 1
transition_matrix["mainstem, ICH to LGR", "mainstem, MCN to ICH or PRA"] <- 1
transition_matrix["mainstem, ICH to LGR", "mainstem, upstream of LGR"] <- 1
# transition_matrix["mainstem, ICH to LGR", "Tucannon River"] <- 1
transition_matrix["mainstem, ICH to LGR", "Tucannon River Mouth"] <- 1
transition_matrix["mainstem, ICH to LGR", "Tucannon River Upstream"] <- 1


# 9: mainstem, upstream of LGR
transition_matrix["mainstem, upstream of LGR", "loss"] <- 1
transition_matrix["mainstem, upstream of LGR", "mainstem, ICH to LGR"] <- 1
# transition_matrix["mainstem, upstream of LGR", "Asotin Creek"] <- 1
transition_matrix["mainstem, upstream of LGR", "Asotin Creek Mouth"] <- 1
transition_matrix["mainstem, upstream of LGR", "Asotin Creek Upstream"] <- 1
transition_matrix["mainstem, upstream of LGR", "Clearwater River"] <- 1
transition_matrix["mainstem, upstream of LGR", "Salmon River"] <- 1
transition_matrix["mainstem, upstream of LGR", "Grande Ronde River"] <- 1
# transition_matrix["mainstem, upstream of LGR", "Imnaha River"] <- 1
transition_matrix["mainstem, upstream of LGR", "Imnaha River Mouth"] <- 1
transition_matrix["mainstem, upstream of LGR", "Imnaha River Upstream"] <- 1


# Deschutes River
# 10: Deschutes River Mouth
# 11: Deschutes River Upstream
transition_matrix["Deschutes River Mouth", "loss"] <- 1
transition_matrix["Deschutes River Mouth", "mainstem, BON to MCN"] <- 1
# transition_matrix["Deschutes River Mouth", "Deschutes River Upstream"] <- 1
transition_matrix["Deschutes River Upstream", "loss"] <- 1
transition_matrix["Deschutes River Upstream", "mainstem, BON to MCN"] <- 1
# transition_matrix["Deschutes River Upstream", "Deschutes River Mouth"] <- 1


# John Day River
# 12: John Day River Mouth
# 13: John Day River Upstream
transition_matrix["John Day River Mouth", "loss"] <- 1
transition_matrix["John Day River Mouth", "mainstem, BON to MCN"] <- 1
# transition_matrix["John Day River Mouth", "John Day River Upstream"] <- 1
transition_matrix["John Day River Upstream", "loss"] <- 1
transition_matrix["John Day River Upstream", "mainstem, BON to MCN"] <- 1
# transition_matrix["John Day River Upstream", "John Day River Mouth"] <- 1


# Hood River
# 14: Hood River Mouth
# 15: Hood River Upstream
# transition_matrix["Hood River Mouth", "loss"] <- 1
# transition_matrix["Hood River Mouth", "mainstem, BON to MCN"] <- 1
# transition_matrix["Hood River Mouth", "Hood River Upstream"] <- 1
# transition_matrix["Hood River Upstream", "loss"] <- 1
# transition_matrix["Hood River Upstream", "mainstem, BON to MCN"] <- 1
# transition_matrix["Hood River Upstream", "Hood River Mouth"] <- 1


# Fifteenmile Creek
# 14: Fifteenmile Creek Mouth
# 15: Fifteenmile Creek Upstream
transition_matrix["Fifteenmile Creek Mouth", "loss"] <- 1
transition_matrix["Fifteenmile Creek Mouth", "mainstem, BON to MCN"] <- 1
# transition_matrix["Fifteenmile Creek Mouth", "Fifteenmile Creek Upstream"] <- 1
transition_matrix["Fifteenmile Creek Upstream", "loss"] <- 1
transition_matrix["Fifteenmile Creek Upstream", "mainstem, BON to MCN"] <- 1
# transition_matrix["Fifteenmile Creek Upstream", "Fifteenmile Creek Mouth"] <- 1


# Umatilla River
# 16: Umatilla River Mouth
# 17: Umatilla River Upstream
transition_matrix["Umatilla River Mouth", "loss"] <- 1
transition_matrix["Umatilla River Mouth", "mainstem, BON to MCN"] <- 1
# transition_matrix["Umatilla River Mouth", "Umatilla River Upstream"] <- 1
transition_matrix["Umatilla River Upstream", "loss"] <- 1
transition_matrix["Umatilla River Upstream", "mainstem, BON to MCN"] <- 1
# transition_matrix["Umatilla River Upstream", "Umatilla River Mouth"] <- 1


# Yakima River
# 18: Yakima River Mouth
# 19: Yakima River Upstream
transition_matrix["Yakima River Mouth", "loss"] <- 1
transition_matrix["Yakima River Mouth", "mainstem, MCN to ICH or PRA"] <- 1
# transition_matrix["Yakima River Mouth", "Yakima River Upstream"] <- 1
transition_matrix["Yakima River Upstream", "loss"] <- 1
transition_matrix["Yakima River Upstream", "mainstem, MCN to ICH or PRA"] <- 1
# transition_matrix["Yakima River Upstream", "Yakima River Mouth"] <- 1


# Walla Walla River
# 20: Walla Walla River Mouth
# 21: Walla Walla River Upstream
transition_matrix["Walla Walla River Mouth", "loss"] <- 1
transition_matrix["Walla Walla River Mouth", "mainstem, MCN to ICH or PRA"] <- 1
# transition_matrix["Walla Walla River Mouth", "Walla Walla River Upstream"] <- 1
transition_matrix["Walla Walla River Upstream", "loss"] <- 1
transition_matrix["Walla Walla River Upstream", "mainstem, MCN to ICH or PRA"] <- 1
# transition_matrix["Walla Walla River Upstream", "Walla Walla River Mouth"] <- 1


# Wenatchee River
# 22: Wenatchee River Mouth
# 23: Wenatchee River Upstream
transition_matrix["Wenatchee River Mouth", "loss"] <- 1
transition_matrix["Wenatchee River Mouth", "mainstem, RIS to RRE"] <- 1
# transition_matrix["Wenatchee River Mouth", "Wenatchee River Upstream"] <- 1
transition_matrix["Wenatchee River Upstream", "loss"] <- 1
transition_matrix["Wenatchee River Upstream", "mainstem, RIS to RRE"] <- 1
# transition_matrix["Wenatchee River Upstream", "Wenatchee River Mouth"] <- 1


# Entiat River
# 24: Entiat River Mouth
# 25: Entiat River Upstream
transition_matrix["Entiat River Mouth", "loss"] <- 1
transition_matrix["Entiat River Mouth", "mainstem, RRE to WEL"] <- 1
# transition_matrix["Entiat River Mouth", "Entiat River Upstream"] <- 1
transition_matrix["Entiat River Upstream", "loss"] <- 1
transition_matrix["Entiat River Upstream", "mainstem, RRE to WEL"] <- 1
# transition_matrix["Entiat River Upstream", "Entiat River Mouth"] <- 1


# Okanogan River
# 26: Okanogan River Mouth
# 27: Okanogan River Upstream
transition_matrix["Okanogan River Mouth", "loss"] <- 1
transition_matrix["Okanogan River Mouth", "mainstem, upstream of WEL"] <- 1
# transition_matrix["Okanogan River Mouth", "Okanogan River Upstream"] <- 1
transition_matrix["Okanogan River Upstream", "loss"] <- 1
transition_matrix["Okanogan River Upstream", "mainstem, upstream of WEL"] <- 1
# transition_matrix["Okanogan River Upstream", "Okanogan River Mouth"] <- 1


# Methow River
# 28: Methow River Mouth
# 29: Methow River Upstream
transition_matrix["Methow River Mouth", "loss"] <- 1
transition_matrix["Methow River Mouth", "mainstem, upstream of WEL"] <- 1
# transition_matrix["Methow River Mouth", "Methow River Upstream"] <- 1
transition_matrix["Methow River Upstream", "loss"] <- 1
transition_matrix["Methow River Upstream", "mainstem, upstream of WEL"] <- 1
# transition_matrix["Methow River Upstream", "Methow River Mouth"] <- 1


# Tucannon River
# 30: Tucannon River Mouth
# 31: Tucannon River Upstream
transition_matrix["Tucannon River Mouth", "loss"] <- 1
transition_matrix["Tucannon River Mouth", "mainstem, ICH to LGR"] <- 1
# transition_matrix["Tucannon River Mouth", "Tucannon River Upstream"] <- 1
transition_matrix["Tucannon River Upstream", "loss"] <- 1
transition_matrix["Tucannon River Upstream", "mainstem, ICH to LGR"] <- 1
# transition_matrix["Tucannon River Upstream", "Tucannon River Mouth"] <- 1


# Asotin Creek
# 32: Asotin Creek Mouth
# 33: Asotin Creek Upstream
transition_matrix["Asotin Creek Mouth", "loss"] <- 1
transition_matrix["Asotin Creek Mouth", "mainstem, upstream of LGR"] <- 1
# transition_matrix["Asotin Creek Mouth", "Asotin Creek Upstream"] <- 1
transition_matrix["Asotin Creek Upstream", "loss"] <- 1
transition_matrix["Asotin Creek Upstream", "mainstem, upstream of LGR"] <- 1
# transition_matrix["Asotin Creek Upstream", "Asotin Creek Mouth"] <- 1


# 34: Clearwater River
transition_matrix["Clearwater River", "loss"] <- 1
transition_matrix["Clearwater River", "mainstem, upstream of LGR"] <- 1


# 35: Salmon River
transition_matrix["Salmon River", "loss"] <- 1
transition_matrix["Salmon River", "mainstem, upstream of LGR"] <- 1


# 36: Grande Ronde River
transition_matrix["Grande Ronde River", "loss"] <- 1
transition_matrix["Grande Ronde River", "mainstem, upstream of LGR"] <- 1


# Imnaha River
# 37: Imnaha River Mouth
# 38: Imnaha River Upstream
transition_matrix["Imnaha River Mouth", "loss"] <- 1
transition_matrix["Imnaha River Mouth", "mainstem, upstream of LGR"] <- 1
# transition_matrix["Imnaha River Mouth", "Imnaha River Upstream"] <- 1
transition_matrix["Imnaha River Upstream", "loss"] <- 1
transition_matrix["Imnaha River Upstream", "mainstem, upstream of LGR"] <- 1
# transition_matrix["Imnaha River Upstream", "Imnaha River Mouth"] <- 1

# 39: BON to MCN other tributaries
transition_matrix["BON to MCN other tributaries", "loss"] <- 1
transition_matrix["BON to MCN other tributaries", "mainstem, BON to MCN"] <- 1

# 40: Upstream WEL other tributaries
transition_matrix["Upstream WEL other tributaries", "loss"] <- 1
transition_matrix["Upstream WEL other tributaries", "mainstem, upstream of WEL"] <- 1


## Note the indices of certain movements
print(paste0(seq(1,41,1), " - ", model_states)) # this just makes it easier to tell state indexing
mainstem_indices <- seq(1,9,1)
upstream_indices <- which(grepl(" Upstream", model_states))
river_mouth_indices <- upstream_indices - 1


#### SECTION 2: PREPARE AND EXPORT DATA FOR MODELING ####
##### Step 1: Load and filter states data #####

# Load states complete
# states_complete <- read.csv(here::here("intermediate_outputs", "adults_states_complete", 
#                                        "upper_columbia_adults_states_complete.csv"), row.names = 1)
states_complete <- read.csv("intermediate_outputs/adults_states_complete/upper_columbia_adults_states_complete.csv", row.names = 1)

# load and join transport data
transport_data <- read.csv(here::here("Data", "covariate_data", "model_inputs", "transport.csv"))

states_complete %>% 
  left_join(., transport_data, by = join_by(tag_code == tag_id)) %>% 
  mutate(transport = ifelse(is.na(transport), 0, transport)) -> states_complete

### KEEP ONLY THE HATCHERY ORIGIN FISH
# tag_code_metadata <- read.csv(here::here("Data", "covariate_data", "tag_code_metadata.csv"))
tag_code_metadata <- read.csv("Data/covariate_data/tag_code_metadata.csv")
# keep only the fish that are in the dataset
tag_code_metadata <- subset(tag_code_metadata, tag_code %in% states_complete$tag_code)
# natal_origin_table <- read.csv(here::here("Data", "covariate_data", "natal_origin_table.csv"))
natal_origin_table <- read.csv("Data/covariate_data/natal_origin_table.csv")
tag_code_metadata %>% 
  left_join(., natal_origin_table, by = "release_site_name") -> tag_code_metadata


states_complete %>% 
  left_join(., dplyr::select(tag_code_metadata, tag_code, rear_type_code, natal_origin), by = "tag_code") %>% 
  subset(!(rear_type_code %in% c("W"))) -> states_complete

# Count how many we have by origin, and drop small sample sizes
states_complete %>% 
  filter(!duplicated(tag_code_2)) %>% 
  count(natal_origin) -> uppcol_hatchery_origin_table

origins_to_drop <- subset(uppcol_hatchery_origin_table, n < 350)$natal_origin
# Dropping hatchery, Entiat fish - there's only a single, lonely fish (and this is an unknown origin fish)
subset(tag_code_metadata, natal_origin %in% origins_to_drop)$tag_code -> drop_tag_codes

states_complete %>% 
  subset(!(tag_code %in% drop_tag_codes)) -> states_complete

# first create the run year df
run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14",
              "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23", "23/24", "24/25")
run_year_start <- seq(ymd_hms("2004-06-01 00:00:00"), ymd_hms("2024-06-01 00:00:00"), by = "years")
run_year_end <- seq(ymd_hms("2005-05-31 23:59:59"), ymd_hms("2025-05-31 23:59:59"), by = "years")
run_year_numeric = seq(4, 24, 1)

run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)

# First, let's join the run year df with the states complete 
# check how long this takes: about 30 seconds
# takes 
print(Sys.time())
states_complete %>% 
  rowwise() %>% # this is apparently crucial to getting this to work
  dplyr::mutate(date_time = ymd_hms(date_time)) %>% 
  dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_time & run_year_end >= date_time)$run_year) -> states_complete
print(Sys.time())

# SAVE THIS AS A MODEL INPUT - to index which matrix to use
run_year_indexing <- data.frame(run_year = run_year[1:20],
                                run_year_index = 1:20)

states_complete %>% 
  left_join(., run_year_indexing, by = "run_year") -> states_complete


# Remove any fish that have detections in the 24/25 run year
run_year_2425_fish <- unique(subset(states_complete, is.na(run_year_index))$tag_code)

states_complete %>% 
  subset(., !(tag_code %in% run_year_2425_fish)) -> states_complete

##### Step 2: Edit states data - remove upstream detections in DE years #####

# Create a data frame for when each tributary has detection efficiency capability

# Make them first each individually, then join them
deschutes_river_trib_det_eff_capability <- data.frame(state = "Deschutes River Upstream",
                                                      run_year = run_year[1:20], # ignore 24/25 to keep consistent
                                                      DE = c(rep(0,9), rep(1,6), rep(0,5)))

john_day_river_trib_det_eff_capability <- data.frame(state = "John Day River Upstream",
                                                     run_year = run_year[1:20], # ignore 24/25 to keep consistent
                                                     DE = c(rep(0,8), rep(1,12)))

# hood_river_trib_det_eff_capability <- data.frame(state = "Hood River Upstream",
#                                                  run_year = run_year[1:18], # ignore 24/25 to keep consistent
#                                                  DE = c(rep(0,9), rep(1,9)))

fifteenmile_creek_trib_det_eff_capability <- data.frame(state = "Fifteenmile Creek Upstream",
                                                        run_year = run_year[1:20], # ignore 24/25 to keep consistent
                                                        DE = c(rep(0,8), rep(1,12)))

umatilla_river_trib_det_eff_capability <- data.frame(state = "Umatilla River Upstream",
                                                     run_year = run_year[1:20], # ignore 24/25 to keep consistent
                                                     DE = c(rep(0,3), rep(1,17)))

yakima_river_trib_det_eff_capability <- data.frame(state = "Yakima River Upstream",     
                                                   run_year = run_year[1:20], # ignore 24/25 to keep consistent
                                                   # remove 04/05 run year
                                                   DE = c(0, rep(1,19)))

walla_walla_river_trib_det_eff_capability <- data.frame(state = "Walla Walla River Upstream",
                                                        run_year = run_year[1:20], # ignore 24/25 to keep consistent
                                                        DE = c(rep(0,1), rep(1,19)))

wenatchee_river_trib_det_eff_capability <- data.frame(state = "Wenatchee River Upstream",
                                                      run_year = run_year[1:18], # ignore 24/25 to keep consistent
                                                      DE = c(rep(0,7), rep(1,11)))

entiat_river_trib_det_eff_capability <- data.frame(state = "Entiat River Upstream",
                                                   run_year = run_year[1:20], # ignore 24/25 to keep consistent
                                                   DE = c(rep(0,4), rep(1,16)))

okanogan_river_trib_det_eff_capability <- data.frame(state = "Okanogan River Upstream",
                                                     run_year = run_year[1:20], # ignore 24/25 to keep consistent
                                                     DE = c(rep(0,9), rep(1,11)))

methow_river_trib_det_eff_capability <- data.frame(state = "Methow River Upstream",
                                                   run_year = run_year[1:20], # ignore 24/25 to keep consistent
                                                   DE = c(rep(0,5), rep(1,15)))

tucannon_river_trib_det_eff_capability <- data.frame(state = "Tucannon River Upstream",
                                                     run_year = run_year[1:20], # ignore 24/25 to keep consistent
                                                     DE = c(rep(0,7), rep(1,13)))

asotin_creek_trib_det_eff_capability <- data.frame(state = "Asotin Creek Upstream",      
                                                   run_year = run_year[1:20], # ignore 24/25 to keep consistent
                                                   DE = c(rep(0,7), rep(1,13)))

imnaha_river_trib_det_eff_capability <- data.frame(state = "Imnaha River Upstream",
                                                   run_year = run_year[1:20], # ignore 24/25 to keep consistent
                                                   DE = c(rep(0,7), rep(1,13)))


deschutes_river_trib_det_eff_capability %>% 
  bind_rows(., john_day_river_trib_det_eff_capability) %>% 
  # bind_rows(., hood_river_trib_det_eff_capability) %>% 
  bind_rows(., fifteenmile_creek_trib_det_eff_capability) %>% 
  bind_rows(., umatilla_river_trib_det_eff_capability) %>% 
  bind_rows(., yakima_river_trib_det_eff_capability) %>% 
  bind_rows(., walla_walla_river_trib_det_eff_capability) %>% 
  bind_rows(., wenatchee_river_trib_det_eff_capability) %>% 
  bind_rows(., entiat_river_trib_det_eff_capability) %>% 
  bind_rows(., okanogan_river_trib_det_eff_capability) %>% 
  bind_rows(., methow_river_trib_det_eff_capability) %>% 
  bind_rows(., tucannon_river_trib_det_eff_capability) %>% 
  bind_rows(., asotin_creek_trib_det_eff_capability) %>% 
  bind_rows(., imnaha_river_trib_det_eff_capability) -> trib_det_eff_capability

# join by state and run year
states_complete %>% 
  left_join(.,trib_det_eff_capability, by = c("state", "run_year")) -> states_complete


# Solution: 1) Remove any transitions that occur between mouth and upstream or upstream and mouth; 2) remove all only upstream transitions
# that occur in years where DE = 0

# KEY STEP - in all run years in upstream states where DE = 1, remove the upstream state, and remove all transitions between mouth-upstream and upstream-mouth
# get all of the upstream states in a vector
upstream_states <- model_states[upstream_indices]

# get all of the river mouth states in a vector
river_mouth_states <- model_states[river_mouth_indices]

# get all mainstem states in a vector
mainstem_states <- model_states[mainstem_indices]


states_complete %>% 
  # step 1): remove states where the tag code and the next tag code are the same (therefore same fish)
  subset(., !(tag_code == lag(tag_code) & 
                # and where the next state is an upstream state and the last state is a river mouth state
                state %in% upstream_states & lag(state) %in% river_mouth_states |
                # OR where the next state is a river mouth state and the last state is an upstream state
                tag_code == lag(tag_code) & 
                state %in% river_mouth_states & lag(state) %in% upstream_states)) %>% 
  # step 1.5): for fish that go mainstem - upstream - mainstem - remove the second two states
  subset(., !(tag_code == lag(tag_code) & 
                tag_code == lead(tag_code) &
                state %in% upstream_states & 
                lag(state) %in% mainstem_states & 
                lead(state) %in% mainstem_states |
           lag(tag_code) == tag_code & 
           lag(tag_code, n = 2) == tag_code & 
           state %in% mainstem_states & 
           lag(state) %in% upstream_states &
           lag(state, 2) %in% mainstem_states)) %>% 
  # step 2): remove only upstream states in years where DE == 1
  # keep NAs and zeros, which will remove 1s
  subset(., is.na(DE) | DE == 0 ) -> states_complete

# this introduces an issue - if a fish goes straight to an upstream state and then comes back to the mainstem in a DE = 1 year,
# we end up with consecutive mainstem states because the intervening state got removed.
# Solution: See step 1.5 above. If a fish goes mainstem - upstream - mainstem, remove the upstream state and the next mainstem state

# Now - we have to check to make sure everything makes sense
states_complete %>% 
  ungroup() %>% 
  mutate(current_tag_code = tag_code_2) %>% 
  mutate(next_tag_code = lead(tag_code_2)) %>% 
  mutate(current_state = state) %>% 
  mutate(next_state = lead(state)) %>% 
  mutate(transition = paste0(current_state, " - ", next_state)) %>% 
  subset(current_tag_code == next_tag_code) -> transitions_df
  # dplyr::select(tag_code_2, state, transition, run_year) -> transitions_df

# inspect output
table(transitions_df$transition)

##### Get values to pass as data for stan model #####

unique_tag_code_2s <- unique(states_complete$tag_code_2)

# first, get the maximum number of visits by any fish
states_complete %>% 
  # group_by(tag_code) %>% 
  # We need to use tag_code_2 because this splits apart repeat spawners
  group_by(tag_code_2) %>% 
  dplyr::count() %>% 
  as.data.frame() -> site_visits_by_tag_code

# + 1 here because we don't yet have a loss state
max_visits <- max(site_visits_by_tag_code$n) + 1

unique_tag_codes <- unique(states_complete$tag_code_2)
nfish <- length(unique_tag_codes)

# second, create an array with the dimensions number of fish, maximum number of site visits, number of states
state_data <- array(data = 0, dim = c(nstates, max_visits, nfish))

# third, populate the state data

# Add a column for the observation
states_complete %>% 
  # group_by(tag_code) %>% 
  group_by(tag_code_2) %>% 
  dplyr::mutate(order = row_number()) -> states_complete

# add a column to index by fish
states_complete %>% 
  # group_by(tag_code) %>% 
  group_by(tag_code_2) %>% 
  dplyr::mutate(tag_code_number = cur_group_id()) -> states_complete


##### Prepare final states data for model #####

# Edit 2023-07-24: Need to create new fields for windows of time - this is 
# necessary to account for implict site visits and create the maximum possible
# window that a fish was in a state

# Edit 2023-07-25: We have to do this before dropping the upstream states in DE
# years, because those detections are still informative for creating the windows,
# even if they aren't for the movement probabilities

## Note that this block takes a while
print(Sys.time())
states_complete %>% 
  mutate(prev_tag_code = lag(tag_code_2)) %>% 
  mutate(next_tag_code = lead(tag_code_2)) %>% 
  # add in next state, including loss, based on whether the next tag code is the same or not
  mutate(next_state = ifelse(next_tag_code != tag_code_2, "loss", lead(state))) %>% 
  # add in the next date (and if it's loss, NA)
  mutate(date_time_2 = ifelse(next_state == "loss", NA, ymd_hms(format(as_datetime(lead(date_time)))))) %>% 
  mutate(implicit_movement_prev_date_time = ymd_hms(ifelse(pathway == "implicit" & lag(pathway != "implicit"), format(as_datetime(lag(date_time))), NA_character_))) %>%
  fill(implicit_movement_prev_date_time, .direction = "down") %>%
  mutate(implicit_movement_next_date_time = ymd_hms(ifelse(pathway == "implicit" & lead(pathway != "implicit"), format(as_datetime(lead(date_time))), NA_character_))) %>%
  fill(implicit_movement_next_date_time, .direction = "up") %>%
  # drop all of the rows that you don't need an implicit movement date tim
  mutate(implicit_movement_prev_date_time = ifelse(pathway != "implicit", NA, ymd_hms(format(as_datetime(implicit_movement_prev_date_time))))) %>%
  mutate(implicit_movement_next_date_time = ifelse(pathway != "implicit", NA, ymd_hms(format(as_datetime(implicit_movement_next_date_time))))) %>% 
  # dplyr::select(tag_code, state, next_state, pathway, date_time, date_time_2, implicit_movement_prev_date_time, implicit_movement_next_date_time)
  mutate(window_date_time_1 = ymd_hms(format(as_datetime(ifelse(pathway == "implicit", implicit_movement_prev_date_time, date_time))))) %>% 
  mutate(window_date_time_2 = ymd_hms(format(as_datetime(ifelse(pathway == "implicit", implicit_movement_next_date_time, date_time_2))))) -> states_complete
print(Sys.time())


# Loop through to convert df to array
# This takes quite a while
print(Sys.time())
for (i in 1:nrow(states_complete)){
  # We need to skip all upstream detections in DE years
  # index by 1) which state, 2) which visit, 3) which fish
  state_data[which(model_states == states_complete[i, "state", drop = TRUE]),states_complete[i, "order", drop = TRUE], states_complete[i, "tag_code_number", drop = TRUE]] <- 1 
  
}
print(Sys.time())

# Now, add the loss state
# First, get the number of non-loss states
n_not_loss_states <- vector(length = nfish)
# Now, get the indices of the fish that already have a loss state
states_complete %>% 
  subset(state == "loss") -> trapped_fish

trapped_fish_tags <- trapped_fish$tag_code_2
# Get the indices of these
which(unique_tag_codes %in% trapped_fish_tags) -> trapped_fish_indices

for (i in 1:nfish){
  # 41 is the loss state
  n_not_loss_states[i] <- sum(state_data[,,i])
  state_data[41, n_not_loss_states[i]+1,i] <- 1
}

# remove the extra loss state for the trapped fish
for (i in 1:length(trapped_fish_indices)){
  j <- trapped_fish_indices[i]
  state_data[41, n_not_loss_states[j]+1,j] <-0
}

##### Load and reformat tributary data for detection efficiency #####

# Load tributary discharge data
# tributary_discharge_data <- read.csv(here::here("Data", "covariate_data", "model_inputs", "tributary_discharge_data_zscore.csv"))
tributary_discharge_data <- read.csv("Data/covariate_data/model_inputs/tributary_discharge_data_zscore.csv")

# Get tributary categorical data (equipment eras)
# this by run year
# we will have multiple eras, then NAs for run years without an ability to calculate detection efficiency

# Create a design matrix for all tributary data
# first columns = intercepts for eras
# last columns = discharge values
# rows = run years

# This creates the master file; we then subset it to make a design matrix for each tributary, where we replace every value except
# the relevant columns with zeros. We will then store this in one big array; each slice of the array will be one tributary,
# then rows will be run years, and columns will be the design matrix for the parameters.
# Then, in the actual model code, we will use the tributary a fish was detected in to select a slice of the array;
# then that slice (which is the design matrix for that tributary) will then be indexed by the run year to get a row vector.
# The row vector will be multiplied by the full parameter vector (all alphas and betas) to get the eta (linear predictor) for the detection efficiency GLM





# initialize the matrix, with dimensions rows = run years, columns = N eras (for alphas) + N tributaries (for betas)
# NOTE: for the Imnaha and for Fifteenmile Creek, we will not have a relationship with discharge. So we will only have an intercept, and put
# in fake data (all zeros, so that the slope term doesn't matter) for the discharge data in order to allow it to run.

# Tributaries that we have no capability to estimate detection capability: Clearwater River, Grande Ronde River, Salmon River

# Note that we are choosing the start run year depending on when in the year the site came online (e.g. 12/2011 would be
# 11/12 start year). If it's active in the spring of that run year then we accept it

# Here's how many eras we have (do it alphabetically)
# Asotin Creek 1: 11/12 - 17/18
# Asotin Creek 2: 18/19 - 23/24
# Deschutes River 1: 13/14 - 18/19
# Entiat River 1: 07/08 - 23/24
# Fifteenmile Creek: 11/12 - 18/19; river mouth site still active after this point but no upstream sites
# Imnaha River 1: 10/11 - 23/24
# John Day River 1: 12/13 - 23/24
# Methow River 1: 09/10 - 16/17
# Methow River 2: 17/18 - 23/24
# Okanogan River 1: 13/14 - 23/24
# Tucannon River 1: 10/11 - 19/20
# Tucannon River 2: 20/21 - 23/24
# Umatilla River 1: 06/07 - 13/14
# Umatilla River 2: 14/15 - 23/24 (see email from Stacy Remple for info on this - not in tributary detection efficiency RMD)
# Walla Walla River 1 (ORB only): 05/06 - 11/12
# Walla Walla River 2 (ORB and PRV): 12/13 - 14/15
# Walla Walla River 3 (PRV only): 16/17-18/19
# Walla Walla River 4 (WWB): 19/20 - 23/24
# Wenatchee River 1: 10/11 (starts January 2011) - 23/24
# Yakima River 1: 04/05 - 23/24

# Total eras: 20
# Total tribs: 13 (17 total - Clearwater, Grande Ronde, and Salmon; also Hood River has now been removed as a state)

# So the first 20 columns are the era/intercept terms, then the next 13 are the slope terms for discharge

# We then have 20 run years (04/05 through 23/24)

# Our parameter vector will also be a column vector of length 33

tributary_design_matrix <- matrix(0, nrow = 20, ncol = 33)

# This is to make it easier to tell indexing
rownames(tributary_design_matrix) <- paste0(run_year_df$run_year[1:(nrow(run_year_df)-1)], "-", seq(1,20,1))

# Now let's populate it with 1s/discharge values

# (1) Asotin Creek 1: 11/12 - 17/18
tributary_design_matrix[8:14,1] <- 1

# (2) Asotin Creek 2: 18/19 - 23/24
tributary_design_matrix[15:20,2] <- 1

# (3) Deschutes River 1: 13/14 - 18/19
tributary_design_matrix[10:15,3] <- 1

# (4) Entiat River 1: 07/08 - 23/24
tributary_design_matrix[4:20,4] <- 1

# (5) Fifteenmile Creek: 11/12 - 18/19; river mouth site still active after this point but no upstream sites
# tributary_design_matrix[8:15,5] <- 1
tributary_design_matrix[8:20,5] <- 1

# # (6) Hood River 1: 12/13 - 23/24
# tributary_design_matrix[9:20,6] <- 1

# (6) Imnaha River 1: 10/11 - 23/24
tributary_design_matrix[7:20,6] <- 1

# (7) John Day River 1: 12/13 - 23/24
tributary_design_matrix[9:20,7] <- 1

# (8) Methow River 1: 09/10 - 16/17
tributary_design_matrix[6:20,8] <- 1

# (9) Methow River 2: 17/18 - 23/24
tributary_design_matrix[14:20,9] <- 1

# (10) Okanogan River 1: 13/14 - 23/24
tributary_design_matrix[10:20,10] <- 1

# (11) Tucannon River 1: 10/11 - 19/20
tributary_design_matrix[7:16,11] <- 1

# (12) Tucannon River 2: 20/21 - 23/24
tributary_design_matrix[17:20,12] <- 1

# (13) Umatilla River 1: 06/07 - 13/14
tributary_design_matrix[3:10,13] <- 1

# (14) Umatilla River 2: 14/15 - 23/24 (see email from Stacy Remple for info on this - not in tributary detection efficiency RMD or on PTAGIS, but there was an antenna installation in 2014)
tributary_design_matrix[11:20,14] <- 1

# (15) Walla Walla River 1 (ORB solo): 05/06 - 11/12
tributary_design_matrix[2:8,15] <- 1

# (16) Walla Walla River 2 (PRV and ORB simultaneously): 12/13 - 14/15
tributary_design_matrix[9:15,16] <- 1

# (17) Walla Walla River 3 (PRV solo): 15/16-18/19
tributary_design_matrix[9:15,17] <- 1

# (18) Walla Walla River 4 (WWB): 19/20 - 23/24
tributary_design_matrix[16:20,18] <- 1

# (19) Wenatchee River 1: 10/11 (starts January 2011) - 23/24
tributary_design_matrix[7:20,19] <- 1

# (20) Yakima River 1: 04/05 - 23/24
# tributary_design_matrix[1:18,20] <- 1
# for now - remove 04/05 run year
tributary_design_matrix[2:20,20] <- 1


# Now, add discharge values in the remaining columns
tributary_design_matrix[2:20,21] <-subset(tributary_discharge_data, tributary == "Asotin Creek")$mean_discharge_zscore

tributary_design_matrix[2:20,22] <-subset(tributary_discharge_data, tributary == "Deschutes River")$mean_discharge_zscore

tributary_design_matrix[2:20,23] <-subset(tributary_discharge_data, tributary == "Entiat River")$mean_discharge_zscore

# tributary_design_matrix[2:20,24] <-subset(tributary_discharge_data, tributary == "Fifteenmile Creek")$mean_discharge_zscore
# No discharge data for Fifteenmile Creek - keep as zeros, since we don't have a complete time series

# tributary_design_matrix[2:20,26] <-subset(tributary_discharge_data, tributary == "Hood River")$mean_discharge_zscore

# tributary_design_matrix[2:10,25] <-subset(tributary_discharge_data, tributary == "Imnaha River")$mean_discharge_zscore
# Indexing is different for the Imnaha because we only have discharge data through 12/13 (13/14 data is not from the complete run year)
# Keep Imnaha as all zeros - we only have one year (12/13) with overlap between discharge data and detection capability

tributary_design_matrix[2:20,26] <-subset(tributary_discharge_data, tributary == "John Day River")$mean_discharge_zscore

tributary_design_matrix[2:20,27] <-subset(tributary_discharge_data, tributary == "Methow River")$mean_discharge_zscore

tributary_design_matrix[2:20,28] <-subset(tributary_discharge_data, tributary == "Okanogan River")$mean_discharge_zscore

tributary_design_matrix[2:20,29] <-subset(tributary_discharge_data, tributary == "Tucannon River")$mean_discharge_zscore

tributary_design_matrix[2:20,30] <-subset(tributary_discharge_data, tributary == "Umatilla River")$mean_discharge_zscore

tributary_design_matrix[2:20,31] <-subset(tributary_discharge_data, tributary == "Walla Walla River")$mean_discharge_zscore

tributary_design_matrix[2:20,32] <-subset(tributary_discharge_data, tributary == "Wenatchee River")$mean_discharge_zscore

tributary_design_matrix[2:20,33] <-subset(tributary_discharge_data, tributary == "Yakima River")$mean_discharge_zscore




# change row names just to run years
rownames(tributary_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]

# Remove the 04/05 run year
# tributary_design_matrix <- tributary_design_matrix[2:18,]
# Don't do this - just leave it in for indexing purposes

# Take the master tributary design matrix, and create design matrices for each of the tributaries individually
asotin_creek_design_matrix <- matrix(0, nrow = 20, ncol = 33)
rownames(asotin_creek_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
asotin_creek_design_matrix[,1] <- tributary_design_matrix[,1]
asotin_creek_design_matrix[,2] <- tributary_design_matrix[,2]
asotin_creek_design_matrix[,21] <- tributary_design_matrix[,21]

deschutes_river_design_matrix <- matrix(0, nrow = 20, ncol = 33)
rownames(deschutes_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
deschutes_river_design_matrix[,3] <- tributary_design_matrix[,3]
deschutes_river_design_matrix[,22] <- tributary_design_matrix[,22]

entiat_river_design_matrix <- matrix(0, nrow = 20, ncol = 33)
rownames(entiat_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
entiat_river_design_matrix[,4] <- tributary_design_matrix[,4]
entiat_river_design_matrix[,23] <- tributary_design_matrix[,23]

fifteenmile_creek_design_matrix <- matrix(0, nrow = 20, ncol = 33)
rownames(fifteenmile_creek_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
fifteenmile_creek_design_matrix[,5] <- tributary_design_matrix[,5]
fifteenmile_creek_design_matrix[,24] <- tributary_design_matrix[,24]

# hood_river_design_matrix <- matrix(0, nrow = 20, ncol = 33)
# rownames(hood_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
# hood_river_design_matrix[,6] <- tributary_design_matrix[,6]
# hood_river_design_matrix[,26] <- tributary_design_matrix[,26]

imnaha_river_design_matrix <- matrix(0, nrow = 20, ncol = 33)
rownames(imnaha_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
imnaha_river_design_matrix[,6] <- tributary_design_matrix[,6]
imnaha_river_design_matrix[,25] <- tributary_design_matrix[,25]

john_day_river_design_matrix <- matrix(0, nrow = 20, ncol = 33)
rownames(john_day_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
john_day_river_design_matrix[,7] <- tributary_design_matrix[,7]
john_day_river_design_matrix[,26] <- tributary_design_matrix[,26]

methow_river_design_matrix <- matrix(0, nrow = 20, ncol = 33)
rownames(methow_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
methow_river_design_matrix[,8] <- tributary_design_matrix[,8]
methow_river_design_matrix[,9] <- tributary_design_matrix[,9]
methow_river_design_matrix[,27] <- tributary_design_matrix[,27]

okanogan_river_design_matrix <- matrix(0, nrow = 20, ncol = 33)
rownames(okanogan_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
okanogan_river_design_matrix[,10] <- tributary_design_matrix[,10]
okanogan_river_design_matrix[,28] <- tributary_design_matrix[,28]

tucannon_river_design_matrix <- matrix(0, nrow = 20, ncol = 33)
rownames(tucannon_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
tucannon_river_design_matrix[,11] <- tributary_design_matrix[,11]
tucannon_river_design_matrix[,12] <- tributary_design_matrix[,12]
tucannon_river_design_matrix[,29] <- tributary_design_matrix[,29]

umatilla_river_design_matrix <- matrix(0, nrow = 20, ncol = 33)
rownames(umatilla_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
umatilla_river_design_matrix[,13] <- tributary_design_matrix[,13]
umatilla_river_design_matrix[,14] <- tributary_design_matrix[,14]
umatilla_river_design_matrix[,30] <- tributary_design_matrix[,30]

walla_walla_river_design_matrix <- matrix(0, nrow = 20, ncol = 33)
rownames(walla_walla_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
walla_walla_river_design_matrix[,15] <- tributary_design_matrix[,15]
walla_walla_river_design_matrix[,16] <- tributary_design_matrix[,16]
walla_walla_river_design_matrix[,17] <- tributary_design_matrix[,17]
walla_walla_river_design_matrix[,18] <- tributary_design_matrix[,18]
walla_walla_river_design_matrix[,31] <- tributary_design_matrix[,31]

wenatchee_river_design_matrix <- matrix(0, nrow = 20, ncol = 33)
rownames(wenatchee_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
wenatchee_river_design_matrix[,19] <- tributary_design_matrix[,19]
wenatchee_river_design_matrix[,32] <- tributary_design_matrix[,32]

yakima_river_design_matrix <- matrix(0, nrow = 20, ncol = 33)
rownames(yakima_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
yakima_river_design_matrix[,20] <- tributary_design_matrix[,20]
yakima_river_design_matrix[,33] <- tributary_design_matrix[,33]

# Now take these matrices and store them in one big array
# indices: rows = run years, columns = parameters/design matrix, slices = tributaries
# Note: this will have the same number of slices as possible state transitions (40 - which is all transitions minus loss).
# But only 13 of these slices (corresponding to the 13 tributaries that have detection efficiency capability)
# will have any non-zero values.
tributary_design_matrices_array <- array(0, dim = c(20,33,40))

# Note that these are indices of only the mouth states (since these are the states that we are calculating detection efficiency at)
tributary_design_matrices_array[,,32] <- asotin_creek_design_matrix
tributary_design_matrices_array[,,10] <- deschutes_river_design_matrix
tributary_design_matrices_array[,,24] <- entiat_river_design_matrix
tributary_design_matrices_array[,,14] <- fifteenmile_creek_design_matrix
# tributary_design_matrices_array[,,14] <- hood_river_design_matrix
tributary_design_matrices_array[,,37] <- imnaha_river_design_matrix
tributary_design_matrices_array[,,12] <- john_day_river_design_matrix
tributary_design_matrices_array[,,28] <- methow_river_design_matrix
tributary_design_matrices_array[,,26] <- okanogan_river_design_matrix
tributary_design_matrices_array[,,30] <- tucannon_river_design_matrix
tributary_design_matrices_array[,,16] <- umatilla_river_design_matrix
tributary_design_matrices_array[,,20] <- walla_walla_river_design_matrix
tributary_design_matrices_array[,,22] <- wenatchee_river_design_matrix
tributary_design_matrices_array[,,18] <- yakima_river_design_matrix

# Create a vector for the order of these tributaries
trib_det_eff_order <- c("Asotin_Creek", 
                        "Deschutes_River", 
                        "Entiat_River", 
                        "Fifteenmile_Creek", 
                        # "Hood_River",
                        "Imnaha_River",
                        "John_Day_River", 
                        "Methow_River", 
                        "Okanogan_River", 
                        "Tucannon_River", 
                        "Umatilla_River",
                        "Walla_Walla_River",
                        "Wenatchee_River", 
                        "Yakima_River")


##### Reformat origin information #####
# tag_code_metadata <- read.csv(here::here("Data", "covariate_data", "tag_code_metadata.csv"))
tag_code_metadata <- read.csv("Data/covariate_data/tag_code_metadata.csv")
# keep only the fish that are in the dataset
# This keeps only hatchery origin fish
tag_code_metadata <- subset(tag_code_metadata, tag_code %in% states_complete$tag_code)

# Convert origins into numbers
# All relative to Yakima, since it's the last one alphabetically
origin_numeric <- data.frame(natal_origin = c("Asotin_Creek", 
                                              "Clearwater_River",
                                              "Deschutes_River", 
                                              "Entiat_River", 
                                              "Fifteenmile_Creek", 
                                              "Grande_Ronde_River", 
                                              # "Hood_River",
                                              "Imnaha_River",
                                              "John_Day_River", 
                                              "Methow_River", 
                                              "Okanogan_River", 
                                              "Salmon_River", 
                                              "Tucannon_River", 
                                              "Umatilla_River",
                                              "Walla_Walla_River",
                                              "Wenatchee_River", 
                                              "Yakima_River"),
                             natal_origin_numeric = seq(1,16,1))

# natal_origin_table <- read.csv(here::here("Data", "covariate_data", "natal_origin_table.csv"))
natal_origin_table <- read.csv("Data/covariate_data/natal_origin_table.csv")
tag_code_metadata %>% 
  left_join(., natal_origin_table, by = "release_site_name") %>% 
  left_join(., origin_numeric, by = "natal_origin") -> tag_code_metadata

# At this point, we need to recreate tag_code_metadata but with the tag_code_2 field
states_complete %>% 
  distinct(tag_code_2, .keep_all = TRUE) %>% 
  dplyr::select(tag_code, tag_code_2) -> tag_codes_2

# reformat this into origin_hatchery info
tag_codes_2 %>%
  left_join(dplyr::select(tag_code_metadata, tag_code, natal_origin_numeric), by = "tag_code") %>% 
  dplyr::rename(natal_origin = natal_origin_numeric) %>% 
  dplyr::select(-tag_code) -> origin_hatchery_actual
  

fish_sim_cat_data_actual <- origin_hatchery_actual
  
  
# Store quantities for loop
# Store the total number of individuals
n.ind <- dim(state_data)[3]
  
  # Store the number of observations per individual
  # -1 because the last state is loss, which isn't actually an observation
  n.obs <- vector(length = n.ind)
  for (i in 1:n.ind){
    n.obs[i] <- sum(state_data[,,i]) - 1
  }
  
  
  
  # Get the state that each fish was in at each n.obs
  
  # First initialize an empty list
  states_list <- list()
  
  for (i in 1:n.ind){
    vec <- vector(length = (n.obs[i]-1))
    
    # Store in list
    states_list[[i]] <- vec
  }
  
  # Store with the states
  for (i in 1:n.ind){
    for (j in 1:(n.obs[i])){
      # states_list[[i]][j] <- rownames(as.data.frame(which(state_data[[i]][,j] == 1)))
      states_list[[i]][j] <- which(state_data[,j,i] == 1) # Get the index of the site instead of the name
    }
  }
  
  # Turn into matrix for stan
  states_mat <- matrix(nrow = n.ind, ncol = max(n.obs))
  for (i in 1:n.ind){
    states_mat[i,1:(n.obs[i])] <- states_list[[i]]
  }
  
  
  # Create the design matrix for origin + run year
  # new for detection efficiency: add a run year field to allow for glm calculation
  # of detection efficiency, as well as to note when certain movements are not possible
  # cat_X_mat_actual <- matrix(0, nrow = n.ind, ncol = 4)
  
  
  # We have 3 origins - so this is going to be three columns
  cat_X_mat_actual <- matrix(0, nrow = n.ind, ncol = 3)
  # Start it so that they're all 0s
  
  # This is for origin
  for (i in 1:n.ind){
    # Natal origin
    if (fish_sim_cat_data_actual$natal_origin[i] == subset(origin_numeric, natal_origin == "Wenatchee_River")$natal_origin_numeric){ # Wenatchee River
      cat_X_mat_actual[i,1] <- 1
    }
    else if (fish_sim_cat_data_actual$natal_origin[i] == subset(origin_numeric, natal_origin == "Okanogan_River")$natal_origin_numeric){ # Okanogan River
      cat_X_mat_actual[i,2] <- 1
    }
    else { # for Methow River
      cat_X_mat_actual[i,3] <- 1
    }
  }
  
  
  # Create a design matrix for temperature - each origin gets an effect; outside of the DPS boundaries, all share an effect
  temp_X_mat_actual <- matrix(0, nrow = n.ind, ncol = 4)
  # The first column everyone gets a 1 (this is temperature effect for all fish, outside of DPS boundaries)
  temp_X_mat_actual[,1] <- 1
  
  for (i in 1:n.ind){
    if (fish_sim_cat_data_actual$natal_origin[i] == subset(origin_numeric, natal_origin == "Wenatchee_River")$natal_origin_numeric){ # Wenatchee River
      temp_X_mat_actual[i,2] <- 1
    }
    else if (fish_sim_cat_data_actual$natal_origin[i] == subset(origin_numeric, natal_origin == "Okanogan_River")$natal_origin_numeric){ # Okanogan River
      temp_X_mat_actual[i,3] <- 1
    }
    else if (fish_sim_cat_data_actual$natal_origin[i] == subset(origin_numeric, natal_origin == "Methow_River")$natal_origin_numeric){ # Methow River
      temp_X_mat_actual[i,4] <- 1
    }
  }
  
  # Create a design matrix for year effects - each origin gets an effect, and these are only inside DPS boundaries
  year_X_mat_actual <- matrix(0, nrow = n.ind, ncol = 3)
  
  for (i in 1:n.ind){
    if (fish_sim_cat_data_actual$natal_origin[i] == subset(origin_numeric, natal_origin == "Wenatchee_River")$natal_origin_numeric){ # Wenatchee River
      year_X_mat_actual[i,1] <- 1
    }
    else if (fish_sim_cat_data_actual$natal_origin[i] == subset(origin_numeric, natal_origin == "Okanogan_River")$natal_origin_numeric){ # Okanogan River
      year_X_mat_actual[i,2] <- 1
    }
    else if (fish_sim_cat_data_actual$natal_origin[i] == subset(origin_numeric, natal_origin == "Methow_River")$natal_origin_numeric){ # Methow River
      year_X_mat_actual[i,3] <- 1
    }
  }
  
  
  # "Tucannon_River", "Asotin_Creek", "Clearwater_River", "Salmon_River", "Grande_Ronde_River", "Imnaha_River"
  # One tweak: We have to replace all NAs in our input data for stan to accept it
  states_mat[is.na(states_mat)] <- -999
  
  
  # We are going to transform our detection histories to instead be vectors containing
  # the number of each site that was visited, instead of a matrix of sites with 0s for
  # not visited and 1s for visited
  
  # First, create an empty matrix to store the newly formatted detection histories
  state_data_2 <- matrix(0, nrow = dim(state_data)[3], ncol = dim(state_data)[2])
  
  # Now fill it in
  for (i in 1:dim(state_data)[3]){
    det_hist <- state_data[,,i]
    
    # Count the number of site visits
    nsite_visits <- sum(det_hist)
    
    for (j in 1:nsite_visits){
      state_data_2[i,j] <- which(det_hist[,j] == 1, arr.ind = TRUE)
    }
    
  }
    
  
  # Get vector of run years of fish
  # need the -1 correction for indices because we're ignoring the 04/05 run year
  fish_run_years <- vector(length = n.ind)
  
  states_complete %>% 
    distinct(tag_code_2, .keep_all = TRUE) %>% 
    mutate(run_year = subset(run_year_df, run_year_start <= date_time & run_year_end >= date_time)$run_year) %>% 
    dplyr::select(tag_code_2, run_year) -> tag_code_2_run_years
  
  for (i in 1:n.ind){
    fish_run_years[i] <- which(run_year_df$run_year == tag_code_2_run_years$run_year[i])-1
  }
  
  
  ##### Create the run year/detection efficiency array #####
  # we need an array that we can index - current (from, rows) x k (to, columns) x run years (slices) - and the values
  # indicate which of the two matrices to use
  # 0 means no detection efficiency correction; 1 means use detection efficiency matrix to pull parameter value from
  # 20 run years: 04/05 through 23/24 (this is the same ad the detection efficiency GLM indexing)
  
  # Also note - we only need to see transitions from mainstem to river mouth
  run_year_DE_array <- array(0, dim = c(nstates, nstates, 20))
  
  # figure out what the indices are of the various tributaries with detection efficiency estimation capability
  print(paste0(seq(1,41,1), " - ", model_states))
  
  # Asotin Creek
  # From Upstream LGR (9) to mouth (32) or upstream (33)
  run_year_DE_array[9,32,8:20] <- 1
  
  # Deschutes River
  # From mainstem BON to MCN (2) to mouth (10) or upstream (11)
  run_year_DE_array[2,10,10:15] <- 1
  
  # Entiat River
  # From RRE to WEL (6) to mouth (24) or upstream (25)
  # CHANGE 2023-07-18 - 07/08 is now an NDE year, not DE
  # previous line: run_year_DE_array[6,26,4:18] <- 1
  run_year_DE_array[6,24,5:20] <- 1
  
  # Fifteenmile Creek
  # From BON to MCN (2) to mouth (16) or upstream (17)
  # Note here that in years 16:18, we no longer have upstream sites, but still have a river mouth site
  # CHANGE 2023-07-18 - 11/12 is now an NDE year, not DE
  # previous line: run_year_DE_array[2,16,8:18] <- 1
  run_year_DE_array[2,14,9:20] <- 1
  
  # Hood River
  # From BON to MCN (2) to mouth (14) or upstream (15)
  # CHANGE 2023-07-18 - 12/13 is now an NDE year, not DE
  # previous line: run_year_DE_array[2,14,9:18] <- 1
  # run_year_DE_array[2,14,10:18] <- 1
  
  # Imnaha River
  # From Upstream LGR (9) to mouth (37) or upstream (38)
  # CHANGE 2023-07-18 - 10/11 is now an NDE year, not DE
  # previous line: run_year_DE_array[9,39,7:18] <- 1
  run_year_DE_array[9,37,8:20] <- 1
  
  # John Day River
  # From BON to MCN (2) to mouth (12) or upstream (13)
  run_year_DE_array[2,12,9:20] <- 1
  
  # Methow River
  # From Upstream WEL (7) to mouth (28) or upstream (29)
  run_year_DE_array[7,28,6:20] <- 1
  
  # Okanogan River
  # From Upstream WEL (7) to mouth (26) or upstream (27)
  run_year_DE_array[7,26,10:20] <- 1
  
  # Tucannon River
  # From ICH to LGR (8) to mouth (30) or upstream (31)
  # CHANGE 2023-07-18 - 10/11 is now an NDE year, not DE
  # previous line: run_year_DE_array[8,32,7:18] <- 1
  run_year_DE_array[8,30,8:20] <- 1
  
  # Umatilla River
  # From BON to MCN (2) to mouth (16) or upstream (17)
  # CHANGE 2023-07-18 - 06/07 is now an NDE year, not DE
  # previous line: run_year_DE_array[2,18,3:18] <- 1
  run_year_DE_array[2,16,4:20] <- 1
  
  # Walla Walla River
  # From MCN to ICH or PRA (3) to mouth (20) or upstream (21)
  run_year_DE_array[3,20,2:20] <- 1
  
  # Wenatchee River
  # From RIS to RRE (5) to mouth (22) or upstream (23)
  # CHANGE 2023-07-18 - 10/11 is now an NDE year, not DE
  # previous line: run_year_DE_array[5,24,7:18] <- 1
  run_year_DE_array[5,22,8:20] <- 1
  
  # Yakima River
  # From MCN to ICH or PRA (3) to mouth (18) or upstream (19)
  # run_year_DE_array[3,20,1:18] <- 1
  # so for now - remove run year 04/05, otherwise it'll all get confused because there's no discharge data for that year
  # also the array doesn't start until late October 2004, so not appropriate to have that be DE anyway
  run_year_DE_array[3,18,2:20] <- 1
  
  
  
  # Load the data from the detection efficiency model to use as priors in this model
  # det_eff_param_posteriors <- as.matrix(read.csv(here::here("Stan", "detection_efficiency", "det_eff_param_posteriors.csv"), row.names = 1))
  det_eff_param_posteriors <- as.matrix(read.csv("Stan/detection_efficiency/det_eff_param_posteriors.csv", row.names = 1))
  
  # get a vector of run years that transitions occurred within - need to do this at the end, once you've made all changes to states_complete
  
  
  
  # Edit 2023-07-21: You also need to drop the ones that have "loss" assigned because they were trapped - 
  # this will make the indexing line up
  states_complete %>% 
    filter(state != "loss") -> states_complete_noloss
  
  
  transition_run_years <- states_complete_noloss$run_year_index
  
  ntransitions = length(transition_run_years)
  
  print("nyears = ", max(transition_run_years))
  
  ##### Convert the dates into a numeric index - with edit from above #####
  
  # # create an empty matrix
  transition_date_matrix <- as.data.frame(matrix(data = NA, nrow = nfish, ncol = max_visits))
  # populate it with the dates
  for (i in 1:nrow(states_complete_noloss)){
    transition_date_matrix[states_complete_noloss[i,"tag_code_number", drop = TRUE], states_complete_noloss[i, "order", drop = TRUE]] <- as.character(date(states_complete_noloss[i, "date_time", drop = TRUE]))
  }
  # deal with the weird numeric issue
  
  
  # 
  # 
  # # We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  date_numeric <- function(x, na.rm = FALSE) as.numeric(ymd(date(x)) - ymd("2005-05-31"))
  transition_date_matrix %>%
    as_tibble() %>%
    mutate_all(date_numeric) %>% 
    as.matrix() -> transition_date_numeric
  
  # Get transition day of year as well
  transition_date_matrix %>%
    as_tibble() %>%
    mutate_all(yday) %>% 
    as.matrix() -> day_of_year_numeric
  
  # Convert day of year to season based on cutoff of June 1
  yday_to_season <- function(x, na.rm = FALSE) ifelse(x < 152, 0, 1)
  
  # 1 = winter/spring; 2 = summer/fall
  transition_date_matrix %>%
    as_tibble() %>%
    mutate_all(yday) %>% 
    mutate_all(yday_to_season) %>% 
    as.matrix() -> transition_seasons_numeric
  
  # replace all NAs with -999s in the above three matrices - they won't be used for anything (and a -999 in this case 
  # is nonsensical), but this keeps Stan from freaking out
  transition_date_numeric[is.na(transition_date_numeric)] <- -999
  day_of_year_numeric[is.na(day_of_year_numeric)] <- -999
  transition_seasons_numeric[is.na(transition_seasons_numeric)] <- -999
  
  # now, convert this to a vector
  as.vector(transition_seasons_numeric)[!(as.vector(transition_seasons_numeric) == -999)] -> transition_seasons_vector
  
  
  
  
  # load in the temperature data
  # temp_data <- read.csv(here::here("Data", "covariate_data", "model_inputs","window_temps_for_stan.csv"), row.names = 1)
  temp_data <- read.csv("Data/covariate_data/model_inputs/window_temps_for_stan.csv", row.names = 1)
  
  # load the spill data
  # spill windows
  # spill_data <- read.csv(here::here("Data", "covariate_data", "model_inputs","window_spill_for_stan.csv"), row.names = 1)
  spill_data <- read.csv("Data/covariate_data/model_inputs/window_spill_for_stan.csv", row.names = 1)
  
  # spill days by winter/spring month
  # jan_spill_days <- read.csv(here::here("Data", "covariate_data", "model_inputs","january_spill_df_for_stan.csv"), row.names = 1)
  # feb_spill_days <- read.csv(here::here("Data", "covariate_data", "model_inputs","february_spill_df_for_stan.csv"), row.names = 1)
  # mar_spill_days <- read.csv(here::here("Data", "covariate_data", "model_inputs","march_spill_df_for_stan.csv"), row.names = 1)
  jan_spill_days <- read.csv("Data/covariate_data/model_inputs/january_spill_df_for_stan.csv", row.names = 1)
  feb_spill_days <- read.csv("Data/covariate_data/model_inputs/february_spill_df_for_stan.csv", row.names = 1)
  mar_spill_days <- read.csv("Data/covariate_data/model_inputs/march_spill_df_for_stan.csv", row.names = 1)
  
  # combine these into one df, for total number of winter spill days in Jan/Feb/Mar
  jan_spill_days + feb_spill_days + mar_spill_days -> winter_spill_days
  
  
  # create a df with combinations of natal origin and state that would qualify
  # as post-overshoot fallback
  # this is structured as a list, where the names of the objects in the list
  # are the states, and then the objects in the list are vectors of the natal
  # origins for which that state qualifies as a post-overshoot state
  
  # Note that this needs to change for the DPS, but not for hatchery vs. wild
  post_overshoot_combos_UC <- list(c("Wenatchee_River"),
                                c("Wenatchee_River","Entiat_River"),
                                c("Wenatchee_River","Entiat_River", "Okanogan River", "Methow River"),
                                c("Wenatchee_River","Entiat_River", "Okanogan River", "Methow River"))
  
  names(post_overshoot_combos_UC) <- c("mainstem, RRE to WEL", "mainstem, upstream of WEL",
                                       "mainstem, ICH to LGR", "mainstem, upstream of LGR")
  
  # a function that takes states_complete and takes post_overshoot combos and a month
  # and returns a vector of movements that qualify as potentially post-overshoot fallback
  
  # now, the function
  # month_numeric can be a vector of months, and this function will check if
  # the fish was in a state in any of those months
  post_overshoot_vector <- function(states_complete, post_overshoot_combos,
                                    month_numeric){
    
    # matrix to index whether or not fish could experienced winter spill conditions
    # one each per month
    # # create an empty matrix
    month_spill_vector <- vector(length = nrow(states_complete))
    
    
    # this happens in two steps:
    # Step 1: check if they're the right state + origin combo to be considered 
    # post-overshoot fallback
    # Step 2: check if the timing is right for them to have experienced 
    # these spill conditions
    
    for (i in 1:nrow(states_complete)){
      
      # Step 1: check if they're the right state + origin combo to be considered 
      # post-overshoot fallback
      if (states_complete[i, "state", drop = TRUE] %in% names(post_overshoot_combos) &
          states_complete[i, "natal_origin", drop = TRUE] %in% post_overshoot_combos[[states_complete[i, "state", drop = TRUE]]]){
        
        # Step 2: check if the timing is right for them to have experienced 
        # these spill conditions
        
        # special case for the last observation in the whole dataset - it has to be 
        # the last time that fish was observed, but we can't do the i + 1 indexing
        if (i == nrow(states_complete)){
          month_spill_vector[i] <- 1
        }

          month_1 <- month(states_complete[i, "window_date_time_1", drop = TRUE])
          month_2 <- month(states_complete[i, "window_date_time_2", drop = TRUE])
          
          # if it is the last observation of a fish, then we don't subset by time
          # because it may have experienced those spill conditions
          if(is.na(month_2)){
            month_spill_vector[i] <- 1
            # if it's not the last observation of a fish, create a sequence of months
            # from the month entering the state to the month exiting the state;
            # if our spill month is in that vector, keep it
          } else if (month_1 > month_2){
            if(any(month_numeric %in% c(seq(month_1, 12), seq(1, month_2)))){
              month_spill_vector[i] <- 1
            } else{
              month_spill_vector[i] <- 0
            }
          }  else {
            if(any(month_numeric %in% seq(month_1, month_2))){
              month_spill_vector[i] <- 1
            } else{
              month_spill_vector[i] <- 0
            }
          }
          

        
        # if it's not in the right state + origin combo, ignore
      } else {
        month_spill_vector[i] <- 0
      }
      
    }
    
    return(month_spill_vector)
    
  }
  
  
  # run it for jan, feb, march - for winter spill days
  
  winter_post_overshoot_vector <- post_overshoot_vector(states_complete = states_complete_noloss, 
                                                     post_overshoot_combos = post_overshoot_combos_UC,
                                    month_numeric = c(1,2,3))

  
  
  
  #### Determine which states are never visited for fish in this dataset
  
  
  
  model_states %>% 
    as.data.frame() %>% 
    dplyr::rename(state = ".") %>% 
    mutate(index = row_number()) -> model_states_df
  
  states_complete %>% 
    left_join(., model_states_df, by = "state") %>% 
    mutate(transition = ifelse(tag_code_2 == lag(tag_code_2), paste0(lag(state), " - ", state), NA)) %>% 
    mutate(transition_numeric = ifelse(tag_code_2 == lag(tag_code_2), paste0(lag(index), " - ", index), NA)) -> state_transitions
  
  # add two more columns: from and to, as numeric
  state_transitions %>% 
    mutate(from = as.numeric(gsub(" -.*", "", transition_numeric))) %>% 
    mutate(to = as.numeric(gsub(".* - ", "", transition_numeric))) -> state_transitions
  
  state_transitions %>% 
    ungroup() %>% 
    count(from) -> from_state_transition_counts
  
  state_transitions %>% 
    ungroup() %>% 
    count(from, to) -> transition_counts
  
  # some data checks
  table(state_transitions$transition)
  table(state_transitions$transition_numeric)
  setdiff(unique(state_transitions$state), model_states)
  
  # EXPORT ALL OF THE STATES THAT WERE NEVER VISITED AS A VECTOR
  setdiff(model_states, unique(state_transitions$state)) -> states_not_visited
  
  ##### Edit the model setup - for this DPS, cut out state transitions for which we have no data, then update parameters based on transition_matrix #####
  # Remove the states that no fish ever enter
  # This is done by setting both from (rows) and to (columns) to zero for that state
  
  # When we remove some of the upstream transitions, we're left with a lot of upstream states
  # that no fish have ever been in.
  # There are also never any upstream to RM or vice versa transitions, because we removed those. I set
  # those all to zeros in the above chunk.
  
  
  for (i in 1:length(states_not_visited)){
    transition_matrix[states_not_visited[i], ] <- 0
    transition_matrix[, states_not_visited[i]] <- 0
  }
  
  ## Get some parameters for model based on updated transition matrix
  # Use the transition matrix to calculate the possible movements
  possible_movements <- rowSums(transition_matrix)
  
  # Get the indices that are 1s (except loss, since that's 1 minus the others)
  movements <- which(transition_matrix[,1:(nstates-1)] == 1, arr.ind = TRUE)
  nmovements <- dim(movements)[1]
  
  # Now get all of the movements which are fixed to zero
  not_movements <- which(transition_matrix[,1:(nstates-1)] == 0, arr.ind = TRUE)
  n_notmovements <- dim(not_movements)[1]
  
  #### Create parameter_indices_matrix, to create a key that allows us to reduce all containers by one dimension ####
  # Create parameter_indices_matrix - the key to match a movement to its spot in the parameter vector
  movements2 <- movements[order(movements[,"row"]),]
  # The data that we are putting in the parameter_indices_matrix is the number of possible movements plus one - this will allow us to index to the last spot in the parameters vector
  parameter_indices_matrix <- matrix(data = nrow(movements2)+1, nrow = 41, ncol = 41)
  
  for(i in 1:nrow(movements2)){
    parameter_indices_matrix[movements2[i, "row"], movements2[i, "col"]] <- i
  }
  

  #### SECTION 3: Write out the model (to be copied to stan model) ####
  
  
  #### Declare parameters to be used in the model ####
  
  
  # Paste movements to name parameters
  movements %>% 
    as.data.frame() %>% 
    mutate(b0_matrix_name = paste0("b0_matrix_", row, "_", col)) -> b0_matrix_names
  
  movements %>% 
    as.data.frame() -> movements_df
  
  
  # Now, note which movements are movements from mainstem into tributaries (these will need two versions of parameters)

  # We also need to note what the indices are of river mouth vs. upstream tributary sites. In years where
  # we can calculate detection efficiency at the river mouth, we need to ignore movements into to upstream sites
  # In years where we can't, we don't ignore them, but the same parameters govern movements into both upstream
  # and river mouth sites.
  
  # We also don't care about any movements from upstream to river mouth states and vice versa

  
  b0_matrix_names %>% 
    mutate(mainstem_upstream_movement = ifelse(row %in% mainstem_indices & col %in% upstream_indices, 1, 0)) -> b0_matrix_names
  
  b0_matrix_names %>% 
    mutate(mainstem_river_mouth_movement = ifelse(row %in% mainstem_indices & col %in% river_mouth_indices, 1, 0)) -> b0_matrix_names
  
  b0_matrix_names %>% 
    mutate(within_trib_movement = ifelse(row %in% upstream_indices & col %in% river_mouth_indices |
                                           col %in% upstream_indices & row %in% river_mouth_indices, 1, 0)) -> b0_matrix_names
  
  b0_matrix_names %>% 
    mutate(upstream_mainstem_movement = ifelse(row %in% upstream_indices & col %in% mainstem_indices, 1, 0)) -> b0_matrix_names
  
  b0_matrix_names %>% 
    mutate(river_mouth_mainstem_movement = ifelse(row %in% river_mouth_indices & col %in% mainstem_indices, 1, 0)) -> b0_matrix_names
  
  
  movements_df %>% 
    mutate(mainstem_upstream_movement = ifelse(row %in% mainstem_indices & col %in% upstream_indices, 1, 0)) -> movements_df
  
  movements_df %>% 
    mutate(mainstem_river_mouth_movement = ifelse(row %in% mainstem_indices & col %in% river_mouth_indices, 1, 0)) -> movements_df
  
  movements_df %>% 
    mutate(within_trib_movement = ifelse(row %in% upstream_indices & col %in% river_mouth_indices |
                                           col %in% upstream_indices & row %in% river_mouth_indices, 1, 0)) -> movements_df
  
  movements_df %>% 
    mutate(upstream_mainstem_movement = ifelse(row %in% upstream_indices & col %in% mainstem_indices, 1, 0)) -> movements_df
  
  movements_df %>% 
    mutate(river_mouth_mainstem_movement = ifelse(row %in% river_mouth_indices & col %in% mainstem_indices, 1, 0)) -> movements_df
  
  upper_columbia_movements <- subset(movements_df, row %in% c(4,5,6,7,seq(22,29,1),40) |
                                       col %in% c(4,5,6,7,seq(22,29,1),40))
  
  non_upper_columbia_movements  <- subset(movements_df, !(row %in% c(4,5,6,7,seq(22,29,1),40)) &
                                            !(col %in% c(4,5,6,7,seq(24,29,1),40)))
  
  # I think these both need to be order by from (rows), to (columns)
  upper_columbia_movements %>% 
    arrange(row, col) -> upper_columbia_movements
  
  # make a version of this for tempxorigin 
  upper_columbia_movements %>% 
    filter(!(col == (row - 1))) %>% 
    filter(!(row == 8 & col == 3)) %>% 
    # remove all movements from tributaries back to the mainstem
    filter(!(col < 10 & row >= 10)) %>% 
    mutate(row_col = paste0(row, "_", col)) -> upper_columbia_tempxorigin
  
  
  b0_matrix_names %>% 
    arrange(row,col) -> b0_matrix_names
  
  # create parameter names for bspillwindow
  # these are only for downstream movements (but all of them)
  movements_df %>% 
    filter(col == (row - 1) | (row == 8 & col == 3)) %>% 
    mutate(row_col = paste0(row, "_", col)) -> fallback_movements
  
  # create parameter names for spill by month parameters
  # post-overshoot fallback movements in the upper columbia DPS: 6 -> 5, 7 -> 6, 8 -> 3, 9 -> 8
  
  fallback_movements %>% 
    filter(row %in% c(6, 7, 8, 9)) -> post_overshoot_fallback_movements
  
  # create parameter name dfs for each of these
  b0_matrix_names %>% 
    mutate(bspillwindow_matrix_name = sub("b0_matrix", "bspillwindow_matrix", b0_matrix_name)) %>% 
    dplyr::select(-b0_matrix_name) %>% 
    mutate(row_col = paste0(row, "_", col)) %>% 
    filter(row_col %in% fallback_movements$row_col)-> bspillwindow_matrix_names
  
  b0_matrix_names %>% 
    mutate(bjanspill_matrix_name = sub("b0_matrix", "bjanspill_matrix", b0_matrix_name)) %>% 
    mutate(bfebspill_matrix_name = sub("b0_matrix", "bfebspill_matrix", b0_matrix_name)) %>% 
    mutate(bmarspill_matrix_name = sub("b0_matrix", "bmarspill_matrix", b0_matrix_name)) %>% 
    mutate(bwinterspill_matrix_name = sub("b0_matrix", "bwinterspill_matrix", b0_matrix_name)) %>% 
    dplyr::select(-b0_matrix_name) %>% 
    mutate(row_col = paste0(row, "_", col)) %>% 
    filter(row_col %in% post_overshoot_fallback_movements$row_col)-> bspill_month_matrix_names
  
  ### 2024-03-07 modification ###
  # change b0 to not have any intercept terms for movements that also have origin terms
  # basically, just remove all upper_columbia_movements from the b0_matrix_names
  
  upper_columbia_movements %>% 
    mutate(movement_abbrev = paste0(row, "_", col)) -> upper_columbia_movements
  
  b0_matrix_names %>% 
    mutate(movement_abbrev = paste0(row, "_", col)) %>% 
    filter(!(movement_abbrev %in% upper_columbia_movements$movement_abbrev)) -> b0_matrix_names
  
  
  # print declaration of b0_matrix parameters - distinguish which need two versions based on mainstem_trib_movement column
  # Note that this block doesn't change between models for different ESUs
  for (i in 1:(nrow(b0_matrix_names))){
    # Paste the numerator
    
    # If it is a movement from the mainstem into the river mouth state, then it gets two versions of the parameter.
    # If it's not, then it just gets one.
    # If its a movement from the upstream to the river mouth state or vice versa, it doesn't get a parameter at all - 
    # we don't care about these movements.
    
    # Movements from mainstem to river mouth: print two versions
    if (b0_matrix_names$mainstem_river_mouth_movement[i] == 1){
      cat("real ", b0_matrix_names$b0_matrix_name[i], "_DE", ";", "\n", sep = "")
      cat("real ", b0_matrix_names$b0_matrix_name[i], "_NDE", ";", "\n", sep = "")
    }
    
    # If it's a within tributary movement - we don't want it. These have been
    # removed from the data (we don't care about them!)
    else if (b0_matrix_names$within_trib_movement[i] == 1){
      # do nothing!
    }
    # If it's a movement from mainstem to upstream - we don't want a parameter for it,
    # because it'll share the same parameter (see transformed parameters block) 
    # as the movement from mainstem to river mouth
    else if (b0_matrix_names$mainstem_upstream_movement[i] == 1){
      # do nothing!
    }
    
    # If we move from the river mouth state to the mainstem, we want a parameter - this is captured in the below
    # BUT - if it's a movement from the upstream state to the mainstem, we don't want a parameter.
    # these transitions will get the same parameters as the river mouth to mainstem
    else if (b0_matrix_names$upstream_mainstem_movement[i] == 1){
      # do nothing!
    } else {
      # Finally - if it's any other movement, just treat it like normal!
      # If it's not, print just the one version
      cat("real ", b0_matrix_names$b0_matrix_name[i],";", "\n", sep = "")
    }
    
  }
  
  # print declaration of window_spill parameters
  # note that because none of these are movements into tributaries, we never need to print DE vs NDE versions
  for (i in 1:(nrow(bspillwindow_matrix_names))){
    cat("real ", bspillwindow_matrix_names$bspillwindow_matrix_name[i],";", "\n", sep = "")
  }
  
  # print declaration of monthly day of spill parameters
  # note that because none of these are movements into tributaries, we never need to print DE vs NDE versions
  for (i in 1:(nrow(bspill_month_matrix_names))){
    cat("real ", bspill_month_matrix_names$bwinterspill_matrix_name[i],";", "\n", sep = "")
  }
  
  
  
  # Print declaration of btemp parameters
  # Paste movements to temp parameters
  movements %>% 
    as.data.frame() %>% 
    # Remove all downstream movements
    filter(!(col == (row - 1))) %>% 
    filter(!(row == 8 & col == 3)) %>% 
    # remove all movements from tributaries back to the mainstem
    filter(!(col < 10 & row >= 10)) %>% 
    mutate(btemp0_matrix_name = paste0("btemp0_matrix_", row, "_", col)) %>% 
    mutate(btemp1_matrix_name = paste0("btemp1_matrix_", row, "_", col)) -> btemp_matrix_names
  
  # add info on what type of movement it is
  btemp_matrix_names %>% 
    mutate(mainstem_upstream_movement = ifelse(row %in% mainstem_indices & col %in% upstream_indices, 1, 0)) -> btemp_matrix_names
  
  btemp_matrix_names %>% 
    mutate(mainstem_river_mouth_movement = ifelse(row %in% mainstem_indices & col %in% river_mouth_indices, 1, 0)) -> btemp_matrix_names
  
  btemp_matrix_names %>% 
    mutate(within_trib_movement = ifelse(row %in% upstream_indices & col %in% river_mouth_indices |
                                           col %in% upstream_indices & row %in% river_mouth_indices, 1, 0)) -> btemp_matrix_names
  
  btemp_matrix_names %>% 
    mutate(upstream_mainstem_movement = ifelse(row %in% upstream_indices & col %in% mainstem_indices, 1, 0)) -> btemp_matrix_names
  
  btemp_matrix_names %>% 
    mutate(river_mouth_mainstem_movement = ifelse(row %in% river_mouth_indices & col %in% mainstem_indices, 1, 0)) -> btemp_matrix_names
  
  btemp_matrix_names %>% 
    arrange(row,col) -> btemp_matrix_names
  
  # REMOVE ALL THE ONES WITHIN THE ESU - these are only origin-specific
  btemp_matrix_names %>% 
    mutate(row_col = paste0(row, "_", col)) %>% 
    filter(!(row_col %in% upper_columbia_tempxorigin$row_col)) -> btemp_matrix_names
  
  
  
  # Note that this block doesn't change between models for different ESUs
  for (i in 1:(nrow(btemp_matrix_names))){
    # Paste the numerator
    
    # If it is a movement from the mainstem into the river mouth state, then it gets two versions of the parameter.
    # If it's not, then it just gets one.
    # If its a movement from the upstream to the river mouth state or vice versa, it doesn't get a parameter at all - 
    # we don't care about these movements.
    
    # Movements from mainstem to river mouth: print two versions
    if (btemp_matrix_names$mainstem_river_mouth_movement[i] == 1){
      cat("real ", btemp_matrix_names$btemp0_matrix_name[i], "_DE", ";", "\n", sep = "")
      cat("real ", btemp_matrix_names$btemp0_matrix_name[i], "_NDE", ";", "\n", sep = "")
      cat("real ", btemp_matrix_names$btemp1_matrix_name[i], "_DE", ";", "\n", sep = "")
      cat("real ", btemp_matrix_names$btemp1_matrix_name[i], "_NDE", ";", "\n", sep = "")
    }
    
    # If it's a within tributary movement - we don't want it. These have been
    # removed from the data (we don't care about them!)
    else if (btemp_matrix_names$within_trib_movement[i] == 1){
      # do nothing!
    }
    # If it's a movement from mainstem to upstream - we don't want a parameter for it,
    # because it'll share the same parameter (see transformed parameters block) 
    # as the movement from mainstem to river mouth
    else if (btemp_matrix_names$mainstem_upstream_movement[i] == 1){
      # do nothing!
    }
    
    # If we move from the river mouth state to the mainstem, we want a parameter - this is captured in the below
    # BUT - if it's a movement from the upstream state to the mainstem, we don't want a parameter.
    # these transitions will get the same parameters as the river mouth to mainstem
    else if (btemp_matrix_names$upstream_mainstem_movement[i] == 1){
      # do nothing!
    } else {
      # Finally - if it's any other movement, just treat it like normal!
      # If it's not, print just the one version
      cat("real ", btemp_matrix_names$btemp0_matrix_name[i],";", "\n", sep = "")
      cat("real ", btemp_matrix_names$btemp1_matrix_name[i],";", "\n", sep = "")
    }
    
  }
  
  
  # Print declaration of btempxorigin (interaction) parameters
  # Paste movements to temp parameters
  movements %>% 
    as.data.frame() %>% 
    # Remove all downstream movements
    filter(!(col == (row - 1))) %>% 
    filter(!(row == 8 & col == 3)) %>% 
    # remove all movements from tributaries back to the mainstem
    filter(!(col < 10 & row >= 10)) %>% 
    mutate(btemp0xorigin_matrix_name = paste0("btemp0xorigin_matrix_", row, "_", col)) %>% 
    mutate(btemp1xorigin_matrix_name = paste0("btemp1xorigin_matrix_", row, "_", col))-> btempxorigin_matrix_names
  
  # add info on what type of movement it is
  btempxorigin_matrix_names %>% 
    mutate(mainstem_upstream_movement = ifelse(row %in% mainstem_indices & col %in% upstream_indices, 1, 0)) -> btempxorigin_matrix_names
  
  btempxorigin_matrix_names %>% 
    mutate(mainstem_river_mouth_movement = ifelse(row %in% mainstem_indices & col %in% river_mouth_indices, 1, 0)) -> btempxorigin_matrix_names
  
  btempxorigin_matrix_names %>% 
    mutate(within_trib_movement = ifelse(row %in% upstream_indices & col %in% river_mouth_indices |
                                           col %in% upstream_indices & row %in% river_mouth_indices, 1, 0)) -> btempxorigin_matrix_names
  
  btempxorigin_matrix_names %>% 
    mutate(upstream_mainstem_movement = ifelse(row %in% upstream_indices & col %in% mainstem_indices, 1, 0)) -> btempxorigin_matrix_names
  
  btempxorigin_matrix_names %>% 
    mutate(river_mouth_mainstem_movement = ifelse(row %in% river_mouth_indices & col %in% mainstem_indices, 1, 0)) -> btempxorigin_matrix_names
  
  btempxorigin_matrix_names %>% 
    arrange(row,col) -> btempxorigin_matrix_names
  
  # print only for the movements that are within the ESU
  for (i in 1:3){ # since we only have three origins for wild
    for (j in 1:nrow(upper_columbia_tempxorigin)){
      
      # Movements from mainstem to river mouth: print two versions
      if (upper_columbia_tempxorigin$mainstem_river_mouth_movement[j] == 1){
        cat("real ", "btemp0xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], "_DE", ";", "\n", sep = "")
        cat("real ", "btemp0xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], "_NDE", ";", "\n", sep = "")
        cat("real ", "btemp1xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], "_DE", ";", "\n", sep = "")
        cat("real ", "btemp1xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], "_NDE", ";", "\n", sep = "")
      }
      
      # If it's a within tributary movement - we don't want it
      else if (upper_columbia_tempxorigin$within_trib_movement[j] == 1){
        # do nothing!
      }
      # If it's a movement from mainstem to upstream - we also don't want it
      else if (upper_columbia_tempxorigin$mainstem_upstream_movement[j] == 1){
        # do nothing!
      } 
      # If we move from the river mouth state to the mainstem, we want a parameter - this is captured in the below
      # BUT - if it's a movement from the upstream state to the mainstem, we don't want a parameter.
      # these transitions will get the same parameters as the river mouth to mainstem
      else if (upper_columbia_tempxorigin$upstream_mainstem_movement[j] == 1){
        # do nothing!
      } else {
        # Finally - if it's any other movement, just treat it like normal!
        # If it's not, print just the one version
        cat("real ", "btemp0xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], ";", "\n", sep = "")
        cat("real ", "btemp1xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], ";", "\n", sep = "")
      }
      
    }
    cat ("\n")
  }
  
  # In this ESU, there are four origins, so we will have three origin parameters
  # We will only allow an origin effect within the ESU (so after PRA). Before, they will only have an intercept term
  # print declaration of borigin parameters
  for (i in 1:3){ # since we only have three origins for wild, we will only go i in 1:3
    for (j in 1:nrow(upper_columbia_movements)){
      
      # Movements from mainstem to river mouth: print two versions
      if (upper_columbia_movements$mainstem_river_mouth_movement[j] == 1){
        cat("real ", "borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], "_DE", ";", "\n", sep = "")
        cat("real ", "borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], "_NDE", ";", "\n", sep = "")
      }
      
      # If it's a within tributary movement - we don't want it
      else if (upper_columbia_movements$within_trib_movement[j] == 1){
        # do nothing!
      }
      # If it's a movement from mainstem to upstream - we also don't want it
      else if (upper_columbia_movements$mainstem_upstream_movement[j] == 1){
        # do nothing!
      } 
      # If we move from the river mouth state to the mainstem, we want a parameter - this is captured in the below
      # BUT - if it's a movement from the upstream state to the mainstem, we don't want a parameter.
      # these transitions will get the same parameters as the river mouth to mainstem
      else if (upper_columbia_movements$upstream_mainstem_movement[j] == 1){
        # do nothing!
      } else {
        # Finally - if it's any other movement, just treat it like normal!
        # If it's not, print just the one version
        cat("real ", "borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], ";", "\n", sep = "")
      }
      
    }
    cat ("\n")
  }
  
  ## Declare year parameters - similar setup as origin effects, but with an interaction with year
  # However, we are now dropping any movements from tributaries back into the mainstem,
  # because those movements don't have enough data to estimate them
  
  # First, filter out the movements that aren't getting a RE of year
  # based on rownames(transition_matrix), tributaries start at state 10; drop
  # all movements from those states
  upper_columbia_movements %>% 
    filter(!(row>=10)) %>% 
    mutate(row_col = paste0(row, "_", col)) -> upper_columbia_RE_year
  
  
  # Now for each origin, print parameter names - first the raw vector
  for (i in 1:3){ # since we only have three origins for wild
    for (j in 1:nrow(upper_columbia_RE_year)){
      
      # Movements from mainstem to river mouth: print two versions
      if (upper_columbia_RE_year$mainstem_river_mouth_movement[j] == 1){
        cat("vector[nyears] ", "byearxorigin", i, "_raw_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_DE", ";", "\n", sep = "")
        cat("vector[nyears] ", "byearxorigin", i, "_raw_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_NDE", ";", "\n", sep = "")
      }
      
      # If it's a within tributary movement - we don't want it
      else if (upper_columbia_RE_year$within_trib_movement[j] == 1){
        # do nothing!
      }
      # If it's a movement from mainstem to upstream - we also don't want it - because this is a shared parameter with movement from mainstem to river mouth
      else if (upper_columbia_RE_year$mainstem_upstream_movement[j] == 1){
        # do nothing!
      } 
      # There are no year effects for movements from tributaries, but this doesn't matter because it should have been removed earlier
      else if (upper_columbia_RE_year$upstream_mainstem_movement[j] == 1){
        # do nothing!
      } else {
        # Finally - if it's any other movement, just treat it like normal!
        # If it's not, print just the one version
        cat("vector[nyears] ", "byearxorigin", i, "_raw_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], ";", "\n", sep = "")
      }
      
    }
    cat ("\n")
  }
  
  ## Then print the sigma
  for (i in 1:3){ # since we only have three origins for wild
    for (j in 1:nrow(upper_columbia_RE_year)){
      
      
      # If this is a post-overshoot fallback movement, print out a commented out version
      # Origin 1: Wenatchee River: 7-6 and 6-5 are both post-overshoot fallback movements within the DPS
      # Origin 2: Okanogan River: no post-overshoot fallback movements within the DPS
      # Origin 3: Methow River: no post-overshoot fallback movements within the DPS
      
      # this is for the Wenatchee
      if(upper_columbia_RE_year$row[j] == 7 & upper_columbia_RE_year$col[j] == 6 & i == 1 |
         upper_columbia_RE_year$row[j] == 6 & upper_columbia_RE_year$col[j] == 5 & i == 1){
        cat("// real<lower=0> ", "sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_DE", ";", "\n", sep = "")
      }
      
      # Movements from mainstem to river mouth: print two versions
      else if (upper_columbia_RE_year$mainstem_river_mouth_movement[j] == 1){
        cat("real<lower=0> ", "sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_DE", ";", "\n", sep = "")
        cat("real<lower=0> ", "sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_NDE", ";", "\n", sep = "")
      }
      
      # If it's a within tributary movement - we don't want it
      else if (upper_columbia_RE_year$within_trib_movement[j] == 1){
        # do nothing!
      }
      # If it's a movement from mainstem to upstream - we also don't want it
      else if (upper_columbia_RE_year$mainstem_upstream_movement[j] == 1){
        # do nothing!
      } 
      # If we move from the river mouth state to the mainstem, we want a parameter - this is captured in the below
      # BUT - if it's a movement from the upstream state to the mainstem, we don't want a parameter.
      # these transitions will get the same parameters as the river mouth to mainstem
      else if (upper_columbia_RE_year$upstream_mainstem_movement[j] == 1){
        # do nothing!
      } else {
        # Finally - if it's any other movement, just treat it like normal!
        # If it's not, print just the one version
        cat("real<lower=0> ", "sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], ";", "\n", sep = "")
      }
      
    }
    cat ("\n")
  }
  
  
  #### write it out for the transformed parameters section ####
  # storing these occurs after you declare all of the vectors/matrices that are used to index the parameters
  
  # Write out two separate matrices for b0 - an NDE (no detection efficiency correction) and a DE (detection efficiency) matrix
  for (i in 1:(nrow(b0_matrix_names))){
    # Movements from mainstem to river mouth: store the two different versions
    if (b0_matrix_names$mainstem_river_mouth_movement[i] == 1){
      cat("b0_vector_DE[", parameter_indices_matrix[b0_matrix_names$row[i], b0_matrix_names$col[i]], "]", " = ", b0_matrix_names$b0_matrix_name[i], "_DE",";", "\n", sep = "")
      cat("b0_vector_NDE[", parameter_indices_matrix[b0_matrix_names$row[i], b0_matrix_names$col[i]], "]", " = ", b0_matrix_names$b0_matrix_name[i], "_NDE",";", "\n", sep = "")
      
    }
    
    # If it's a within tributary movement - we don't want it
    else if (b0_matrix_names$within_trib_movement[i] == 1){
      # do nothing!
    }
    # If it's a movement from mainstem to upstream:
    # For NDE - give it the same parameter as the one for the mainstem to the river mouth site (which by index, is the site before)
    # For DE - we need to give these all -100000 values, so that in the logit they come out as zeros (since they're not transitions that are allowed)
    else if (b0_matrix_names$mainstem_upstream_movement[i] == 1){
      # cat("b0_matrix_DE[", b0_matrix_names$row[i], ",", b0_matrix_names$col[i], "]", " = ", b0_matrix_names$b0_matrix_name[i-1], "_DE",";", "\n", sep = "")
      cat("b0_vector_DE[", parameter_indices_matrix[b0_matrix_names$row[i], b0_matrix_names$col[i]], "]", " = ", "-100000;", "\n", sep = "")
      cat("b0_vector_NDE[", parameter_indices_matrix[b0_matrix_names$row[i], b0_matrix_names$col[i]], "]", " = ", b0_matrix_names$b0_matrix_name[i-1], "_NDE",";", "\n", sep = "")
    }
    # If it's a movement from the upstream state back to the mainstem, give it the same parameter as movement from the river mouth to the upstream
    # And note - that these are not DE or NDE parameters
    # Because of order, we need -2 here instead of -1, to pick the right from state
    # NO WE DON'T, just use -1
    else if (b0_matrix_names$upstream_mainstem_movement[i] == 1){
      cat("b0_vector_DE[", parameter_indices_matrix[b0_matrix_names$row[i], b0_matrix_names$col[i]], "]", " = ", b0_matrix_names$b0_matrix_name[i-1],";", "\n", sep = "")
      cat("b0_vector_NDE[", parameter_indices_matrix[b0_matrix_names$row[i], b0_matrix_names$col[i]], "]", " = ", b0_matrix_names$b0_matrix_name[i-1],";", "\n", sep = "")
      # If it's not, store the same parameter in both matrices
    } else {
      cat("b0_vector_DE[", parameter_indices_matrix[b0_matrix_names$row[i], b0_matrix_names$col[i]], "]", " = ", b0_matrix_names$b0_matrix_name[i], ";", "\n", sep = "")
      cat("b0_vector_NDE[", parameter_indices_matrix[b0_matrix_names$row[i], b0_matrix_names$col[i]], "]", " = ", b0_matrix_names$b0_matrix_name[i], ";", "\n", sep = "")
    }
    
  } 
  
  ### Edit 2024-03-09 ###
  # We need to have our intercept vector contain zeros for any allowed movements that
  # are now only getting an origin term.
  
  ### 2024-03-09 edits - add a zero for everything in the intercept vector
  for (i in 1:nrow(upper_columbia_movements)){
    cat("b0_vector_DE[", parameter_indices_matrix[upper_columbia_movements$row[i], upper_columbia_movements$col[i]], "]", " = 0;", "\n", sep = "")
    cat("b0_vector_NDE[", parameter_indices_matrix[upper_columbia_movements$row[i], upper_columbia_movements$col[i]], "]", " = 0;", "\n", sep = "")
  }
  
  # We only need one matrix for spill, because none of these are DE corrected
  # print declaration of window_spill parameters
  # note that because none of these are movements into tributaries, we never need to print DE vs NDE versions
  for (i in 1:(nrow(bspillwindow_matrix_names))){
    cat("bspillwindow_vector[", parameter_indices_matrix[bspillwindow_matrix_names$row[i], bspillwindow_matrix_names$col[i]], "]", " = ", bspillwindow_matrix_names$bspillwindow_matrix_name[i], ";", "\n", sep = "")
  }
  
  # print declaration of monthly day of spill parameters
  # note that because none of these are movements into tributaries, we never need to print DE vs NDE versions
  for (i in 1:(nrow(bspill_month_matrix_names))){
    cat("bwinterspill_vector[",  parameter_indices_matrix[bspill_month_matrix_names$row[i], bspill_month_matrix_names$col[i]], "]", " = ", bspill_month_matrix_names$bwinterspill_matrix_name[i], ";", "\n", sep = "")
    
  }
  
  
  ## Year effects
  # 1. store the individual byear_raw vectors into the matrix that contains all of them
  for (i in 1:3){ # since we only have three origins for wild
    for (j in 1:nrow(upper_columbia_RE_year)){
      
      # Movements from mainstem to river mouth: print two versions
      if (upper_columbia_RE_year$mainstem_river_mouth_movement[j] == 1){
        cat("byearxorigin", i, "_raw_parameters_matrix_DE[", parameter_indices_matrix[upper_columbia_RE_year$row[j], upper_columbia_RE_year$col[j]], ", ] = to_array_1d(byearxorigin", i, "_raw_vector_",
            upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_DE);", "\n", sep = "")
        cat("byearxorigin", i, "_raw_parameters_matrix_NDE[", parameter_indices_matrix[upper_columbia_RE_year$row[j], upper_columbia_RE_year$col[j]], ", ] = to_array_1d(byearxorigin", i, "_raw_vector_",
            upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_NDE);", "\n", sep = "")
      }
      
      # If it's a within tributary movement - we don't want it
      else if (upper_columbia_RE_year$within_trib_movement[j] == 1){
        # do nothing!
      }
      # If it's a movement from mainstem to upstream:
      # For NDE - give it the same parameter as the one for the mainstem to the river mouth site (which by index, is the site before)
      # For DE - do nothing - it'll just keep the zero from before
      else if (upper_columbia_RE_year$mainstem_upstream_movement[j] == 1){
        cat("byearxorigin", i, "_raw_parameters_matrix_NDE[", parameter_indices_matrix[upper_columbia_RE_year$row[j], upper_columbia_RE_year$col[j]], ", ] = to_array_1d(byearxorigin", i, "_raw_vector_",
            upper_columbia_RE_year$row[j-1], "_", upper_columbia_RE_year$col[j-1], "_NDE);", "\n", sep = "")
      } 
      # For year effects, there is no RE of year when moving from a tributary to the mainstem. So we will not do anything here.
      # These should not be in the movement df though, so it shouldn't matter.
      else if (upper_columbia_RE_year$upstream_mainstem_movement[j] == 1){
        # do nothing!
      } else {
        # Finally - if it's any other movement, just treat it like normal!
        # If it's not, print just the one version
        cat("byearxorigin", i, "_raw_parameters_matrix_DE[", parameter_indices_matrix[upper_columbia_RE_year$row[j], upper_columbia_RE_year$col[j]], ", ] = to_array_1d(byearxorigin", i, "_raw_vector_",
            upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], ");", "\n", sep = "")
        cat("byearxorigin", i, "_raw_parameters_matrix_NDE[", parameter_indices_matrix[upper_columbia_RE_year$row[j], upper_columbia_RE_year$col[j]], ", ] = to_array_1d(byearxorigin", i, "_raw_vector_",
            upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], ");", "\n", sep = "")
      }
      
    }
    cat ("\n")
  }
  
  # 2. Write out the conversion from raw * sigma into actual
  for (i in 1:3){ # since we only have three origins for wild
    for (j in 1:nrow(upper_columbia_RE_year)){
      
      # Movements from mainstem to river mouth: print two versions
      if (upper_columbia_RE_year$mainstem_river_mouth_movement[j] == 1){
        cat("byearxorigin", i, "_actual_parameters_matrix_DE[", parameter_indices_matrix[upper_columbia_RE_year$row[j], upper_columbia_RE_year$col[j]], ", ] = to_array_1d(byearxorigin", i, "_raw_vector_",
            upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_DE * sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_DE);", "\n", sep = "")
        cat("byearxorigin", i, "_actual_parameters_matrix_NDE[", parameter_indices_matrix[upper_columbia_RE_year$row[j], upper_columbia_RE_year$col[j]], ", ] = to_array_1d(byearxorigin", i, "_raw_vector_",
            upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_NDE * sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_NDE);", "\n", sep = "")
      }
      
      # If it's a within tributary movement - we don't want it
      else if (upper_columbia_RE_year$within_trib_movement[j] == 1){
        # do nothing!
      }
      # If it's a movement from mainstem to upstream:
      # For NDE - give it the same parameter as the one for the mainstem to the river mouth site (which by index, is the site before)
      # For DE - do nothing - it'll just keep the zero from before
      else if (upper_columbia_RE_year$mainstem_upstream_movement[j] == 1){
        cat("byearxorigin", i, "_actual_parameters_matrix_NDE[", parameter_indices_matrix[upper_columbia_RE_year$row[j], upper_columbia_RE_year$col[j]], ", ] = to_array_1d(byearxorigin", i, "_raw_vector_",
            upper_columbia_RE_year$row[j-1], "_", upper_columbia_RE_year$col[j-1], "_NDE * sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j-1], "_", upper_columbia_RE_year$col[j-1], "_NDE);", "\n", sep = "")
      } 
      # For year effects, there is no RE of year when moving from a tributary to the mainstem. So we will not do anything here.
      # These should not be in the movement df though, so it shouldn't matter.
      else if (upper_columbia_RE_year$upstream_mainstem_movement[j] == 1){
        # do nothing!
      } else {
        # Finally - if it's any other movement, just treat it like normal!
        # If it's not, print just the one version
        cat("byearxorigin", i, "_actual_parameters_matrix_DE[", parameter_indices_matrix[upper_columbia_RE_year$row[j], upper_columbia_RE_year$col[j]], ", ] = to_array_1d(byearxorigin", i, "_raw_vector_",
            upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], " * sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], ");", "\n", sep = "")
        cat("byearxorigin", i, "_actual_parameters_matrix_NDE[", parameter_indices_matrix[upper_columbia_RE_year$row[j], upper_columbia_RE_year$col[j]], ", ] = to_array_1d(byearxorigin", i, "_raw_vector_",
            upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], " * sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], ");", "\n", sep = "")
      }
      
    }
    cat ("\n")
  }
  
  # 3. Store the sigma parameters in the vector that contains all of them
  for (i in 1:3){ # since we only have three origins for wild
    for (j in 1:nrow(upper_columbia_RE_year)){
      
      # Movements from mainstem to river mouth: print two versions
      if (upper_columbia_RE_year$mainstem_river_mouth_movement[j] == 1){
        cat("sigma_yearxorigin", i, "_vector_DE[", parameter_indices_matrix[upper_columbia_RE_year$row[j], upper_columbia_RE_year$col[j]], 
            "] = sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_DE;", "\n", sep = "")
        cat("sigma_yearxorigin", i, "_vector_NDE[", parameter_indices_matrix[upper_columbia_RE_year$row[j], upper_columbia_RE_year$col[j]], 
            "] = sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_NDE;", "\n", sep = "")
      }
      
      # If it's a within tributary movement - we don't want it
      else if (upper_columbia_RE_year$within_trib_movement[j] == 1){
        # do nothing!
      }
      # If it's a movement from mainstem to upstream:
      # For NDE - give it the same parameter as the one for the mainstem to the river mouth site (which by index, is the site before)
      # For DE - do nothing - it'll just keep the zero from before
      else if (upper_columbia_RE_year$mainstem_upstream_movement[j] == 1){
        cat("sigma_yearxorigin", i, "_vector_NDE[", parameter_indices_matrix[upper_columbia_RE_year$row[j], upper_columbia_RE_year$col[j]], 
            "] = sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j-1], "_", upper_columbia_RE_year$col[j-1], "_NDE;", "\n", sep = "")
      } 
      # For year effects, there is no RE of year when moving from a tributary to the mainstem. So we will not do anything here.
      # These should not be in the movement df though, so it shouldn't matter.
      else if (upper_columbia_RE_year$upstream_mainstem_movement[j] == 1){
        # do nothing!
      } else {
        # Finally - if it's any other movement, just treat it like normal!
        # If it's not, print just the one version
        cat("sigma_yearxorigin", i, "_vector_DE[", parameter_indices_matrix[upper_columbia_RE_year$row[j], upper_columbia_RE_year$col[j]], 
            "] = sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], ";", "\n", sep = "")
        cat("sigma_yearxorigin", i, "_vector_NDE[", parameter_indices_matrix[upper_columbia_RE_year$row[j], upper_columbia_RE_year$col[j]], 
            "] = sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], ";", "\n", sep = "")
      }
      
    }
    cat ("\n")
  }
  
  
  
  # Write out the temp matrix - DE and NDE
  for (i in 1:(nrow(btemp_matrix_names))){
    # Movements from mainstem to river mouth: store the two different versions
    if (btemp_matrix_names$mainstem_river_mouth_movement[i] == 1){
      cat("btemp0_vector_DE[", parameter_indices_matrix[btemp_matrix_names$row[i], btemp_matrix_names$col[i]], "]", " = ", btemp_matrix_names$btemp0_matrix_name[i], "_DE",";", "\n", sep = "")
      cat("btemp0_vector_NDE[", parameter_indices_matrix[btemp_matrix_names$row[i], btemp_matrix_names$col[i]], "]", " = ", btemp_matrix_names$btemp0_matrix_name[i], "_NDE",";", "\n", sep = "")
      cat("btemp1_vector_DE[", parameter_indices_matrix[btemp_matrix_names$row[i],  btemp_matrix_names$col[i]], "]", " = ", btemp_matrix_names$btemp1_matrix_name[i], "_DE",";", "\n", sep = "")
      cat("btemp1_vector_NDE[", parameter_indices_matrix[btemp_matrix_names$row[i], btemp_matrix_names$col[i]], "]", " = ", btemp_matrix_names$btemp1_matrix_name[i], "_NDE",";", "\n", sep = "")
      
    }
    
    # If it's a within tributary movement - we don't want it
    else if (btemp_matrix_names$within_trib_movement[i] == 1){
      # do nothing!
    }
    # If it's a movement from mainstem to upstream:
    # For NDE - give it the same parameter as the one for the mainstem to the river mouth site (which by index, is the site before)
    # For DE - do nothing, because these transition are not allowed in DE
    else if (btemp_matrix_names$mainstem_upstream_movement[i] == 1){
      cat("btemp0_vector_NDE[", parameter_indices_matrix[btemp_matrix_names$row[i], btemp_matrix_names$col[i]], "]", " = ", btemp_matrix_names$btemp0_matrix_name[i-1], "_NDE",";", "\n", sep = "")
      cat("btemp1_vector_NDE[", parameter_indices_matrix[btemp_matrix_names$row[i], btemp_matrix_names$col[i]], "]", " = ", btemp_matrix_names$btemp1_matrix_name[i-1], "_NDE",";", "\n", sep = "")
    }
    # If it's a movement from the upstream state back to the mainstem, give it the same parameter as movement from the river mouth to the upstream
    # And note - that these are not DE or NDE parameters
    # Because of order, we need -2 here instead of -1, to pick the right from state
    # NO WE DON'T, just use -1
    else if (btemp_matrix_names$upstream_mainstem_movement[i] == 1){
      cat("btemp0_vector_DE[", parameter_indices_matrix[btemp_matrix_names$row[i], btemp_matrix_names$col[i]], "]", " = ", btemp_matrix_names$btemp0_matrix_name[i-1],";", "\n", sep = "")
      cat("btemp0_vector_NDE[", parameter_indices_matrix[btemp_matrix_names$row[i], btemp_matrix_names$col[i]], "]", " = ", btemp_matrix_names$btemp0_matrix_name[i-1],";", "\n", sep = "")
      cat("btemp1_vector_DE[", parameter_indices_matrix[btemp_matrix_names$row[i], btemp_matrix_names$col[i]], "]", " = ", btemp_matrix_names$btemp1_matrix_name[i-1],";", "\n", sep = "")
      cat("btemp1_vector_NDE[", parameter_indices_matrix[btemp_matrix_names$row[i], btemp_matrix_names$col[i]], "]", " = ", btemp_matrix_names$btemp1_matrix_name[i-1],";", "\n", sep = "")
      # If it's not, store the same parameter in both matrices
    } else {
      cat("btemp0_vector_DE[", parameter_indices_matrix[btemp_matrix_names$row[i], btemp_matrix_names$col[i]], "]", " = ", btemp_matrix_names$btemp0_matrix_name[i], ";", "\n", sep = "")
      cat("btemp0_vector_NDE[", parameter_indices_matrix[btemp_matrix_names$row[i], btemp_matrix_names$col[i]], "]", " = ", btemp_matrix_names$btemp0_matrix_name[i], ";", "\n", sep = "")
      cat("btemp1_vector_DE[", parameter_indices_matrix[btemp_matrix_names$row[i], btemp_matrix_names$col[i]], "]", " = ", btemp_matrix_names$btemp1_matrix_name[i], ";", "\n", sep = "")
      cat("btemp1_vector_NDE[", parameter_indices_matrix[btemp_matrix_names$row[i], btemp_matrix_names$col[i]], "]", " = ", btemp_matrix_names$btemp1_matrix_name[i], ";", "\n", sep = "")
    }
    
  } 
  
  
  # Write out the tempxorigin matrix - DE
  for (i in 1:3){
    for (j in 1:nrow(upper_columbia_tempxorigin)){
      # Movements from mainstem to river mouth: store the two different versions
      if (upper_columbia_tempxorigin$mainstem_river_mouth_movement[j] == 1){
        cat("btemp0xorigin", i, "_vector_DE[", parameter_indices_matrix[upper_columbia_tempxorigin$row[j], upper_columbia_tempxorigin$col[j]], "]", " = ", "btemp0xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], "_DE", ";", "\n", sep = "")
        cat("btemp0xorigin", i, "_vector_NDE[", parameter_indices_matrix[upper_columbia_tempxorigin$row[j], upper_columbia_tempxorigin$col[j]], "]", " = ", "btemp0xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], "_NDE", ";", "\n", sep = "")
        cat("btemp1xorigin", i, "_vector_DE[", parameter_indices_matrix[upper_columbia_tempxorigin$row[j], upper_columbia_tempxorigin$col[j]], "]", " = ", "btemp1xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], "_DE", ";", "\n", sep = "")
        cat("btemp1xorigin", i, "_vector_NDE[", parameter_indices_matrix[upper_columbia_tempxorigin$row[j], upper_columbia_tempxorigin$col[j]], "]", " = ", "btemp1xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], "_NDE", ";", "\n", sep = "")
        
      }
      
      # If it's a within tributary movement - we don't want it
      else if (upper_columbia_tempxorigin$within_trib_movement[j] == 1){
        # do nothing!
      }
      # If it's a movement from mainstem to upstream:
      # For NDE - give it the same parameter as the one for the mainstem to the river mouth site (which by index, is the site before)
      # For DE - we don't care about this, because we're removing the upstream states
      else if (upper_columbia_tempxorigin$mainstem_upstream_movement[j] == 1){
        cat("btemp0xorigin", i, "_vector_NDE[", parameter_indices_matrix[upper_columbia_tempxorigin$row[j], upper_columbia_tempxorigin$col[j]], "]", " = ", "btemp0xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j-1], "_", upper_columbia_tempxorigin$col[j-1], "_NDE", ";", "\n", sep = "")
        cat("btemp1xorigin", i, "_vector_NDE[", parameter_indices_matrix[upper_columbia_tempxorigin$row[j], upper_columbia_tempxorigin$col[j]], "]", " = ", "btemp1xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j-1], "_", upper_columbia_tempxorigin$col[j-1], "_NDE", ";", "\n", sep = "")
        
      }
      # If it's a movement from the upstream state back to the mainstem, give it the same parameter as movement from the river mouth to the upstream
      # And note - that these are not DE or NDE parameters
      else if (upper_columbia_tempxorigin$upstream_mainstem_movement[j] == 1){
        cat("btemp0xorigin", i, "_vector_DE[", parameter_indices_matrix[upper_columbia_tempxorigin$row[j], upper_columbia_tempxorigin$col[j]], "]", " = ", "btemp0xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j-1], "_", upper_columbia_tempxorigin$col[j],  ";", "\n", sep = "")
        cat("btemp0xorigin", i, "_vector_NDE[", parameter_indices_matrix[upper_columbia_tempxorigin$row[j], upper_columbia_tempxorigin$col[j]], "]", " = ", "btemp0xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j-1], "_", upper_columbia_tempxorigin$col[j], ";", "\n", sep = "")
        cat("btemp1xorigin", i, "_vector_DE[", parameter_indices_matrix[upper_columbia_tempxorigin$row[j], upper_columbia_tempxorigin$col[j]], "]", " = ", "btemp1xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j-1], "_", upper_columbia_tempxorigin$col[j],  ";", "\n", sep = "")
        cat("btemp1xorigin", i, "_vector_NDE[", parameter_indices_matrix[upper_columbia_tempxorigin$row[j], upper_columbia_tempxorigin$col[j]], "]", " = ", "btemp1xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j-1], "_", upper_columbia_tempxorigin$col[j], ";", "\n", sep = "")
        # If it's not, store the same parameter in both matrices
      } else {
        cat("btemp0xorigin", i, "_vector_DE[", parameter_indices_matrix[upper_columbia_tempxorigin$row[j], upper_columbia_tempxorigin$col[j]], "]", " = ", "btemp0xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], ";", "\n", sep = "")
        cat("btemp0xorigin", i, "_vector_NDE[", parameter_indices_matrix[upper_columbia_tempxorigin$row[j], upper_columbia_tempxorigin$col[j]], "]", " = ", "btemp0xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j],  ";", "\n", sep = "")
        cat("btemp1xorigin", i, "_vector_DE[", parameter_indices_matrix[upper_columbia_tempxorigin$row[j], upper_columbia_tempxorigin$col[j]], "]", " = ", "btemp1xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], ";", "\n", sep = "")
        cat("btemp1xorigin", i, "_vector_NDE[", parameter_indices_matrix[upper_columbia_tempxorigin$row[j], upper_columbia_tempxorigin$col[j]], "]", " = ", "btemp1xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j],  ";", "\n", sep = "")
      }
      
    }
    
    cat ("\n")
  }
  
  
  # There are three wild origins, we will only go for i in 1:3
  for (i in 1:3){
    for (j in 1:nrow(upper_columbia_movements)){
      # Movements from mainstem to river mouth: store the two different versions
      if (upper_columbia_movements$mainstem_river_mouth_movement[j] == 1){
        cat("borigin", i, "_vector_DE[", parameter_indices_matrix[upper_columbia_movements$row[j], upper_columbia_movements$col[j]], "]", " = ", "borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], "_DE", ";", "\n", sep = "")
        cat("borigin", i, "_vector_NDE[", parameter_indices_matrix[upper_columbia_movements$row[j], upper_columbia_movements$col[j]], "]", " = ", "borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], "_NDE", ";", "\n", sep = "")
        
      }
      
      # If it's a within tributary movement - we don't want it
      else if (upper_columbia_movements$within_trib_movement[j] == 1){
        # do nothing!
      }
      # If it's a movement from mainstem to upstream:
      # For NDE - give it the same parameter as the one for the mainstem to the river mouth site (which by index, is the site before)
      # For DE - we don't care about this, because we're removing the upstream states
      else if (upper_columbia_movements$mainstem_upstream_movement[j] == 1){
        cat("borigin", i, "_vector_NDE[", parameter_indices_matrix[upper_columbia_movements$row[j], upper_columbia_movements$col[j]], "]", " = ", "borigin", i, "_matrix_", upper_columbia_movements$row[j-1], "_", upper_columbia_movements$col[j-1], "_NDE", ";", "\n", sep = "")
        
      }
      # If it's a movement from the upstream state back to the mainstem, give it the same parameter as movement from the river mouth to the upstream
      # However, this movement isn't possible in DE, so don't include those
      # And note - that these are not DE or NDE parameters
      else if (upper_columbia_movements$upstream_mainstem_movement[j] == 1){
        # cat("borigin", i, "_vector_DE[", parameter_indices_matrix[upper_columbia_movements$row[j], upper_columbia_movements$col[j]], "]", " = ", "borigin", i, "_matrix_", upper_columbia_movements$row[j-1], "_", upper_columbia_movements$col[j],  ";", "\n", sep = "")
        cat("borigin", i, "_vector_NDE[", parameter_indices_matrix[upper_columbia_movements$row[j], upper_columbia_movements$col[j]], "]", " = ", "borigin", i, "_matrix_", upper_columbia_movements$row[j-1], "_", upper_columbia_movements$col[j], ";", "\n", sep = "")
        # If it's not, store the same parameter in both matrices
      } else {
        cat("borigin", i, "_vector_DE[", parameter_indices_matrix[upper_columbia_movements$row[j], upper_columbia_movements$col[j]], "]", " = ", "borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], ";", "\n", sep = "")
        cat("borigin", i, "_vector_NDE[", parameter_indices_matrix[upper_columbia_movements$row[j], upper_columbia_movements$col[j]], "]", " = ", "borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j],  ";", "\n", sep = "")
      }
      
    }
    
    cat ("\n")
  }
  
  
  ##### write out the priors ####
  
  # b0 (intercept) priors
  for (i in 1:(nrow(b0_matrix_names))){
    # Paste the numerator
    
    # If it is a movement from the mainstem into the river mouth state, then it gets two versions of the parameter.
    # If it's not, then it just gets one.
    # If its a movement from the upstream to the river mouth state or vice versa, it doesn't get a parameter at all - 
    # we don't care about these movements.
    
    # Movements from mainstem to river mouth: print two versions
    if (b0_matrix_names$mainstem_river_mouth_movement[i] == 1){
      cat(b0_matrix_names$b0_matrix_name[i], "_DE", " ~ normal(0,5);", "\n", sep = "")
      cat(b0_matrix_names$b0_matrix_name[i], "_NDE", " ~ normal(0,5);", "\n", sep = "")
    }
    
    # If it's a within tributary movement - we don't want it. These have been
    # removed from the data (we don't care about them!)
    else if (b0_matrix_names$within_trib_movement[i] == 1){
      # do nothing!
    }
    # If it's a movement from mainstem to upstream - we don't want a parameter for it,
    # because it'll share the same parameter (see transformed parameters block) 
    # as the movement from mainstem to river mouth
    else if (b0_matrix_names$mainstem_upstream_movement[i] == 1){
      # do nothing!
    }
    
    # If we move from the river mouth state to the mainstem, we want a parameter - this is captured in the below
    # BUT - if it's a movement from the upstream state to the mainstem, we don't want a parameter.
    # these transitions will get the same parameters as the river mouth to mainstem
    else if (b0_matrix_names$upstream_mainstem_movement[i] == 1){
      # do nothing!
    } else {
      # Finally - if it's any other movement, just treat it like normal!
      # If it's not, print just the one version
      cat(b0_matrix_names$b0_matrix_name[i]," ~ normal(0,5);", "\n", sep = "")
    }
    
  }
  
  # window_spill priors
  # note that because none of these are movements into tributaries, we never need to print DE vs NDE versions
  for (i in 1:(nrow(bspillwindow_matrix_names))){
    cat(bspillwindow_matrix_names$bspillwindow_matrix_name[i]," ~ normal(0,5);", "\n", sep = "")
  }
  
  # winter spill priors
  # note that because none of these are movements into tributaries, we never need to print DE vs NDE versions
  for (i in 1:(nrow(bspill_month_matrix_names))){
    cat(bspill_month_matrix_names$bwinterspill_matrix_name[i]," ~ normal(0,5);", "\n", sep = "")
  }
  
  
  
  # btemp priors

  # Note that this block doesn't change between models for different ESUs
  for (i in 1:(nrow(btemp_matrix_names))){
    # Paste the numerator
    
    # If it is a movement from the mainstem into the river mouth state, then it gets two versions of the parameter.
    # If it's not, then it just gets one.
    # If its a movement from the upstream to the river mouth state or vice versa, it doesn't get a parameter at all - 
    # we don't care about these movements.
    
    # Movements from mainstem to river mouth: print two versions
    if (btemp_matrix_names$mainstem_river_mouth_movement[i] == 1){
      cat(btemp_matrix_names$btemp0_matrix_name[i], "_DE", " ~ normal(0,5);", "\n", sep = "")
      cat(btemp_matrix_names$btemp0_matrix_name[i], "_NDE", " ~ normal(0,5);", "\n", sep = "")
      cat(btemp_matrix_names$btemp1_matrix_name[i], "_DE", " ~ normal(0,5);", "\n", sep = "")
      cat(btemp_matrix_names$btemp1_matrix_name[i], "_NDE", " ~ normal(0,5);", "\n", sep = "")
    }
    
    # If it's a within tributary movement - we don't want it. These have been
    # removed from the data (we don't care about them!)
    else if (btemp_matrix_names$within_trib_movement[i] == 1){
      # do nothing!
    }
    # If it's a movement from mainstem to upstream - we don't want a parameter for it,
    # because it'll share the same parameter (see transformed parameters block) 
    # as the movement from mainstem to river mouth
    else if (btemp_matrix_names$mainstem_upstream_movement[i] == 1){
      # do nothing!
    }
    
    # If we move from the river mouth state to the mainstem, we want a parameter - this is captured in the below
    # BUT - if it's a movement from the upstream state to the mainstem, we don't want a parameter.
    # these transitions will get the same parameters as the river mouth to mainstem
    else if (btemp_matrix_names$upstream_mainstem_movement[i] == 1){
      # do nothing!
    } else {
      # Finally - if it's any other movement, just treat it like normal!
      # If it's not, print just the one version
      cat(btemp_matrix_names$btemp0_matrix_name[i]," ~ normal(0,5);", "\n", sep = "")
      cat(btemp_matrix_names$btemp1_matrix_name[i]," ~ normal(0,5);", "\n", sep = "")
    }
    
  }
  
  
  # Priors for btempxorigin (interaction) parameters
  # print only for the movements that are within the ESU
  for (i in 1:3){ # since we only have three origins for wild
    for (j in 1:nrow(upper_columbia_tempxorigin)){
      
      # Movements from mainstem to river mouth: print two versions
      if (upper_columbia_tempxorigin$mainstem_river_mouth_movement[j] == 1){
        cat("btemp0xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], "_DE", " ~ normal(0,5);", "\n", sep = "")
        cat("btemp0xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], "_NDE", " ~ normal(0,5);", "\n", sep = "")
        cat("btemp1xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], "_DE", " ~ normal(0,5);", "\n", sep = "")
        cat("btemp1xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], "_NDE", " ~ normal(0,5);", "\n", sep = "")
      }
      
      # If it's a within tributary movement - we don't want it
      else if (upper_columbia_tempxorigin$within_trib_movement[j] == 1){
        # do nothing!
      }
      # If it's a movement from mainstem to upstream - we also don't want it
      else if (upper_columbia_tempxorigin$mainstem_upstream_movement[j] == 1){
        # do nothing!
      } 
      # If we move from the river mouth state to the mainstem, we want a parameter - this is captured in the below
      # BUT - if it's a movement from the upstream state to the mainstem, we don't want a parameter.
      # these transitions will get the same parameters as the river mouth to mainstem
      else if (upper_columbia_tempxorigin$upstream_mainstem_movement[j] == 1){
        # do nothing!
      } else {
        # Finally - if it's any other movement, just treat it like normal!
        # If it's not, print just the one version
        cat("btemp0xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], " ~ normal(0,5);", "\n", sep = "")
        cat("btemp1xorigin", i, "_matrix_", upper_columbia_tempxorigin$row[j], "_", upper_columbia_tempxorigin$col[j], " ~ normal(0,5);", "\n", sep = "")
      }
      
    }
    cat ("\n")
  }
  
  # In this ESU, there are three origins, so we will have three origin parameters
  # We will only allow an origin effect within the ESU (so after PRA). Before, they will only have an intercept term
  # print declaration of borigin parameters
  for (i in 1:3){ # since we only have three origins for wild, we will only go i in 1:3
    for (j in 1:nrow(upper_columbia_movements)){
      
      # Movements from mainstem to river mouth: print two versions
      if (upper_columbia_movements$mainstem_river_mouth_movement[j] == 1){
        cat("borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], "_DE", " ~ normal(0,5);", "\n", sep = "")
        cat("borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], "_NDE", " ~ normal(0,5);", "\n", sep = "")
      }
      
      # If it's a within tributary movement - we don't want it
      else if (upper_columbia_movements$within_trib_movement[j] == 1){
        # do nothing!
      }
      # If it's a movement from mainstem to upstream - we also don't want it
      else if (upper_columbia_movements$mainstem_upstream_movement[j] == 1){
        # do nothing!
      } 
      # If we move from the river mouth state to the mainstem, we want a parameter - this is captured in the below
      # BUT - if it's a movement from the upstream state to the mainstem, we don't want a parameter.
      # these transitions will get the same parameters as the river mouth to mainstem
      else if (upper_columbia_movements$upstream_mainstem_movement[j] == 1){
        # do nothing!
      } else {
        # Finally - if it's any other movement, just treat it like normal!
        # If it's not, print just the one version
        cat("borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], " ~ normal(0,5);", "\n", sep = "")
      }
      
    }
    cat ("\n")
  }
  
  ## Priors for year parameters - similar setup as origin effects, but with an interaction with year
  # However, we are now dropping any movements from tributaries back into the mainstem,
  # because those movements don't have enough data to estimate them
  
  
  # Now for each origin, first the priors for the raw vector
  for (i in 1:3){ # since we only have three origins for wild
    for (j in 1:nrow(upper_columbia_RE_year)){
      
      # Movements from mainstem to river mouth: print two versions
      if (upper_columbia_RE_year$mainstem_river_mouth_movement[j] == 1){
        cat("byearxorigin", i, "_raw_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_DE", " ~ normal(0,1);", "\n", sep = "")
        cat("byearxorigin", i, "_raw_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_NDE", " ~ normal(0,1);", "\n", sep = "")
      }
      
      # If it's a within tributary movement - we don't want it
      else if (upper_columbia_RE_year$within_trib_movement[j] == 1){
        # do nothing!
      }
      # If it's a movement from mainstem to upstream - we also don't want it - because this is a shared parameter with movement from mainstem to river mouth
      else if (upper_columbia_RE_year$mainstem_upstream_movement[j] == 1){
        # do nothing!
      } 
      # There are no year effects for movements from tributaries, but this doesn't matter because it should have been removed earlier
      else if (upper_columbia_RE_year$upstream_mainstem_movement[j] == 1){
        # do nothing!
      } else {
        # Finally - if it's any other movement, just treat it like normal!
        # If it's not, print just the one version
        cat("byearxorigin", i, "_raw_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], " ~ normal(0,1);", "\n", sep = "")
      }
      
    }
    cat ("\n")
  }
  
  ## Then priors on the sigma parameters
  for (i in 1:3){ # since we only have three origins for wild
    for (j in 1:nrow(upper_columbia_RE_year)){
      
      # If this is a post-overshoot fallback movement, print out a commented out version
      # Origin 1: Wenatchee River: 7-6 and 6-5 are both post-overshoot fallback movements within the DPS
      # Origin 2: Okanogan River: no post-overshoot fallback movements within the DPS
      # Origin 3: Methow River: no post-overshoot fallback movements within the DPS
      
      # this is for the Wenatchee
      if(upper_columbia_RE_year$row[j] == 7 & upper_columbia_RE_year$col[j] == 6 & i == 1 |
         upper_columbia_RE_year$row[j] == 6 & upper_columbia_RE_year$col[j] == 5 & i == 1){
        cat("// sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_DE", " ~ cauchy(0,1);", "\n", sep = "")
      }
      
      # Movements from mainstem to river mouth: print two versions
      else if (upper_columbia_RE_year$mainstem_river_mouth_movement[j] == 1){
        cat("sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_DE", " ~ cauchy(0,1);", "\n", sep = "")
        cat("sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], "_NDE", " ~ cauchy(0,1);", "\n", sep = "")
      }
      
      # If it's a within tributary movement - we don't want it
      else if (upper_columbia_RE_year$within_trib_movement[j] == 1){
        # do nothing!
      }
      # If it's a movement from mainstem to upstream - we also don't want it
      else if (upper_columbia_RE_year$mainstem_upstream_movement[j] == 1){
        # do nothing!
      } 
      # If we move from the river mouth state to the mainstem, we want a parameter - this is captured in the below
      # BUT - if it's a movement from the upstream state to the mainstem, we don't want a parameter.
      # these transitions will get the same parameters as the river mouth to mainstem
      else if (upper_columbia_RE_year$upstream_mainstem_movement[j] == 1){
        # do nothing!
      } else {
        # Finally - if it's any other movement, just treat it like normal!
        # If it's not, print just the one version
        cat("sigma_yearxorigin", i, "_vector_", upper_columbia_RE_year$row[j], "_", upper_columbia_RE_year$col[j], " ~ cauchy(0,1);", "\n", sep = "")
      }
      
    }
    cat ("\n")
  }
  
  
  
  
  ##### Take all outputs from above and export as a single data file for model fitting #####
  
  # step 0: data in a list #
  data <- list(y = state_data_2, n_ind = n.ind, n_obs = n.obs, possible_movements = possible_movements,
               states_mat = states_mat, max_visits = dim(state_data_2)[2],
               movements = movements, not_movements = not_movements,
               nmovements = nmovements, parameter_indices_matrix = parameter_indices_matrix,
               n_states = nrow(transition_matrix), n_temp_days = nrow(temp_data),
               transition_dates = transition_date_numeric, transition_seasons_vector = transition_seasons_vector, temperature_data = as.matrix(temp_data),
               spill_window_data = as.matrix(spill_data), winter_spill_days_data = as.matrix(winter_spill_days),
               winter_post_overshoot_vector = winter_post_overshoot_vector,
               n_notmovements = n_notmovements, possible_states = transition_matrix, cat_X_mat = cat_X_mat_actual, temp_X_mat = temp_X_mat_actual, year_X_mat = year_X_mat_actual,
               grainsize = 1, N = dim(state_data_2)[1],
               # New data for detection efficiency
               tributary_design_matrices_array = tributary_design_matrices_array,
               ntransitions = ntransitions,
               transition_run_years = transition_run_years,
               nyears = max(transition_run_years),
               run_year_DE_array = run_year_DE_array,
               det_eff_param_posteriors = det_eff_param_posteriors)
  
  
  
  save(data, file = "Stan/upper_columbia_hatchery/UCH_model_data.rda")