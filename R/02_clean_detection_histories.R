#### 02: clean detection histories

# load libraries
library(tidyverse)
library(lubridate)
library(here)

#### Step 1: Run complete tag history files through create_detection_history.R script ####
# This step can either be run locally or through a HPC

# Did you run the scripts in R/hyak/02_clean_detection_histories ?
HPC <- TRUE

if (HPC == TRUE) {
  print("Please consult the readme file and run the jobs found in R/hyak to generate the following files")
  det_hist_1 <- read.csv("intermediate_outputs/CTH_1_4_complete_det_hist.csv")
  det_hist_2 <- read.csv("intermediate_outputs/CTH_5_8_complete_det_hist.csv")
  det_hist_3 <- read.csv("intermediate_outputs/CTH_9_13_complete_det_hist.csv")
  
  event_site_metadata_1 <- read.csv("intermediate_outputs/CTH_1_4_event_site_metadata.csv")
  event_site_metadata_2 <- read.csv("intermediate_outputs/CTH_5_8_event_site_metadata.csv")
  event_site_metadata_3 <- read.csv("intermediate_outputs/CTH_9_13_event_site_metadata.csv")

  event_site_metadata_1 %>%
    bind_rows(., event_site_metadata_2) %>%
    bind_rows(., event_site_metadata_3) -> complete_event_site_metadata

  write.csv(complete_event_site_metadata, here::here("intermediate_outputs", "complete_event_site_metadata.csv"))
} else {
  # source function
  source("R/functions/create_detection_history.R")
  
  file_paths_1 <- c("Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_1_2025-01-28.csv",
                    "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_2_2025-01-28.csv",
                    "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_3_2025-01-28.csv",
                    "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_4_2025-01-28.csv")
  
  file_paths_2 <- c("Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_5_2025-01-28.csv",
                    "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_6_2025-01-28.csv",
                    "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_7_2025-01-28.csv",
                    "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_8_2025-01-28.csv")
  
  file_paths_3 <- c("Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_9_2025-01-28.csv",
                    "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_10_2025-01-28.csv",
                    "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_11_2025-01-28.csv",
                    "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_12_2025-01-28.csv",
                    "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_13_2025-01-28.csv")
  
  output_1 <- create_detection_history(file_paths = file_paths_1)
  output_2 <- create_detection_history(file_paths = file_paths_2)
  output_3 <- create_detection_history(file_paths = file_paths_3)
  
  det_hist_1 <- output1$det_hist
  det_hist_2 <- output2$det_hist
  det_hist_3 <- output3$det_hist
  
  event_site_metadata_1 <- output1$event_site_metadata
  event_site_metadata_2 <- output2$event_site_metadata
  event_site_metadata_3 <- output3$event_site_metadata
  
  # join complete event site metadata and export
  
  event_site_metadata_1 %>% 
    bind_rows(., event_site_metadata_2) %>% 
    bind_rows(., event_site_metadata_3) -> complete_event_site_metadata
  
  write.csv(complete_event_site_metadata, here::here("intermediate_outputs", "complete_event_site_metadata.csv"))
}

#### Step 2: Join outputs and reformat

# Join detection histories
det_hist_1 %>% 
  bind_rows(., det_hist_2) %>% 
  bind_rows(., det_hist_3) -> det_hist_complete

# fix and convert datetimes
# some are missing a time, which we will assume to be 00:00:00
det_hist_complete %>% 
  mutate(start_time = ifelse(nchar(start_time) == 10, paste0(start_time, " 00:00:00"), start_time)) %>% 
  mutate(start_time = ymd_hms(start_time)) %>% 
  mutate(end_time = ifelse(nchar(end_time) == 10, paste0(end_time, " 00:00:00"), end_time)) %>% 
  mutate(end_time = ymd_hms(end_time)) -> complete_det_hist

#### Post-hoc processing fixes ####

# Fix 1: Fish that are seen in the exit antennas at BO2, but then skip BO4, are actually not aborting, but just missed detections at BO4
complete_det_hist %>% 
  mutate(non_ascent = ifelse(end_ant_group %in% c("BO2-WEIR 52", "BO2-WEIR 51", "BO2-UMT Entrance"), NA, non_ascent)) -> complete_det_hist

# Fix 2: Fish that are seen in the exit antennas at BO3, but then skip BO4, are actually not aborting, but just missed detections at BO4
# Use the exit coils, same as Susannah
complete_det_hist %>% 
  mutate(non_ascent = ifelse(end_ant_group %in% c("BO3-WEIR 59", "BO3-WEIR 58"), NA, non_ascent)) -> complete_det_hist

# Fixes 3-5: BO1 issues
# note sites that are between BON and MCN
BON_MCN_inriver <- c("COLR4 - Columbia River - Bonneville Dam to John Day Dam (km 234-347)",
                     "The Dalles Adult Fishways (combined)", "JDJ - John Day Dam Juvenile",
                     "LMILIS - Little Miller Island, Columbia River", "TDLPI - Lone Pine Island and associated unnamed islands near The Dalles Dam",
                     "COLR5 - Columbia River - John Day Dam to Snake River (km 347-522)",
                     # These are all new since 2017
                     "JO2 - John Day North Fish Ladder", "JDALD1 - JDA - Release into south fish ladder", 
                     "JO1 - John Day South Fish Ladder", "JDALD2 - JDA - Release into north fish ladder",
                     "John Day Dam Adult Fishways (combined)")


complete_det_hist %>% 
  mutate(to_remove = "no") -> complete_det_hist

# First, let's look for all of the cases where there are more than two events, and remove all but the first and last ones
# For example: 3D9.1C2C8BAB16 - 7 events that ultimately are all part of the same aborted attempt
complete_det_hist %>% 
  mutate(to_remove = ifelse(event_site_name == "BO1 - Bonneville Bradford Is. Ladder" & # make sure all three events are at BO1
                              lag(event_site_name) == "BO1 - Bonneville Bradford Is. Ladder" &
                              lead(event_site_name) == "BO1 - Bonneville Bradford Is. Ladder" &
                              lag(tag_code) == tag_code & lead(tag_code) == tag_code & #  make sure it's the same fish))
                              !(end_antenna_id %in% c("01", "02")), # don't select it if it's seen exiting the ladder
                            "yes", "no")) -> complete_det_hist

# Manually inspect all of these
complete_det_hist %>% 
  group_by(tag_code) %>% 
  filter(any(to_remove == "yes")) -> complete_det_hist_test_removes
# these look good!

# drop them
complete_det_hist %>% 
  subset(to_remove != "yes") -> complete_det_hist


for (i in 1:(nrow(complete_det_hist)-1)){
  # Fix 3: Fish that are taking a very long time to pass through BO1
  # if:
  if (complete_det_hist[i, 'tag_code'] == complete_det_hist[i+1, 'tag_code'] &  # 1) two rows are from the same fish 
      complete_det_hist[i, 'event_site_name'] == "BO1 - Bonneville Bradford Is. Ladder" & # 2) and are both at BO1
      complete_det_hist[i+1, 'event_site_name'] == "BO1 - Bonneville Bradford Is. Ladder" &
      complete_det_hist[i, 'end_ant_group'] %in% c("A-BRANCH WEIR 50", "A-BRANCH WEIR 49", "A-BRANCH WEIR 48", # 3) and the first event is in the lower weirs
                                                   "B-BRANCH WEIR 50", "B-BRANCH WEIR 49", "B-BRANCH WEIR 48") &
      complete_det_hist[i+1, 'end_ant_group'] == "VERTICAL SLOT DETECTORS") { # 4) and the second is in the upper weirs
    
    # Change non-ascent to not aborted
    complete_det_hist[i, "non_ascent"] <- NA
    
    # Move the end information from the subsequent row to the current row
    complete_det_hist[i, 'end_time'] <- complete_det_hist[i+1, 'end_time']
    complete_det_hist[i, 'end_antenna_id'] <- complete_det_hist[i+1, 'end_antenna_id']
    complete_det_hist[i, 'end_ant_group'] <- complete_det_hist[i+1, 'end_ant_group']
    
    # then, identify rows to remove - these are the rows following the current, from
    # which we have pulled the necessary information
    complete_det_hist[i+1, "to_remove"] <- "yes"
    
  }
  
  # Fix 4: Fish that are missed detections in vertical slot detectors at exit of BO1
  else if (complete_det_hist[i, 'tag_code'] == complete_det_hist[i+1, 'tag_code'] &  # 1) two rows are from the same fish 
           complete_det_hist[i, 'event_site_name'] == "BO1 - Bonneville Bradford Is. Ladder" & # the first is at BO1 and 
           complete_det_hist[i, 'end_ant_group'] %in% c("A-BRANCH WEIR 50", "A-BRANCH WEIR 49", "A-BRANCH WEIR 48", # the detection is in the upper portion of the lower weir
                                                        "B-BRANCH WEIR 50", "B-BRANCH WEIR 49", "B-BRANCH WEIR 48") &
           complete_det_hist[i+1, 'event_site_name'] %in% c(BON_MCN_inriver, "MC1 - McNary Oregon Shore Ladder",  "MC2 - McNary Washington Shore Ladder")){ # and the next row is upstream of BON
    
    # Change non-ascent to not aborted
    complete_det_hist[i, "non_ascent"] <- NA
    
  }
  
  # Fix 5: Fish that go part way up the ladder to the junction pool, then hang out and turn around a while later
  else if (complete_det_hist[i, 'tag_code'] == complete_det_hist[i+1, 'tag_code'] &  # 1) two rows are from the same fish 
           complete_det_hist[i, 'event_site_name'] == "BO1 - Bonneville Bradford Is. Ladder" & # 2) and are both at BO1
           complete_det_hist[i+1, 'event_site_name'] == "BO1 - Bonneville Bradford Is. Ladder" &
           complete_det_hist[i, 'end_ant_group'] %in% c("A-BRANCH WEIR 50", "A-BRANCH WEIR 49", "A-BRANCH WEIR 48", 
                                                        "B-BRANCH WEIR 50", "B-BRANCH WEIR 49", "B-BRANCH WEIR 48") & # 3) the first row ends in the upper portion of the lower weir
           complete_det_hist[i+1, 'start_ant_group'] %in% c("A-BRANCH WEIR 50", "A-BRANCH WEIR 49", "A-BRANCH WEIR 48", 
                                                            "B-BRANCH WEIR 50", "B-BRANCH WEIR 49", "B-BRANCH WEIR 48") & # 4) the second detection starts in the upper portion of the lower weir and ends in the exit
           complete_det_hist[i+1, 'end_ant_group'] %in% c("A-BRANCH WEIR 43", "A-BRANCH WEIR 44", "A-BRANCH WEIR 45", 
                                                          "B-BRANCH WEIR 43", "B-BRANCH WEIR 44", "B-BRANCH WEIR 45")){
    
    # Make sure it's aborted
    complete_det_hist[i, "non_ascent"] <- "aborted"
    
    # Move the end information from the subsequent row to the current row
    complete_det_hist[i, 'end_time'] <- complete_det_hist[i+1, 'end_time']
    complete_det_hist[i, 'end_antenna_id'] <- complete_det_hist[i+1, 'end_antenna_id']
    complete_det_hist[i, 'end_ant_group'] <- complete_det_hist[i+1, 'end_ant_group']
    
    # then, identify rows to remove - these are the rows following the current, from
    # which we have pulled the necessary information
    complete_det_hist[i+1, "to_remove"] <- "yes"
    
    
  }
}

# Fix 7: Fix rows that are incorrectly assigned at BO1 descents, that are really aborts
complete_det_hist %>% 
  mutate(non_ascent = ifelse(event_site_name == "BO1 - Bonneville Bradford Is. Ladder" & non_ascent == "descent" & 
                               start_antenna_id == "04" & end_ant_group %in% c("A-BRANCH WEIR 43", "A-BRANCH WEIR 44", "A-BRANCH WEIR 45", 
                                                                               "B-BRANCH WEIR 43", "B-BRANCH WEIR 44", "B-BRANCH WEIR 45"),
                             "aborted", non_ascent)) -> complete_det_hist








# check if fix for BO1 worked
# Look into BO1 fish
complete_det_hist %>% 
  group_by(tag_code) %>% 
  filter(any(end_ant_group %in% c("A-BRANCH WEIR 50", "A-BRANCH WEIR 50",
                                  "B-BRANCH WEIR 50", "B-BRANCH WEIR 50"))) -> BO1_branch_ends

# Now, remove those rows
complete_det_hist %>% 
  subset(to_remove == "no") %>% 
  dplyr::select(-to_remove) -> complete_det_hist_postprocessed

# Export this file

write.csv(complete_det_hist_postprocessed, here::here("intermediate_outputs", "complete_det_hist_postprocessed.csv"))

# Export a record of how many fish were detected at each site

complete_det_hist_postprocessed %>% 
  # Make a correction for JDA - shouldn't be necessary later once the script is re-run
  mutate(event_site_name = ifelse(event_site_name %in% c("JO1 - John Day South Fish Ladder", "JO2 - John Day North Fish Ladder"),
                                  "John Day Dam Adult Fishways (combined)", event_site_name)) %>% 
  # mutate(event_site_latitude = ifelse(event_site_name == "John Day Dam Adult Fishways (combined)", 45.71866, event_site_latitude)) %>% 
  # mutate(event_site_longitude = ifelse(event_site_name == "John Day Dam Adult Fishways (combined)", -120.6978, event_site_longitude)) %>% 
  group_by(event_site_name) %>% 
  summarise(n()) %>% 
  dplyr::rename(count = `n()`)-> event_det_counts

event_det_counts %>% 
  left_join(., event_site_metadata, by = "event_site_name") -> event_det_counts

# Get another df for dams
# Note: this is outdated, but also not important to the output
event_det_counts %>% 
  mutate(dam = ifelse(event_site_name %in% c("Bonneville Adult Fishways (combined)",
                                             "McNary Adult Fishways (combined)",
                                             "Priest Rapids Adult Fishways (combined)",
                                             "Rock Island Adult Fishways (combined)", 
                                             "RRF - Rocky Reach Fishway", 
                                             "Wells Dam Adult Fishways (combined)",
                                             "Ice Harbor Adult Fishways (combined)",  
                                             "Lower Granite Dam Adult Fishways (combined)",
                                             "John Day Dam Adult Fishways (combined)",
                                             # Dams without consistent PIT tag detectors
                                             # Missing John Day for this dataset (installed 2017)
                                             "The Dalles Adult Fishways (combined)",
                                             "LMA - Lower Monumental Adult Ladders",
                                             "GOA - Little Goose Fish Ladder"), "dam",
                      "in stream")) %>% 
  # Get a field for dam abbreviations
  mutate(dam_abbr = ifelse(event_site_name == "Bonneville Adult Fishways (combined)", "BON",
                           ifelse(event_site_name == "McNary Adult Fishways (combined)", "MCN",
                                  ifelse(event_site_name == "Priest Rapids Adult Fishways (combined)", "PRA",
                                         ifelse(event_site_name == "Rock Island Adult Fishways (combined)", "RIS",
                                                ifelse(event_site_name == "RRF - Rocky Reach Fishway", "RRE",
                                                       ifelse(event_site_name == "Wells Dam Adult Fishways (combined)", "WEL",
                                                              ifelse(event_site_name == "Ice Harbor Adult Fishways (combined)", "ICH",
                                                                     ifelse(event_site_name == "Lower Granite Dam Adult Fishways (combined)", "LGR",
                                                                            # Dams without consistent PIT tag detectors
                                                                            ifelse(event_site_name == "John Day Dam Adult Fishways (combined)", "JDA",
                                                                                   ifelse(event_site_name == "The Dalles Adult Fishways (combined)", "TDA",
                                                                                          ifelse(event_site_name == "LMA - Lower Monumental Adult Ladders", "LMO",
                                                                                                 ifelse(event_site_name == "GOA - Little Goose Fish Ladder", "LGO", NA))))))))))))) -> event_det_counts

# the old version somehow has lat/lons for the sites and this doesn't, so let's
# leave the old one in
# write.csv(event_det_counts,  "intermediate_outputs/complete_event_det_counts.csv", row.names = FALSE)


