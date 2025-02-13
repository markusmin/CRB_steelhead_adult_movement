#### create_detection_history

# Author: Markus Min
# Last updated: 2025-01-29

## Description:
# This function takes as input the raw PTAGIS query from 01_query_data.Rmd
# and outputs a detection history, where detections at multiple antennas
# at a single site are collapsed into one observation at a site, along with
# a start and end time at the site and the direction of movement

# this function takes as an input a list of file paths and runs each of them
# through the processing script

### Load libraries
library(tidyverse)
library(lubridate)
# library(janitor)

# for testing:
# file_paths <- c("Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_1_2025-01-28.csv",
#                 "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_2_2025-01-28.csv",
#                 "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_3_2025-01-28.csv",
#                 "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_4_2025-01-28.csv")

create_detection_history <- function(file_paths){
  ##### Store the antenna configurations #####
  
  BO1_110_upper <- c("01", "02", "03", "04")
  
  BO1_110_lower <- c("05", "06", "07", "08", "09", "0A",
                     "0B", "0C", "0D", "0E", "0F", "10",
                     "11", "12", "13", "14",
                     "15", "16", "17", "18", "19", 
                     "1A","1B", "1C", "1D", "1E", "1F",
                     "20", "21", "22", "23", "24")
  
  BO1_120_upper <- c("01", "02", "03", "04")
  
  BO1_120_lower <- c("05", "06", "07", "08", "09", "0A",
                     "0B", "0C", "0D", "0E", "0F", "10",
                     "11", "12", "13", "14",
                     "15", "16", "17", "18", "19", 
                     "1A","1B", "1C", "1D", "1E", "1F",
                     "20", "21", "22", "23", "24")
  
  BO3_lower <- as.character(c(11, 12, 13, 14, 15, 16, 17, 18))
  
  BO3_upper <- c("01", "02", "03", "04", "05", "06", "07", "08" , "09", "0A",
                 "0B", "0C", "0D", "0E", "0F", 10)
  
  BO3_AFF <- as.character(c(22, 24, 26, 28, 32, 34, 36, 38))
  
  MC1_upper <- c("01", "02")
  MC1_lower <- c("03", "04", "05", "06", "07", "08", "09",
                 "0A", "0B", "0C", "0D", "0E", "0F", 10, 11, 12)
  MC1_entrance <- c("0D", "0E", "0F", 10, 11, 12)
  
  MC2_120_upper <- c("01", "02", "03")
  MC2_120_lower <- c("04", "05", "06", "07", "08", "09",
                     "0A", "0B", "0C", "0D", "0E", "0F", 10, 11, 12, 13)
  MC2_120_entrance <- c("0E", "0F", 10, 11, 12, 13)
  
  PRA_100_east <- c("01", "02", "03", "04")
  PRA_100_west <- c("05", "06", "07", "08")
  
  PRA_110_east <- c("01", "02", "03", "04")
  PRA_110_west <- c("05", "06", "07", "08")
  PRA_110_AFF <- c("A1", "A2", "A3")
  
  RIA_100_left <- c("01", "02", "03", "04")
  RIA_100_middle <- c("05", "06", "07", "08")
  RIA_100_right <- c("09", "0A", "0B", "0C")
  
  RIA_110_left <- c("01", "02", "03", "04")
  RIA_110_middle <- c("05", "06", "07", "08")
  RIA_110_right <- c("09", "0A", "0B", "0C", "A1")
  
  RIA_120_left <- c("01", "02", "03", "04")
  RIA_120_middle <- c("05", "06", "07", "08")
  RIA_120_right <- c("09", "0A", "0B")
  
  RIA_130_left <- c("01", "02", "03", "04")
  RIA_130_middle <- c("05", "06", "07", "08")
  RIA_130_right <- c("09", "0A", "0B", "2C", "2D")
  
  RIA_140_left <- c("01", "02", "03", "04")
  RIA_140_middle <- c("05", "06", "07", "08")
  RIA_140_right <- c("09", "0A", "0B", "2A", "2B", "2C", "2D",
                     31, 32, 33, 34, 35, 36, 37, 38)
  
  RRF_100 <- c("01", "02", "03", "04")
  
  RRF_110 <- c("A1", "A2", "A3", "A4")
  RRF_110_trap <- c("A0")
  
  WEA_110_left <- c("01", "02", "03", "04")
  WEA_110_right <- c("05", "06", "07", "08")
  WEA_110_left_trap <- c("09", "0A")
  WEA_110_right_trap <- c("0B", "0C")
  
  WEA_120_left <- c("01", "02", "03", "04")
  WEA_120_right <- c("05", "06", "07", "08")
  WEA_120_left_trap <- c("09", "0A")
  WEA_120_right_trap <- c("B1", "B2", "B3")
  
  WEA_130_left_upper <- c("01", "02", "03", "04")
  WEA_130_right_upper <- c("05", "06", "07", "08")
  WEA_130_right_lower <- c("F3", "F4")
  WEA_130_left_trap <- c("09", "0A")
  WEA_130_right_trap <- c("B1", "B2", "B3")
  
  WEA_140_left_upper <- c("01", "02", "03", "04")
  WEA_140_left_lower <- c("F1", "F2")
  WEA_140_right_upper <- c("05", "06", "07", "08")
  WEA_140_right_lower <- c("F3", "F4")
  WEA_140_left_trap <- c("09", "0A")
  WEA_140_right_trap <- c("B1", "B2", "B3")
  
  WEA_150_left_upper <- c("01", "02", "03", "04")
  WEA_150_left_lower <- c("A1", "A2")
  WEA_150_right_upper <- c("05", "06", "07", "08")
  WEA_150_right_lower <- c("A3", "A4")
  WEA_150_left_trap <- c("09", "0A")
  WEA_150_right_trap <- c("B1", "B2", "B3")
  
  WEA_160_left_upper <- c("01", "02", "03", "04")
  WEA_160_left_lower <- c("A1", "A2")
  WEA_160_right_upper <- c("05", "06", "07", "08")
  WEA_160_right_lower <- c("A3", "A4")
  WEA_160_left_trap <- c("09", "0A")
  WEA_160_right_trap <- c("B1", "B2", "B3")
  WEA_160_right_AFF <- c("C1", "C2")
  
  WEA_170_left_upper <- c("01", "02", "03", "04")
  WEA_170_left_lower <- c("A1", "A2")
  WEA_170_right_upper <- c("05", "06", "07", "08")
  WEA_170_right_lower <- c("A3", "A4")
  WEA_170_left_trap <- c("09", "0A")
  WEA_170_right_trap <- c("B1", "B2", "B3")
  
  WEA_180_left_upper <- c("01", "02", "03", "04")
  WEA_180_left_lower <- c("A1", "A2")
  WEA_180_right_upper <- c("05", "06", "07", "08")
  WEA_180_right_lower <- c("A3", "A4")
  WEA_180_left_trap <- c("09", "0A")
  WEA_180_right_trap <- c("B1", "B2", "B3")
  
  ICH_100_north_ladder <- c("09", "0A", "0B", "0C", "0D", "0E", "0F", "10")
  ICH_100_south_ladder <- c("01", "02", "03", "04", "05", "06", "07", "08")
  ICH_100_bypass <- c("A1", "A2", "A3", "A4")
  
  ICH_110_north_ladder <- c("09", "0A", "0B", "0C", "0D", "0E", "0F", "10")
  ICH_110_south_ladder <- c("01", "02", "03", "04", "05", "06", "07", "08")
  ICH_110_bypass <- c("A1", "A2", "A3", "A4")
  
  
  ICH_110_south_trap <- c("F1")
  
  GRA_140_lower <- as.character(c(12, 14, 16, 18, 22, 24, 26, 28))
  GRA_140_upper <- c("01", "02", "03", "04", "05", "06", "07", "08")
  
  GRA_150_entrance <- c("B1", "B2", "B3", "B4")
  GRA_150_lower <- as.character(c(12, 14, 16, 18, 22, 24, 26, 28))
  GRA_150_upper <- c("01", "02", "03", "04", "05", "06", "07", "08")
  GRA_150_exit <- c("A1", "A2")
  
  # concatenate some of them
  
  RIA_left <- c(RIA_100_left, RIA_110_left, RIA_120_left,
                RIA_130_left, RIA_140_left)
  RIA_middle <- c(RIA_100_middle, RIA_110_middle, RIA_120_middle,
                  RIA_130_middle, RIA_140_middle)
  RIA_right <- c(RIA_100_right, RIA_110_right, RIA_120_right,
                 RIA_130_right, RIA_140_right)
  
  WEL_left <- c(WEA_110_left, WEA_110_left_trap,
                WEA_120_left, WEA_120_left_trap,
                WEA_130_left_upper, WEA_130_left_trap,
                WEA_140_left_upper, WEA_140_left_lower, WEA_140_left_trap,
                WEA_150_left_upper, WEA_150_left_lower, WEA_150_left_trap,
                WEA_160_left_upper, WEA_160_left_lower, WEA_160_left_trap,
                WEA_170_left_upper, WEA_170_left_lower, WEA_170_left_trap,
                WEA_180_left_upper, WEA_180_left_lower, WEA_180_left_trap)
  
  WEL_right <- c(WEA_110_right, WEA_110_right_trap,
                 WEA_120_right, WEA_120_right_trap,
                 WEA_130_right_upper, WEA_130_right_lower, WEA_130_right_trap,
                 WEA_140_right_upper, WEA_140_right_lower, WEA_140_right_trap,
                 WEA_150_right_upper, WEA_150_right_lower, WEA_150_right_trap,
                 WEA_160_right_upper, WEA_160_right_lower, WEA_160_right_trap, WEA_160_right_AFF,
                 WEA_170_right_upper, WEA_170_right_lower, WEA_170_right_trap,
                 WEA_180_right_upper, WEA_180_right_lower, WEA_180_right_trap)
  
  WEA_traps <- unique(c(WEA_110_left_trap, WEA_120_left_trap, WEA_130_left_trap, 
                        WEA_140_left_trap, WEA_150_left_trap, WEA_160_left_trap, 
                        WEA_170_left_trap, WEA_180_left_trap,
                        WEA_110_right_trap, WEA_120_right_trap, WEA_130_right_trap, 
                        WEA_140_right_trap, WEA_150_right_trap, WEA_160_right_trap, 
                        WEA_160_right_AFF, WEA_170_right_trap, WEA_180_right_trap))
  
  
  
  
  ### Combine files, fix some column data types, rename some columns with long names
  CTH_list <- vector(mode = "list", length = length(file_paths))
  for (i in 1:length(file_paths)){
    # CTH_list[[i]] <- clean_names(read.csv(file_paths[i]))
    CTH_list[[i]] <- read.csv(file_paths[i])
    CTH_list[[i]] %>% 
      dplyr::rename("tag_code" = "Tag.Code", 
                    "event_type_name" = "Event.Type.Name", 
                    "event_site_name" = "Event.Site.Name",
                    "event_date_mmddyyyy" = "Event.Date.MMDDYYYY", 
                    "event_release_site_name" = "Event.Release.Site.Name", 
                    "event_release_date_mmddyyyy" = "Event.Release.Date.MMDDYYYY",
                    "antenna_id" = "Antenna.ID", 
                    "antenna_group_name" = "Antenna.Group.Name", 
                    "antenna_group_configuration_value" = "Antenna.Group.Configuration.Value",
                    "event_date_time_value" = "Event.Date.Time.Value", 
                    "cth_count" = "CTH.Count") -> CTH_list[[i]]
  }
  
  bind_rows(CTH_list) %>% 
    mutate(event_date_time_value = mdy_hms(event_date_time_value)) %>%
    dplyr::rename(ant_config = antenna_group_configuration_value) -> CTH_complete

  # Get the date of arrival at BON - take min in case it fell back over Bonneville
  CTH_complete %>% 
    group_by(tag_code) %>% 
    subset(event_site_name %in% c("BHL - Adult Fishway at BONH", "BO1 - Bonneville Bradford Is. Ladder", 
                                  "BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF", 
                                  "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility")) %>% 
    filter(event_date_time_value == min(event_date_time_value)) %>% 
    dplyr::select(tag_code, event_date_time_value) %>% 
    dplyr::rename(BON_arrival = event_date_time_value) -> BON_arrival
  
  # Add run year info
  # Add run year info
  run_year <- c("05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21", "21/22", "22/23", "23/24")
  run_year_start <- seq(ymd("2005-06-01"), ymd("2023-06-01"), by = "years")
  run_year_end <- seq(ymd("2006-05-31"), ymd("2024-05-31"), by = "years")
  
  run_year_df <- data.frame(run_year, run_year_start, run_year_end)
  
  BON_arrival %>% 
    mutate(dummy =TRUE) %>% 
    left_join(run_year_df %>% mutate(dummy=TRUE), by = "dummy") %>% 
    filter(run_year_start <= BON_arrival, run_year_end >= BON_arrival) %>% 
    select(-c(dummy, run_year_start, run_year_end)) -> BON_arrival

  
  
  # Subset only the adult migration history
  CTH_complete %>% 
    left_join(., BON_arrival, by = "tag_code") %>% 
    # Keep only observations, no mark or recapture events (there is some weird stuff going on with those)
    subset(., event_type_name == "Observation") %>% 
    group_by(tag_code) %>% 
    filter(event_date_time_value >= BON_arrival) %>% 
    dplyr::select(-BON_arrival) %>% 
    # Make sure that they are in chronological order
    arrange(tag_code, event_date_time_value) -> CTH_adult
  
  CTH_adult %>% 
    group_by(tag_code) %>% 
    filter(any(event_site_name %in% c("BHL - Adult Fishway at BONH", "BON - Bonneville Dam Complex"))) -> BON_non_ladders
  
  # So it would appear that there wasn't a single detection at "BON - Bonneville Dam Complex" in our detection histories - unsure of where this came from, then
  
  # EDIT CTH_adult
  # We have to edit this so that the code recognizes that different detectors
  # are at the same site
  
  # HOWEVER: We want to keep the original event site names as well so that we don't lose that detail in later files
  CTH_adult %>% 
    mutate(event_site_name_original = event_site_name) -> CTH_adult
  
  # Let's check where the BO4 fish were previously
  CTH_adult %>% 
    mutate(last_event_site_name = lag(event_site_name)) %>% 
    subset(event_site_name == "BO4 - Bonneville WA Ladder Slots") -> BO4_fish
  
  BO4_fish %>% 
    subset(last_event_site_name != "BO4 - Bonneville WA Ladder Slots") -> BO4_fish
  
  # table(BO4_fish$last_event_site_name)
  
  # Look at fish that weren't previously at BO2 or BO3
  subset(BO4_fish, !(last_event_site_name %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF")))$tag_code -> nondirect_BO4s
  subset(BO4_fish, last_event_site_name == "BO1 - Bonneville Bradford Is. Ladder")$tag_code -> BO1_BO4s
  
  subset(BO4_fish, last_event_site_name == "DRM - Deschutes River mouth")$tag_code -> DRM_BO4s
  
  # Let's look at fish that were released into the adult fish ladder
  subset(CTH_adult, event_site_name == "LGRLDR - LGR - Release into the Adult Fish Ladder")$tag_code -> LGR_ladder_held_fish
  
  # Okay, some edits:
  # BHL is the Bonneville Hatchery ladder, which was only active 2012-2014, and is actually before the dam itself (and is not part of the adult fishways)
  # Combine the arrays for the dams that don't separate states in our model (because they weren't operational for the full study period)
  
  # ALSO: rename the ICH bypass as such, instead of leaving it in with ICH (combined)
  CTH_adult %>%
    #   mutate(event_site_name = ifelse(event_site_name %in% c("BHL - Adult Fishway at BONH", "BO1 - Bonneville Bradford Is. Ladder",
    #                                                          "BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
    #                                                          "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility",
    #                                                          "BON - Bonneville Dam Complex"),
    #                                   "Bonneville Adult Fishways (combined)", event_site_name)) %>%
    #   mutate(event_site_name = ifelse(event_site_name %in% c("MC1 - McNary Oregon Shore Ladder", "MC2 - McNary Washington Shore Ladder"),
    #                                   "McNary Adult Fishways (combined)", event_site_name)) %>% 
    mutate(event_site_name = ifelse(event_site_name %in% c("TD1 - The Dalles East Fish Ladder", "TD2 - The Dalles North Fish Ladder"),
                                    "The Dalles Adult Fishways (combined)", event_site_name)) %>%
    mutate(event_site_name = ifelse(event_site_name %in% c("JO1 - John Day South Fish Ladder", "JO2 - John Day North Fish Ladder"),
                                    "John Day Dam Adult Fishways (combined)", event_site_name)) %>% 
    mutate(event_site_name = ifelse(event_site_name == "ICH - Ice Harbor Dam (Combined)" & antenna_id %in% ICH_100_bypass, 
                                    "ICB - Ice Harbor Bypass", event_site_name)) -> CTH_adult
  #   mutate(event_site_name = ifelse(event_site_name %in% c("ICH - Ice Harbor Dam (Combined)", "IHR - Ice Harbor Dam"),
  #                                   "Ice Harbor Adult Fishways (combined)", event_site_name)) %>% 
  #   mutate(event_site_name = ifelse(event_site_name %in% c("LGRLDR - LGR - Release into the Adult Fish Ladder", "GRA - Lower Granite Dam Adult",
  #                                                          "LGR - Lower Granite Dam"),
  #                                   "Lower Granite Dam Adult Fishways (combined)", event_site_name)) %>% 
  #   mutate(event_site_name = ifelse(event_site_name %in% c("RIA - Rock Island Adult", "RIS - Rock Island Dam"),
  #                                   "Rock Island Adult Fishways (combined)", event_site_name)) %>% 
  #   mutate(event_site_name = ifelse(event_site_name %in% c("WEA - Wells Dam, DCPUD Adult Ladders", 
  #                                                          "WELLD2 - WEL - Release into the West Adult Fish Ladder",
  #                                                          "WELLD1 - WEL - Release into the East Adult Fish Ladder",
  #                                                          "WEL - Wells Dam"), "Wells Dam Adult Fishways (combined)", event_site_name)) %>% 
  #   mutate(event_site_name = ifelse(event_site_name %in% c("PRA - Priest Rapids Adult", "PRD - Priest Rapids Dam",
  #                                                          "PRDLD1 - PRD - Release into the Left Bank (facing downstream) Adult Fish Ladder"), 
  #                                   "Priest Rapids Adult Fishways (combined)", event_site_name)) -> CTH_adult
  
  
  # BON array names
  BON_arrays <- c("BHL - Adult Fishway at BONH", "BO1 - Bonneville Bradford Is. Ladder", 
                  "BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                  "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility",
                  "BON - Bonneville Dam Complex")
  
  # McNary array names
  MCN_arrays <- c("MC1 - McNary Oregon Shore Ladder", "MC2 - McNary Washington Shore Ladder")
  
  # ICH array names
  ICH_arrays <- c("ICH - Ice Harbor Dam (Combined)", "IHR - Ice Harbor Dam")
  
  # LGR array names
  LGR_arrays <- c("LGRLDR - LGR - Release into the Adult Fish Ladder", "GRA - Lower Granite Dam Adult", "LGR - Lower Granite Dam")
  
  # RIA array names
  RIS_arrays <- c("RIA - Rock Island Adult", "RIS - Rock Island Dam")
  
  # RRE array names
  RRE_arrays <- c("RRF - Rocky Reach Fishway")
  
  # WEL array names
  WEL_arrays <- c("WEA - Wells Dam, DCPUD Adult Ladders", "WELLD2 - WEL - Release into the West Adult Fish Ladder",
                  "WELLD1 - WEL - Release into the East Adult Fish Ladder", "WEL - Wells Dam")
  
  # Note that the only adult route is WEA - the others are only for juveniles
  
  # PRA array names
  PRA_arrays <- c("PRA - Priest Rapids Adult", "PRD - Priest Rapids Dam", 
                  "PRDLD1 - PRD - Release into the Left Bank (facing downstream) Adult Fish Ladder")
  
  
  ##### Convert complete tag histories into history of detection events #####
  
  
  # Get the first and last time of an event, store those
  
  
  # Get a list of the unique tag IDs
  unique_tag_IDs <- unique(CTH_adult$tag_code) 

  det_hist <- data.frame(tag_code = character(), event_type_name = character(), event_site_name = character(), 
                         start_time = as.POSIXct(character()), end_time = as.POSIXct(character()),  
                         start_antenna_id = character(), end_antenna_id = character(), 
                         start_ant_group = character(), end_ant_group = character(),
                         # antenna_group_name = character(), antenna_id = character(), ant_config = numeric())
                         ant_config = numeric(), non_ascent = character())
  
  
  # Loop through the unique tags
  print("For loop starting")
  for (i in 1:length(unique_tag_IDs)){
    # Need to run the loop for longer to actual get observations at certain dams
    # Get the start time
    if (i == 1){
      start_time <- Sys.time() 
    }
    
    # Make a new dataframe to store the history for each fish
    ind_det_hist <- data.frame(tag_code = character(), event_type_name = character(), event_site_name = character(), 
                               start_time = as.POSIXct(character()), end_time = as.POSIXct(character()), 
                               # antenna_group_name = character(), antenna_id = character(), ant_config = numeric())
                               start_ant_group = character(), end_ant_group = character(),
                               start_antenna_id = character(), end_antenna_id = character(), 
                               ant_config = numeric(), non_ascent = character())
    
    # subset the complete dataset to only this fish
    tag_hist <- subset(CTH_adult, tag_code == unique_tag_IDs[i])
    
    
    # Loop through the rows of the tag history
    for (j in 1:nrow(tag_hist)){
      # for (j in 1:36){
      
      ### For the first entry, just store these values plus the time as the start time
      if (j == 1){
        # store the tag code
        ind_det_hist[1,'tag_code'] <- unique_tag_IDs[i]
        
        # store the event type name
        ind_det_hist[1,'event_type_name'] <- tag_hist[j,'event_type_name']
        
        
        # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
        ind_det_hist[1,'ant_config'] <- tag_hist[j,'ant_config']
        
        
        # Store the start time
        ind_det_hist[1, 'start_time'] <- tag_hist[j,'event_date_time_value']
        # Store the start antenna group
        ind_det_hist[1,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
        # store the start antenna ID
        ind_det_hist[1,'start_antenna_id'] <- tag_hist[j,'antenna_id']
        
        # If it's BO2-BO3-BO4, store the site name that way, as well as the antennas
        if (tag_hist[j,'event_site_name'] %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                 "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility")){
          # site name
          ind_det_hist[1,'event_site_name'] <- 'BO2-BO3-BO4' 
          
          if (tag_hist[j,'event_site_name'] == "BO2 - Bonneville Cascades Is. Ladder"){
            # Store the start antenna group
            ind_det_hist[1,'start_ant_group'] <- paste0("BO2-", tag_hist[j,'antenna_group_name'])
            # store the start antenna ID
            ind_det_hist[1,'start_antenna_id'] <- paste0("BO2-", tag_hist[j,'antenna_id'])
          }
          else if (tag_hist[j,'event_site_name'] %in% c("BO3 - Bonneville WA Shore Ladder/AFF", "BONAFF - BON - Adult Fish Facility")){
            # Store the start antenna group
            ind_det_hist[1,'start_ant_group'] <- paste0("BO3-", tag_hist[j,'antenna_group_name'])
            # store the start antenna ID
            ind_det_hist[1,'start_antenna_id'] <- paste0("BO3-", tag_hist[j,'antenna_id'])
          }
          else if (tag_hist[j,'event_site_name'] == "BO4 - Bonneville WA Ladder Slots"){
            # Store the start antenna group
            ind_det_hist[1,'start_ant_group'] <- paste0("BO4-", tag_hist[j,'antenna_group_name'])
            # store the start antenna ID
            ind_det_hist[1,'start_antenna_id'] <- paste0("BO4-", tag_hist[j,'antenna_id'])
          }
          
        }
        # otherwise, store the site name like normal
        else {
          ind_det_hist[1,'event_site_name'] <- tag_hist[j,'event_site_name']
        }  
        
        # START THE COUNTER
        counter <- 1
        
        # SPECIAL CASE: If there is only a single detection in the entire history, then end it
        if(nrow(tag_hist) == 1){
          ind_det_hist[1,'end_time'] <- tag_hist[j,'event_date_time_value']
          ind_det_hist[1,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
          ind_det_hist[1,'end_antenna_id'] <- tag_hist[j,'antenna_id']
        }
        
        # SPECIAL CASE: If there is only one detection at the first site, store the 
        # end time as well
        
        # else if (tag_hist[j+1, 'event_site_name'] != tag_hist[j, 'event_site_name'] &
        #          # Need to make sure they're not in the BO2-BO3-BO4 complex
        #          !(tag_hist[j, 'event_site_name']  %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
        #                                                   "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility"))){
        else if (tag_hist[j+1, 'event_site_name'] != tag_hist[j, 'event_site_name']){
          
          # edit 2022-08-10: Excluding BO2-BO3-BO4 here entirely leads to problems. need to instead expand the if statement
          if (tag_hist[j, 'event_site_name']  %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                     "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility")){
            
            # So if it's seen at BO2-BO3-BO4, then we only store the end time and update the counter IF the next site is not at BO2-BO3-BO4
            if (!(tag_hist[j+1, 'event_site_name'] %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                          "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility"))){
              # Store the end time
              ind_det_hist[counter, 'end_time'] <- tag_hist[j,'event_date_time_value']
              # Store the end antenna group and ID
              ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
              ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
              
              # UPDATE THE COUNTER
              # every time we store an end time, we update the counter. This allows
              # us to move through the detection history df
              counter <- counter + 1
              
              
            }
            
          } 
          
          else {
            # Store the end time
            ind_det_hist[counter, 'end_time'] <- tag_hist[j,'event_date_time_value']
            # Store the end antenna group and ID
            ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
            ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
            
            # UPDATE THE COUNTER
            # every time we store an end time, we update the counter. This allows
            # us to move through the detection history df
            counter <- counter + 1
            
          }
          
        }
        
        
      }
      
      ### If it's the last entry, store those values
      
      else if (j == nrow(tag_hist)){
        # Store the tag code
        ind_det_hist[counter, 'tag_code'] <- tag_hist[j,'tag_code']
        # store the event type name
        ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
        # Store the end time
        ind_det_hist[counter, 'end_time'] <- tag_hist[j,'event_date_time_value']
        # Store the end antenna group
        ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
        # store the end antenna ID
        ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
        # store the antenna config
        ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
        
        
        # don't need to update the counter here because it's the last observation for the fish
        
        # Change the BO2, BO3, BO4 antenna names
        if (tag_hist[j,'event_site_name'] == "BO2 - Bonneville Cascades Is. Ladder"){
          # Store the end antenna group
          ind_det_hist[counter,'end_ant_group'] <- paste0("BO2-", tag_hist[j,'antenna_group_name'])
          # store the end antenna ID
          ind_det_hist[counter,'end_antenna_id'] <- paste0("BO2-", tag_hist[j,'antenna_id'])
        }
        else if (tag_hist[j,'event_site_name'] %in% c("BO3 - Bonneville WA Shore Ladder/AFF", "BONAFF - BON - Adult Fish Facility")){
          # Store the end antenna group
          ind_det_hist[counter,'end_ant_group'] <- paste0("BO3-", tag_hist[j,'antenna_group_name'])
          # store the end antenna ID
          ind_det_hist[counter,'end_antenna_id'] <- paste0("BO3-", tag_hist[j,'antenna_id'])
        }
        else if (tag_hist[j,'event_site_name'] == "BO4 - Bonneville WA Ladder Slots"){
          # Store the end antenna group
          ind_det_hist[counter,'end_ant_group'] <- paste0("BO4-", tag_hist[j,'antenna_group_name'])
          # store the end antenna ID
          ind_det_hist[counter,'end_antenna_id'] <- paste0("BO4-", tag_hist[j,'antenna_id'])
        }
        
        
        # SPECIAL CASE: Terminal trap detections
        # note in the non-ascent category that these were not ascents
        
        # Traps at Wells
        if (tag_hist[j, 'event_site_name'] == "WEA - Wells Dam, DCPUD Adult Ladders"){
          if (tag_hist[j, 'antenna_id'] %in% c(WEA_traps)){
            ind_det_hist[counter, 'non_ascent'] <- "trapped"
          }
          
        }
        
        
        # SPECIAL CASE: If it's the last entry AND the only detection at a site
        # Because we're changing the names of some sites in
        # the ind_det_hist file (e.g., ICH1 and ICH2, WEA1 and WEA2), we need
        # to note multiple routes of passage here. So it's just going to expand the
        # if statement quite a bit
        if (tag_hist[j-1, 'event_site_name'] != tag_hist[j, 'event_site_name'] |
            # Wells Dam
            tag_hist[j-1, 'event_site_name'] == "WEA - Wells Dam, DCPUD Adult Ladders" & tag_hist[j-1, 'antenna_id'] %in% WEL_left &
            tag_hist[j, 'event_site_name'] == "WEA - Wells Dam, DCPUD Adult Ladders" & tag_hist[j, 'antenna_id'] %in% WEL_right |
            tag_hist[j-1, 'event_site_name'] == "WEA - Wells Dam, DCPUD Adult Ladders" & tag_hist[j-1, 'antenna_id'] %in% WEL_right &
            tag_hist[j, 'event_site_name'] == "WEA - Wells Dam, DCPUD Adult Ladders" & tag_hist[j, 'antenna_id'] %in% WEL_left |
            # ICH
            tag_hist[j-1, 'event_site_name'] == "ICH - Ice Harbor Dam (Combined)" & tag_hist[j-1, 'antenna_id'] %in% ICH_110_north_ladder &
            tag_hist[j, 'event_site_name'] == "ICH - Ice Harbor Dam (Combined)" & tag_hist[j, 'antenna_id'] %in% c(ICH_110_south_ladder, ICH_110_south_trap) |
            tag_hist[j-1, 'event_site_name'] == "ICH - Ice Harbor Dam (Combined)" & tag_hist[j-1, 'antenna_id'] %in% c(ICH_110_south_ladder, ICH_110_south_trap) &
            tag_hist[j, 'event_site_name'] == "ICH - Ice Harbor Dam (Combined)" & tag_hist[j, 'antenna_id'] %in% ICH_110_north_ladder |
            # PRA
            tag_hist[j-1, 'event_site_name'] == "PRA - Priest Rapids Adult" & tag_hist[j-1, 'antenna_id'] %in% PRA_110_west &
            tag_hist[j, 'event_site_name'] == "PRA - Priest Rapids Adult" & tag_hist[j, 'antenna_id'] %in% c(PRA_110_AFF, PRA_110_east) |
            tag_hist[j-1, 'event_site_name'] == "PRA - Priest Rapids Adult" & tag_hist[j-1, 'antenna_id'] %in% c(PRA_110_AFF, PRA_110_east) &
            tag_hist[j, 'event_site_name'] == "PRA - Priest Rapids Adult" & tag_hist[j, 'antenna_id'] %in% PRA_110_west |
            # RIA
            tag_hist[j-1, 'event_site_name'] == "RIA - Rock Island Adult" & tag_hist[j-1, 'antenna_id'] %in% RIA_140_left &
            tag_hist[j, 'event_site_name'] == "RIA - Rock Island Adult" & tag_hist[j, 'antenna_id'] %in% c(RIA_140_middle, RIA_140_right) |
            tag_hist[j-1, 'event_site_name'] == "RIA - Rock Island Adult" & tag_hist[j-1, 'antenna_id'] %in% RIA_140_middle &
            tag_hist[j, 'event_site_name'] == "RIA - Rock Island Adult" & tag_hist[j, 'antenna_id'] %in% c(RIA_140_left, RIA_140_right) |
            tag_hist[j-1, 'event_site_name'] == "RIA - Rock Island Adult" & tag_hist[j-1, 'antenna_id'] %in% RIA_140_right &
            tag_hist[j, 'event_site_name'] == "RIA - Rock Island Adult"  & tag_hist[j, 'antenna_id'] %in% c(RIA_140_left, RIA_140_middle) |
            
            tag_hist[j-1, 'event_site_name'] == tag_hist[j, 'event_site_name'] & tag_hist[j, 'event_date_time_value'] -
            tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)){
          
          
          # FIRST: If it's in the BO2-BO3-BO4 complex, we need to treat it differently
          
          # But - we need to treat BO4 differently from BO2 and BO3, since BO4 basically
          # functions as antennas at the top of a ladder
          if (tag_hist[j, 'event_site_name'] %in% c("BO4 - Bonneville WA Ladder Slots")){
            
            # store the tag code
            ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
            
            # store the event type name
            ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
            
            
            # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
            ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
            
            # store the location field
            ind_det_hist[counter,'event_site_name'] <- 'BO2-BO3-BO4'
            
            # Don't overwrite what was already added - this way we keep the start time and antenna
            # that we previously recorded
            
            # Store the start time
            # ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
            # Store the start antenna
            # ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
            # store the start antenna info
            # ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
            
            # now, apply the same non-descent rules as below
            # If the last antenna seen is not in BO4, then note it was an aborted attempt
            # This now can't be true, since we took these out using the above if statement
            if (tag_hist[j, 'event_site_name'] %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                      "BONAFF - BON - Adult Fish Facility")){
              
              ind_det_hist[counter,'non_ascent'] <- "aborted"
              
            }
            
            
          }
          
          else if (tag_hist[j, 'event_site_name'] %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                         "BONAFF - BON - Adult Fish Facility")){
            
            # store the tag code
            ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
            
            # store the event type name
            ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
            
            
            # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
            ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
            
            # store the location field
            ind_det_hist[counter,'event_site_name'] <- 'BO2-BO3-BO4'
            
            # In this case, we do want to store the start time, because if we're back
            # at BO2 or BO3, then it means we've started a new ascent attempt
            
            # Store the start time
            ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
            
            if (tag_hist[j,'event_site_name'] == "BO2 - Bonneville Cascades Is. Ladder"){
              # Store the start antenna group
              ind_det_hist[counter,'start_ant_group'] <- paste0("BO2-", tag_hist[j,'antenna_group_name'])
              # store the start antenna ID
              ind_det_hist[counter,'start_antenna_id'] <- paste0("BO2-", tag_hist[j,'antenna_id'])
            }
            else if (tag_hist[j,'event_site_name'] %in% c("BO3 - Bonneville WA Shore Ladder/AFF", "BONAFF - BON - Adult Fish Facility")){
              # Store the start antenna group
              ind_det_hist[counter,'start_ant_group'] <- paste0("BO3-", tag_hist[j,'antenna_group_name'])
              # store the start antenna ID
              ind_det_hist[counter,'start_antenna_id'] <- paste0("BO3-", tag_hist[j,'antenna_id'])
            }
            
            # now, apply the same non-descent rules as below
            # If the last antenna seen is not in BO4, then note it was an aborted attempt
            # This is by default true for this section because of the above else if statement
            if (tag_hist[j, 'event_site_name'] %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                      "BONAFF - BON - Adult Fish Facility")){
              
              ind_det_hist[counter,'non_ascent'] <- "aborted"
              
            }
            
            
          }
          
          else {
            
            # Here we are using 48 hours to be in line with the dams.
            
            # store the tag code
            ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
            
            # store the event type name
            ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
            
            
            # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
            ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
            
            # store the location field
            ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
            
            # Store the start time
            ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
            # Store the start antenna
            ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
            # store the start antenna info
            ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
            
            # If it's at RIA, then replace the event site name with the correct ladder
            if(tag_hist[j, 'event_site_name'] %in% RIS_arrays){
              if(tag_hist[j, 'antenna_id'] %in% RIA_left){
                ind_det_hist[counter,'event_site_name'] <- 'RIA1 - Rock Island Adult Left Ladder'
              }
              else if(tag_hist[j, 'antenna_id'] %in% RIA_middle){
                ind_det_hist[counter,'event_site_name'] <- 'RIA2 - Rock Island Adult Middle Ladder'
              }
              else if(tag_hist[j, 'antenna_id'] %in% RIA_right){
                ind_det_hist[counter,'event_site_name'] <- 'RIA3 - Rock Island Adult Right Ladder'
              }
              
            }
            
            # If it's at ICH, then replace the event site name with the correct ladder
            if(tag_hist[j, 'event_site_name'] == "ICH - Ice Harbor Dam (Combined)"){
              if(tag_hist[j, 'antenna_id'] %in% ICH_110_north_ladder){
                ind_det_hist[counter,'event_site_name'] <- 'ICH1 - Ice Harbor Dam North Ladder'
              }
              else if(tag_hist[j, 'antenna_id'] %in% c(ICH_110_south_ladder, ICH_110_south_trap)){
                ind_det_hist[counter,'event_site_name'] <- 'ICH2 - Ice Harbor Dam South Ladder'
              }
              
            }
            
            # If it's at PRA, then replace the event site name with the correct ladder
            if(tag_hist[j, 'event_site_name'] == "PRA - Priest Rapids Adult"){
              if(tag_hist[j, 'antenna_id'] %in% PRA_110_west){
                ind_det_hist[counter,'event_site_name'] <- 'PRA1 - Priest Rapids Adult West Ladder'
              }
              else if(tag_hist[j, 'antenna_id'] %in% c(PRA_110_AFF, PRA_110_east)){
                ind_det_hist[counter,'event_site_name'] <- 'PRA2 - Priest Rapids Adult East Ladder'
              }
              
            }
            
            # If it's at WEA, then replace the event site name with the correct ladder
            if(tag_hist[j, 'event_site_name'] == "WEA - Wells Dam, DCPUD Adult Ladders"){
              if(tag_hist[j, 'antenna_id'] %in% WEL_left){
                ind_det_hist[counter,'event_site_name'] <- 'WEA1 - Wells Dam, DCPUD Left Adult Ladder'
              }
              else if(tag_hist[j, 'antenna_id'] %in% WEL_right){
                ind_det_hist[counter,'event_site_name'] <- 'WEA2 - Wells Dam, DCPUD Right Adult Ladder'
              }
              
            }
            
          }
        }
      }
      
      ### If there was only one detection at a site, then store the start and end times at the same time
      
      # See if this event site name is different than the one ahead and behind
      else if (tag_hist[j-1, 'event_site_name'] != tag_hist[j, 'event_site_name'] &
               tag_hist[j+1, 'event_site_name'] != tag_hist[j, 'event_site_name'] &
               # Need to make sure they're not in the BO2-BO3-BO4 complex
               !(tag_hist[j, 'event_site_name']  %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                        "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility"))){
        # Excluding BO2-BO3-BO4 here is problematic if we have a single detection in this whole complex, which I feel like
        # must have a very low probability of happening. I think with the old script it happened a single time, but that'll
        # be fixed in another part of the script.
        
        
        # store the tag code
        ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
        
        # store the event type name
        ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
        
        
        # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
        ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
        
        # store the location fields
        ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
        
        # Store the start time
        ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
        # Store the start antenna
        ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
        # store the start antenna info
        ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
        
        # Store the end time
        ind_det_hist[counter, 'end_time'] <- tag_hist[j,'event_date_time_value']
        # Store the end antenna
        ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
        # Store the end antenna ID
        ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
        
        # If it's at RIA, then replace the event site name with the correct ladder
        if(tag_hist[j, 'event_site_name'] %in% RIS_arrays){
          if(tag_hist[j, 'antenna_id'] %in% RIA_left){
            ind_det_hist[counter,'event_site_name'] <- 'RIA1 - Rock Island Adult Left Ladder'
          }
          else if(tag_hist[j, 'antenna_id'] %in% RIA_middle){
            ind_det_hist[counter,'event_site_name'] <- 'RIA2 - Rock Island Adult Middle Ladder'
          }
          else if(tag_hist[j, 'antenna_id'] %in% RIA_right){
            ind_det_hist[counter,'event_site_name'] <- 'RIA3 - Rock Island Adult Right Ladder'
          }
          
        }
        
        # If it's at PRA, then replace the event site name with the correct ladder
        if(tag_hist[j, 'event_site_name'] == "PRA - Priest Rapids Adult"){
          if(tag_hist[j, 'antenna_id'] %in% c(PRA_110_east, PRA_110_AFF)){
            ind_det_hist[counter,'event_site_name'] <- 'PRA2 - Priest Rapids Adult East Ladder'
          }
          else if(tag_hist[j, 'antenna_id'] %in% PRA_110_west){
            ind_det_hist[counter,'event_site_name'] <- 'PRA1 - Priest Rapids Adult West Ladder'
          }
          
        }
        
        # If it's at ICH, then replace the event site name with the correct ladder
        if(tag_hist[j, 'event_site_name'] == "ICH - Ice Harbor Dam (Combined)"){
          if(tag_hist[j, 'antenna_id'] %in% ICH_110_north_ladder){
            ind_det_hist[counter,'event_site_name'] <- 'ICH1 - Ice Harbor Dam North Ladder'
          }
          else if(tag_hist[j, 'antenna_id'] %in% c(ICH_110_south_ladder, ICH_110_south_trap)){
            ind_det_hist[counter,'event_site_name'] <- 'ICH2 - Ice Harbor Dam South Ladder'
          }
          
        }
        
        # If it's at WEA, then replace the event site name with the correct ladder
        if(tag_hist[j, 'event_site_name'] == "WEA - Wells Dam, DCPUD Adult Ladders"){
          if(tag_hist[j, 'antenna_id'] %in% WEL_left){
            ind_det_hist[counter,'event_site_name'] <- 'WEA1 - Wells Dam, DCPUD Left Adult Ladder'
          }
          else if(tag_hist[j, 'antenna_id'] %in% WEL_right){
            ind_det_hist[counter,'event_site_name'] <- 'WEA2 - Wells Dam, DCPUD Right Adult Ladder'
          }
          
        }
        
        # UPDATE THE COUNTER
        # every time we store an end time, we update the counter. This allows
        # us to move through the detection history df
        counter <- counter + 1
        
      }
      
      
      
      
      # Adult ladders
      
      ##### BONNEVILLE #####
      
      else if (tag_hist[j, 'event_site_name'] %in% BON_arrays){
        
        # BON ROUTE 1 - BO1
        
        if(tag_hist[j, 'event_site_name'] == "BO1 - Bonneville Bradford Is. Ladder"){
          
          # If it's the first detection in this route, or it hasn't been seen in this route for at least 48 hours, store the start time and antenna
          if (tag_hist[j-1, 'event_site_name'] != "BO1 - Bonneville Bradford Is. Ladder" |
              tag_hist[j-1, 'event_site_name'] == "BO1 - Bonneville Bradford Is. Ladder" &
              tag_hist[j, 'event_date_time_value'] - 
              tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)) {
            # Store the event site name
            ind_det_hist[counter, 'event_site_name'] <- tag_hist[j,'event_site_name']
            # Store the start time
            ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
            # Store the start antenna group
            ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
            # Store the start antenna ID
            ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
          }
          
          
          # Different if statements for varying antenna configurations
          
          # BO1 100
          # This one only has a single array
          else if(tag_hist[j, 'ant_config'] == 100){
            
            # If the next detection is at a different site, or is more than 48 hours later, store the current time as the end time
            if (tag_hist[j+1, 'event_site_name'] != "BO1 - Bonneville Bradford Is. Ladder" | # If it's at a different site
                tag_hist[j+1, 'event_date_time_value'] - 
                tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's within 48 hours
              
              # store the info
              # store the tag code
              ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
              
              # store the event type name
              ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
              
              
              # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
              ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
              
              # store the location fields
              ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
              # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
              # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
              # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
              # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
              
              # store the end time
              ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
              # Store the end antenna
              ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
              # store the end antenna ID
              ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
              
              # UPDATE THE COUNTER
              # every time we store an end time, we update the counter. This allows
              # us to move through the detection history df
              counter <- counter + 1
              
              
            }
            
            # if it's within 48 hours at BO1, then don't record it
            else {
              
            }
            
          }
          
          # BO1 110 & 120
          # these now have exit antennas
          # These are the same for our purposes, since the only difference is the installation of lamprey monitors for BO1 120
          else if(tag_hist[j, 'ant_config'] %in% c(110, 120)){
            
            
            # If the next detection is at a different site, or is more than 48 hours later, store the current time as the end time
            if (tag_hist[j+1, 'event_site_name'] != "BO1 - Bonneville Bradford Is. Ladder" | # If it's at a different site
                tag_hist[j+1, 'event_date_time_value'] - 
                tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's more than 48 hours later
              
              # store the info
              # store the tag code
              ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
              
              # store the event type name
              ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
              
              
              # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
              ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
              
              # store the location fields
              ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
              # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
              # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
              # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
              # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
              
              # store the end time
              ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
              # Store the end antenna
              ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
              # store the end antenna ID
              ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
              
              # However, if the last antenna seen is not an exit coil, note that it was an aborted attempt
              if (tag_hist[j, 'antenna_id'] %in% BO1_110_lower & !(ind_det_hist[counter,'start_antenna_id'] %in% BO1_110_upper)){
                
                ind_det_hist[counter,'non_ascent'] <- "aborted"
                
              }
              
              # Other alternative - if the first detection was at an exit coil, and the last at an entrance coil, then it went downstream through the ladder
              else if (tag_hist[j, 'antenna_id'] %in% BO1_110_lower & ind_det_hist[counter,'start_antenna_id'] %in% BO1_110_upper){
                ind_det_hist[counter,'non_ascent'] <- "descent"
              }
              
              # UPDATE THE COUNTER
              # every time we store an end time, we update the counter. This allows
              # us to move through the detection history df
              counter <- counter + 1
              
              
            }
            
            # if it's within 48 hours at BO1, then don't record it
            
            else {
              
            }
            
          }
          
        } 
        # ROUTE 2 - BO2-BO3-BO4 complex
        # For this route, we need to make sure that the fish came out through BO4, otherwise it's an aborted attempt.
        # We will use the same code as for BO1, but use BO2 and BO3 as the "lower" portions of the ladder
        else if(tag_hist[j, 'event_site_name'] %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                      "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility")){
          
          # if it's the first detection in this complex, or it hasn't been seen in more than 48 hours, store the start time 
          # edit 07-29-22 - because we're observing fish spending a lot of time in this complex (in some cases, over 100 days), 
          # it needs to have last been seen in BO4 in order for it to be a new event based on time difference (BUT NOT CONSECUTIVE BO4 events - see tag code 384.3B23AD9276 for example)
          
          # Need specific rule for each
          if (!(tag_hist[j-1, 'event_site_name']  %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                         "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility")) |
              # If it's currently at BO2, then call it a new event if the last attempt was any of BO2, BO3, or BO4
              # If it's currently at BO3, then call it a new event if the last attempt was any of BO2, BO3, or BO4
              tag_hist[j, 'event_site_name'] %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                    "BONAFF - BON - Adult Fish Facility") &
              tag_hist[j, 'event_date_time_value'] -
              tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)){
            # BUT - if it's currently at BO4, don't call it a new event if any of the last detections were at any of BO2, BO3, or BO4
            # this is covered by the beginning of the if statement, so we don't need to re-write it

            # Store the event site name
            ind_det_hist[counter, 'event_site_name'] <- 'BO2-BO3-BO4'
            # Store the start time
            ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
            # # Store the start antenna
            # ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
            # # Store the start antenna id
            # ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
            
            if (tag_hist[j,'event_site_name'] == "BO2 - Bonneville Cascades Is. Ladder"){
              # Store the start antenna group
              ind_det_hist[counter,'start_ant_group'] <- paste0("BO2-", tag_hist[j,'antenna_group_name'])
              # store the start antenna ID
              ind_det_hist[counter,'start_antenna_id'] <- paste0("BO2-", tag_hist[j,'antenna_id'])
            }
            else if (tag_hist[j,'event_site_name'] %in% c("BO3 - Bonneville WA Shore Ladder/AFF", "BONAFF - BON - Adult Fish Facility")){
              # Store the start antenna group
              ind_det_hist[counter,'start_ant_group'] <- paste0("BO3-", tag_hist[j,'antenna_group_name'])
              # store the start antenna ID
              ind_det_hist[counter,'start_antenna_id'] <- paste0("BO3-", tag_hist[j,'antenna_id'])
            }
            else if (tag_hist[j,'event_site_name'] == "BO4 - Bonneville WA Ladder Slots"){
              # Store the start antenna group
              ind_det_hist[counter,'start_ant_group'] <- paste0("BO4-", tag_hist[j,'antenna_group_name'])
              # store the start antenna ID
              ind_det_hist[counter,'start_antenna_id'] <- paste0("BO4-", tag_hist[j,'antenna_id'])
            }
          }
          
          # If the next detection is at a site not in BO2, BO3, or BO4,
          # or is more than 48 hours later, store the current time as the end time
          # EDIT 2022-07-29 #
          # again, we need to remove the time component in the movement to BO4 because we're seeing them spent 100+ days in here
          else if (!(tag_hist[j+1, 'event_site_name']  %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                              "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility")) | # If it's at a different site
                   tag_hist[j+1, 'event_site_name']  %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF") & # exclude BO4 from the time component here
                   tag_hist[j+1, 'event_date_time_value'] -
                   tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's over 48 hours later
            
            # old else if
            # else if (!(tag_hist[j+1, 'event_site_name']  %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
            #                                                "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility")) | # If it's at a different site
            #     tag_hist[j+1, 'event_date_time_value'] - 
            #     tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's over 48 hours later
            
            # store the info
            # store the tag code
            ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
            
            # store the event type name
            ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
            
            
            # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
            ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
            
            # store the location fields
            # ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
            # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
            # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
            # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
            # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
            
            # store the end time
            ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
            # # Store the end antenna
            # ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
            # # store the end antenna ID
            # ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
            
            if (tag_hist[j,'event_site_name'] == "BO2 - Bonneville Cascades Is. Ladder"){
              # Store the end antenna group
              ind_det_hist[counter,'end_ant_group'] <- paste0("BO2-", tag_hist[j,'antenna_group_name'])
              # store the end antenna ID
              ind_det_hist[counter,'end_antenna_id'] <- paste0("BO2-", tag_hist[j,'antenna_id'])
            }
            else if (tag_hist[j,'event_site_name'] %in% c("BO3 - Bonneville WA Shore Ladder/AFF", "BONAFF - BON - Adult Fish Facility")){
              # Store the end antenna group
              ind_det_hist[counter,'end_ant_group'] <- paste0("BO3-", tag_hist[j,'antenna_group_name'])
              # store the end antenna ID
              ind_det_hist[counter,'end_antenna_id'] <- paste0("BO3-", tag_hist[j,'antenna_id'])
            }
            else if (tag_hist[j,'event_site_name'] == "BO4 - Bonneville WA Ladder Slots"){
              # Store the end antenna group
              ind_det_hist[counter,'end_ant_group'] <- paste0("BO4-", tag_hist[j,'antenna_group_name'])
              # store the end antenna ID
              ind_det_hist[counter,'end_antenna_id'] <- paste0("BO4-", tag_hist[j,'antenna_id'])
            }
            
            # Save event site name
            ind_det_hist[counter,'event_site_name'] <- 'BO2-BO3-BO4'
            
            
            
            # However, if the last antenna seen is not in BO4, then note it was an aborted attempt
            if (tag_hist[j, 'event_site_name'] %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                      "BONAFF - BON - Adult Fish Facility")){
              
              ind_det_hist[counter,'non_ascent'] <- "aborted"
              
            }
            
            # we're not going to allow fish to go downstream from BO4 to BO3 or BO2.
            
            
            
            # UPDATE THE COUNTER
            # every time we store an end time, we update the counter. This allows
            # us to move through the detection history df
            counter <- counter + 1
            
            
          }
          
          
          
        }
        
      }
      
      
      ##### McNary #####
      
      else if (tag_hist[j, 'event_site_name'] %in% MCN_arrays){
        
        # MCN route 1 - MC1, Oregon Shore Ladder
        
        if(tag_hist[j, 'event_site_name'] == "MC1 - McNary Oregon Shore Ladder"){
          
          # If it's the first detection in this route, or it hasn't been seen at this route in at least 48 hours store the start time and antenna
          if (tag_hist[j-1, 'event_site_name'] != "MC1 - McNary Oregon Shore Ladder" |
              tag_hist[j-1, 'event_site_name'] == "MC1 - McNary Oregon Shore Ladder" &
              tag_hist[j, 'event_date_time_value'] - 
              tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)) {
            # Store the event site name
            ind_det_hist[counter, 'event_site_name'] <- tag_hist[j,'event_site_name']
            # Store the start time
            ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
            # Store the start antenna
            ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
            # Store the start antenna id
            ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
          }
          
          
          # Antenna configuration has not changed at MC1, so no need for different statements
          
          # MC1 has antennas in the lower weirs, and 2 at the counting windows
          # However, a quick look reveals fish that weren't seen at the counting window antennas but clearly ascended (384.1B796A06FF)
          
          # If the next detection is at a different site, or is more than 48 hours later, store the current time as the end time
          if (tag_hist[j+1, 'event_site_name'] != "MC1 - McNary Oregon Shore Ladder" | # If it's at a different site
              tag_hist[j+1, 'event_date_time_value'] - 
              tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's more than 48 hours later
            
            # store the info
            # store the tag code
            ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
            
            # store the event type name
            ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
            
            
            # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
            ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
            
            # store the location fields
            ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
            # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
            # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
            # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
            # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
            
            # store the end time
            ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
            # Store the end antenna
            ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
            # store the end antenna ID
            ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
            
            
            
            # However, if the last antenna seen is an entrance coil (and we didn't see it enter at an exit coil), note that it was an aborted attempt
            if (tag_hist[j, 'antenna_id'] %in% MC1_entrance & !(ind_det_hist[counter,'start_antenna_id'] %in% MC1_upper)){
              
              ind_det_hist[counter,'non_ascent'] <- "aborted"
              
            }
            
            # Other alternative - if the first detection was at an exit coil, and the last at an entrance coil, then it went downstream through the ladder
            else if (tag_hist[j, 'antenna_id'] %in% MC1_entrance & ind_det_hist[counter,'start_antenna_id'] %in% MC1_upper){
              ind_det_hist[counter,'non_ascent'] <- "descent"
            }
            
            # UPDATE THE COUNTER
            # every time we store an end time, we update the counter. This allows
            # us to move through the detection history df
            counter <- counter + 1
            
            
          }
          
          # if it's within 48 hours at MC1, then don't record it
          
          else {
            
          }
          
          
        } 
        # MCN route 2 - MC2 - McNary Washington Shore Ladder
        else if(tag_hist[j, 'event_site_name'] == "MC2 - McNary Washington Shore Ladder"){
          
          # if it's the first detection in this route, or it hasn't been seen in this route for more than 48 hours, store the start time 
          if (!(tag_hist[j-1, 'event_site_name']  %in% c("MC2 - McNary Washington Shore Ladder")) |
              tag_hist[j-1, 'event_site_name'] == "MC2 - McNary Washington Shore Ladder" &
              tag_hist[j, 'event_date_time_value'] - 
              tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)) {
            # Store the event site name
            ind_det_hist[counter, 'event_site_name'] <- tag_hist[j,'event_site_name']
            # Store the start time
            ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
            # Store the start antenna
            ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
            # Store the start antenna id
            ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
          }
          
          # Antenna configuration 110
          # there are no counting window coils in this configuration, so we'll just 
          
          else if(tag_hist[j, 'ant_config'] == 110){
            
            # If the next detection is not at MC2
            # or is more than 48 hours later, store the current time as the end time
            if (!(tag_hist[j+1, 'event_site_name']  == "MC2 - McNary Washington Shore Ladder") | # If it's at a different site
                tag_hist[j+1, 'event_date_time_value'] - 
                tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's over 48 hours later
              
              # store the info
              # store the tag code
              ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
              
              # store the event type name
              ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
              
              
              # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
              ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
              
              # store the location fields
              # ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
              # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
              # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
              # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
              # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
              
              # store the end time
              ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
              # Store the end antenna
              ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
              # store the end antenna ID
              ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
              
              # Save event site name
              ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
              
              
              
              # With this configuration of coils, we're not going to look for any strange behavior; it's also only in use for the first 9 months of our study
              
              
              # UPDATE THE COUNTER
              # every time we store an end time, we update the counter. This allows
              # us to move through the detection history df
              counter <- counter + 1
              
              
            }
            
            
          }
          
          # Antenna configuration 120
          
          else if(tag_hist[j, 'ant_config'] == 120){
            
            # If the next detection is not at MC2
            # or is more than 48 hours later, store the current time as the end time
            if (!(tag_hist[j+1, 'event_site_name']  == "MC2 - McNary Washington Shore Ladder") | # If it's at a different site
                tag_hist[j+1, 'event_date_time_value'] - 
                tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's over 48 hours later
              
              # store the info
              # store the tag code
              ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
              
              # store the event type name
              ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
              
              
              # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
              ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
              
              # store the location fields
              # ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
              # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
              # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
              # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
              # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
              
              # store the end time
              ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
              # Store the end antenna
              ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
              # store the end antenna ID
              ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
              
              # Save event site name
              ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
              
              # However, if the last antenna seen is an entrance coil (but it was never seen at the exit coil), note that it was an aborted attempt
              if (tag_hist[j, 'antenna_id'] %in% MC2_120_entrance & !(ind_det_hist[counter,'start_antenna_id'] %in% MC2_120_upper)){
                
                ind_det_hist[counter,'non_ascent'] <- "aborted"
                
              }
              
              # Other alternative - if the first detection was at an exit coil, and the last at an entrance coil, then it went downstream through the ladder
              else if (tag_hist[j, 'antenna_id'] %in% MC2_120_entrance & ind_det_hist[counter,'start_antenna_id'] %in% MC2_120_upper){
                ind_det_hist[counter,'non_ascent'] <- "descent"
              }
              
              
              
              # UPDATE THE COUNTER
              # every time we store an end time, we update the counter. This allows
              # us to move through the detection history df
              counter <- counter + 1
              
              
            }
            
            
          }
          
          
          
          
          
        }
        
      }
      
      ##### Priest Rapids #####
      
      else if (tag_hist[j, 'event_site_name'] %in% PRA_arrays){
        
        # For priest rapids, we are just going to treat the AFF antennas as antennas in the lower ladder (for the east/left ladder).
        # configuration 110, which started in June 2007, introduced these arrays
        
        # PRA route 1 - the west (right) ladder
        # In this route, the configuration has not changed through time, so no ifelse based on configuration
        
        if(tag_hist[j, 'antenna_id'] %in% PRA_110_west){
          
          # If it's the first detection in this route, or it hasn't been seen at this route in at least 48 hours store the start time and antenna
          if (!(tag_hist[j-1, 'antenna_id'] %in% PRA_110_west & tag_hist[j-1, 'event_site_name'] == "PRA - Priest Rapids Adult") |
              tag_hist[j-1, 'antenna_id'] %in% PRA_110_west & 
              tag_hist[j-1, 'event_site_name'] == "PRA - Priest Rapids Adult" &
              tag_hist[j, 'event_date_time_value'] - 
              tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)) {
            # Store the event site name
            ind_det_hist[counter, 'event_site_name'] <- "PRA1 - Priest Rapids Adult West Ladder"
            # Store the start time
            ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
            # Store the start antenna
            ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
            # Store the start antenna id
            ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
          }
          
          
          
          # If the next detection is at a different site, or is more than 48 hours later, store the current time as the end time
          if (!(tag_hist[j+1, 'antenna_id'] %in% PRA_110_west) | # If it's at a different site
              tag_hist[j+1, 'event_date_time_value'] - 
              tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's more than 48 hours later
            
            # store the info
            # store the tag code
            ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
            
            # store the event type name
            ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
            
            
            # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
            ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
            
            # store the location fields
            ind_det_hist[counter,'event_site_name'] <- "PRA1 - Priest Rapids Adult West Ladder"
            # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
            # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
            # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
            # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
            
            # store the end time
            ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
            # Store the end antenna
            ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
            # store the end antenna ID
            ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
            
            # With only four antennas at the exit here, we don't have the ability to identify any non-ascent movements
            
            # UPDATE THE COUNTER
            # every time we store an end time, we update the counter. This allows
            # us to move through the detection history df
            counter <- counter + 1
            
            
          }
          
          # if it's within 48 hours at PRA west, then don't record it
          
          else {
            
          }
          
          
        } 
        # PRA route 2 - the east (left) ladder
        # config 110 (June 2007) introduced antennas in the AFF
        else if(tag_hist[j, 'antenna_id'] %in% c(PRA_110_AFF, PRA_110_east)){
          
          
          # If it's the first detection in this route in at least 48 hours,
          if (!(tag_hist[j-1, 'antenna_id'] %in% c(PRA_110_east, PRA_110_AFF) & tag_hist[j-1, 'event_site_name'] == "PRA - Priest Rapids Adult") |
              tag_hist[j-1, 'antenna_id'] %in% c(PRA_110_east, PRA_110_AFF) &
              tag_hist[j-1, 'event_site_name'] == "PRA - Priest Rapids Adult" &
              tag_hist[j, 'event_date_time_value'] - 
              tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)) {
            # Store the event site name
            ind_det_hist[counter, 'event_site_name'] <- "PRA2 - Priest Rapids Adult East Ladder"
            # Store the start time
            ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
            # Store the start antenna
            ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
            # Store the start antenna id
            ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
          }
          
          
          # If it's not the first detection in this route
          
          # config 100 - nothing in the AFF, so we can't identify non-ascents
          else if (tag_hist[j, 'ant_config'] == 100){ 
            
            # If the next detection is not at PRA, east
            # or is more than 48 hours later, store the current time as the end time
            if (!(tag_hist[j+1, 'event_site_name']  == "PRA - Priest Rapids Adult") & !(tag_hist[j+1, 'antenna_id'] %in% c(PRA_110_east, PRA_110_AFF)) | # If it's at a different site
                tag_hist[j+1, 'event_date_time_value'] - 
                tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's over 48 hours later
              
              # store the info
              # store the tag code
              ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
              
              # store the event type name
              ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
              
              
              # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
              ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
              
              # store the location fields
              ind_det_hist[counter,'event_site_name'] <- "PRA2 - Priest Rapids Adult East Ladder"
              # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
              # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
              # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
              # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
              
              # store the end time
              ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
              # Store the end antenna
              ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
              # store the end antenna ID
              ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
              
              # With only four antennas at the exit here, we don't have the ability to identify any non-ascent movements
              
              # UPDATE THE COUNTER
              # every time we store an end time, we update the counter. This allows
              # us to move through the detection history df
              counter <- counter + 1
              
              
            }
            
          }
          
          else if (tag_hist[j, 'ant_config'] == 110){ 
            
            # If the next detection is not at PRA, east
            # or is more than 48 hours later, store the current time as the end time
            if (!(tag_hist[j+1, 'event_site_name']  == "PRA - Priest Rapids Adult") & !(tag_hist[j+1, 'antenna_id'] %in% c(PRA_110_east, PRA_110_AFF)) | # If it's at a different site
                tag_hist[j+1, 'event_date_time_value'] - 
                tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's over 48 hours later
              
              # store the info
              # store the tag code
              ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
              
              # store the event type name
              ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
              
              
              # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
              ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
              
              # store the location fields
              # ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
              # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
              # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
              # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
              # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
              
              # store the end time
              ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
              # Store the end antenna
              ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
              # store the end antenna ID
              ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
              
              # Save event site name
              ind_det_hist[counter,'event_site_name'] <- "PRA2 - Priest Rapids Adult East Ladder"
              
              
              
              # Now, if it was only seen in the AFF but not the exit coils, call it aborted
              if (tag_hist[j, 'antenna_id'] %in% PRA_110_AFF){
                
                ind_det_hist[counter,'non_ascent'] <- "aborted"
                
              }
              
              
              # UPDATE THE COUNTER
              # every time we store an end time, we update the counter. This allows
              # us to move through the detection history df
              counter <- counter + 1
              
              
            }
            
            
          }
          
          
        }
        
      }
      
      
      
      
      
      
      
      
      ##### Rock Island #####
      
      # Rock Island has there separate routes - a left, middle, and right ladder
      # In each of these routes, antennas are only by the exits, in each configuration,
      # but the names of the specific antennas have changed. However, if we just concatenate
      # the antennas for each configuration, we should be fine.
      
      # All antennas are only at the exit, so we cannot distinguish direction in the ladders or aborted attempts.
      
      else if (tag_hist[j, 'event_site_name'] %in% RIS_arrays){
        
        # For RIA, we are not going to worry about configurations, because they have not changed much.
        
        # RIA route 1 - the left bank ladder
        if(tag_hist[j, 'antenna_id'] %in% RIA_left){
          
          # If it's the first detection in this route, or it hasn't been seen at this route in at least 48 hours store the start time and antenna
          if (!(tag_hist[j-1, 'antenna_id'] %in% RIA_left & tag_hist[j-1, 'event_site_name'] == "RIA - Rock Island Adult") |
              tag_hist[j-1, 'antenna_id'] %in% RIA_left & 
              tag_hist[j-1, 'event_site_name'] == "RIA - Rock Island Adult" &
              tag_hist[j, 'event_date_time_value'] - 
              tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)) {
            # Store the event site name (RIA1 - Rock Island Adult Left Ladder)
            ind_det_hist[counter, 'event_site_name'] <- 'RIA1 - Rock Island Adult Left Ladder'
            # Store the start time
            ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
            # Store the start antenna
            ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
            # Store the start antenna id
            ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
          }
          
          
          
          # If the next detection is at a different site, or is more than 48 hours later, store the current time as the end time
          if (!(tag_hist[j+1, 'antenna_id'] %in% RIA_left) | # If it's at a different site
              tag_hist[j+1, 'event_date_time_value'] - 
              tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's more than 48 hours later
            
            # store the info
            # store the tag code
            ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
            
            # store the event type name
            ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
            
            
            # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
            ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
            
            # store the location fields
            ind_det_hist[counter,'event_site_name'] <- "RIA1 - Rock Island Adult Left Ladder"
            # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
            # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
            # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
            # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
            
            # store the end time
            ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
            # Store the end antenna
            ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
            # store the end antenna ID
            ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
            
            # With only four antennas at the exit here, we don't have the ability to identify any non-ascent movements
            
            # UPDATE THE COUNTER
            # every time we store an end time, we update the counter. This allows
            # us to move through the detection history df
            counter <- counter + 1
            
            
          }
          
          # if it's within 48 hours at RIA left, then don't record it
          
          else {
            
          }
          
          
        } 
        # RIA route 2 - the middle bank ladder
        else if(tag_hist[j, 'antenna_id'] %in% RIA_middle){
          
          # If it's the first detection in this route, or it hasn't been seen at this route in at least 48 hours store the start time and antenna
          if (!(tag_hist[j-1, 'antenna_id'] %in% RIA_middle & tag_hist[j-1, 'event_site_name'] == "RIA - Rock Island Adult") |
              tag_hist[j-1, 'antenna_id'] %in% RIA_middle & 
              tag_hist[j-1, 'event_site_name'] == "RIA - Rock Island Adult" &
              tag_hist[j, 'event_date_time_value'] - 
              tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)) {
            # Store the event site name
            ind_det_hist[counter, 'event_site_name'] <- "RIA2 - Rock Island Adult Middle Ladder"
            # Store the start time
            ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
            # Store the start antenna
            ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
            # Store the start antenna id
            ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
          }
          
          
          
          # If the next detection is at a different site, or is more than 48 hours later, store the current time as the end time
          if (!(tag_hist[j+1, 'antenna_id'] %in% RIA_middle) | # If it's at a different site
              tag_hist[j+1, 'event_date_time_value'] - 
              tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's more than 48 hours later
            
            # store the info
            # store the tag code
            ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
            
            # store the event type name
            ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
            
            
            # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
            ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
            
            # store the location fields
            ind_det_hist[counter,'event_site_name'] <- "RIA2 - Rock Island Adult Middle Ladder"
            # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
            # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
            # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
            # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
            
            # store the end time
            ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
            # Store the end antenna
            ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
            # store the end antenna ID
            ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
            
            # With only four antennas at the exit here, we don't have the ability to identify any non-ascent movements
            
            # UPDATE THE COUNTER
            # every time we store an end time, we update the counter. This allows
            # us to move through the detection history df
            counter <- counter + 1
            
            
          }
          
          # if it's within 48 hours at RIA middle, then don't record it
          
          else {
            
          }
          
          
        } 
        
        # RIA route 3 - the right bank ladder
        
        else if(tag_hist[j, 'antenna_id'] %in% RIA_right){
          
          # If it's the first detection in this route, or it hasn't been seen at this route in at least 48 hours store the start time and antenna
          if (!(tag_hist[j-1, 'antenna_id'] %in% RIA_right & tag_hist[j-1, 'event_site_name'] == "RIA - Rock Island Adult") |
              tag_hist[j-1, 'antenna_id'] %in% RIA_right & 
              tag_hist[j-1, 'event_site_name'] == "RIA - Rock Island Adult" &
              tag_hist[j, 'event_date_time_value'] - 
              tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)) {
            # Store the event site name
            ind_det_hist[counter, 'event_site_name'] <- "RIA3 - Rock Island Adult Right Ladder"
            # Store the start time
            ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
            # Store the start antenna
            ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
            # Store the start antenna id
            ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
          }
          
          
          
          # If the next detection is at a different site, or is more than 48 hours later, store the current time as the end time
          if (!(tag_hist[j+1, 'antenna_id'] %in% RIA_right) | # If it's at a different site
              tag_hist[j+1, 'event_date_time_value'] - 
              tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's more than 48 hours later
            
            # store the info
            # store the tag code
            ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
            
            # store the event type name
            ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
            
            
            # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
            ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
            
            # store the location fields
            ind_det_hist[counter,'event_site_name'] <- "RIA3 - Rock Island Adult Right Ladder"
            # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
            # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
            # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
            # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
            
            # store the end time
            ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
            # Store the end antenna
            ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
            # store the end antenna ID
            ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
            
            # With only four antennas at the exit here, we don't have the ability to identify any non-ascent movements
            
            # UPDATE THE COUNTER
            # every time we store an end time, we update the counter. This allows
            # us to move through the detection history df
            counter <- counter + 1
            
            
          }
          
          # if it's within 48 hours at RIA right, then don't record it
          
          else {
            
          }
          
          
        } 
        
      }
      
      ##### Rocky Reach Dam #####
      
      # Rocky Reach Dam has only one ladder, with a set of 4 antennas (a fifth was added in April 2014, config 110) at the exit.
      # Given this, we will not be able to identify aborted ascension attempts or descent through the ladder.
      
      # All antennas are only at the exit, so we cannot distinguish direction in the ladders or aborted attempts.
      
      else if (tag_hist[j, 'event_site_name'] %in% RRE_arrays){
        
        # For RRE, we are not going to do anything special, except that we are going to have a 48 hour window.
        
        # If it's the first detection at RRE, or it hasn't been seen at this route in at least 48 hours store the start time and antenna
        if (!(tag_hist[j-1, 'event_site_name'] == "RRF - Rocky Reach Fishway") |
            tag_hist[j-1, 'event_site_name'] == "RRF - Rocky Reach Fishway" &
            tag_hist[j, 'event_date_time_value'] - 
            tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)) {
          # Store the event site name
          # ind_det_hist[counter, 'event_site_name'] <- "RIA3 - Rock Island Adult Right Ladder"
          ind_det_hist[counter, 'event_site_name'] <- tag_hist[j,'event_site_name']
          # Store the start time
          ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
          # Store the start antenna
          ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
          # Store the start antenna id
          ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
        }
        
        
        
        # If the next detection is at a different site, or is more than 48 hours later, store the current time as the end time
        if (!(tag_hist[j+1, 'event_site_name'] == "RRF - Rocky Reach Fishway") | # If it's at a different site
            tag_hist[j+1, 'event_date_time_value'] - 
            tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's more than 48 hours later
          
          # store the info
          # store the tag code
          ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
          
          # store the event type name
          ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
          
          
          # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
          ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
          
          # store the location fields
          ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
          # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
          # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
          # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
          # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
          
          # store the end time
          ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
          # Store the end antenna
          ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
          # store the end antenna ID
          ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
          
          # With only four antennas at the exit here, we don't have the ability to identify any non-ascent movements
          
          # UPDATE THE COUNTER
          # every time we store an end time, we update the counter. This allows
          # us to move through the detection history df
          counter <- counter + 1
          
          
        }
        
        # if it's within 48 hours at RRF, then don't record it
        
        else {
          
        }
        
        
      }
      
      
      
      ##### Wells Dam #####
      
      # Wells Dam has a left and a right ladder. Each of these has a trap, with a note
      # that "Trap fish are removed to the hatchery or trucked off site"
      
      # If the last antenna a fish is seen at at Wells is a trap antenna, then we
      # note that it was trapped in the non-ascent field
      # This occurs earlier in the script, in the j == nrow(tag_hist)
      
      # Wells has gone through many configurations; 110 and 120 have detections near the exit and in the traps.
      # 130 additionally has antennas in the lower part of the right ladder;
      # 140 and onwards have antennas in the lower parts of the left and right ladders
      
      else if (tag_hist[j, 'event_site_name'] %in% WEL_arrays){
        
        # Wells route 1: Left ladder
        if(tag_hist[j, 'antenna_id'] %in% WEL_left){
          
          # If it's the first detection in this route, or it hasn't been seen at this route in at least 48 hours store the start time and antenna
          if (!(tag_hist[j-1, 'antenna_id'] %in% WEL_left & tag_hist[j-1, 'event_site_name'] == "WEA - Wells Dam, DCPUD Adult Ladders") |
              tag_hist[j-1, 'antenna_id'] %in% WEL_left & 
              tag_hist[j-1, 'event_site_name'] == "WEA - Wells Dam, DCPUD Adult Ladders" &
              tag_hist[j, 'event_date_time_value'] - 
              tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)) {
            # Store the event site name
            ind_det_hist[counter, 'event_site_name'] <- 'WEA1 - Wells Dam, DCPUD Left Adult Ladder'  # tag_hist[j,'event_site_name']
            # Store the start time
            ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
            # Store the start antenna
            ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
            # Store the start antenna id
            ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
          }
          
          # If it's not the first detection in this route
          
          # config 110, 120, and 130 are functionally the same for the left ladder - antennas in the trap and in the upper ladder
          # We don't have to deal with the trap sites here because they're dealt with earlier in the script
          else if (tag_hist[j, "ant_config"] %in% c(110, 120, 130)){
            # If the next detection is at a different site, or is more than 48 hours later, store the current time as the end time
            if (!(tag_hist[j+1, 'antenna_id'] %in% WEL_left) | # If it's at a different site
                tag_hist[j+1, 'event_date_time_value'] - 
                tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's more than 48 hours later
              
              # store the info
              # store the tag code
              ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
              
              # store the event type name
              ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
              
              
              # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
              ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
              
              # store the location fields
              ind_det_hist[counter,'event_site_name'] <- 'WEA1 - Wells Dam, DCPUD Left Adult Ladder' # tag_hist[j,'event_site_name']
              # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
              # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
              # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
              # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
              
              # store the end time
              ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
              # Store the end antenna
              ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
              # store the end antenna ID
              ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
              
              # With only four antennas at the exit here, we don't have the ability to identify any non-ascent movements
              
              # UPDATE THE COUNTER
              # every time we store an end time, we update the counter. This allows
              # us to move through the detection history df
              counter <- counter + 1
              
              
            }
            
            # if it's within 48 hours at PRA west, then don't record it
            
            else {
              
            }
            
            
          }
          
          # antenna configs 140, 150, 160, 170, and 180 all have antennas in the lower ladder
          # they're functionally the same
          else if (tag_hist[j, "ant_config"] %in% c(140, 150, 160, 170, 180)){
            # If the next detection is at a different site, or is more than 48 hours later, store the current time as the end time
            if (!(tag_hist[j+1, 'antenna_id'] %in% WEL_left) | # If it's at a different site
                tag_hist[j+1, 'event_date_time_value'] - 
                tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's more than 48 hours later
              
              # store the info
              # store the tag code
              ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
              
              # store the event type name
              ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
              
              
              # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
              ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
              
              # store the location fields
              ind_det_hist[counter,'event_site_name'] <- 'WEA1 - Wells Dam, DCPUD Left Adult Ladder' # tag_hist[j,'event_site_name']
              # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
              # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
              # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
              # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
              
              # store the end time
              ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
              # Store the end antenna
              ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
              # store the end antenna ID
              ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
              
              # Now, if it was only seen in the lower ladder but not the exit coils, (and wasn't first seen in the exit coils) call it aborted
              if (tag_hist[j, 'antenna_id'] %in% c(WEA_140_left_lower, WEA_150_left_lower,
                                                   WEA_160_left_lower, WEA_170_left_lower,
                                                   WEA_180_left_lower) &
                  !(ind_det_hist[counter,'start_antenna_id'] %in% c(WEA_170_left_upper))){ # the names of these coils didn't change, so just need one config
                
                ind_det_hist[counter,'non_ascent'] <- "aborted"
              }
              
              # Other alternative - if the first detection was at an exit coil, and the last at an entrance coil, then it went downstream through the ladder
              else if (tag_hist[j, 'antenna_id'] %in% c(WEA_140_left_lower, WEA_150_left_lower,
                                                        WEA_160_left_lower, WEA_170_left_lower,
                                                        WEA_180_left_lower) & 
                       ind_det_hist[counter,'start_antenna_id'] %in% WEA_170_left_upper){
                ind_det_hist[counter,'non_ascent'] <- "descent"
              }
              
              # UPDATE THE COUNTER
              # every time we store an end time, we update the counter. This allows
              # us to move through the detection history df
              counter <- counter + 1
              
              
            }
            
            # if it's within 48 hours at PRA west, then don't record it
            
            else {
              
            }
            
            
          }
          
          
          
        } 
        # WEL route 2 - the right ladder
        else if(tag_hist[j, 'antenna_id'] %in% WEL_right){
          
          # If it's the first detection in this route, or it hasn't been seen at this route in at least 48 hours store the start time and antenna
          if (!(tag_hist[j-1, 'antenna_id'] %in% WEL_right & tag_hist[j-1, 'event_site_name'] == "WEA - Wells Dam, DCPUD Adult Ladders") |
              tag_hist[j-1, 'antenna_id'] %in% WEL_right & 
              tag_hist[j-1, 'event_site_name'] == "WEA - Wells Dam, DCPUD Adult Ladders" &
              tag_hist[j, 'event_date_time_value'] - 
              tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)) {
            # Store the event site name
            ind_det_hist[counter, 'event_site_name'] <- 'WEA2 - Wells Dam, DCPUD Right Adult Ladder' # tag_hist[j,'event_site_name']
            # Store the start time
            ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
            # Store the start antenna
            ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
            # Store the start antenna id
            ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
          }
          
          # If it's not the first detection in this route
          
          # config 110 and 120 are functionally the same for the right ladder - antennas in the trap and in the upper ladder
          # We don't have to deal with the trap sites here because they're dealt with earlier in the script
          else if (tag_hist[j, "ant_config"] %in% c(110, 120)){
            # If the next detection is at a different site, or is more than 48 hours later, store the current time as the end time
            if (!(tag_hist[j+1, 'antenna_id'] %in% WEL_right) | # If it's at a different site
                tag_hist[j+1, 'event_date_time_value'] - 
                tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's more than 48 hours later
              
              # store the info
              # store the tag code
              ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
              
              # store the event type name
              ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
              
              
              # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
              ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
              
              # store the location fields
              ind_det_hist[counter,'event_site_name'] <- 'WEA2 - Wells Dam, DCPUD Right Adult Ladder' # tag_hist[j,'event_site_name']
              # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
              # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
              # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
              # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
              
              # store the end time
              ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
              # Store the end antenna
              ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
              # store the end antenna ID
              ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
              
              # With only four antennas at the exit here, we don't have the ability to identify any non-ascent movements
              
              # UPDATE THE COUNTER
              # every time we store an end time, we update the counter. This allows
              # us to move through the detection history df
              counter <- counter + 1
              
              
            }
            
            # if it's within 48 hours at PRA west, then don't record it
            
            else {
              
            }
            
            
          }
          
          # antenna configs 130, 140, 150, 160, 170, and 180 all have antennas in the lower ladder
          # they're functionally the same
          else if (tag_hist[j, "ant_config"] %in% c(130, 140, 150, 160, 170, 180)){
            # If the next detection is at a different site, or is more than 48 hours later, store the current time as the end time
            if (!(tag_hist[j+1, 'antenna_id'] %in% WEL_right) | # If it's at a different site
                tag_hist[j+1, 'event_date_time_value'] - 
                tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's more than 48 hours later
              
              # store the info
              # store the tag code
              ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
              
              # store the event type name
              ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
              
              
              # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
              ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
              
              # store the location fields
              ind_det_hist[counter,'event_site_name'] <- 'WEA2 - Wells Dam, DCPUD Right Adult Ladder' # tag_hist[j,'event_site_name']
              # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
              # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
              # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
              # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
              
              # store the end time
              ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
              # Store the end antenna
              ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
              # store the end antenna ID
              ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
              
              # Now, if it was only seen in the lower ladder but not the exit coils, (and wasn't first seen in the exit coils) call it aborted
              if (tag_hist[j, 'antenna_id'] %in% c(WEA_130_right_lower, WEA_140_right_lower, WEA_150_right_lower,
                                                   WEA_160_right_lower, WEA_170_right_lower, WEA_180_right_lower) &
                  !(ind_det_hist[counter,'start_antenna_id'] %in% c(WEA_170_right_upper))){ # the names of these coils didn't change, so just need one config
                
                ind_det_hist[counter,'non_ascent'] <- "aborted"
              }
              
              # Other alternative - if the first detection was at an exit coil, and the last at an entrance coil, then it went downstream through the ladder
              else if (tag_hist[j, 'antenna_id'] %in% c(WEA_130_right_lower, WEA_140_right_lower, WEA_150_right_lower,
                                                        WEA_160_right_lower, WEA_170_right_lower, WEA_180_right_lower) & 
                       ind_det_hist[counter,'start_antenna_id'] %in% WEA_170_right_upper){
                ind_det_hist[counter,'non_ascent'] <- "descent"
              }
              
              # UPDATE THE COUNTER
              # every time we store an end time, we update the counter. This allows
              # us to move through the detection history df
              counter <- counter + 1
              
              
            }
            
            # if it's within 48 hours at PRA west, then don't record it
            
            else {
              
            }
            
            
          }
          
          
          
        } 
        
      }
      
      
      
      
      
      
      
      
      
      
      
      ##### Ice Harbor Dam #####
      
      # Ice Harbor dam has a north and a south shore ladder. 
      # Two configurations, 100 and 110. 110 only added a trap entrance detector in the south shore, which doesn't really matter.
      
      # ICH (combined) originally included the juvenile bypass as well. I have made an edit
      # earlier in the script to rename all of these as the bypass, so they're not dealt with here
      
      # In each case, we only have antennas at the exit, so we can't monitor non-ascensions
      
      else if (tag_hist[j, 'event_site_name'] %in% ICH_arrays){
        
        # ICH route 1 - the north shore ladder
        # In this route, the configuration has not changed through time, so no ifelse based on configuration
        
        if(tag_hist[j, 'antenna_id'] %in% ICH_110_north_ladder){
          
          # If it's the first detection in this route, or it hasn't been seen at this route in at least 48 hours store the start time and antenna
          if (!(tag_hist[j-1, 'antenna_id'] %in% ICH_110_north_ladder & tag_hist[j-1, 'event_site_name'] == "ICH - Ice Harbor Dam (Combined)") |
              tag_hist[j-1, 'antenna_id'] %in% ICH_110_north_ladder & 
              tag_hist[j-1, 'event_site_name'] == "ICH - Ice Harbor Dam (Combined)" &
              tag_hist[j, 'event_date_time_value'] - 
              tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)) {
            # Store the event site name
            ind_det_hist[counter, 'event_site_name'] <- "ICH1 - Ice Harbor Dam North Ladder"
            # Store the start time
            ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
            # Store the start antenna
            ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
            # Store the start antenna id
            ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
          }
          
          
          
          # If the next detection is at a different site, or is more than 48 hours later, store the current time as the end time
          if (!(tag_hist[j+1, 'antenna_id'] %in% ICH_110_north_ladder) | # If it's at a different site
              tag_hist[j+1, 'event_date_time_value'] - 
              tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's more than 48 hours later
            
            # store the info
            # store the tag code
            ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
            
            # store the event type name
            ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
            
            
            # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
            ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
            
            # store the location fields
            ind_det_hist[counter,'event_site_name'] <- "ICH1 - Ice Harbor Dam North Ladder"
            # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
            # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
            # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
            # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
            
            # store the end time
            ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
            # Store the end antenna
            ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
            # store the end antenna ID
            ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
            
            # With only four antennas at the exit here, we don't have the ability to identify any non-ascent movements
            
            # UPDATE THE COUNTER
            # every time we store an end time, we update the counter. This allows
            # us to move through the detection history df
            counter <- counter + 1
            
            
          }
          
          # if it's within 48 hours at ICH1 north, then don't record it
          
          else {
            
          }
          
          
        } 
        # ICH route 2 - the south ladder
        else if(tag_hist[j, 'antenna_id'] %in% c(ICH_110_south_ladder, ICH_110_south_trap)){
          
          
          # If it's the first detection in this route in at least 48 hours,
          if (!(tag_hist[j-1, 'antenna_id'] %in% c(ICH_110_south_ladder, ICH_110_south_trap) & tag_hist[j-1, 'event_site_name'] == "ICH - Ice Harbor Dam (Combined)") |
              tag_hist[j-1, 'antenna_id'] %in% c(ICH_110_south_ladder, ICH_110_south_trap) &
              tag_hist[j-1, 'event_site_name'] == "ICH - Ice Harbor Dam (Combined)" &
              tag_hist[j, 'event_date_time_value'] - 
              tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)) {
            # Store the event site name
            ind_det_hist[counter, 'event_site_name'] <- "ICH2 - Ice Harbor Dam South Ladder"
            # Store the start time
            ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
            # Store the start antenna
            ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
            # Store the start antenna id
            ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
          }
          
          
          # If it's not the first detection in this route
          
          # If the next detection is not at ICH2 south
          # or is more than 48 hours later, store the current time as the end time
          if (!(tag_hist[j+1, 'event_site_name']  == "ICH - Ice Harbor Dam (Combined)") & !(tag_hist[j+1, 'antenna_id'] %in% c(ICH_110_south_ladder, ICH_110_south_trap)) | # If it's at a different site
              tag_hist[j+1, 'event_date_time_value'] - 
              tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's over 48 hours later
            
            # store the info
            # store the tag code
            ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
            
            # store the event type name
            ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
            
            
            # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
            ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
            
            # store the location fields
            ind_det_hist[counter,'event_site_name'] <- "ICH2 - Ice Harbor Dam South Ladder"
            # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
            # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
            # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
            # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
            
            # store the end time
            ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
            # Store the end antenna
            ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
            # store the end antenna ID
            ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
            
            # With only four antennas at the exit here, we don't have the ability to identify any non-ascent movements
            
            # UPDATE THE COUNTER
            # every time we store an end time, we update the counter. This allows
            # us to move through the detection history df
            counter <- counter + 1
            
            
          }
          
          
          
        }
        
      }
      
      
      
      
      
      
      
      
      
      
      
      ##### Lower Granite Dam #####
      
      # Lower Granite Dam has a single ladder, with arrays in the bottom of the ladder before the AFF,
      # and arrays near the ladder exit.
      
      # Also note that there are recapture/release events in the ladder, where event_type_name == recapture &
      # event_site_name == LGRLDR - LGR - Release into the Adult Fish Ladder
      
      # It has two configurations, 140 and 150. 150 added additional arrays in
      # the entrance and exit of the ladder.
      
      # Based on a quick inspection, I don't really see any evidence of aborted ascension
      # attempts or descents through the ladder. But the code will still allow for it.
      # detection probability at the entrance seems far from 100%. I think we'll just group 
      # the entrance antennas with the lower ones and the exit antennas with the upper ones
      
      else if(tag_hist[j, 'event_site_name'] %in% LGR_arrays){
        
        # If it's the first detection in this route, or if the last detection was in
        # the upper antennas at least 48 hours ago
        if (!(tag_hist[j-1, 'event_site_name'] %in% LGR_arrays) |
            tag_hist[j-1, 'event_site_name'] %in% LGR_arrays &
            # tag_hist[j-1, 'antenna_id'] %in% c(GRA_150_exit, GRA_150_upper) &
            tag_hist[j, 'event_date_time_value'] - 
            tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)) {
          # Store the event site name
          ind_det_hist[counter, 'event_site_name'] <- tag_hist[j,'event_site_name']
          # Store the start time
          ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
          # Store the start antenna group
          ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
          # Store the start antenna ID
          ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
        }
        
        
        # Different if statements for varying antenna configurations
        
        # GRA 140 - upper and lower
        else if(tag_hist[j, 'ant_config'] == 140 ){
          
          
          # If the next detection is at a different site, or is more than 48 hours later
          # and is in the lower ladder, store the current time as the end time
          if (!(tag_hist[j+1, 'event_site_name'] %in% LGR_arrays) | # If it's at a different site
              tag_hist[j+1, 'antenna_id'] %in% c(GRA_150_entrance, GRA_150_lower) &
              tag_hist[j+1, 'event_date_time_value'] - 
              tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's within 48 hours
            
            # store the info
            # store the tag code
            ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
            
            # store the event type name
            ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
            
            
            # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
            ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
            
            # store the location fields
            ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
            # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
            # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
            # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
            # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
            
            # store the end time
            ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
            # Store the end antenna
            ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
            # store the end antenna ID
            ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
            
            
            
            # However, if the last antenna seen is not an exit coil, note that it was an aborted attempt
            if (tag_hist[j, 'antenna_id'] %in% GRA_140_lower & !(ind_det_hist[counter,'start_antenna_id'] %in% GRA_140_upper)){
              
              ind_det_hist[counter,'non_ascent'] <- "aborted"
              
            }
            
            # Other alternative - if the first detection was at an exit coil, and the last at an entrance coil, then it went downstream through the ladder
            else if (tag_hist[j, 'antenna_id'] %in% GRA_140_lower & ind_det_hist[counter,'start_antenna_id'] %in% GRA_140_upper){
              ind_det_hist[counter,'non_ascent'] <- "descent"
            }
            
            # UPDATE THE COUNTER
            # every time we store an end time, we update the counter. This allows
            # us to move through the detection history df
            counter <- counter + 1
            
            
          }
          
          # if it's within 48 hours at GRA, then don't record it
          
          else {
            
          }
          
        }
        
        
        # GRA 150 - upper, lower, entrance, exit
        # Include the recapture possibility here
        else if(tag_hist[j, 'ant_config'] == 150 | tag_hist[j, 'event_site_name'] == "LGRLDR - LGR - Release into the Adult Fish Ladder"){
          
          
          # If the next detection is at a different site, or is more than 48 hours later
          # and is in the lower ladder, store the current time as the end time
          if (!(tag_hist[j+1, 'event_site_name'] %in% LGR_arrays) | # If it's at a different site
              tag_hist[j+1, 'antenna_id'] %in% c(GRA_150_entrance, GRA_150_lower) &
              tag_hist[j+1, 'event_date_time_value'] - 
              tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's within 48 hours
            
            # store the info
            # store the tag code
            ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
            
            # store the event type name
            ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
            
            
            # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
            ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
            
            # store the location fields
            ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
            # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
            # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
            # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
            # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
            
            # store the end time
            ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
            # Store the end antenna
            ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
            # store the end antenna ID
            ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
            
            
            # However, if the last antenna seen is not an exit coil, note that it was an aborted attempt
            if (tag_hist[j, 'antenna_id'] %in% c(GRA_150_lower, GRA_150_entrance) & !(ind_det_hist[counter,'start_antenna_id'] %in% c(GRA_150_upper, GRA_150_exit))){
              
              ind_det_hist[counter,'non_ascent'] <- "aborted"
              
            }
            
            # Other alternative - if the first detection was at an exit coil, and the last at an entrance coil, then it went downstream through the ladder
            else if (tag_hist[j, 'antenna_id'] %in% c(GRA_150_lower, GRA_150_entrance) & ind_det_hist[counter,'start_antenna_id'] %in% c(GRA_150_upper, GRA_150_exit)){
              ind_det_hist[counter,'non_ascent'] <- "descent"
            }
            
            # UPDATE THE COUNTER
            # every time we store an end time, we update the counter. This allows
            # us to move through the detection history df
            counter <- counter + 1
            
            
          }
          
          # if it's within 48 hours at GRA, then don't record it
          
          else {
            
          }
          
        }
        
      } 
      
      ##### non-dams #####   
      
      # For every other entry that is not at a dam, look at the previous entry to see if it
      # was the same site < 48 hours ago
      else {
        
        # If the next entry isn't the same site OR is >=48 hours ahead, store
        # the current time as the end time
        if (tag_hist[j+1, 'event_site_name'] != tag_hist[j, 'event_site_name'] |
            tag_hist[j+1, 'event_date_time_value'] -
            tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){
          
          # store the end time
          ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
          # Store the end antenna
          ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
          # store the end antenna ID
          ind_det_hist[counter,'end_antenna_id'] <- tag_hist[j,'antenna_id']
          
          # SPECIAL CASE: If there is only one detection at a site, store the 
          # start time as well
          # Here we will look at both the previous and next entry
          # If the previous detection is at a different site OR is more than 6 hours ago,
          # store the start time
          if (tag_hist[j-1, 'event_site_name'] != tag_hist[j, 'event_site_name'] |
              tag_hist[j, 'event_date_time_value'] -
              tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)){
            
            # store the tag code
            ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
            
            # store the event type name
            ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
            
            
            # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
            ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
            
            # store the location fields
            ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
            
            # Store the start time
            ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
            # Store the start antenna
            ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
            # store the end antenna ID
            ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
            
          }
          
          # UPDATE THE COUNTER
          # every time we store an end time, we update the counter. This allows
          # us to move through the detection history df
          counter <- counter + 1
          
        }
        
        # If the previous site was a different site OR was >6 hours ago, 
        # start a new entry and start time
        else if (tag_hist[j-1, 'event_site_name'] != tag_hist[j, 'event_site_name'] |
                 tag_hist[j, 'event_date_time_value'] -
                 tag_hist[j-1, 'event_date_time_value'] >= hours(x = 48)){
          
          
          # store the tag code
          ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
          
          # store the event type name
          ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
          
          
          # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
          ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
          
          # store the location fields
          ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
          
          # store the start time
          ind_det_hist[counter,'start_time'] <- tag_hist[j,'event_date_time_value']
          # Store the start antenna
          ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
          # store the start antenna ID
          ind_det_hist[counter,'start_antenna_id'] <- tag_hist[j,'antenna_id']
          
        }
        
        # If the previous entry was the same site <=6 hours ago AND
        # The next entry is the same site <=6 hours ago, skip it
        else {
          
        }
        
      }
      
      
    }
    
    # Cut out extra rows
    ind_det_hist <- ind_det_hist[!(is.na(ind_det_hist$tag_code)),]
    
    # Append the individual tag history to the complete tag history
    det_hist %>% 
      bind_rows(., ind_det_hist) -> det_hist
  }
  # Print the run time
  end_time <- Sys.time()
  
  print(paste0("Total tags: ", length(unique_tag_IDs)))
  print(paste0("Total time: ", end_time - start_time))
  print("For loop finished")
  # Change the timezone back to PST, with a correction for DST
  det_hist %>% 
    mutate(start_time = force_tz(start_time, "PST") - dhours(1)) %>% 
    mutate(end_time = force_tz(end_time, "PST") - dhours(1)) -> det_hist_PST
  
  # Create a data frame that contains the metadata for each site
  
  det_hist %>% 
    # Make a correction for JDA - shouldn't be necessary later once the script is re-run
    mutate(event_site_name = ifelse(event_site_name %in% c("JO1 - John Day South Fish Ladder", "JO2 - John Day North Fish Ladder"),
                                    "John Day Dam Adult Fishways (combined)", event_site_name)) %>% 
    # mutate(event_site_latitude = ifelse(event_site_name == "John Day Dam Adult Fishways (combined)", 45.71866, event_site_latitude)) %>% 
    # mutate(event_site_longitude = ifelse(event_site_name == "John Day Dam Adult Fishways (combined)", -120.6978, event_site_longitude)) %>% 
    group_by(event_site_name) %>% 
    summarise(n()) %>% 
    dplyr::rename(count = `n()`)-> event_det_counts
  
  # Get the metadata
  det_hist %>% 
    dplyr::select(-c(tag_code, start_time, end_time)) %>% 
    distinct(event_site_name, .keep_all = TRUE) -> event_site_metadata
  
  output <- list(det_hist = det_hist_PST, event_site_metadata = event_site_metadata)
  return(output)
}


