#### 03: create state histories

# load libraries
library(tidyverse)
library(lubridate)
library(here)

#### Step 1: Run detection histories through functions/create_state_history.R ####
# This step can either be run locally or through a HPC

# Did you run the scripts in R/hyak/03_create_state_histories ?
HPC <- TRUE

if (HPC == TRUE) {
  print("Please consult the readme file and run the jobs found in R/hyak to generate the following files")
  states_complete_part1 <- read.csv("intermediate_outputs/states_complete_part1.csv")
  states_complete_part2 <- read.csv("intermediate_outputs/states_complete_part2.csv")
  states_complete_part3 <- read.csv("intermediate_outputs/states_complete_part3.csv")
  states_complete_part4 <- read.csv("intermediate_outputs/states_complete_part4.csv")
} else {
  # source function
  source("R/functions/create_state_history.R")
  
  states_complete_part1 <- create_state_history(data_quartile = 1)
  states_complete_part2 <- create_state_history(data_quartile = 2)
  states_complete_part3 <- create_state_history(data_quartile = 3)
  states_complete_part4 <- create_state_history(data_quartile = 4)
}

#### Step 2: Join outputs and clean up state histories

states_complete_part1 %>% 
  bind_rows(., states_complete_part2) %>% 
  bind_rows(., states_complete_part3) %>% 
  bind_rows(., states_complete_part4) -> states_complete

# Change all Hood River states to BON to MCN other tributaries - we don't care
# about this tributary individually because we don't have any Hood River fish in our final model
states_complete %>% 
  mutate(state = ifelse(grepl("Hood River", state), 
                        "BON to MCN other tributaries", state)) -> states_complete


# Read in tag code metadata
read.csv(here::here("Data", "covariate_data", "tag_code_metadata.csv")) -> tag_code_metadata

# Fix the release dates
tag_code_metadata %>% 
  mutate(release_date = mdy(release_date_mmddyyyy)) -> tag_code_metadata

# Join release date with the states complete
states_complete %>% 
  left_join(dplyr::select(tag_code_metadata, tag_code, release_date), by = "tag_code") %>% 
  # fix date times: some are missing a time, which we will assume to be 00:00:00
  # Apparently having a time be 00:00:00 causes downstream issues, so make it 00:00:01
  mutate(date_time = ifelse(nchar(date_time) == 10, paste0(date_time, " 00:00:01"), date_time)) %>% 
  mutate(date_time = ymd_hms(date_time)) -> states_complete

# We're going to have to interpolate the times for implicit visits for the code to work

# First, figure out indices of missing date_time
missing_date_time <- is.na(states_complete$date_time)
date_times <- states_complete$date_time

# Use the next time that the fish was detected as the next time
for (i in 1:nrow(states_complete)){
  
  if (is.na(states_complete[i,"date_time"])){
    # Find the next known time
    # Truncate the missing date times to only ones after current state
    missing_date_time_subset <- missing_date_time[(i+1):nrow(states_complete)]
    next_time_index <- min(which(missing_date_time_subset %in% FALSE)) + i
    
    # Now, interpolate the missing time
    # next_time <- ymd_hms(date_times[next_time_index])
    next_time <- format(as.POSIXct(date_times[next_time_index], tz = "UTC"), "%Y-%m-%d %H:%M:%S")
    
    # Get the missing time - add the time difference divided by the number
    # of missing steps plus 1, multiply by which number step it is
    missing_time <- ymd_hms(next_time)
    
    # populate the missing date_time
    states_complete[i, "date_time"] <- missing_time
    
  }
  # If it's not, leave it alone (do nothing)
  else {
    
  }
}



# Put this checkpoint in so we don't have to re-run the for loop
# write.csv(states_complete, here::here("intermediate_outputs", "states_complete_times_interpolated.csv"))

# read.csv(here::here("intermediate_outputs", "states_complete_times_interpolated.csv")) %>%
#   dplyr::select(-X) -> states_complete

# Fix date time
states_complete %>% 
  mutate(date_time = ymd_hms(date_time)) -> states_complete

# Figure out how long detections were after release time
states_complete %>% 
  mutate(time_after_release = date_time - as.POSIXct(release_date)) -> states_complete

# Take apart the date time into y, m, d for subsetting
states_complete %>% 
  mutate(release_year = year(release_date)) %>% 
  mutate(release_month = month(release_date)) %>% 
  mutate(release_day = day(release_date)) %>% 
  mutate(event_year = year(date_time)) %>% 
  mutate(event_month = month(date_time)) %>% 
  mutate(event_day = day(date_time)) -> states_complete

##### Identify juvenile movements #####

# Identify juvenile movements as:
# 1) within 90 days of release; or
# 2) On or before June 15 of the release year, or if the release was on or after July 1 of the previous year, then any movements before June 15 of this year

states_complete %>% 
  mutate(life_stage = "Adult") %>% #first, classify all movement as adult (we will then change this to kelt, repeat spawner, and juvenile later if certain conditions are met)
  mutate(life_stage = ifelse(time_after_release <= days(x = 90) | # fish that were seen right after release
                               # mutate(life_stage = ifelse(
                               # or fish that were seen the same year as release and before June 15
                               release_year == event_year & event_month <= 5 | 
                               release_year == event_year & event_month == 6 & event_day <= 15 |
                               # or fish that were released the year before (on or after July 1) and seen the subsequent year on or before June 15
                               release_year == event_year-1 & release_month >= 7 & event_month <= 5 | 
                               release_year == event_year-1 & release_month >= 7 & event_month == 6 & event_day <= 15, "Juvenile",  
                             life_stage)) -> states_complete

# Remove all detection histories that are only juvenile (keep only those detection histories with at least one adult movement)
states_complete %>% 
  group_by(tag_code) %>% 
  filter(any(life_stage == "Adult")) -> states_complete

# Let's figure out the first ADULT arrival at BON
states_complete %>% 
  group_by(tag_code) %>% 
  subset(pathway %in% c("BON (adult)", "BON_aborted")) %>% 
  subset(life_stage != "Juvenile") %>% 
  filter(date_time == min(date_time)) %>% 
  dplyr::select(tag_code, date_time) %>% 
  dplyr::rename(BON_arrival = date_time) -> BON_arrival_df

# merge that back in
states_complete %>% 
  left_join(., BON_arrival_df, by = "tag_code") -> states_complete

# calculate time after adult arrival
# this line takes quite a while
states_complete %>% 
  mutate(time_after_arrival = date_time - as.POSIXct(BON_arrival)) -> states_complete



# One correction: if an implicit state follows a juvenile movement, it also has to be a juvenile movement. This way, we will remove
# the implicit mouth to BON state that we get when we see a juvenile in the adult ladder and return as an adult. That state
# visit in between shouldn't be modeled
states_complete %>% 
  mutate(life_stage = ifelse(pathway == "implicit" & lag(life_stage) == "Juvenile", "Juvenile", life_stage)) -> states_complete

checkpoint_states_complete <- states_complete

# Another correction: if mouth to BON follows a juvenile movement, it also has to be a juvenile movement. We want all fish starting when they reach BON (adult)
# Edit 2022-08-11: exclude mouth to BON that shows up because of an aborted BON attempt
# Second fix - make sure there are no NAs, so keep all juveniles as juveniles here
states_complete %>% 
  mutate(life_stage = ifelse(life_stage == "Juvenile", "Juvenile", 
                             # Make a fix here, so that any detections prior to arrival date at BON are juvenile. This will fix individuals seen at BHL
                             ifelse(time_after_arrival < 0, "Juvenile",
                                    ifelse(state == "mainstem, mouth to BON" & lag(life_stage) == "Juvenile" & pathway != "BON_aborted", "Juvenile", life_stage)))) -> states_complete



# Let's look for juveniles
states_complete %>% 
  group_by(tag_code) %>% 
  filter(any(life_stage == "Juvenile")) -> juv_movements

# Let's look at the adult histories; do they all start with BON (adult)? They should
states_complete %>% 
  subset(life_stage == "Adult") %>% 
  group_by(tag_code) %>% 
  dplyr::mutate(order_2 = row_number()) %>% 
  subset(order_2 == 1) -> first_adult_states

table(first_adult_states$pathway)
# This looks good! They're all starting either with BON (adult) or BON (aborted)

##### Identify kelt movements #####

# First, we need to identify the end of spawning movements.

# Identify which sites are tributaries and which are mainstem sites
# so we never even use this in this script!
# We jsut say that everything that isn't mainstem is a tributary. Which makes sense. So this isn't affected by division into upstream and river mouth interrogation sites.
tributary_mainstem <- data.frame(tributary = c("Asotin Creek", 
                                               "Clearwater River",
                                               "Deschutes River", 
                                               "Entiat River", 
                                               "Fifteenmile Creek", 
                                               "Grande Ronde River", 
                                               # "Hood River",
                                               "Imnaha River",
                                               "John Day River", 
                                               "Methow River", 
                                               "Okanogan River", 
                                               "Salmon River", 
                                               "Tucannon River", 
                                               "Umatilla River",
                                               "Walla Walla River",
                                               "Wenatchee River", 
                                               "Yakima River",
                                               "Upstream LGR other tributaries",
                                               "ICH to LGR other tributaries",
                                               "BON to MCN other tributaries",
                                               "Upstream WEL other tributaries"
),
mainstem = c("mainstem, upstream of LGR",
             "mainstem, upstream of LGR",
             "mainstem, BON to MCN",
             "mainstem, RRE to WEL",
             "mainstem, BON to MCN",
             "mainstem, upstream of LGR",
             # "mainstem, BON to MCN",
             "mainstem, upstream of LGR",
             "mainstem, BON to MCN",
             "mainstem, upstream of WEL",
             "mainstem, upstream of WEL",
             "mainstem, upstream of LGR",
             "mainstem, ICH to LGR",
             "mainstem, BON to MCN",
             "mainstem, MCN to ICH or PRA",
             "mainstem, RIS to RRE",
             "mainstem, MCN to ICH or PRA",
             "mainstem, upstream of LGR",
             "mainstem, ICH to LGR",
             "mainstem, BON to MCN",
             "mainstem, upstream of WEL"))

# Get the order of sites
# Order of sites (no tributaries), Columbia River
site_order_notrib_columbia <- c("mainstem, mouth to BON", "mainstem, BON to MCN", 
                                "mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS",
                                "mainstem, RIS to RRE", "mainstem, RRE to WEL", 
                                "mainstem, upstream of WEL")

# Order of sites (no tributaries), Snake
# site_order_notrib_snake <- c("mainstem, ICH to LGR", "mainstem, upstream of LGR")
# Contains all non-tributary sites from upstream of LGR to downstream of BON
site_order_notrib_snake <- c("mainstem, mouth to BON", "mainstem, BON to MCN", 
                             "mainstem, MCN to ICH or PRA",
                             "mainstem, ICH to LGR",
                             "mainstem, upstream of LGR")


# Site order, no tributaries, Snake directly to upper Columbia
site_order_notrib_columbia_snake <- c("mainstem, upstream of WEL",
                                      "mainstem, RRE to WEL", 
                                      "mainstem, RIS to RRE",
                                      "mainstem, PRA to RIS",
                                      "mainstem, MCN to ICH or PRA",
                                      "mainstem, ICH to LGR",
                                      "mainstem, upstream of LGR")


# Site order, both:
site_order_combined <- data.frame(state = c("mainstem, mouth to BON", "mainstem, BON to MCN", 
                                            "mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS",
                                            "mainstem, RIS to RRE", "mainstem, RRE to WEL", 
                                            "mainstem, upstream of WEL", "mainstem, ICH to LGR",
                                            "mainstem, upstream of LGR"), state_order = c(1,2,3,4,5,6,7,4,5))


# Final spawning movements:
# Seen at BON adult (repeat spawner); then all downstream movement before that and after last being seen is kelt movement

states_complete %>% 
  group_by(tag_code) %>% 
  mutate(order = row_number()) %>% 
  mutate(max_order = max(order)) -> states_complete

# Determine if movement is occurring upstream or downstream
states_complete %>% 
  left_join(., site_order_combined, by = "state") %>% 
  # now add a line so that every tributary gets a 99, which makes it so that any movement 
  # out of a tributary is considered downstream movement and any movement into a tributary is considered an upstream movement
  mutate(state_order = ifelse(is.na(state_order), 99, state_order)) -> states_complete

# Now classify movements as upstream or downstream
states_complete %>% 
  group_by(tag_code) %>% 
  mutate(direction = ifelse(lag(state_order) > state_order, "downstream", "upstream")) %>% 
  mutate(direction = ifelse(order == 1, "upstream", direction))-> states_complete

states_complete %>% 
  mutate(life_stage = ifelse(tag_code == lag(tag_code) & # if it's the same fish
                               pathway == "BON (adult)" & # if it was seen in the adult ladder at BON
                               time_after_arrival > days(x = 180) &  # and a long time after arrival
                               event_month >= 6 & event_month <= 11 & # and it was seen at BON in a month consistent with adult run timing, based on https://www.cbr.washington.edu/dart/wrapper?type=php&fname=hrt_adult_1657840462_690.php
                               lag(direction) == "downstream", # & it was just seen moving downstream (kelt movement)
                             "repeat",  life_stage)) -> states_complete  # then call it a repeat spawner

# Now, call all downstream movement prior to a return to BON kelt movement

# Make a new field - order of direction, based on tag code and direction
# every time it changes, change direction index
states_complete$direction_index <- c(0,cumsum(as.numeric(with(states_complete,direction[1:(length(direction)-1)] != direction[2:length(direction)]))))

# Now change the downstream movements prior to being a repeat spawner to kelt movement
# First, find the direction indices preceding repeat
states_complete %>% 
  subset(life_stage == "repeat") %>% 
  mutate(repeat_start = order) %>% 
  distinct(tag_code, .keep_all = TRUE)-> repeat_spawners

# Get the direction indices of when individuals are repeat spawners
repeat_spawners$direction_index -> repeat_direction_indices

pre_repeat_movement <- repeat_direction_indices - 1

# Now call the downstream movement prior to repeat spawning kelt movement
states_complete %>% 
  mutate(life_stage = ifelse(direction_index %in% pre_repeat_movement & direction == "downstream", "kelt", life_stage)) -> states_complete

# Also find kelt movement that does not precede a return
# if a fish does not return, it's really hard to determine where it spawned - you can't really say
# for sure when fallback is downstream kelt movement and when it's just trying to get to trib

# steelhead are known to spawn in the spring, so kelts would be moving in spring/summer back downstream

# Let's look first at downstream movements at kelt timing
# states_complete %>% 
#   mutate(life_stage = ifelse(event_month %in% c(5,6,7) & # select months we might expect kelts
#                                direction == "downstream" & # make sure that they're going downstream and (next line) make sure it's a long time after
#                                time_after_arrival > days(180), "kelt", life_stage)) -> states_complete

# so by month doesn't work, because most of the kelt movement is implicit movements
# that aren't actually detected (due to bad detection in downstream movements)
# instead, we'll want to choose based on just time after arrival
# states_complete %>% 
#   mutate(life_stage = ifelse(direction == "downstream" & # make sure that they're going downstream and (next line) make sure it's a long time after
#                                time_after_arrival > days(180), "kelt", life_stage)) -> states_complete

# not going to work, again because of implicit move timing.
# Let's try lag(direction)
# also only allow it if it's been at least 90 days and movement is occurring when we'd expect it to (following spring spawning)
# states_complete_orig <- states_complete

states_complete %>%
  mutate(life_stage = ifelse(order == max_order & event_month >= 3 & event_month <= 7 & time_after_arrival > days(x = 90) & direction == "downstream", "kelt",
                             ifelse(event_month >= 3 & event_month <= 7 & time_after_arrival > days(x = 90) & direction == "downstream" & lag(direction) == "downstream" | # if two consecutive downstream movements, then they're possible kelt movements
                                      event_month >= 3 & event_month <= 7 & time_after_arrival > days(x = 90) & direction == "downstream" & lead(direction) == "downstream", "kelt", 
                                    life_stage))) -> states_complete


# TEMPORARY CALL: Any terminal downstream movement is kelt movement

# This sort of works, but also IDs lots of movement, especially between upper snake and upper columbia or vice versa, as kelt.
# We'll fix up the time interpolation, then subset by only spring/late winter (probably around March to July or something like that)
# Okay, so the problem with subsetting only spring/late winter downstream movement
# is that again, with implicit site movements, we don't have that time.
# SO: I think we need to identify repeat spawners, then identify all consecutive 
# downstream movement prior to it arriving back at BON as kelt movement


# Here we are figuring out what the next direction is - this is important so that we can ID when downstream movement is NOT kelt movement
# e.g., when we have fallback events in the middle of an adult history
states_complete %>% 
  ungroup() %>% 
  mutate(direction_index = direction_index - 1) %>% 
  subset(direction_index >= 0) %>%  # remove the first one
  # mutate(index = row_number()) %>% 
  dplyr::rename(next_direction = direction) %>% 
  dplyr::select(next_direction, direction_index, tag_code) %>% 
  # dplyr::rename(next_direction_tag_code = tag_code) %>% 
  distinct(., next_direction, direction_index) -> direction_index_lead_df

states_complete %>% 
  ungroup() %>% 
  # mutate(index = row_number()) %>% 
  left_join(., direction_index_lead_df, by = c("direction_index")) -> states_complete

# Okay, I think we can use group_by(tag_code) and max(direction index) to see if it still changes direction for the same tag code
# BUT - we need to make sure it's not a repeat spawner. But then theoretically this would exclude fallback by repeat spawners
# So it definitely does exclude fallback by repeat spawners, which is an issue, but it looks like only for tag code 3D9.1C2D7EBE8B, which seems to be the only affected individual
states_complete %>% 
  group_by(tag_code) %>% 
  mutate(max_direction_index = max(direction_index)) -> states_complete

# First, create a new field for if at any point it's a repeat spawner
states_complete %>% 
  mutate(eventual_repeat = ifelse(tag_code %in% repeat_spawners$tag_code, "eventual_repeat", "single")) %>% 
  left_join(., dplyr::select(repeat_spawners, c(tag_code, repeat_start)), by = "tag_code")-> states_complete

# Now, change every movement after repeat start to repeat
states_complete %>% 
  mutate(life_stage = ifelse(eventual_repeat == "eventual_repeat" & order >= repeat_start, "repeat", life_stage)) -> states_complete

# Now change those kelt life stage events that are not actually kelts (because they're followed by upstream movement) back to Adult
states_complete %>% 
  mutate(life_stage = ifelse(life_stage == "kelt" & next_direction == "upstream" & direction_index != max_direction_index & eventual_repeat == "single",
                             "Adult", life_stage)) -> states_complete

## Breakpoint: Let's look into this and see if our kelts make sense
subset(states_complete, life_stage == "kelt")$tag_code -> kelt_tag_codes

subset(states_complete, tag_code %in% kelt_tag_codes) -> current_kelts
## This currently does not look right at all...? Let's comment
# out the fixes for individual fish, run the rest of the code, and see how it looks.

# edit 2022-08-12 - make the fix for tag code 3D9.1C2D7EBE8B, which is one repeat spawner with fallback that's not kelt movement
# states_complete %>% 
  # mutate(life_stage = ifelse(tag_code == "3D9.1C2D7EBE8B" & event_year == 2012 & event_month == 4 & event_day == 11, "Adult", life_stage)) -> states_complete

# edit 2022-12-01 - tag code 3D9.1C2D9303CC is also a repeat spawner with fallback that's not kelt movement
# states_complete %>% 
  # mutate(life_stage = ifelse(tag_code == "3D9.1C2D9303CC" & event_year == 2014 & event_month == 3 & event_day == 14, "Adult", life_stage)) -> states_complete

# states_complete_ckpt <- states_complete



# Let's check on our repeat spawners
states_complete %>% 
  subset(., eventual_repeat == "eventual_repeat") -> eventual_repeat_df


# Let's check on our kelts
states_complete %>% 
  group_by(tag_code) %>% 
  filter(any(life_stage == "kelt")) %>% 
  dplyr::select(tag_code, state, date_time, life_stage, direction, pathway, order, next_direction, max_direction_index) -> kelts


# Check some individuals
# check on our guy
states_complete %>% 
  subset(tag_code == "384.36F2B442DB") -> test_kelt

# check how the triple spawners are doing
subset(eventual_repeat_df, tag_code == "3D9.1C2C430C8D") -> test_triple

# find possible repeats
kelts %>% 
  group_by(tag_code) %>% 
  slice(tail(row_number(),1)) %>% 
  subset(life_stage != "kelt") -> non_terminal_kelts

non_terminal_kelts$tag_code -> non_terminal_kelts_tags

non_terminal_kelts_hist <- subset(kelts, tag_code %in% non_terminal_kelts_tags)


# This guy (384.1B796A8EE1) is weird - one random detection in BON adult. Likely same as with juveniles, where it was seen there but going downstream?
# Need to make it so that all movement following kelt detection is kelt, until we see them again as a repeat spawner
# 384.36F2B32373 also odd - definitely looks like kelt movement, but it would mean that it spawned by March. I guess that's not impossible
# 384.36F2B442DB: This is not kelt movement - some downstream movement in April, but then back upstream to Clearwater. 
# Need to not call things kelt movement if they then go back upstream after
# 	384.36F2B48AEC - same story as above - two downstream movements in the middle, but back upstream after
# So kelt movement either has to be terminal in the detection history, or has to be followed by repeat spawner movement.
# we can make this fix I think
# not a kelt - 384.3B2397826C. Same as before (downstream in the middle, switching from upper columbia to snake)
# 384.3B239A8BD6 - individual where we're not identifying some of the kelt movements. Would be fixed by making the change to call everything prior
# to repeat spawning as kelt movements
# mis IDed fallback as kelt movement :	384.3B239CF5F8

# IDENTIFY REPEAT KELT MOVEMENT

# This guy definitely has kelt movement after being a repeat spawner and returning
# 3D9.1BF18C45D3

# Use the same code again, but only for repeat spawners

states_complete %>% 
  group_by(tag_code) %>% 
  filter(any(is.na(life_stage))) -> NA_life_stages

states_complete %>% 
  mutate(life_stage = ifelse(life_stage == "repeat" & event_month >= 3 & event_month <= 7 & time_after_arrival > days(x = 90) & 
                               direction == "downstream" & order == max_order, "repeat_kelt",
                             ifelse(life_stage == "repeat" & event_month >= 3 & event_month <= 7 & time_after_arrival > days(x = 90) & 
                                      direction == "downstream" & lag(direction) == "downstream" |
                                      life_stage == "repeat" & event_month >= 3 & event_month <= 7 & time_after_arrival > days(x = 90) & 
                                      direction == "downstream" & lead(direction) == "downstream", "repeat_kelt", life_stage))) -> states_complete

# states_complete_ckpt <- states_complete

# Quick fix for individuals (e.g., 3D9.1C2CCDA88C) with an upstream movement that's clearly a kelt (sandwiched between downstream movements)
# This is introducing a lot of NAs because if it's the last line in the history, then there's no lead, and if it's the first, there's no lag
# fixed!
states_complete %>% 
  mutate(life_stage = ifelse(eventual_repeat == "eventual_repeat" & order == max_order | eventual_repeat == "eventual_repeat" & order == 1, life_stage,
                             ifelse(eventual_repeat == "eventual_repeat" & lag(life_stage == "kelt") & lead(life_stage == "kelt"), "kelt", life_stage))) -> states_complete

known_trouble_tags <- c("3D9.1C2D7EBE8B", "3D9.1C2D9303CC")
subset(states_complete, tag_code %in% known_trouble_tags) -> trouble_tags_df
# Ok, these still need to be fixed. We need to handle downstream movements
# after reaching a tributary that are actually not kelt movements, but are
# post-overshoot fallback movements after staging in a non-natal tributary.

# Solution: if there's kelt movement sandwiched
# in adult movements, then they're just adults. Because kelt should always either
# be the last life history stage observed, or be followed by repeat spawner

# step 1: Create a next life stage column
states_complete$life_stage_index <- c(0,cumsum(as.numeric(with(states_complete,life_stage[1:(length(life_stage)-1)] != life_stage[2:length(life_stage)]))))

states_complete %>% 
  ungroup() %>% 
  mutate(life_stage_index = life_stage_index - 1) %>% 
  subset(life_stage_index >= 0) %>%  # remove the first one
  # mutate(index = row_number()) %>% 
  dplyr::rename(next_life_stage = life_stage) %>% 
  dplyr::select(next_life_stage, life_stage_index, tag_code) %>% 
  # dplyr::rename(next_life_stage_tag_code = tag_code) %>% 
  distinct(., next_life_stage, life_stage_index) -> life_stage_index_lead_df

states_complete %>% 
  ungroup() %>% 
  left_join(., life_stage_index_lead_df, by = c("life_stage_index")) -> states_complete

# now, change it so that if they will eventually repeat spawning, they can't have adult/kelt/adult movement
states_complete %>% 
  mutate(life_stage = ifelse(eventual_repeat == "eventual_repeat" &
                               life_stage == "kelt" & 
                               next_life_stage == "Adult", "Adult", life_stage)) -> states_complete

subset(states_complete, tag_code %in% known_trouble_tags) -> test
# great - we fixed these two tag codes!


states_complete %>% 
  group_by(tag_code) %>% 
  filter(any(is.na(life_stage))) -> NA_life_stages
# there aren't any fish in here, which is good!

# Look into our repeat kelts
states_complete %>% 
  group_by(tag_code) %>% 
  filter(any(life_stage == "repeat_kelt")) -> repeat_kelts

# Ok - lets check on all kelts
subset(states_complete, life_stage == "kelt")$tag_code -> kelt_tag_codes
subset(states_complete, tag_code %in% kelt_tag_codes) -> current_kelts
# Looks good!

##### Truncate and export final history, without juveniles and kelts #####

# We want to model adult movement, but not kelt movement, and note when an adult has become a repeat spawner.
# For that, we will create a new field based on tag code, with repeat appended for those individuals

# Edit 2022-08-12 - sometimes if there are repeat repeats, then we need to append a number.
# This might be the only fish that's a repeat repeat: 3D9.1BF27BDC8D
# So we can see if it was a repeat kelt, and if it was, then the next repeat after it is another repeat

states_complete %>% 
  group_by(tag_code) %>% 
  subset(life_stage == "repeat_kelt") %>% 
  mutate(repeat_kelt_end = max(order)) %>% 
  distinct(tag_code, .keep_all = TRUE)-> repeat_kelts


states_complete %>% 
  left_join(., dplyr::select(repeat_kelts, c(tag_code, repeat_kelt_end)), by = "tag_code") %>% 
  # need to make this edit so that we aren't left with NA tag_code_2
  mutate(repeat_kelt_end = ifelse(is.na(repeat_kelt_end), 999, repeat_kelt_end)) %>% 
  mutate(tag_code_2 = ifelse(eventual_repeat == "eventual_repeat" & order > repeat_kelt_end, paste0(tag_code, "_repeat_2"),
                             ifelse(eventual_repeat == "eventual_repeat" & order >= repeat_start, paste0(tag_code, "_repeat"), tag_code))) -> states_complete

# Check to make sure our new tag codes make sense
states_complete %>% 
  subset(., eventual_repeat == "eventual_repeat") -> eventual_repeat_df_2


# FILTER FINAL FILE FOR EXPORT
states_complete %>% 
  subset(., !(life_stage %in% c("kelt", "repeat_kelt", "Juvenile"))) -> adults_only_states

# Let's make sure we're starting upstream of BON
adults_only_states %>% 
  group_by(tag_code_2) %>% 
  filter(row_number() == 1) -> adults_first_states

length(unique(adults_first_states$tag_code_2))
# 58,233 migrations

##### Remove ghost tags - use criteria from DART (provided by Susannah Iltis) #####

# DART applies these criteria for identifying ghost/shed tags for our tributary queries:
# (1) detected 3 or more years after last detection for all species tagged as adults
# (2) detected 5 years after last detection for all stages and species
# (3) steelhead tagged as juvenile with first detection 4 or more years after tagging
# (4) salmon tagged as juvenile with first detection 2 or more years after tagging
# (5) all tagged as adults with first detection 2 or more years after tagging.

# for our script, we've already removed the entire juvenile history at this point. So we can apply criteria (1)
# and remove everything that's occurred >= 3 years after the last detection

adults_only_states %>% 
  # subset(tag_code %in% bad_tags$tag_code) %>%  # for testing
  mutate(date_time = ymd_hms(date_time)) %>% 
  mutate(time_between_detections = ifelse(tag_code == lag(tag_code), time_after_release - lag(time_after_release), NA)) %>%
  # mutate(time_between_detections = ifelse(tag_code == lag(tag_code), interval(date_time, lag(date_time)), NA)) %>% 
  mutate(ghost_tag = ifelse(time_between_detections >= 365*3, "ghost", 
                            ifelse(is.na(time_between_detections) | time_between_detections < 365*3, "not_ghost", NA))) -> adults_only_states

# now note that everything that occurs after when we first note it's a ghost is also ghost activity
eventual_ghost_tags <- unique(subset(adults_only_states, ghost_tag == "ghost")$tag_code)
adults_only_states %>% 
  mutate(ghost_start = ifelse(ghost_tag == "ghost", order, NA)) %>% 
  subset(ghost_tag == "ghost") %>% 
  dplyr::select(tag_code, ghost_start) -> ghost_tag_start

adults_only_states %>% 
  ungroup() %>% 
  left_join(., ghost_tag_start, by = "tag_code") %>% 
  # group_by(tag_code) %>% 
  mutate(ghost_tag = ifelse(tag_code %in% eventual_ghost_tags & order >= ghost_start, "ghost", "not_ghost")) -> adults_only_states

# check on our ghost tags
adults_only_states %>% 
  group_by(tag_code) %>% 
  filter(any(ghost_tag == "ghost")) -> ghost_tags

# remove all of the ghost tag activity
adults_only_states %>% 
  subset(ghost_tag != "ghost") -> adults_only_states


##### Now, export this file, removing irrelevant columns #####
adults_only_states %>% 
  dplyr::select(tag_code, state, date_time, pathway, life_stage, tag_code_2) -> adults_only_states_for_export
write.csv(adults_only_states_for_export, here::here("intermediate_outputs", "adults_states_complete.csv"))


# Export three more files, into the appropriate modeling folder
# Snake: everything above ICH
# Upper Columbia: Everything above PRA
# Middle Columbia: Everything above BON, except Hood River

# load natal origins
natal_origin_table <- read.csv(here::here("Data", "covariate_data", "natal_origin_table.csv"))

# match natal origin to tag codes
tag_code_metadata %>% 
  left_join(natal_origin_table, by = "release_site_name") -> tag_code_metadata

tag_code_metadata %>% 
  mutate(ESU = ifelse(natal_origin %in% c("Tucannon_River", "Asotin_Creek", "Clearwater_River", "Salmon_River", "Grande_Ronde_River", "Imnaha_River"), "snake",
                      ifelse(natal_origin %in% c("Wenatchee_River", "Entiat_River", "Okanogan_River","Methow_River"), "upper_columbia",
                             ifelse(natal_origin %in% c("Deschutes_River", "Fifteenmile_Creek", "John_Day_River", "Umatilla_River", "Yakima_River", "Walla_Walla_River"), "middle_columbia",
                                    ifelse(natal_origin %in% c("Hood_River"), "lower_columbia", "ERROR"))))) -> tag_code_metadata

### Export natal origin information
# keep only the fish that are in the dataset
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

tag_code_metadata %>% 
  left_join(., origin_numeric, by = "natal_origin") %>% 
  mutate(rear_type_numeric = ifelse(rear_type_code %in% c("H", "U"), 2, 1))-> tag_code_metadata

# At this point, we need to recreate tag_code_metadata but with the tag_code_2 field
states_complete %>% 
  distinct(tag_code_2, .keep_all = TRUE) %>% 
  dplyr::select(tag_code, tag_code_2) -> tag_codes_2

# reformat this into origin_rear info
tag_codes_2 %>%
  left_join(dplyr::select(tag_code_metadata, tag_code, natal_origin_numeric, rear_type_numeric), by = "tag_code") %>% 
  dplyr::rename(natal_origin = natal_origin_numeric, rear_type = rear_type_numeric) %>% 
  dplyr::select(-tag_code) -> origin_rear_actual

write.csv(origin_rear_actual,here::here("Data", "covariate_data", "origin_rear_actual.csv"))

### Export movement files

snake_tags <- subset(tag_code_metadata, ESU == "snake")$tag_code
upper_columbia_tags <- subset(tag_code_metadata, ESU == "upper_columbia")$tag_code
middle_columbia_tags <- subset(tag_code_metadata, ESU == "middle_columbia")$tag_code

subset(adults_only_states_for_export, tag_code %in% snake_tags) -> snake_adults_only_states
subset(adults_only_states_for_export, tag_code %in% upper_columbia_tags) -> upper_columbia_adults_only_states
subset(adults_only_states_for_export, tag_code %in% middle_columbia_tags) -> middle_columbia_adults_only_states

# write them to the the intermediate outputs folders
write.csv(snake_adults_only_states, here::here("intermediate_outputs", "adults_states_complete", "snake_adults_states_complete.csv"))
write.csv(upper_columbia_adults_only_states, here::here("intermediate_outputs", "adults_states_complete", "upper_columbia_adults_states_complete.csv"))
write.csv(middle_columbia_adults_only_states, here::here("intermediate_outputs", "adults_states_complete", "middle_columbia_adults_states_complete.csv"))

# write them to the modeling folders
# write.csv(snake_adults_only_states, here::here("Stan", "snake_river_wild", "snake_adults_states_complete.csv"))
# write.csv(snake_adults_only_states, here::here("Stan", "snake_river_hatchery", "snake_adults_states_complete.csv"))
# write.csv(upper_columbia_adults_only_states, here::here("Stan", "upper_columbia_wild", "upper_columbia_adults_states_complete.csv"))
# write.csv(upper_columbia_adults_only_states, here::here("Stan", "upper_columbia_hatchery", "upper_columbia_adults_states_complete.csv"))
# write.csv(middle_columbia_adults_only_states, here::here("Stan", "middle_columbia_wild", "middle_columbia_adults_states_complete.csv"))
# write.csv(middle_columbia_adults_only_states, here::here("Stan", "middle_columbia_hatchery", "middle_columbia_adults_states_complete.csv"))