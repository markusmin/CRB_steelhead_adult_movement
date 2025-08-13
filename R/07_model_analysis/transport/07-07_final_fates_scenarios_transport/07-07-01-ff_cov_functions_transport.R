# 07-07-01-ff_cov_functions.R

# This script contains all of the functions for running the simulations for final fates and covariates.
# This script will be sourced by the other scripts that run individual populations through these simulations.

##### Load the model runs and packages #####


# If we are running this on hyak, we need to change the source path so that the remaining file paths will work
setwd("/gscratch/scrubbed/mmin/")

# load the model runs and packages
source("R/07_model_analysis/transport/07-01_load_stan_models_transport.R")

#### Select representative warm vs. cold years ####

run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14",
              "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23", "23/24", "24/25")
run_year_start <- seq(ymd_hms("2004-06-01 00:00:00"), ymd_hms("2024-06-01 00:00:00"), by = "years")
run_year_end <- seq(ymd_hms("2005-05-31 23:59:59"), ymd_hms("2025-05-31 23:59:59"), by = "years")
run_year_numeric = seq(4, 24, 1)

run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)

# subset states_dates for origin and population to pull out just month day, and then
# use those months and days to index temps from that year
dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2024-12-31"), by = "days"),
                           date = 1:length(seq(ymd("2005-06-01"), ymd("2024-12-31"), by = "days")))
dates_actual %>% 
  mutate(month_day = format(as.Date(date_actual), "%m-%d")) -> dates_actual

temp_ts <- as.data.frame(SRW_T_envir$data$temperature_data)
dates <- seq(ymd("2005-06-01"), ymd("2024-12-31"), by = "days")
years <- year(dates)
temp_ts$date <- dates
temp_ts$year <- years

# Get the median for each run year for both winter/spring (temp0) and summer/fall (temp1)

temp_ts %>% 
  rowwise() %>%
  dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date & run_year_end >= date)$run_year) %>% 
  dplyr::mutate(season = ifelse(yday(date)<152, 0, 1))  %>% 
  # drop 24/25 because it's not a complete year (and it isn't modeled)
  filter(!(run_year == "24/25")) %>% 
  subset(season == 0) %>% 
  group_by(run_year) %>% 
  summarise_all(median) %>% 
  arrange(BON) -> season0_annual_temp_medians

temp_ts %>% 
  rowwise() %>%
  dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date & run_year_end >= date)$run_year) %>% 
  dplyr::mutate(season = ifelse(yday(date)<152, 0, 1))  %>% 
  # drop 24/25 because it's not a complete year (and it isn't modeled)
  filter(!(run_year == "24/25")) %>% 
  subset(season == 1) %>% 
  group_by(run_year) %>% 
  summarise_all(median) %>% 
  arrange(BON) -> season1_annual_temp_medians

spill_ts <- as.data.frame(SRW_T_envir$data$spill_window_data)
dates <- seq(ymd("2005-06-01"), ymd("2024-12-31"), by = "days")
years <- year(dates)
spill_ts$date <- dates
spill_ts$year <- years

spill_ts %>% 
  rowwise() %>%
  dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date & run_year_end >= date)$run_year) %>% 
  dplyr::mutate(season = ifelse(yday(date)<152, 0, 1))  %>% 
  # drop 24/25 because it's not a complete year (and it isn't modeled)
  filter(!(run_year == "24/25")) %>% 
  subset(season == 0) %>% 
  group_by(run_year) %>% 
  summarise_all(median) %>% 
  arrange(BON) -> season0_annual_spill_medians

spill_ts %>% 
  rowwise() %>%
  dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date & run_year_end >= date)$run_year) %>% 
  dplyr::mutate(season = ifelse(yday(date)<152, 0, 1))  %>% 
  # drop 24/25 because it's not a complete year (and it isn't modeled)
  filter(!(run_year == "24/25")) %>% 
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
midpoint_year <- round(median(1:nrow(season1_annual_temp_medians)))
min_year <- 1
max_year <- nrow(season1_annual_temp_medians)
quantile_25 <- round(quantile(1:nrow(season1_annual_temp_medians), 0.25))
quantile_75 <- round(quantile(1:nrow(season1_annual_temp_medians), 0.75))

rep_years <- data.frame(fix_run_year = season1_annual_temp_medians[c(min_year, quantile_25, midpoint_year, quantile_75, max_year),]$run_year,
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


final_fates_simulation_fixcov_SRW_T <- function(nsim,
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
      sample_date[i] <- sample(subset(SRW_T_states_dates, state == i)$date, 1)
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
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRW_T_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRW_T_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRW_T_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRW_T_states_dates, state == i  & direction == "downstream")$season)
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
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
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
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    
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
    
    
    
    
    
    SRW_T_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> SRW_T_states_month_day
    
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
        if(nrow(subset(SRW_T_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRW_T_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_T_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRW_T_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_T_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRW_T_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_T_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRW_T_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
        if(nrow(subset(SRW_T_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRW_T_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_T_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRW_T_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_T_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRW_T_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_T_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRW_T_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
    possible_movements <- SRW_T_envir$data$movements[, "col"][SRW_T_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[j], season_upstream[j]))
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_T[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRW_T[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRW_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    btemp1xorigin6_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin6 +
                                                                    borigin1_array_SRW_T[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW_T[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW_T[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW_T[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW_T[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW_T[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_T[i,possible_movements,iter] +
                    btemp1_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRW_T[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRW_T[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin6 +
                    borigin1_array_SRW_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_T[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_T[i,possible_movements,iter]*origin6))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_T[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRW_T[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRW_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    btemp0xorigin6_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin6 +
                                                                    borigin1_array_SRW_T[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW_T[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW_T[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW_T[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW_T[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW_T[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_T[i,possible_movements,iter] +
                    btemp0_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRW_T[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRW_T[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin6 +
                    borigin1_array_SRW_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_T[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_T[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRW_T_envir$data$movements[, "col"][SRW_T_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[j], season_downstream[j]))
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_T[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRW_T[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRW_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      btemp1xorigin6_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin6 +
                                                                      borigin1_array_SRW_T[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW_T[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW_T[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW_T[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW_T[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW_T[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_T[i,possible_movements,iter] +
                    btemp1_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRW_T[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRW_T[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin6 +
                    borigin1_array_SRW_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_T[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_T[i,possible_movements,iter]*origin6))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_T[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRW_T[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRW_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      btemp0xorigin6_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin6 +
                                                                      borigin1_array_SRW_T[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW_T[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW_T[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW_T[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW_T[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW_T[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_T[i,possible_movements,iter] +
                    btemp0_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRW_T[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRW_T[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin6 +
                    borigin1_array_SRW_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_T[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_T[i,possible_movements,iter]*origin6))
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

final_fates_simulation_fixcov_SRW_NT <- function(nsim,
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
      sample_date[i] <- sample(subset(SRW_NT_states_dates, state == i)$date, 1)
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
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRW_NT_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRW_NT_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRW_NT_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRW_NT_states_dates, state == i  & direction == "downstream")$season)
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
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
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
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    
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
    
    
    
    
    
    SRW_NT_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> SRW_NT_states_month_day
    
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
        if(nrow(subset(SRW_NT_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRW_NT_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_NT_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRW_NT_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_NT_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRW_NT_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_NT_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRW_NT_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
        if(nrow(subset(SRW_NT_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRW_NT_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_NT_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRW_NT_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_NT_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRW_NT_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_NT_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRW_NT_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
    possible_movements <- SRW_NT_envir$data$movements[, "col"][SRW_NT_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[j], season_upstream[j]))
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_NT[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRW_NT[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRW_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    btemp1xorigin6_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin6 +
                                                                    borigin1_array_SRW_NT[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW_NT[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW_NT[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW_NT[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW_NT[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW_NT[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_NT[i,possible_movements,iter] +
                    btemp1_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRW_NT[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRW_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin6 +
                    borigin1_array_SRW_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_NT[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_NT[i,possible_movements,iter]*origin6))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_NT[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRW_NT[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRW_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    btemp0xorigin6_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin6 +
                                                                    borigin1_array_SRW_NT[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW_NT[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW_NT[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW_NT[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW_NT[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW_NT[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_NT[i,possible_movements,iter] +
                    btemp0_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRW_NT[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRW_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin6 +
                    borigin1_array_SRW_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_NT[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_NT[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRW_NT_envir$data$movements[, "col"][SRW_NT_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[j], season_downstream[j]))
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_NT[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRW_NT[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRW_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      btemp1xorigin6_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin6 +
                                                                      borigin1_array_SRW_NT[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW_NT[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW_NT[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW_NT[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW_NT[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW_NT[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_NT[i,possible_movements,iter] +
                    btemp1_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRW_NT[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRW_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin6 +
                    borigin1_array_SRW_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_NT[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_NT[i,possible_movements,iter]*origin6))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_NT[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRW_NT[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRW_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      btemp0xorigin6_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin6 +
                                                                      borigin1_array_SRW_NT[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW_NT[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW_NT[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW_NT[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW_NT[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW_NT[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_NT[i,possible_movements,iter] +
                    btemp0_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRW_NT[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRW_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin6 +
                    borigin1_array_SRW_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_NT[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_NT[i,possible_movements,iter]*origin6))
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

final_fates_simulation_fixcov_SRH_T <- function(nsim,
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
      sample_date[i] <- sample(subset(SRH_T_states_dates, state == i)$date, 1)
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
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRH_T_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRH_T_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRH_T_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRH_T_states_dates, state == i  & direction == "downstream")$season)
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
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
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
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    
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
    
    
    
    
    
    SRH_T_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> SRH_T_states_month_day
    
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
        if(nrow(subset(SRH_T_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRH_T_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_T_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRH_T_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_T_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRH_T_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_T_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRH_T_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
        if(nrow(subset(SRH_T_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRH_T_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_T_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRH_T_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_T_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRH_T_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_T_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRH_T_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
    possible_movements <- SRH_T_envir$data$movements[, "col"][SRH_T_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[j], season_upstream[j]))
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_T[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRH_T[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRH_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    borigin1_array_SRH_T[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH_T[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH_T[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH_T[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH_T[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_T[i,possible_movements,iter] +
                    btemp1_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRH_T[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRH_T[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    borigin1_array_SRH_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_T[i,possible_movements,iter]*origin5))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_T[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRH_T[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRH_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    borigin1_array_SRH_T[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH_T[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH_T[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH_T[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH_T[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_T[i,possible_movements,iter] +
                    btemp0_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRH_T[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRH_T[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    borigin1_array_SRH_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_T[i,possible_movements,iter]*origin5))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRH_T_envir$data$movements[, "col"][SRH_T_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[j], season_downstream[j]))
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_T[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRH_T[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRH_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      borigin1_array_SRH_T[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH_T[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH_T[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH_T[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH_T[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_T[i,possible_movements,iter] +
                    btemp1_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRH_T[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRH_T[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    borigin1_array_SRH_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_T[i,possible_movements,iter]*origin5))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_T[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRH_T[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRH_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      borigin1_array_SRH_T[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH_T[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH_T[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH_T[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH_T[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_T[i,possible_movements,iter] +
                    btemp0_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRH_T[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRH_T[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    borigin1_array_SRH_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_T[i,possible_movements,iter]*origin5))
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

final_fates_simulation_fixcov_SRH_NT <- function(nsim,
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
      sample_date[i] <- sample(subset(SRH_NT_states_dates, state == i)$date, 1)
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
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRH_NT_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRH_NT_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRH_NT_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRH_NT_states_dates, state == i  & direction == "downstream")$season)
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
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
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
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    
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
    
    
    
    
    
    SRH_NT_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> SRH_NT_states_month_day
    
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
        if(nrow(subset(SRH_NT_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRH_NT_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_NT_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRH_NT_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_NT_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRH_NT_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_NT_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRH_NT_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
        if(nrow(subset(SRH_NT_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRH_NT_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_NT_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRH_NT_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_NT_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRH_NT_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_NT_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRH_NT_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
    possible_movements <- SRH_NT_envir$data$movements[, "col"][SRH_NT_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[j], season_upstream[j]))
      
      if (upstream_season_sim == 1){
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_NT[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRH_NT[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRH_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    borigin1_array_SRH_NT[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH_NT[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH_NT[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH_NT[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH_NT[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_NT[i,possible_movements,iter] +
                    btemp1_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRH_NT[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRH_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    borigin1_array_SRH_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_NT[i,possible_movements,iter]*origin5))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_NT[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRH_NT[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRH_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    borigin1_array_SRH_NT[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH_NT[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH_NT[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH_NT[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH_NT[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_NT[i,possible_movements,iter] +
                    btemp0_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRH_NT[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRH_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    borigin1_array_SRH_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_NT[i,possible_movements,iter]*origin5))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRH_NT_envir$data$movements[, "col"][SRH_NT_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    for (j in 1:length(possible_movements)){
      # toggle temp0 v temp1 based on resampling season
      downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[j], season_downstream[j]))
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_NT[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRH_NT[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRH_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      borigin1_array_SRH_NT[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH_NT[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH_NT[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH_NT[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH_NT[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_NT[i,possible_movements,iter] +
                    btemp1_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRH_NT[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRH_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    borigin1_array_SRH_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_NT[i,possible_movements,iter]*origin5))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_NT[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRH_NT[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRH_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      borigin1_array_SRH_NT[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH_NT[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH_NT[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH_NT[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH_NT[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_NT[i,possible_movements,iter] +
                    btemp0_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRH_NT[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRH_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    borigin1_array_SRH_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_NT[i,possible_movements,iter]*origin5))
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
    ntransitions <- length(transitions[i,][which(!(transitions[i,] %in% c(0,41)))])
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
      transitions_df$previous_state[i] <- ifelse(state_history[i-1]==41, 0, state_history[i-1]) # all first detections come from state 0 (fish purgatory)
      transitions_df$date[i] <- transition_date_history[i]
      transitions_df$season[i] <- transition_seasons_history[i]
    }
  }
  
  # drop the non-data
  transitions_df %>% 
    filter(!(state %in% c(0,41))) -> transitions_df
  
  # determine if the fish is going upstream or downstream
  transitions_df %>% 
    mutate(direction = ifelse(state > previous_state, "upstream", "downstream")) -> transitions_df
  
  return(transitions_df)
}



SRW_T_states_dates <- get_states_dates_direction(envir = SRW_T_envir)
SRW_NT_states_dates <- get_states_dates_direction(envir = SRW_NT_envir)
SRH_T_states_dates <- get_states_dates_direction(envir = SRH_T_envir)
SRH_NT_states_dates <- get_states_dates_direction(envir = SRH_NT_envir)

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
    ntransitions <- length(transitions[i,][which(!(transitions[i,] %in% c(0,41)))])
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
      transitions_df$previous_state[i] <- ifelse(state_history[i-1]==41, 0, state_history[i-1]) # all first detections come from state 0 (fish purgatory)
      transitions_df$date[i] <- transition_date_history[i]
      transitions_df$season[i] <- transition_seasons_history[i]
    }
  }
  
  # drop the non-data
  transitions_df %>% 
    filter(!(state %in% c(0,41))) -> transitions_df
  
  # determine if the fish is going upstream or downstream
  transitions_df %>% 
    mutate(direction = ifelse(state > previous_state, "upstream", "downstream")) -> transitions_df
  
  return(transitions_df)
  
}


compare_final_fate_fixcov_rear_type_SR_T <- function(niter, nsim,
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
      sim_wild <- final_fates_simulation_fixcov_SRW_T(nsim = nsim,
                                                      start_state = 2, states_dates = get_origin_states_dates(envir = SRW_T_envir, origin_select = origin_select, rear = "wild"),
                                                      origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                      origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                      temp_data = SRW_T_envir$data$temperature_data, spillwindow_data = SRW_T_envir$data$spill_window_data, 
                                                      winterspill_data = SRW_T_envir$data$winter_spill_days_data,
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
      sim_hatchery <- final_fates_simulation_fixcov_SRH_T(nsim = nsim,
                                                          start_state = 2, states_dates = get_origin_states_dates(envir = SRH_T_envir, origin_select = origin_select, rear = "hatchery"),
                                                          origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], origin3 = hatchery_origin_params[3], 
                                                          origin4 = hatchery_origin_params[4], origin5 = hatchery_origin_params[5], 
                                                          temp_data = SRH_T_envir$data$temperature_data, spillwindow_data = SRH_T_envir$data$spill_window_data, 
                                                          winterspill_data = SRH_T_envir$data$winter_spill_days_data,
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
      sim_wild <- final_fates_simulation_fixcov_SRW_T(nsim = nsim,
                                                      start_state = 2, states_dates = get_origin_states_dates(envir = SRW_T_envir, origin_select = origin_select, rear = "wild"),
                                                      origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                      origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                      temp_data = SRW_T_envir$data$temperature_data, spillwindow_data = SRW_T_envir$data$spill_window_data, 
                                                      winterspill_data = SRW_T_envir$data$winter_spill_days_data,
                                                      fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                      fix_temp_season1_value = fix_temp_season1_value, 
                                                      fix_spillwindow_season0_value = fix_spillwindow_season0_value, fix_spillwindow_season1_value = fix_spillwindow_season1_value,
                                                      fix_winterspill_value = fix_winterspill_value,
                                                      fix_run_year = fix_run_year)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_fixcov_SRH_T(nsim = nsim,
                                                          start_state = 2, states_dates = get_origin_states_dates(envir = SRH_T_envir, origin_select = origin_select, rear = "hatchery"),
                                                          origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], origin3 = hatchery_origin_params[3], 
                                                          origin4 = hatchery_origin_params[4], origin5 = hatchery_origin_params[5], 
                                                          temp_data = SRH_T_envir$data$temperature_data, spillwindow_data = SRH_T_envir$data$spill_window_data, 
                                                          winterspill_data = SRH_T_envir$data$winter_spill_days_data,
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

compare_final_fate_fixcov_rear_type_SR_NT <- function(niter, nsim,
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
      sim_wild <- final_fates_simulation_fixcov_SRW_NT(nsim = nsim,
                                                       start_state = 2, states_dates = get_origin_states_dates(envir = SRW_NT_envir, origin_select = origin_select, rear = "wild"),
                                                       origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                       origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                       temp_data = SRW_NT_envir$data$temperature_data, spillwindow_data = SRW_NT_envir$data$spill_window_data, 
                                                       winterspill_data = SRW_NT_envir$data$winter_spill_days_data,
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
      sim_hatchery <- final_fates_simulation_fixcov_SRH_NT(nsim = nsim,
                                                           start_state = 2, states_dates = get_origin_states_dates(envir = SRH_NT_envir, origin_select = origin_select, rear = "hatchery"),
                                                           origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], origin3 = hatchery_origin_params[3], 
                                                           origin4 = hatchery_origin_params[4], origin5 = hatchery_origin_params[5], 
                                                           temp_data = SRH_NT_envir$data$temperature_data, spillwindow_data = SRH_NT_envir$data$spill_window_data, 
                                                           winterspill_data = SRH_NT_envir$data$winter_spill_days_data,
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
      sim_wild <- final_fates_simulation_fixcov_SRW_NT(nsim = nsim,
                                                       start_state = 2, states_dates = get_origin_states_dates(envir = SRW_NT_envir, origin_select = origin_select, rear = "wild"),
                                                       origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                       origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                       temp_data = SRW_NT_envir$data$temperature_data, spillwindow_data = SRW_NT_envir$data$spill_window_data, 
                                                       winterspill_data = SRW_NT_envir$data$winter_spill_days_data,
                                                       fix_state = fix_state, fix_temp_season0_value = fix_temp_season0_value, 
                                                       fix_temp_season1_value = fix_temp_season1_value, 
                                                       fix_spillwindow_season0_value = fix_spillwindow_season0_value, fix_spillwindow_season1_value = fix_spillwindow_season1_value,
                                                       fix_winterspill_value = fix_winterspill_value,
                                                       fix_run_year = fix_run_year)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_fixcov_SRH_NT(nsim = nsim,
                                                           start_state = 2, states_dates = get_origin_states_dates(envir = SRH_NT_envir, origin_select = origin_select, rear = "hatchery"),
                                                           origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], origin3 = hatchery_origin_params[3], 
                                                           origin4 = hatchery_origin_params[4], origin5 = hatchery_origin_params[5], 
                                                           temp_data = SRH_NT_envir$data$temperature_data, spillwindow_data = SRH_NT_envir$data$spill_window_data, 
                                                           winterspill_data = SRH_NT_envir$data$winter_spill_days_data,
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

final_fates_simulation_fixcov_times_SRW_T <- function(nsim,
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
      sample_date[i] <- sample(subset(SRW_T_states_dates, state == i)$date, 1)
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
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRW_T_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRW_T_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRW_T_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRW_T_states_dates, state == i  & direction == "downstream")$season)
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
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
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
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRW_T_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRW_T_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    
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
    
    
    
    
    
    SRW_T_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> SRW_T_states_month_day
    
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
        if(nrow(subset(SRW_T_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRW_T_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_T_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRW_T_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_T_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRW_T_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_T_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRW_T_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
        if(nrow(subset(SRW_T_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRW_T_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_T_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRW_T_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_T_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRW_T_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_T_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRW_T_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
    
    SRW_T_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(fix = ifelse(yday(date_actual) > fix_start_date &
                            yday(date_actual) < fix_end_date, 1, 0)) %>% 
      mutate(fix = ifelse(is.na(fix), 0, fix)) -> SRW_T_states_dates_fix
    
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
        if(nrow(subset(SRW_T_states_dates_fix, state == i & direction == "upstream" & season == 0)) == 0){
          season0_upstream_fixprop[i] <- 0
        } else {
          season0_upstream_fixprop[i] <- sum(subset(SRW_T_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
            length(subset(SRW_T_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
        }
      } else {
        season0_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
      }
      # proportion of season1, upstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_T_states_dates_fix, state == i & direction == "upstream" & season == 1)) == 0){
          season1_upstream_fixprop[i] <- 0
        } else {
          season1_upstream_fixprop[i] <- sum(subset(SRW_T_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
            length(subset(SRW_T_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)
        }
      } else {
        season1_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1 )$fix)
      }
      # proportion of season0, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_T_states_dates_fix, state == i & direction == "downstream" & season == 0)) == 0){
          season0_downstream_fixprop[i] <- 0
        } else {
          season0_downstream_fixprop[i] <- sum(subset(SRW_T_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
            length(subset(SRW_T_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
        }
      } else {
        season0_downstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
      }
      # proportion of season1, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_T_states_dates_fix, state == i & direction == "downstream" & season == 1)) == 0){
          season1_downstream_fixprop[i] <- 0
        } else {
          season1_downstream_fixprop[i] <- sum(subset(SRW_T_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)/
            length(subset(SRW_T_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)
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
    
    SRW_T_states_month_day %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(drop_date = ifelse(yday(date_actual.x) > fix_start_date &
                                  yday(date_actual.x) < fix_end_date, 1, 0)) %>% 
      mutate(drop_date = ifelse(is.na(drop_date), 0, drop_date)) %>% 
      filter(drop_date == 0) -> SRW_T_states_month_day_trimmed
    
    
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRW_T_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(SRW_T_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_T_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(SRW_T_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_T_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(SRW_T_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_T_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season1_notfixed[i] <- median(spillwindow_data[subset(SRW_T_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
    possible_movements <- SRW_T_envir$data$movements[, "col"][SRW_T_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      
      if (upstream_season_sim == 1){
        
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_T[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRW_T[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRW_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    btemp1xorigin6_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin6 +
                                                                    borigin1_array_SRW_T[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW_T[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW_T[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW_T[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW_T[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW_T[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_T[i,possible_movements,iter] +
                    btemp1_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRW_T[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRW_T[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin6 +
                    borigin1_array_SRW_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_T[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_T[i,possible_movements,iter]*origin6))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_T[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRW_T[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRW_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    btemp0xorigin6_array_SRW_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin6 +
                                                                    borigin1_array_SRW_T[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW_T[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW_T[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW_T[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW_T[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW_T[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_T[i,possible_movements,iter] +
                    btemp0_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRW_T[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRW_T[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin6 +
                    borigin1_array_SRW_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_T[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_T[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRW_T_envir$data$movements[, "col"][SRW_T_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_T[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRW_T[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRW_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      btemp1xorigin6_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin6 +
                                                                      borigin1_array_SRW_T[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW_T[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW_T[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW_T[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW_T[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW_T[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_T[i,possible_movements,iter] +
                    btemp1_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRW_T[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRW_T[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin6 +
                    borigin1_array_SRW_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_T[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_T[i,possible_movements,iter]*origin6))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_T[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRW_T[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRW_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      btemp0xorigin6_array_SRW_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin6 +
                                                                      borigin1_array_SRW_T[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW_T[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW_T[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW_T[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW_T[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW_T[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_T[i,possible_movements,iter] +
                    btemp0_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRW_T[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRW_T[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin6 +
                    borigin1_array_SRW_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_T[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_T[i,possible_movements,iter]*origin6))
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

final_fates_simulation_fixcov_times_SRW_NT <- function(nsim,
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
      sample_date[i] <- sample(subset(SRW_NT_states_dates, state == i)$date, 1)
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
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRW_NT_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRW_NT_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRW_NT_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRW_NT_states_dates, state == i  & direction == "downstream")$season)
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
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
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
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRW_NT_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRW_NT_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    
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
    
    
    
    
    
    SRW_NT_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> SRW_NT_states_month_day
    
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
        if(nrow(subset(SRW_NT_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRW_NT_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_NT_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRW_NT_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_NT_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRW_NT_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_NT_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRW_NT_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
        if(nrow(subset(SRW_NT_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRW_NT_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_NT_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRW_NT_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_NT_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRW_NT_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_NT_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRW_NT_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
    
    SRW_NT_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(fix = ifelse(yday(date_actual) > fix_start_date &
                            yday(date_actual) < fix_end_date, 1, 0)) %>% 
      mutate(fix = ifelse(is.na(fix), 0, fix)) -> SRW_NT_states_dates_fix
    
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
        if(nrow(subset(SRW_NT_states_dates_fix, state == i & direction == "upstream" & season == 0)) == 0){
          season0_upstream_fixprop[i] <- 0
        } else {
          season0_upstream_fixprop[i] <- sum(subset(SRW_NT_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
            length(subset(SRW_NT_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
        }
      } else {
        season0_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
      }
      # proportion of season1, upstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_NT_states_dates_fix, state == i & direction == "upstream" & season == 1)) == 0){
          season1_upstream_fixprop[i] <- 0
        } else {
          season1_upstream_fixprop[i] <- sum(subset(SRW_NT_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
            length(subset(SRW_NT_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)
        }
      } else {
        season1_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1 )$fix)
      }
      # proportion of season0, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_NT_states_dates_fix, state == i & direction == "downstream" & season == 0)) == 0){
          season0_downstream_fixprop[i] <- 0
        } else {
          season0_downstream_fixprop[i] <- sum(subset(SRW_NT_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
            length(subset(SRW_NT_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
        }
      } else {
        season0_downstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
      }
      # proportion of season1, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_NT_states_dates_fix, state == i & direction == "downstream" & season == 1)) == 0){
          season1_downstream_fixprop[i] <- 0
        } else {
          season1_downstream_fixprop[i] <- sum(subset(SRW_NT_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)/
            length(subset(SRW_NT_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)
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
    
    SRW_NT_states_month_day %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(drop_date = ifelse(yday(date_actual.x) > fix_start_date &
                                  yday(date_actual.x) < fix_end_date, 1, 0)) %>% 
      mutate(drop_date = ifelse(is.na(drop_date), 0, drop_date)) %>% 
      filter(drop_date == 0) -> SRW_NT_states_month_day_trimmed
    
    
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRW_NT_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(SRW_NT_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRW_NT_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(SRW_NT_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRW_NT_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(SRW_NT_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRW_NT_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season1_notfixed[i] <- median(spillwindow_data[subset(SRW_NT_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
    possible_movements <- SRW_NT_envir$data$movements[, "col"][SRW_NT_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      
      if (upstream_season_sim == 1){
        
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_NT[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRW_NT[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRW_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    btemp1xorigin6_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin6 +
                                                                    borigin1_array_SRW_NT[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW_NT[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW_NT[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW_NT[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW_NT[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW_NT[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_NT[i,possible_movements,iter] +
                    btemp1_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRW_NT[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRW_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin6 +
                    borigin1_array_SRW_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_NT[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_NT[i,possible_movements,iter]*origin6))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_NT[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRW_NT[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRW_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    btemp0xorigin6_array_SRW_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin6 +
                                                                    borigin1_array_SRW_NT[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW_NT[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW_NT[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW_NT[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW_NT[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW_NT[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_NT[i,possible_movements,iter] +
                    btemp0_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRW_NT[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRW_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin6 +
                    borigin1_array_SRW_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_NT[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_NT[i,possible_movements,iter]*origin6))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRW_NT_envir$data$movements[, "col"][SRW_NT_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_NT[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRW_NT[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRW_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      btemp1xorigin6_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin6 +
                                                                      borigin1_array_SRW_NT[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW_NT[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW_NT[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW_NT[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW_NT[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW_NT[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_NT[i,possible_movements,iter] +
                    btemp1_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRW_NT[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRW_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    btemp1xorigin6_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin6 +
                    borigin1_array_SRW_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_NT[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_NT[i,possible_movements,iter]*origin6))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW_NT[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRW_NT[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRW_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      btemp0xorigin6_array_SRW_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin6 +
                                                                      borigin1_array_SRW_NT[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRW_NT[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRW_NT[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRW_NT[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRW_NT[i,possible_movements[j],iter]*origin5 +
                                                                      borigin6_array_SRW_NT[i,possible_movements[j],iter]*origin6)/
          sum(exp(b0_array_SRW_NT[i,possible_movements,iter] +
                    btemp0_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRW_NT[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRW_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    btemp0xorigin6_array_SRW_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin6 +
                    borigin1_array_SRW_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRW_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRW_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRW_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRW_NT[i,possible_movements,iter]*origin5 +
                    borigin6_array_SRW_NT[i,possible_movements,iter]*origin6))
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

final_fates_simulation_fixcov_times_SRH_T <- function(nsim,
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
      sample_date[i] <- sample(subset(SRH_T_states_dates, state == i)$date, 1)
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
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRH_T_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRH_T_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRH_T_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRH_T_states_dates, state == i  & direction == "downstream")$season)
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
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
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
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRH_T_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRH_T_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    
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
    
    
    
    
    
    SRH_T_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> SRH_T_states_month_day
    
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
        if(nrow(subset(SRH_T_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRH_T_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_T_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRH_T_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_T_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRH_T_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_T_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRH_T_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
        if(nrow(subset(SRH_T_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRH_T_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_T_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRH_T_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_T_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRH_T_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_T_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRH_T_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
    
    SRH_T_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(fix = ifelse(yday(date_actual) > fix_start_date &
                            yday(date_actual) < fix_end_date, 1, 0)) %>% 
      mutate(fix = ifelse(is.na(fix), 0, fix)) -> SRH_T_states_dates_fix
    
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
        if(nrow(subset(SRH_T_states_dates_fix, state == i & direction == "upstream" & season == 0)) == 0){
          season0_upstream_fixprop[i] <- 0
        } else {
          season0_upstream_fixprop[i] <- sum(subset(SRH_T_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
            length(subset(SRH_T_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
        }
      } else {
        season0_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
      }
      # proportion of season1, upstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_T_states_dates_fix, state == i & direction == "upstream" & season == 1)) == 0){
          season1_upstream_fixprop[i] <- 0
        } else {
          season1_upstream_fixprop[i] <- sum(subset(SRH_T_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
            length(subset(SRH_T_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)
        }
      } else {
        season1_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1 )$fix)
      }
      # proportion of season0, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_T_states_dates_fix, state == i & direction == "downstream" & season == 0)) == 0){
          season0_downstream_fixprop[i] <- 0
        } else {
          season0_downstream_fixprop[i] <- sum(subset(SRH_T_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
            length(subset(SRH_T_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
        }
      } else {
        season0_downstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
      }
      # proportion of season1, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_T_states_dates_fix, state == i & direction == "downstream" & season == 1)) == 0){
          season1_downstream_fixprop[i] <- 0
        } else {
          season1_downstream_fixprop[i] <- sum(subset(SRH_T_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)/
            length(subset(SRH_T_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)
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
    
    SRH_T_states_month_day %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(drop_date = ifelse(yday(date_actual.x) > fix_start_date &
                                  yday(date_actual.x) < fix_end_date, 1, 0)) %>% 
      mutate(drop_date = ifelse(is.na(drop_date), 0, drop_date)) %>% 
      filter(drop_date == 0) -> SRH_T_states_month_day_trimmed
    
    
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRH_T_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(SRH_T_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_T_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(SRH_T_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_T_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(SRH_T_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_T_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season1_notfixed[i] <- median(spillwindow_data[subset(SRH_T_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
    possible_movements <- SRH_T_envir$data$movements[, "col"][SRH_T_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      
      if (upstream_season_sim == 1){
        
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_T[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRH_T[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRH_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    borigin1_array_SRH_T[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH_T[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH_T[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH_T[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH_T[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_T[i,possible_movements,iter] +
                    btemp1_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRH_T[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRH_T[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH_T[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    borigin1_array_SRH_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_T[i,possible_movements,iter]*origin5))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_T[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRH_T[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRH_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRH_T[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    borigin1_array_SRH_T[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH_T[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH_T[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH_T[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH_T[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_T[i,possible_movements,iter] +
                    btemp0_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRH_T[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRH_T[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH_T[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    borigin1_array_SRH_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_T[i,possible_movements,iter]*origin5))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRH_T_envir$data$movements[, "col"][SRH_T_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_T[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRH_T[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRH_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      borigin1_array_SRH_T[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH_T[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH_T[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH_T[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH_T[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_T[i,possible_movements,iter] +
                    btemp1_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRH_T[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRH_T[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH_T[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    borigin1_array_SRH_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_T[i,possible_movements,iter]*origin5))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_T[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRH_T[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRH_T[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRH_T[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      borigin1_array_SRH_T[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH_T[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH_T[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH_T[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH_T[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_T[i,possible_movements,iter] +
                    btemp0_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRH_T[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRH_T[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH_T[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    borigin1_array_SRH_T[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_T[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_T[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_T[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_T[i,possible_movements,iter]*origin5))
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

final_fates_simulation_fixcov_times_SRH_NT <- function(nsim,
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
      sample_date[i] <- sample(subset(SRH_NT_states_dates, state == i)$date, 1)
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
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "upstream")) == 0){
        season_upstream[i] <- 0
      } else {
        season_upstream[i] <- sum(subset(SRH_NT_states_dates, state == i  & direction == "upstream")$season)/length(subset(SRH_NT_states_dates, state == i  & direction == "upstream")$season)
      }
    } else {
      season_upstream[i] <- sum(subset(states_dates, state == i  & direction == "upstream")$season)/length(subset(states_dates, state == i  & direction == "upstream")$season)
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "downstream")) == 0){
        season_downstream[i] <- 0
      } else {
        season_downstream[i] <- sum(subset(SRH_NT_states_dates, state == i  & direction == "downstream")$season)/length(subset(SRH_NT_states_dates, state == i  & direction == "downstream")$season)
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
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        temp_upstream_season0[i] <- 0
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      temp_upstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        temp_downstream_season0[i] <- 0
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      temp_downstream_season0[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        temp_upstream_season1[i] <- 0
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      temp_upstream_season1[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        temp_downstream_season1[i] <- 0
      } else {
        temp_downstream_season1[i] <- median(temp_data[subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
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
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 0)) == 0){
        spillwindow_upstream_season0[i] <- 0
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 0)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 0)) == 0){
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 0)) == 0){
        spillwindow_downstream_season0[i] <- 0
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 0)$date,i])
      }
    } else {
      spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 0)$date,i])
    }
    # season 1
    if(nrow(subset(states_dates, state == i  & direction == "upstream" & season == 1)) == 0){
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 1)) == 0){
        spillwindow_upstream_season1[i] <- 0
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRH_NT_states_dates, state == i & direction == "upstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream" & season == 1)$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream" & season == 1)) == 0){
      if(nrow(subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 1)) == 0){
        spillwindow_downstream_season1[i] <- 0
      } else {
        spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRH_NT_states_dates, state == i & direction == "downstream" & season == 1)$date,i])
      }
    } else {
      spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream" & season == 1)$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    
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
    
    
    
    
    
    SRH_NT_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> SRH_NT_states_month_day
    
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
        if(nrow(subset(SRH_NT_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          temp_upstream_season0[i] <- 0
        } else {
          temp_upstream_season0[i] <- median(temp_data[subset(SRH_NT_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_upstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_NT_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          temp_downstream_season0[i] <- 0
        } else {
          temp_downstream_season0[i] <- median(temp_data[subset(SRH_NT_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        temp_downstream_season0[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_NT_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          temp_upstream_season1[i] <- 0
        } else {
          temp_upstream_season1[i] <- median(temp_data[subset(SRH_NT_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        temp_upstream_season1[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_NT_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          temp_downstream_season1[i] <- 0
        } else {
          temp_downstream_season1[i] <- median(temp_data[subset(SRH_NT_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
        if(nrow(subset(SRH_NT_states_month_day, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0[i] <- 0
        } else {
          spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(SRH_NT_states_month_day, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_NT_states_month_day, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0[i] <- 0
        } else {
          spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(SRH_NT_states_month_day, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_NT_states_month_day, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1[i] <- 0
        } else {
          spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(SRH_NT_states_month_day, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_NT_states_month_day, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1[i] <- 0
        } else {
          spillwindow_downstream_season1[i] <- median(spillwindow_data[subset(SRH_NT_states_month_day, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
    
    SRH_NT_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(fix = ifelse(yday(date_actual) > fix_start_date &
                            yday(date_actual) < fix_end_date, 1, 0)) %>% 
      mutate(fix = ifelse(is.na(fix), 0, fix)) -> SRH_NT_states_dates_fix
    
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
        if(nrow(subset(SRH_NT_states_dates_fix, state == i & direction == "upstream" & season == 0)) == 0){
          season0_upstream_fixprop[i] <- 0
        } else {
          season0_upstream_fixprop[i] <- sum(subset(SRH_NT_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
            length(subset(SRH_NT_states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
        }
      } else {
        season0_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 0)$fix)
      }
      # proportion of season1, upstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_NT_states_dates_fix, state == i & direction == "upstream" & season == 1)) == 0){
          season1_upstream_fixprop[i] <- 0
        } else {
          season1_upstream_fixprop[i] <- sum(subset(SRH_NT_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
            length(subset(SRH_NT_states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)
        }
      } else {
        season1_upstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "upstream" & season == 1 )$fix)
      }
      # proportion of season0, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_NT_states_dates_fix, state == i & direction == "downstream" & season == 0)) == 0){
          season0_downstream_fixprop[i] <- 0
        } else {
          season0_downstream_fixprop[i] <- sum(subset(SRH_NT_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
            length(subset(SRH_NT_states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
        }
      } else {
        season0_downstream_fixprop[i] <- sum(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)/
          length(subset(states_dates_fix, state == i  & direction == "downstream" & season == 0)$fix)
      }
      # proportion of season1, downstream that should be fixed to new conditions
      if(nrow(subset(states_dates_fix, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_NT_states_dates_fix, state == i & direction == "downstream" & season == 1)) == 0){
          season1_downstream_fixprop[i] <- 0
        } else {
          season1_downstream_fixprop[i] <- sum(subset(SRH_NT_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)/
            length(subset(SRH_NT_states_dates_fix, state == i  & direction == "downstream" & season == 1)$fix)
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
    
    SRH_NT_states_month_day %>% 
      left_join(to_fix, by = "state") %>% 
      mutate(drop_date = ifelse(yday(date_actual.x) > fix_start_date &
                                  yday(date_actual.x) < fix_end_date, 1, 0)) %>% 
      mutate(drop_date = ifelse(is.na(drop_date), 0, drop_date)) %>% 
      filter(drop_date == 0) -> SRH_NT_states_month_day_trimmed
    
    
    
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state, by direction and season
      # season 0
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)) == 0){
        if(nrow(subset(SRH_NT_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)) == 0){
          spillwindow_upstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(SRH_NT_states_month_day_trimmed, state == i & direction == "upstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 0)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)) == 0){
        if(nrow(subset(SRH_NT_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)) == 0){
          spillwindow_downstream_season0_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(SRH_NT_states_month_day_trimmed, state == i & direction == "downstream" & season == 0)$date_run_year,i])
        }
      } else {
        spillwindow_downstream_season0_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 0)$date_run_year,i])
      }
      # season 1
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)) == 0){
        if(nrow(subset(SRH_NT_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)) == 0){
          spillwindow_upstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(SRH_NT_states_month_day_trimmed, state == i & direction == "upstream" & season == 1)$date_run_year,i])
        }
      } else {
        spillwindow_upstream_season1_notfixed[i] <- median(spillwindow_data[subset(states_month_day_trimmed, state == i  & direction == "upstream" & season == 1)$date_run_year,i])
      }
      if(nrow(subset(states_month_day_trimmed, state == i  & direction == "downstream" & season == 1)) == 0){
        if(nrow(subset(SRH_NT_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)) == 0){
          spillwindow_downstream_season1_notfixed[i] <- 0
        } else {
          spillwindow_downstream_season1_notfixed[i] <- median(spillwindow_data[subset(SRH_NT_states_month_day_trimmed, state == i & direction == "downstream" & season == 1)$date_run_year,i])
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
    possible_movements <- SRH_NT_envir$data$movements[, "col"][SRH_NT_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    upstream_season_sim <- sample(c(0,1),1, prob = c(1-season_upstream[i], season_upstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      
      if (upstream_season_sim == 1){
        
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_NT[i,possible_movements[j],iter] +
                                                                    btemp1_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i] + 
                                                                    bspillwindow_array_SRH_NT[i,possible_movements[j],iter]*spillwindow_upstream_season1[i] + 
                                                                    bwinterspill_array_SRH_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp1xorigin1_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin1 +
                                                                    btemp1xorigin2_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin2 + 
                                                                    btemp1xorigin3_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin3 +
                                                                    btemp1xorigin4_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin4 +
                                                                    btemp1xorigin5_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season1[i]*origin5 +
                                                                    borigin1_array_SRH_NT[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH_NT[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH_NT[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH_NT[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH_NT[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_NT[i,possible_movements,iter] +
                    btemp1_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i] + 
                    bspillwindow_array_SRH_NT[i,possible_movements,iter]*spillwindow_upstream_season1[i] + 
                    bwinterspill_array_SRH_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season1[i]*origin5 +
                    borigin1_array_SRH_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_NT[i,possible_movements,iter]*origin5))
      } else {
        upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_NT[i,possible_movements[j],iter] +
                                                                    btemp0_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i] +
                                                                    bspillwindow_array_SRH_NT[i,possible_movements[j],iter]*spillwindow_upstream_season0[i] + 
                                                                    bwinterspill_array_SRH_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                    btemp0xorigin1_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin1 +
                                                                    btemp0xorigin2_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin2 +
                                                                    btemp0xorigin3_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin3 +
                                                                    btemp0xorigin4_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin4 +
                                                                    btemp0xorigin5_array_SRH_NT[i,possible_movements[j],iter]*temp_upstream_season0[i]*origin5 +
                                                                    borigin1_array_SRH_NT[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH_NT[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH_NT[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH_NT[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH_NT[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_NT[i,possible_movements,iter] +
                    btemp0_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i] +
                    bspillwindow_array_SRH_NT[i,possible_movements,iter]*spillwindow_upstream_season0[i] + 
                    bwinterspill_array_SRH_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH_NT[i,possible_movements,iter]*temp_upstream_season0[i]*origin5 +
                    borigin1_array_SRH_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_NT[i,possible_movements,iter]*origin5))
      }
      
      
      
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRH_NT_envir$data$movements[, "col"][SRH_NT_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 41)
    
    # make sure to drop all movements to downstream states (in DE years, these aren't possible)
    grep(" downstream", model_states) -> downstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% downstream_indices)]
    
    # toggle temp0 v temp1 based on resampling season
    downstream_season_sim <- sample(c(0,1),1, prob = c(1-season_downstream[i], season_downstream[i]))
    
    for (j in 1:length(possible_movements)){
      
      if (downstream_season_sim == 1){
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_NT[i,possible_movements[j],iter] +
                                                                      btemp1_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i] + 
                                                                      bspillwindow_array_SRH_NT[i,possible_movements[j],iter]*spillwindow_downstream_season1[i] + 
                                                                      bwinterspill_array_SRH_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp1xorigin1_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin1 +
                                                                      btemp1xorigin2_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin2 + 
                                                                      btemp1xorigin3_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin3 +
                                                                      btemp1xorigin4_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin4 +
                                                                      btemp1xorigin5_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season1[i]*origin5 +
                                                                      borigin1_array_SRH_NT[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH_NT[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH_NT[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH_NT[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH_NT[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_NT[i,possible_movements,iter] +
                    btemp1_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i] + 
                    bspillwindow_array_SRH_NT[i,possible_movements,iter]*spillwindow_downstream_season1[i] + 
                    bwinterspill_array_SRH_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp1xorigin1_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin1 +
                    btemp1xorigin2_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin2 + 
                    btemp1xorigin3_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin3 +
                    btemp1xorigin4_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin4 +
                    btemp1xorigin5_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season1[i]*origin5 +
                    borigin1_array_SRH_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_NT[i,possible_movements,iter]*origin5))
      } else {
        downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH_NT[i,possible_movements[j],iter] +
                                                                      btemp0_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i] +
                                                                      bspillwindow_array_SRH_NT[i,possible_movements[j],iter]*spillwindow_downstream_season0[i] + 
                                                                      bwinterspill_array_SRH_NT[i,possible_movements[j],iter]*winterspill[i] +
                                                                      btemp0xorigin1_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin1 +
                                                                      btemp0xorigin2_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin2 +
                                                                      btemp0xorigin3_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin3 +
                                                                      btemp0xorigin4_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin4 +
                                                                      btemp0xorigin5_array_SRH_NT[i,possible_movements[j],iter]*temp_downstream_season0[i]*origin5 +
                                                                      borigin1_array_SRH_NT[i,possible_movements[j],iter]*origin1 +
                                                                      borigin2_array_SRH_NT[i,possible_movements[j],iter]*origin2 +
                                                                      borigin3_array_SRH_NT[i,possible_movements[j],iter]*origin3 +
                                                                      borigin4_array_SRH_NT[i,possible_movements[j],iter]*origin4 +
                                                                      borigin5_array_SRH_NT[i,possible_movements[j],iter]*origin5)/
          sum(exp(b0_array_SRH_NT[i,possible_movements,iter] +
                    btemp0_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i] +
                    bspillwindow_array_SRH_NT[i,possible_movements,iter]*spillwindow_downstream_season0[i] + 
                    bwinterspill_array_SRH_NT[i,possible_movements,iter]*winterspill[i] +
                    btemp0xorigin1_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin1 +
                    btemp0xorigin2_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin2 +
                    btemp0xorigin3_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin3 +
                    btemp0xorigin4_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin4 +
                    btemp0xorigin5_array_SRH_NT[i,possible_movements,iter]*temp_downstream_season0[i]*origin5 +
                    borigin1_array_SRH_NT[i,possible_movements,iter]*origin1 +
                    borigin2_array_SRH_NT[i,possible_movements,iter]*origin2 +
                    borigin3_array_SRH_NT[i,possible_movements,iter]*origin3 +
                    borigin4_array_SRH_NT[i,possible_movements,iter]*origin4 +
                    borigin5_array_SRH_NT[i,possible_movements,iter]*origin5))
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



# ORDER THE STATES FOR PLOTTING
states_order_for_plot <- gsub(" Mouth| Upstream", "", model_states)
states_order_for_plot <- states_order_for_plot[!(duplicated(states_order_for_plot))]
# Make a couple of changes to make them be in the order from most downstream to most upstream
states_order_for_plot[10] <- "Fifteenmile Creek"
states_order_for_plot[11] <- "Deschutes River"
states_order_for_plot[12] <- "John Day River"
states_order_for_plot[14] <- "Walla Walla River"
states_order_for_plot[15] <- "Yakima River"
states_order_for_plot[18] <- "Methow River"
states_order_for_plot[19] <- "Okanogan River"

states_order_for_plot[15:28] <- states_order_for_plot[14:27]
states_order_for_plot[14] <- "BON to MCN other tributaries"
states_order_for_plot[22:28] <- states_order_for_plot[21:27]
states_order_for_plot[21] <- "Upstream WEL other tributaries"
states_order_for_plot[28] <- "loss"

states_order_for_plot[23] <- "Clearwater River"
states_order_for_plot[24] <- "Asotin Creek"
states_order_for_plot[25] <- "Grande Ronde River"
states_order_for_plot[26] <- "Salmon River"


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