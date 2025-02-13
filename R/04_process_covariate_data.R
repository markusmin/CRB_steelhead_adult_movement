#### 04: process covariate data

# load libraries
library(here)
library(janitor)
library(tidyverse)
library(lubridate)
library(car)
library(ggthemes)
library(zoo)
library(ggpubr)
library(MARSS)
library(mixtools)

#### define temperature functions ####

# Function 1: reformat covariate data
cov_reformat <- function(input_data, variable_name){
  
  # need to deal with leap years
  leap_years <- seq(1960, 2024, by = 4)
  
  input_data %>% 
    pivot_longer(., cols = colnames(.)[2:ncol(.)]) %>% 
    dplyr::rename(year = name,
                  !!quo_name(variable_name) := value) %>% 
    mutate(year = gsub("X","",year)) %>% 
    mutate(year = as.numeric(year), day = as.numeric(day)) %>% 
    # drop day 366 in non-leap years
    dplyr::filter(!(day == 366 & !(year %in% leap_years))) %>%
    # Add a new column for date (combine year and day)
    mutate(date = format(as.Date(day, origin = paste0(year-1, "-12-31")))) %>% 
    arrange(year, day) %>% 
    # drop all years before 2005
    filter(year >= 2005) -> output_data
  
  
  return(output_data)
}

# Function 2: Generate plot and regression of forebay vs. tailrace temperature  
forebay_v_tailrace <- function(forebay_t, tailrace_t, dam){
  # reformat into one df
  dplyr::select(forebay_t, c(temp, date)) %>% 
    dplyr::rename(forebay_temp = temp) -> forebay_t_forjoin
  
  tailrace_t %>% 
    dplyr::rename(tailrace_temp = temp) %>% 
    dplyr::select(tailrace_temp, date) %>% 
    left_join(., forebay_t_forjoin, by = "date") -> temp_comparison
  
  # get the lm fit
  lm_fit <- lm(forebay_temp~tailrace_temp, data = temp_comparison)
  
  # generate plot
  forebay_v_tailrace_plot <- ggplot(subset(temp_comparison, !is.na(tailrace_temp) & !(is.na(forebay_temp))), 
                                    aes(y = forebay_temp, x = tailrace_temp)) +
    geom_point() +
    ggtitle(dam) +
    scale_x_continuous(lim = c(0, 28)) +
    scale_y_continuous(lim = c(0, 28)) +
    annotate(geom = "text", x = 2, y = 20, hjust = 0,
             label = paste0("y = ", round(lm_fit$coefficients[1],3), " + ", round(lm_fit$coefficients[2],3)," * x"))
  
  # return plot and equation
  # return(list(forebay_v_tailrace_plot, lm_fit))
  return(forebay_v_tailrace_plot)
}

# Function 3: Plot ts of forebay and tailrace temperatures
plot_temp_ts <- function(forebay_t, tailrace_t, dam){
  forebay_t %>% 
    mutate(date = ymd(date)) %>% 
    mutate(day = yday(date)) %>% 
    mutate(month = as.factor(month(date))) -> forebay_t
  
  tailrace_t %>% 
    mutate(date = ymd(date)) %>% 
    mutate(day = yday(date)) %>% 
    mutate(month = as.factor(month(date))) -> tailrace_t
  
  forebay_t_ts <- ggplot(forebay_t, aes(x = day, y = temp, color = year)) +
    geom_point() +
    ggtitle(dam) +
    ylab("Forebay T")
  
  tailrace_t_ts <- ggplot(tailrace_t, aes(x = day, y = temp, color = year)) +
    geom_point() +
    ggtitle(dam) +
    ylab("Tailrace T")
  
  two_plots <- ggarrange(forebay_t_ts, tailrace_t_ts, ncol = 1)
  
  return(two_plots)
}

# Function 4: Filter out the outliers based on +/- 4 degrees, as well as based on 
# a moving average of seven days, and if it's more than 2 degrees outside that moving average, drop it
# issue with this is that we get runs of outliers that aren't dropped by that alone, which is why we need the first step to take out major outliers
filter_outliers_temp <- function(forebay_t, tailrace_t, dam){
  forebay_t %>% 
    mutate(date = ymd(date)) %>% 
    mutate(day = yday(date)) %>% 
    mutate(month = as.factor(month(date))) -> forebay_t
  
  # get a moving average using rollapply()
  # forebay_t %>% 
  #   mutate(MA = rollapply(temp, 7 , mean, na.rm = TRUE, align='center', fill = NA)) %>% 
  #   # if it's the beginning or end of the TS, then MA will be NA, so replace those with the temperature (so these won't be dropped)
  #   mutate(MA = ifelse(is.na(MA), temp, MA)) %>% 
  #   filter(abs(temp - MA) <= 2) %>% 
  #   dplyr::select(date, temp) %>%
  #   dplyr::rename(forebay_temp = temp)-> forebay_t_filtered
  
  forebay_t %>%
    group_by(day) %>%
    summarise(day_mean = mean(temp, na.rm = TRUE)) -> forebay_t_daily_means
  
  forebay_t %>%
    # filter based on major outliers
    left_join(., forebay_t_daily_means, by = "day") %>%
    filter(abs(temp - day_mean) <= 4) %>%
    # filter based on moving average
    mutate(MA = rollapply(temp, 7 , mean, na.rm = TRUE, align='center', fill = NA)) %>%
    # if it's the beginning or end of the TS, then MA will be NA, so replace those with the temperature (so these won't be dropped)
    mutate(MA = ifelse(is.na(MA), temp, MA)) %>%
    filter(abs(temp - MA) <= 2) %>%
    dplyr::select(date, temp) %>%
    dplyr::rename(forebay_temp = temp)-> forebay_t_filtered
  
  tailrace_t %>% 
    mutate(date = ymd(date)) %>% 
    mutate(day = yday(date)) %>% 
    mutate(month = as.factor(month(date))) -> tailrace_t
  
  # tailrace_t %>% 
  #   mutate(MA = rollapply(temp, 7 , mean, na.rm = TRUE, align='center', fill = NA)) %>% 
  #   # if it's the beginning or end of the TS, then MA will be NA, so replace those with the temperature (so these won't be dropped)
  #   mutate(MA = ifelse(is.na(MA), temp, MA)) %>% 
  #   filter(abs(temp - MA) <= 2) %>% 
  #       dplyr::select(date, temp) %>%
  #   dplyr::rename(tailrace_temp = temp)-> tailrace_t_filtered
  
  tailrace_t %>%
    group_by(day) %>%
    summarise(day_mean = mean(temp, na.rm = TRUE)) -> tailrace_t_daily_means
  
  tailrace_t %>%
    left_join(., tailrace_t_daily_means, by = "day") %>%
    filter(abs(temp - day_mean) <= 4) %>%
    # filter based on moving average
    mutate(MA = rollapply(temp, 7 , mean, na.rm = TRUE, align='center', fill = NA)) %>%
    # if it's the beginning or end of the TS, then MA will be NA, so replace those with the temperature (so these won't be dropped)
    mutate(MA = ifelse(is.na(MA), temp, MA)) %>%
    filter(abs(temp - MA) <= 2) %>%
    dplyr::select(date, temp) %>%
    dplyr::rename(tailrace_temp = temp)-> tailrace_t_filtered
  
  return(list(forebay_t_filtered, tailrace_t_filtered))
}

#### Load and reformat the temperature data ####

# BON
BON_temp_t <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "BON", "BON_tailrace_temp_2025-02-03.csv"))[1:366,]
BON_temp_t_long <- cov_reformat(input_data = BON_temp_t, variable_name = "temp")
BON_temp_f <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "BON", "BON_temp_2025-02-03.csv"))[1:366,]
BON_temp_f_long <- cov_reformat(input_data = BON_temp_f, variable_name = "temp")

# ICH
ICH_temp_t <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "ICH", "ICH_tailrace_temp_2025-02-03.csv"))[1:366,]
ICH_temp_t_long <- cov_reformat(input_data = ICH_temp_t, variable_name = "temp")
ICH_temp_f <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "ICH", "ICH_temp_2025-02-03.csv"))[1:366,]
ICH_temp_f_long <- cov_reformat(input_data = ICH_temp_f, variable_name = "temp")

# LGR
LGR_temp_t <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "LGR", "LGR_tailrace_temp_2025-02-03.csv"))[1:366,]
LGR_temp_t_long <- cov_reformat(input_data = LGR_temp_t, variable_name = "temp")
LGR_temp_f <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "LGR", "LGR_temp_2025-02-03.csv"))[1:366,]
LGR_temp_f_long <- cov_reformat(input_data = LGR_temp_f, variable_name = "temp")

# MCN
MCN_temp_t <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "MCN", "MCN_tailrace_temp_2025-02-03.csv"))[1:366,]
MCN_temp_t_long <- cov_reformat(input_data = MCN_temp_t, variable_name = "temp")
MCN_temp_f <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "MCN", "MCN_temp_2025-02-03.csv"))[1:366,]
MCN_temp_f_long <- cov_reformat(input_data = MCN_temp_f, variable_name = "temp")


# PRA
PRA_temp_t <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "PRA", "PRA_tailrace_temp_2025-02-03.csv"))[1:366,]
PRA_temp_t_long <- cov_reformat(input_data = PRA_temp_t, variable_name = "temp")
PRA_temp_f <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "PRA", "PRA_temp_2025-02-03.csv"))[1:366,]
PRA_temp_f_long <- cov_reformat(input_data = PRA_temp_f, variable_name = "temp")


# RIS
RIS_temp_t <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "RIS", "RIS_tailrace_temp_2025-02-03.csv"))[1:366,]
RIS_temp_t_long <- cov_reformat(input_data = RIS_temp_t, variable_name = "temp")
RIS_temp_f <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "RIS", "RIS_temp_2025-02-03.csv"))[1:366,]
RIS_temp_f_long <- cov_reformat(input_data = RIS_temp_f, variable_name = "temp")


# RRE
RRE_temp_t <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "RRE", "RRE_tailrace_temp_2025-02-03.csv"))[1:366,]
RRE_temp_t_long <- cov_reformat(input_data = RRE_temp_t, variable_name = "temp")
RRE_temp_f <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "RRE", "RRE_temp_2025-02-03.csv"))[1:366,]
RRE_temp_f_long <- cov_reformat(input_data = RRE_temp_f, variable_name = "temp")

# WEL
WEL_temp_t <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "WEL", "WEL_tailrace_temp_2025-02-03.csv"))[1:366,]
WEL_temp_t_long <- cov_reformat(input_data = WEL_temp_t, variable_name = "temp")
WEL_temp_f <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "WEL", "WEL_temp_2025-02-03.csv"))[1:366,]
WEL_temp_f_long <- cov_reformat(input_data = WEL_temp_f, variable_name = "temp")






#### Visualize the data ####

# Plot forebay vs. tailrace temperatures collected on the same day
forebay_v_tailrace(forebay_t = BON_temp_f_long, tailrace_t = BON_temp_t_long, dam = "BON")
forebay_v_tailrace(forebay_t = ICH_temp_f_long, tailrace_t = ICH_temp_t_long, dam = "ICH")
forebay_v_tailrace(forebay_t = LGR_temp_f_long, tailrace_t = LGR_temp_t_long, dam = "LGR")
forebay_v_tailrace(forebay_t = MCN_temp_f_long, tailrace_t = MCN_temp_t_long, dam = "MCN")
forebay_v_tailrace(forebay_t = PRA_temp_f_long, tailrace_t = PRA_temp_t_long, dam = "PRA")
forebay_v_tailrace(forebay_t = RIS_temp_f_long, tailrace_t = RIS_temp_t_long, dam = "RIS")
forebay_v_tailrace(forebay_t = RRE_temp_f_long, tailrace_t = RRE_temp_t_long, dam = "RRE")
forebay_v_tailrace(forebay_t = WEL_temp_f_long, tailrace_t = WEL_temp_t_long, dam = "WEL")

# there are some complete nonsense values that it looks like we will have to filter. 
# Worst offenders are RIS, RRE, and WEL. Is it the forebay or tailrace temperatures that are bad? 

# Plot time series of temperature throughout the year
plot_temp_ts(forebay_t = BON_temp_f_long, tailrace_t = BON_temp_t_long, dam = "BON")
plot_temp_ts(forebay_t = ICH_temp_f_long, tailrace_t = ICH_temp_t_long, dam = "ICH")
plot_temp_ts(forebay_t = LGR_temp_f_long, tailrace_t = LGR_temp_t_long, dam = "LGR")
plot_temp_ts(forebay_t = MCN_temp_f_long, tailrace_t = MCN_temp_t_long, dam = "MCN")
plot_temp_ts(forebay_t = PRA_temp_f_long, tailrace_t = PRA_temp_t_long, dam = "PRA")
plot_temp_ts(forebay_t = RIS_temp_f_long, tailrace_t = RIS_temp_t_long, dam = "RIS")
plot_temp_ts(forebay_t = RRE_temp_f_long, tailrace_t = RRE_temp_t_long, dam = "RRE")
plot_temp_ts(forebay_t = WEL_temp_f_long, tailrace_t = WEL_temp_t_long, dam = "WEL")

# There are some clear outliers in both the forebay and tailrace temperatures that need to be filtered

#### Filter outliers, based on manual inspection ####

### For total nonsense runs of values (RIS, RRE, and WEL) - manually select and drop those days

## Wells Dam
# WEL, forebay
WEL_temp_f_long %>% 
  mutate(date = ymd(date)) %>% 
  mutate(day = yday(date)) %>% 
  mutate(year = as.factor(year)) -> WEL_temp_f_long_2

ggplot(WEL_temp_f_long_2, aes(x = day, y = temp, color = year)) +
  geom_point() +
  scale_color_tableau(palette = "Tableau 20") +
  ggtitle("Wells Forebay T, Original")

# run of values in 2009
# subset(WEL_temp_f_long, year == 2009)
# nonsense starts on day 270 and end on day 300

WEL_temp_f_long %>% 
  filter(!(year == 2009 & day >= 270 & day <=300)) -> WEL_temp_f_long_clean

ggplot(WEL_temp_f_long_clean, aes(x = day, y = temp, color = as.factor(year))) +
  geom_point() +
  scale_color_tableau(palette = "Tableau 20") +
  ggtitle("Wells Forebay T, Clean")

# WEL, tailrace
WEL_temp_t_long %>% 
  mutate(date = ymd(date)) %>% 
  mutate(day = yday(date)) %>% 
  mutate(year = as.factor(year)) -> WEL_temp_t_long_2

ggplot(WEL_temp_t_long_2, aes(x = day, y = temp, color = year)) +
  geom_point() +
  scale_color_tableau(palette = "Tableau 20") +
  ggtitle("Wells Tailrace T, Original")

# run of values in 2007 and run of values in 2008 and run of values in 2009 and run of values in 2017
# subset(WEL_temp_t_long, year == 2007)
# 2007: nonsense starts on day 161 and ends on day 260

# subset(WEL_temp_t_long, year == 2008)
# 2008: nonsense starts on day 79 and ends on day 182

# subset(WEL_temp_t_long, year == 2009)
# 2009: nonsense starts on day 272 and ends on day 300

# subset(WEL_temp_t_long, year == 2017)
# 2017: nonsense starts on day 177 and ends on day 197

WEL_temp_t_long %>% 
  filter(!(year == 2007 & day >= 161 & day <=260)) %>% 
  filter(!(year == 2008 & day >= 79 & day <=182)) %>% 
  filter(!(year == 2009 & day >= 272 & day <=300)) %>% 
  filter(!(year == 2017 & day >= 177 & day <=197)) -> WEL_temp_t_long_clean

ggplot(WEL_temp_t_long_clean, aes(x = day, y = temp, color = as.factor(year))) +
  geom_point() +
  scale_color_tableau(palette = "Tableau 20") +
  ggtitle("Wells Tailrace T, Clean")

## Rock Island Dam

# RIS, forebay
RIS_temp_f_long %>% 
  mutate(date = ymd(date)) %>% 
  mutate(day = yday(date)) %>% 
  mutate(year = as.factor(year)) -> RIS_temp_f_long_2

ggplot(RIS_temp_f_long_2, aes(x = day, y = temp, color = year)) +
  geom_point() +
  scale_color_tableau(palette = "Tableau 20") +
  ggtitle("RIS Forebay T, Original")

# so there are three flat runs, two in 2005 and one in 2006. Let's filter that out based on lag/lead of at least three of the same value in a row
# there's also a stupid run in 2007, drop that: day 161-281
# subset(RIS_temp_f_long, year == 2007)

RIS_temp_f_long %>% 
  filter(!(year == 2007 & day >= 161 & day <=281)) %>% 
  filter(!(temp == lag(temp) & temp == lag(temp, 2) |
             temp == lead(temp) & temp == lead(temp, 2))) -> RIS_temp_f_long_clean

ggplot(RIS_temp_f_long_clean, aes(x = day, y = temp, color = as.factor(year))) +
  geom_point() +
  scale_color_tableau(palette = "Tableau 20") +
  ggtitle("RIS Forebay T, Clean")

# RIS, tailrace
RIS_temp_t_long %>% 
  mutate(date = ymd(date)) %>% 
  mutate(day = yday(date)) %>% 
  mutate(year = as.factor(year)) -> RIS_temp_t_long_2

ggplot(RIS_temp_t_long_2, aes(x = day, y = temp, color = year)) +
  geom_point() +
  scale_color_tableau(palette = "Tableau 20") +
  ggtitle("RIS Tailrace T, Original")
# just looks like two flat runs, one in 2005 and one in 2006

# there's also a single value in 2006 that's totally messing up the linear interpolation that we need to drop on day 101
# subset(RIS_temp_t_long, year == 2006)

RIS_temp_t_long %>% 
  filter(!(year == 2006 & day == 101)) %>% 
  filter(!(temp == lag(temp) & temp == lag(temp, 2) |
             temp == lead(temp) & temp == lead(temp, 2))) -> RIS_temp_t_long_clean

ggplot(RIS_temp_t_long_clean, aes(x = day, y = temp, color = as.factor(year))) +
  geom_point() +
  scale_color_tableau(palette = "Tableau 20") +
  ggtitle("RIS Tailrace T, Clean")

## Rocky Reach Dam

# RRE, forebay
RRE_temp_f_long %>% 
  mutate(date = ymd(date)) %>% 
  mutate(day = yday(date)) %>% 
  mutate(year = as.factor(year)) -> RRE_temp_f_long_2

ggplot(RRE_temp_f_long_2, aes(x = day, y = temp, color = year)) +
  geom_point() +
  scale_color_tableau(palette = "Tableau 20") +
  ggtitle("RRE Forebay T, Original")

# so there are four flat runs, two in 2005 and two in 2006. Let's filter that out based on lag/lead of at least three of the same value in a row

RRE_temp_f_long %>% 
  filter(!(temp == lag(temp) & temp == lag(temp, 2) |
             temp == lead(temp) & temp == lead(temp, 2))) -> RRE_temp_f_long_clean

ggplot(RRE_temp_f_long_clean, aes(x = day, y = temp, color = as.factor(year))) +
  geom_point() +
  scale_color_tableau(palette = "Tableau 20") +
  ggtitle("RRE Forebay T, Clean")

# RRE, tailrace
RRE_temp_t_long %>% 
  mutate(date = ymd(date)) %>% 
  mutate(day = yday(date)) %>% 
  mutate(year = as.factor(year)) -> RRE_temp_t_long_2

ggplot(RRE_temp_t_long_2, aes(x = day, y = temp, color = year)) +
  geom_point() +
  scale_color_tableau(palette = "Tableau 20") +
  ggtitle("RRE Tailrace T, Original")
# just looks like three flat runs, two in 2005 and one in 2006

RRE_temp_t_long %>% 
  filter(!(temp == lag(temp) & temp == lag(temp, 2) |
             temp == lead(temp) & temp == lead(temp, 2))) -> RRE_temp_t_long_clean

ggplot(RRE_temp_t_long_clean, aes(x = day, y = temp, color = as.factor(year))) +
  geom_point() +
  scale_color_tableau(palette = "Tableau 20") +
  ggtitle("RRE Tailrace T, Clean")

#### Filter outliers, based on filtering algorithm ####
BON_temps_clean <- filter_outliers_temp(forebay_t = BON_temp_f_long, tailrace_t = BON_temp_t_long, dam = "BON")
ICH_temps_clean <- filter_outliers_temp(forebay_t = ICH_temp_f_long, tailrace_t = ICH_temp_t_long, dam = "ICH")
LGR_temps_clean <- filter_outliers_temp(forebay_t = LGR_temp_f_long, tailrace_t = LGR_temp_t_long, dam = "LGR")
MCN_temps_clean <- filter_outliers_temp(forebay_t = MCN_temp_f_long, tailrace_t = MCN_temp_t_long, dam = "MCN")
PRA_temps_clean <- filter_outliers_temp(forebay_t = PRA_temp_f_long, tailrace_t = PRA_temp_t_long, dam = "PRA")
RIS_temps_clean <- filter_outliers_temp(forebay_t = RIS_temp_f_long_clean, tailrace_t = RIS_temp_t_long_clean, dam = "RIS")
RRE_temps_clean <- filter_outliers_temp(forebay_t = RRE_temp_f_long_clean, tailrace_t = RRE_temp_t_long_clean, dam = "RRE")
WEL_temps_clean <- filter_outliers_temp(forebay_t = WEL_temp_f_long_clean, tailrace_t = WEL_temp_t_long_clean, dam = "WEL")

#### Reformat and export temperature ####
reformat_temp_for_export <- function(temps_clean){
  # Generate the full sequence of dates for which we need temperatures
  temperature_df_skeleton <- data.frame(date = seq(ymd('2005-06-01'),ymd('2024-12-31'), by = 'days'))
  
  temperature_df_skeleton %>% 
    left_join(., temps_clean[[2]], by = "date") %>% 
    left_join(., temps_clean[[1]], by = "date") -> temp_ts
  
  return(temp_ts)
}


# Keep only 2005-06-01 and onwards

# TO DO: Need to complete these to make sure we aren't missing any date/times

BON_temp_ts <- reformat_temp_for_export(BON_temps_clean)
ICH_temp_ts <- reformat_temp_for_export(ICH_temps_clean)
LGR_temp_ts <- reformat_temp_for_export(LGR_temps_clean)
MCN_temp_ts <- reformat_temp_for_export(MCN_temps_clean)
PRA_temp_ts <- reformat_temp_for_export(PRA_temps_clean)
RIS_temp_ts <- reformat_temp_for_export(RIS_temps_clean)
RRE_temp_ts <- reformat_temp_for_export(RRE_temps_clean)
WEL_temp_ts <- reformat_temp_for_export(WEL_temps_clean)





write.csv(BON_temp_ts, here::here("Data", "covariate_data", "temp_processed", "BON_temp_processed.csv"), row.names = FALSE)
write.csv(MCN_temp_ts, here::here("Data", "covariate_data", "temp_processed", "MCN_temp_processed.csv"), row.names = FALSE)
write.csv(PRA_temp_ts, here::here("Data", "covariate_data", "temp_processed", "PRA_temp_processed.csv"), row.names = FALSE)
write.csv(RIS_temp_ts, here::here("Data", "covariate_data", "temp_processed", "RIS_temp_processed.csv"), row.names = FALSE)
write.csv(RRE_temp_ts, here::here("Data", "covariate_data", "temp_processed", "RRE_temp_processed.csv"), row.names = FALSE)
write.csv(ICH_temp_ts, here::here("Data", "covariate_data", "temp_processed", "ICH_temp_processed.csv"), row.names = FALSE)
write.csv(WEL_temp_ts, here::here("Data", "covariate_data", "temp_processed", "WEL_temp_processed.csv"), row.names = FALSE)
write.csv(LGR_temp_ts, here::here("Data", "covariate_data", "temp_processed", "LGR_temp_processed.csv"), row.names = FALSE)






#### Fit MARSS temperature model ####
# This section of code fits a MAR model to temperature data from different dams
# to estimate temperature for every date in our time series.

# The MAR model that we will fit is a state-space model, with forebay and tailrace
# temperatures as observations of a "Columbia River Basin" temperature, with
# each dam estimated as an offset of it
# We will have sixteen observations of the same process (forebay and tailrace temperatures for eight dams)
BON_temp <- BON_temp_ts
MCN_temp <- MCN_temp_ts
PRA_temp <- PRA_temp_ts
RIS_temp <- RIS_temp_ts
RRE_temp <- RRE_temp_ts
ICH_temp <- ICH_temp_ts
WEL_temp <- WEL_temp_ts
LGR_temp <- LGR_temp_ts

# make the column data names more informative by including dam
colnames(BON_temp) <- paste0("BON_", colnames(BON_temp))
colnames(MCN_temp) <- paste0("MCN_", colnames(MCN_temp))
colnames(PRA_temp) <- paste0("PRA_", colnames(PRA_temp))
colnames(RIS_temp) <- paste0("RIS_", colnames(RIS_temp))
colnames(RRE_temp) <- paste0("RRE_", colnames(RRE_temp))
colnames(ICH_temp) <- paste0("ICH_", colnames(ICH_temp))
colnames(WEL_temp) <- paste0("WEL_", colnames(WEL_temp))
colnames(LGR_temp) <- paste0("LGR_", colnames(LGR_temp))


# reformat data for model - drop date and transpose, bind all dams together
t(BON_temp[,2:3]) %>% 
  rbind(., t(MCN_temp[,2:3])) %>% 
  rbind(., t(PRA_temp[,2:3])) %>% 
  rbind(., t(RIS_temp[,2:3])) %>% 
  rbind(., t(RRE_temp[,2:3])) %>% 
  rbind(., t(ICH_temp[,2:3])) %>% 
  rbind(., t(WEL_temp[,2:3])) %>% 
  rbind(., t(LGR_temp[,2:3])) -> temp_for_MAR

dat = temp_for_MAR
# fit the model


# A are the offsets/biases for each dam - eight elements, each repeated twice
# BON is the reference
a <- matrix(data = as.list(rep(0, 16)),
            nrow = 16, ncol = 1)
a[3:16] <- c("a_MCN", "a_MCN", "a_PRA", "a_PRA", "a_RIS", "a_RIS",
             "a_RRE", "a_RRE", "a_ICH", "a_ICH", "a_WEL", "a_WEL", "a_LGR", "a_LGR")
A <- a

# no bias term for Columbia river temp
U <- "zero"

# just one process and no interactions
B <- matrix(1)

# Z is sixteen observations of the same process
Z <- matrix(1, 16, 1)

# we need eight elements to correspond to the 
# eight dams, and each element will be represented twice (once for forebay and once for tailrace)
r <- matrix(data = as.list(rep(0, 256)), nrow = 16, ncol = 16)
diag(r) <- c("r_BON", "r_BON", "r_MCN", "r_MCN", "r_PRA", "r_PRA", "r_RIS", "r_RIS",
             "r_RRE", "r_RRE", "r_ICH", "r_ICH", "r_WEL", "r_WEL", "r_LGR", "r_LGR")

# variance among dams is probably the same, given all the same equipment
R <- "diagonal and equal"
Q <- matrix("q") # this matrix is a 1x1 matrix, because there's only one trend we're estimating
model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)
kem <- MARSS(dat, model = model.list)


# some have considerably larger r (variance) values. Why?
# is it because of the difference between forebay and tailrace temperatures?
BON_temp[,2] - BON_temp[,3] -> BON_fore_tail_diff
summary(BON_fore_tail_diff)
MCN_temp[,2] - MCN_temp[,3] -> MCN_fore_tail_diff
summary(MCN_fore_tail_diff)
PRA_temp[,2] - PRA_temp[,3] -> PRA_fore_tail_diff
summary(PRA_fore_tail_diff)
RIS_temp[,2] - RIS_temp[,3] -> RIS_fore_tail_diff
summary(RIS_fore_tail_diff)
RRE_temp[,2] - RRE_temp[,3] -> RRE_fore_tail_diff
summary(RRE_fore_tail_diff)
ICH_temp[,2] - ICH_temp[,3] -> ICH_fore_tail_diff
summary(ICH_fore_tail_diff)
WEL_temp[,2] - WEL_temp[,3] -> WEL_fore_tail_diff
summary(WEL_fore_tail_diff)
LGR_temp[,2] - LGR_temp[,3] -> LGR_fore_tail_diff
summary(LGR_fore_tail_diff)

# plot Columbia River temperature
dates <- seq(ymd('2005-06-01'),ymd('2024-12-31'), by = 'days')

par(mfrow = c(1,1))


plot(dates, kem$states, ylab = "Columbia River Temperature", 
     xlab = "", type = "l")
lines(dates, kem$states - 1.96 * kem$states.se, type = "l", lwd = 1, lty = 2, col = "red")
lines(dates, kem$states + 1.96 * kem$states.se, type = "l", lwd = 1, lty = 2, col = "red")
title("Columbia River Temperature")


# extract the y hat values
hat_yt <- MARSShatyt(kem)

# plot the yhat for each dam
# add points for actual data

par(mfrow = c(2,2))
for (i in 1:16) {
  plot(dates, kem$ytT[i, ], ylab = "Temperature", 
       xlab = "", type = "l")
  lines(dates, kem$ytT[i, ] - 1.96 * kem$ytT.se[i, 
  ], type = "l", lwd = 1, lty = 2, col = "red")
  lines(dates, kem$ytT[i, ] + 1.96 * kem$ytT.se[i, 
  ], type = "l", lwd = 1, lty = 2, col = "red")
  # points(dates, dat[i, ], cex = 0.2)
  title(rownames(dat)[i])
}

# plot just the data
par(mfrow = c(1,1))
for (i in 1:16) {
  plot(dates, dat[i, ], cex = 0.4)
  title(rownames(dat)[i])
}

# export the yhat values for tailrace temperatures for modeling
temp_yhat <- data.frame(date = dates,
                        BON = kem$ytT[1,],
                        MCN = kem$ytT[3,],
                        PRA = kem$ytT[5,],
                        RIS = kem$ytT[7,],
                        RRE = kem$ytT[9,],
                        ICH = kem$ytT[11,],
                        WEL = kem$ytT[13,],
                        LGR = kem$ytT[15,])


# extract A values and add to the state
a_est <- as.data.frame(coef(kem)$A)
a_MCN <- a_est["a_MCN",]
a_PRA <- a_est["a_PRA",]
a_RIS <- a_est["a_RIS",]
a_RRE <- a_est["a_RRE",]
a_ICH <- a_est["a_ICH",]
a_WEL <- a_est["a_WEL",]
a_LGR <- a_est["a_LGR",]


# BON estimate is just the state
BON_temp_model_est <- kem$states

# every other is just plus the bias term
MCN_temp_model_est <- kem$states + a_MCN
PRA_temp_model_est <- kem$states + a_PRA
RIS_temp_model_est <- kem$states + a_RIS
RRE_temp_model_est <- kem$states + a_RRE
ICH_temp_model_est <- kem$states + a_ICH
WEL_temp_model_est <- kem$states + a_WEL
LGR_temp_model_est <- kem$states + a_LGR

temp_model_est <- data.frame(date = dates,
                             MOUTH = rep(NA, 7154),
                             BON = BON_temp_model_est[1,],
                             MCN = MCN_temp_model_est[1,],
                             PRA = PRA_temp_model_est[1,],
                             RIS = RIS_temp_model_est[1,],
                             RRE = RRE_temp_model_est[1,],
                             ICH = ICH_temp_model_est[1,],
                             WEL = WEL_temp_model_est[1,],
                             LGR = LGR_temp_model_est[1,])

# plot fits to data
BON_temp %>% 
  pivot_longer(cols = c(BON_tailrace_temp, BON_forebay_temp)) -> BON_temp_2

ggplot(temp_model_est, aes(x = dates, y = BON)) +
  geom_point(data = BON_temp_2, aes(x = ymd(BON_date), y = value, color = name)) +
  geom_line() +
  ggtitle("Modeled temperature at Bonneville Dam")

MCN_temp %>% 
  pivot_longer(cols = c(MCN_tailrace_temp, MCN_forebay_temp)) -> MCN_temp_2

ggplot(temp_model_est, aes(x = dates, y = MCN)) +
  geom_point(data = MCN_temp_2, aes(x = ymd(MCN_date), y = value, color = name)) +
  geom_line() +
  ggtitle("Modeled temperature at McNary Dam")

PRA_temp %>% 
  pivot_longer(cols = c(PRA_tailrace_temp, PRA_forebay_temp)) -> PRA_temp_2

ggplot(temp_model_est, aes(x = dates, y = PRA)) +
  geom_point(data = PRA_temp_2, aes(x = ymd(PRA_date), y = value, color = name)) +
  geom_line() +
  ggtitle("Modeled temperature at Priest Rapids Dam")

RIS_temp %>% 
  pivot_longer(cols = c(RIS_tailrace_temp, RIS_forebay_temp)) -> RIS_temp_2

ggplot(temp_model_est, aes(x = dates, y = RIS)) +
  geom_point(data = RIS_temp_2, aes(x = ymd(RIS_date), y = value, color = name)) +
  geom_line() +
  ggtitle("Modeled temperature at Rock Island Dam")

RRE_temp %>% 
  pivot_longer(cols = c(RRE_tailrace_temp, RRE_forebay_temp)) -> RRE_temp_2

ggplot(temp_model_est, aes(x = dates, y = RRE)) +
  geom_point(data = RRE_temp_2, aes(x = ymd(RRE_date), y = value, color = name)) +
  geom_line() +
  ggtitle("Modeled temperature at Rocky Reach Dam")

ICH_temp %>% 
  pivot_longer(cols = c(ICH_tailrace_temp, ICH_forebay_temp)) -> ICH_temp_2

ggplot(temp_model_est, aes(x = dates, y = ICH)) +
  geom_point(data = ICH_temp_2, aes(x = ymd(ICH_date), y = value, color = name)) +
  geom_line() +
  ggtitle("Modeled temperature at Ice Harbor Dam")

WEL_temp %>% 
  pivot_longer(cols = c(WEL_tailrace_temp, WEL_forebay_temp)) -> WEL_temp_2

ggplot(temp_model_est, aes(x = dates, y = WEL)) +
  geom_point(data = WEL_temp_2, aes(x = ymd(WEL_date), y = value, color = name)) +
  geom_line() +
  ggtitle("Modeled temperature at Wells Dam")

LGR_temp %>% 
  pivot_longer(cols = c(LGR_tailrace_temp, LGR_forebay_temp)) -> LGR_temp_2

ggplot(temp_model_est, aes(x = dates, y = LGR)) +
  geom_point(data = LGR_temp_2, aes(x = ymd(LGR_date), y = value, color = name)) +
  geom_line() +
  ggtitle("Modeled temperature at Lower Granite Dam")


write.csv(temp_model_est, here::here("Data", "covariate_data", "temp_processed", "temp_mod_est.csv"), row.names = FALSE)










#### Calculate temperature across windows experienced by fish ####
temp_mod_est <- temp_model_est
# Get the temperature data
temp_mod_est$date <- ymd(temp_mod_est$date)

# extract each dam
BON_temp_ts <- dplyr::select(temp_mod_est, date, BON)
MCN_temp_ts <- dplyr::select(temp_mod_est, date, MCN)
PRA_temp_ts <- dplyr::select(temp_mod_est, date, PRA)
RIS_temp_ts <- dplyr::select(temp_mod_est, date, RIS)
RRE_temp_ts <- dplyr::select(temp_mod_est, date, RRE)
WEL_temp_ts <- dplyr::select(temp_mod_est, date, WEL)
ICH_temp_ts <- dplyr::select(temp_mod_est, date, ICH)
LGR_temp_ts <- dplyr::select(temp_mod_est, date, LGR)

## Load fish movement data to calculate windows
snake_adults_states_complete <- read.csv(here::here("intermediate_outputs", "adults_states_complete", "snake_adults_states_complete.csv"))
midcol_adults_states_complete <- read.csv(here::here("intermediate_outputs", "adults_states_complete", "middle_columbia_adults_states_complete.csv"))
uppcol_adults_states_complete <- read.csv(here::here("intermediate_outputs", "adults_states_complete", "upper_columbia_adults_states_complete.csv"))

# combine all
snake_adults_states_complete %>% 
  bind_rows(., midcol_adults_states_complete) %>% 
  bind_rows(., uppcol_adults_states_complete) -> ASC

# now add tag code metadata, for natal origins
origin_numeric <- data.frame(natal_origin = c("Asotin_Creek", 
                                              "Clearwater_River",
                                              "Deschutes_River", 
                                              "Entiat_River", 
                                              "Fifteenmile_Creek", 
                                              "Grande_Ronde_River", 
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

origin_rear_actual <- read.csv(here::here("Data", "covariate_data", "origin_rear_actual.csv"), row.names = 1)
origin_rear_actual %>% 
  dplyr::rename(natal_origin_numeric = natal_origin) %>% 
  left_join(origin_numeric, by = "natal_origin_numeric") %>% 
  dplyr::select(tag_code_2, natal_origin) -> tag_code_origins

ASC %>% 
  left_join(., tag_code_origins, by = "tag_code_2") -> ASC

# also note which ESU they're from, for plotting
ESU_origins <- data.frame(natal_origin = unique(ASC$natal_origin),
                          ESU = c(rep("Snake", 6),
                                  rep("Middle Columbia", 6),
                                  rep("Upper Columbia", 4)))

ASC %>% 
  left_join(., ESU_origins, by = "natal_origin") -> ASC

# define a function that takes just one state, and finds the average amount of time to move out of that state

residence_time <- function(residence_state){
  
  # don't keep any with implicit time interpolated
  ASC %>% 
    filter(state == residence_state & tag_code_2 == lead(tag_code_2) & pathway != "implicit" & lead(pathway) != "implicit" |
             lag(state) == residence_state & tag_code_2 == lag(tag_code_2) & pathway != "implicit" & lag(pathway) != "implicit") %>% 
    mutate(date_time = ymd_hms(date_time)) -> one_state_df
  
  # make a data frame to record all of the transitions
  n_transitions <- nrow(one_state_df)/2
  passage_df <- data.frame(tag_code_2 = one_state_df$tag_code_2[seq(1,(nrow(one_state_df)-1),2)],
                           ESU = one_state_df$ESU[seq(1,(nrow(one_state_df)-1),2)],
                           natal_origin = one_state_df$natal_origin[seq(1,(nrow(one_state_df)-1),2)],
                           passage_time = NA)
  
  
  for (i in 1:nrow(passage_df)){
    passage_df$passage_time[i] <- one_state_df$date_time[(i*2)] - one_state_df$date_time[(i*2-1)]
    
  }
  
  return(passage_df)
  
}

## Calculate median residence time in each mainstem state

main_mainstem_states <- c(
  "mainstem, BON to MCN",
  "mainstem, MCN to ICH or PRA",
  "mainstem, PRA to RIS",
  "mainstem, RIS to RRE",
  "mainstem, RRE to WEL",
  "mainstem, ICH to LGR",
  "mainstem, upstream of WEL",
  "mainstem, upstream of LGR")

list_quantiles <- list()

for (i in 1:length(main_mainstem_states)) {
  residence_df <- residence_time(residence_state = main_mainstem_states[i])
  
  # Look into windows - show: 1) median of the whole distribution, 
  # and what cutoff you'd need to capture 50%, 80%, and 95% of fish residence time
  quantiles <- quantile(residence_df$passage_time, probs = c(0.25, 0.5, 0.75, 0.8, 0.95))
  list_quantiles[[i]] <- quantiles
  
  
}

# quantiles isn't really going to work for upstream LGR and WEL, need a mixture model for those
library(mixtools)

# Upstream of Lower Granite
LGR_residence_df <- residence_time(residence_state = "mainstem, upstream of LGR")
LGR_mix <- normalmixEM(LGR_residence_df$passage_time, k = 2)

# most fish are longer residence time (81% vs. shorter residence time)
LGR_mix$lambda
# shorter residence time is 24 days, longer is 168 days
LGR_mix$mu

# show the fit
plot(LGR_mix, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8,
     main2="Residence time, upstream of LGR", xlab2="Days")


# Upstream of Wells
WEL_residence_df <- residence_time(residence_state = "mainstem, upstream of WEL")
WEL_mix <- normalmixEM(WEL_residence_df$passage_time, k = 2)

# 51%/49% for short vs. long residence time at Wells, much different than at LGR
WEL_mix$lambda
# shorter residence time is 21 days, longer is 187 days
WEL_mix$mu

# show the fit
plot(WEL_mix, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8,
     main2="Residence time, upstream of WEL", xlab2="Days")


# make an actual table
quantiles_table <- data.frame(reach = main_mainstem_states,
                              t25 = c(list_quantiles[[1]][1],
                                      list_quantiles[[2]][1],
                                      list_quantiles[[3]][1],
                                      list_quantiles[[4]][1],
                                      list_quantiles[[5]][1],
                                      list_quantiles[[6]][1],
                                      list_quantiles[[7]][1],
                                      list_quantiles[[8]][1]),
                              t50 = c(list_quantiles[[1]][2],
                                      list_quantiles[[2]][2],
                                      list_quantiles[[3]][2],
                                      list_quantiles[[4]][2],
                                      list_quantiles[[5]][2],
                                      list_quantiles[[6]][2],
                                      list_quantiles[[7]][2],
                                      list_quantiles[[8]][2]),
                              t75 = c(list_quantiles[[1]][3],
                                      list_quantiles[[2]][3],
                                      list_quantiles[[3]][3],
                                      list_quantiles[[4]][3],
                                      list_quantiles[[5]][3],
                                      list_quantiles[[6]][3],
                                      list_quantiles[[7]][3],
                                      list_quantiles[[8]][3]),
                              t80 = c(list_quantiles[[1]][4],
                                      list_quantiles[[2]][4],
                                      list_quantiles[[3]][4],
                                      list_quantiles[[4]][4],
                                      list_quantiles[[5]][4],
                                      list_quantiles[[6]][4],
                                      list_quantiles[[7]][4],
                                      list_quantiles[[8]][4]),
                              t95 = c(list_quantiles[[1]][5],
                                      list_quantiles[[2]][5],
                                      list_quantiles[[3]][5],
                                      list_quantiles[[4]][5],
                                      list_quantiles[[5]][5],
                                      list_quantiles[[6]][5],
                                      list_quantiles[[7]][5],
                                      list_quantiles[[8]][5]))

## Use median residence time in each state to calculate a temperature window for each dam

quantiles_table %>% 
  mutate(median = round(t50)) -> quantiles_table

# match them using outflow

# use a function
# optional argument if you want to start window at a different time
# also update to make it flexible at the end of the time series to shorten window at the end of the time series (but this shouldn't really matter, we're not really seeing any movements in December 2022); but it'll just take the temp at that day
window_temp <- function(temp_data, start_window_days = 0, end_window_days){
  colnames(temp_data) <- c("date", "temp")
  
  temp_data %>% 
    mutate(window_temp = NA) -> temp_data
  
  # loop to get all windows
  for (i in 1:nrow(temp_data)){
    if (i > nrow(temp_data)-start_window_days) {
      temp_data$window_temp[i] <- temp_data$temp[i]
      
    } else {
      temp_data$window_temp[i] <- mean(subset(temp_data, date >= temp_data$date[i] + days(x = start_window_days) & date <= temp_data$date[i] + days(x = end_window_days))$temp)    
      
    }
    
    
  }
  
  return(temp_data)
}

# Bonneville Dam, for mainstem, BON to MCN reach
BON_window_temp <- window_temp(temp_data = BON_temp_ts, end_window_days = subset(quantiles_table, reach == "mainstem, BON to MCN")$median)

# McNary Dam, for mainstem, MCN to ICH or PRA reach
MCN_window_temp <- window_temp(temp_data = MCN_temp_ts, end_window_days = subset(quantiles_table, reach == "mainstem, MCN to ICH or PRA")$median)

# PRA Dam, for PRA to RIS reach
PRA_window_temp <- window_temp(temp_data = PRA_temp_ts, end_window_days = subset(quantiles_table, reach == "mainstem, PRA to RIS")$median)

# RIS Dam, for RIS to RRE reach
RIS_window_temp <- window_temp(temp_data = RIS_temp_ts, end_window_days = subset(quantiles_table, reach == "mainstem, RIS to RRE")$median)

# RRE Dam, for RRE to WEL reach
RRE_window_temp <- window_temp(temp_data = RRE_temp_ts, end_window_days = subset(quantiles_table, reach == "mainstem, RRE to WEL")$median)

# Ice Harbor Dam, for ICH to LGR reach
ICH_window_temp <- window_temp(temp_data = ICH_temp_ts, end_window_days = subset(quantiles_table, reach == "mainstem, ICH to LGR")$median)

# Wells Dam, for upstream of WEL reach - quick fish
WEL_quick_index <- which.min(WEL_mix$mu)
WEL_slow_index <- which.max(WEL_mix$mu)


WEL_quick_window_temp <- window_temp(temp_data = WEL_temp_ts, 
                                     start_window_days = round(WEL_mix$mu[WEL_quick_index] - WEL_mix$sigma[WEL_quick_index]), 
                                     end_window_days = round(WEL_mix$mu[WEL_quick_index] + WEL_mix$sigma[WEL_quick_index]))

# Wells Dam, for upstream of WEL reach - slow fish
WEL_slow_window_temp <- window_temp(temp_data = WEL_temp_ts, 
                                    start_window_days = round(WEL_mix$mu[WEL_slow_index] - WEL_mix$sigma[WEL_slow_index]), 
                                    end_window_days = round(WEL_mix$mu[WEL_slow_index] + WEL_mix$sigma[WEL_slow_index]))

# Lower Granite Dam, for upstream of LGR reach - quick fish
LGR_quick_index <- which.min(LGR_mix$mu)
LGR_slow_index <- which.max(LGR_mix$mu)

LGR_quick_window_temp <- window_temp(temp_data = LGR_temp_ts, 
                                     start_window_days = round(LGR_mix$mu[LGR_quick_index] - LGR_mix$sigma[LGR_quick_index]), 
                                     end_window_days = round(LGR_mix$mu[LGR_quick_index] + LGR_mix$sigma[LGR_quick_index]))

# Lower Granite Dam, for upstream of LGR reach - slow fish
LGR_slow_window_temp <- window_temp(temp_data = LGR_temp_ts, 
                                    start_window_days = round(LGR_mix$mu[LGR_slow_index] - LGR_mix$sigma[LGR_slow_index]), 
                                    end_window_days = round(LGR_mix$mu[LGR_slow_index] + LGR_mix$sigma[LGR_slow_index]))

## Reformat and export for Stan model 

# Trim to only the run years we're modeling: 05/06-23/24
dplyr::rename(dplyr::select(BON_window_temp, date, window_temp), BON = window_temp) %>% 
  left_join(dplyr::rename(dplyr::select(MCN_window_temp, date, window_temp), MCN = window_temp), by = "date") %>% 
  left_join(dplyr::rename(dplyr::select(PRA_window_temp, date, window_temp), PRA = window_temp), by = "date") %>% 
  left_join(dplyr::rename(dplyr::select(RIS_window_temp, date, window_temp), RIS = window_temp), by = "date") %>% 
  left_join(dplyr::rename(dplyr::select(RRE_window_temp, date, window_temp), RRE = window_temp), by = "date") %>% 
  left_join(dplyr::rename(dplyr::select(WEL_quick_window_temp, date, window_temp), WEL_quick = window_temp), by = "date") %>% 
  left_join(dplyr::rename(dplyr::select(ICH_window_temp, date, window_temp), ICH = window_temp), by = "date") %>% 
  left_join(dplyr::rename(dplyr::select(LGR_quick_window_temp, date, window_temp), LGR_quick = window_temp), by = "date") %>% 
  left_join(dplyr::rename(dplyr::select(WEL_slow_window_temp, date, window_temp), WEL_slow = window_temp), by = "date") %>%
  left_join(dplyr::rename(dplyr::select(LGR_slow_window_temp, date, window_temp), LGR_slow = window_temp), by = "date") %>%
  filter(date <= ymd("2024-05-31")) %>% 
  mutate(index = row_number()) %>% 
  relocate(index) %>% 
  dplyr::select(-date) %>% 
  column_to_rownames("index") -> window_temps

# before z-scoring, export the mean and SD of each set of window temps
window_temps %>% 
  mutate(BON_mean = mean(BON)) %>% 
  mutate(BON_sd = sd(BON)) %>% 
  mutate(MCN_mean = mean(MCN)) %>% 
  mutate(MCN_sd = sd(MCN)) %>% 
  mutate(PRA_mean = mean(PRA)) %>% 
  mutate(PRA_sd = sd(PRA)) %>% 
  mutate(RIS_mean = mean(RIS)) %>% 
  mutate(RIS_sd = sd(RIS)) %>% 
  mutate(RRE_mean = mean(RRE)) %>% 
  mutate(RRE_sd = sd(RRE)) %>% 
  mutate(WEL_quick_mean = mean(WEL_quick)) %>% 
  mutate(WEL_quick_sd = sd(WEL_quick)) %>% 
  mutate(WEL_slow_mean = mean(WEL_slow)) %>% 
  mutate(WEL_slow_sd = sd(WEL_slow)) %>% 
  mutate(ICH_mean = mean(ICH)) %>% 
  mutate(ICH_sd = sd(ICH)) %>% 
  mutate(LGR_quick_mean = mean(LGR_quick)) %>% 
  mutate(LGR_quick_sd = sd(LGR_quick)) %>% 
  mutate(LGR_slow_mean = mean(LGR_slow)) %>% 
  mutate(LGR_slow_sd = sd(LGR_slow)) %>% 
  dplyr::select(-c(BON, MCN, PRA, RIS, RRE, WEL_quick,
                   WEL_slow, ICH, LGR_quick, LGR_slow)) %>% 
  filter(!duplicated(BON_mean)) -> window_temps_summary


# now, z-score every column
window_temps %>% 
  mutate(BON = (BON - mean(BON))/sd(BON)) %>% 
  mutate(MCN = (MCN - mean(MCN))/sd(MCN)) %>% 
  mutate(PRA = (PRA - mean(PRA))/sd(PRA)) %>% 
  mutate(RIS = (RIS - mean(RIS))/sd(RIS)) %>% 
  mutate(RRE = (RRE - mean(RRE))/sd(RRE)) %>% 
  mutate(WEL_quick = (WEL_quick - mean(WEL_quick))/sd(WEL_quick)) %>% 
  mutate(WEL_slow = (WEL_slow - mean(WEL_slow))/sd(WEL_slow)) %>% 
  mutate(ICH = (ICH - mean(ICH))/sd(ICH)) %>% 
  mutate(LGR_quick = (LGR_quick - mean(LGR_quick))/sd(LGR_quick)) %>% 
  mutate(LGR_slow = (LGR_slow - mean(LGR_slow))/sd(LGR_slow)) -> window_temps

# add an empty column for the first state (that doesn't get a temperature effect)
# make sure that the order is the same is the order of states
window_temps %>% 
  mutate(MOUTH = 0) %>% 
  relocate(MOUTH) -> window_temps

# add rows of -999s for all of the dates that we aren't modeling to keep the dimensions the same
# for all the other data
extra_days <- data.frame(MOUTH = rep(-999, 7154-nrow(window_temps)),
                         BON = rep(-999, 7154-nrow(window_temps)),
                         MCN = rep(-999, 7154-nrow(window_temps)),
                         PRA = rep(-999, 7154-nrow(window_temps)),
                         RIS = rep(-999, 7154-nrow(window_temps)),
                         RRE = rep(-999, 7154-nrow(window_temps)),
                         WEL_quick = rep(-999, 7154-nrow(window_temps)),
                         ICH = rep(-999, 7154-nrow(window_temps)),
                         LGR_quick = rep(-999, 7154-nrow(window_temps)),
                         WEL_slow = rep(-999, 7154-nrow(window_temps)),
                         LGR_slow = rep(-999, 7154-nrow(window_temps)))

window_temps %>% 
  bind_rows(., extra_days) -> window_temps


## Export the files
write.csv(window_temps, here::here("Data", "covariate_data", "model_inputs", "window_temps_for_stan.csv"))
write.csv(window_temps_summary, here::here("Data", "covariate_data", "model_inputs", "window_temps_summary.csv"))

## Export another copy in the folders where the Stan models are being run
write.csv(window_temps, here::here("Stan", "middle_columbia_wild", "window_temps_for_stan.csv"))
write.csv(window_temps, here::here("Stan", "middle_columbia_hatchery", "window_temps_for_stan.csv"))
write.csv(window_temps, here::here("Stan", "upper_columbia_wild", "window_temps_for_stan.csv"))
write.csv(window_temps, here::here("Stan", "upper_columbia_hatchery", "window_temps_for_stan.csv"))
write.csv(window_temps, here::here("Stan", "snake_river_wild", "window_temps_for_stan.csv"))
write.csv(window_temps, here::here("Stan", "snake_river_hatchery", "window_temps_for_stan.csv"))

## Visualize temperatuer data that is going into stan model

temp_ts <- window_temps[!(window_temps$MOUTH == -999),]
dates <- seq(ymd("2005-06-01"), ymd("2024-05-31"), by = "days")
years <- year(dates)
temp_ts$date <- dates
temp_ts$year <- years

run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", 
              "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", 
              "20/21","21/22", "22/23", "23/24")
run_year_start <- seq(ymd("2004-06-01"), ymd("2023-06-01"), by = "years")
run_year_end <- seq(ymd("2005-05-31"), ymd("2024-05-31"), by = "years")
run_year_numeric = seq(4, 23, 1)

run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)


temp_ts %>% 
  rowwise() %>%
  dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date & run_year_end >= date)$run_year) %>% 
  group_by(run_year) %>% 
  summarise_all(median) %>% 
  arrange(BON) -> annual_temp_medians

temp_ts %>% 
  rowwise() %>%
  dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date & run_year_end >= date)$run_year) %>% 
  group_by(run_year) %>% 
  summarise_all(mean) %>% 
  arrange(BON) -> annual_temp_means


summary(annual_temp_medians$BON)
summary(annual_temp_means$BON)
summary(temp_ts$BON)

# so yeah, the year medians are slightly negative; the full dataset is centered, and the year means are essentially centered. That makes more sense

median_temp_by_run_year_BON <- ggplot(annual_temp_medians, aes(x = year, y = BON, label = run_year)) +
  geom_point() +
  ylab("Z-scored basin temperature") +
  geom_text(aes(y = BON + 0.01))

ggsave(here::here("figures", "data_plots", "median_temp_by_run_year_BON.png"), median_temp_by_run_year_BON, height = 8, width = 8)





#### define spill volume functions ####
# Function 1: reformat covariate data

# Spill is in kcfs (1000 cu ft/s), but the scale of this is still too high. To get it approx
# on the same scale as temp (which is z-scored), divide by 100. This means that the units
# now are in 100s of kcfs
spill_reformat <- function(input_data, variable_name){
  
  # need to deal with leap years
  leap_years <- seq(1960, 2024, by = 4)
  
  input_data %>% 
    pivot_longer(., cols = colnames(.)[2:ncol(.)]) %>% 
    dplyr::rename(year = name,
                  !!quo_name(variable_name) := value) %>% 
    mutate(year = gsub("X","",year)) %>% 
    mutate(spill = spill/100) %>% 
    mutate(year = as.numeric(year), day = as.numeric(day)) %>% 
    # drop day 366 in non-leap years
    dplyr::filter(!(day == 366 & !(year %in% leap_years))) %>%
    # Add a new column for date (combine year and day)
    mutate(date = format(as.Date(day, origin = paste0(year-1, "-12-31")))) %>% 
    arrange(year, day) %>% 
    # don't include 2025 data
    filter(year <= 2024) -> output_data
  
  return(output_data)
}

# Function 2: Calculate spill across a window of time
window_spill <- function(spill_data, start_window_days = 0, end_window_days){
  spill_data %>% 
    mutate(date = ymd(date)) %>% 
    filter(date >= ymd("2005-06-01") & date <= ymd("2024-12-31")) %>% 
    dplyr::select(date, spill) %>% 
    relocate(date) -> spill_data
  
  spill_data %>% 
    mutate(window_spill = NA) -> spill_data
  
  # loop to get all windows
  for (i in 1:nrow(spill_data)){
    if (i > nrow(spill_data)-start_window_days) {
      spill_data$window_spill[i] <- spill_data$spill[i]
      
    } else {
      spill_data$window_spill[i] <- mean(subset(spill_data, date >= spill_data$date[i] + days(x = start_window_days) & date <= spill_data$date[i] + days(x = end_window_days))$spill)    
      
    }
    
    
  }
  
  return(spill_data)
}

#### Load and reformat the spill data ####

# BON
BON_spill <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "BON", "BON_spill_2025-02-03.csv"))[1:366,]
BON_spill_long <- spill_reformat(input_data = BON_spill, variable_name = "spill")

# ICH
ICH_spill <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "ICH", "ICH_spill_2025-02-03.csv"))[1:366,]
ICH_spill_long <- spill_reformat(input_data = ICH_spill, variable_name = "spill")

# LGR
LGR_spill <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "LGR", "LGR_spill_2025-02-03.csv"))[1:366,]
LGR_spill_long <- spill_reformat(input_data = LGR_spill, variable_name = "spill")

# MCN
MCN_spill <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "MCN", "MCN_spill_2025-02-03.csv"))[1:366,]
MCN_spill_long <- spill_reformat(input_data = MCN_spill, variable_name = "spill")

# PRA
PRA_spill <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "PRA", "PRA_spill_2025-02-03.csv"))[1:366,]
PRA_spill_long <- spill_reformat(input_data = PRA_spill, variable_name = "spill")

# RIS
RIS_spill <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "RIS", "RIS_spill_2025-02-03.csv"))[1:366,]
RIS_spill_long <- spill_reformat(input_data = RIS_spill, variable_name = "spill")

# RRE
RRE_spill <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "RRE", "RRE_spill_2025-02-03.csv"))[1:366,]
RRE_spill_long <- spill_reformat(input_data = RRE_spill, variable_name = "spill")

# WEL
WEL_spill <- read.csv(here::here("Data", "covariate_data", "CBR_DART", "WEL", "WEL_spill_2025-02-03.csv"))[1:366,]
WEL_spill_long <- spill_reformat(input_data = WEL_spill, variable_name = "spill")

# so, some gaps in WEL spill data, on Jan 1 and 2 2020, and in December 2021.
# I think it's safe to assume that these are zero, especially looking at typical spill during these months
# let's look at that data
ggplot(filter(WEL_spill_long, year %in% c(2020, 2021)), aes(x = ymd(date), y = spill)) + geom_point()
# yeah, most likely zeros. will interpolate as such

WEL_spill_long %>% 
  mutate(spill = ifelse(is.na(spill) & year < 2023, 0, spill)) -> WEL_spill_long







#### Use the residence time windows calculated previously to calculate spill windows ####

# Bonneville Dam, for mainstem, BON to MCN reach
BON_window_spill <- window_spill(spill_data = BON_spill_long, end_window_days = subset(quantiles_table, reach == "mainstem, BON to MCN")$median)

# McNary Dam, for mainstem, MCN to ICH or PRA reach
MCN_window_spill <- window_spill(spill_data = MCN_spill_long, end_window_days = subset(quantiles_table, reach == "mainstem, MCN to ICH or PRA")$median)

# PRA Dam, for PRA to RIS reach
PRA_window_spill <- window_spill(spill_data = PRA_spill_long, end_window_days = subset(quantiles_table, reach == "mainstem, PRA to RIS")$median)

# RIS Dam, for RIS to RRE reach
RIS_window_spill <- window_spill(spill_data = RIS_spill_long, end_window_days = subset(quantiles_table, reach == "mainstem, RIS to RRE")$median)

# RRE Dam, for RRE to WEL reach
RRE_window_spill <- window_spill(spill_data = RRE_spill_long, end_window_days = subset(quantiles_table, reach == "mainstem, RRE to WEL")$median)

# Ice Harbor Dam, for ICH to LGR reach
ICH_window_spill <- window_spill(spill_data = ICH_spill_long, end_window_days = subset(quantiles_table, reach == "mainstem, ICH to LGR")$median)

# Wells Dam, for upstream of WEL reach - quick fish
WEL_quick_index <- which.min(WEL_mix$mu)
WEL_slow_index <- which.max(WEL_mix$mu)


WEL_quick_window_spill <- window_spill(spill_data = WEL_spill_long, 
                                       start_window_days = round(WEL_mix$mu[WEL_quick_index] - WEL_mix$sigma[WEL_quick_index]), 
                                       end_window_days = round(WEL_mix$mu[WEL_quick_index] + WEL_mix$sigma[WEL_quick_index]))

# Wells Dam, for upstream of WEL reach - slow fish
WEL_slow_window_spill <- window_spill(spill_data = WEL_spill_long, 
                                      start_window_days = round(WEL_mix$mu[WEL_slow_index] - WEL_mix$sigma[WEL_slow_index]), 
                                      end_window_days = round(WEL_mix$mu[WEL_slow_index] + WEL_mix$sigma[WEL_slow_index]))

# Lower Granite Dam, for upstream of LGR reach - quick fish
LGR_quick_index <- which.min(LGR_mix$mu)
LGR_slow_index <- which.max(LGR_mix$mu)

LGR_quick_window_spill <- window_spill(spill_data = LGR_spill_long, 
                                       start_window_days = round(LGR_mix$mu[LGR_quick_index] - LGR_mix$sigma[LGR_quick_index]), 
                                       end_window_days = round(LGR_mix$mu[LGR_quick_index] + LGR_mix$sigma[LGR_quick_index]))

# Lower Granite Dam, for upstream of LGR reach - slow fish
LGR_slow_window_spill <- window_spill(spill_data = LGR_spill_long, 
                                      start_window_days = round(LGR_mix$mu[LGR_slow_index] - LGR_mix$sigma[LGR_slow_index]), 
                                      end_window_days = round(LGR_mix$mu[LGR_slow_index] + LGR_mix$sigma[LGR_slow_index]))

## Reformat and export for Stan model 

# Trim to only the run years we're modeling: 05/06-23/24
dplyr::rename(dplyr::select(BON_window_spill, date, window_spill), BON = window_spill) %>% 
  left_join(dplyr::rename(dplyr::select(MCN_window_spill, date, window_spill), MCN = window_spill), by = "date") %>% 
  left_join(dplyr::rename(dplyr::select(PRA_window_spill, date, window_spill), PRA = window_spill), by = "date") %>% 
  left_join(dplyr::rename(dplyr::select(RIS_window_spill, date, window_spill), RIS = window_spill), by = "date") %>% 
  left_join(dplyr::rename(dplyr::select(RRE_window_spill, date, window_spill), RRE = window_spill), by = "date") %>% 
  left_join(dplyr::rename(dplyr::select(WEL_quick_window_spill, date, window_spill), WEL_quick = window_spill), by = "date") %>% 
  left_join(dplyr::rename(dplyr::select(ICH_window_spill, date, window_spill), ICH = window_spill), by = "date") %>% 
  left_join(dplyr::rename(dplyr::select(LGR_quick_window_spill, date, window_spill), LGR_quick = window_spill), by = "date") %>% 
  left_join(dplyr::rename(dplyr::select(WEL_slow_window_spill, date, window_spill), WEL_slow = window_spill), by = "date") %>%
  left_join(dplyr::rename(dplyr::select(LGR_slow_window_spill, date, window_spill), LGR_slow = window_spill), by = "date") %>%
  filter(date <= ymd("2024-05-31")) %>% 
  mutate(index = row_number()) %>% 
  relocate(index) %>% 
  dplyr::select(-date) %>% 
  column_to_rownames("index") -> window_spill

# add an empty column for the first state (that doesn't get a spill effect)
# make sure that the order is the same is the order of states
window_spill %>% 
  mutate(MOUTH = 0) %>% 
  relocate(MOUTH) -> window_spill

# add rows of -999s for all of the dates that we aren't modeling to keep the dimensions the same
# for all the other data
extra_days <- data.frame(MOUTH = rep(-999, 7154-nrow(window_spill)),
                         BON = rep(-999, 7154-nrow(window_spill)),
                         MCN = rep(-999, 7154-nrow(window_spill)),
                         PRA = rep(-999, 7154-nrow(window_spill)),
                         RIS = rep(-999, 7154-nrow(window_spill)),
                         RRE = rep(-999, 7154-nrow(window_spill)),
                         WEL_quick = rep(-999, 7154-nrow(window_spill)),
                         ICH = rep(-999, 7154-nrow(window_spill)),
                         LGR_quick = rep(-999, 7154-nrow(window_spill)),
                         WEL_slow = rep(-999, 7154-nrow(window_spill)),
                         LGR_slow = rep(-999, 7154-nrow(window_spill)))

window_spill %>% 
  bind_rows(., extra_days) -> window_spill

## Export the files
write.csv(window_spill, here::here("Data", "covariate_data", "model_inputs", "window_spill_for_stan.csv"))

## Export another copy in the folders where the Stan models are being run
write.csv(window_spill, here::here("Stan", "middle_columbia_wild", "window_spill_for_stan.csv"))
write.csv(window_spill, here::here("Stan", "middle_columbia_hatchery", "window_spill_for_stan.csv"))
write.csv(window_spill, here::here("Stan", "upper_columbia_wild", "window_spill_for_stan.csv"))
write.csv(window_spill, here::here("Stan", "upper_columbia_hatchery", "window_spill_for_stan.csv"))
write.csv(window_spill, here::here("Stan", "snake_river_wild", "window_spill_for_stan.csv"))
write.csv(window_spill, here::here("Stan", "snake_river_hatchery", "window_spill_for_stan.csv"))

#### define spill days functions ####
april_spill_days <- function(spill_data){
  spill_data %>% 
    # keep only 2005-2024
    subset(year >= 2005 & year <= 2024) %>% 
    mutate(date = ymd(date)) %>% 
    mutate(month = month(date)) %>% 
    mutate(any_spill = ifelse(spill > 0, "yes", "no")) %>% 
    group_by(year, month) %>% 
    count(any_spill == "yes") %>% 
    dplyr::rename(positive_spill = `any_spill == "yes"`) %>% 
    subset(positive_spill == TRUE & month == 4) %>% 
    dplyr::rename(april_spill_days = n) %>% 
    mutate(april_spill_days = april_spill_days/100) %>% 
    ungroup() %>% 
    complete(year = seq(2005,2024,1), fill = list(april_spill_days = 0)) %>% 
    dplyr::select(year, april_spill_days) -> april_spill_days
  
  return(april_spill_days)
}

march_spill_days <- function(spill_data){
  spill_data %>% 
    # keep only 2005-2024
    subset(year >= 2005 & year <= 2024) %>% 
    mutate(date = ymd(date)) %>% 
    mutate(month = month(date)) %>% 
    mutate(any_spill = ifelse(spill > 0, "yes", "no")) %>% 
    group_by(year, month) %>% 
    count(any_spill == "yes") %>% 
    dplyr::rename(positive_spill = `any_spill == "yes"`) %>% 
    subset(positive_spill == TRUE & month == 3) %>% 
    dplyr::rename(march_spill_days = n) %>% 
    mutate(march_spill_days = march_spill_days/100) %>% 
    ungroup() %>% 
    complete(year = seq(2005,2024,1), fill = list(march_spill_days = 0)) %>% 
    dplyr::select(year, march_spill_days) -> march_spill_days
  
  return(march_spill_days)
}

february_spill_days <- function(spill_data){
  spill_data %>% 
    # keep only 2005-2024
    subset(year >= 2005 & year <= 2024) %>% 
    mutate(date = ymd(date)) %>% 
    mutate(month = month(date)) %>% 
    mutate(any_spill = ifelse(spill > 0, "yes", "no")) %>% 
    group_by(year, month) %>% 
    count(any_spill == "yes") %>% 
    dplyr::rename(positive_spill = `any_spill == "yes"`) %>% 
    subset(positive_spill == TRUE & month == 2) %>% 
    dplyr::rename(february_spill_days = n) %>% 
    mutate(february_spill_days = february_spill_days/100) %>% 
    ungroup() %>% 
    complete(year = seq(2005,2024,1), fill = list(february_spill_days = 0)) %>% 
    dplyr::select(year, february_spill_days) -> february_spill_days
  
  return(february_spill_days)
}

january_spill_days <- function(spill_data){
  spill_data %>% 
    # keep only 2005-2024
    subset(year >= 2005 & year <= 2024) %>% 
    mutate(date = ymd(date)) %>% 
    mutate(month = month(date)) %>% 
    mutate(any_spill = ifelse(spill > 0, "yes", "no")) %>% 
    group_by(year, month) %>% 
    count(any_spill == "yes") %>% 
    dplyr::rename(positive_spill = `any_spill == "yes"`) %>% 
    subset(positive_spill == TRUE & month == 1) %>% 
    dplyr::rename(january_spill_days = n) %>% 
    mutate(january_spill_days = january_spill_days/100) %>% 
    ungroup() %>% 
    complete(year = seq(2005,2024,1), fill = list(january_spill_days = 0)) %>% 
    dplyr::select(year, january_spill_days) -> january_spill_days
  
  return(january_spill_days)
}

BON_april_spill <- april_spill_days(spill_data = BON_spill_long)
MCN_april_spill <- april_spill_days(spill_data = MCN_spill_long)
PRA_april_spill <- april_spill_days(spill_data = PRA_spill_long)
RIS_april_spill <- april_spill_days(spill_data = RIS_spill_long)
RRE_april_spill <- april_spill_days(spill_data = RRE_spill_long)
WEL_april_spill <- april_spill_days(spill_data = WEL_spill_long)
ICH_april_spill <- april_spill_days(spill_data = ICH_spill_long)
LGR_april_spill <- april_spill_days(spill_data = LGR_spill_long)


BON_march_spill <- march_spill_days(spill_data = BON_spill_long)
MCN_march_spill <- march_spill_days(spill_data = MCN_spill_long)
PRA_march_spill <- march_spill_days(spill_data = PRA_spill_long)
RIS_march_spill <- march_spill_days(spill_data = RIS_spill_long)
RRE_march_spill <- march_spill_days(spill_data = RRE_spill_long)
WEL_march_spill <- march_spill_days(spill_data = WEL_spill_long)
ICH_march_spill <- march_spill_days(spill_data = ICH_spill_long)
LGR_march_spill <- march_spill_days(spill_data = LGR_spill_long)

BON_february_spill <- february_spill_days(spill_data = BON_spill_long)
MCN_february_spill <- february_spill_days(spill_data = MCN_spill_long)
PRA_february_spill <- february_spill_days(spill_data = PRA_spill_long)
RIS_february_spill <- february_spill_days(spill_data = RIS_spill_long)
RRE_february_spill <- february_spill_days(spill_data = RRE_spill_long)
WEL_february_spill <- february_spill_days(spill_data = WEL_spill_long)
ICH_february_spill <- february_spill_days(spill_data = ICH_spill_long)
LGR_february_spill <- february_spill_days(spill_data = LGR_spill_long)

BON_january_spill <- january_spill_days(spill_data = BON_spill_long)
MCN_january_spill <- january_spill_days(spill_data = MCN_spill_long)
PRA_january_spill <- january_spill_days(spill_data = PRA_spill_long)
RIS_january_spill <- january_spill_days(spill_data = RIS_spill_long)
RRE_january_spill <- january_spill_days(spill_data = RRE_spill_long)
WEL_january_spill <- january_spill_days(spill_data = WEL_spill_long)
ICH_january_spill <- january_spill_days(spill_data = ICH_spill_long)
LGR_january_spill <- january_spill_days(spill_data = LGR_spill_long)

#### Export spill as number of days of spill for post-overshoot fallback covariate ####

dplyr::rename(BON_january_spill, BON = january_spill_days) %>% 
  left_join(dplyr::rename(MCN_january_spill, MCN = january_spill_days), by = "year") %>% 
  left_join(dplyr::rename(PRA_january_spill, PRA = january_spill_days), by = "year") %>% 
  left_join(dplyr::rename(RIS_january_spill, RIS = january_spill_days), by = "year") %>% 
  left_join(dplyr::rename(RRE_january_spill, RRE = january_spill_days), by = "year") %>% 
  left_join(dplyr::rename(WEL_january_spill, WEL = january_spill_days), by = "year") %>% 
  left_join(dplyr::rename(ICH_january_spill, ICH = january_spill_days), by = "year") %>% 
  left_join(dplyr::rename(LGR_january_spill, LGR = january_spill_days), by = "year") %>% 
  dplyr::select(-year) %>% 
  mutate(index = row_number()) %>% 
  relocate(index) %>% 
  column_to_rownames("index") -> january_spill_df

# add an empty column for the first state (that doesn't get a temperature effect)
# make sure that the order is the same is the order of states
january_spill_df %>% 
  mutate(MOUTH = 0) %>% 
  relocate(MOUTH) -> january_spill_df

dplyr::rename(BON_february_spill, BON = february_spill_days) %>% 
  left_join(dplyr::rename(MCN_february_spill, MCN = february_spill_days), by = "year") %>% 
  left_join(dplyr::rename(PRA_february_spill, PRA = february_spill_days), by = "year") %>% 
  left_join(dplyr::rename(RIS_february_spill, RIS = february_spill_days), by = "year") %>% 
  left_join(dplyr::rename(RRE_february_spill, RRE = february_spill_days), by = "year") %>% 
  left_join(dplyr::rename(WEL_february_spill, WEL = february_spill_days), by = "year") %>% 
  left_join(dplyr::rename(ICH_february_spill, ICH = february_spill_days), by = "year") %>% 
  left_join(dplyr::rename(LGR_february_spill, LGR = february_spill_days), by = "year") %>% 
  dplyr::select(-year) %>% 
  mutate(index = row_number()) %>% 
  relocate(index) %>% 
  column_to_rownames("index") -> february_spill_df

# add an empty column for the first state (that doesn't get a temperature effect)
# make sure that the order is the same is the order of states
february_spill_df %>% 
  mutate(MOUTH = 0) %>% 
  relocate(MOUTH) -> february_spill_df

dplyr::rename(BON_march_spill, BON = march_spill_days) %>% 
  left_join(dplyr::rename(MCN_march_spill, MCN = march_spill_days), by = "year") %>% 
  left_join(dplyr::rename(PRA_march_spill, PRA = march_spill_days), by = "year") %>% 
  left_join(dplyr::rename(RIS_march_spill, RIS = march_spill_days), by = "year") %>% 
  left_join(dplyr::rename(RRE_march_spill, RRE = march_spill_days), by = "year") %>% 
  left_join(dplyr::rename(WEL_march_spill, WEL = march_spill_days), by = "year") %>% 
  left_join(dplyr::rename(ICH_march_spill, ICH = march_spill_days), by = "year") %>% 
  left_join(dplyr::rename(LGR_march_spill, LGR = march_spill_days), by = "year") %>% 
  dplyr::select(-year) %>% 
  mutate(index = row_number()) %>% 
  relocate(index) %>% 
  column_to_rownames("index") -> march_spill_df

# add an empty column for the first state (that doesn't get a temperature effect)
# make sure that the order is the same is the order of states
march_spill_df %>% 
  mutate(MOUTH = 0) %>% 
  relocate(MOUTH) -> march_spill_df

dplyr::rename(BON_april_spill, BON = april_spill_days) %>% 
  left_join(dplyr::rename(MCN_april_spill, MCN = april_spill_days), by = "year") %>% 
  left_join(dplyr::rename(PRA_april_spill, PRA = april_spill_days), by = "year") %>% 
  left_join(dplyr::rename(RIS_april_spill, RIS = april_spill_days), by = "year") %>% 
  left_join(dplyr::rename(RRE_april_spill, RRE = april_spill_days), by = "year") %>% 
  left_join(dplyr::rename(WEL_april_spill, WEL = april_spill_days), by = "year") %>% 
  left_join(dplyr::rename(ICH_april_spill, ICH = april_spill_days), by = "year") %>% 
  left_join(dplyr::rename(LGR_april_spill, LGR = april_spill_days), by = "year") %>% 
  dplyr::select(-year) %>% 
  mutate(index = row_number()) %>% 
  relocate(index) %>% 
  column_to_rownames("index") -> april_spill_df

# add an empty column for the first state (that doesn't get a temperature effect)
# make sure that the order is the same is the order of states
april_spill_df %>% 
  mutate(MOUTH = 0) %>% 
  relocate(MOUTH) -> april_spill_df

write.csv(january_spill_df, here::here("Data", "covariate_data", "model_inputs", "january_spill_df_for_stan.csv"))
write.csv(february_spill_df, here::here("Data", "covariate_data", "model_inputs", "february_spill_df_for_stan.csv"))
write.csv(march_spill_df, here::here("Data", "covariate_data", "model_inputs", "march_spill_df_for_stan.csv"))
write.csv(april_spill_df, here::here("Data", "covariate_data", "model_inputs", "april_spill_df_for_stan.csv"))



