# CRB_steelhead_adult_movement
Workflow for using PIT tag detections to fit a multistate movement model to adult Steelhead returning to the Columbia River Basin.

*last updated: 2025-01-29*

Note: The data files are very large and therefore are not uploaded to GitHub. If you would like access to the data files, please contact me directly.

To reproduce this analysis, you must run the following scripts **sequentially**:

The directories are set up as follows:
-- Data ------- all input data files for project
  |-- covariate_data ------- temperature, flow, and spill data queried from CBR DART
  |-- PTAGIS_queries ------- tag histories for individual fish queried from PTAGIS
-- R/ ------- all scripts necessary to reproduce analysis
  |-- functions ------- individual function definitions, which are sourced by other scripts in the R/ folder
  |-- hyak ------- jobs and scripts for running the workflow on UW' Hyak HPC
-- intermediate_outputs/ ------- outputs from scripts in the R/ folder that are inputs to later scripts



## Part 1: Prepare the fish movement data for inclusion in the model

#### Step 1: Query from PTAGIS

The code to setup/explain the queries is can be found in `R/01_data_querying.Rmd`.

All of the PTAGIS files that were queried can be found in `Data/PTAGIS_queries`

#### Step 2: Convert raw detections into detections at specific sites

This step processes the raw complete tag histories from the first step into a detection
history for each fish. It distills the detections at individual antennas down to
detections at each site and uses the sequence of detections at different antennas
to determine the directionality of movement.

This step can be run either entirely locally or partially using a high-performance computing (HPC) cluster.

The primary script to run this is `R/02_clean_detection_histories.R`. The script notes steps
that can be run separately on a HPC cluster. Code to run these steps on UW's Hyak HPC
cluster can be found in `R/hyak/02_clean_detection_histories`. The three jobs that need to be run are the following:

1. job_ckpt_02_clean_detection_histories_CTH_1_4.slurm
2. job_ckpt_02_clean_detection_histories_CTH_5_8.slurm
3. job_ckpt_02_clean_detection_histories_CTH_9_13.slurm

The primary output from this step is the following file: `intermediate_outputs/complete_det_hist_postprocessed.csv`.
If you wish to run step 3 via HPC, you must move this file to the cluster before proceeding.

#### Step 3: Convert detections at specific sites into state occupancy/transitions

This script takes the `intermediate_outputs/complete_det_hist_postprocessed.csv` file that was exported 
during step 2 and converts this series of detections at different sites into a 
history of transitions between different states in our model. This step can be run 
either entirely locally or partially using a high-performance computing (HPC) cluster. 
It takes as input a file from PTAGIS that contains information on all interrogation 
sites. This file was created using the following steps:

1. From the [PTAGIS Dashboard](https://dashboard.ptagis.org): Advanced Reporting Home Page > Standard Reports > MRR and Interrogation Site Info > All Interrogation Sites Table
2. Add (drag and drop from REPORT OBJECTS to the table view) the following columns: Interrogation Site Basin
3. Export as CSV

This file was stored as  Data/All Interrogation Sites Table 2025-01-30.csv

The primary script to run this is `R/03_create_state_histories.R`. The script notes steps
that can be run separately on a HPC cluster. Code to run these steps on UW's Hyak HPC cluster 
can be found in `R/hyak/03_create_state_histories`. The four jobs that need to be run are the following:

1. job_ckpt_03_create_state_histories_part1.slurm
2. job_ckpt_03_create_state_histories_part2.slurm
3. job_ckpt_03_create_state_histories_part3.slurm
4. job_ckpt_03_create_state_histories_part4.slurm

The script `R/hyak/03_create_state_histories` exports the complete state history 
for fish from each DPS in the appropriate folders within the `Stan/` folder. These
are the input files for the R scripts that also call and run the Stan model
for each DPS.


## Part 2: Prepare the covariate data for inclusion in the model

#### Step 1: Download the data from CBR DART

The four data types were downloaded:
1. Temperature C, Daily Average -> "temp"
2. Temperature C, Tailrace, Daily Average -> "tailrace_temp"
3. Outflow KCFS, Daily Average -> "flow"
4. Spill KCFS, Daily Average -> "spill"

These files are stored in `Data/covariate_data/CBR_DART`. All files have the suffix 
"2025-02-03" to indicate the date that they were downloaded.

#### Step 2: Process the data for inclusion in the model

The script `04_process_covariate_data.R` processes both the temperature and spill data
for input into the model, using the following steps.

To process the temperature data:

1. Filters out outlier temperatures (a) manually (based on a visual inspection of 
temperature data), (b) based on interannual averages (data points >4 degrees from 
the average value for that day of year are removed), and (c) based on a moving average 
(data points >2 degrees from the 7-day moving average are removed).
2. Fits a MARSS model to estimate temperature throughout the basin, which addresses
the gaps in the data. This model uses as inputs the forebay and tailrace temperatures 
from the previous step and uses them to fit a MARSS model to estimate temperature.
3. Estimates temperature within windows of time experienced by individual fish. This step
uses as inputs the MARSS-estimated temperatures from the previous step, calculates 
the median residence time for in each state, and then uses that median residence 
time to calculate the mean temperature using a window for each date (defined as the 
window of time from the date plus the median residence time in that state). 
These window temperatures are exported in `Data/covariate_data/model_inputs/window_temps_for_stan.csv`,
as well as in the individual modeling folders in `Stan/`.

To process the spill data:

The script takes the daily average spill at each dam and processes it in two ways: 
by calculating a window of the volume of spill around a date (the same approach as 
is used for temperature) and by calculating the total number of days of spill per 
year in the winter months (January, February, and March). The spill window data 
is exported in `Data/covariate_data/model_inputs/window_spill_for_stan.csv`
and the spill days data is exported in `Data/covariate_data/model_inputs/january_spill_df_for_stan.csv`, 
`Data/covariate_data/model_inputs/february_spill_df_for_stan.csv`, and 
`Data/covariate_data/model_inputs/march_spill_df_for_stan.csv`. Note that 
although the spill days data is exported separately for each of the winter months, 
in the primary version of the model it is combined into a single value for the 
total number of days with any amount of spill across all months. These spill
files are also exported in the individual modeling folders in `Stan/`.

#### Step 4: Estimate detection efficiency in tributaries

Estimating detection efficiency in tributaries was done in three steps, all
of which are documented in `R/05_estimate_tributary_detection_efficiency.Rmd`.
1. Query discharge data from the [USGS dashboard](https://dashboard.waterdata.usgs.gov/app/nwd/en/)
2. Reformat discharge data
3. Fit discharge model in Stan





#### Step 5: Combine movement and covariate data for input into Stan model

In the folder `R/06_prepare_stan_inputs`, there is a script that prepares the
data for each DPS (e.g., `R/06_prepare_stan_inputs/MCH_stan_data_prep.R`).









