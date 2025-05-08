# CRB_steelhead_adult_movement
Workflow for using PIT tag detections to fit a multistate movement model to adult Steelhead returning to the Columbia River Basin.

*last updated: 2025-05-07*

Note: Some files are not uploaded because they are very large. These include the raw data files are very large and the modeling outputs.
For the raw data, please contact me directly or follow the instructions below to query PTAGIS and obtain the files. For
the modeling outputs, please contact me directly or follow the steps below to fit the models.

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

#### Step 6: Run the Stan model

This step was run on UW's Hyak high perfomance computing cluster.

The scripts can all be found in the `Stan/` folder:
- the *_model.stan file is the stan model
- the *_chain[1-4].R files each run one chain of the Stan model
- the job_*_chain[1-4].slurm are jobs that run the files above on hyak.
- the *_submit_chains.sh scripts are shell scripts that submit all of the above jobs.

#### Step 7: Analyze model results

This step is handled by multiple scripts stored in the `R/07_model_analysis` folder
that examine model diagnostics and create different results for the paper. Because 
the model output files are quite large, many of these scripts were run on Hyak.

The scripts within this folder are as follows:
- `07-00_install_dependencies.R`
  - This script installs all dependices needed to run the following scripts in this folder.
- `07-01_load_stan_models.R`
  - This script loads the Stan model outputs, combines the chains into one model
  object, and does some processing of the chains to give parameter estimates.
- `07-02_stan_model_diagnostics.R`
  - This script generates diagnostic figures for the model, including histograms
  of rhat and effective sample size.
- `07-03_final_fates_simulation.R`
  - This script takes the output from the stan model runs 
  generates estimates of final fate distributions
- `07-04_Fig3_final_fate_overshoot.R`
  - This script takes the final fate distributions estimated in `07-03_final_fates.R`
  and generates figure 3 in the manuscript, which shows the probability of homing 
  as a final fate vs. the probability of overshoot under median conditions.
- `07-05_Fig4_temperature_effects_plot.R`
  - This script generates figure 4 in the manuscript, which shows the model 
  estimated temperature response when fish reach the mainstem state that connects 
  to their natal tributary.
- `07-06_Fig5_spilldays_effects_plot.R`
  - This script generates figure 5 in the manuscript, which shows the model 
  estimated spill days response when fish have overshot their natal tributary.
- `/07-07_final_fates_scenarios/`
  - This is a folder that contains multiple jobs and scripts to generate estimates
  of final fates under different scenarios of temperature and spill, using the
  model-estimated parameter posteriors
  - `07-07_submit_FF_cov.sh` is a Bash script that runs all of the below jobs/scripts
  - `07-07-01-ff_cov_functions.R` contains many of the functions that are necessary
  for estimating final fates under different scenarios
  - A separate R script and associated job to run the script on Hyak exists for
  each population of interest:
    - `07-07-02_jdr_ff_cov.R`
    - `07-07-02_jdr_ff_cov_job.slurm`
    - `07-07-03_uma_ff_cov.R`
    - `07-07-03_uma_ff_cov_job.slurm`
    - `07-07-04_wawa_ff_cov.R`
    - `07-07-04_wawa_ff_cov_job.slurm`
    - `07-07-05_wen_ff_cov.R`
    - `07-07-05_wen_ff_cov_job.slurm`
    - `07-07-06_tuc_ff_cov.R`
    - `07-07-06_tuc_ff_cov_job.slurm`
    - `07-07-07_ent_ff_cov.R`
    - `07-07-07_ent_ff_cov_job.slurm`
    - `07-07-08_yak_ff_cov.R`
    - `07-07-08_yak_ff_cov_job.slurm`
- `07-08_Fig6_final_fates_scenarios.R`
  - This script uses the output generated by the scripts in folder
  `/07-07_final_fates_scenarios/` and generates figure 6 in the manuscript,
  which shows how natal homing varies under different temperature and spill scenarios.
- `07-09_S8_spillvolume_effects_plot.R`
  - This script generates generates figures for the supplement which shows the 
  estimated effect of the volume of spill at four dams in the Columbia River Basin
  on en-route fallback.
- `/alternative_spill_treatment/`
  - This folder contains scripts that generates all of the core manuscript figures
  when the effect of winter spill on post-overshoot fallback is modeled as only 
  the days of spill in March, rather than as all days in January, February, and
  March as is done in the base model.

  


#### Step 8: Generate additional figures and tables

This step generates other results for the paper that are not related to the Stan
model outputs.

- `08_generate_map.R`
  - This generates the map of the study region, which is Fig. 1A. Fig. 1B and 
  Fig. 2 are not generated by a script.












