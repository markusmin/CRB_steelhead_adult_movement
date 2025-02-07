#!/bin/bash
## clean detection histories for CTH files 5_8

## Job name
#SBATCH --job-name=det_hist_2

## Partition and Allocation
#SBATCH -p ckpt
#SBATCH -A stf

## Resources
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --ntasks=3
#SBATCH --mem=100G

## Specify the working directory for this job
#SBATCH --chdir=/gscratch/stf/mmin/CRB_steelhead_adult_movement

## We don't need any modules beyond the default loaded ones

##  Scripts to be executed here
apptainer exec --bind ../../../:/gscratch/stf/mmin ../cmdstanr_latest.sif Rscript R/hyak/02_clean_detection_histories/02_clean_detection_histories_CTH_5_8.R
