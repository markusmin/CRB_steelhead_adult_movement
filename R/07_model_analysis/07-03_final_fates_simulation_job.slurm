#!/bin/bash

## Job name
#SBATCH --job-name=FF_sim

## Partition and Allocation
#SBATCH -p ckpt
#SBATCH -A stf

## Resources
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --ntasks=3
#SBATCH --mem=100G

## Specify the working directory for this job
#SBATCH --chdir=.


## We don't need any modules beyond the default loaded ones

##  Scripts to be executed here
apptainer exec --bind ../../:/gscratch/scrubbed/mmin ../../../../../gscratch/stf/mmin/cmdstanr_latest.sif Rscript 07-03_final_fates_simulation.R
