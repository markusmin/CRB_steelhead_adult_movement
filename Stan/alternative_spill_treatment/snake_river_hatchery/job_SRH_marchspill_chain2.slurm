#!/bin/bash
## Snake River, Hatchery chain 2


## Job name
#SBATCH --job-name=SRH2
#SBATCH --mail-type=ALL

## Partition and Allocation
#SBATCH -p compute
#SBATCH -A stf

## Resources
#SBATCH --nodes=1
#SBATCH --time=240:00:00
#SBATCH --ntasks=28
#SBATCH --mem=120G

## Specify the working directory for this job
#SBATCH --chdir=.

## We don't need any modules beyond the default loaded ones

##  Scripts to be executed here
apptainer exec --bind ../../../../:/gscratch/stf/mmin ../../../../cmdstanr_latest.sif Rscript SRH_marchspill_chain2.R

