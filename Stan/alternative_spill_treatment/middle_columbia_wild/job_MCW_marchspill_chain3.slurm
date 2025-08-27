#!/bin/bash
## Middle Columbia, Wild chain 3


## Job name
#SBATCH --job-name=MCW3
#SBATCH --mail-type=ALL

## Partition and Allocation
#SBATCH -p compute
#SBATCH -A stf

## Resources
#SBATCH --nodes=1
#SBATCH --time=84:00:00
#SBATCH --ntasks=28
#SBATCH --mem=120G

## Specify the working directory for this job
#SBATCH --chdir=.

## We don't need any modules beyond the default loaded ones

##  Scripts to be executed here
apptainer exec --bind ../../../../:/gscratch/stf/mmin ../../../../cmdstanr_latest.sif Rscript MCW_marchspill_chain3.R

