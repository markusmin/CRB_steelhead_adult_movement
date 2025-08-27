#!/bin/bash
## Upper Columbia, Wild chain 2


## Job name
#SBATCH --job-name=UCW2
#SBATCH --mail-type=ALL

## Partition and Allocation
#SBATCH -p compute
#SBATCH -A stf

## Resources
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --ntasks=28
#SBATCH --mem=120G

## Specify the working directory for this job
#SBATCH --chdir=.

## We don't need any modules beyond the default loaded ones

##  Scripts to be executed here
apptainer exec --bind ../../../../:/gscratch/stf/mmin ../../../../cmdstanr_latest.sif Rscript UCW_marchspill_chain2.R

