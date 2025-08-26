#!/bin/bash
# submit the four different jobs for four different chains

sbatch job_SRH_chain1.slurm —mpi=pmix_v1
sbatch job_SRH_chain2.slurm —mpi=pmix_v1
sbatch job_SRH_chain3.slurm —mpi=pmix_v1
sbatch job_SRH_chain4.slurm —mpi=pmix_v1