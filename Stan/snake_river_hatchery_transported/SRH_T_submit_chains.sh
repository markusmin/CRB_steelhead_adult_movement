#!/bin/bash
# submit the four different jobs for four different chains

sbatch job_SRH_T_chain1.slurm —mpi=pmix_v1
sbatch job_SRH_T_chain2.slurm —mpi=pmix_v1
sbatch job_SRH_T_chain3.slurm —mpi=pmix_v1
sbatch job_SRH_T_chain4.slurm —mpi=pmix_v1