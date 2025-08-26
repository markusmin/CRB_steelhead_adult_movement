#!/bin/bash
# submit the four different jobs for four different chains

sbatch job_SRH_NT_chain1.slurm —mpi=pmix_v1
sbatch job_SRH_NT_chain2.slurm —mpi=pmix_v1
sbatch job_SRH_NT_chain3.slurm —mpi=pmix_v1
sbatch job_SRH_NT_chain4.slurm —mpi=pmix_v1