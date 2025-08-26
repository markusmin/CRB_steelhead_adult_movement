#!/bin/bash
# submit the four different jobs for four different chains

sbatch job_MCW_marchspill_chain1.slurm —mpi=pmix_v1
sbatch job_MCW_marchspill_chain2.slurm —mpi=pmix_v1
sbatch job_MCW_marchspill_chain3.slurm —mpi=pmix_v1
sbatch job_MCW_marchspill_chain4.slurm —mpi=pmix_v1