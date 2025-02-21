#!/bin/bash
# submit the different FF_cov jobs

sbatch 07-07-02_jdr_ff_cov_job.slurm —mpi=pmix_v1
sbatch 07-07-03_uma_ff_cov_job.slurm —mpi=pmix_v1
sbatch 07-07-04_wawa_ff_cov_job.slurm —mpi=pmix_v1
sbatch 07-07-05_wen_ff_cov_job.slurm —mpi=pmix_v1
sbatch 07-07-06_tuc_ff_cov_job.slurm —mpi=pmix_v1
sbatch 07-07-07_ent_ff_cov_job.slurm —mpi=pmix_v1
sbatch 07-07-08_yak_ff_cov_job.slurm —mpi=pmix_v1
