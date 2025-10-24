#!/bin/bash
#SBATCH --array=1-5%5
#SBATCH --job-name=simulation_job
#SBATCH --partition=wrobel
#SBATCH --output=main_simulation_fttm_time.out
#SBATCH --error=main_simulation_fttm_time.err

module purge
module load R/4.4.0

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript main_simulation_fttm_time.R $JOBID


