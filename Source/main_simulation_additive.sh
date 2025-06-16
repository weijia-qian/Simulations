#!/bin/bash
#SBATCH --array=1-18%18
#SBATCH --job-name=simulation_job
#SBATCH --partition=wrobel
#SBATCH --output=main_simulation_additive.out
#SBATCH --error=main_simulation_additive.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript main_simulation_additive.R $JOBID


