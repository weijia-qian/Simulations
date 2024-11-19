#!/bin/bash
#SBATCH --array=1-27%20
#SBATCH --job-name=simulation_job
#SBATCH --partition=chang
#SBATCH --output=simulation.out
#SBATCH --error=simulation.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript simulation.R $JOBID


