#!/bin/bash

#SBATCH --job-name=expr
#SBATCH --output=../log/expr/%a.out
#SBATCH --array=1-N_JOBS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00

#SBATCH --mail-user=yixiangluo@berkeley.edu
#SBATCH --mail-type=ALL

EXPERIMENT=expr

# run experiments

ml R

Rscript ../R/cluster/experiment.R $EXPERIMENT $SLURM_ARRAY_TASK_ID


