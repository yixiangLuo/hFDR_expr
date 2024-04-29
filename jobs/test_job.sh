#!/bin/bash

#SBATCH --job-name=test
#SBATCH --output=../log/test/%a.out
#SBATCH --array=1-14
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00

#SBATCH --mail-user=yixiangluo@berkeley.edu
#SBATCH --mail-type=ALL

EXPERIMENT=test

# run experiments

ml R

Rscript ../R/cluster/experiment.R $EXPERIMENT $SLURM_ARRAY_TASK_ID


