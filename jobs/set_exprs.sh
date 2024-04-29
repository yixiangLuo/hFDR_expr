#!/bin/bash

#SBATCH --job-name=set_expr
#SBATCH --output=../log/%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=06:00:00

#SBATCH --mail-user=yixiangluo@berkeley.edu
#SBATCH --mail-type=ALL

N_JOBS=200
EXPR_NAMES=$(./expr_names.sh)

ml R

for EXPR_NAME in $EXPR_NAMES; do
    Rscript ../R/cluster/set_expr.R $EXPR_NAME $N_JOBS
done

