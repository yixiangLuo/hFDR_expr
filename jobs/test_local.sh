#!/bin/bash

N_JOBS=14
EXPERIMENT=test

# set experiment data
Rscript ../R/cluster/set_expr.R $EXPERIMENT $N_JOBS

# run experiments
for JOB_ID in $(seq 1 1 $N_JOBS)
do
    Rscript ../R/cluster/experiment.R $EXPERIMENT $JOB_ID &
done

wait

Rscript ../R/cluster/post_process.R $EXPERIMENT
