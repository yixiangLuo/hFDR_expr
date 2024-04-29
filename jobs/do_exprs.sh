#!/bin/bash

EXPR_NAMES=$(./expr_names.sh)

for EXPR_NAME in $EXPR_NAMES; do
    sbatch "${EXPR_NAME}_job.sh"
done

