#!/usr/bin/env bash

DIR=$1
CONF=$(find $DIR -name "*_obiwizard_*ini" | sed -n ${SLURM_ARRAY_TASK_ID}p)

03_ObiWizard_main.sh -obiwizard-config $CONF
