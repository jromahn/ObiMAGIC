#!/usr/bin/env bash

CONFIGS=("$@")
CONF=${CONFIGS[${SLURM_ARRAY_TASK_ID}-1]}

03_ObiWizard_main.sh -obiwizard-config $CONF
