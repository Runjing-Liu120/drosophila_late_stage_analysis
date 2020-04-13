#!/bin/bash
seed=$(($SLURM_ARRAY_TASK_ID + 454534))
sample=$SLURM_ARRAY_TASK_ID
Rscript run_iRF.R 0.5 $seed ensemble_iRF/irfSpatialFit_gut_set12pp_sample$sample
