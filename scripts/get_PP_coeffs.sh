#!/bin/bash

K=$(($SLURM_ARRAY_TASK_ID))

Rscript get_PP_coeffs.R $K
