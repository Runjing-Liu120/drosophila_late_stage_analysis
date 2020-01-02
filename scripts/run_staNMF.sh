#!/bin/bash

source activate staNMF-py27

seed=$((492870 + $SLURM_ARRAY_TASK_ID))
K=$(($SLURM_ARRAY_TASK_ID))

python run_staNMF.py --seed $seed --k1 $K --k2 $K
