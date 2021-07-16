#!/bin/bash

# Increase as desired
NREPS=1
NPROCESSES=1

ARGS_VARIED_RHO="--p_dgp 500 
        --method_dgp [ar1]
        --max_corr_dgp start0.5end1numvals6
        --a_dgp [3]
        --fstat_filter [lasso,ridge]
        --coeff_size_dgp 2
        --reps $NREPS
        --num_processes $NPROCESSES 
        --cv_score_fstat 0
        --seed_start 0
        --coeff_dist_dgp uniform
        --corr_signals_dgp [True,False]
        --sparsity_dgp 0.1
        --description Linear stats on ar1 with varied rho + maximum correlation
"

python3.9 ../mrcrep.py $ARGS_VARIED_RHO
