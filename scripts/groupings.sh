#!/bin/bash

# Increase these as desired
NREPS=1
NPROCESSES=1

COMMON_ARGS="
        --p_dgp 500
        --fstat_filter [lasso,ridge]
        --sparsity_dgp 0.1
        --cv_score_fstat 0
        --seed_start 0
        --num_processes $NPROCESSES 
        --reps $NREPS
        --coeff_dist_dgp uniform
        --n_dgp [125,250,333]
"

ARGS_AR1="${COMMON_ARGS}
        --p_dgp 500 
        --method_dgp [ar1]
        --cutoff_dgp [1,0.9,0.8,0.7,0.6]
        --a_dgp [3]
        --coeff_size_dgp 2
        --corr_signals_dgp [True,False]
        --description Group knockoffs + linear statistics on ar1 
"

### Actually run
python3.9 ../mrcrep.py $ARGS_AR1

