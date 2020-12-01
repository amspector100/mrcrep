#!/bin/bash

# Set nreps>300
NREPS=1
NPROCESSES=1

# -- Linear stuff, in particular AR1
COMMON_ARGS="
	--p_dgp 500
	--n_sample [250]
	--method_dgp [blockequi]
	--gamma_dgp 1
	--fstat_filter [lasso,ridge]
	--sparsity_dgp 0.1
	--cv_score_fstat 0
	--seed_start 1000
	--num_processes $NPROCESSES
	--reps $NREPS
	--coeff_dist_dgp uniform
	--S_curve True
	--description Vary S in equiplot
"

# All rho except rho = 0.9
ARGS1="${COMMON_ARGS}
	--rho_dgp start0.1end0.7numvals4
	--coeff_size_dgp 1
"
python3 ../mrcrep.py $ARGS1


# rho = 0.9 gets a higher signal size
ARGS2="${COMMON_ARGS}
	--rho_dgp 0.9
	--coeff_size_dgp 2
"
python3 ../mrcrep.py $ARGS2