#!/bin/bash

# Set nreps>300
NREPS=1
NPROCESSES=1

# -- Linear stuff, in particular AR1
ARGS="
	--p_dgp 500
	--n_sample [250]
	--method_dgp [blockequi]
	--rho_dgp start0.1end0.9numvals5
	--gamma_dgp 1
	--fstat_filter [lasso,ridge]
	--sparsity_dgp 0.1
	--cv_score_fstat 0
	--seed_start 1000
	--num_processes $NPROCESSES
	--reps $NREPS
	--coeff_size_dgp 1
	--coeff_dist_dgp uniform
	--S_curve True
	--description Vary S in equiplot
"
python3 ../mrcrep.py $ARGS

