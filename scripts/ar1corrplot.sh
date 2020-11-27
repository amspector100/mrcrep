#!/bin/bash

# Change this to > 128 to replicate
NREPS=1
NPROCESSES=1

ARGS_CONST_RHO="--p_dgp 500 
	--method_dgp [ar1] 
	--rho_dgp start0.1end0.9numvals9
	--fstat_filter [lasso,ridge,dlasso]
	--coeff_size_dgp 1 
	--reps $NREPS
	--num_processes $NPROCESSES
	--cv_score_fstat 0
	--seed_start 0
	--coeff_dist_dgp uniform
	--corr_signals_dgp [True,False]
	--sparsity_dgp 0.1
	--description Linear statistics on ar1 with constant rho
"
python3 ../mrcrep.py $ARGS_CONST_RHO

ARGS_VARIED_RHO="--p_dgp 500 
	--method_dgp [ar1] 
	--a_dgp [0.5,1,2,3]
	--fstat_filter [lasso,ridge,dlasso]
	--coeff_size_dgp 2
	--reps $NREPS
	--num_processes $NPROCESSES 
	--cv_score_fstat 0
	--seed_start 0
	--coeff_dist_dgp uniform
	--corr_signals_dgp [True,False]
	--sparsity_dgp 0.1
	--description Linear statistics on ar1 with constant rho
"
python3 ../mrcrep.py $ARGS_VARIED_RHO