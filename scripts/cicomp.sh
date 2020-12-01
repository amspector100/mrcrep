#!/bin/bash

# Set to > 512 to replicate
NREPS=35
NPROCESSES=1

ARGS_EQUI="
	--p_dgp 500
	--n_dgp [375]
	--coeff_size_dgp [1]
	--method_dgp blockequi
	--rho_dgp start0.75end0.95numvals5
	--gamma_dgp 0
	--block_size_dgp 2
	--fstat_filter [lasso,ridge]
	--sparsity_dgp 0.1
	--cv_score_fstat 0 
	--seed_start 0
	--num_processes $NPROCESSES
	--reps $NREPS
	--coeff_dist_dgp none
	--comp_ci True
	--description CICOMP1: Lasso,Ridge on Block equicorr with p500
"

# Block-equicorrelated panel
python3 ../mrcrep.py $ARGS_EQUI

ARGS_ER="
	--p_dgp 500
	--method_dgp ver
	--delta_dgp [0.2,0.4,0.6,0.8]
	--fstat_filter [lasso,ridge]
	--sparsity_dgp 0.1
	--cv_score_fstat 0
	--seed_start 0
	--num_processes $NPROCESSES
	--reps $NREPS
	--n_dgp 333
	--coeff_size_dgp [1]
	--coeff_dist_dgp uniform
	--comp_ci True
	--description Lasso,Ridge on VER with p500
"

# ErdosRenyi (cov)
python3 ../mrcrep.py $ARGS_ER