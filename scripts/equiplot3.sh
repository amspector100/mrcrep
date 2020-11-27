#!/bin/bash

# Set to 128 or greater
NREPS=1
NPROCESSES=1

ARGS="
	--p_dgp 500
	--n_sample 2000
	--sign_prob_dgp start0.0end0.5numvals6
	--method_dgp blockequi
	--rho_dgp 0.6
	--gamma_dgp 1
	--fstat_filter lasso
	--y_dist_sample gaussian
	--sparsity_dgp 0.5
	--cv_score_fstat 0
	--seed_start 0
	--num_processes $NPROCESSES
	--reps $NREPS
	--use_lars_fstat True
	--coeff_size_dgp 5
	--coeff_dist_dgp none
	--iid_signs_dgp False
	--description Q conjecture for equicorrelated graph, vary q, sign_prob
"
python3 ../mrcrep.py $ARGS
