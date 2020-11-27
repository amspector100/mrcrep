#!/bin/bash

# Set nreps=128
NREPS=1
NPROCESSES=1

LOW_DIM_ARGS="
	--p_dgp 500
	--n_dgp 500
	--n_sample 500
	--method_dgp [blockequi]
	--rho_dgp [0,0.8]
	--gamma_dgp 1
	--fstat_filter lasso
	--use_lars_fstat True
	--y_dist_sample gaussian
	--sparsity_sample 0.1
	--cv_score_fstat 0
	--seed_start 0
	--num_processes $NPROCESSES 
	--reps $NREPS 
	--coeff_dist_sample dsliu2020 
	--resample_beta True
	--description DSLIU2020 REPLICATION p500
"

# -- Low dimensional case
python3 ../mrcrep.py $LOW_DIM_ARGS

HIGH_DIM_ARGS="
	--p_dgp 2000
	--n_dgp 500
	--n_sample 500
	--method_dgp [daibarber2016]
	--rho_dgp [0,0.8] 
	--gamma_dgp 1
	--fstat_filter lasso
	--use_lars_fstat True
	--y_dist_sample gaussian
	--sparsity_sample 0.025
	--cv_score_fstat 0
	--seed_start 0
	--num_processes $NPROCESSES
	--reps $NPROCESSES
	--coeff_dist_sample dsliu2020 
	--resample_beta True 
	--description DSLIU2020 REPLICATION p2000
"

python3 ../mrcrep.py $HIGH_DIM_ARGS
