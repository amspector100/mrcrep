#!/bin/bash

#Set nreps=256
NREPS=1
NPROCESSES=1

ARGS="
	--p_dgp 500
	--n_sample [250,500,750,1000]
	 --method_dgp [blockequi]
	 --rho_dgp start0end0.9numvals10
	 --gamma_dgp 1 
	 --fstat_filter [lasso,ridge]
	 --y_dist_sample [gaussian,binomial]
	 --sparsity_dgp 0.1
	 --cv_score_fstat 0
	 --seed_start 0 
	 --num_processes $NPROCESSES
	 --reps $NREPS
	 --coeff_size_dgp 1
	 --coeff_dist_dgp uniform
	 --description Main equicorrelated plot
"
python3 ../mrcrep.py $ARGS