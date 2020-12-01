#!/bin/bash

NREPS=1
NPROCESSES=1

# Heavy-tailed Markov Chain
AR1T_ARGS="
	--p_dgp 500
	--method_dgp [ar1]
	--a_dgp [3]
	--fstat_filter [lasso,dlasso,ridge]
	--sparsity_dgp 0.1
	--cv_score_fstat 0
	--seed_start 0
	--num_processes $NPROCESSES
	--reps $NREPS
	--coeff_size_dgp 0.3
	--coeff_dist_dgp uniform
	--corr_signals_dgp [True,False]
	--n_sample [350,500,750,1000]
	--x_dist_sample ar1t
	--df_t_dgp 3
	--ksampler_filter [artk,gaussian]
	--description Lasso,Ridge on AR1T with p500, linear cond_mean, t_dist design 
"
python3 ../mrcrep.py $AR1T_ARGS

# Block-equicorrelated t distribution (independent blocks)
BLOCKT_ARGS="
	--p_dgp 500
	--method_dgp blockequi
	--gamma_dgp 0
	--rho_dgp [0.5]
	--fstat_filter [lasso,dlasso,ridge]
	--sparsity_dgp 0.1
	--cv_score_fstat 0
	--seed_start 0
	--num_processes $NPROCESSES
	--reps $NREPS 
	--coeff_size_dgp 0.4
	--coeff_dist_dgp uniform
	--n_sample [350,500,750,1000]
	--x_dist_sample blockt 
	--df_t_dgp 3
	--df_t_sample 3
	--ksampler_filter [blockt,gaussian]
	--description Lasso,Ridge on block-equi t with p500, linear cond_mean
"
python3 ../mrcrep.py $BLOCKT_ARGS

# Gibbs grid --- this one is the most computationally intensive
GIBBS_GRID_ARGS="
	--p_dgp 625
	--x_dist_dgp gibbs
	--x_dist_sample gibbs
	--fstat_filter [lasso,ridge,dlasso]
	--sparsity_dgp 0.1 
	--cv_score_fstat 0
	--seed_start 0
	--num_processes $NPROCESSES
	--reps $NREPS
	--coeff_size_dgp 0.4
	--coeff_dist_dgp uniform
	--n_sample [350,500,750,1000]
	--ksampler_filter [gibbs_grid,gaussian]
	--resample_beta True
	--resample_sigma False
	--description Linear statistics on gibbs grid, linear cond mean
"
python3 ../mrcrep.py $GIBBS_GRID_ARGS
