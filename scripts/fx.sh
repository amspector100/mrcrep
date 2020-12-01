#!/bin/bash

# ---- FX knockoffs
NREPS=1
NPROCESSES=1

# 1. Equicorrelated and block equicorrelated
EQUI_ARGS="
	--p_dgp 500
	--method_dgp [blockequi]
	--rho_dgp [0.5]
	--gamma_dgp [0,1]
	--fstat_filter [lasso]
	--sparsity_dgp 0.1
	--cv_score_fstat 0
	--seed_start 0
	--zstat_fstat lars_path
	--num_processes $NPROCESSES
	--reps $NREPS
	--coeff_size_dgp 0.15
	--coeff_dist_dgp uniform
	--n_sample [1005,1250,1500,1750,2000]
	--ksampler_filter fx
	--description lasso signed max on equicorr\block equicorr with p500, linear cond_mean, FX knockoffs and low coeff_size
"
python3 ../mrcrep.py $EQUI_ARGS

# 2. AR1
AR1_ARGS="
	--p_dgp 500
	--method_dgp [ar1]
	--a_dgp [3] 
	--fstat_filter [lasso]
	--sparsity_dgp 0.1
	--cv_score_fstat 0
	--seed_start 0
	--zstat_fstat lars_path
	--num_processes $NPROCESSES
	--reps $NREPS
	--coeff_size_dgp 0.45
	--coeff_dist_dgp uniform
	--corr_signals_dgp [True,False]
	--n_sample [1005,1250,1500,1750,2000]
	--ksampler_filter fx
	--description Lasso on AR1 with p500, linear cond_mean, FX knockoffs, low coeff size
"
python3 ../mrcrep.py $AR1_ARGS

# 3. ER Covariance runs
VER_ARGS="
	--p_dgp 500
	--method_dgp [ver]
	--delta_dgp [0.2]
	--fstat_filter [lasso]
	--sparsity_dgp 0.1
	--cv_score_fstat 0
	--seed_start 0
	--zstat_fstat lars_path
	--num_processes $NPROCESSES
	--reps $NREPS
	--coeff_size_dgp 0.5
	--coeff_dist_dgp uniform
	--n_sample [1005,1250,1500,1750,2000]
	--ksampler_filter fx
	--description Linear stats on VER with p500, linear cond_mean, FX knockoffs
"
python3 ../mrcrep.py $VER_ARGS

# 4. ER Precision runs
QER_ARGS="
	--p_dgp 500
	--method_dgp [qer]
	--delta_dgp [0.2]
	--fstat_filter [lasso]
	--sparsity_dgp 0.1
	--cv_score_fstat 0
	--seed_start 0
	--zstat_fstat lars_path
	--num_processes $NPROCESSES
	--reps $NREPS
	--coeff_size_dgp 0.25
	--coeff_dist_dgp uniform
	--n_sample [1005,1250,1500,1750,2000]
	--ksampler_filter fx
	--description Linear stats on QER with p500, linear cond_mean, FX knockoffs
"
python3 ../mrcrep.py $QER_ARGS