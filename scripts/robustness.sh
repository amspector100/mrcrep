#!/bin/bash

NREPS=1
NPROCESSES=1

COMMON_ARGS="
	--p_dgp 500
	--fstat_filter [lasso]
	--sparsity_dgp 0.1
	--cv_score_fstat 0
	--seed_start 0
	--num_processes $NPROCESSES
	--reps $NREPS
	--coeff_dist_dgp uniform
	--n_sample [200,350,500,650,875,1000]
	--infer_sigma_filter True 
"

BLOCK_EQUI_ARGS="${COMMON_ARGS}
	--method_dgp [blockequi]
	--rho_dgp [0.6]
	--gamma_dgp 0
	--coeff_size_dgp 0.5
	--shrinkage_filter [none,graphicallasso,ledoitwolf]
	--description Lasso on Block-equicorr with p500, linear cond_mean, but infer Cov matrix, low coeff size
"
python3 ../mrcrep.py $BLOCK_EQUI_ARGS

# NO graphical lasso for equicorrelated
EQUI_ARGS="${COMMON_ARGS}
	--method_dgp [blockequi]
	--rho_dgp [0.6]
	--gamma_dgp 1
	--coeff_size_dgp 0.5
	--shrinkage_filter [none,ledoitwolf]
	--description Lasso on equicorr with p500, linear cond_mean, but infer Cov matrix, low coeff size
"
python3 ../mrcrep.py $EQUI_ARGS

# AR1 case
AR1_ARGS="${COMMON_ARGS}
	--method_dgp [ar1]
	--a_dgp [3]
	--coeff_size_dgp 1
	--corr_signals_dgp [True,False]
	--shrinkage_filter [none,graphicallasso,ledoitwolf]
	--description Lasso on AR1 with p500, linear cond_mean, but infer Cov matrix, low coeff size
"
python3 ../mrcrep.py $AR1_ARGS


# ER case
ER_ARGS="${COMMON_ARGS}
	--p_dgp 500
	--method_dgp [ver,qer]
	--delta_dgp [0.2]
	--coeff_size_dgp 0.5
	--shrinkage_filter [none,graphicallasso,ledoitwolf]
	--description Lasso on ER/VER with p500, linear cond_mean, but infer Cov matrix, low coeff size
"
python3 ../mrcrep.py $ER_ARGS

