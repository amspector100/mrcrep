#!/bin/bash

# Set to >128
NREPS=1
NPROCESSES=1

COMMON_ARGS="
	--p_dgp 200
	--fstat_filter [lasso,randomforest,deeppink]
	--sparsity_dgp 0.15
	--seed_start 0
	--cv_score_fstat 0
	--num_processes $NPROCESSES
	--reps $NREPS
	--coeff_size_dgp 5
	--cond_mean_sample [cos,quadratic,pairint,cubic,trunclinear]
	--n_sample [200,300,500,750,1000,2000,3000]
"
# --- Blockequi and equicorrelated (gamma=1 means equicorrelated)
EQUI_ARGS="${COMMON_ARGS}
	--method_dgp [blockequi]
	--rho_dgp [0.5]
	--gamma_dgp [0,1]
	--description Lasso + nonlinear fstats on equi p200 large n, all nonlinear cond means
"
python3 ../mrcrep.py $EQUI_ARGS

# --- AR1
AR1_ARGS="${COMMON_ARGS}
	--method_dgp ar1
	--a_dgp 3
	--corr_signals_dgp [True,False]
	--description Lasso/nonlinear stats on AR1 p200 large n, all nonlinear cond means
"
python3 ../mrcrep.py $AR1_ARGS

# --- VER/QER
ER_ARGS="${COMMON_ARGS}
	--method_dgp [ver,qer]
	--delta_dgp [0.2]	
"
python3 ../mrcrep.py $ER_ARGS

module purge

