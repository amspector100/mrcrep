#!/bin/bash

# Set to > 128
NPROCESSES=1
NREPS=1

# -- Linear stuff, in particular AR1
COMMON_ARGS="
	--p_dgp 500
	--fstat_filter [lasso,dlasso,ridge]
	--sparsity_dgp 0.1
	--cv_score_fstat 0
	--seed_start 0
	--num_processes $NPROCESSES 
	--reps $NREPS
	--coeff_dist_dgp uniform
	--n_dgp [125,250,333,500,1000]
"
EQUI_ARGS="${COMMON_ARGS}
	--method_dgp [blockequi]
	--rho_dgp [0.5]
	--gamma_dgp [0,1]
	--coeff_size_dgp 1
	--description Block equicorrelated/equicorrelated gaussian plots
"
python3 ../mrcrep.py $EQUI_ARGS

ER_ARGS="${COMMON_ARGS}
	--method_dgp [ver,qer]
	--delta_dgp [0.2]
	--coeff_size_dgp 1
	--description VER/QER gaussian plots
"
python3 ../mrcrep.py $ER_ARGS

AR1_ARGS="${COMMON_ARGS}
	--method_dgp [ar1]
	--a_dgp [3]
	--coeff_size_dgp 2
	--corr_signals_dgp [True,False]
	--description AR1 gaussian plots
"
python3 ../mrcrep.py $AR1_ARGS