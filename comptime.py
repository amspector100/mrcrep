# To profile, run python3.9 -m cprofilev comptime.py

import os
import sys
import time
import datetime


import warnings
import numpy as np
import pandas as pd

# Ensure we are using the right version of knockpy
import knockpy
print(f"Using version {knockpy.__version__}, should be using v1.1.0")


#### GLOBALS
NUM_PROCESSES = 12 # change to whatever you want
SPACING = 200
REPS = 12
MAX_BLOCK = 100
NUM_FACTORS = 25


def calc_mac(V, S):
	""" Mac loss function """
	return S.shape[0] - np.mean(np.diag(S))

def test_timing(seed, loss_fn, sample_kwargs, smatrix_kwargs):
	"""
	Tests how long it creates to compute s-matrix
	"""
	np.random.seed(seed)
	dgprocess = knockpy.dgp.DGP()
	dgprocess.sample_data(**sample_kwargs)
	Sigma = dgprocess.Sigma
	time0 = time.time()
	S = knockpy.smatrix.compute_smatrix(
		Sigma=Sigma, **smatrix_kwargs
	)
	tdelta = time.time() - time0
	loss = loss_fn(Sigma, S)
	return (tdelta, seed, loss)

def main():


	# Create output paths
	today = str(datetime.date.today())
	hour = str(datetime.datetime.today().time())
	hour = hour.replace(':','-').split('.')[0]
	output_path = f'data/{today}/{hour}/comptime.csv'
	output_dir = os.path.dirname(output_path)
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	time_0 = time.time()
	ps = SPACING * np.arange(1, 10)
	outputs = []
	#output = pd.DataFrame(columns = ['method', 'p', 'seed', 'time'])
	for p in ps:
		for (method, solver, loss_fn) in [
			('mvr', 'cd', knockpy.mrc.mvr_loss), 
			('maxent', 'cd', knockpy.mrc.maxent_loss),
			('sdp', 'sdp', calc_mac),
			('sdp', 'cd', calc_mac),
		]:
			for (how_approx, max_block) in [
				#('blockdiag', 1001 if p <= 1000 else 10),
				('blockdiag', MAX_BLOCK),
				#('factor', MAX_BLOCK)
			]:
				if how_approx == 'factor' and method == 'sdp' and solver == 'sdp':
					continue
				if max_block == 10:
					continue
				total_time = np.around(time.time() - time_0, 2)
				print(f"At p={p}, method={method}, solver={solver}, how_approx={how_approx}, max_block={max_block}, total_time={total_time}")
				tdeltas = knockpy.utilities.apply_pool(
					func=test_timing,
					constant_inputs={
						'sample_kwargs':{'p':p, 'method':'ar1', 'a':1, 'b':1},
						'smatrix_kwargs':{
							'method':method, 
							'solver':solver,
							'how_approx':how_approx,
							'max_block':max_block,
							'num_factors':min(NUM_FACTORS, int(p/2)),
							'line_search':False,
						},
						'loss_fn':loss_fn,
					},
					seed=list(range(REPS)),
					num_processes=NUM_PROCESSES
				)
				for (tdelta, seed, loss) in tdeltas:
					outputs.append(
						[p, method, solver, how_approx, max_block, seed, tdelta, loss]
					)

				output_df = pd.DataFrame(outputs, columns=[
					'p', 'method', 'solver', 'how_approx', 'max_block', 'seed', 'time', 'loss'
				])

				output_df.to_csv(output_path)

if __name__ == '__main__':
	main()