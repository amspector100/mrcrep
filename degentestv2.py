import os
import sys
import time
import datetime
import knockpy
from knockadapt import knockoff_stats as kstats
from knockadapt import utilities
from knockadapt.knockoff_filter import KnockoffFilter

import warnings
import numpy as np
import pandas as pd

import itertools
from multiprocessing import Pool
from functools import partial

# Global: the set of antisymmetric functions we use
PAIR_AGGS = ['cd', 'sm']
# q-values to evaluate for
q_values = 0.05*np.arange(1, 21)
# Degrees of freedom for t-distributions
DEFAULT_DF_T = knockadapt.graphs.DEFAULT_DF_T

def fetch_kwarg(kwargs, key, default=None):
	""" Utility function for parsing """
	if key in kwargs:
		return kwargs.pop(key)
	else:
		return default

def str2bool(v):
	""" Helper function, converts strings to boolean vals""" 
	if isinstance(v, bool):
		return v
	elif isinstance(v, str):
		if v.lower() in ['yes', 'true', 't', 'y', '1']:
			return True
		elif v.lower() in ['no', 'false', 'f', 'n', '0']:
			return False
	return v

def obj2int(v):
	try:
		v = float(v)
		if v.is_integer():
			v = int(v)
		return v
	except: 
		return str2bool(v)

def val2list(val):
	""" Turns discrete values into lists, otherwise returns """
	if not isinstance(val, list):
		if isinstance(val, np.ndarray):
			if len(val.shape) == 1:
				return val
		return [val] 
	return val

def dict2keyproduct(dictionary):
	""" Takes a dictionary mapping to lists
	and returns the sorted list of keys and 
	the cartesian product of each list."""
	keys = sorted([key for key in dictionary])
	components = []
	for key in keys:
		dictionary[key] = val2list(dictionary[key])
		components.append(dictionary[key])
	product = list(itertools.product(*components))
	return keys, list(product)

def apply_pool(func, all_inputs, num_processes):
	""" Utility function"""

	# Don't use the pool object if num_processes=1
	if num_processes == 1:
		all_outputs = []
		for inp in all_inputs:
			all_outputs.append(func(inp))
	else:
		with Pool(num_processes) as thepool:
			all_outputs = thepool.map(func, all_inputs)

	return all_outputs

def fetch_competitor_S(
	Sigma,
	groups,
	time0,
	S_curve=False,
	verbose=False,
	**kwargs
):
	"""
	Creates various different S-matrices.
	"""
	# If S_curve is true, just do gamma * I for many gammas
	p = Sigma.shape[0]
	fetch_kwarg(kwargs, 'sdp_tol', None) # Prevent double arg errors
	if S_curve:
		mineig = np.linalg.eigh(Sigma)[0].min()
		gammas = 2*mineig*np.arange(0, 11, 1)/10
		gammas[0] += 0.001
		gammas[-1] -= 0.001
		S_matrices = {
			f'S{np.around(gamma, 3)}':gamma*np.eye(p) for gamma in gammas
		}
		S_matrices['mvr'] = knockadapt.knockoffs.compute_S_matrix(
			Sigma=Sigma,
			groups=groups,
			method='mvr',
			solver='cd',
		)
		S_matrices['maxent'] = knockadapt.knockoffs.compute_S_matrix(
			Sigma=Sigma,
			groups=groups,
			method='maxent',
			solver='cd'
		)
		S_matrices['sdp'] = knockadapt.knockoffs.compute_S_matrix(
			Sigma=Sigma,
			groups=groups,
			method='sdp',
			tol=1e-10,
		)
		S_matrices['ci'] = knockadapt.knockoffs.compute_S_matrix(
			Sigma=Sigma,
			groups=groups,
			method='ciknock',
		)
		# Create perturbed option when G_SDP is low rank
		if Sigma[0,1] >= 0.5:
			S_matrices['sdp_perturbed'] = 0.99*S_matrices['sdp']
		return S_matrices

	### Special case: detect if Sigma is equicorrelated,
	# in which case we can calculate the solution analytically
	# for ungrouped knockoffs.
	if Sigma is not None:
		if np.unique(groups).shape[0] == p:
			rho = Sigma[0, 1]
			if rho >= 0:
				equicorr = rho*np.ones((p, p)) + (1-rho)*np.eye(p)
				if np.all(np.abs(Sigma-equicorr) < 1e-10):
					print(f"Sigma is equicorr (rho={rho}), using analytical solution")
					S_SDP = min(1, 2-2*rho)*np.eye(p)
					S_matrices = {'sdp':S_SDP}
					# Compute exact solution for smaller p
					if p < 1000:
						S_matrices['mvr'] = knockadapt.knockoffs.compute_S_matrix(
							Sigma=Sigma,
							groups=groups,
							method='mvr',
							solver='cd',
						)
						S_matrices['mmi'] = knockadapt.knockoffs.compute_S_matrix(
							Sigma=Sigma,
							groups=groups,
							method='maxent',
							solver='cd'
						)
						S_matrices['ci'] = knockadapt.knockoffs.compute_S_matrix(
							Sigma=Sigma,
							method='ciknock',
						)
					# Else, asymptotically these methods are the same
					else:
						S_matrices['mvr'] = (1-rho)*np.eye(p)
					if rho < 0.5:
						return S_matrices
					else:
						S_matrices['sdp_perturbed'] = S_SDP*(0.99)
						return S_matrices

	### Calculate (A)SDP S-matrix
	if time0 is None:
		time0 = time.time() 

	# Compute SDP matrix
	S_matrices = {}
	if verbose:
		(f'Now computing SDP S matrices, time is {time.time() - time0}')
	S_SDP = knockadapt.knockoffs.compute_S_matrix(
		Sigma=Sigma,
		groups=groups,
		max_block=501, 
		sdp_tol=1e-5,
		method='sdp' if p <= 500 else 'asdp',
		**kwargs
	)
	S_matrices['sdp'] = S_SDP
	if verbose:
		print(f'Finished computing SDP matrices, time is {time.time() - time0}')

	### Calculate MCR matrices
	rej_rate = fetch_kwarg(kwargs, 'rej_rate', 0)
	if rej_rate != 0:
		kwargs['rej_rate'] = rej_rate
	for new_method in ['mvr','maxent']:

		# Choose between projected gradient decent and coord descent
		# solvers
		if new_method == 'maxent' and rej_rate != 0:
			solver = 'psgd'
		else:
			solver = 'cd'
		S_MRC = knockadapt.knockoffs.compute_S_matrix(
			Sigma=Sigma,
			groups=groups,
			method=new_method,
			solver=solver,
			**kwargs
		)
		S_matrices[new_method] = S_MRC
	S_matrices['ci'] = knockadapt.knockoffs.compute_S_matrix(
		Sigma=Sigma,
		method='ciknock',
	)
	if verbose:
		print(f'Finished computing MRC matrices, time is {time.time() - time0}')

	return S_matrices

def Z2selections(Z, groups, q, **kwargs):
	
	# Calculate W statistics
	W = kstats.combine_Z_stats(Z, groups, **kwargs)

	# Calculate selections 
	T = kstats.data_dependent_threshhold(W=W, fdr=q)
	selected_flags = (W >= T).astype("float32")
	return selected_flags, W

def selection2power(selections, group_nonnulls):
	# Calculate fdp, power
	fdp = np.sum(selections*(1-group_nonnulls))
	power = np.sum(selections*group_nonnulls)

	# Normalize
	fdp = fdp/max(1, np.sum(selections))
	power = power/max(1, np.sum(group_nonnulls))

	return [power, fdp]

def single_dataset_power_fdr(
	seed,
	Sigma,
	beta,
	groups,
	normalize=True,
	sample_kwargs={},
	filter_kwargs={
		'feature_stat_kwargs':{},
		'knockoff_kwargs':{},
	},
	S_matrices=None,
	time0=None,
	S_curve=False,
):
	""" 
	Knockoff kwargs should be included in filter_kwargs
	"""

	# Fetch groups, infer q
	if Sigma is not None:
		p = Sigma.shape[0]
	else:
		p = sample_kwargs['p']
	if groups is None:
		groups = np.arange(1, p+1, 1)

	# Prevent global changes to kwargs
	filter_kwargs = filter_kwargs.copy()
	sample_kwargs = sample_kwargs.copy()

	# Sample data, record time
	localtime = time.time()
	np.random.seed(seed)
	X, y, beta, _, Sigma = knockadapt.graphs.sample_data(
		corr_matrix=Sigma,
		beta=beta,
		**sample_kwargs
	)

	# Parse nonnulls (this allows for group knockoffs although
	# we do not actually use this in our experiments)
	group_nonnulls = knockadapt.utilities.fetch_group_nonnulls(
		beta, groups
	)

	# Some updates for fixedX knockoffs
	# and MX knockoffs when Sigma must be inferred
	fixedX = fetch_kwarg(filter_kwargs, 'fixedx', default=False)
	infer_sigma = fetch_kwarg(filter_kwargs, 'infer_sigma', default=False)

	# For the metro sampler, we can compute better S matrices if we have a 
	# guess of the rejection rate. We do not use this in our experiments
	# by default.
	rej_rate = fetch_kwarg(filter_kwargs, 'rej_rate', default=0)

	# We calculate S-matrices if we do not already know them.
	if 'knockoff_kwargs' in filter_kwargs:
		kwargs = filter_kwargs['knockoff_kwargs']
	else:
		kwargs = {}
	# Case 1: We have to estimate Sigma from the data
	if infer_sigma:
		shrinkage = fetch_kwarg(filter_kwargs, 'shrinkage', default='ledoitwolf')
		Sigma, _ = knockadapt.utilities.estimate_covariance(X, shrinkage=shrinkage)
		Sigma = utilities.cov2corr(Sigma)
		invSigma = utilities.chol2inv(Sigma)
	# Case 2: We are running FX knockoffs
	if fixedX:
		Sigma = np.dot(X.T, X)
		invSigma = None
	# Calculate S-matrices from Sigma
	if infer_sigma or fixedX or S_matrices is None:
		verbose = fetch_kwarg(kwargs, 'verbose', default=False)
		S_matrices = fetch_competitor_S(
			Sigma=Sigma, 
			groups=groups,
			time0=time0,
			rej_rate=rej_rate,
			verbose=verbose,
			S_curve=S_curve,
			**kwargs
		)

	# Now we loop through the S matrices
	degen_flag = 'sdp_perturbed' in S_matrices 
	output = []
	for S_method in S_matrices:

		# Pull S matrix
		S = S_matrices[S_method]

		# If the method produces fully degenerate knockoffs,
		# signal this as part of the filter kwargs.
		_sdp_degen = (degen_flag and S_method == 'sdp')

		# Do NOT run OLS or debiased lasso for degenerate case,
		# because this will lead to linear algebra errors.
		if _sdp_degen and 'feature_stat' in filter_kwargs:
			if filter_kwargs['feature_stat'] == 'dlasso':
				for pair_agg in PAIR_AGGS:
					for q in q_values:
						output.append([
							S_method,
							0,
							0,
							q,
							0,
							"NULL",
							pair_agg,
							np.zeros(p),
							np.zeros(p),
							np.zeros(p),
							np.zeros(p),
							seed
						])
				continue

		# Create knockoff_kwargs
		if 'knockoff_kwargs' not in filter_kwargs:
			filter_kwargs['knockoff_kwargs'] = {}
		filter_kwargs['knockoff_kwargs']['S'] = S
		filter_kwargs['knockoff_kwargs']['method'] = S_method
		filter_kwargs['knockoff_kwargs']['_sdp_degen'] = _sdp_degen
		if 'verbose' not in filter_kwargs['knockoff_kwargs']:
			filter_kwargs['knockoff_kwargs']['verbose'] = False

		# Possibly pass a few parameters to metro sampler
		if 'x_dist' in sample_kwargs:
			# For Ising
			if 'gibbs_graph' in sample_kwargs:
				filter_kwargs['knockoff_kwargs']['gibbs_graph'] = sample_kwargs['gibbs_graph']
			# For t-distributions
			if str(sample_kwargs['x_dist']).lower() in ['ar1t', 'blockt']:
				if 'df_t' in sample_kwargs:
					filter_kwargs['knockoff_kwargs']['df_t'] = sample_kwargs['df_t']
				else:
					filter_kwargs['knockoff_kwargs']['df_t'] = DEFAULT_DF_T # This matters

		# Run MX knockoff filter to obtain
		# Z statistics
		knockoff_filter = KnockoffFilter(fixedX=fixedX)
		knockoff_filter.forward(
			X=X, 
			y=y, 
			mu=np.zeros(p),
			Sigma=Sigma, 
			groups=groups,
			**filter_kwargs
		)
		Z = knockoff_filter.Z
		score = knockoff_filter.score
		score_type = knockoff_filter.score_type

		# Calculate power/fdp/score for a variety of 
		# antisymmetric functions
		for pair_agg in PAIR_AGGS:
			for q in q_values:
				# Start by creating selections
				selections, W = Z2selections(
					Z=Z,
					groups=groups,
					q=q,
					pair_agg=pair_agg
				)
				# Then create power/fdp
				power, fdp = selection2power(
					selections, group_nonnulls
				)
				output.append([
					S_method,
					power,
					fdp,
					q,
					score, 
					score_type,
					pair_agg,
					np.around(W, 4),
					np.around(Z[0:p], 4),
					np.around(Z[p:], 4),
					selections,
					seed,
				])

	# Possibly log progress
	try:
		if seed % 10 == 0 and sample_kwargs['n'] == MAXIMUM_N:
			overall_cost = time.time() - time0
			local_cost = time.time() - localtime
			print(f"Finished one seed {seed}, took {local_cost} per seed, {overall_cost} total")
	except:
		# In notebooks this will error
		pass

	# Output: list of 
	# [S_method, power, fdp, q, score, score_type, pair_agg, W, Z, tildeZ, selection, seed]
	return output

def calc_power_and_fdr(
	time0,
	Sigma,
	beta,
	groups=None,
	q=0.2,
	reps=100,
	num_processes=1,
	sample_kwargs={},
	filter_kwargs={},
	seed_start=0,
	S_matrices={'sdp':None, 'mvr':None},
	S_curve=False,
):

	# Fetch nonnulls
	if Sigma is not None:
		p = Sigma.shape[0]
	else:
		p = sample_kwargs['p']

	# Sample data reps times and calc power/fdp
	partial_sd_power_fdr = partial(
		single_dataset_power_fdr, 
		Sigma=Sigma,
		beta=beta,
		groups=groups,
		sample_kwargs=sample_kwargs,
		filter_kwargs=filter_kwargs,
		S_matrices=S_matrices,
		time0=time0,
		S_curve=S_curve,
	)
	# (Possibly) apply multiprocessing
	all_inputs = list(range(seed_start, seed_start+reps))
	num_processes = min(len(all_inputs), num_processes)
	all_outputs = apply_pool(
		func = partial_sd_power_fdr,
		all_inputs = all_inputs,
		num_processes=num_processes
	)

	# Extract output and return
	final_out = []
	for output in all_outputs:
		final_out.extend(output)

	# Final out: list of arrays: 
	# S_method, power, fdp, q, score, score_type, pair_agg, W, Z, tildeZ, selection
	return final_out

def analyze_degen_solns(
	Sigma,
	beta,
	dgp_number,
	groups=None,
	sample_kwargs={},
	filter_kwargs={},
	fstat_kwargs={},
	reps=50,
	num_processes=5,
	seed_start=0,
	time0=None,
	S_curve=False,
	storew=True,
	):
	"""
	:param dgp_number: A number corresponding to which dgp
	each row corresponds to.
	:param sample_kwargs: 
	A dictionary. Each key is a sample parameter name, with the value
	being either a single value or a list or array of values. 
	:param filter_kwargs: 
	A dictionary. Each key is a mx_filter parameter name, with the value
	being either a single value or a list or array of values. 
	"""
	global MAXIMUM_N # A hack to allow for better logging
	
	# Infer p and set n defaults
	if Sigma is not None:
		p = Sigma.shape[0]
	else:
		p = sample_kwargs['p']
		if isinstance(p, list):
			p = p[0]
	sample_kwargs['p'] = [p]
	if 'n' not in sample_kwargs:
		sample_kwargs['n'] = [
			int(p/4), int(p/2), int(p/1.5), int(p), int(2*p), int(4*p)
		]
	if not isinstance(sample_kwargs['n'], list):
		sample_kwargs['n'] = [sample_kwargs['n']]
	MAXIMUM_N = max(sample_kwargs['n']) # Helpful for logging
	if groups is None:
		groups = np.arange(1, p+1, 1)

	# Construct/iterate cartesian product of sample, filter, fstat kwargs
	sample_keys, sample_product = dict2keyproduct(sample_kwargs)
	filter_keys, filter_product = dict2keyproduct(filter_kwargs)
	fstat_keys, fstat_product = dict2keyproduct(fstat_kwargs)

	# Initialize final output
	counter = 0
	columns = ['power', 'fdp', 'q', 'S_method', 'antisym', 'score', 'score_type']
	if storew:
		columns += [f'W{i}' for i in range(1, p+1)]
		columns += [f'Z{i}' for i in range(1, p+1)]
		columns += [f'tildeZ{i}' for i in range(1, p+1)]
		columns += [f'selection{i}' for i in range(1, p+1)] 
	columns += ['dgp_number']
	columns += [x for x in sample_keys if x != 'gibbs_graph']
	columns += filter_keys + fstat_keys + ['seed']
	result_df = pd.DataFrame(columns=columns)
	
	# Check if we are going to ever fit MX knockoffs on the
	# "ground truth" covariance matrix. If so, we'll memoize
	# the SDP/MVR/MMI results.
	if 'fixedx' in filter_kwargs:
		fixedX_vals = filter_kwargs['fixedx']
		if False in fixedX_vals:
			MX_flag = True
		else:
			MX_flag = False
	else:
		MX_flag = True
	if 'infer_sigma' in filter_kwargs:
		infer_sigma_vals = filter_kwargs['infer_sigma']
		if False in infer_sigma_vals:
			ground_truth = True
		else:
			ground_truth = False
	else:
		ground_truth = True
	# Possibly save S-matrices
	if ground_truth and MX_flag and Sigma is not None:
		rej_rate = fetch_kwarg(filter_kwargs, 'rej_rate', default=[0])[0]
		print(f"Storing SDP/MVR/MMI results with rej_rate={rej_rate}")
		S_matrices = fetch_competitor_S(
			Sigma=Sigma,
			groups=groups,
			time0=time0,
			S_curve=S_curve,
			rej_rate=rej_rate,
			verbose=True,
		)
	else:
		print(f"Not storing SDP/MVR results")
		S_matrices = None

	### Calculate power of knockoffs for the different methods
	for filter_vals in filter_product:
		filter_vals = list(filter_vals)
		new_filter_kwargs = {
			key:val for key, val in zip(filter_keys, filter_vals)
		}

		# In high dimensional cases or binomial cases,
		# don't fit OLS.
		if 'feature_stat' in new_filter_kwargs:
			fstat = new_filter_kwargs['feature_stat']
		else:
			fstat = 'lasso'

		for sample_vals in sample_product:
			sample_vals = list(sample_vals)
			new_sample_kwargs = {
				key:val for key, val in zip(sample_keys, sample_vals)
			}

			# Don't run OLS in certain cases
			if fstat == 'ols':
				if new_sample_kwargs['n'] < 2*new_sample_kwargs['p']:
					continue
				if 'y_dist' in new_sample_kwargs:
					if new_sample_kwargs['y_dist'] == 'binomial':
						continue
			if fstat == 'dlasso':
				if 'y_dist' in new_sample_kwargs:
					if new_sample_kwargs['y_dist'] == 'binomial':
						continue

			for fstat_vals in fstat_product:
				# Extract feature-statistic kwargs
				# and place them properly (as a dictionary value
				# of the filter_kwargs)
				fstat_vals = list(fstat_vals)
				new_fstat_kwargs = {
					key:val for key, val in zip(fstat_keys, fstat_vals)
				}
				if 'feature_stat_kwargs' in new_filter_kwargs:
					new_filter_kwargs['feature_stat_kwargs'] = dict(
						new_filter_kwargs['feature_stat_kwargs'],
						**new_fstat_kwargs
					)
				else:
					new_filter_kwargs['feature_stat_kwargs'] = new_fstat_kwargs

				# Power/FDP for the S methods
				out = calc_power_and_fdr(
					Sigma=Sigma,
					beta=beta,
					groups=groups,
					reps=reps,
					num_processes=num_processes,
					sample_kwargs=new_sample_kwargs,
					filter_kwargs=new_filter_kwargs,
					seed_start=seed_start,
					time0=time0,
					S_matrices=S_matrices,
					S_curve=S_curve,
				)

				# Loop through antisymmetric functions and S matrices
				for vals in out:
					S_method, power, fdp, q, score, score_type, agg, W, Z, tildeZ, selections, seed = vals
					row = [power, fdp, q, S_method, agg, score, score_type] 
					if storew:
						row.extend(W.tolist())
						row.extend(Z.tolist())
						row.extend(tildeZ.tolist())
						row.extend(selections.astype(np.int32).tolist())
					row += [dgp_number]
					sample_vals_to_add = [
						sample_vals[i] for i in range(len(sample_vals)) if sample_keys[i] != 'gibbs_graph'
					] 
					row += sample_vals_to_add + filter_vals + fstat_vals + [seed]
					result_df.loc[counter] = row 
					counter += 1

	return result_df, S_matrices

def parse_args(args):
	""" Homemade argument parser 
	Usage: --argname_{dgp/sample/filter} value
	Value should be one of:
		- integer/float
		- arbitrary string
		- Alternatively, string of the form
		"start{num}end{num}numvals{num}" 
		which indicates that the parameter ought 
		to be varied along a linear interpolation
		from the specified start/end with the specified
		number of values.
		- Alternatively, a string like
		"[str1, str2, str3]" which will be parsed as
		["str1", "str2", "str3"]
	
	Note that "reps" and "num_processes" are special kwargs.
	They get placed as sample_kwargs by default but will be
	extracted.

	The --description argument must be the last argument, as 
	all arguments after the --description argument will be ignored.
	"""
	
	# Get rid of script name
	args = args[1:]

	# Initialize kwargs constructor
	special_keys = [
		'reps', 'num_processes', 'seed_start', 'description',
		's_curve', 'resample_beta', 'resample_sigma', 'storew',
	]
	key_types = ['dgp', 'sample', 'filter', 'fstat']
	all_kwargs = {ktype:{} for ktype in key_types}
	key = None
	key_type = None # One of dgp, sample, filter
	description_index = None # At the end, can write description

	#Parse
	for i, arg in enumerate(args):
		arg = str(arg).lower()

		# Even placements should be keyword
		if i % 2 == 0:
			if str(arg)[0:2] != '--':
				raise ValueError(
					f'{i}th arg ({arg}) should be a kwarg name preceded by "--"'
				)
			key = arg[2:]

			# Check what type of keyword
			if key in special_keys:
				key_type = 'sample'
			else:
				key_split = key.split('_')
				key_type = key_split[-1]
				key = ('_').join(key_split[0:-1])

			# Raise error if unexpected key type
			if key_type not in key_types:
				raise ValueError(
					f'{i}th arg ({arg}) has key_type {key_type}, must be one of {key_types}'
				)

			# Friendly reminder
			if key == 'feature_stat_fn':
				raise ValueError("feature_stat_fn is depreciated, use feature_stat")

			# Description
			if key == 'description':
				description_index = i
				break

		# Parse values
		if i % 2 == 1:

			# Check if it's a list written out
			if arg[0] == '[':
				# Strip brackets, whitspace, and split by commas
				value = arg.replace('[', '').replace(']', '')
				value = value.replace(' ', '').split(',')
				# Try to process
				value = [obj2int(v) for v in value]

			# Check if it's start{num}end{num}numvals{num}
			elif arg[0:5] == 'start':
				start = float(arg.replace('start', '').split('end')[0])
				end = float(arg.split('end')[-1].split('numvals')[0])
				numvals = int(arg.split('numvals')[-1])
				value = np.linspace(
					start, end, numvals
				)

			# Apply obj2int (preserves strings, infers floats, bools, ints)
			else:
				value = obj2int(arg)


			all_kwargs[key_type][key] = val2list(value)

	# Parse description 
	description = ''
	if description_index is not None:
		description += (' ').join(args[description_index+1:])
	description += f'\n \n Other arguments were: {args[0:description_index]}'
	all_kwargs['sample']['description'] = description

	out = (all_kwargs[key_type] for key_type in key_types)
	return out

def main(args):
	""" Layers of keyword arguments are as follows.
	- dgp_kwargs: kwargs needed to create Sigma (e.g. p, method)
	- sample_kwargs: kwargs needed to sample data (e.g. n)
	- filter_kwargs: kwargs for the knockoff filter. (e.g. recycling)
	- fstat_kwargs: kwargs for the feature statistic. (e.g. pair_agg)
	
	For each of these, keys mapping to iterables like lists or 
	numpy arrays will be varied. See parse_args comments.

	The MAIN constraint is that the same sample_kwargs must be used
	for all dgp_kwargs. E.g., it will be difficult to vary both
	a for an AR1 covariance matrix and delta for an ErdosRenyi covariance
	matrix. The same goes for filter_kwargs abd fstat_kwargs.
	"""

	# Create kwargs
	dgp_kwargs, sample_kwargs, filter_kwargs, fstat_kwargs = parse_args(args)
	print(f"Args were {args}")

	# Make sure pair_aggs is not being duplicated
	if 'pair_agg' in fstat_kwargs:
		raise ValueError("Many pair_aggs will be analyzed anyway. Do not add this as a fstat_kwarg.")
	# Some common errors I make
	if 'y_dist' in dgp_kwargs:
		raise ValueError("y_dist ought to be in sample_kwargs")

	# Parse some special non-graph kwargs
	reps = fetch_kwarg(sample_kwargs, 'reps', default=[1])[0]
	num_processes = fetch_kwarg(sample_kwargs, 'num_processes', default=[12])[0]
	seed_start = fetch_kwarg(sample_kwargs, 'seed_start', default=[0])[0]
	description = fetch_kwarg(sample_kwargs, 'description', default='')
	S_curve = fetch_kwarg(sample_kwargs, 's_curve', default=[False])[0]
	resample_beta = fetch_kwarg(sample_kwargs, 'resample_beta', default=[True])[0]
	resample_sigma = fetch_kwarg(sample_kwargs, 'resample_sigma', default=[True])[0]
	storew = fetch_kwarg(sample_kwargs, 'storew', default=[False])[0]
	print(f"S_curve flag = {S_curve}")
	print(f"DGP kwargs are {dgp_kwargs}")
	print(f"Sample kwargs are {sample_kwargs}")
	print(f"Filter kwargs are {filter_kwargs}")
	print(f"Ftat kwargs are {fstat_kwargs}")
	print(f"Description is {description.split('Other arguments were:')[0]}")


	# Create output paths
	today = str(datetime.date.today())
	hour = str(datetime.datetime.today().time())
	hour = hour.replace(':','-').split('.')[0]
	output_path = f'data/mrcfinal/{today}/{hour}'
	all_key_types = ['dgp', 'sample', 'filter', 'fstat']
	all_kwargs = [dgp_kwargs,sample_kwargs, filter_kwargs, fstat_kwargs]

	for key_type,kwargs in zip(all_key_types, all_kwargs):
		output_path += f'/{key_type}_'
		keys = sorted([k for k in kwargs])
		for k in keys:
			path_val = ('').join(str(kwargs[k]).split(' '))
			output_path += f'{k}{path_val}_'

	# Put it all together and ensure directory exists
	output_path += f'seedstart{seed_start}_reps{reps}_results.csv'
	beta_path = output_path.split('.csv')[0] + '_betas.csv'
	S_path = output_path.split('.csv')[0] + '_S.csv'
	description_path = f'data/mrcfinal/{today}/{hour}/' + 'description.txt'
	output_dir = os.path.dirname(output_path)
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	print(f"Output path is {output_path}")

	# Save description
	with open(description_path, 'w') as thefile:
		thefile.write(description)	


	# Initialize final final output
	all_results = pd.DataFrame()

	# Loop through DGP parameters
	dgp_keys = sorted([key for key in dgp_kwargs])
	dgp_components = []
	for key in dgp_keys:
		dgp_kwargs[key] = val2list(dgp_kwargs[key])
		dgp_components.append(dgp_kwargs[key])
	dgp_product = itertools.product(*dgp_components)

	# p sadly does need to have 1 value for a variety of reasons 
	if 'p' in dgp_kwargs:
		if len(dgp_kwargs['p']) > 1:
			raise ValueError(f"Must have only one value of p, not {dgp_kwarys[p]}")
		else:
			p = dgp_kwargs['p'][0]
	else:
		p = 100

	# Initialize ways to save beta, dgp
	dgp_number = 0
	beta_df = pd.DataFrame(columns=np.arange(1, p+1, 1))
	beta_df.index.name = 'dgp_number'
	# Stores diagonal of S
	S_diags_df = pd.DataFrame(
		columns = ['dgp_number', 'S_method'] + [i for i in range(1, p+1)]
	) 
	S_counter = 0


	for dgp_vals in dgp_product:

		# Create DGP using dgp kwargs
		dgp_vals = list(dgp_vals)
		new_dgp_kwargs = {key:val for key,val in zip(dgp_keys, dgp_vals)}
		print(f"DGP kwargs are now: {new_dgp_kwargs}")
		np.random.seed(110)
		_, _, beta, invSigma, Sigma = knockadapt.graphs.sample_data(
			**new_dgp_kwargs
		)

		# Use precomputed Sigma for ising model
		if 'x_dist' in new_dgp_kwargs:
			if new_dgp_kwargs['x_dist'] == 'gibbs':
				if 'method' not in new_dgp_kwargs:
					raise ValueError(f"Method must be supplied for x_dist == gibbs")
				if new_dgp_kwargs['method'] == 'ising':
					p = new_dgp_kwargs['p']
					file_dir = os.path.dirname(os.path.abspath(__file__))
					v_file = f'{file_dir}/qcache/vout{p}.txt'
					if os.path.exists(v_file):
						print(f"Loading custom Sigma for gibbs ising model")
						Sigma = np.loadtxt(f'{file_dir}/qcache/vout{p}.txt')
					else:
						print(f"No custom Sigma available for gibbs ising model: using default")
			# Save gibbs_graph
			sample_kwargs['gibbs_graph'] = [invSigma]

		# Cache beta
		if not resample_beta:
			beta_df.loc[dgp_number] = beta
			beta_df.to_csv(beta_path)

		# Push sample_kwargs to full output
		if resample_beta or resample_sigma:
			for key in new_dgp_kwargs:
				val = new_dgp_kwargs[key]
				if key in sample_kwargs:
					old_val = sample_kwargs[key]
					if len(old_val) != 1:
						raise ValueError(f"DGP / sample keys {key} conflict when resampling Sigma / beta")
					old_val = old_val[0]
					if old_val not in dgp_kwargs[key]:
						raise ValueError(f"DGP / sample keys {key} conflict when resampling Sigma / beta")
				sample_kwargs[key] = [val]

		# Create results
		result, S_matrices = analyze_degen_solns(
			Sigma=Sigma if not resample_sigma else None,
			beta = beta if not resample_beta else None,
			dgp_number=dgp_number,
			groups=None,
			sample_kwargs=sample_kwargs,
			filter_kwargs=filter_kwargs,
			fstat_kwargs=fstat_kwargs,
			reps=reps,
			num_processes=num_processes,
			seed_start=seed_start,
			time0=time0,
			S_curve=S_curve,
			storew=storew,
		)
		for key in dgp_keys:
			if key in sample_kwargs:
				continue
			result[key] = new_dgp_kwargs[key]
		all_results = all_results.append(
			result, 
			ignore_index = True
		)
		all_results.to_csv(output_path)

		# Cache S outputs
		if S_matrices is not None:
			for S_method in S_matrices:
				S = S_matrices[S_method]
				if S is None:
					continue
				S_diag = np.diag(S)
				S_diags_df.loc[S_counter] = [dgp_number, S_method] + S_diag.tolist()
				S_counter += 1
			S_diags_df.to_csv(S_path)

		# Increment dgp number
		dgp_number += 1

	return all_results

if __name__ == '__main__':

	time0 = time.time()
	MAXIMUM_N = 0
	main(sys.argv)


