import os
import sys
import time
import datetime

# Ensure we are using the right version of knockpy
import knockpy
print(f"Using version {knockpy.__version__} of knockpy")
from knockpy import knockoff_stats as kstats
from knockpy import utilities
from knockpy.knockoff_filter import KnockoffFilter

import warnings
import numpy as np
import pandas as pd

import itertools
from multiprocessing import Pool
from functools import partial

# Global: the set of antisymmetric functions we use
antisyms = ['cd', 'sm']
# Levels of FDR control values to evaluate for
q_values = 0.05*np.arange(1, 21)
# Degrees of freedom for t-distributions
DEFAULT_DF_T = knockpy.dgp.DEFAULT_DF_T

# Defaults 
DEFAULT_FSTAT = 'lasso'
DEFAULT_KSAMPLER = 'gaussian'
DEFAULT_SHRINKAGE = 'ledoitwolf'
SHRINKAGE_TOL = 1e-4
S_METHODS = ['mvr', 'maxent', 'sdp']
DEFAULT_INFER_SIGMA = False
SDPTOL = 0.01

# Defaults for printing output
DEFAULT_ANTISYM_4PRINT = 'cd'
DEFAULT_Q_4PRINT = 0.1

### Utility functions which improve parsing
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

def check_equi_degen(Sigma):
    """
    Checks if Sigma is equicorr. with rho \ge 0.5
    """
    if Sigma is None:
        return False
    p = Sigma.shape[0]
    rho = Sigma[0, 1]
    if rho >= 0:
        equicorr = rho*np.ones((p, p)) + (1-rho)*np.eye(p)
        if np.all(np.abs(Sigma-equicorr) < 1e-10):
            if rho >= 0.5:
                return True
    return False


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
    kwargs.pop('sdp_tol', None) # Prevent double arg errors
    if S_curve:
        mineig = np.linalg.eigh(Sigma)[0].min()
        gammas = 2*mineig*np.arange(0, 11, 1)/10
        gammas[0] += 0.001
        gammas[-1] -= 0.001
        S_matrices = {
            f'S{np.around(gamma, 3)}':gamma*np.eye(p) for gamma in gammas
        }
        for S_method in S_METHODS:
            S_matrices[S_method] = knockpy.smatrix.compute_smatrix(
                Sigma=Sigma,
                groups=groups,
                method=S_method,
                solver='cd',
                tol=1e-10,
            )
        # Create perturbed option when G_SDP is low rank
        if Sigma[0,1] >= 0.5:
            S_matrices['sdp_perturbed'] = (1-SDPTOL)*S_matrices['sdp']
        return S_matrices

    ### Special case: detect if Sigma is equicorrelated,
    # in which case we can calculate the solution analytically
    # for ungrouped knockoffs.
    if Sigma is not None:
        if np.unique(groups).shape[0] == p:
            if check_equi_degen(Sigma):
                rho = Sigma[0,1]
                print(f"Sigma is equicorr (rho={rho}), using analytical solution")
                S_SDP = min(1, 2-2*rho)*np.eye(p)
                S_matrices = {'sdp':S_SDP}
                # Compute exact solution for smaller p
                if p < 1000:
                    for method in S_METHODS:
                        if method == 'sdp':
                            continue
                        S_matrices[method] = knockpy.smatrix.compute_smatrix(
                            Sigma=Sigma,
                            groups=groups,
                            method=method,
                            solver='cd',
                        )

                    # Else, asymptotically these methods are the same
                    else:
                        S_matrices['mvr'] = (1-rho)*np.eye(p)
                    # For high rho, check non-degenerate version of SDP
                    if rho < 0.5:
                        return S_matrices
                    else:
                        S_matrices['sdp_perturbed'] = S_SDP*(0.99)
                        return S_matrices

    ### Generically calculate SDP/MCR matrices
    S_matrices = {}
    for new_method in S_METHODS:
        S_MRC = knockpy.smatrix.compute_smatrix(
            Sigma=Sigma,
            groups=groups,
            method=new_method,
            solver='cd',
            **kwargs
        )
        S_matrices[new_method] = S_MRC

    if verbose:
        print(f'Finished computing MRC matrices, time is {time.time() - time0}')

    return S_matrices

class ExactSDPEquiKnock(knockpy.knockoffs.GaussianSampler):
    """
    Samples SDP knockoffs exactly, with no
    tolerance added (so no numerical errors)
    for equicorrelated covariance matrices with
    rho >= 0.5.
    """

    def sample_knockoffs(self):
        """ 
        ensures GSDP is exactly rank p+1 with the right low-rank structure
        """

        super().sample_knockoffs()
        sumcols = self.X[:, 0] + self.Xk[:, 0]
        self.Xk = sumcols.reshape(-1, 1) - self.X
        return self.Xk

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
    cutoff,
    groups,
    normalize=True,
    sample_kwargs={},
    filter_kwargs={
        'fstat_kwargs':{},
        'knockoff_kwargs':{},
    },
    S_matrices=None,
    time0=None,
    S_curve=False,
):
    """ 
    Knockoff kwargs should be included in filter_kwargs
    """

    # Prevent global changes to kwargs
    filter_kwargs = filter_kwargs.copy()
    sample_kwargs = sample_kwargs.copy()

    ##### Step 1: Sample the data
    if Sigma is not None:
        p = Sigma.shape[0]
    else:
        p = sample_kwargs['p']

    # If DGP parameters are previously generated, 
    # add them in here. The gibbs_graph parameter
    # will be `None` unless we are doing gibbs sampling.
    np.random.seed(seed)
    gibbs_graph = sample_kwargs.pop('gibbs_graph', None) 
    dgprocess = knockpy.dgp.DGP(
        Sigma=Sigma,
        beta=beta,
        gibbs_graph=gibbs_graph
    )
    dgprocess.sample_data(
        **sample_kwargs
    )
    X = dgprocess.X
    n = X.shape[0]
    y = dgprocess.y
    beta = dgprocess.beta
    Sigma = dgprocess.Sigma

    # Create groups, almost always trivial groupsin our experiments
    if groups is None:
        if cutoff==1:
            groups = np.arange(1, p + 1)
        else:
            groups = knockpy.dgp.create_grouping(Sigma, cutoff=cutoff)

    # Parse nonnulls 
    group_nonnulls = knockpy.utilities.fetch_group_nonnulls(
        beta, groups
    )

    ### Step 2: Misc tasks before we loop through S matrices. 
    # Step 2A: we extract keyword arguments.
    ksampler = filter_kwargs.pop('ksampler', DEFAULT_KSAMPLER)
    fstat = filter_kwargs.pop('fstat', DEFAULT_FSTAT)
    # This means we do not pass the (true) Sigma to the filter
    infer_sigma = filter_kwargs.pop('infer_sigma', DEFAULT_INFER_SIGMA)

    # Step 2B: check if we're in the equicorrelated case ---
    # we do a a few extra things: (i) we compute exact SDP knockoffs,
    # (ii) we compute a perturbed version of SDP knockoffs which is not
    # exact.
    equi_degen_flag = check_equi_degen(dgprocess.Sigma)
    if equi_degen_flag and not infer_sigma and ksampler == 'gaussian':
        rho = dgprocess.Sigma[0, 1]
        S_matrices['sdp_perturbed'] = (1-SDPTOL) * (2 - 2*rho) * np.eye(p)

    # Step 2C: if we are inferring Sigma, do this now to share computation
    # for all S-matrix generation.
    if infer_sigma:
        shrinkage = filter_kwargs.pop('shrinkage', DEFAULT_SHRINKAGE)
        if shrinkage == 'none' and n < p: # Shrinkage is required in the high-dim case
            return []
        inferred_Sigma, _ = knockpy.utilities.estimate_covariance(X, tol=SHRINKAGE_TOL, shrinkage=shrinkage)

    ### Step 3: Now we loop through the S matrices
    output = []
    for S_method in S_matrices:
        n = X.shape[0]
        print(f"Now computing for S_method={S_method} with n={n} seed={seed} fstat={fstat} ksampler={ksampler}...")

        # Pull S matrix
        S = S_matrices[S_method]

        # Create knockoff_kwargs
        if 'knockoff_kwargs' not in filter_kwargs:
            filter_kwargs['knockoff_kwargs'] = {}
        filter_kwargs['knockoff_kwargs']['S'] = S # possibly None
        filter_kwargs['knockoff_kwargs']['method'] = S_method

        # Possibly pass a few parameters to metro sampler
        if 'x_dist' in sample_kwargs:
            # For gibbs graph
            if gibbs_graph is not None and ksampler == 'gibbs_grid':
                filter_kwargs['knockoff_kwargs']['gibbs_graph'] = gibbs_graph
            # For t-distributions
            if str(sample_kwargs['x_dist']).lower() in ['ar1t', 'blockt']:
                if 'df_t' in sample_kwargs and ksampler != 'gaussian':
                    filter_kwargs['knockoff_kwargs']['df_t'] = sample_kwargs['df_t']
                elif ksampler != 'gaussian':
                    filter_kwargs['knockoff_kwargs']['df_t'] = DEFAULT_DF_T # This matters
                else:
                    filter_kwargs['knockoff_kwargs'].pop('df_t', None)
            if ksampler != 'gaussian' and ksampler != 'gibbs_grid':
                filter_kwargs['knockoff_kwargs']['tol'] = min(1e-2, np.linalg.eigh(Sigma)[0].min()/10)

        # In the equicorrelated case
        # we overwride the KnockoffFilter's attempt to prevent
        # numerical errors and sample the knockoffs ourselves.
        # This flag tells us when to do that.
        # We don't want to run OLS/debiased lasso in this case,
        # since everything is low rank.
        _sdp_degen = (equi_degen_flag and S_method == 'sdp' and ksampler=='gaussian' and not infer_sigma)
        if _sdp_degen:
            if filter_kwargs.get('fstat', 'lasso') in ['dlasso', 'ols']:
                continue
            print(f"Using exact SDP knockoffs (no numerical errors) for equi cov matrix")
            sdp_ksampler = ExactSDPEquiKnock(
                X=X, Sigma=Sigma, **filter_kwargs['knockoff_kwargs']
            )

        # Initialize knockoff filter
        kfilter = KnockoffFilter(
            ksampler=ksampler if not _sdp_degen else sdp_ksampler,
            fstat=fstat,
        )
        kfilter.forward(
            X=X, 
            y=y, 
            mu=np.zeros(p), # common to all experiments
            Sigma=inferred_Sigma if infer_sigma else Sigma, # ignored by fx 
            groups=groups,
            **filter_kwargs
        )
        Z = kfilter.Z
        score = kfilter.score
        score_type = kfilter.score_type

        # Calculate power/fdp/score for a variety of 
        # antisymmetric functions
        for antisym in antisyms:
            for q in q_values:
                # Start by creating selections
                selections, W = Z2selections(
                    Z=Z,
                    groups=groups,
                    q=q,
                    antisym=antisym
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
                    antisym,
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
    # [S_method, power, fdp, q, score, score_type, antisym, seed]
    return output

def calc_power_and_fdr(
    time0,
    Sigma,
    beta,
    cutoff=1,
    groups=None,
    q=0.2,
    reps=100,
    num_processes=1,
    sample_kwargs={},
    filter_kwargs={},
    seed_start=0,
    S_matrices=None,
    S_curve=False,
):

    # Sample data reps times and calc power/fdp
    num_processes = min(reps, num_processes)
    all_outputs = knockpy.utilities.apply_pool(
        func=single_dataset_power_fdr,
        seed=list(range(seed_start, seed_start+reps)),
        num_processes=num_processes,
        constant_inputs={
            'Sigma':Sigma,
            'beta':beta,
            'cutoff':cutoff,
            'groups':groups,
            'sample_kwargs':sample_kwargs,
            'filter_kwargs':filter_kwargs,
            'S_matrices':S_matrices,
            'time0':time0,
            'S_curve':S_curve,
        }
    )

    # Rormat slightly
    final_out = []
    for output in all_outputs:
        final_out.extend(output)

    # Final out: list of arrays: 
    # S_method, power, fdp, q, score, score_type, antisym
    return final_out

def compare_S_methods(
    Sigma,
    beta,
    dgp_number,
    cutoff=1,
    sample_kwargs={},
    filter_kwargs={},
    fstat_kwargs={},
    reps=50,
    num_processes=5,
    seed_start=0,
    time0=None,
    S_curve=False,
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
    
    # Infer p, set defaults for n
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

    # Helpful for logging
    MAXIMUM_N = max(sample_kwargs['n']) 

    # Construct/iterate cartesian product of sample, filter, fstat kwargs
    sample_keys, sample_product = dict2keyproduct(sample_kwargs)
    filter_keys, filter_product = dict2keyproduct(filter_kwargs)
    fstat_keys, fstat_product = dict2keyproduct(fstat_kwargs)

    # Initialize final output
    counter = 0
    columns = ['power', 'fdp', 'q', 'S_method', 'antisym', 'score', 'score_type']
    columns += ['dgp_number']
    columns += [x for x in sample_keys if x != 'gibbs_graph']
    columns += filter_keys + fstat_keys + ['seed']
    result_df = pd.DataFrame(columns=columns)
    
    # Check if we are going to ever fit MX knockoffs on the
    # "ground truth" covariance matrix. If so, we'll memoize
    # the SDP/MRC results.
    ksampler_vals = filter_kwargs.get('ksampler', ['gaussian'])
    if ksampler_vals == ['fx']:
        MX_flag = False
    else:
        MX_flag = True
    infer_sigma_vals = filter_kwargs.get(
        'infer_sigma', [DEFAULT_INFER_SIGMA]
    )
    if False in infer_sigma_vals:
        ground_truth = True
    else:
        ground_truth = False

    # Possibly save S-matrices
    if ground_truth and MX_flag and Sigma is not None:
        print(f"Storing SDP/MRC results")
        if cutoff != 1:
            groups = knockpy.dgp.create_grouping(Sigma, cutoff=cutoff)
        else:
            groups = np.arange(1, p + 1)
        S_matrices = fetch_competitor_S(
            Sigma=Sigma,
            groups=groups,
            time0=time0,
            S_curve=S_curve,
            verbose=True,
        )
    else:
        print(f"Not storing SDP/MVR results")
        S_matrices = {method:None for method in S_METHODS}
        groups = None

    ### Calculate power of knockoffs for the different methods
    for filter_vals in filter_product:
        filter_vals = list(filter_vals)
        new_filter_kwargs = {
            key:val for key, val in zip(filter_keys, filter_vals)
        }

        # In high dimensional cases or binomial cases,
        # don't fit OLS.
        if 'fstat' in new_filter_kwargs:
            fstat = new_filter_kwargs['fstat']
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
                if 'fstat_kwargs' in new_filter_kwargs:
                    new_filter_kwargs['fstat_kwargs'] = dict(
                        new_filter_kwargs['fstat_kwargs'],
                        **new_fstat_kwargs
                    )
                else:
                    new_filter_kwargs['fstat_kwargs'] = new_fstat_kwargs

                # Power/FDP for the S methods
                out = calc_power_and_fdr(
                    Sigma=Sigma,
                    beta=beta,
                    cutoff=cutoff,
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
                    S_method, power, fdp, q, score, score_type, agg, seed = vals
                    row = [power, fdp, q, S_method, agg, score, score_type] 
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
    Usage: --argname_{dgp/sample/filter/fstat} value
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
        's_curve', 'resample_beta', 'resample_sigma', 'comp_ci',
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
            if key == 'feature_stat':
                raise ValueError("fstat_fn is depreciated, use fstat")

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
        description += '\n'
    description = 'Arguments were: \n'
    i = 0
    arg_msg = ''
    max_ind = len(args) if description_index is None else description_index
    while i < max_ind:
        description += f'{args[i]}={args[i+1]}\n'
        i += 2
    all_kwargs['sample']['description'] = description

    out = (all_kwargs[key_type] for key_type in key_types)
    return out

def print_results(results):
    """
    Helper function which prints power/FDR grouped by relevant variables
    """
    # Subset
    results = results.loc[
        (results['antisym'] == DEFAULT_ANTISYM_4PRINT) & 
        (results['q'] == DEFAULT_Q_4PRINT)
    ]

    # Find active variables, meaning params which vary
    MEASUREMENTS = ['power', 'fdp', 'score', 'mac']
    EXCLUDED = ['seed', 'dgp_number', 'p'] 
    ACTIVE_VARIABLES = []
    for col in results.columns:
        if col in MEASUREMENTS or col in EXCLUDED:
            continue
        try:
            if np.unique(results[col]).shape[0] == 1:
                continue
        except TypeError:
            results[col] = results[col].astype(str)
            if np.unique(results[col].shape[0]) == 1:
                continue
        if 'score' not in col:
            ACTIVE_VARIABLES.append(col)

    # Take mean, std
    meas_subset = [m for m in MEASUREMENTS if m in results.columns]
    num_samples = np.unique(results['seed']).shape[0]
    mean_meas = results.groupby(
        ACTIVE_VARIABLES
    )[meas_subset].mean().reset_index()
    std_meas = (results.groupby(
        ACTIVE_VARIABLES
    )[meas_subset].std()/np.sqrt(num_samples)).reset_index()

    # Merge measurements together
    print(ACTIVE_VARIABLES)
    all_meas = pd.merge(
        mean_meas,
        std_meas, 
        on=ACTIVE_VARIABLES,
        suffixes = ('_mean', '_se')
    )
    all_meas = all_meas.drop(['score_mean', 'score_se'], 'columns')

    # Add n/p information
    all_meas['p'] = results['p'].unique().tolist()[0]
    if 'n' not in all_meas.columns:
        assert(results['n'].unique().shape[0] == 1)
        all_meas['n'] = results['n'].unique().tolist()[0]

    # print!
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 30)
    print(all_meas)
    pd.reset_option('display.max_rows|display.max_columns|display.width')



def main(args):
    """ Layers of keyword arguments are as follows.
    - dgp_kwargs: kwargs needed to create Sigma (e.g. p, method)
    - sample_kwargs: kwargs needed to sample data (e.g. n)
    - filter_kwargs: kwargs for the knockoff filter. (e.g. recycling)
    - fstat_kwargs: kwargs for the feature statistic. (e.g. antisym)
    
    For each of these, keys mapping to iterables like lists or 
    numpy arrays will be varied. See parse_args comments.

    The MAIN constraint is that the same sample_kwargs must be used
    for all dgp_kwargs. E.g., it will be difficult to vary both
    a for an AR1 covariance matrix and delta for an ErdosRenyi covariance
    matrix. The same goes for filter_kwargs abd fstat_kwargs.
    """
    global S_METHODS

    # Create kwargs
    dgp_kwargs, sample_kwargs, filter_kwargs, fstat_kwargs = parse_args(args)
    print(f"Args were {args}")

    # Make sure antisyms is not being duplicated
    if 'antisym' in fstat_kwargs:
        raise ValueError("Many antisyms will be analyzed anyway. Do not add this as a fstat_kwarg.")
    # Some common errors I make
    if 'y_dist' in dgp_kwargs:
        raise ValueError("y_dist ought to be in sample_kwargs")

    # Parse some special kwargs
    reps = sample_kwargs.pop('reps', [1])[0]
    num_processes = sample_kwargs.pop('num_processes', [12])[0]
    seed_start = sample_kwargs.pop('seed_start', [0])[0]
    description = sample_kwargs.pop('description', '')
    S_curve = sample_kwargs.pop('s_curve', [False])[0]
    resample_beta = sample_kwargs.pop('resample_beta', [True])[0]
    resample_sigma = sample_kwargs.pop('resample_sigma', [True])[0]
    comp_ci = sample_kwargs.pop('comp_ci', [False])[0]
    if comp_ci:
        S_METHODS.append('ciknock')

    # Log what's going on
    print(f"S_curve flag = {S_curve}")
    print(f"DGP kwargs are {dgp_kwargs}")
    print(f"Sample kwargs are {sample_kwargs}")
    print(f"Filter kwargs are {filter_kwargs}")
    print(f"Ftat kwargs are {fstat_kwargs}")
    print(f"Description is {description.split('Arguments were:')[0]}")

    # Create output paths
    today = str(datetime.date.today())
    hour = str(datetime.datetime.today().time())
    hour = hour.replace(':','-').split('.')[0]
    output_path = f'data/{today}/{hour}/'

    # Put it all together and ensure directory exists
    output_path += f'seedstart{seed_start}_reps{reps}_results.csv'
    description_path = f'data/{today}/{hour}/' + 'description.txt'
    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    print(f"Output path is {output_path}")

    # Save description
    with open(description_path, 'w') as thefile:
        thefile.write(description)  

    # Initialize final final output
    all_results = pd.DataFrame()

    # p sadly does need to have 1 value for a variety of reasons 
    if 'p' in dgp_kwargs:
        if len(dgp_kwargs['p']) > 1:
            raise ValueError(f"Must have only one value of p, not {dgp_kwargs[p]}")
        else:
            p = dgp_kwargs['p'][0]
    else:
        p = 100
        dgp_kwargs['p'] = [p]

    # Loop through DGP parameters
    dgp_keys = sorted([key for key in dgp_kwargs])
    dgp_components = []
    for key in dgp_keys:
        dgp_kwargs[key] = val2list(dgp_kwargs[key])
        dgp_components.append(dgp_kwargs[key])
    dgp_product = itertools.product(*dgp_components)

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

        # Create DGP using dgp kwargs. These will be ignored
        # if resample_sigma and resample_beta is True
        dgp_vals = list(dgp_vals)
        new_dgp_kwargs = {key:val for key,val in zip(dgp_keys, dgp_vals)}
        print(f"DGP kwargs are now: {new_dgp_kwargs}")
        cutoff = new_dgp_kwargs.pop("cutoff", 1)
        np.random.seed(110)
        dgprocess = knockpy.dgp.DGP()
        _, _, beta, invSigma, Sigma = dgprocess.sample_data(
            **new_dgp_kwargs
        )

        # Use precomputed Sigma for ising model
        if 'x_dist' in new_dgp_kwargs:
            if new_dgp_kwargs['x_dist'] == 'gibbs':
                p = new_dgp_kwargs['p']
                file_dir = os.path.dirname(os.path.abspath(__file__))
                v_file = f'{file_dir}/covcache/vout{p}.txt'
                if os.path.exists(v_file):
                    print(f"Loading custom Sigma for gibbs model at {v_file}")
                    Sigma = np.loadtxt(f'{file_dir}/covcache/vout{p}.txt')
                else:
                    print(f"No custom Sigma available for gibbs model at {v_file}: using default (expect low power)")
            # Save gibbs_graph
            sample_kwargs['gibbs_graph'] = [dgprocess.gibbs_graph]

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
        result, S_matrices = compare_S_methods(
            Sigma=Sigma if not resample_sigma else None,
            beta = beta if not resample_beta else None,
            dgp_number=dgp_number,
            cutoff=cutoff,
            sample_kwargs=sample_kwargs,
            filter_kwargs=filter_kwargs,
            fstat_kwargs=fstat_kwargs,
            reps=reps,
            num_processes=num_processes,
            seed_start=seed_start,
            time0=time0,
            S_curve=S_curve,
        )
        new_dgp_kwargs['cutoff'] = cutoff # put this back for logging
        for key in dgp_keys:
            if key in sample_kwargs:
                continue
            result[key] = new_dgp_kwargs[key]
        all_results = all_results.append(
            result, 
            ignore_index = True
        )
    
        all_results.to_csv(output_path)

        # Reset printing options and print, then reset to defaults
        print_results(all_results)

        # Increment dgp number
        dgp_number += 1

    return all_results

if __name__ == '__main__':

    time0 = time.time()
    MAXIMUM_N = 0
    main(sys.argv)


