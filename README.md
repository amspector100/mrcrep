## MRC Replications

This repository contains all the codes to replicate the experiments from ___. Running these scripts requires a python environment with ``knockpy`` installed: please see https://github.com/amspector100/knockpy for installation. These scripts run correctly with knockpy=1.0.1. A few figures were generated using the attached .ipynb notebook instead,  

### Overview

The mrcrep.py file contains an experimentation framework as well as a custom argument parser which makes it easy to compare different types of knockoffs (namely, methods based on minimizing covariace vs. methods which minimize the reconstructability of the features).

mrcrep.py broadly takes four types of arguments:
1. Arguments which affect the data generating process. All arguments of this form should follow the convention `--{argument_name}_dgp {value}`. For example, to set the coeff_size parameter to 1, you should add ``--coeff_size_dgp 500``. 
2. Arguments to pass to the KnockoffFilter class from knockpy. All arguments of this form should follow the convention `--{argument_name}_filter {value}`. For example, to set the feature statistic to deeppink statistics, you should add ``--fstat_filter deeppink``.
3. Arguments to pass to your chosen feature statistic class. All arguments of this form should follow the convention `--{argument_name}_fstat {value}`. For example, to set get the (default) lasso statistic to use the LARS solver, you should add ``--_filter deeppink``.
4. A `--description` argument which must be the last argument.

For all arguments, you can add multiple values for an argument by writing brackets around the values. For example, the following script will run knockoffs on a linear Gaussian model with p equal to 500, 2000 and n equal to 100, 500.

```
python3 mrcrep.py --p_dgp [500,2000] --n_dgp [100,500]
```

The output file will be saved as a .csv in the data directory.

### How to replicate specific figures

The figures correspond to:

- Figure 1: see scripts/gaussianplot1.sh
- Figure 2: see the .ipynb notebook
- FIgure 3: see scripts/equiplot1.sh
- Figure 4: see scripts/gaussianplot1.sh
- Figure 5: see scripts/gaussian_nonlinear.sh
- Figure 6: see scripts/robustness.sh
- Figure 7: see scripts/robustness.sh
- Figure 8: see scripts/fx.sh
- Figure 9: see scripts/nongaussian_linear.sh
- Figure 10: see scripts/equiplot3.sh
- Figure 11: see the .ipynb notebook
- Figure 12: see the .ipynb notebook
- Figure 13: see scripts/equiplot2.sh
- Figure 14: see scripts/cicomp.sh
- Figure 15: see the .ipynb notebook
- Figure 16: see scripts/dsliu2020.sh and xing2019rep.R
- Figure 17: see scripts/dsliu2020.sh
- Figure 18: see scripts/ar1corrplot.sh
- Figure 19: see scripts/ar1corrplot.sh
- Figure 20: see scripts/gaussian_nonlinear.sh
- Figure 21: see scripts/gaussian_nonlinear.sh
- Figure 22: see scripts/nongaussian_linear.sh