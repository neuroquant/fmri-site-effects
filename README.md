# A toolbox to quantify and adjust for session and site effects in fMRI time-series

This toolbox implements a number of approaches to correct fMRI time-series or volumes for session specific artifacts and site effects, aka "batch Effects". 


## Contents

This toolbox has 3 types of functions. 

- Primary functions to measure reliability, implement alternatives to the standard correlation matrix that adjust for nuisance variation or artifacts. 

- Demo functions that illustrate the function API and usage of primary functions. These functions start with `demo_`

- Analysis functions that perform a detailed case study of a particular method or pipeline on a dataset. These functions start with `analyze_`

## Primary functions

- `siteEffects`: A class containing useful methods for diagnosis of site effects and other metrics of improvement. 
	- `siteEffects.within_group_error()` 
	- `siteEffects.test_retest_mccc()`

- `projpca`, `demo_projpca`: Multivariate correction on time-series data tensor. Provides residual time-series after regressing out site effects can be used to infer correlation matrices. *Disclaimer: This function is still being tested, unfinished.*

## Demo functions

- `demo_successive_normalization`: Demonstrates examples of calling ggmClass for data cleaning and correlation estimation. Is a wrapper around local functions that combine different forms of standardization of data matrix with a variety of covariance estimators. Exported figures can be found in `tmp/<date>/demo_successive_norm_*.png`

- `demo_conditional_correlation.m`: Demonstrates examples of calling ggmClass for data cleaning and correlation estimation. Is a wrapper around local functions that perform factor model based decomposition of the observed covariance by conditioning on covariates. 

- `demo_robust_covariance.m`: Demonstrates examples of calling ggmClass to call robust estimators of the covariance or correlation matrix. 


## Installation

Use the `--recursive` setting to make sure the `ggmClass` and `connectivity-diagnostics` submodule are added. *Disclaimer: You will need to have keyless ssh login to github to access private repositories.*

```bash
git clone --recursive git@github.com:TheEtkinLab/fmri-site-effects.git
```

Modify `setup.m` to specificy paths to `matlab-library` or your local equivalent. To add package dependencies, run setup using 

```MATLAB
run('setup.m')
```

