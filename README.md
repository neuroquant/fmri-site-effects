# Adjust fMRI time-series for session and site effects

This toolbox implements a number of approaches to correct fMRI time-series or volumes for session specific artifacts and site effects. 


## Contents

- `siteEffects`: A class containing useful methods for diagnosis of site effects and other metrics of improvement. 
	- `siteEffects.within_group_error()` 

- `projpca`, `demo_projpca`: Multivariate correction on time-series data tensor. Provides residual time-series after regressing out site effects can be used to infer correlation matrices. *Disclaimer: This function is still being tested, unfinished.*

- `demo_successive_normalization`: Demonstrates examples of calling ggmClass for data cleaning and correlation estimation. Is a wrapper around local functions that combine different forms of standardization of data matrix with a variety of covariance estimators. Exported figures can be found in `tmp/<date>/demo_successive_norm_*.png`

## Installation

Use the `--recursive` setting to make sure the `ggmClass` and `connectivity-diagnostics` submodule are added. *Disclaimer: You will need to have keyless ssh login to github to access private repositories.*

```bash
git clone --recursive git@github.com:TheEtkinLab/fmri-site-effects.git
```

Modify `setup.m` to specificy paths to `matlab-library` or your local equivalent. To add package dependencies, run setup using 

```MATLAB
run('setup.m')
```

