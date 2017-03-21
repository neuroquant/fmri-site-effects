# Adjust fMRI time-series for session and site effects

This toolbox implements a number of approaches to correct fMRI time-series or volumes for session specific artifacts and site effects. 

## Contents

- `siteEffects`: A class containing useful methods for diagnosis of site effects and other metrics of improvement. 
	- `siteEffects.within_group_error()` 

- `projPCA`: Multivariate correction on time-series data tensor. Provides residual time-series after regressing out site effects can be used to infer correlation matrices. *Disclaimer: This function is still being tested, unfinished.*