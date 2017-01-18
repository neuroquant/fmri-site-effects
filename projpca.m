function [output] = projpca(Y, X, varargin)
	% Heterogeneity Test and Projected PCA
	% 
	% An example usecase would be to clean site effects per fMRI volume.
	% Here,
	%	 Y - subjects X ROIs 
	% 	 X - subjects X covariates
	% and output.u_it + output.mu_t would be the site corrected version of Y
	% 
	% Model Example: 
	% y_it = x'_i beta_t + mu_t + u_it
	% where 
	% 	- i=(1,...,p) variables
	% 	- t=(1,...,T) samples/time-varying observations
	% 	- x'_i is a d-dimensional covariate fixed for each i
	%   - beta_t is a regression coefficient common to all i
	% 
	% The above model can be re-written as a factor model 
	% y_it = g(X_i)'f_t + u_it
	% where
	% 	- g(X_i) is d+1 dimensional loadings (1,X_i')'
	%   - f_t is d+1 dimensional factors (mu_t,beta_t')' 
	% 
	% INPUTS
	% Data matrix Y should be n_subjects x n_rois
	% 
	% Covariate matrix X should be n_subjects x n_covariates. 
	% If you have a categorical site covariate, then columns of X should contain dummy variables encoding each site. 
	% 
	% 
	% OUPUTS
	% output.U 	contains residuals
	% output.Y  contains output.U + output.mu
	% output.beta contains the covariate or site effects
	% 
	% NOTE:
	% This function does not yet implement the projected PCA method 
	% but employs regression over observed categorical covariates
	% This function does not yet support corrections using continuous demographic (age,education) covariates
	% 
	% References: 
	% Fan, J., Liao, Y., & Wang, W. (2016). PROJECTED PRINCIPAL COMPONENT ANALYSIS IN FACTOR MODELS. Annals of Statistics, 44(1), 219–254. http://doi.org/10.1214/15-AOS1364
			

	% Elementary implementation given observed covariates and no latent site effects
	
	[p t] = size(Y);
	d = size(X,2); 
	
	% Center rows
	mu = mean(Y,1); 
	Yc = bsxfun(@minus,Y,mu); 
	
	% Solve for beta in Y = X*beta
	SigX = pinv(X'*X);  
	B = SigX*X'*Yc;

	% Residuals
	% Either Yc-XB = (I-P)Yc
	P = X*SigX*X';		
	U = (eye(p)-P)*Yc;
	
	
	output.Y = bsxfun(@plus,U,mu); 
	output.U = U; 
	output.mu = mu;
	output.beta = B;
end