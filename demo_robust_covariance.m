function results =  demo_robust_covariance(X,varargin)
	
	opts = struct();
	opts.exportfig = true;
	opts.exportfun = @(fname)(print('-dpng','-r150',fname));	
	% if(exist('process_options'))
	% 	[exportfig exportfun opts] = process_options(varargin,  ...
	% 									'exportfig', opts.exportfig, ...
	% 									'exportfun', opts.exportfun ...
	% 									);
	% else
	% 	warning('KPMTools from matlab-library needed to process optional arguments');
	% end
	exportfig = opts.exportfig;
	exportfun = opts.exportfun;
	
	results = {};
	
	results{1}.method = 'Sample Correlation';
	results{1}.output = standard_correlation(X);
	
	results{2}.method = 'RC, Sample Correlation';
	results{2}.output = standard_correlation_sn(X);

	results{3}.method = 'Robust Correlation';
	results{3}.output = robust_correlation(X);
	
	results{4}.method = 'RC, Robust Correlation';
	results{4}.output = robust_correlation_sn(X);
		
	disp(sprintf('Frob. MSE: Sigma_std - Sigma_rob = %.3f',  ...
		norm(abs(results{2}.output.corr-results{4}.output.corr),'fro')));
	
	disp(sprintf('Frob. SSE: Sigma_std - Sigma_rob = %.3f', ...	
		sum(sum((results{2}.output.corr-results{4}.output.corr).^2)) ));
		
	
	if(exist('brewermap'))	
		colormapfun = @()(brewermap(length(colormap),'RdYlBu'));
		close all;
	else
		colormapfun = @winter;
	end
	
	figure('Position',[50 100 800 650]);
	subplot(2,4,1); 
	imagesc(results{1}.output.corr); axis equal image;
	colormap(colormapfun()); 
	title(results{1}.method);
	subplot(2,4,2); 
	imagesc(results{2}.output.corr); axis equal image;
	colormap(colormapfun())
	title(results{2}.method);
	subplot(2,4,3); 
	imagesc(results{3}.output.corr); axis equal image;
	colormap(colormapfun())
	title(results{3}.method);
	subplot(2,4,4); 
	imagesc(results{4}.output.corr); axis equal image;
	colormap(colormapfun())
	title(results{4}.method);
	subplot(2,4,5); 	
	histogram(results{1}.output.corr(:),'Normalization','probability'); 
	ylim([-1.2 1.2]);axis tight; 
	title(results{1}.method); xlabel('correlation'); ylabel('pdf')
	subplot(2,4,6); 
	histogram(results{2}.output.corr(:),'Normalization','probability'); 
	ylim([-1.2 1.2]);axis tight; 
	title(results{2}.method);xlabel('correlation'); ylabel('pdf')
	subplot(2,4,7); 
	histogram(results{3}.output.corr(:),'Normalization','probability'); 
	ylim([-1.2 1.2]);axis tight; 
	title(results{3}.method);xlabel('correlation'); ylabel('pdf')
	subplot(2,4,8); 
	histogram(results{4}.output.corr(:),'Normalization','probability'); 
	ylim([-1.2 1.2]);axis tight; 
	title(results{4}.method);xlabel('correlation'); ylabel('pdf')
	
	fname = ['tmp' filesep datestr(now,'dd-mmm-yyyy-HHMMSS')]; 
	if(~exist(fname,'dir'))
		mkdir(fname)
	end
	exportfun(fullfile(fname,mfilename));	
	
end


function output =  standard_correlation(X)
	% Only uses usual column standardize (i.e. correlation)
	output = struct()
		
	[Sigma results] = covariance.mle_sample_covariance(X,struct('standardize',false));
	 	
	output.corr = Sigma;
	
end

function output = standard_correlation_sn(X)
	% Automatically applies row-first successive norm
	%standardize.successive_normalize(X');
	
	output = struct();
	
	
	[Sigma results] = covariance.mle_sample_covariance(X,struct('standardize',true));
	
	output.corr = Sigma;
	
end


function output =  robust_correlation(X)
	% Only uses usual column standardize (i.e. correlation)
	output = struct()
	
	%X = standardize.standardize_cols(X);
	
	[Sigma results] = covariance.tylermle_sample_covariance(X,struct('standardize',false));
	 	
	output.corr = Sigma;
	
end

function output = robust_correlation_sn(X)
	% Automatically applies row-first successive norm
	%standardize.successive_normalize(X');
	
	output = struct();
	
	
	[Sigma results] = covariance.tylermle_sample_covariance(X,struct('standardize',true));
	
	output.corr = Sigma;
	
end

