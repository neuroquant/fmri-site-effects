function results =  demo_conditional_covariance(X,Y,varargin)


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

	results{3}.method = 'Conditional Correlation';
	results{3}.output = conditional_correlation(X,Y);

	results{4}.method = 'NSR plus Sample Correlation';
	results{4}.output = results{3}.output;
	results{4}.output.corr = results{3}.output.corr2;

	disp(sprintf('Frob. MSE: Sigma_std - Sigma_cond = %.3f',  ...
		norm(abs(results{1}.output.corr-results{3}.output.corr),'fro')));
	
	disp(sprintf('Frob. SSE: Sigma_std - Sigma_cond = %.3f', ...	
		sum(sum((results{1}.output.corr-results{3}.output.corr).^2)) ));
		
	
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



function output = standard_correlation(X)
	% Usual standard correlation matrix
	
	output = struct();
	
	[Sigma results] = covariance.mle_sample_covariance(X, ...
												struct('standardize',false));
	
	output.corr = Sigma;
	
end


function output = standard_correlation_sn(X)
	% Automatically applies row-first successive norm
	%standardize.successive_normalize(X');
	
	output = struct();
	
	[Sigma results] = covariance.mle_sample_covariance(X, ...
												struct('standardize',true));
	
	output.corr = Sigma;
	
end


function output =  conditional_correlation(X,Y)
	% Only uses usual column standardize (i.e. correlation)
	output = struct()
		
	%[Sigma results] = covariance.conditional_sample_covariance_separate(X);
	
		
	[Sigma results] = covariance.conditional_sample_covariance_separate(X, ...
									struct('verbose',false,...
											'nuisance',Y) ...	
											);
	 	
	output.corr = Sigma;
	output.nuisance = results.nCov;
	output.corr2 = covariance.mle_sample_covariance(results.X_perpY, ...
												struct('standardize',true));
	
end