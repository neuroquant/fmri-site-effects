function results =  demo_conditional_covariance(X,Y,varargin)


	opts = struct();
	opts.exportfig = true;
	opts.exportfun = @(fname)(print('-dpng','-r300',fname));	
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

	fname = ['tmp' filesep datestr(now,'dd-mmm-yyyy-HHMM')]; 
	if(~exist(fname,'dir'))
		mkdir(fname)
	end
    opts.outputdir = fname;

	results = {};
	
	results{1}.method = 'Sample Correlation';
	results{1}.output = standard_correlation(X);
	
	results{2}.method = 'RC, Sample Correlation';
	results{2}.output = standard_correlation_sn(X);

	results{3}.method = 'Nuisance Correlation';
	results{3}.output = conditional_correlation(X,Y);

	results{4}.method = 'Denoised Sample Correlation';
	results{4}.output = results{3}.output;
	results{4}.output.corr = results{3}.output.corr;

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
	
    addpath(genpath('../MATLAB/packages/spreadFigures/tightfig/'));
    figure(1);
	set(gcf,'Position',[10 150  1200 650]);
	subplot(2,2,1); 
	imagesc(results{1}.output.corr); 
    colormap(colormapfun()); axis image; colorbar;;
	title(results{1}.method);
    set(gca,'fontsize',16);
	subplot(2,2,2); 
	imagesc(results{3}.output.nuisance); 
    colormap(colormapfun()); axis image; colorbar;
	title(results{3}.method);
    set(gca,'fontsize',16);
	subplot(2,2,3); 
	imagesc(results{4}.output.corr); 
    colormap(colormapfun()); axis image; colorbar;;
	title(results{4}.method);
    set(gca,'fontsize',16);
	subplot(2,2,4); 
	imagesc(results{2}.output.corr); 
    colormap(colormapfun()); axis image; colorbar;
	title(results{2}.method);
    set(gca,'fontsize',16);
	
    exportfun(fullfile(fname,mfilename));
    % if(exist('tightfig'))
    %     tightfig;
    % end
    
    figure(2)
	set(gcf,'Position',[10 150  1200 650]);
	subplot(2,2,1); 
    upper_idx = find(reshape(triu(ones(size(results{1}.output.corr)),1),  ...
                                [1 numel(results{1}.output.corr)])); 	
	histogram(results{1}.output.corr(upper_idx),'Normalization','probability'); 
	ylim([-1.2 1.2]);axis tight; 
	title(results{1}.method); xlabel('correlation'); ylabel('pdf')
    set(gca,'fontsize',16);
	subplot(2,2,2); 
	histogram(results{3}.output.nuisance(upper_idx),'Normalization','probability'); 
	ylim([-1.2 1.2]);axis tight; 
	title(results{3}.method);xlabel('correlation'); ylabel('pdf')
    set(gca,'fontsize',16);
	subplot(2,2,3); 
	histogram(results{4}.output.corr(upper_idx),'Normalization','probability'); 
	ylim([-1.2 1.2]);axis tight; 
	title(results{4}.method);xlabel('correlation'); ylabel('pdf')
    set(gca,'fontsize',16);
	subplot(2,2,4); 
	histogram(results{2}.output.corr(upper_idx),'Normalization','probability'); 
	ylim([-1.2 1.2]);axis tight; 
	title(results{2}.method);xlabel('correlation'); ylabel('pdf')
    set(gca,'fontsize',16);
    % if(exist('tightfig'))
    %     tightfig;
    % end
	exportfun(fullfile(fname,[mfilename '2']));


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
	output = struct();
		
	%[Sigma results] = covariance.conditional_sample_covariance_separate(X);
	
		
	[~, results] = covariance.conditional_sample_covariance_separate(X, ...
									struct('verbose',true,...
											'nuisance',Y) ...	
											);
	output.corr = results.Sigma;
	output.nuisance = results.nCov;
	output.corr2 = covariance.mle_sample_covariance(results.X_perpY, ...
												struct('standardize',false));
	
end