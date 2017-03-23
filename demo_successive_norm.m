function results =  demo_successive_norm(X,varargin)
	
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
	
	results{1}.method = 'Column Standardize';
	results{1}.output = standard_correlation(X)
	
	results{2}.method = 'Row-Column Standardize';
	results{2}.output = standard_correlation_sn(X)
	
	results{3}.method = 'Row-Column Sym Standardize';
	results{3}.output = standard_correlation_sym_sn(X)
		
	disp(sprintf('Frob. MSE: Sigma_c - Sigma_sym = %.3f',  ...
		norm(abs(results{2}.output.corr-results{3}.output.corr),'fro')));
	
	disp(sprintf('Frob. SSE: Sigma_c - Sigma_sym = %.3f', ...	
		sum(sum((results{2}.output.corr-results{3}.output.corr).^2)) ));
		
	figure('Position',[50 100 800 650]);
	subplot(2,3,1); 
	imagesc(results{1}.output.corr); axis equal image;
	colormap('winter'); 
	title(results{1}.method);
	subplot(2,3,2); 
	imagesc(results{2}.output.corr); axis equal image;
	colormap('winter')
	title(results{2}.method);
	subplot(2,3,3); 
	imagesc(results{3}.output.corr); axis equal image;
	colormap('winter')
	title(results{3}.method);
	subplot(2,3,4); 	
	histogram(results{1}.output.corr(:),'Normalization','pdf'); 
	axis equal image;
	title(results{1}.method); xlabel('correlation'); ylabel('pdf')
	subplot(2,3,5); 
	histogram(results{2}.output.corr(:),'Normalization','pdf'); 
	axis equal image;
	title(results{2}.method);xlabel('correlation'); ylabel('pdf')
	subplot(2,3,6); 
	histogram(results{3}.output.corr(:),'Normalization','pdf'); 
	axis equal image;
	colormap('winter')
	title(results{3}.method);xlabel('correlation'); ylabel('pdf')
	
	fname = ['tmp' filesep datestr(now,'dd-mmm-yyyy-HHMMSS')]; 
	if(~exist(fname,'dir'))
		mkdir(fname)
	end
	exportfun(fullfile(fname,mfilename));	
	
end


function output =  standard_correlation(X)
	
	output = struct()
	
	X = standardize.standardize_cols(X);
	
	ggmobj = GGM(X);
	ggmobj.MLECovEstimate();	 	
	output.corr = ggmobj.Sigma;
	
end

function output = standard_correlation_sn(X)
	
	output = struct();
	
	X = standardize.successive_normalize(X);
	
	ggmobj = GGM(X);
	ggmobj.MLECovEstimate();
	output.corr = ggmobj.Sigma;
	
end

function output = standard_correlation_sym_sn(X)
	
	output = struct();
	
	X = standardize.successive_normalize(X, ...
						struct('method','sym'));
	
	ggmobj = GGM(X);
	ggmobj.MLECovEstimate();
	output.corr = ggmobj.Sigma;
	
end

