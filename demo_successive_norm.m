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
	
	
	figure;
	subplot(1,2,1); 
	imagesc(results{1}.output.corr); axis equal image;
	colormap('winter'); 
	title(results{1}.method);
	subplot(1,2,2); 
	imagesc(results{2}.output.corr); axis equal image;
	colormap('winter')
	title(results{2}.method);
	
	fname = ['tmp' filesep datestr(now,'DD-MM-YYYY-hhmmss')]; 
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
	
	X = standardize.successive_normalize(X, ...
						struct('method','sym'));
	
	ggmobj = GGM(X);
	ggmobj.MLECovEstimate();
	output.corr = ggmobj.Sigma;
	
end
