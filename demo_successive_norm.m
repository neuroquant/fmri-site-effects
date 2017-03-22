function demo_successive_norm(X)
	
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
