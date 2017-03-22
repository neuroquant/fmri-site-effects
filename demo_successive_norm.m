function demo_successive_norm(X)
	
	results = {};
	
	result{1}.method = 'standard';
	result{1}.output = standard_correlation(X)
	
	result{2}.method = 'standard';
	result{2}.output = standard_correlation_sn(X)
	
	
	figure;
	subplot(1,2,1); 
	imagesc(result{1}.output.corr); axis equal image;
	colormap('winter'); 
	title(result{1}.method);
	subplot(1,2,2); 
	imagesc(result{2}.output.corr); axis equal image;
	colormap('winter')
	title(result{2}.method);
		
	
end


function output =  standard_correlation(X)
	
	output = struct()
	ggmobj = GGM(X); 	
	output.corr = ggmobj.Sigma;
	
end

function output = standard_correlation_sn(X)
	
	output = struct();
	ggmobj = GGM(X);
	ggmobj.MLECovEstimate();
	output.corr = ggmobj.Sigma;
	
end
