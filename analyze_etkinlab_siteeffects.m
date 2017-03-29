function results = analyze_etkinlab_siteeffects()
% Compute and return conditional correlations, nuisance correlations, and NSR on all subjects x sites. 
	
	studydata = load('~/COMET/adim-preprocessed/others/HCData4SiteEfectTests.mat'); 
	studydata.noise = load('~/COMET/adim-preprocessed/others/noisePCs4SiteEffectData_03222017.mat')	

	studydata.compcor1 = studydata.noise.noisePCs.noisePCs(1:20,:,:); 
	studydata.compcor2 = studydata.noise.noisePCs.noisePCs(21:end,:,:); 

	results = struct();
	results.Site1 = {};
	results.Site2 = {};
	for ii=1:length(studydata.subIDsSite1)
		X = squeeze(studydata.signalsSite1(ii,:,:))';
		Y = squeeze(studydata.compcor1(ii,:,:)); 
		try
			results.Site1{ii} = conditional_correlation(X,Y);
		catch me
			disp(me)
			results.Site1{ii} = {};
		end 
	end
	
	for ii=1:length(studydata.subIDsSite2)
		X = squeeze(studydata.signalsSite2(ii,:,:))';
		Y = squeeze(studydata.compcor2(ii,:,:)); 
		try
			results.Site2{ii} = conditional_correlation(X,Y); 
		catch me
			disp(me)
			results.Site2{ii} = {};
		end
	end
	
	
	% CAUTION
	% These diagnostics will not work until signal has been cleaned. 
	% 
	
	X1 = permute(studydata.signalsSite1,[3 2 1]); 
	X2 = permute(studydata.signalsSite2,[3 2 1]); 
	results.diagnostic1 = get_influence_diagnostic(X1); 								
	results.diagnostic2 = get_influence_diagnostic(X2);
	
	results.diagnostic3 = get_influence_diagnostic(cat(3,X1,X2));
										
	
	save('~/COMET/adim-preprocessed/others/HCData4SiteEffect_Decomposition_CompCorNuisance','results');
	
end

function output =  conditional_correlation(X,Y)
	% Only uses usual column standardize (i.e. correlation)
	output = struct()
		
	[Sigma results] = covariance.conditional_sample_covariance_separate(X, ...
									struct('verbose',false,...
											'nuisance',Y) ...	
											);
	 	
	output.corr = Sigma;
	output.nuisance = results.nCov;
	output.NSR = results.NSR
	
end

function output = get_influence_diagnostic(X)
	
	output = struct();
    
    [Sigmasb Sigmas] = covjackknife(X,[1 2 3]);
	
	output.nansubjects = isnan(squeeze(sum(sum(Sigmas,1),2))); 	
	Sigmasb(isnan(Sigmasb)) = 0; 
	Sigmas(isnan(Sigmas)) = 0;
	
	[output.influence output.norminfluence] = influence(Sigmasb,Sigmas);
	
end