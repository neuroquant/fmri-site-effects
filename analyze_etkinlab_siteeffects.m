function results = analyze_etkinlab_siteeffects()
% Compute and return conditional correlations, nuisance correlations, and NSR on all subjects x sites. 
    datadir = 'data';
	studydata = load(fullfile(datadir,'HCData4SiteEfectTests.mat')); 
	studydata.noise = load(fullfile(datadir,'noisePCs.mat'));	
    studydata.yeo_order = load(fullfile(datadir,'community_ordering'));
    reorder_community = studydata.yeo_order.orderbycommunities;

	studydata.compcor1 = studydata.noise.noisePCs.noisePCs(1:20,:,:); 
	studydata.compcor2 = studydata.noise.noisePCs.noisePCs(21:end,:,:); 

	results = struct();
	results.Site1 = {};
	results.Site2 = {};
    observed_correlation = [];
    denoised_correlation = []; 
    denoised_correlation2 = []; 
    nuisance_correlation = [];
    nsr = [];
    sitelabels = zeros(44,1); sitelabels(1:20) = 1; sitelabels(21:44) = 2; 
    includeROIs1 = find(sum(squeeze(abs(std(studydata.signalsSite1,[],3))>1e-2),1)==length(studydata.subIDsSite1));
    includeROIs2 = find(sum(squeeze(abs(std(studydata.signalsSite2,[],3))>1e-2),1)==length(studydata.subIDsSite2));
    includeROIs = intersect(includeROIs1,includeROIs2); 
    [loca locb] = ismember(reorder_community,includeROIs); 
    choose_community = reorder_community<=3;
    
    compcor_dims = size(studydata.compcor1); 
    %Y = get_shared_nuisance(squeeze(mean(studydata.signalsSite1,2))');
    Y = [];
	for ii=1:length(studydata.subIDsSite1)
        ii
		X = squeeze(studydata.signalsSite1(ii,reorder_community(loca(:)),:))';  
		Y = cat(2,Y,squeeze(studydata.compcor1(ii,:,:))); 
        try
            observed_correlation(ii,:,:) = standard_correlation(X); 
            results.Site1{ii} = conditional_correlation(X,Y);
            denoised_correlation(ii,:,:) = results.Site1{ii}.corr;
            denoised_correlation2(ii,:,:) = standard_correlation_sn(X);
            nuisance_correlation(ii,:,:) = results.Site1{ii}.nuisance;                      nsr(ii) = results.Site1{ii}.NSR;
        catch me
            disp(me)
            results.Site1{ii} = {};
        end
        % if(ii==1)
        %     demo_conditional_correlation(X,Y);
        % end
	end
	
    baseidx = size(denoised_correlation,1); 
    %Y = get_shared_nuisance(squeeze(mean(studydata.signalsSite2,2))');
    Y = [];
	for ii=1:length(studydata.subIDsSite2)
		X = squeeze(studydata.signalsSite2(ii,reorder_community(loca(:)),:))';
		Y = cat(2,Y,squeeze(studydata.compcor2(ii,:,:))); 
		try
            observed_correlation(baseidx+ii,:,:) = standard_correlation(X); 
			results.Site2{ii} = conditional_correlation(X,Y); 
            denoised_correlation(baseidx+ii,:,:) = results.Site2{ii}.corr;      
            denoised_correlation2(baseidx+ii,:,:) = standard_correlation_sn(X); 
            nuisance_correlation(baseidx+ii,:,:) = results.Site2{ii}.nuisance; 
            nsr(baseidx+ii) = results.Site2{ii}.NSR;
		catch me
			disp(me)
			results.Site2{ii} = {};
		end
        % if(ii==1)
       %      demo_conditional_correlation(X,Y);
       %  end
	end
	    
	
    results.denoised = denoised_correlation; 
    results.denoised2 = denoised_correlation2; 
    results.nuisance = nuisance_correlation; 
    results.observed = observed_correlation;
    results.nsr = nsr; 
    
    
	if(exist('brewermap'))	
		colormapfun = @()(flipud(brewermap(length(colormap),'RdYlBu')));
		close all;
	else
		colormapfun = @winter;
	end
	
    
    addpath(genpath('~/MATLAB/packages/spreadFigures/tightfit/'));
    figure(1);
	%set(gcf,'Position',[10 450  1200  300]);
	set(gcf,'Position',[10 450  850 850]);    
    subplot(2,2,1); 
    site_effect{1} = detect_site_effect(results.observed,sitelabels); 
    imagesc(real(site_effect{1}.similarity)); 
    colormap(colormapfun()); colorbar; axis image; 
    title('Similarity Matrix (Observed)','fontsize',24);
    xlabel(sprintf('(wit,bet,rat) = (%.2f, %.2f, %.2f)', ...
                    site_effect{1}.within,  ...
                    site_effect{1}.between, ...
                    site_effect{1}.ratio));
    set(gca,'fontsize',16);
    
    subplot(2,2,2); 
    site_effect{2} = detect_site_effect(results.nuisance,sitelabels); 
    imagesc(real(site_effect{2}.similarity)); 
    colormap(colormapfun()); colorbar; axis image; 
    title('Similarity Matrix (Nuisance)','fontsize',24);    
    xlabel(sprintf('(wit,bet,rat) = (%.2f, %.2f, %.2f)', ...
                    site_effect{2}.within,  ...
                    site_effect{2}.between, ...
                    site_effect{2}.ratio));
    set(gca,'fontsize',16);
    
    
    subplot(2,2,3); 
    site_effect{3} = detect_site_effect(results.denoised,sitelabels); 
    imagesc(site_effect{3}.similarity); 
    colormap(colormapfun()); colorbar; axis image; 
    title('Similarity Matrix (Denoised)','fontsize',24);    
    xlabel(sprintf('(wit,bet,rat) = (%.2f, %.2f, %.2f)', ...
                    site_effect{3}.within,  ...
                    site_effect{3}.between, ...
                    site_effect{3}.ratio));
    set(gca,'fontsize',16);
    
    subplot(2,2,4); 
    site_effect{4} = detect_site_effect(results.denoised2,sitelabels); 
    imagesc(site_effect{4}.similarity); 
    colormap(colormapfun()); colorbar; axis image;
    title('Similarity Matrix (SN)','fontsize',24);    
    xlabel(sprintf('(wit,bet,rat) = (%.2f, %.2f, %.2f)', ...
                    site_effect{4}.within,  ...
                    site_effect{4}.between, ...
                    site_effect{4}.ratio));
    set(gca,'fontsize',16);
    
    % if(exist('tightfig'))
    %     tightfig;
    % end
    
    results.site_effect = site_effect;

    exportfun = @(fname)(print('-dpng','-r600',fname));
	fname = ['tmp' filesep datestr(now,'dd-mmm-yyyy-HHMM')]; 
	if(~exist(fname,'dir'))
		mkdir(fname)
	end
	exportfun(fullfile(fname,[mfilename '1']));
    
    %     figure;
    %     % stem([1:length(nsr)],nsr); hold on;
    %     % scatter([1:length(nsr)],nsr,25); hold off;
    %     % ylim([0 1.1]); xlim([0 length(nsr)+2]);
    %     % set(gca,'fontsize',24);
    %     % colormap(colormapfun()); axis image; colorbar;
    %     boxplot(nsr,'group',sitelabels,'color',sitelabels)
    %     title('Nuisance to Signal Ratio by Site')
    % exportfun(fullfile(fname,[mfilename '2']));
    %
    %     results(1).site_effect = site_effect;
    
	% CAUTION
	% These diagnostics will not work until signal has been cleaned. 
	% 
	
    % X1 = permute(studydata.signalsSite1,[3 2 1]);
    % X2 = permute(studydata.signalsSite2,[3 2 1]);
    % results.diagnostic1 = get_influence_diagnostic(X1);
    % results.diagnostic2 = get_influence_diagnostic(X2);
    %
    % results.diagnostic3 = get_influence_diagnostic(cat(3,X1,X2));
    %
    %
    % save('~/COMET/adim-preprocessed/others/HCData4SiteEffect_Decomposition_CompCorNuisance','results');
	
end

function output =  conditional_correlation(X,Y)
	% Only uses usual column standardize (i.e. correlation)
	output = struct();
		
	[Sigma results] = covariance.conditional_sample_covariance_separate(X, ...
									struct('verbose',false,...
											'nuisance',Y) ...	
											);
	 	
	output.corr = results.Sigma;
	output.nuisance = results.nCorr; 
	output.corr2 = covariance.mle_sample_covariance(results.X_perpY, ...
												struct('standardize', false));
	output.NSR = results.NSR;
	
end


function Sigma = standard_correlation(X)
	% Usual standard correlation matrix
	
	
	[Sigma results] = covariance.mle_sample_covariance(X, ...
												struct('standardize',false));
	
	
end


function Sigma = standard_correlation_sn(X)
	% Automatically applies row-first successive norm
	%standardize.successive_normalize(X');
	
	[Xnew] = standardize.successive_normalize(X');    
	[Sigma results] = corr(Xnew');
    
    % [Sigma results] = covariance.mle_sample_covariance(X,
%                                                 struct('standardize',true));
%
end


function site_effect = detect_site_effect(X,sitelabels)

    corrfun = @(X)(corr(X)); %@(X)(covariance.rank_sample_covariance(X(:,1:min(50,size(X,2)))));
    upper_idx = find(reshape(triu(ones(size(X,2), size(X,3)),1), [1 size(X,2)^2])); 
    X = reshape(X,[size(X,1) size(X,2)*size(X,3)]);
    tmpSimilarity = corrfun(X(:,upper_idx)');    
    [hom sep mw] = compareWithinAndBetweenGroupsSim(tmpSimilarity,sitelabels);
    
    site_effect.similarity = tmpSimilarity;           
    site_effect.within = hom; 
    site_effect.between = sep; 
    site_effect.stat = mw;  
    site_effect.ratio = (hom-sep)/(hom+sep); 
    
    
end

function [nuisance] = get_shared_nuisance(Y); 
   
   % Y is time-series x subjects
   Yz = zscore(Y')'; 
   [U S V] = svd(Y);  
   nfactors = min(size(Y,2),5);   
   nuisance = U(1:nfactors,:); 
   nuisance = reshape(nuisance,[size(Y,1) nfactors]); 
    
end

function output = get_influence_diagnostic(X)
	
	output = struct();
    
    [Sigmasb Sigmas] = covjackknife(X,[1 2 3]);
	output.nansubjects = isnan(squeeze(sum(sum(Sigmas,1),2))); 	
	Sigmasb(isnan(Sigmasb)) = 0; 
	Sigmas(isnan(Sigmas)) = 0;
	
	[output.influence output.norminfluence] = influence(Sigmasb,Sigmas);
	
end