function results = analyze_etkinlab_sess_length()
    
    studyname = 'best-eegfmri-cpac'
    addpath(fullfile(getenv('HOME'),'MATLAB','ggmClass'));

    
	if(exist('brewermap'))	
		colormapfun = @()(flipud(brewermap(length(colormap),'RdYlBu')));
		close all;
	else
		colormapfun = @winter;
	end

    % test re-test mccc 
    exportfun = @(filename)(print('-dpng','-r150',filename));
    %fname = ['tmp' filesep datestr(now,'dd-mmm-yyyy')];
    fname = ['tmp' filesep '15-Oct-2018'];    
    if(~exist(fname,'dir'))
        mkdir(fname)
    end
    fname = [fname filesep studyname]; 
    if(~exist(fname,'dir'))
        mkdir(fname)
    end
    
    nmethods = 1;
    methodnames = {'sn_denoised'};    
     load(['tmp' filesep '15-Oct-2018' filesep 'best-eegfmri-cpac-trt-t840/etkinlab_analysis/data_etkinlab_analysis_best-eegfmri-cpac_rfmri.mat']);
    
    nsubjects = length(studydata.session1.X);
    p = size(studydata.session1.X{1},2);
    
    
    % Get community labels
    nodeidx = 1:size(studydata.session1.X{1},2);
    tmp_community_parts = ...
         regexp(studydata.roilabels.table.Var2,{'_'}, 'split');
    for ii=1:length(tmp_community_parts)     
        communitylabels{ii} = tmp_community_parts{ii}{3};
    end
    communitylabel = grp2idx(communitylabels);
    communitycolor = communitylabel;
    
    
    results = {};
    results.cpac = {};
    results.community = communitylabel;
    results.communitycolor = communitycolor;
    
    
    % Get network estimate in increments of 80 measurements 
    % (80*.7 = 56 seconds)
    
    n_step = 80;
    n_minutes = 10;
    n_trials = 1;

    for minuteno=1:1:n_minutes
        tmp_kendallW = zeros(1,nsubjects);
        tmp_mccc = zeros(p,1);
        for trialno=1:n_trials
            corr_mat = helper_corr_session(studydata,minuteno*n_step);
            results.networks = corr_mat;
            graphs = results.networks(nodeidx,:,:,:);
            tmp_kendallW = tmp_kendallW + helper_kendallW(graphs,nodeidx);
            trt_reliability = helper_node_mccc(graphs);
            tmp_mccc = tmp_mccc + real(trt_reliability.mccc);
        end
        results.kendallW(:,minuteno) = tmp_kendallW/n_trials;
        results.mccc(:,minuteno) = tmp_mccc/n_trials;
    end
    
    dlmwrite([fname filesep 'trt_' studyname '_mccc.txt'],results.mccc,'precision',6);
    dlmwrite([fname filesep 'trt_' studyname '_kendallw.txt'],results.kendallW,'precision',6);
    
    mccc_summary(:,1) = mean(results.mccc,1);
    mccc_summary(:,2) = prctile(results.mccc,[5]);
    mccc_summary(:,3) = prctile(results.mccc,[95]);
    
    kw_summary(:,1) = mean(results.kendallW,1);
    kw_summary(:,2) = prctile(results.kendallW,[5]);
    kw_summary(:,3) = prctile(results.kendallW,[95]);

    dlmwrite([fname filesep 'trt_' studyname '_mccc_summary.txt'],mccc_summary,'precision',4);
    dlmwrite([fname filesep 'trt_' studyname '_kendallw_summary.txt'],kw_summary,'precision',4);


    %   save([fname filesep 'trt_' studyname '_rfmri.mat'],'results','-append');
    % end
    %
    % try
    %     mean(results.cpac.sn_denoised.trt_reliability.kendallsw-...
    %         results.cpac.observed.trt_reliability.kendallsw,2);
    % catch
    %
    % end
    
end



function corr_mat = helper_corr_session(studydata,n_time)
    
    [n p] = size(studydata.(['session1']).X{1});
    nsubjects = length(studydata.(['session1']).X);
    step_t = n_time
    verbose = false;
    
    nsessions = 2;
    corr_mat = zeros(p,p,nsubjects,nsessions);
    for sessionno=1:nsessions
    	for ii=1:nsubjects
            if(verbose)
                ii
                ['session' num2str(sessionno)]
            end
            X = studydata.(['session' num2str(sessionno)]).X{ii};
            start_idx = randi([9 max(9,(n-9-step_t))]);
            X = X(start_idx+1:start_idx+n_time,:);
            
            try
                corr_mat(:,:,ii,sessionno) = ...
                            standard_correlation_sn(X);
            catch me
                disp(me)
                results.Site1{ii,sessionno} = {};
                disp('Session')
                sessionno
                disp('Subject')
                ii
            end
            % if(ii==1)
            %     demo_conditional_correlation(X,Y);
            % end
    	end
    end
    
    
end




function output = helper_kendallW(graphs,nodeidx)
    
    
    p = size(graphs,1);
    nsubjects = size(graphs,3);
    tmp_kendallW = nan(1,nsubjects);
    
    
    ba_nodes = nodeidx;
    trid_idx = find(reshape(triu(ones(p,p),1),[p^2 1]));
    for subjectno=1:nsubjects
        rater1 = graphs(ba_nodes,ba_nodes,subjectno,1);
        rater2 = graphs(ba_nodes,ba_nodes,subjectno,2);
        rater1 = rater1(trid_idx);
        rater2 = rater2(trid_idx);

        tmp_kendallW(subjectno) = reliability.kendallsW(cat(2,rater1,rater2));
    end
    
    output = tmp_kendallW;
    
end



function output = helper_mccc()
    
    
    
    
    
    
    
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

function trt_reliability = helper_node_mccc(graphs)
    % 
    % Inputs
    %   - graphs is p x p x subjects x {rater|session}    
    
    p = size(graphs,1);
    nsubjects = size(graphs,3);
    if(p<=50)
        nodeidx = 1:p;
    end
    verbose = false;
    
    nboot = 50;
    mccc = zeros(p,1); 
    mccc_ci = zeros(p,2);
    mccc_var = zeros(p,1);
    mccc_wvar = zeros(p,1);
    mccc_boot = zeros(p,nboot);
    similarity = zeros(2*nsubjects,2*nsubjects);
    
    for nodeno=1:p;
        jj = nodeno;
        %jj = nodeidx(nodeno);
        features = squeeze(graphs(jj,:,:,:));
        similarity = similarity + cov(cat(2,features(:,:,1),features(:,:,2)));
        
    end
    similarity = similarity./p;    
    for nodeno=1:p
        %jj = nodeidx(nodeno)
        if(verbose)
            jj = nodeno
        else
            jj = nodeno;
            if(mod(jj,10)==0)
                disp(sprintf('MCCC %2.2f Complete',jj/p));
            end
        end
        features = squeeze(graphs(jj,:,:,:));
        [mccc(jj), MCCC, mccc_ci(jj,:), mccc_boot(jj,:)] = ...
             reliability.mccc(features(:,:,1)',features(:,:,2)');
        mccc_var(jj) = MCCC.normVind;
        mccc_wvar(jj) = MCCC.normVdep;
        if(verbose)
            disp(sprintf('MCCC: %2.4f, CI: (%2.4f,%2.4f), WithinVar: %.4f, TotalVar: %.4f', ...
        mccc(jj), mccc_ci(jj,1), mccc_ci(jj,2), mccc_wvar(jj), mccc_var(jj)))
        end
    end
    trt_reliability.mccc = real(mccc); 
    trt_reliability.mccc_tvar = mccc_var;
    trt_reliability.mccc_wvar = mccc_wvar;
    trt_reliability.mccc_ci = real(mccc_ci);
    trt_reliability.mccc_boot = real(mccc_boot);
    trt_reliability.similarity = similarity;
end