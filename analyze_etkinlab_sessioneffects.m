function results = analyze_etkinlab_sessioneffects(varargin)
% Compute and return conditional correlations, nuisance correlations, and NSR on all subjects x sites. 

    % import ggmClass/standardize.*
    % import ggmClass/covariance.*
    addpath(fullfile(getenv('HOME'),'MATLAB','ggmClass'));

    % % DANA healthy controls
    % tmpdata = load(['../prelim2016-resting-predicts-tms/' ...
    %                     'tmp_workspace_1.mat'],'obj');
    % tmpnames = regexp(tmpdata.obj.subjlabels,'_','split');
    % for ii=1:length(tmpnames)
    %     tmpdata.obj.subjlabels{ii} = tmpnames{ii}{2};
    % end 
    
    studyname = 'pnas';
    verbose = false;
    
    if(nargin==0)
        plotOnly = false;
    elseif(nargin>=1)
        plotOnly = varargin{1};
    end
    
	if(exist('brewermap'))	
		colormapfun = @()(flipud(brewermap(length(colormap),'RdYlBu')));
		close all;
	else
		colormapfun = @winter;
	end

    % test re-test mccc 
    exportfun = @(filename)(print('-dpng','-r150',filename));
    fname = ['tmp' filesep datestr(now,'dd-mmm-yyyy')];

    if(~exist(fname,'dir'))
        mkdir(fname)
    end

    nmethods = 2;
    methodnames = {'observed','sn_denoised'};
    

    switch studyname 
    
    case 'pnas'
        datadir = fullfile(getenv('HOME'),'MATLAB','kggm2016-paper','data');
        studydata.cpac = load(fullfile(datadir,'dataset_pnas_rfmri'));                  
        subjlia = cellfun(@(x)(~isempty(x)),studydata.cpac.confounds);
        include_subjects = zeros(1,length(subjlia)); 
        subject_idx = find(subjlia);
        fd_thresh = .2;
        for subjno=1:length(subject_idx)
            subject_no = subject_idx(subjno);
            tmpconfounds = studydata.cpac.confounds{1};
            sessionno1 = ...
                 find(strcmp(studydata.cpac.tms{subject_no},'rest'));
            sessionno2 = ...
                 find(strcmp(studydata.cpac.tms{subject_no},'rest1Hzpre'));
            if(...
            sqrt(mean((...
            tmpconfounds{sessionno1}.framewise_displacement).^2))<fd_thresh ...
            | ...
            sqrt(mean((...
            tmpconfounds{sessionno2}.framewise_displacement).^2))<fd_thresh ...
            )
                include_subjects(subject_no) = 1;
            end
        end
        if(verbose)
            disp('Checking Motion Inclusion')
            sum(subjlia)
            sum(subjlia & include_subjects)
        end
        
        studydata.cpac.confounds = studydata.cpac.confounds(...
                                        subjlia & include_subjects);
        studydata.cpac.X = studydata.cpac.X(subjlia & include_subjects);
        studydata.cpac.subjects = studydata.cpac.subjects(...
                                        subjlia & include_subjects);
        studydata.cpac.files = studydata.cpac.files(subjlia & include_subjects);

        studydata.session1.compcor = [];
        studydata.session2.compcor = [];
        studydata.session1.motion = [];
        studydata.session2.motion = [];
        nsubjects = length(studydata.cpac.X)
        for ii=1:nsubjects
            try
                sessionno = find(strcmp(studydata.cpac.tms{ii},'rest'));
                studydata.session1.X{ii} = studydata.cpac.X{ii}(:,:,sessionno);
                % studydata.session1.compcor(ii,:,:) = studydata.cpac.confounds{ii}{sessionno}.compcor;
                % studydata.session1.motion(ii,:,:) = studydata.cpac.confounds{ii}{sessionno}.motion;
    
                sessionno = find(strcmp(studydata.cpac.tms{ii},'rest1Hzpre'));
                studydata.session2.X{ii} = studydata.cpac.X{ii}(:,:,sessionno);
                % studydata.session2.compcor(ii,:,:) = studydata.cpac.confounds{ii}{sessionno}.compcor;
                % studydata.session2.motion(ii,:,:) = studydata.cpac.confounds{ii}{sessionno}.motion;
            catch me 
                studydata.cpac.subjects{ii}
                disp(me)
            end
        end
    case 'corr-bnu1'
        datadir = fullfile(getenv('HOME'),'MATLAB','rcubed','data');
        studydata.cpac = load(fullfile(datadir,'dataset_corr_rfmri'));
    
        subjlia = cellfun(@(x)(~isempty(x)),studydata.cpac.confounds);
        studydata.cpac.confounds = studydata.cpac.confounds(subjlia);
        studydata.cpac.X = studydata.cpac.X(subjlia);
        studydata.cpac.subjects = studydata.cpac.subjects(subjlia);
        studydata.cpac.files = studydata.cpac.files(subjlia);

        nsubjects = length(studydata.cpac.X)
    
        for ii=1:nsubjects
            ii
            sessionno = find(strcmp(studydata.cpac.tms{ii},'ses-1'))
            studydata.session1.X{ii} = studydata.cpac.X{ii}(:,:,sessionno);        
            sessionno = find(strcmp(studydata.cpac.tms{ii},'ses-2'))
            studydata.session2.X{ii} = studydata.cpac.X{ii}(:,:,sessionno);
        end
    
    end

    communityno = 3;
    ncommunities = length(unique(studydata.cpac.roilabels.Yeo7));
    nodeidx = find(studydata.cpac.roilabels.Yeo7==communityno);
    nodeidx = cat(1,nodeidx,find(studydata.cpac.roilabels.Yeo7==4));
    nodeidx = cat(1,nodeidx,find(studydata.cpac.roilabels.Yeo7==6));
    communitylabel = studydata.cpac.roilabels.Yeo7Label(nodeidx);
    if(length(unique(communitylabel))>2)
        communitylabel = {'AllNetworks'};
    end
    communitycolor = studydata.cpac.roilabels.Yeo7(nodeidx);

	results = struct();
    observed_correlation = [];
    denoised_correlation = []; 
    denoised_correlation2 = []; 
    nuisance_correlation = [];
    nsr = [];

    if(~plotOnly)
        %Y = get_shared_nuisance(squeeze(mean(studydata.signalsSite1,2))');
        Y = [];
        for sessionno=1:2
        	for ii=1:nsubjects
                if(verbose)
                    ii
                    ['session' num2str(sessionno)]
                end
                X = studydata.(['session' num2str(sessionno)]).X{ii};
                Y = mean(X,2);
        		%Y = squeeze(studydata.compcor(ii,:,:));
        		%Y = cat(2,Y,squeeze(studydata.motion(ii,:,:))); 
                try
                    observed_correlation(:,:,ii,sessionno) = ...
                                standard_correlation(X);
                    % results.Site1{ii,sessionno} = conditional_correlation(X,Y);
                    % denoised_correlation(:,:,ii,sessionno) = ...
                    %             results.Site1{ii,sessionno}.corr;
                    % nuisance_correlation(:,:,ii,sessionno) = ...
                    %              results.Site1{ii,sessionno}.nuisance;
                    % nsr(ii,sessionno) = results.Site1{ii,sessionno}.NSR;
                    denoised_correlation2(:,:,ii,sessionno) = ...
                                standard_correlation_sn(X);

                catch me
                    disp(me)
                    results.Site1{ii,sessionno} = {};
                end
                % if(ii==1)
                %     demo_conditional_correlation(X,Y);
                % end
        	end
        end
    
        % results.denoised2 = denoised_correlation;
        % results.nuisance = nuisance_correlation;
        % results.nsr = nsr; 
    
        results.denoised2 = denoised_correlation2;
        results.sn_denoised = denoised_correlation2;
        results.observed = observed_correlation;
        
        sitelabels = ones(1,2*nsubjects);
    
        % if(exist('brewermap'))
        %     colormapfun = @()(flipud(brewermap(length(colormap),'RdYlBu')));
        %     close all;
        % else
        %     colormapfun = @winter;
        % end
        %
        %         % test re-test mccc
        %         exportfun = @(filename)(print('-dpng','-r150',filename));
        %         fname = ['tmp' filesep datestr(now,'dd-mmm-yyyy')];
        %
        %         if(~exist(fname,'dir'))
        %             mkdir(fname)
        %         end
        %
        %         % upper_idx = find(reshape(triu(ones(p,p),1), [p*p 1]));
        %         % features = reshape(results.denoised(upper_idx),[p*p nsubjects 2]);
        %
        %         nmethods = 2;
        %         methodnames = {'observed','sn_denoised'};
    
        for methodno=1:nmethods
            methodname = methodnames{methodno};
            p = size(results.(methodname),1);
            graphs = results.(methodname)(nodeidx,:,:,:);
            results.cpac.(methodname).trt_reliability = ...
                                 helper_node_mccc(graphs);
            results.cpac.(methodname).trt_reliability.nodeidx = nodeidx;
            tmp_kendallW = [];
            for nodeno=1:size(graphs,1)
                for subjectno=1:nsubjects
                    tmp_kendallW(nodeno,subjectno) = ...
                         reliability.kendallsW(...
                             cat(1,squeeze(graphs(nodeno,:,subjectno,1)),...
                                  squeeze(graphs(nodeno,:,subjectno,2)))' ...
                                  );
                end
            end
            results.cpac.(methodname).trt_reliability.kendallsw = tmp_kendallW;
            mean(results.cpac.sn_denoised.trt_reliability.kendallsw-results.cpac.observed.trt_reliability.kendallsw,2);
            save(['tmp/trt_' studyname '_rfmri.mat'],'results');    
        end
    else
        load(['tmp/trt_' studyname '_rfmri.mat'],'results');
    end
    
    % Plotting code only 
    figure;
    set(gcf,'Position',[50 75 1200 600]);  
    for methodno=1:nmethods
        methodname = methodnames{methodno};
        y_methods = results.cpac.(methodname).trt_reliability.mccc;
        y_lo = results.cpac.(methodname).trt_reliability.mccc_ci(:,1);
        y_hi = results.cpac.(methodname).trt_reliability.mccc_ci(:,2);
        
        grammobj(methodno,1) = plot_reliability_dotplot(y_methods, ...
                                        y_lo, ...
                                        y_hi, ...
                                        methodnames,  ...
                                        communitycolor  ...
                                        );
        switch methodname
        case 'sn_denoised'
            grammobj(methodno,1).set_title(['ICC ' 'after centering']);
        case 'observed'
            grammobj(methodno,1).set_title(['ICC ' 'before centering']);
        otherwise
            grammobj(methodno,1).set_title(['ICC for ' methodname]);
        end
    end
    grammobj.draw();
    savefilename = fullfile(fname, ['sessioneffects_reliability' '_' ...
                     'before_after_sn' '_' communitylabel{1}]);
    savefig(savefilename);
    grammobj.export('file_name', ...
        savefilename,'file_type','pdf');
    exportfun(savefilename);
    
    % Plot full mccc node distributions on the same plot
    y_methods = cat(2,...
                results.cpac.(methodnames{1}).trt_reliability.mccc,...
                results.cpac.(methodnames{2}).trt_reliability.mccc ...
                );                
    % y_methods_lo = cat(2,...
    %         results.cpac.(methodnames{1}).trt_reliability.mccc_ci(nodeidx,1),...
    %         results.cpac.(methodnames{2}).trt_reliability.mccc_ci(nodeidx,1) ...
    %         );
    % y_methods_hi = cat(2,...
    %         results.cpac.(methodnames{1}).trt_reliability.mccc_ci(nodeidx,2),...
    %         results.cpac.(methodnames{2}).trt_reliability.mccc_ci(nodeidx,2) ...
    %         );
    % y_groups = {};
    % y_groups{1} = repmat(methodnames(1),[1 size(y_methods,1)]);
    % y_groups{2} = repmat(methodnames(2),[1 size(y_methods,1)]);
    % 
                
    % Plot cornerhistogram
    grammobj = plot_reliability_cornerhist(y_methods, ...
                                    {'Before centering','After centering'}, ...
                                    communitycolor  ...
                                    );
    figure;    
    set(gcf,'Position',[25 75 850 750]);                  
    grammobj.draw();
    savefilename = fullfile(fname, ...
                        ['sessioneffects_reliability' '_' ...
                        'difference' '_' methodnames{1} '_' methodnames{2} ...
                        communitylabel{1} ...
                        ] ...
                    );
    savefig(savefilename);                
    grammobj.export('file_name',savefilename,'file_type','pdf');
    exportfun(savefilename);

    % figure;
    % alternate_color = zeros(length(nodeidx),2);
    % switch color_type
    % case 'odd-even'
    %     alternate_color(1:2:end) = 1;
    %     alternate_color(2:2:end) = 2;
    %     mycolormap = brewermap(length(unique(alternate_color))*...
    %                             ncommunities,...
    %                             'Paired' ...
    %                             );
    % otherwise
    %     alternate_color(:,1) = communityno;
    %     alternate_color(:,2) = 2*communityno;
    %     mycolormap = brewermap(2*ncommunities,'Dark2');
    %     %mycolormap = mycolormap([communityno 2*communityno],:);
    % end
    %
    % grammobj(1,1) =  gramm('x', [1:length(nodeidx)]', ...
    %                     'y', [{y_methods(:,1)}, {y_methods(:,2)}], ...
    %                     'ymin', [{y_methods_lo(:,1)}, {y_methods_lo(:,2)}], ...
    %                     'ymax', [{y_methods_hi(:,1)}, {y_methods_hi(:,2)}], ...
    %                     'label', methodnames, ...
    %                     'color', [{alternate_color(:,1)}, ...
    %                                  {alternate_color(:,2)}], ...
    %                     'marker',{'x','o'} ...
    %                     );
    % switch geom_type
    % case 'geom_bar'
    %     grammobj(1,1).geom_bar('dodge',0.8,'width',0.6);
    % case 'geom_point'
    %     grammobj(1,1).geom_point();
    %     grammobj(1,1).set_point_options('base_size',15);
    % otherwise
    %     grammobj(1,1).stat_density();
    % end
    % % grammobj(1,1).geom_interval('geom','black_errorbar',...
    % %                 'dodge',0.8,'width',1);
    % grammobj(1,1).set_color_options('map', mycolormap);
    % grammobj(1,1).set_names('x', [communitylabel ' ROIs'] ,...
    %                  'y','MCCC',....
    %                  'marker','method');
    % grammobj(1,1).set_text_options('base_size',18);
    % ylim_upper = round(10*max(y_methods(:)))/10;
    % grammobj(1,1).axe_property('ylim',[0 max([.6 ylim_upper])]);
    % grammobj(1,1).no_legend();
    % grammobj(1,1).draw()
    %
    % exportfun(fullfile(fname, ['sessioneffects_reliability' '_' ...
    %                  'density' '_' methodnames{1} '_' methodnames{2} ...
    %                   communitylabel{1}]) ...
    %                  );
    
    
    nboot = 10;
    mccc = zeros(p,1); 
    mccc_ci = zeros(p,2);
    mccc_boot = zeros(p,nboot);
    similarity = zeros(2*nsubjects,2*nsubjects);
    for nodeno=1:length(nodeidx)
        jj = nodeidx(nodeno)
        features1 = ...
                reshape(permute(squeeze(results.(methodnames{1})(jj,:,:,:)),...
                        [1 3 2]),[p*2 nsubjects]);
        features2 = ...
                reshape(permute(squeeze(results.(methodnames{2})(jj,:,:,:)),...
                        [1 3 2]),[p*2 nsubjects]);
        similarity = similarity + cov(cat(2,features1,features2));
        [mccc(jj),~,mccc_ci(jj,:),mccc_boot(jj,:)] = ...
             reliability.mccc(features1(1:p,:),features2(1:p,:));
         mccc(jj)
         mccc_ci(jj,:)
    end
    similarity = similarity./length(nodeidx);
    results.cpac.sn_gm_agreement.mccc = real(mccc);
    results.cpac.sn_gm_agreement.mccc_ci = real(mccc_ci);
    results.cpac.sn_gm_agreement.similarity = similarity;
    results.cpac.sn_gm_agreement.nodeidx = nodeidx;
    
    save(['tmp/trt_' studyname '_rfmri.mat'],'results','-append');


    % Compare consistency of successive normalization to global signal regression
    geom_type = 'geom_point';
    color_type = 'odd-even';
    alternate_color = zeros(1,length(nodeidx)); 
    switch color_type
    case 'odd-even'        
        alternate_color(1:2:end) = 1;
        alternate_color(2:2:end) = 2;
        mycolormap = brewermap(length(unique(alternate_color))*...
                                ncommunities,...
                                'Paired' ...
                                );
    otherwise
        alternate_color(:) = communityno;
        mycolormap = brewermap(ncommunities,'Spectral');
    end
    figure;
    set(gcf,'Position',[50 75 850 600]);
    grammobj = gramm('x',1:length(nodeidx),...
                    'y',results.cpac.sn_gm_agreement.mccc(nodeidx),...
                    'ymin',results.cpac.sn_gm_agreement.mccc_ci(nodeidx,1),...
                    'ymax',results.cpac.sn_gm_agreement.mccc_ci(nodeidx,2),...
                    'color',alternate_color);
    switch geom_type
    case 'geom_bar'
        grammobj.geom_bar('dodge',0.8,'width',0.6);
    case 'geom_point'
        grammobj.geom_point();
        grammobj.set_point_options('base_size',15);
    otherwise
        grammobj.geom_point();
    end
    grammobj.geom_interval('geom','black_errorbar',...
    'dodge',0.8,'width',1); 
    grammobj.set_line_options('base_size',2);
    grammobj.set_color_options('map', mycolormap);
    grammobj.set_names('x', [communitylabel ' ROIs'] ,...
                        'y','MCCC',....
                        'color','ROI');
    grammobj.set_text_options('base_size',18);
    ylim_upper = round(10*max(...
                        results.cpac.sn_gm_agreement.mccc_ci(nodeidx,2)))/10;
    grammobj.axe_property('ylim',[0 max([.6,ylim_upper])]);
    grammobj.no_legend();
    grammobj.draw();
    exportfun(fullfile(fname,['sessioneffects_reliability' ...
         '2_' communitylabel{1} ]));
     
end

function output =  conditional_correlation(X,Y)
	% Only uses usual column standardize (i.e. correlation)
	output = struct();
	
    n1 = size(X,1);
    n2 = size(Y,1); 
    if(n2>n1)
        Y = Y(end-n1+1:end,:);
    elseif(n1>n2)
        X = X(end-n2+1:end,:);
    end
    	
	[Sigma results] = covariance.conditional_sample_covariance_separate(X, ...
									struct('verbose',true,...
											'nuisance',Y) ...	
											);
	 	
	output.corr = results.Sigma;
	output.nuisance = results.nCorr; 
	output.corr2 = covariance.mle_sample_covariance(results.X_perpY, ...
												struct('standardize', 'cols'));
	output.NSR = results.NSR;
	
end


function Sigma = standard_correlation(X)
	% Usual standard correlation matrix
	
	
	[Sigma results] = covariance.mle_sample_covariance(X, ...
												struct('standardize','cols'));
	
	
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


function trt_reliability = helper_node_mccc(graphs)
    % 
    % Inputs
    %   - graphs is p x p x subjects x {rater|session}    
    
    p = size(graphs,1);
    nsubjects = size(graphs,3);
    if(p<=50)
        nodeidx = 1:p;
    end
    
    nboot = 10;
    mccc = zeros(p,1); 
    mccc_ci = zeros(p,2);
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
        jj = nodeno
        features = squeeze(graphs(jj,:,:,:));
        [mccc(jj),~,mccc_ci(jj,:),mccc_boot(jj,:)] = ...
             reliability.mccc(features(:,:,1)',features(:,:,2)');
        mccc(jj)
        mccc_ci(jj,:)
    end
    trt_reliability.mccc = real(mccc); 
    trt_reliability.mccc_ci = real(mccc_ci);
    trt_reliability.mccc_boot = real(mccc_boot);
    trt_reliability.similarity = similarity;
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


function grammobj = plot_reliability_dotplot(y_methods,y_lo,y_hi,methodnames,community)
    
    geom_type = 'geom_point';
    color_type = 'common';
    ncommunities = 7;
    
    alternate_color = zeros(1,length(y_methods)); 
    switch color_type
    case 'odd-even'        
        alternate_color(1:2:end) = 1;
        alternate_color(2:2:end) = 2;
        mycolormap = brewermap(length(unique(alternate_color))*...
                                ncommunities,...
                                'Paired' ...
                                );
    otherwise
        alternate_color = community;
        mycolormap = brewermap(ncommunities,'Spectral');
    end
  
    grammobj = gramm('x',1:length(y_methods),...
                     'y', y_methods,...
                     'ymin', y_lo,...
                     'ymax', y_hi,...
                     'color', alternate_color);
    switch geom_type
    case 'geom_bar'
        grammobj.geom_bar('dodge',0.8,'width',0.6);
    case 'geom_point'
        grammobj.geom_point();
        grammobj.set_point_options('base_size',15);
    otherwise
        grammobj.geom_point();
    end
    grammobj.geom_interval('geom','black_errorbar',...
    'dodge',0.8,'width',1); 
    grammobj.set_color_options('map', mycolormap);
    grammobj.set_names('x', 'ROIs',...
                        'y','MCCC',....
                        'color','ROI');
    grammobj.set_text_options('base_size',18);
    ylim_upper = round(10*max(y_methods(:)))/10;
    grammobj.axe_property('ylim',[0 max([.45 ylim_upper])]);
    grammobj.no_legend();
    
end


function grammobj = plot_reliability_cornerhist(y_methods,methodnames,community)
    % y_methods is features | regions | nodes x methods
    
    ncommunities = 7;
    color_type = 'common';
    geom_type = 'cornerhist';
    [nnodes nmethods] = size(y_methods);
    alternate_color = zeros(1,nnodes); 
    switch color_type
    case 'odd-even'        
        alternate_color(1:2:end) = 1;
        alternate_color(2:2:end) = 2;
        mycolormap = brewermap(length(unique(alternate_color))*...
                                ncommunities,...
                                'Paired' ...
                                );
    otherwise
        alternate_color = community; 
        mycolormap = brewermap(ncommunities,'Spectral');
    end
    grammobj = gramm('x', y_methods(:,1), ...
                     'y', y_methods(:,2), ...
                     'color', alternate_color ...
                     );                 
    grammobj.geom_point();
    grammobj.set_point_options('base_size',12);
    switch geom_type
    case 'qq'
        grammobj.stat_qq();
    case 'cornerhist'
        grammobj.stat_cornerhist('edges',-.3:0.02:.3,'aspect',0.45);
        grammobj.set_color_options('map', mycolormap);
        grammobj.geom_abline();
    end
    grammobj.set_line_options('base_size',3);
    grammobj.set_text_options('base_size',24);
    grammobj.set_title('Difference in ICCs');
    ylim_upper = round(10*max(y_methods(:)))/10;
    grammobj.axe_property('ylim', [0 max([.5 ylim_upper])],...
                          'xlim', [0 max([.5 ylim_upper])] ...
                          );
    grammobj.set_names('x',regexprep(methodnames{1},'_','-'),...
                        'y', regexprep(methodnames{2},'_','-'));
    grammobj.no_legend();     
end