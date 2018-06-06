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
    
    studyname = 'best-eegfmri-cpac';
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
    %fname = ['tmp' filesep datestr(now,'dd-mmm-yyyy')];
    fname = ['tmp' filesep '22-May-2018'];    
    if(~exist(fname,'dir'))
        mkdir(fname)
    end
    fname = [fname filesep studyname]; 
    if(~exist(fname,'dir'))
        mkdir(fname)
    end
    
    nmethods = 1;
    methodnames = {'sn_denoised'};

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
        pipelinename = 'pipeline_global00'; 
        atlasname = 'CC200'; 
        
        fname = [fname filesep pipelinename]; 
        if(~exist(fname,'dir'))
            mkdir(fname)
        end
        
        fname = [fname filesep atlasname]; 
        if(~exist(fname,'dir'))
            mkdir(fname)
        end
    
        studydata.cpac = load(fullfile(datadir,...
                        ['corr-bnu1_task-rest_'  ...
                        pipelinename '_' atlasname ] ...
                        ));
        %subjlia = cellfun(@(x)(~isempty(x)),studydata.cpac.confounds);
        subjlia = studydata.cpac.trt_idx;
        studydata.cpac.confounds = studydata.cpac.confounds(subjlia);
        studydata.cpac.X = studydata.cpac.X(subjlia);
        studydata.cpac.subjects = studydata.cpac.subjects(subjlia);
        studydata.cpac.files = studydata.cpac.subjects;

        nsubjects = length(studydata.cpac.X)
    
        for ii=1:nsubjects
            ii
            sessionno = find(strcmp(studydata.cpac.sessions,'ses-1'));
            studydata.session1.X{ii} = studydata.cpac.X{ii}{sessionno};        
            sessionno = find(strcmp(studydata.cpac.sessions,'ses-2'));
            studydata.session2.X{ii} = studydata.cpac.X{ii}{sessionno};
        end
    
    case 'best-eegfmri-cpac'
        
        pipeline = 'etkinlab_analysis';
        resttype = '';
        n_time = 800;
        fname = [fname '-' resttype 'trt' ['-t' num2str(n_time)] ];
        if(~exist(fname))
            mkdir(fname)
        end
        fname = fullfile(fname,pipeline);
        if(~exist(fname))
            mkdir(fname)
        end
        
        datadir = fullfile(getenv('HOME'),'MATLAB','kggm2016-paper','data');
        studydata = ...
             load(fullfile(datadir,'dataset_SchaeferYeo100_best_rfmri')); 

        
        session1_prefix = ['tp1_efmri_' 'rest2'];
        session1_idx = 1:length(studydata.files);
        studydata.session1.subjects = studydata.files;
        
        session2_prefix = ['tp2_efmri_' 'rest1'];
        session2_idx = 1:length(studydata.files)
        studydata.session2.subjects = studydata.files;
        
        
        all_subnames = ...
              cellfun(@num2str,num2cell(1:14),'UniformOutput',false);
        
        for ii=1:length(all_subnames)

            if(ii<=length(studydata.X) && ...
                 ~isempty(studydata.X{ii}{1}))
                studydata.session1.X{ii} = ...
                 squeeze(studydata.X{ii}{1}(11:810,1:100));
            else
                 studydata.session1.X{ii} = [];
                 studydata.session1.subjects{ii} = [];
            end

            if(ii<=length(studydata.X) && length(studydata.X{ii})>2 && ...
                 ~isempty(studydata.X{ii}{3}))
                 studydata.session2.X{ii} = ...
                 squeeze(studydata.X{ii}{3}(11:810,1:100));
            else
                  studydata.session2.X{ii} = [];
                  studydata.session2.subjects{ii} = [];                  
            end
              
             tmp_subj_no = [];  
        end
        
        both_subjects_idx = ...
                     ~(cellfun(@isempty,studydata.session2.subjects) | ...
                                cellfun(@isempty,studydata.session1.subjects));
        studydata.session1.X = studydata.session1.X(both_subjects_idx);
        studydata.session2.X = studydata.session2.X(both_subjects_idx);
        nsubjects = sum(both_subjects_idx)
        studydata.session1
        studydata.roilabels.table = readtable('/Volumes/MACBOOKUSB/Datasets/Schaefer2018_LocalGlobal/Parcellations/MNI/Schaefer2018_100Parcels_7Networks_order.txt','Delimiter','\t','ReadVariableNames',0);
        studydata.subjects.labels = all_subnames;
        studydata.subjects.index = both_subjects_idx;
        
        save([fname filesep 'data_' pipeline '_' ...
                 studyname '_rfmri.mat'],'studydata');
                
    
    case 'best-eegfmri'
        
        pipeline = 'etkinlab_analysis';
        resttype = '';
        n_time = 400;
        fname = [fname '-' resttype 'tp1' ['-t' num2str(n_time)] ];
        if(~exist(fname))
            mkdir(fname)
        end
        fname = fullfile(fname,pipeline);
        if(~exist(fname))
            mkdir(fname)
        end
        
        datadir = fullfile(getenv('HOME'),'Downloads');
        tmpdata = load(fullfile(datadir,'BEST_Schaefer100PlusSubcortical_BoldSignalsHomAndSNR_w1_05282018.mat'));
        % tmpdata = load(fullfile(datadir,'BEST_Schaefer100PlusSubcortical_BoldSignalsHomAndSNR_w1_03152018.mat'));
        
        subNames = ...
        regexprep(tmpdata.subNames,{'.*best_','.connectivity.*'},{'',''});
        
        session1_prefix = ['tp1_efmri_' 'rest1'];
        session1_idx = ...
         find(~cellfun(@isempty,strfind(subNames,session1_prefix)));
        studydata.session1.subjects = subNames(session1_idx);
        
        session2_prefix = ['tp1_efmri_' 'rest2'];
        session2_idx = ...
            find(~cellfun(@isempty,strfind(subNames,session2_prefix)));
        studydata.session2.subjects = subNames(session2_idx);
        
        all_subnames = ...
              cellfun(@num2str,num2cell(3001:3020),'UniformOutput',false);
        
        for ii=1:length(all_subnames)
            tmp_subj_no = ...
                find(~cellfun(@isempty, ...
                strfind(subNames(session1_idx),all_subnames(ii))));
             if(~isempty(tmp_subj_no))
                studydata.session1.X{ii} = ...
                 squeeze(tmpdata.signals(session1_idx(tmp_subj_no),1:100,11:410))';
                 studydata.session1.subjects{ii} = subNames(session1_idx(tmp_subj_no));
                 
             else
                 studydata.session1.X{ii} = [];
                 studydata.session1.subjects{ii} = [];
                 
             end
             tmp_subj_no = [];     
             
             tmp_subj_no = ...
                 find(~cellfun(@isempty, ...
                 strfind(subNames(session2_idx),all_subnames(ii))));
              if(~isempty(tmp_subj_no))
                 studydata.session2.X{ii} = ...
                  squeeze(tmpdata.signals(session2_idx(tmp_subj_no),1:100,11:410))';
                 studydata.session2.subjects{ii} = ...
                  subNames(session2_idx(tmp_subj_no));
              else
                  studydata.session2.X{ii} = [];
                  studydata.session2.subjects{ii} = [];                  
              end
              
             tmp_subj_no = [];  
        end
        
        both_subjects_idx = ...
                     ~(cellfun(@isempty,studydata.session2.subjects) | ...
                                cellfun(@isempty,studydata.session1.subjects));
        studydata.session1.X = studydata.session1.X(both_subjects_idx);
        studydata.session2.X = studydata.session2.X(both_subjects_idx);
        nsubjects = sum(both_subjects_idx)
        studydata.session1
        studydata.roilabels.table = readtable('/Volumes/MACBOOKUSB/Datasets/Schaefer2018_LocalGlobal/Parcellations/MNI/Schaefer2018_100Parcels_7Networks_order.txt','Delimiter','\t','ReadVariableNames',0);
        studydata.subjects.labels = all_subnames;
        studydata.subjects.index = both_subjects_idx;
        
        save([fname filesep 'data_' pipeline '_' ...
                 studyname '_rfmri.mat'],'studydata');
        
    end
    

    if(isfield(studydata,'cpac'))
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
    else
        nodeidx = 1:size(studydata.session1.X{1},2);
        tmp_community_parts = ...
             regexp(studydata.roilabels.table.Var2,{'_'}, 'split');
        for ii=1:length(tmp_community_parts)     
            communitylabels{ii} = tmp_community_parts{ii}{3};
        end
        communitylabel = grp2idx(communitylabels);
        communitycolor = communitylabel;
    end
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
                X = X(1:n_time,:);
                %Y = mean(X,2);
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
        
        if(isfield(studydata,'subjects'))
            subject_labels = ...
                 studydata.subjects.labels(studydata.subjects.index)
        end
        % results.denoised2 = denoised_correlation;
        % results.nuisance = nuisance_correlation;
        % results.nsr = nsr; 
    
        results.denoised2 = denoised_correlation2;
        results.sn_denoised = denoised_correlation2;
        results.observed = observed_correlation;
        if(exist('subject_labels'))
            results.subject_labels = subject_labels;
        end
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
            save([fname filesep 'trt_' studyname '_rfmri.mat'],'results');

            tmp_kendallW = [];
            if(~exist(fullfile(fname,methodname)))
                mkdir(fullfile(fname,methodname));
            end
            
            ba_nodes = nodeidx;
            trid_idx = find(reshape(triu(ones(p,p),1),[p^2 1]));
            all_y_diff = [];
            all_x_mean = [];
            tmp_propLOA = [];
            for subjectno=1:nsubjects                
         %  [altmanplot, kendallcoef, propLOA, y_diff, x_mean] = ...
         %  plot_bland_altman_agreement(graphs(ba_nodes, ba_nodes,subjectno,1),...
         %  graphs(ba_nodes, ba_nodes,subjectno,2));
                rater1 = graphs(ba_nodes,ba_nodes,subjectno,1); 
                rater2 = graphs(ba_nodes,ba_nodes,subjectno,2);
                rater1 = rater1(trid_idx);
                rater2 = rater2(trid_idx);
                baStats = reliability.bland_altman_MD(rater1,rater2); 
                y_diff = baStats.D;
                x_mean = baStats.M;
                propLOA = nan;
                kendallcoef = reliability.kendallsW(cat(2,rater1,rater2));
                % savefilename = fullfile(fname, methodname, ...
                %         ['agreement_subj' num2str(subjectno)]);
                % altmanplot.export('file_name', ...
                %    savefilename,'file_type','pdf');
                tmp_kendallW(subjectno) = kendallcoef;
                tmp_propLOA(subjectno) = propLOA;
                all_y_diff(:,subjectno) = y_diff;
                all_x_mean(:,subjectno) = x_mean;
                close;
            end
            
            loa_matrix = zeros(p^2,3);
            for ii=1:size(all_y_diff,1)
                bastats = reliability.bland_altman_loa(...
                                    all_x_mean(ii,:)', ...
                                    all_y_diff(ii,:)' ...
                                    );
                loa_matrix(trid_idx(ii),1) = bastats.stats.mu;
                loa_matrix(trid_idx(ii),2) = bastats.stats.ci(1);                               loa_matrix(trid_idx(ii),3) = bastats.stats.ci(2);    
            end
            figure;
            %scatter(mean(all_x_mean,2),mean(all_y_diff,2))
            loa_int = [];
            loa_int(1,1) = mean(loa_matrix(trid_idx,2)); 
            loa_int(1,2) = mean(loa_matrix(trid_idx,3)); 
            loa_int(2:3,1) = prctile(loa_matrix(trid_idx,2),[2.5 97.5]);
            loa_int(2:3,2) = prctile(loa_matrix(trid_idx,3),[2.5 97.5]);
            
            [altmanplot tmp_propLOA] = ...
                             plot_bland_altman_agreement(mean(all_x_mean,2), ...
                                        mean(all_y_diff,2),  ...
                                        loa_int ...
                                        );
            savefilename = fullfile(fname, methodname, ...
                         ['group_edgelevel_mdplot']);
            altmanplot.export('file_name', ...
               savefilename,'file_type','pdf');
            results.cpac.(methodname).trt_reliability.edgeidx = trid_idx;
            results.cpac.(methodname).trt_reliability.y_diff = all_y_diff;
            results.cpac.(methodname).trt_reliability.x_mean = all_x_mean;
            results.cpac.(methodname).trt_reliability.loa_matrix = ...
                                         sparse(loa_matrix); 
            results.cpac.(methodname).trt_reliability.propLOA = tmp_propLOA;
            results.cpac.(methodname).trt_reliability.global_kendallsw = tmp_kendallW;
            save([fname filesep 'trt_' studyname '_rfmri.mat'],'results','-append');
        end
    
        try
            mean(results.cpac.sn_denoised.trt_reliability.kendallsw-...
                results.cpac.observed.trt_reliability.kendallsw,2);
        catch
            
        end
        

        % Compare Agreement With and Without Successive Normalization
        % nboot = 10;
        % mccc = zeros(p,1);
        % mccc_ci = zeros(p,2);
        % mccc_boot = zeros(p,nboot);
        % similarity = zeros(2*nsubjects,2*nsubjects);
        % for nodeno=1:length(nodeidx)
        %     jj = nodeidx(nodeno)
        %     features1 = ...
        %             reshape(permute(squeeze(results.(methodnames{1})(jj,:,:,:)),...
        %                     [1 3 2]),[p*2 nsubjects]);
        %     features2 = ...
        %             reshape(permute(squeeze(results.(methodnames{2})(jj,:,:,:)),...
        %                     [1 3 2]),[p*2 nsubjects]);
        %     similarity = similarity + cov(cat(2,features1,features2));
        %     [mccc(jj),~,mccc_ci(jj,:),mccc_boot(jj,:)] = ...
        %          reliability.mccc(features1(1:p,:),features2(1:p,:));
        %      mccc(jj)
        %      mccc_ci(jj,:)
        % end
        % similarity = similarity./length(nodeidx);
        % results.cpac.sn_gm_agreement.mccc = real(mccc);
        % results.cpac.sn_gm_agreement.mccc_ci = real(mccc_ci);
        % results.cpac.sn_gm_agreement.similarity = similarity;
        % results.cpac.sn_gm_agreement.nodeidx = nodeidx;
        %
        % results = setfield(results,'sn_denoised',[]);
        % results = setfield(results,'observed',[]);
        % results = setfield(results,'denoised2',[]);
        % save([fname filesep 'trt_' studyname '_rfmri.mat'],'results','-append');

    else
        load([fname filesep 'trt_' studyname '_rfmri.mat'],'results');
        
        if(~isempty(results.observed))
            results = setfield(results,'sn_denoised',[]);
            results = setfield(results,'observed',[]);
            results = setfield(results,'denoised2',[]);
            save([fname filesep 'trt_' studyname ...
                                     '_rfmri.mat'],'results','-append');
        end
        
        if(isfield(results,'subject_labels'))
            subject_labels = ...
                 results.subject_labels;
        end
        
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
    if(strcmp(studyname,'best-eegfmri')|| strcmp(studyname,'best-eegfmri-cpac'))
        savefilename = fullfile(fname, ['sessioneffects_reliability' '_' ...
                         'before_after_sn']);
    else
        savefilename = fullfile(fname, ['sessioneffects_reliability' '_' ...
                     'before_after_sn' '_' communitylabel{1}]);
    end
    savefig(savefilename);
    grammobj.export('file_name', ...
        savefilename,'file_type','pdf');
    exportfun(savefilename);
    
    figure; 
    addpath('~/MATLAB/packages/Violinplot-Matlab'); 
    set(gcf,'Position',[25 75 400 600]); 
    vplotobj = violinplot(results.cpac.('sn_denoised').trt_reliability.mccc); 
    title('ROI ICC Distribution')
    set(gcf,'Color',[.9 .9 .9]); 
    set(gcf,'Color',[.9 .9 .9]);
    set(gca,'fontsize',20)
    savefilename = fullfile(fname, ['MCCC_Summary.png']);
    export_fig(savefilename,'-transparent');
    pause(3)
    
    figure; 
    addpath('~/MATLAB/packages/Violinplot-Matlab'); 
    set(gcf,'Position',[25 75 400 600]); 
    vplotobj = violinplot(results.cpac.('sn_denoised').trt_reliability.global_kendallsw); 
    title('Subjects Kendall (W) Distribution')
    set(gcf,'Color',[.9 .9 .9]); 
    set(gcf,'Color',[.9 .9 .9]);
    set(gca,'fontsize',20)
    savefilename = fullfile(fname, ['KendallW_Summary.png']);
    export_fig(savefilename,'-transparent');
    pause(3)
    
    
    % Plot full mccc node distributions on the same plot
    for methodno=1:nmethods
        if(nmethods>1)
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
        
        if(strcmp(studyname,'best-eegfmri'))
            savefilename = fullfile(fname, ['sessioneffects_reliability' '_' ...
                        'difference' '_' methodnames{1} '_' methodnames{2}] ...
                        );
        else
            savefilename = fullfile(fname, ...
                        ['sessioneffects_reliability' '_' ...
                        'difference' '_' methodnames{1} '_' methodnames{2} ...
                        communitylabel{1} ...
                        ] ...
                        );
        end
        

        savefig(savefilename);                
        grammobj.export('file_name',savefilename,'file_type','pdf');
        exportfun(savefilename);

        y_methods = cat(2,...
             mean(results.cpac.(methodnames{1}).trt_reliability.kendallsw,2),...
             mean(results.cpac.(methodnames{2}).trt_reliability.kendallsw,2) ...
            );
            % Plot cornerhistogram
            grammobj = plot_reliability_cornerhist(y_methods, ...
                                    {'Before centering','After centering'}, ...
                                    communitycolor  ...
                                    );
            figure;
            set(gcf,'Position',[25 75 850 750]);
            grammobj.draw();
            savefilename = fullfile(fname, ...
                        ['sessioneffects_concordance' '_' ...
                        'difference' '_' methodnames{1} '_' methodnames{2} ...
                        communitylabel{1} ...
                        ] ...
                    );
            savefig(savefilename);
            grammobj.export('file_name',savefilename,'file_type','pdf');
            exportfun(savefilename);

        else
            y_methods = ...
                 results.cpac.(methodnames{1}).trt_reliability.mccc;
        end
    end


    addpath('~/MATLAB/packages/Violinplot-Matlab'); 
    figure;
    set(gcf,'Position',[25 75 1200 650]);
    %vplotobj = violinplot(results.cpac.(methodname).trt_reliability.kendallsw');
    vplotobj = plot_violinplot(...
            results.cpac.(methodname).trt_reliability.kendallsw', ...
            {methodname}, ...
            communitycolor ...
            );
    if(iscell(communitylabel))
        communitylabelname = communitylabel{1}
    else
        communitylabelname = ''
    end
    savefilename = fullfile(fname, ...
                        ['sessioneffects_violinplot' '_' ...
                        methodname ...
                        communitylabelname ...
                        ] ...
                    );
    savefig(savefilename);
    exportfun(savefilename);   



    % Export kendallsW with powers parameters into a single table. 
    reliability_table = table();
    methodname = 'sn_denoised';
    
    for subjno=1:nsubjects
        
        if(isfield(studydata,'cpac') && isfield(studydata.cpac,'confounds'))
        
            power_table = studydata.cpac.confounds{subjno}{1}.powers_params(:,...
                                            [1 3 6]);     
            tmp_table = power_table;
            tmp_table.Properties.VariableNames(2:3) = ...
                 strcat(tmp_table.Properties.VariableNames(2:3),{'_ses1'});
             
        
            power_table = studydata.cpac.confounds{subjno}{2}.powers_params(:,...
                                                 [1 3 6]); 
            tmp_table2 = power_table;
            tmp_table2.Properties.VariableNames(2:3) = ...
                 strcat(tmp_table2.Properties.VariableNames(2:3),{'_ses2'});
        
            tmp_table = horzcat(tmp_table,tmp_table2(:,2:3));
            
            flag_reliability_vs_motion = 1;
        else
            tmp_table = table();
            try
                tmp_table.subject = subjno;
                tmp_table.subjectlabel = subject_labels(subjno);
            catch me
                disp(me)
            end
        end
        tmp_table.kendallw = ...
         [results.cpac.(methodname).trt_reliability.kendallsw(:,subjno)'];
        reliability_table = vertcat(reliability_table,tmp_table);
        
        flag_reliability_vs_motion = 0;
        
    end
    
    save([fname filesep 'table_kendallw_' methodname ],'reliability_table');
    writetable(reliability_table, [fname filesep 'table_kendallw_' methodname '.csv'],'Delimiter','\t');
    

    if(flag_reliability_vs_motion)
        % Scatterplot    
        communityno = 3;
        y_methods = cat(2,(reliability_table.rootMeanSquareFD_ses1+ ...
                            reliability_table.rootMeanSquareFD_ses2)/2, ...
                            mean(reliability_table.kendallw(:, ...
                                   communitycolor==communityno),2) ...
                        );

        grammobj = plot_reliability_glm(y_methods, ...
                                        {'RMSE FD','Kendall W'}, ...
                                        communityno  ...
                                        );
        figure;
        set(gcf,'Position',[25 75 850 750]);
        grammobj.draw();                                    
        % scatter( (reliability_table.rootMeanSquareFD_ses1+ ...
        %             reliability_table.rootMeanSquareFD_ses2)/2,  ...
        %             mean(reliability_table.kendallw,2), 50);
        savefilename = fullfile(fname, ...
                            ['sessioneffects_kendallw' '_' ...
                            methodname ...
                            communitylabel{1} ...
                            ] ...
                        );    
        savefig(savefilename);
        grammobj.export('file_name',savefilename,'file_type','pdf');
        exportfun(savefilename); 
    end
    % % Compare consistency of successive normalization to global signal regression
    %
    % y_methods = results.cpac.sn_gm_agreement.mccc(nodeidx);
    % y_lo = results.cpac.sn_gm_agreement.mccc_ci(nodeidx,1);
    % y_hi = results.cpac.sn_gm_agreement.mccc_ci(nodeidx,2);
    % figure;
    % set(gcf,'Position',[50 75 1200 600]);
    % grammobj = plot_reliability_dotplot(y_methods, ...
    %                                 y_lo, ...
    %                                 y_hi, ...
    %                                 'Before/After SN Agreement',  ...
    %                                 communitycolor  ...
    %                                 );
    % grammobj.set_title(['ICC Method Agreement']);
    % grammobj.draw();
    % exportfun(fullfile(fname,['sessioneffects_reliability' ...
    %      '2_' communitylabel{1} ]));
    %
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
        jj = nodeno
        features = squeeze(graphs(jj,:,:,:));
        [mccc(jj), MCCC, mccc_ci(jj,:), mccc_boot(jj,:)] = ...
             reliability.mccc(features(:,:,1)',features(:,:,2)');
        mccc_var(jj) = MCCC.normVind;
        mccc_wvar(jj) = MCCC.normVdep;
        disp(sprintf('MCCC: %2.4f, CI: (%2.4f,%2.4f), WithinVar: %.4f, TotalVar: %.4f', ...
        mccc(jj), mccc_ci(jj,1), mccc_ci(jj,2), mccc_wvar(jj), mccc_var(jj)))
    end
    trt_reliability.mccc = real(mccc); 
    trt_reliability.mccc_tvar = mccc_var;
    trt_reliability.mccc_wvar = mccc_wvar;
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


function vplotobj = plot_violinplot(datamatrix,methodnames,community)
    
    
    addpath('~/MATLAB/packages/Violinplot-Matlab'); 
    
    
    ncommunities = 7;
    vcolormap = brewermap(ncommunities,'Spectral');
    
    vplotobj = violinplot(datamatrix, ...
                'ShowNotches');
    
    
    for plotno=1:length(vplotobj)
        set(vplotobj(plotno).ViolinPlot, ...
            'FaceColor',vcolormap(community(plotno),:) ...
            );
        set(vplotobj(plotno).ScatterPlot, ...
            'MarkerEdgeColor', vcolormap(community(plotno),:) ...
            );
        set(vplotobj(plotno).ScatterPlot, ...
            'CData', vcolormap(community(plotno),:) ...
            );
        set(vplotobj(plotno).ScatterPlot, ...
            'MarkerFaceColor', vcolormap(community(plotno),:) ...
            );
        set(vplotobj(plotno).MeanPlot, ...
            'Color', vcolormap(community(plotno),:) ...
            );
    end  
    xlabel('ROIs','fontsize',20,'fontweight','bold'); 
    ylabel('Reliability', 'fontsize',20,'fontweight','bold');
    title(regexprep(methodnames{1},'_','-'));
    xtickvals = get(gca,'XTick');
    xticklabels = get(gca,'XTickLabels');         
    set(gca,'FontSize',18, ...
            'XTick',xtickvals(1:4:end), ...
            'XTickLabels',xticklabels(1:4:end) ...
            )
            %'XTickLabelRotation',90)
    
    
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
    grammobj.axe_property('ylim',[0 max([.6 ylim_upper])]);
    grammobj.no_legend();
    
end


function grammobj = plot_reliability_glm(y_methods,methodnames,community)
    % y_methods is features | regions | nodes x methods
    
    ncommunities = 7;
    color_type = 'common';
    geom_type = 'glm';
    [nsubjects nmethods] = size(y_methods);
    alternate_color = zeros(1,nsubjects); 
    switch color_type
    case 'odd-even'        
        alternate_color(1:2:end) = 1;
        alternate_color(2:2:end) = 2;
        mycolormap = brewermap(length(unique(alternate_color))*...
                                ncommunities,...
                                'Paired' ...
                                );
    otherwise
        alternate_color(:) = community;
        mycolormap = brewermap(ncommunities,'Spectral');
    end
    ylim_upper = round(10*max(y_methods(:)))/10;
    ylim_lower = round(10*min(y_methods(:)))/10;
    if(ylim_lower<.2)
        ylim_lower = 0;
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
    case 'glm'
        grammobj.stat_glm();
    case 'cornerhist'
        if(ylim_upper>.8)
            edgelim = -.1:.005:.1;
        else
            edgelim = -.3:.01:.3;
        end
        grammobj.stat_cornerhist('edges', edgelim,'aspect',.5);
        %grammobj.stat_cornerhist('edges',-.3:01:.3,'aspect',0.6);
    end
    grammobj.set_color_options('map', mycolormap(unique(community),:));
    grammobj.set_line_options('base_size',3);
    grammobj.set_text_options('base_size',24);
    grammobj.set_title('QC covariates vs. Concordance');
    % grammobj.axe_property('ylim', [min([ylim_lower,.6]) ...
    %                                 max([.75 ylim_upper])],...
    %                       'xlim', [min([ylim_lower,.6]) ...
    %                                 max([.75 ylim_upper])] ...
    %                       );
    grammobj.set_names('x',regexprep(methodnames{1},'_','-'),...
                        'y', regexprep(methodnames{2},'_','-'));
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
    ylim_upper = round(10*max(y_methods(:)))/10;
    ylim_lower = round(10*min(y_methods(:)))/10;
    if(ylim_lower<.2)
        ylim_lower = 0;
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
    case 'glm'
        grammobj.stat_glm();
    case 'cornerhist'
        if(ylim_upper>.8)
            edgelim = -.1:.005:.1;
        else
            edgelim = -.3:.01:.3;
        end
        grammobj.stat_cornerhist('edges', edgelim,'aspect',.5);
        %grammobj.stat_cornerhist('edges',-.3:01:.3,'aspect',0.6);
    end
    grammobj.set_color_options('map', mycolormap);
    grammobj.geom_abline();
    grammobj.set_line_options('base_size',3);
    grammobj.set_text_options('base_size',24);
    grammobj.set_title('Difference in ICCs');
    grammobj.axe_property('ylim', [min([ylim_lower,.6]) ...
                                    max([.75 ylim_upper])],...
                          'xlim', [min([ylim_lower,.6]) ...
                                    max([.75 ylim_upper])] ...
                          );
    grammobj.set_names('x',regexprep(methodnames{1},'_','-'),...
                        'y', regexprep(methodnames{2},'_','-'));
    grammobj.no_legend();     
end