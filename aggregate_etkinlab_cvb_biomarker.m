function [roi_table results] = aggregate_etkinlab_cvb_biomarker()
    
    roi_table = table();
    
    methodname = 'sn_denoised';
    metricname = 'mccc';
    studyname = 'best-eegfmri';
    resttype = 'rest2';
    sessiontype = 'tp1';
    sessionlabel = 'within';
    runlabel = 'within';
    samplesize = 't400';
    samplesizename = '5min';
    results_name = [studyname '-' resttype sessiontype '-' samplesize];
    disp(['Results: ' results_name])
    
    pipeline_name = 'etkinlab_analysis';
    
    resultsdir = fullfile('tmp','15-Oct-2018',results_name,pipeline_name);
    ls(resultsdir)
    load(fullfile(resultsdir,'trt_best-eegfmri_rfmri.mat'))
    
    metric = results.cpac.(methodname).trt_reliability.(metricname);
    metric = round(metric,4);
    sm_roi_table = table();
    sm_roi_table.(metricname) = reshape(metric,[length(metric) 1]);
    sm_roi_table.restName = repmat({resttype},[length(metric) 1]);
    sm_roi_table.sessionType = repmat({sessionlabel},[length(metric) 1]);
    sm_roi_table.runType = repmat({runlabel},[length(metric) 1]);
    sm_roi_table.n_time = repmat({samplesize},[length(metric) 1]);
    sm_roi_table.method = repmat({methodname},[length(metric) 1]);
    roi_table = vertcat(roi_table,sm_roi_table);
    
    
    resttype = 'rest1';
    sessiontype = 'tp2';
    sessionlabel = 'within';
    runlabel = 'within';
    samplesize = 't400';
    samplesizename = '5min';
    results_name = [studyname '-' resttype sessiontype '-' samplesize];
    disp(['Results: ' results_name])
    
    pipeline_name = 'etkinlab_analysis';
    
    resultsdir = fullfile('tmp','22-Mar-2018',results_name,pipeline_name);
    ls(resultsdir)
    load(fullfile(resultsdir,'trt_best-eegfmri_rfmri.mat'))
    
    metric = results.cpac.(methodname).trt_reliability.(metricname);
    metric = round(metric,4);
    sm_roi_table = table();
    sm_roi_table.(metricname) = reshape(metric,[length(metric) 1]);
    sm_roi_table.restName = repmat({resttype},[length(metric) 1]);
    sm_roi_table.sessionType = repmat({sessionlabel},[length(metric) 1]);
    sm_roi_table.runType = repmat({runlabel},[length(metric) 1]);
    sm_roi_table.n_time = repmat({samplesize},[length(metric) 1]);
    sm_roi_table.method = repmat({methodname},[length(metric) 1]);
    roi_table = vertcat(roi_table,sm_roi_table);
    
    
    resttype = 'rest2';
    sessiontype = 'tp1';
    sessionlabel = 'within';
    runlabel = 'within';
    samplesize = 't200';
    samplesizename = '2min';
    results_name = [studyname '-' resttype sessiontype '-' samplesize];
    disp(['Results: ' results_name])
    
    pipeline_name = 'etkinlab_analysis';
    
    resultsdir = fullfile('tmp','22-Mar-2018',results_name,pipeline_name);
    ls(resultsdir)
    load(fullfile(resultsdir,'trt_best-eegfmri_rfmri.mat'))
    
    metric = results.cpac.(methodname).trt_reliability.(metricname);
    metric = round(metric,4);
    sm_roi_table = table();
    sm_roi_table.(metricname) = reshape(metric,[length(metric) 1]);
    sm_roi_table.restName = repmat({resttype},[length(metric) 1]);
    sm_roi_table.sessionType = repmat({sessionlabel},[length(metric) 1]);
    sm_roi_table.runType = repmat({runlabel},[length(metric) 1]);
    sm_roi_table.n_time = repmat({samplesize},[length(metric) 1]);
    sm_roi_table.method = repmat({methodname},[length(metric) 1]);
    roi_table = vertcat(roi_table,sm_roi_table);
    
    resttype = 'rest1';
    sessiontype = 'tp2';
    sessionlabel = 'within';
    runlabel = 'within';
    samplesize = 't200';
    samplesizename = '5min';
    results_name = [studyname '-' resttype sessiontype '-' samplesize];
    disp(['Results: ' results_name])
    
    pipeline_name = 'etkinlab_analysis';
    
    resultsdir = fullfile('tmp','22-Mar-2018',results_name,pipeline_name);
    ls(resultsdir)
    load(fullfile(resultsdir,'trt_best-eegfmri_rfmri.mat'))
    
    metric = results.cpac.(methodname).trt_reliability.(metricname);
    metric = round(metric,4);
    sm_roi_table = table();
    sm_roi_table.(metricname) = reshape(metric,[length(metric) 1]);
    sm_roi_table.restName = repmat({resttype},[length(metric) 1]);
    sm_roi_table.sessionType = repmat({sessionlabel},[length(metric) 1]);
    sm_roi_table.runType = repmat({runlabel},[length(metric) 1]);
    sm_roi_table.n_time = repmat({samplesize},[length(metric) 1]);
    sm_roi_table.method = repmat({methodname},[length(metric) 1]);
    roi_table = vertcat(roi_table,sm_roi_table);

    resttype = '';
    sessiontype = 'tp1';
    sessionlabel = 'within';
    runlabel = 'between';
    samplesize = 't200';
    samplesizename = '2min';
    results_name = [studyname '-' resttype sessiontype '-' samplesize];
    disp(['Results: ' results_name])
    
    pipeline_name = 'etkinlab_analysis';
    
    resultsdir = fullfile('tmp','22-Mar-2018',results_name,pipeline_name);
    ls(resultsdir)
    load(fullfile(resultsdir,'trt_best-eegfmri_rfmri.mat'))
    
    metric = results.cpac.(methodname).trt_reliability.(metricname);
    metric = round(metric,4);
    sm_roi_table = table();
    sm_roi_table.(metricname) = reshape(metric,[length(metric) 1]);
    sm_roi_table.restName = repmat({resttype},[length(metric) 1]);
    sm_roi_table.sessionType = repmat({sessionlabel},[length(metric) 1]);
    sm_roi_table.runType = repmat({runlabel},[length(metric) 1]);
    sm_roi_table.n_time = repmat({samplesize},[length(metric) 1]);
    sm_roi_table.method = repmat({methodname},[length(metric) 1]);
    roi_table = vertcat(roi_table,sm_roi_table);


    resttype = '';
    sessiontype = 'tp1';
    sessionlabel = 'within';
    runlabel = 'between';
    samplesize = 't400';
    samplesizename = '5min';
    results_name = [studyname '-' resttype sessiontype '-' samplesize];
    disp(['Results: ' results_name])
    
    pipeline_name = 'etkinlab_analysis';
    
    resultsdir = fullfile('tmp','22-Mar-2018',results_name,pipeline_name);
    ls(resultsdir)
    load(fullfile(resultsdir,'trt_best-eegfmri_rfmri.mat'))
    
    metric = results.cpac.(methodname).trt_reliability.(metricname);
    metric = round(metric,4);
    sm_roi_table = table();
    sm_roi_table.(metricname) = reshape(metric,[length(metric) 1]);
    sm_roi_table.restName = repmat({resttype},[length(metric) 1]);
    sm_roi_table.sessionType = repmat({sessionlabel},[length(metric) 1]);
    sm_roi_table.runType = repmat({runlabel},[length(metric) 1]);
    sm_roi_table.n_time = repmat({samplesize},[length(metric) 1]);
    sm_roi_table.method = repmat({methodname},[length(metric) 1]);
    roi_table = vertcat(roi_table,sm_roi_table);
    
    resttype = '';
    sessiontype = 'tp1';
    sessionlabel = 'within';
    runlabel = 'between';
    samplesize = 't800';
    samplesizename = '10min';
    results_name = [studyname '-' resttype sessiontype '-' samplesize];
    disp(['Results: ' results_name])
    
    pipeline_name = 'etkinlab_analysis';
    
    resultsdir = fullfile('tmp','22-Mar-2018',results_name,pipeline_name);
    ls(resultsdir)
    load(fullfile(resultsdir,'trt_best-eegfmri_rfmri.mat'))

    metric = results.cpac.(methodname).trt_reliability.(metricname);
    metric = round(metric,4);
    sm_roi_table = table();
    sm_roi_table.(metricname) = reshape(metric,[length(metric) 1]);
    sm_roi_table.restName = repmat({resttype},[length(metric) 1]);
    sm_roi_table.sessionType = repmat({sessionlabel},[length(metric) 1]);
    sm_roi_table.runType = repmat({runlabel},[length(metric) 1]);
    sm_roi_table.n_time = repmat({samplesize},[length(metric) 1]);
    sm_roi_table.method = repmat({methodname},[length(metric) 1]);
    roi_table = vertcat(roi_table,sm_roi_table);

    resttype = '';
    sessiontype = '';
    sessionlabel = 'between';
    runlabel = 'between';
    samplesize = 't200';
    samplesizename = '2min';
    results_name = [studyname '-' resttype sessiontype '-' samplesize];
    disp(['Results: ' results_name])
    
    pipeline_name = 'etkinlab_analysis';
    
    resultsdir = fullfile('tmp','22-Mar-2018',results_name,pipeline_name);
    ls(resultsdir)
    load(fullfile(resultsdir,'trt_best-eegfmri_rfmri.mat'))
    
    metric = results.cpac.(methodname).trt_reliability.(metricname);
    metric = round(metric,4);
    sm_roi_table = table();
    sm_roi_table.(metricname) = reshape(metric,[length(metric) 1]);
    sm_roi_table.restName = repmat({resttype},[length(metric) 1]);
    sm_roi_table.sessionType = repmat({sessionlabel},[length(metric) 1]);
    sm_roi_table.runType = repmat({runlabel},[length(metric) 1]);
    sm_roi_table.n_time = repmat({samplesize},[length(metric) 1]);
    sm_roi_table.method = repmat({methodname},[length(metric) 1]);
    roi_table = vertcat(roi_table,sm_roi_table);
    
    resttype = '';
    sessiontype = '';
    sessionlabel = 'between';
    runlabel = 'between';
    samplesize = 't400';
    samplesizename = '5min';
    results_name = [studyname '-' resttype sessiontype '-' samplesize];
    disp(['Results: ' results_name])
    
    pipeline_name = 'etkinlab_analysis';
    
    resultsdir = fullfile('tmp','22-Mar-2018',results_name,pipeline_name);
    ls(resultsdir)
    load(fullfile(resultsdir,'trt_best-eegfmri_rfmri.mat'))
    
    metric = results.cpac.(methodname).trt_reliability.(metricname);
    metric = round(metric,4);
    sm_roi_table = table();
    sm_roi_table.(metricname) = reshape(metric,[length(metric) 1]);
    sm_roi_table.restName = repmat({resttype},[length(metric) 1]);
    sm_roi_table.sessionType = repmat({sessionlabel},[length(metric) 1]);
    sm_roi_table.runType = repmat({runlabel},[length(metric) 1]);
    sm_roi_table.n_time = repmat({samplesize},[length(metric) 1]);
    sm_roi_table.method = repmat({methodname},[length(metric) 1]);
    roi_table = vertcat(roi_table,sm_roi_table);
    
    resttype = '';
    sessiontype = '';
    sessionlabel = 'between';
    runlabel = 'between';
    samplesize = 't800';
    samplesizename = '10min';
    results_name = [studyname '-' resttype sessiontype '-' samplesize];
    disp(['Results: ' results_name])
    
    pipeline_name = 'etkinlab_analysis';
    
    resultsdir = fullfile('tmp','22-Mar-2018',results_name,pipeline_name);
    ls(resultsdir)
    load(fullfile(resultsdir,'trt_best-eegfmri_rfmri.mat'))
    
    metric = results.cpac.(methodname).trt_reliability.(metricname);
    metric = round(metric,4);
    sm_roi_table = table();
    sm_roi_table.(metricname) = reshape(metric,[length(metric) 1]);
    sm_roi_table.restName = repmat({resttype},[length(metric) 1]);
    sm_roi_table.sessionType = repmat({sessionlabel},[length(metric) 1]);
    sm_roi_table.runType = repmat({runlabel},[length(metric) 1]);
    sm_roi_table.n_time = repmat({samplesize},[length(metric) 1]);
    sm_roi_table.method = repmat({methodname},[length(metric) 1]);
    roi_table = vertcat(roi_table,sm_roi_table);
    
    savefilename = ['tmp/22-Mar-2018/all_' metricname ];
    writetable(roi_table,savefilename,'WriteVariableNames',1,'Delimiter','\t');
    
    within_run_idx = strcmp(roi_table.sessionType,'within') & ...
             strcmp(roi_table.runType,'within');
    t400_run_idx = strcmp(roi_table.n_time,'t200');
    y_methods = roi_table(within_run_idx,:).(metricname)
    x_methods = roi_table(within_run_idx,:).n_time;
%    x_methods = roi_table(within_run_idx,:).runType;
    [y_mean y_ci] = grpstats(y_methods,x_methods,{'mean','meanci'})
    grammobj(1,1) = plot_reliability_dotplot(unique(x_methods),fliplr(y_mean'),...
                 fliplr(y_ci(:,1)'),fliplr(y_ci(:,2)'),(unique(x_methods)),[]);
    
     within_run_idx = strcmp(roi_table.sessionType,'between') & ...
              strcmp(roi_table.runType,'between');
     y_methods = roi_table(within_run_idx,:).(metricname)
     x_methods = roi_table(within_run_idx,:).n_time;
    %    x_methods = roi_table(within_run_idx,:).runType;
     [y_mean y_ci] = grpstats(y_methods,x_methods,{'mean','meanci'})
     grammobj(2,1) = plot_reliability_dotplot(unique(x_methods),(y_mean'),...
                 (y_ci(:,1)'),(y_ci(:,2)'),(unique(x_methods)),[]);
    grammobj.draw()
    resultsdir
    print('-dpng','-r150', ...
            fullfile('tmp','22-Mar-2018','all_mccc_t200_t400_mccc.png'));   %grammobj.export('file_name','tmp/22-Mar-2018/between_within_t400_mccc.png');
    
end


function export_mean_roi_reliability(roi_table)
    
    
    
    
end


function grammobj = plot_reliability_dotplot(x_methods,  y_methods,y_lo,y_hi,methodnames,community)
    
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
  
    grammobj = gramm('x', x_methods,...
                     'y', y_methods, ...
                      'ymin', y_lo,...
                      'ymax', y_hi ...
                     );
    switch geom_type
    case 'geom_bar'
        grammobj.geom_bar('dodge',0.8,'width',0.6);
    case 'geom_point'
        grammobj.geom_point();
        grammobj.set_point_options('base_size',10);
    otherwise
        grammobj.geom_point();
    end
    grammobj.geom_interval('geom','black_errorbar',...
    'dodge',0.8,'width',.5); 
    grammobj.set_color_options('map', mycolormap);
    grammobj.set_names( 'x', 'Method',...
                        'y','Mean Reliability' ...
                        );
    grammobj.set_text_options('base_size',22);
    ylim_upper = round(10*max(y_methods(:)))/10;
    grammobj.axe_property('ylim',[0 max([.6 ylim_upper])]);
    grammobj.no_legend();
    
end