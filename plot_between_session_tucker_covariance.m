function plot_between_session_tucker_covariance(X,varargin)
    % 
    % INPUT
    %   X - is a cell array;
    
    nsubjects = length(X);
    
    
    savefilename = 'tmp/Mar-22-2018/best-eegfmri/cpac/tucker_covariances';
    mkdir('tmp/Mar-22-2018/best-eegfmri/');
    mkdir('tmp/Mar-22-2018/best-eegfmri/cpac');
    
    for subject=1:nsubjects
        Y = cell2mat(X{subject});
        [t p s] = size(Y); 
        Y = reshape(Y,[t p*s]);
        imagesc(corr(standardize.successive_normalize(Y)));
        axis image;
        title(sprintf('Tucker ROI Correlation in Subject %d',subject));
        if(subject==1)
            export_fig(savefilename,'-pdf','-q101');
        else
            export_fig(savefilename,'-pdf','-q101','-append');
        end
        clear s;
    end
    
    
end