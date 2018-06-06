function [grammobj prop_wLOA y_diff x_mean] = plot_altman_agreement(rater1,rater2,varargin)
% PLOT_ALTMAN_AGREEMENT
% Plots differences vs. mean
% USAGE: 
%   plot_altman_agreement(rater1,rater2)
%   plot_altman_agreement(rater1,rater2,loa);
% 
% INPUTS
%   - rater1 is n_samples x 1 or n_features x 1
%   - rater2 is n_samples x 1 or n_features x 1
%   (optional)
%   - LOA is a 3 x 2 where each column (1) provides mean/ci for lower agreement interval and column (2) provides mean/ci for upper agreement interval


    p1 = length(rater1);
    p2 = length(rater2);
    
    diff_method = 'MD';
    
	if(exist('brewermap'))	
		colormapfun = @()(flipud(brewermap(length(colormap),'RdYlBu')));
		close all;
	else
		colormapfun = @winter;
	end
    
    if(ndims(rater1)==2 && size(rater1,1)==size(rater1,2))
        trid_idx = find(triu(ones(p1,p1),1));
        rater1 = rater1(trid_idx);
        rater2 = rater2(trid_idx);
    end


    switch diff_method
    case 'raw'
        x_mean = (rater1 + rater2)/2;
        y_diff = rater1-rater2;
    case 'proportion'
        x_mean = (rater1 + rater2)/2;
        y_diff = rater1-rater2;
        y_diff = y_diff./x_mean;
    case 'MD'
        disp('MD')
        y_diff = rater2;
        x_mean = rater1;
        size(y_diff)
        
    end
    
    y_diff_std = std(y_diff); 
    
    switch nargin
        
    case 3
        loa_int = varargin{1};
        
        disp('Proportion of Features within Limits');
        x_strong = 0;
        idx_large = find(abs(x_mean)>x_strong);
        prop_wLOA = 1 - sum((y_diff(idx_large) >  loa_int(1,2)) | ...
                        (y_diff(idx_large) < loa_int(1,1)))/length(idx_large);
        prop_wLOA
        
    otherwise
        disp('Assumes two arguments by default.')
        loa_int = y_diff_std;
        
        disp('Proportion of Strong Features within Limits');
        x_strong = .3;
        idx_large = find(abs(x_mean)>x_strong);
        prop_wLOA = 1 - sum((y_diff(idx_large) >  y_diff_std) | ...
                        (y_diff(idx_large) < - y_diff_std))/length(idx_large);
        prop_wLOA
    end
    

    % kendallcoef = reliability.kendallsW(cat(2,rater1,rater2));
    % disp('Kendalls Agreement (W)')
    % kendallcoef
    
    
    geom_type = 'geom_point';
    
    switch geom_type
    case 'qq'
        
        grammobj = gramm('x',y_diff,...
                         'y',x_mean,'group',abs(x_mean)>x_strong);
        
        grammobj.stat_qq();
        
        grammobj.set_names('y', 'Mean',...
                           'x','Difference' ...
                            );
        
    case 'geom_bar'
        
        grammobj = gramm('y',y_diff,...
                         'x',x_mean);
        
        grammobj.geom_bar('dodge',0.8,'width',0.6);
        
        grammobj.set_names('x', 'Mean',...
                           'y', 'Difference' ...
                            );
        
    case 'geom_point'
        
        grammobj = gramm('y',y_diff,...
                         'x',x_mean);
        
        grammobj.geom_point();
        grammobj.geom_jitter('width',.05);
        grammobj.set_point_options('base_size',5);
        
        %
        grammobj.geom_hline('yintercept',0,'style','k-');
        grammobj.geom_hline('yintercept',mean(y_diff),'style','k-')
        if(isscalar(loa_int))
            grammobj.geom_hline('yintercept',y_diff_std,'style','k--');
            grammobj.geom_hline('yintercept',-y_diff_std,'style','k--');
        else
            loa_int
            grammobj.geom_hline('yintercept',loa_int(1,1),'style','r--');
            grammobj.geom_hline('yintercept',loa_int(1,2),'style','r--');
            % grammobj.geom_hline('yintercept',loa_int(2,1),'style','k-.');
            % grammobj.geom_hline('yintercept',loa_int(3,1),'style','k-.');
            % grammobj.geom_hline('yintercept',loa_int(2,2),'style','k-.');
            % grammobj.geom_hline('yintercept',loa_int(3,2),'style','k-.');
        end
        grammobj.set_names('x', 'Mean',...
                           'y', 'Difference' ...
                            );
        
    otherwise
        grammobj.geom_point();
    end
    
    
    % grammobj.geom_interval('geom','black_errorbar',...
    % 'dodge',0.8,'width',1);
    grammobj.set_color_options('map', colormapfun());
    grammobj.set_title(sprintf('Bland Altman MD (LOA)'));
    grammobj.set_text_options('base_size',20,'title_scaling',1.1);
    switch geom_type
    case 'geom_point'
        grammobj.axe_property('xlim',[-0.3 1],'ylim',[-.4 .4]);
    end
    grammobj.no_legend();
    
    grammobj.draw()

    
end