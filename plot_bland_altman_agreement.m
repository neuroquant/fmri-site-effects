function [grammobj kendallcoef] = plot_altman_agreement(rater1,rater2)
% PLOT_ALTMAN_AGREEMENT
% Plots differences vs. mean
%     

    p1 = length(rater1);
    p2 = length(rater2);
    
	if(exist('brewermap'))	
		colormapfun = @()(flipud(brewermap(length(colormap),'RdYlBu')));
		close all;
	else
		colormapfun = @winter;
	end
    
    if(ndims(rater1)==2)
        trid_idx = find(triu(ones(p1,p1),1));
        rater1 = rater1(trid_idx);
        rater2 = rater2(trid_idx);
    end

    y_diff = rater1-rater2;
    x_mean = (rater1 + rater2)/2;

    kendallcoef = reliability.kendallsW(cat(2,rater1,rater2));
    disp('Kendalls Agreement (W)')
    kendallcoef

    geom_type = 'geom_point';
    
    grammobj = gramm('x',x_mean,...
                     'y',y_diff);
    switch geom_type
    case 'qq'
        grammobj.stat_qq();
    case 'geom_bar'
        grammobj.geom_bar('dodge',0.8,'width',0.6);
    case 'geom_point'
        grammobj.geom_point();
        grammobj.set_point_options('base_size',5);
    otherwise
        grammobj.geom_point();
    end
    grammobj.geom_hline('yintercept',.1,'style','k--');
    grammobj.geom_hline('yintercept',-.1,'style','k--');
     
    % grammobj.geom_interval('geom','black_errorbar',...
    % 'dodge',0.8,'width',1);
    grammobj.set_color_options('map', colormapfun());
    grammobj.set_names('x', 'Mean',...
                        'y','Difference' ...
                        );
    grammobj.set_title(sprintf('Agreement (Kendalls W): %2.3f',kendallcoef));
    grammobj.set_text_options('base_size',20,'title_scaling',1.1);
    grammobj.axe_property('xlim',[-0.3 1],'ylim',[-.65 .65]);
    grammobj.no_legend();
    
    grammobj.draw()

    
end