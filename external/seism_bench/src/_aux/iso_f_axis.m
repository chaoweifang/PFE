function iso_f_axis(figure_name,ftsizet)
% iso_f_axis(figure_name)
% ------------------------------------------------------------------------ 
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
% 
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  March 2011
% ------------------------------------------------------------------------
% This function plots a precision-recall axis with iso-f curves overlaid
% ------------------------------------------------------------------------
    if nargin==0
        figure_name = 'iso_f';
    end

    curr_fig = findobj('type','figure','name',figure_name);
    if isempty(curr_fig)
        figure('name',figure_name)
        width = 300; height = 300;
        ml = 45; mr = 3;
        mb = 40; mt = 20;
        set(gcf,'Units','points','position',[100 100 width+ml+mr height+mt+mb])
        axes('Units','points','position',[ml mb width height])
        ylabel('Precision')
        xlabel('Recall')
        axis equal
        axis([0 1 0 1])
        set(gca,'XTick',0:0.1:1)
        set(gca,'YTick',0:0.1:1)
        set(gca,'FontName','times','FontSize',ftsizet) 
       
        grid on
        box on
        hold on

        for F = 0.1:0.1:0.9
            iso_f_plot(F)
        end

    else
        figure(curr_fig)
    end
end

