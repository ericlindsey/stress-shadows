function stretch_fig_no_whitespace(fig,scale)
    % function by Eric Lindsey, June 2019 derived from this mathworks example:
    % https://www.mathworks.com/help/matlab/creating_plots/save-figure-with-minimal-white-space.html
%%

    % stretch all subplot axes 
    axes=get(fig,'children');
    axestypes=get(get(fig,'children'),'type');
    for i=1:length(axes)
        if strcmp(axestypes{i}, 'axes')
            ax=axes(i);
            outerpos = ax.OuterPosition;
            ti = ax.TightInset; 
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - scale*ti(3);
            ax_height = outerpos(4) - ti(2) - scale*ti(4);
            ax.Position = [left bottom ax_width ax_height];
        end
    end

    % set figure to fill paper
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];

end