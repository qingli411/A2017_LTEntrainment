function ph = plotDataScatter(xdat,ydat,...
                            xlabel_str, ylabel_str,...
                            xlog,ylog,xlims,ylims,...
                            lMColor,lMarker,l_filled,...
                            varargin)
%
    % by default no error bar
    l_errbar = 0;
    % count arguments
    nArgs = length(varargin);
    if nArgs==3
        ydat1 = varargin{1};
        err_n = varargin{2};
        err_p = varargin{3};
        l_errbar = 1;
    elseif nArgs>0
        error(['h = plotScatterAll(xdat0,ydat0,'...
               'xlabel_str, ylabel_str,'...
               'xlog,ylog,xlims,ylims,'...
               'lColor,lMarker,lgname,'...
               '[ydat1],[err0_n],[err0_p])']);
    end

    hold on;
    % constants
    dsize = size(xdat);
    nc = dsize(1);

    % plot
    for ic = 1:nc
        p = plot(xdat(ic),ydat(ic));
        p.LineWidth = 1.0;
        p.Marker = lMarker{ic};
        p.MarkerEdgeColor = lMColor{ic};
        if l_filled{ic}
            p.MarkerFaceColor = lMColor{ic};
        else
            p.MarkerFaceColor = 'none';
        end
        if l_errbar
            ep = errorbar(xdat(ic),ydat1(ic),err_n(ic),err_p(ic));
            ep.LineWidth = 0.5;
            ep.Color = lMColor{ic};
            ep.Marker = '+';
        end
    end

    xlabel(xlabel_str,'Interpreter','Latex');
    ylabel(ylabel_str,'Interpreter','Latex');

    if xlog
        set(gca,'xscale','log');
    end
    if ylog
        set(gca,'yscale','log');
    end
    
    if ~isempty(xlims)
        xlim(xlims);
    end
    if ~isempty(ylims)
        ylim(ylims);
    end
    ph = p;
end