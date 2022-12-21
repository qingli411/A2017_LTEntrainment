function h = newFigure(l_save_fig)
% newFigure() creates new figure object.
%   Visible if l_save_fig = 0, invisible otherwise.

    if (l_save_fig)
        figure('visible','off');
    else
        figure;
    end
    h = gcf;
end
