function h = plotDataProfile(dat,wgt,label,zz,h,lColor,...
                            lgname,l_save_fig,figname)

dsize = size(dat);
nCase = dsize(1);
newFigure(l_save_fig);
% plot object
p = gobjects(nCase,1);
hold on;
for ic = 1:nCase
    p(ic) = plot(dat(ic,:)./wgt(ic),zz(ic,:)./h(ic),lColor{ic});
    p(ic).LineWidth = 1.5;
end
xlabel(label,'Interpreter','latex');
ylabel('$z/h_\mathrm{b}$','Interpreter','latex');
ylim([-1.2,0]);
if ~isempty(lgname)
    lg = legend(p,char(lgname));
    lg.Interpreter = 'latex';
    set(lg,'Location','SouthEast');
end
% reference line
xx = xlim;
yy = ylim;
line('XData', [xx(1) xx(2)], 'YData', [-1 -1],'LineStyle','--',...
    'LineWidth', 1, 'Color', 'k');
if (xx(1)<0 && xx(2)>0)
    line('XData', [0 0], 'YData', [yy(1) yy(2)], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color', 'k');
end
hold off;
if (l_save_fig)
    figName = [figname '.fig'];
    saveas(gcf,figName,'fig');
    postProcessFig(figName,5,5);
end