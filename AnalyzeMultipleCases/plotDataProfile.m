function ph = plotDataProfile(dat,wgt,xlab,ylab,xlims,ylims,xlog,...
                            zz,h,lColor,lgname,varargin)

nArg = length(varargin);
if nArg == 2
    l_aux = 1;
    xaux = varargin{1};
    yaux = varargin{2};
elseif nArg == 4
    l_aux = 2;
    xaux = varargin{1};
    yaux = varargin{2};
    xaux2 = varargin{3};
    yaux2 = varargin{4};
elseif nArg == 0
    l_aux =0;
else
    error('12 or 14 arguments required, got %i',12+nArg);
end
dsize = size(dat);
nCase = dsize(1);

% color gray
c_gray = [0.5 0.5 0.5];
% plot object
p = gobjects(nCase,1);
hold on;
% plot data
for ic = 1:nCase
    p(ic) = plot(dat(ic,:)./wgt(ic),zz(ic,:)./h(ic),lColor{ic});
    p(ic).LineWidth = 1.5;
end
xlabel(xlab,'Interpreter','latex');
ylabel(ylab,'Interpreter','latex');
if ~isempty(xlims)
    xlim(xlims);
end
ylim(ylims);
if xlog
    set(gca,'xscale','log');
end
if ~isempty(lgname)
    lg = legend(p,char(lgname));
    lg.Interpreter = 'latex';
    set(lg,'Location','SouthEast');
end
% reference line
xx = xlim;
yy = ylim;
rl1 = line('XData', [xx(1) xx(2)], 'YData', [-1 -1],'LineStyle','-',...
    'Color', c_gray);
uistack(rl1,'bottom');
if (xx(1)<0 && xx(2)>0)
    rl2 = line('XData', [0 0], 'YData', [yy(1) yy(2)], 'LineStyle', '-', ...
    'Color', c_gray);
    uistack(rl2,'bottom');
end
% plot auxiliary lines if requested
if l_aux >= 1
    for i = 1:nCase
        plot(xaux(i,:),yaux(i,:),lColor{i},'LineWidth', 1);
    end
    if l_aux >= 2
        for i = 1:nCase
            plot(xaux2(i,:),yaux2(i,:),lColor{i},'LineWidth', 1);
        end
    end 
end
ph = gcf;

end