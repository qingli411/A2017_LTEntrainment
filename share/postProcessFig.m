function postProcessFig(filename,varargin)
% postProcessFig resizes a figure, convert it to eps format
%   postProcessFig(figname,[width, hight])

% Qing Li, 150225
%          160505, support width and hight input
%          170208, support fontsize for inlets

% counts arguments
nArgs = length(varargin);
if nArgs==0
    wsize = 6;
    hsize = 4;
elseif nArgs==2
    if isnumeric(varargin{1}) && isnumeric(varargin{2})
        wsize = varargin{1};
        hsize = varargin{2};
    else
        error('Numeric arguments required...');
    end
else
    error('No arguments or pass in width and hight');
end
mleft = 0.25;
mdown = 0.25;
mright = wsize-0.25;
mup = hsize-0.25;

% load figure
[dir, name, ~] = fileparts(filename);
if isempty(dir)
    dir = './';
end
open(filename);
% Get x, y limits
xx = xlim;
yy = ylim;
% Change size
set(gcf,'Units','inches');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPositionMode', 'manual');
%set(gcf,'PaperOrientation','landscape');
paperposition = [mleft mdown mright mup];
set(gcf,'PaperPosition',paperposition);
size = [wsize hsize];
set(gcf,'PaperSize',size);
position = [mleft mdown mright mup];
set(gcf,'Position',position);
% Change renderer, 'painters' works well for vector formats while
% 'opengl' works well for bitmap format
set(gcf, 'renderer', 'painters');
% Change font
a = findall(gcf,'-property','FontSize');
for i=1:numel(a)
    if strcmp(a(i).Type,'axes')
        a(i).FontSize = 12;
    elseif strcmp(a(i).Type,'text')
        a(i).FontSize = 14;
    else
        a(i).FontSize = 12;
    end
end
% set(findall(gcf,'-property','FontSize'),'FontSize',12)
% set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')
% Export file
% check if contains inlet
h = findall(gcf,'Type','Axes');
if length(h)==2
    h(1).FontSize = 9;
end
xlim(xx);
ylim(yy);
print('-depsc2',[dir '/' name]);
end
