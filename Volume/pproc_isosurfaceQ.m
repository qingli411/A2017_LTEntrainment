close all; clear variables;
% case name
casenames = {'R8_BF05WD05WV00_ST01_ens03',...
             'R8_BF05WD05WV12_ST01_ens03'};
figs = {'Q_dom1',...
        'Q_dom2',...
        'Q_dom3',...
        'VM_dom1',...
        'VM_dom2',...
        'VM_dom3',...
        };
wfig = 'w_level';
nf = numel(figs);
nc = numel(casenames);

wsize = 6;
hsize = [2.4, 3.2, 3];

get_dataRootDir;    % get dataRootDir
figdir = [outRootDir '/volume/'];
tttttt = '032401';

angle_view1 = [20,20,20];
angle_view2 = [-10,20,-10];
xdom = 80;
ydom = 80;
x0 = [65,4];
y0 = [90,190];

for j=1:nc
    for i=1:nf
        % process isosurface plots
        filename = [figdir casenames{j} '/' tttttt '/' figs{i} '.fig'];
        % load figure
        [dir, name, ~] = fileparts(filename);
        open(filename);

        % change xlim, ylim
        xmin = x0(j);
        xmax = x0(j)+xdom;
        xlim([xmin, xmax]);
        ymin = y0(j);
        ymax = y0(j)+ydom;
        ylim([ymin, ymax]);
        zlims = zlim;
        hold on;
        p = addArrows(x0(j)+5,y0(j)+5,zlims(1));
        % change view angle
        ni = mod(i-1,3)+1;
        view(angle_view1(ni),angle_view2(ni))

        % Change size
        mleft = 0.15;
        mdown = 0.15;
        mright = wsize-0.15;
        mup = hsize(ni)-0.15;
        set(gcf,'Units','inches');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperPositionMode', 'manual');
        paperposition = [mleft mdown mright mup];
        set(gcf,'PaperPosition',paperposition);
        size = [wsize hsize(ni)];
        set(gcf,'PaperSize',size);
        position = [mleft mdown mright mup];
        set(gcf,'Position',position);
        a = findall(gcf,'-property','FontSize');
        for k=1:numel(a)
            if strcmp(a(k).Type,'axes')
                a(k).FontSize = 12;
            elseif strcmp(a(k).Type,'text')
                a(k).FontSize = 14;
            else
                a(k).FontSize = 12;
            end
        end
        print('-depsc2',[dir '/' name]);
    end
    % process w plot
%     filename = [figdir casenames{j} '/' tttttt '/' wfig '.fig'];
%     % load figure
%     [dir, name, ~] = fileparts(filename);
%     open(filename);
%     h = gcf;
%     for k=1:3
%         hc = h.Children(4-k);
%         hold on;
%         addBox(hc,x0(j),y0(j),xdom,ydom);
%     end
end

function addBox(h, x0, y0, xl, yl)
    bcolor = 'w';
    rectangle(h,'Position',[x0,y0,xl,yl],'EdgeColor',bcolor);
end

function [p1, p2, p3] = addArrows(x0,y0,z0)
    u = 10;
    v = 10;
    w = 10;
    p1 = addArrow(x0, y0, z0, u, 0, 0, 'x');
    p2 = addArrow(x0, y0, z0, 0, v, 0, 'y');
    p3 = addArrow(x0, y0, z0, 0, 0, w, 'z');

end

function p = addArrow(x0, y0, z0, u, v, w, ax)
    p = quiver3(x0, y0, z0, u, v, w);
    p.LineWidth = 1.5;
    p.Color = 'k';
    p.AutoScaleFactor = 0.8;
    p.MaxHeadSize = 0.5;
    if strcmp(ax,'x')
        xx = x0+u;
        yy = y0;
        zz = z0;
    elseif strcmp(ax, 'y')
        xx = x0;
        yy = y0+v;
        zz = z0;
    elseif strcmp(ax, 'z')
        xx = x0;
        yy = y0;
        zz = z0+w;
    else
        error('Unknown ax label');
    end
    tx = text(xx,yy,zz,['$' ax '$']);
    tx.FontSize = 9;
    tx.HorizontalAlignment = 'center';
    tx.Interpreter = 'latex';
end
