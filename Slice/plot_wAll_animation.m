close all; clear variables;
% load shared scripts
path(path,'../share');

%% flags
l_save_fig = 1;
l_update_data = 1;    % 1 if update data, 0 if read from saved workspace
l_swap = 0; % 1 if saved x-z and y-z slices are at the middle,
            % 0 if those indices are 1
c_gray = [0.5, 0.5, 0.5];

%% Parameters
namelist = {'R8_BF05WD05WV00_ST01_ens03',...
            'R8_BF05WD05WV12_ST01_ens03',...
            'R8_BF05WD05WV00_fL00_ens03',...
            'R8_BF05WD05WV12_fL00_ens03',...
            'R8_BF1hWD00WV00_ST00_ens03'};
tstartlist   = {'028800','028800','021600','021600','021600'};
tendlist = {'032400','032400','025200','025200','028800'};
izilist = [2, 6, 40, 63, 64, 65, 66, 67, 68;...
           2, 6, 40, 63, 64, 65, 66, 67, 68;...
           2, 6, 41, 68, 69, 70, 71, 72, 73;...
           2, 6, 40, 68, 69, 70, 71, 72, 73;...
           2, 6, 44, 77, 78, 79, 80, 81, 82];

ic = 5;
casename = namelist{ic};
tstart = tstartlist{ic};
tend = tendlist{ic};
izi = izilist(ic,:);
nizi = numel(izi);

% set path
get_dataRootDir;
dataDir = [dataRootDir '/viz/'];
outDir = [outRootDir '/slice/' casename '/' tstart '-' tend '/png'];
system(['mkdir -p ' outDir]);

% filename
filename_xy = ['viz.vis.' tstart '.' tend '.xy.nc'];
filename_xz = ['viz.vis.' tstart '.' tend '.xz.nc'];
filename_yz = ['viz.vis.' tstart '.' tend '.yz.nc'];

%% load data
loadSliceData;

Lx = x(1)+x(end);
Ly = y(1)+y(end);
zw = z;

% profile data
fpfl = ncarlesPflData('caseName', casename,...
                      'dataDir',  [dataRootDir '/hist/'],...
                      'fPrefix',  'his.mp.vis',...
                      'tStart',   str2double(tstart)+1,...
                      'tEnd',     str2double(tend)+1);

% boundary layer depth
hb = -mean(fpfl.getBLDepth);
[~,indhb] = min(abs(hb+z));
indentr = (0.9*indhb);

%% plot w
% all figures together
fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 5 7];

ind1 = 2;
ind2 = 3;
[~, ind3] = min(abs(izi-0.9*indhb));
hold on;
[xx,yy,zz] = meshgrid(x./Lx,y./Ly,zw(1:indhb)./hb);
zs = [zw(izi(ind1)), zw(izi(ind2)), zw(izi(ind3))];
vol = zeros(nx,ny,indhb);

% loop over time
for it=1:nt
for i=1:nizi
    tmp = w_xy(:,:,i,it)';
    tmpstd = std(tmp(:));
    vol(:,:,izi(i)) = tmp./tmpstd;
end
% vol_std = std(vol(:));
% vol = vol./vol_std;
xx = double(xx);
yy = double(yy);
zz = double(zz);
vol = double(vol);
zlev = double(zs./hb);
ps = slice(xx,yy,zz,vol,[],[],zlev); shading flat;
colorbar;
vmax = 3;
colormap(db2dr(-vmax,vmax));
xlabel('$x/L_x$', 'Interpreter', 'latex');
ylabel('$y/L_y$', 'Interpreter', 'latex');
zlabel('$z/h_b$',  'Interpreter', 'latex')
view(-25, 35);
xlim([0,1]);
ylim([0,1]);
zlim([-1,0]); 
daspect([1,1,0.6]);
colorbar('off');
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

if it == 1
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end
    
if l_save_fig
    it_str = sprintf('%04d', it);
    figname = [outDir '/w_all_' it_str '.fig'];
%     saveas(gcf, figname, 'fig');
    [dir, name, ~] = fileparts(figname);
%     print('-depsc2',[dir '/' name]);
    print('-dpng', [dir '/' name],'-r200');
end
end


%% functions
function f = plotW(x, y, dat)
    figure;
    f = pcolor(x, y, dat);
    shading flat;
    colorbar;
    vmax = std(dat(:))*5;
    colormap(db2dr(-vmax,vmax));
    daspect([1, 1, 1]);
end

