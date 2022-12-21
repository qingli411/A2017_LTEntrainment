close all; clear variables;
% add shared code
path(path,'../share');
path(path,'../share/cbrewer')

casename = 'R9_BF05WD05WV12_ST00_ens01';
tttttt = '014401';

get_dataRootDir;    % get dataRootDir
dataDir = [dataRootDir '/rest/'];
inDir = [outRootDir '/volume/' casename '/' tttttt];

% get grid
srcname = [dataDir casename '/les.F' ];
nnx = getValue(srcname, 'nnx');
nny = getValue(srcname, 'nny');
Lx  = getValue(srcname, 'xl');   % x-domain size in meter
Ly  = getValue(srcname, 'yl');   % y-domain size in meter

dx = Lx/nnx;
x = dx/2:dx:Lx;     % center of the grid, x
dy = Ly/nny;
y = dy/2:dy:Ly;     % center of the grid, y


load([inDir '/w_maxwps.mat']);

[xx2d,yy2d] = meshgrid(x, y);
figure;
wdat = permute(w_maxwps, [2, 1]);
% cmax = max(abs(wdat(:)))*0.3;
% w at the level of maximum wps
ps1 = pcolor(xx2d,yy2d,wdat);
shading flat;
% caxis([-cmax, cmax]);
colorbar;
vmax = std(wdat(:))*3;
colormap(db2dr(-vmax,vmax));
colorbar('off');
daspect([1 1 1]);
xlabel('x (m)');
ylabel('y (m)');
title(['z = ' sprintf('%4.2f', z_maxwps) ' m']);
axis tight

figname = 'w_maxwps_rb';
saveas(gcf, [figname '.fig'], 'fig');
print('-depsc2',[figname '.eps'],'-r400');
