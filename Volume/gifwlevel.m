close all; clear all;
% load shared script
path(path,'../share');
path(path,'../share/cbrewer')

%% data info
% casename = 'R2_lap30_ds4p8_st00_ens02';
casename = 'R2_lainf_dsnan_st00_ens02';
tstart = '028800';
tend = '032400';
ns = 1;

get_dataRootDir;    % get dataRootDir
dataDir = [dataRootDir '/viz/'];
outDir = [outRootDir '/volume'];

filename_xy = ['viz.vis.' tstart '.' tend '.xy.nc'];

h0 = 42;
rl = [0.1, 0.2, 0.5, 0.8, 1.0];
ldepth = -h0.*rl;

l_savegif = 1;

disp('Reading data...');
%% read data
w_xy = ncread([dataDir casename '/' filename_xy],'w');
x = ncread([dataDir casename '/' filename_xy],'x');
y = ncread([dataDir casename '/' filename_xy],'y');
time = ncread([dataDir casename '/' filename_xy],'time');
time_mm = time./60;

nsize = size(w_xy);
nx = nsize(1);
ny = nsize(2);
nl = nsize(3);
nt = nsize(4);

%% loop over levels
dlt = 0.1;
for i=1:nl
    lev_str = sprintf('lev%d',i);
    figname = [casename '_wSlice_' lev_str '.gif'];
    ii = 1;
    fig = figure;
    fig.Color = 'white';
%% time loop
cmax = 0.02;
% tmpdir = [outDir '/wlevelSTLT/' lev_str];
tmpdir = [outDir '/wlevelSTLT'];
system(['mkdir -p ' tmpdir]);
gifname = [tmpdir '/' casename '_wlevel_' lev_str '.gif'];
for it = ns:nt
    w_dat = squeeze(w_xy(:,:,i,it))';
    if it==ns
        cmax = max(abs(w_dat(:)))*0.5;
    end
    pcolor(x,y,w_dat); shading flat;colorbar;
    caxis([-cmax,cmax]);
    daspect([1,1,1]);
    xlim([0,x(nx)]);
    ylim([0,y(ny)]);
    xlabel('x (m)');
    ylabel('y (m)');
    title(sprintf('z = %3.1f m',ldepth(i)));

%     ttag = sprintf('%04d',it);
%     figname = [tmpdir '/tmp_' ttag '.eps'];
%     print('-depsc2',figname);
if l_savegif
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if ii == 1
        imwrite(imind,cm,gifname,'gif','DelayTime',dlt,'Loopcount',inf);
    else
        imwrite(imind,cm,gifname,'gif','DelayTime',dlt,'WriteMode','append');
    end
    ii = ii+1;
end
end
% system(['/opt/local/bin/convert -density 150 -delay 10 '...
%     tmpdir '/tmp_*.eps ' gifname]);
end
