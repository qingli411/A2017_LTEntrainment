close all; clear variables;
% load shared script
path(path,'../share');

%% data info
tstart = '032400';
tend = '036000';
leg = 1;    % 0 if first leg; 1 if second leg
casename = 'R8_BF05WD05WV00_ST01_ens02';

h0 = 42;    % initial mixed layer depth
izi = [2, 6, 40, 65, 72];   % indices of saved x-y slices
ns = 1; % time step to start

l_swap = 0; % 1 if saved x-z and y-z slices are at the middle,
            % 0 if those indices are 1
l_update_data = 1;    % 1 if update data, 0 if read from saved workspace
l_savegif = 0;  % 1 if save gif
l_savepng = 1;  % 1 if save png for all frames

% filename
filename_xy = ['viz.vis.' tstart '.' tend '.xy.nc'];
filename_xz = ['viz.vis.' tstart '.' tend '.xz.nc'];
filename_yz = ['viz.vis.' tstart '.' tend '.yz.nc'];

% set up directories
get_dataRootDir;    % get dataRootDir and outRootDir
dataDir = [dataRootDir '/viz/'];
outDir = [outRootDir '/slice/' casename '/' tstart '-' tend];
system(['mkdir -p ' outDir]);

% output figure name
figname = [outDir '/wSliceAll.gif'];
[fdir, fname, ~] = fileparts(figname);
pngdir = [fdir '/png'];
system(['mkdir -p ' pngdir]);

% set up profile data parameters
pflDir = [dataRootDir '/hist'];
fPrefix = 'his.mp.vis';
% creating object ncarlesPflData
f = ncarlesPflData('caseName', casename,...
                   'dataDir',  pflDir,...
                   'fPrefix',  fPrefix,...
                   'tStart',   1,...
                   'tEnd',     36001);

%% load data
loadSliceData;
its = f.getTimeIndex(time(1));
ite = f.getTimeIndex(time(end));

tmp = f.getVar('ups');
ups = tmp(:,its:ite);
tmp = f.getVar('vps');
vps = tmp(:,its:ite);
tmp = f.getVar('wps');
wps = tmp(:,its:ite);
tmp = f.getVar('uvle');
upvp = tmp(:,its:ite);
tmp = f.getVar('uwle');
upwp = tmp(:,its:ite);
tmp = f.getVar('vwle');
vpwp = tmp(:,its:ite);
clear tmp;


hb_ts = f.getBLDepth('max_Nsquared');
hb = mean(hb_ts(its:ite));

[~,ind_hb] = min(abs(z-hb));

%% time loop
cmax = 0.02;    % colorbar axis range
dlt = 0.1;  % delay time
ii = 1; % frame counter
nzp = nz/2; % depth limit
ldepth = z(izi);    % depth of x-y slices
il = 2; % level index

% setup figure
fig=figure;
fig.Units = 'inches';
fig.Position = [1 1 12 10];
fig.Color = 'white';

% locations for subplots
% ll = [0.3, 0.03, 0.03, 0.73, 0.73, 0.35];
% bb = [0.5, 0.58, 0.08, 0.58, 0.08, 0.05];
% ww = [0.4, 0.25, 0.25, 0.25, 0.25, 0.3];
% hh = [0.4, 0.33, 0.33, 0.33, 0.33, 0.4];

ll = [0.35, 0.03, 0.03, 0.03, 0.33, 0.64];
bb = [0.37, 0.66, 0.36, 0.06, 0.06, 0.06];
ww = [0.59, 0.30, 0.30, 0.30, 0.30, 0.30];
hh = [0.60, 0.26, 0.26, 0.26, 0.30, 0.30];

for it = ns:nt
    % figure 1: 3D contour, z = 0.1h0, x, y
    ax1 = subplot('position',[ll(1) bb(1) ww(1) hh(1)]);
    vol = zeros(nx,ny,nzp);
    [xx,yy,zz] = meshgrid(x,y,z(1:nzp));
    vol(:,:,1) = w_xy(:,:,il,it)';
    vol(:,1,:) = w_yz(:,1:nzp,1,it);
    vol(1,:,:) = w_xz(:,1:nzp,1,it);
    ps = slice(xx,yy,zz,vol,x(1),y(1),z(1)); shading flat;
    cb_label = '$w$ (m/s)';
    colormap(ax1, 'parula');
    cb = colorbar;
    cb.FontSize = 12;
    ylabel(cb,cb_label,'FontSize',14,...
        'Interpreter','latex');
    daspect([1,1,1]);
    xlim([0,x(nx)]);
    ylim([0,y(ny)]);
    zlim([z(nzp),0]);
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');
    if leg == 0
        title(sprintf('z = %4.1f m;   t = %5.1f min',...
            ldepth(il),time_mm(it)-time_mm(end)));
    else
        title(sprintf('z = %4.1f m;   t = %5.1f min',...
            ldepth(il),time_mm(it)-time_mm(1)));
    end
    caxis([-cmax,cmax]);

    % figure 2-4
    lv = [1,3,4];
    for i=1:3
        k = lv(i);
        ip1 = i+1;
        axx = subplot('Position',[ll(ip1) bb(ip1) ww(ip1) hh(ip1)]);
        pcolor(x,y,w_xy(:,:,k,it)'); shading flat;
        colormap(axx, 'parula');
        caxis([-cmax,cmax]);
        daspect([1,1,1]);
        xlim([0,x(nx)]);
        ylim([0,y(ny)]);
        if i==3
            xlabel('x (m)');
        end
        ylabel('y (m)');
        title(sprintf('z = %3.1f m',ldepth(k)));
    end

    c = zeros([ind_hb,3]);
    vcmin = zeros([ind_hb,3]);
    vcmax = zeros([ind_hb,3]);
    lambda = zeros([ind_hb,3]);
    for i=1:ind_hb
        a = anisotropyTensor( ups(i,it),  vps(i,it),  wps(i,it),...
                              upvp(i,it), upwp(i,it), vpwp(i,it));
        c(i,:) = barycentricCoord(a);
        [lambda(i,:),vcmax(i,:),vcmin(i,:)] = eigMaxMin3(a);
    end
    
    % figure 5: barycentric map
    ax5 = subplot('position',[ll(5) bb(5) ww(5) hh(5)]);
    cb_label = '$z/h_\mathrm{b}$';
    plotAnisotropicBarycentricMap(c,-z(1:ind_hb)./z(ind_hb),0,cb_label);
    colormap(ax5, 'jet');
    hold off;
        
    % figure 6: directional map
    ax6 = subplot('position',[ll(6) bb(6) ww(6) hh(6)]);
    plotEigenVectorDirectionMaxMin(lambda, vcmax, vcmin,...
        -z(1:ind_hb)./z(ind_hb),1,cb_label);
    colormap(ax6, 'jet');
    hold off;
    
    if l_savepng
        it_str = sprintf('%04d', it);
        print('-dpng', [pngdir '/' fname '_' it_str],'-r300');
    end
    
    if l_savegif
        drawnow;
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if ii == 1
            imwrite(imind,cm,figname,'gif','DelayTime',dlt,...
                'Loopcount',inf);
        else
            imwrite(imind,cm,figname,'gif','DelayTime',dlt,...
                'WriteMode','append');
        end
        ii = ii+1;
    end
end
