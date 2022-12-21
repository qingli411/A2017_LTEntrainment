close all; clear variables;
% load shared script
path(path,'../share');

%% data info
tstart = '028800';
tend = '032400';
casename = 'R8_BF05WD05WV12_ST01_ens03';

h0 = 42;    % initial mixed layer depth
izi = [2, 6, 40, 65, 72];   % indices of saved x-y slices
ns = 0; % time step to start

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
figname = [outDir '/wSlice.gif'];
[fdir, fname, ~] = fileparts(figname);
pngdir = [fdir '/png'];
system(['mkdir -p ' pngdir]);

%% load data
loadSliceData;

%% time loop
cmax = 0.02;    % colorbar axis range
dlt = 0.1;  % delay time
ii = 1; % frame counter
nzp = nz/2; % depth limit
ldepth = z(izi);    % depth of x-y slices
il = 2; % level index

% setup figure
f=figure;
f.Units = 'inches';
f.Position = [1 1 6 4];
f.Color = 'white';

for it = ns:nt
    % figure 1: 3D contour, z = 0.1h0, x, y
    vol = zeros(nx,ny,nzp);
    [xx,yy,zz] = meshgrid(x,y,z(1:nzp));
    vol(:,:,1) = w_xy(:,:,il,it)';
    vol(:,1,:) = w_yz(:,1:nzp,1,it);
    vol(1,:,:) = w_xz(:,1:nzp,1,it);
    ps = slice(xx,yy,zz,vol,x(1),y(1),z(1)); shading flat;
    daspect([1,1,1]);
    xlim([0,x(nx)]);
    ylim([0,y(ny)]);
    zlim([z(nz/2),0]);
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');
    title(sprintf('z = %3.1f m; t = %4.1f min',ldepth(il),time_mm(it)));
    cb = colorbar;
    caxis([-cmax,cmax]);

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
