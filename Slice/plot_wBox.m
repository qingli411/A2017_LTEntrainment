close all; clear variables;
% load shared script
path(path,'../share');

%% data info
ns = 1; % time step to start

l_swap = 0; % 1 if saved x-z and y-z slices are at the middle,
            % 0 if those indices are 1
l_update_data = 1;    % 1 if update data, 0 if read from saved workspace
l_savepng = 1;  % 1 if save png for all frames
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

ic = 4;
casename = namelist{ic};
tstart = tstartlist{ic};
tend = tendlist{ic};
izi = izilist(ic,:);
nizi = numel(izi);

% set path
get_dataRootDir;
dataDir = [dataRootDir '/viz/'];
if l_savepng
    outDir = [outRootDir '/slice/' casename '/' tstart '-' tend '/png'];
    system(['mkdir -p ' outDir]);
    fname = 'wBox';
end

% filename
filename_xy = ['viz.vis.' tstart '.' tend '.xy.nc'];
filename_xz = ['viz.vis.' tstart '.' tend '.xz.nc'];
filename_yz = ['viz.vis.' tstart '.' tend '.yz.nc'];

%% load data
loadSliceData;
Lx = x(1)+x(end);
Ly = y(1)+y(end);

%% time loop
cmax = 0.02;    % colorbar axis range
dlt = 0.1;  % delay time
ii = 1; % frame counter
nzp = nz/2; % depth limit
ldepth = z(izi(2));    % depth of x-y slices
il = 2; % level index

% setup figure
f=figure;
f.Units = 'inches';
f.Position = [1 1 6 4];
f.Color = 'white';

vol = zeros(nx,ny,nzp);
[xx,yy,zz] = meshgrid(x,y,z(1:nzp));
dz = z(1)-z(2);
for it = ns:nt
    % figure 1: 3D contour, z = 0.1h0, x, y
    vol(:,:,1) = w_xy(:,:,il,it)';
    vol(:,end,:) = w_yz(:,1:nzp,1,it);
    vol(1,:,:) = w_xz(:,1:nzp,1,it);
    ps = slice(xx,yy,zz,vol,x(end),y(1),z(1)); shading flat;
    daspect([1,1,1]);
    hold on;
    xlim([0,x(nx)]);
    ylim([0,y(ny)]);
    zlim([z(nz/2),0]);
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');
    title(sprintf('t = %4.1f min',time_mm(it)));
    cb = colorbar;
    colormap(db2dr(-cmax,cmax));
%     colorbar('off');
    view(30, 30)
    
    % draw boundary layer base
    temp_yz = t_yz(:,1:nzp,1,it);
    N2_yz = (temp_yz(:,1:end-1)-temp_yz(:,2:end))/dz;
    [~,ind_yz] = max(N2_yz,[],2);
    temp_xz = t_xz(:,1:nzp,1,it);
    N2_xz = (temp_xz(:,1:end-1)-temp_xz(:,2:end))/dz;
    [~,ind_xz] = max(N2_xz,[],2);
    ind_xz = idx_smooth(ind_xz);
    ind_yz = idx_smooth(ind_yz);
    x3 = x(end).*ones(size(x));
    y3 = y(1).*ones(size(y));
    
    plot3(x,y3,z(ind_xz),'-','Color',c_gray,'LineWidth',1);
    plot3(x3,y,z(ind_yz),'-','Color',c_gray,'LineWidth',1);
    hold off;
    
    % save png
    if l_savepng
        it_str = sprintf('%04d', it);
        print('-dpng', [outDir '/' fname '_' it_str],'-r200');
    end

end

function idx_new = idx_smooth(idx)
    ni = numel(idx);
    idx_new = idx;
    idx_ref = floor(median(idx(:)));
    for i = 1:ni
        didx = idx(i)-idx_ref;
        if didx > 20
            idx_new(i) = idx_ref;
        end
    end
end