close all; clear variables;
% load shared scripts
path(path,'../share');

%% flags
l_save_fig = 1;

%% Parameters
get_dataRootDir;
dataDir = [dataRootDir '/rest/'];
casename = 'R8_BF05WD05WV00_ST01_ens02';
[s0, r0] = system(['ls ' dataDir casename ]);
if ~s0
    tttttt = r0(1:6);
end
% u file name
pattern = 'u\.[a-z]{2}\.[a-z]{3}[0-9]{3}';
dir = [dataDir casename '/' tttttt '/'];
[s0, r0] = system(['ls ' dir]);
if ~s0
    r1 = regexp(r0, pattern, 'match');
    fileIn = r1{1};
end
% output directory
outDir = [outRootDir '/volume/' casename '/' tttttt];

% get grid
srcname = [dataDir casename '/les.F' ];
nnz = getValue(srcname, 'nnz');
nnx = getValue(srcname, 'nnx');
nny = getValue(srcname, 'nny');
Lx  = getValue(srcname, 'xl');   % x-domain size in meter
Ly  = getValue(srcname, 'yl');   % y-domain size in meter
Lz  = getValue(srcname, 'zl');   % z-domain size in meter
nscl = getValue(srcname, 'nscl'); % number of scalers
nvar = 4+nscl; % number of variables
dx = Lx/nnx;
x = dx/2:dx:Lx;     % center of the grid, x
dy = Ly/nny;
y = dy/2:dy:Ly;     % center of the grid, y
dz = Lz/nnz;
zw  = dz:dz:Lz;     % z for w
zu  = dz/2:dz:Lz;   % z for u, v, t
batag = 9.81/5000;

% read file
fid = fopen([dir fileIn]);
tmp = fread(fid, nvar*nnx*nny*nnz, 'double', 'l');
fclose(fid);
clear var;
var = reshape(tmp, nvar, nnx, nny, nnz);
clear tmp;

%% u, v, w at maximum wps depth
% get maximum wps depth

w = squeeze(var(3,:,:,:));
wps = zeros(1,nnz);
for i=1:nnz
    tmp = w(:,:,i);
    tmp2 = (tmp-mean(tmp(:))).^2;
    wps(i) = mean(tmp2(:));
end
[~,mi] = max(wps(:));

pre_str = 'maxwps';
depth1 = zw(mi);
ind1 = mi;
fprintf('Depth: %g\n', depth1);

h1 = uvwSpectrum(var, mi, x, y, Lx, pre_str, l_save_fig);


%% u, v, w at hb/2
% get boundary layer depth

N2 = zeros(1,nnz);
txym = zeros(1,nnz);
for i=1:nnz
    tmp = squeeze(var(4,:,:,i));
    txym(i) = mean(tmp(:));
end
N2(1:end-1) = (txym(2:end)-txym(1:end-1))./dz;
N2(end) = N2(end-1);
[~,mi0] = max(N2(:));

mi = round(mi0/2);
pre_str = '0p5h';
depth2 = zw(mi);
ind2 = mi;
fprintf('Depth: %g\n', depth2);

h2 = uvwSpectrum(var, mi, x, y, Lx, pre_str, l_save_fig);


%% u, v, w near hb: z = 0.9*hb

mi = round(mi0*0.9);
pre_str = '0p9h';
depth3 = zw(mi);
ind3 = mi;
fprintf('Depth: %g\n', depth3);

h3 = uvwSpectrum(var, mi, x, y, Lx, pre_str, l_save_fig);

%% u, v, w near surface: z = 1.28 m

mi = 2;
pre_str = 'lev2';
depth4 = zw(mi);
ind4 = mi;
fprintf('Depth: %g\n', depth4);

h4 = uvwSpectrum(var, mi, x, y, Lx, pre_str, l_save_fig);

%% wps profile
newFigure(l_save_fig);
plot(wps./wps(ind1), -zw./zw(mi0), '-k', 'LineWidth', 1.5);
hold on;
plot(wps([ind1, ind2, ind3, ind4])./wps(ind1),...
    -[depth1, depth2, depth3, depth4]./zw(mi0), '+k', 'LineWidth', 1.5);
xlabel('$\overline{w''^2}/{\overline{w''^2}}_\mathrm{max}$',...
    'Interpreter', 'latex');
ylabel('$z/h_\mathrm{b}$', 'Interpreter', 'latex');
xlim([0, 1.2]);
ylim([-1.2, 0]);
if l_save_fig
    figname = './profile_wps.fig';
    saveas(gcf, figname, 'fig');
    postProcessFig(figname, 5, 5);
end

%% save data
if l_save_fig

    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    system(['mv snapshot* ' outDir]);
    system(['mv spectrum* ' outDir]);
    system(['mv profile* ' outDir]);
end


%% functions

function h = uvwSpectrum(var, mi, x, y, Lx, pre_str, l_save_fig)
    % mi: level of zw

    c_gray = [0.5, 0.5, 0.5];
    uli = squeeze(var(1, :, :, mi));
    ulp = squeeze(var(1, :, :, mi+1));
    ul  = 0.5.*(uli+ulp);
    vli = squeeze(var(2, :, :, mi));
    vlp = squeeze(var(2, :, :, mi+1));
    vl  = 0.5.*(vli+vlp);
    wl  = squeeze(var(3, :, :, mi));
    up  = ul - mean(ul(:));
    vp  = vl - mean(vl(:));
    wp  = wl - mean(wl(:));
    wps = wp.^2;
    wgt = sqrt(mean(wps(:)));
    % plot figure w
    newFigure(l_save_fig);
    pcolor(x,y,wp'./wgt); shading flat; colorbar;
    caxis([-3, 3]);
    xlabel('$x$ (m)', 'Interpreter', 'latex');
    ylabel('$y$ (m)', 'Interpreter', 'latex');
    daspect([1 1 1]);

    if l_save_fig
        figname = ['./snapshot_wp_' pre_str '.fig'];
        saveas(gcf, figname, 'fig');
        postProcessFig(figname, 6, 4);
    end

    % Calculate spectra
    [ew_x, ew_y, ~] = spectra2D(wp, Lx);
    [eu_x, eu_y, ~] = spectra2D(up, Lx);
    [ev_x, ev_y, ~] = spectra2D(vp, Lx);

    newFigure(l_save_fig);
    hold on;
    wgt = eu_y(2)+ev_y(2)+ew_y(2);
    p1 = plot(eu_x,eu_y./wgt, '-k',  'LineWidth', 1.5);
    p2 = plot(ev_x,ev_y./wgt, '--k', 'LineWidth', 1.5);
    p3 = plot(ew_x,ew_y./wgt, '-.k', 'LineWidth', 1.5);
    ff = 10./(eu_x(2)).^(-5/3);
    plot(eu_x, ff*eu_x.^(-5/3),'-',...
        'Color', c_gray, 'LineWidth', 0.75);
    nnx = numel(x);
    dx = Lx/nnx;
    xlim([2*pi/Lx, 2*pi/dx/3]);
    ylim([1e-3, 1e1]);
    xlabel('$k$ (m$^{-1}$)', 'Interpreter', 'latex');
    ylabel('$E_i(k)/E(k_0)$', 'Interpreter', 'latex');
    set(gca,'xscale', 'log');
    set(gca,'yscale', 'log');
    lg = legend([p1, p2, p3], '$E_u$','$E_v$','$E_w$');
    lg.Location = 'NorthEast';
    lg.Interpreter = 'latex';

    if l_save_fig
        figname = ['./spectrum_uvw_' pre_str '.fig'];
        saveas(gcf, figname, 'fig');
        postProcessFig(figname, 6, 4);
    end

    h = gcf;
end

function [e_x, e_y, km] = spectra2D(dat, Lx)
% spectra2D calculates the 2D spectrum

    nsize = size(dat);
    nnx = nsize(1);
    nny = nsize(2);
    fw = fftshift(fft2(dat));
    pw = abs(fw).^2;
    ew = zeros(nnx,1);
    km = 0;
    kxc=floor(nnx/2)+1;
    kyc=floor(nny/2)+1;
    for ky=1:nny
        for kx=1:nnx
            ki = sqrt((ky-kyc).^2+(kx-kxc).^2);
            k = round(ki)+1;
            ew(k) = ew(k) + pw(kx,ky);
            km = km + pw(kx,ky).*ki;
        end
    end
    kx = 1:1:nnx/2;
    sca = 2*pi/Lx;
    e_x = sca.*(kx-1);
    e_y = ew(1:nnx/2);
    km = km./sum(pw(:)).*sca;
end
