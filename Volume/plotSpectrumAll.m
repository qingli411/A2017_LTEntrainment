close all; clear variables;
% load shared scripts
path(path,'../share');

%% flags
l_save_fig = 0;
global cgray;
cgray = [0.5, 0.5, 0.5];

%% Parameters
get_dataRootDir;
dataDir = [dataRootDir '/rest/'];

namelist = {'R8_BF05WD05WV00_ST01_ens03',...
            'R8_BF05WD05WV12_ST01_ens03',...
            'R8_BF05WD05WV00_fL00_ens03',...
            'R8_BF05WD05WV12_fL00_ens03'};
restlist   = {'032401','032401','025201','025201'};
izilist = [64, 66, 69, 71];
% casename = 'R8_BF05WD05WV00_fL00_ens03';
% [s0, r0] = system(['ls ' dataDir casename ]);
% if ~s0
%     tttttt = r0(1:6);
% end
% tttttt = '025201';
ic = 1;
casename = namelist{ic};
tttttt = restlist{ic};
izi = izilist(ic);
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
system(['mkdir -p ' outDir]);

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

%% spectrum
% u, v, w at maximum wps depth
% get maximum wps depth
w = squeeze(var(3,:,:,:));
wps = zeros(1,nnz);
for i=1:nnz
    tmp = w(:,:,i);
    tmp2 = (tmp-mean(tmp(:))).^2;
    wps(i) = mean(tmp2(:));
end
[~,mi] = max(wps(:));

depth1 = zw(mi);
ind1 = mi;
fprintf('Depth: %g\n', depth1);
[xx1, yy1] = calSpectrum(var, mi, Lx);

% u, v, w at hb/2
% get boundary layer depth
N2 = zeros(1,nnz);
txym = zeros(1,nnz);
for i=1:nnz
    tmp = squeeze(var(4,:,:,i));
    txym(i) = mean(tmp(:));
end
N2(1:end-1) = (txym(2:end)-txym(1:end-1))./dz;
N2(end) = N2(end-1);
[~,indhb] = max(N2(:));
hb = zw(indhb);

mi = round(indhb/2);
depth2 = zw(mi);
ind2 = mi;
fprintf('Depth: %g\n', depth2);
[xx2, yy2] = calSpectrum(var, mi, Lx);

% u, v, w near hb: z = 0.9*hb
% mi = round(indhb*0.9);
mi = izi;
depth3 = zw(mi);
ind3 = mi;
fprintf('Depth: %g\n', depth3);
[xx3, yy3] = calSpectrum(var, mi, Lx);

% u, v, w near surface: z = 1.28 m
mi = 2;
depth4 = zw(mi);
ind4 = mi;
fprintf('Depth: %g\n', depth4);
[xx4, yy4] = calSpectrum(var, mi, Lx);

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 4 7];
hold on;
xx = [xx3, xx2, xx1, xx4];
yy = [yy3, yy2, yy1, yy4];

ylims = [1e-1, 2e11];
refk(-2.*pi./hb, ylims);
refk(-2.*pi./hb*4, ylims);
% rl2=refk(-pi./sqrt(2)/depth2, ylims);
% rl3=refk(-pi./sqrt(2)/depth3, ylims);
% rl4=refk(-pi./sqrt(2)/depth4, ylims);

for i = 1:4
    eu_y = yy(i).eu;
    ev_y = yy(i).ev;
    ew_y = yy(i).ew;
    ek_y = yy(i).ek;
    eu_x = xx(i).eu;
    ev_x = xx(i).ev;
    ew_x = xx(i).ew;
    ek_x = xx(i).ek;
    factor = 10^((i-1)*3);
    nx = numel(ew_x);
    wgt = ek_y(nx/4)/factor;

    ff = 1.*factor./(ek_x(nx/4)).^(-5/3);
    plot(eu_x, ff*eu_x.^(-5/3),'-',...
        'Color', cgray, 'LineWidth', 0.75);
%     p0 = plot(ek_x,ek_y./wgt, '-k', 'LineWidth', 1.5);
    p1 = plot(eu_x,eu_y./wgt, '-k', 'LineWidth', 1.5);
    p2 = plot(ev_x,ev_y./wgt, '-b', 'LineWidth', 1.5);
    p3 = plot(ew_x,ew_y./wgt, '-r', 'LineWidth', 1.5);
end
nnx = numel(x);
dx = Lx/nnx;
% xlim([2*pi/Lx, 2*pi/dx/3]);
xlim([2*pi/Lx, 2*pi/dx/2/sqrt(2)]);
ylim(ylims);
xticks([1e-1, 1]);
xlabel('$k$ (m$^{-1}$)', 'Interpreter', 'latex');
ylabel('$E_i(k)$', 'Interpreter', 'latex');
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
% lg = legend([p1, p2, p3],'$E_u$','$E_v$','$E_w$');
% lg.Location = 'NorthEast';
% lg.Interpreter = 'latex';
pbaspect([1, 2.5, 1]);


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

if l_save_fig
    figname = [outDir '/spectrum_uvw_all.fig'];
    saveas(gcf, figname, 'fig');
    [dir, name, ~] = fileparts(figname);
    print('-depsc2',[dir '/' name]);
end

%% ups, vps, wps profile

% get ups and vps
u = squeeze(var(1,:,:,:));
ups = zeros(1,nnz);
v = squeeze(var(2,:,:,:));
vps = zeros(1,nnz);
for i=1:nnz-1
    tmp = 0.5.*(u(:,:,i)+u(:,:,i+1));
    tmp2 = (tmp-mean(tmp(:))).^2;
    ups(i) = mean(tmp2(:));
    tmp = 0.5.*(v(:,:,i)+v(:,:,i+1));
    tmp2 = (tmp-mean(tmp(:))).^2;
    vps(i) = mean(tmp2(:));
end
% save data
matfile = [outDir '/wpsdata.mat'];
save(matfile, 'wps', 'ups', 'vps', 'zw',...
    'indhb', 'ind1', 'ind2', 'ind3', 'ind4');


%% functions

function [xx, yy] = calSpectrum(var, mi, Lx)
    % mi: level of zw

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
    q = sqrt(0.5.*(up.^2+vp.^2+wp.^2));

    % Calculate spectra
    [xx.ew, yy.ew, ~] = spectra2D(wp, Lx);
    [xx.eu, yy.eu, ~] = spectra2D(up, Lx);
    [xx.ev, yy.ev, ~] = spectra2D(vp, Lx);
    [xx.ek, yy.ek, ~] = spectra2D(q, Lx);
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

function rl = refk(k,ylims)
    global cgray
    rl = line([k,k], ylims);
    rl.Color = cgray;
    rl.LineStyle = '--';
    rl.LineWidth = 0.75;
end