close all; clear variables;

%% flags
l_save_fig = 0;

%% Parameters
path(path,'../share');
get_dataRootDir;
dataDir = [dataRootDir '/rest/'];
casename = 'R8_BF05WD05WV12_ST01_ens02';
% casename = 'R8_BF1hWD00WV00_ST00_ens03';
% [s0, r0] = system(['ls ' dataDir casename ]);
% if ~s0
%     tttttt = r0(1:6);
% end
tttttt = '028801';
% tttttt = '025201';
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

% get boundary layer depth

N2 = zeros(1,nnz);
txym = zeros(1,nnz);
for i=1:nnz
    tmp = squeeze(var(4,:,:,i));
    txym(i) = mean(tmp(:));
end
N2(1:end-1) = (txym(2:end)-txym(1:end-1))./dz;
N2(end) = N2(end-1);
[~,idx_hb] = max(N2(:));

% set up profile data parameters
pflDir = [dataRootDir '/hist'];
fPrefix = 'his.mp.vis';
% creating object ncarlesPflData
f = ncarlesPflData('caseName', casename,...
                   'dataDir',  pflDir,...
                   'fPrefix',  fPrefix,...
                   'tStart',   1,...
                   'tEnd',     28801);

% resolved buoyancy flux profile
tmp = f.getVar('wtle');
wbpflts = squeeze(tmp(:,:,end)).*batag;

% vertical velocity
w = squeeze(var(3,:,:,:));
wxym = mean(mean(w,1),2);
wp = w-wxym;

% buoyancy and buoyancy flux with flux limiter
b = squeeze(var(4,:,:,:)).*batag;
bxym = mean(mean(b,1),2);
bp = b-bxym;
bpw = bp;
bpw(:,:,1:end-1) = 0.5.*(bp(:,:,1:end-1)+bp(:,:,2:end));
wpbp = wp.*bpw;
wpbpm = squeeze(mean(mean(wp.*bpw,1),2));
idx = 2:nnz-2;
idxp1 = idx+1;
idxp2 = idx+2;
idxm1 = idx-1;
tmp = zeros(size(w));
tmp2 = zeros(size(w));
tmp(:,:,idx) = max(-w(:,:,idx), 0).*(b(:,:,idx)...
    +fluxLimiter(b(:,:,idxp1), b(:,:,idx), b(:,:,idxm1)));
tmp2(:,:,idx) = min(-w(:,:,idx), 0).*(b(:,:,idxp1)...
    +fluxLimiter(b(:,:,idx), b(:,:,idxp1), b(:,:,idxp2)));
wpbp_rlim = -tmp-tmp2;
wpbp_rlim(:,:,1) = wpbp(:,:,1);
wpbp_rlim(:,:,end-2:end) = wpbp(:,:,end-2:end);
wpbpm_rlim = squeeze(mean(mean(wpbp_rlim,1),2));

% implied buoyancy perturbation
bp_imp = wpbp_rlim./wp;
bp_imp = bp_imp-mean(mean(bp_imp,1),2);

%% u, v, w near hb: z = 0.9*hb

idx_p9hb = round(idx_hb*0.9);
pre_str = '0p9h';
depth3 = zw(idx_p9hb);
fprintf('Depth: %g\n', depth3);

c_gray = [0.5, 0.5, 0.5];
dat = wpbp_rlim(:,:,idx_p9hb);
wgt = std(dat(:));
% plot figure w
newFigure(l_save_fig);
pcolor(x,y,dat'./wgt); shading flat; colorbar;
caxis([-5, 5]);
xlabel('$x$ (m)', 'Interpreter', 'latex');
ylabel('$y$ (m)', 'Interpreter', 'latex');
daspect([1 1 1]);

if l_save_fig
    figname = 'snapshot_wb.fig';
    saveas(gcf, figname, 'fig');
    postProcessFig(figname, 6, 4);
end

% Calculate spectra
% [ew_x, ew_y, ~] = coSpectra2D(wp(:,:,idx_p9hb), bp_imp(:,:,idx_p9hb), Lx);
dat = squeeze(wpbp_rlim(:,:,idx_p9hb));
[ew_x, ew_y, ~] = spectra2D(dat, Lx);

newFigure(l_save_fig);
hold on;
wgt = ew_y(2);
p3 = plot(ew_x,ew_y./wgt, '-.k', 'LineWidth', 1.5);
ff = 10./(ew_x(2)).^(-5/3);
plot(ew_x, ff*ew_x.^(-5/3),'-',...
    'Color', c_gray, 'LineWidth', 0.75);
nnx = numel(x);
dx = Lx/nnx;
xlim([2*pi/Lx, 2*pi/dx/3]);
% ylim([1e-3, 1e1]);
xlabel('$k$ (m$^{-1}$)', 'Interpreter', 'latex');
ylabel('$E_i(k)/E(k_0)$', 'Interpreter', 'latex');
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
% lg = legend([p1, p2, p3], '$E_u$','$E_v$','$E_w$');
% lg.Location = 'NorthEast';
% lg.Interpreter = 'latex';

if l_save_fig
    figname = 'cospectrum_wb.fig';
    saveas(gcf, figname, 'fig');
    postProcessFig(figname, 6, 4);
end

%
figure;
dat1 = wp(:,:,idx_p9hb);
dat2 = bp_imp(:,:,idx_p9hb);
wgt1 = std(dat1(:));
wgt2 = std(dat2(:));
% plot(dat1(:),dat2(:),'.');
[n, xedge, yedge] = histcounts2(dat1./wgt1,dat2./wgt2,100,'Normalization','pdf');
xx = 0.5.*(xedge(1:end-1)+xedge(2:end));
yy = 0.5.*(yedge(1:end-1)+yedge(2:end));
pcolor(xx, yy, log10(n'));
shading flat;
cb = colorbar;
cb.Label.String = 'log(PDF)';
% cb.Label.Interpreter = 'latex';
c = pink;
c = flipud(c);
colormap(c);
xlabel('wp');
ylabel('bp');
xlims = [-4, 4];
ylims = [-4, 4];
xlim(xlims);
ylim(ylims);
rl = line(xlims,[0,0]);
rl.Color = 'k';
rl = line([0,0],ylims);
rl.Color = 'k';

figure;
dat1 = wp(:,:,idx_p9hb);
dat2 = bp_imp(:,:,idx_p9hb);
wgt1 = std(dat1(:));
wgt2 = std(dat2(:));
xlims = [-4, 4];
ylims = [-4, 4];
[hst, xi, yi] = jointDist(dat1./wgt1, xlims(1), xlims(2), dat2./wgt2, ylims(1), ylims(2));
plot_dist_4p(hst,xi,yi);
xlabel('wp');
ylabel('bp');
xlim(xlims);
ylim(ylims);
rl = line(xlims,[0,0]);
rl.Color = 'k';
rl = line([0,0],ylims);
rl.Color = 'k';
%% save data
% if l_save_fig
% 
%     if ~exist(outDir, 'dir')
%         mkdir(outDir);
%     end
%     system(['mv snapshot* ' outDir]);
%     system(['mv spectrum* ' outDir]);
%     system(['mv profile* ' outDir]);
% end


%% functions

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

function [e_x, e_y, km] = coSpectra2D(dat1, dat2, Lx)
% coSpectra2D calculates the 2D coincident spectrum of dat1 and dat2

    if ~all(size(dat1)==size(dat2))
        error('Dimension of two dateset should be the same\n');
    end
    nsize = size(dat1);
    nnx = nsize(1);
    nny = nsize(2);
    fw1 = fftshift(fft2(dat1));
    fw2 = fftshift(fft2(dat2));
    pw = fw1.*conj(fw2);
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
    e_y = real(ew(1:nnx/2));
    km = km./sum(pw(:)).*sca;
end

function h = plot_dist_4p(hst,xi,yi)
% plot_dist_4p plot the highest 50%, 75%, 90% and 95%
%   centered distribution

    % find the isolines of pdf that enclose the area in which the
    % total probability is 50%, 75%, 90% and 95%
    hsum = sum(hst(:));
    hlist = sort(hst(:),'descend')./hsum;
    hcum = cumsum(hlist);
    vl = [0.5, 0.75, 0.9, 0.95];
    nv = numel(vl);
    vlev = zeros(1,nv);
    for i=1:nv
        [~,ind] = min(abs(hcum-vl(i)));
        vlev(i) = hlist(ind);
    end
    pdfData = hst./hsum;
    % plot log10(pdfData) to make use of the full colarbar
    pdfData(pdfData==0) = 1e-12;    % avoid -inf for log10(0)
    [~,h] = contourf(xi,yi,log10(pdfData'));
    h.LevelListMode = 'manual';
    h.LevelList = log10(vlev);
    h.ShowText = 'on';
    h.TextList = vl;

    caxis([log10(vlev(end)) log10(vlev(1))]);
    tmp = colormap;
    inds = [2, 2, 15, 15, 15, 15, 15, 15, 35, 35, 35, 35, 35, 35, 60, 60];
    my_colorm = tmp(inds,:);
    colormap(my_colorm);
end
