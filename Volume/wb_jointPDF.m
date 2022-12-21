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
%
figure;
plot(wbpflts, zw, 'r-')
hold on;
plot(wpbpm_rlim, zw, 'k--')

figure;
dat1 = wp(:,:,idx_p9hb);
dat2 = bp_imp(:,:,idx_p9hb);
wgt1 = std(dat1(:));
wgt2 = std(dat2(:));
xlims = [-8, 8];
ylims = [-8, 8];
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
