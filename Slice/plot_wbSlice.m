close all; clear variables;
% load shared script
path(path,'../share');

%% data info
leg = 1;    % 0 if first leg; 1 if second leg
ic = 5;
namelist = {'R8_BF05WD05WV00_ST01_ens03',...
            'R8_BF05WD05WV12_ST01_ens03',...
            'R8_BF1hWD00WV00_ST00_ens03',...
            'R8_BF05WD05WV00_fL00_ens03',...
            'R8_BF05WD05WV12_fL00_ens03'};
tstartleg = {'028800','032400'};
tendleg   = {'032400','036000'};
startlist = {tstartleg{leg},tstartleg{leg},'021600','021600','021600'};
endlist   = {tendleg{leg},tendleg{leg},'028800','025200','025200'};
tstart_hist = [28801, 28801, 1, 21601, 21601];
tend_hist = [36001, 36001, 28801, 25201, 25201];
vmaxlist = [5, 5, 5, 5, 5];
vmax = vmaxlist(ic);
casename = namelist{ic};
tstart = startlist{ic};
tend = endlist{ic};
% initial mixed layer depth
h0 = 42;
% indices of saved x-y slices
if ic == 3
    izi = [2, 6, 44, 77, 78, 79, 80, 81, 82];
elseif ic == 4 || ic == 5
    izi = [2, 6, 41, 68, 69, 70, 71, 72, 73];
else
    izi = [2, 6, 40, 63, 64, 65, 66, 67, 68];
end
idxentr = [5, 7, 7, 5, 7];
idx0 = idxentr(ic); % entrainment level
ns = 1; % time step to start
l_update_data = 1;    % 1 if update data, 0 if read from saved workspace
l_savegif = 0;  % 1 if save gif
l_savepng = 1;  % 1 if save png for all frames
c_gray = [0.5, 0.5, 0.5];
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
figname = [outDir '/wbSlice.gif'];
[fdir, fname, ~] = fileparts(figname);
if l_savepng
    pngdir = [outDir '/png'];
    system(['mkdir -p ' pngdir]);
end

% set up profile data parameters
pflDir = [dataRootDir '/hist'];
fPrefix = 'his.mp.vis';
% creating object ncarlesPflData
f = ncarlesPflData('caseName', casename,...
                   'dataDir',  pflDir,...
                   'fPrefix',  fPrefix,...
                   'tStart',   tstart_hist(ic),...
                   'tEnd',     tend_hist(ic));

%% load data
% also read the buoyancy above an below the entrainment
% layer for the flux limiter
idx0m1 = idx0-1;
idx0p1 = idx0+1;
idx0p2 = idx0+2;
tmp = ncread([dataDir casename '/' filename_xy],'w');
wxy0 = squeeze(tmp(:,:,idx0,:));
tmp = ncread([dataDir casename '/' filename_xy],'t').*f.alpha.*f.g;
bxy0 = squeeze(tmp(:,:,idx0,:));
bxy0m1 = squeeze(tmp(:,:,idx0m1,:));
bxy0p1 = squeeze(tmp(:,:,idx0p1,:));
bxy0p2 = squeeze(tmp(:,:,idx0p2,:));
x = ncread([dataDir casename '/' filename_xy],'x');
y = ncread([dataDir casename '/' filename_xy],'y');
zw = ncread([dataDir casename '/' filename_xz],'zw');
time = ncread([dataDir casename '/' filename_xy],'time');
[nx, ny, nt] = size(wxy0);
Lx = x(1)+x(end);
dx = Lx/nx;
% perturbation
wp = wxy0-mean(mean(wxy0,1),2);
bp = bxy0-mean(mean(bxy0,1),2);
% calculate the entrainment buoyancy fluc, apply flux limiter
tmp = max(-wxy0, 0).*(bxy0+fluxLimiter(bxy0p1, bxy0, bxy0m1));
tmp2 = min(-wxy0, 0).*(bxy0p1+fluxLimiter(bxy0, bxy0p1, bxy0p2));
wpbp = -tmp-tmp2;
wpbpm = squeeze(mean(mean(wpbp,1),2)); % mean entrainment flux

% profile
its = f.getTimeIndex(time(1));
ite = f.getTimeIndex(time(end));
tmp = f.getBLDepth('max_Nsquared');
hb_ts = tmp(its:ite);
hb = mean(hb_ts);
[~,ind_hb] = min(abs(zw-hb));

% mean profile
tmpv = f.getProfiles('second',[time(1), time(end)]);
wbmean = tmpv.wtle.*f.g.*f.alpha;

% wb profile time series
ind_wb = floor(ind_hb*1.2);
tmp = f.getVar('wtle');
wbpflts = squeeze(tmp(1:ind_wb,:,its:ite)).*f.g.*f.alpha;

tmp = f.getVar('utau');
utau = tmp(end);
wgt = utau.^3./abs(hb_ts);
% wgtm = utau.^3./abs(hb);
wgtm = std(wpbp(:));
% wgt = f.getVar('wtsfc')*ones(size(hb_ts)).*f.g.*f.alpha;

%% time loop
dlt = 0.1;  % delay time
ii = 1; % frame counter

% setup figure
fig=figure;
fig.Units = 'inches';
fig.Position = [1 1 8 8];
fig.Color = 'white';

ll = [0.05, 0.80, 0.08, 0.57];
bb = [0.40, 0.42, 0.08, 0.08];
ww = [0.70, 0.15, 0.38, 0.38];
hh = [0.55, 0.53, 0.24, 0.24];

% plot figure
for it = ns:nt
    % entrainment buoyancy flux map
    ax1 = subplot('position',[ll(1) bb(1) ww(1) hh(1)]);
    wpbpnm = wpbp(:,:,it)'./wgtm;
    pcolor(x,y,wpbpnm); shading flat;
%     pcolor(x,y,-wpbp(:,:,it)'); shading flat;
    colorbar;
%     colormap(b2r(-1e6,1e6));
%     caxis([-1e6,1e6]);
    colormap(b2r(-vmax,vmax));
    caxis([-vmax,vmax]);
    daspect([1,1,1]);
    xlim([0,Lx]);
    ylim([0,Lx]);
    xlabel('$x$ (m)', 'Interpreter', 'latex');
    ylabel('$y$ (m)', 'Interpreter', 'latex');

    % vertical profile
    ax2 = subplot('position',[ll(2) bb(2) ww(2) hh(2)]);
    plot(wbmean(1:ind_wb)./wgt(it), zw(1:ind_wb)./abs(hb), '-',...
        'Color', c_gray, 'LineWidth', 2);
    hold on;
    plot(wbpflts(:,it)./wgt(it), zw(1:ind_wb)./abs(hb_ts(it)), '-k',...
        'LineWidth', 1);
    plot(wpbpm(it)./wgt(it), zw(izi(idx0))./abs(hb_ts(it)), '+k',...
        'LineWidth', 1);
    xlims = [-1.5,0.5];
%     xlims = [-0.5,1];
    ylims = [-1.2,0];
    xlim(xlims);
    ylim(ylims);
    rl=line([0, 0], ylims);
    rl.Color = 'k';
    rl.LineStyle = ':';
    xlabel('$\overline{w''b''} h_\mathrm{b}/u_*^3$','Interpreter', 'latex');
    ylabel('$z/h_b$','Interpreter', 'latex');
    hold off;
    
    % pdf
    ax3 = subplot('position',[ll(3) bb(3) ww(3) hh(3)]);
    [histy, histedge] = histcounts(wpbpnm(:));
    histx = 0.5.*(histedge(1:end-1)+histedge(2:end));
    plot(histx,log10(histy),'-k','Linewidth', 1.5);
    xlabel('${w''b''} h_\mathrm{b}/u_*^3$', 'Interpreter', 'latex');
    ylabel('$\log(N)$', 'Interpreter', 'latex');
    xlim([-3*vmax, 3*vmax]);
    ylim([0, 4]);
    
    % spectrum
    ax4 = subplot('position',[ll(4) bb(4) ww(4) hh(4)]);
    [specx, specy, km] = spectra2D(wpbpnm, Lx);
    specwgt = specy(2);
    p4 = plot(specx,specy./specwgt, '-k', 'LineWidth', 1.5);
    hold on;
    xlim([2*pi/Lx, 2*pi/dx/3]);
    ylim([1e-1, 1e2]);
    ylims = ylim;
    plot([km, km], ylim, '--', 'Color', c_gray, 'LineWidth', 1);
    xlabel('$k$ (m$^{-1}$)', 'Interpreter', 'latex');
    ylabel('$E(k)/E(k_0)$', 'Interpreter', 'latex');
    set(gca,'xscale', 'log');
    set(gca,'yscale', 'log');
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
