 close all; clear variables;
% This script plots scatters of all cases in Li and Fox-Kemper, 2017.
% Figures 2a, 5-8, and 9a.

% flag for saving figures. 0: do not save, 1: save
l_save_fig = 1;

%% cases
BF_str = ['05';'10';'25';'50';'1h';'2h';'3h';'5h'];
WD_str = ['05';'08';'10'];
WV_str = ['00';'01';'02';'03';'04';'12';'13';'14'];
nsize = size(BF_str);
nBF = nsize(1);
nsize = size(WD_str);
nWD = nsize(1);
nsize = size(WV_str);
nWV = nsize(1);

%% load data
dat = load('LF17_Data_Scatter.mat');
utau = dat.utau;
us = dat.us;
b0 = dat.b0;
hb = dat.hb;
wpsm = dat.wpsm;
tkem = dat.tkem;
lasl = dat.lasl;
stat_wb_mean = dat.stat_wb_mean;
stat_wb_median = dat.stat_wb_median;
stat_wb_p25 = dat.stat_wb_p25;
stat_wb_p75 = dat.stat_wb_p75;

% Langmuir number
la = sqrt(utau./us);
% h/LL
hL = b0.*hb./utau.^2./us;

% fitting for entrainment buoyancy flux
wgt = utau.^3./hb;
xf1 = b0./wgt;
xf2 = lasl.^(-2);
ydat0 = -(stat_wb_mean)./wgt;
yf1 = ydat0;
xf1(isnan(xf1)) = [];
xf2(isnan(xf2)) = [];
yf1(isnan(yf1)) = [];
[fg, rg] = fit([xf1(:),xf2(:)],yf1(:),'poly11');

% color gray
c_gray = [0.5 0.5 0.5];

%% Figure 2a: parameter space
newFigure(l_save_fig);

% Reproduce Fig. 3 in Belcher et al., 2012
xlims = [0.1,10];   % La
ylims = [1e-3,1e3]; % h/LL

% color contour (background)
nx = 500;
ny = 500;
xx = logspace(log10(xlims(1)),log10(xlims(2)),nx);
yy = logspace(log10(ylims(1)),log10(ylims(2)),ny);
zz = zeros(nx,ny);
zz1 = zz;
zz2 = zz;
zz3 = zz;
for i = 1:nx
for j = 1:ny
    zz1(j,i) = 2.*(1-exp(-0.5.*xx(i)));
    zz2(j,i) = 0.22.*xx(i).^(-2);
    zz3(j,i) = 0.3.*xx(i).^(-2).*yy(j);
end
end
zz = zz1 + zz2 + zz3;

hold on;
colormap('summer');
[~,p] = contourf(xx,yy,log10(zz));
p.LevelList = [-0.1, 0, 0.1, 0.25, 0.5, 1, 2, 3];
p.ShowText = 'on';
p.LineColor = [0.3 0.3 0.3];

[~,p] = contour(xx,yy,zz1./zz);
p.LevelList = 0.9;
p.LineWidth = 1.5;
p.LineColor = 'k';
[~,p] = contour(xx,yy,zz2./zz);
p.LevelList = 0.9;
p.LineWidth = 1.5;
p.LineColor = 'k';
[~,p] = contour(xx,yy,zz3./zz);
p.LevelList = 0.9;
p.LineWidth = 1.5;
p.LineColor = 'k';

% white contour (parameter space in CESM)
hstdat = load('wstar3ustar2uS_LaTurb.mat');
plot_dist_4p_l(hstdat.hst,hstdat.xi,hstdat.yi);

% scatter data points
ids3w = 2:5;
ids3s = 6:nWV;
xdat0 = la;
ydat0 = hL;
zdat0 = log10(eps./utau.^3.*hb);
% group 1, wind waves
xdat = xdat0(:,1,ids3w);
ydat = ydat0(:,1,ids3w);
plot(xdat(:),ydat(:),'ow','LineWidth',1.0);
xdat = xdat0(:,2,ids3w);
ydat = ydat0(:,2,ids3w);
plot(xdat(:),ydat(:),'sw','LineWidth',1.0);
xdat = xdat0(:,3,ids3w);
ydat = ydat0(:,3,ids3w);
plot(xdat(:),ydat(:),'dw','LineWidth',1.0);
% group 2, swells
xdat = xdat0(:,1,ids3s);
ydat = ydat0(:,1,ids3s);
scatter(xdat(:),ydat(:),'ow','filled');
xdat = xdat0(:,2,ids3s);
ydat = ydat0(:,2,ids3s);
scatter(xdat(:),ydat(:),'sw','filled');
xdat = xdat0(:,3,ids3s);
ydat = ydat0(:,3,ids3s);
scatter(xdat(:),ydat(:),'dw','filled');
set(gca,'xscale','log');
set(gca,'yscale','log');
xlabel('$La_\mathrm{t}$','Interpreter','latex');
ylabel('$h_\mathrm{b}/L_\mathrm{L}$','Interpreter','latex');
xlim([0.1 10]);
ylim([1e-3,1e3]);

tx = text(0.11, 2e-3,'Langmuir','Interpreter','Latex');
tx.BackgroundColor = 'w';
tx = text(3, 2e-2,'Wind','Interpreter','Latex');
tx.BackgroundColor = 'w';
tx = text(0.13, 1e2,'Convection','Interpreter','Latex');
tx.BackgroundColor = 'w';

if (l_save_fig)
    figname = 'fig2a_sct_ParSpace.pdf';
    postProcessFig(gcf, figname, 5, 5);
end


%% Figure 5a: wps vs surface buoyancy flux (w*^2)
newFigure(l_save_fig);

% set data
xdat0 = (b0.*hb).^(2/3)./utau.^2;
ydat0 = wpsm(:,:,:,1)./utau.^2;
xlabel_str = '${w^*}^2/{u^*}^2$';
ylabel_str = '$\langle \overline{w''^2}\rangle_{h_\mathrm{m}}/{u^*}^2$';
xlog = 1;
ylog = 0;
xlims = [0.1, 30];
ylims = [0, 6];

% --- ref line
hold on;
nx = 10000;
dx = (xlims(2)-xlims(1))/nx;
x = xlims(1):dx:xlims(2)*2;
y = 0.3.*x;
y2 = x.*0+0.6;
plot(x,y,'-','LineWidth',1.0,'Color',c_gray);
plot(x,y2,'-','LineWidth',1.0,'Color',c_gray);
text(6, 5,'$y=0.3x$','Interpreter','Latex');
text(10, 1,'$y=0.6$','Interpreter','Latex');

% plot scatter 
plotScatterAll(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims);

% save figure
if (l_save_fig)
    figname = 'fig5a_sct_wps_wstar.pdf';
    postProcessFig(gcf, figname);
end


%% Figure 5b: wps vs La^{-2} (us0)
newFigure(l_save_fig);

xdat0 = la.^-2;
ydat0 = wpsm(:,:,:,1)./utau.^2;
xlabel_str = '$La_\mathrm{t}^{-2}=u^\mathrm{S}_0/{u^*}$';
ylabel_str = '$\langle \overline{w''^2}\rangle_{h_\mathrm{m}}/{u^*}^2$';
xlog = 0;
ylog = 0;
xlims = [0, 16];
ylims = [0, 6];

% --- ref line
hold on;
nx = 10000;
dx = (xlims(2)-xlims(1))/nx;
x = xlims(1):dx:xlims(2)*2;
y = 0.6.*(1+(3.1).^(-2).*x+(5.7).^(-4).*x.^2);
y2 = 0.64.*(1+0.098.*x);
p1 = plot(x,y,'-','LineWidth',0.8,'Color',c_gray);
p2 = plot(x,y2,'--','LineWidth',0.8,'Color',c_gray);

xcl = squeeze(xdat0(1,:,:));
ycl = squeeze(ydat0(1,:,:));
ycl1 = 0.6.*(1+(3.1).^(-2).*xcl+(5.7).^(-4).*xcl.^2);
ycl2 = 0.64.*(1+0.098.*xcl);
rmse11 = sqrt(mean((ycl1(:)-ycl(:)).^2));
rmse12 = sqrt(mean((ycl2(:)-ycl(:)).^2));
str_rmse11 = sprintf('%4.2f',rmse11);
str_rmse12 = sprintf('%4.2f',rmse12);
eq_str = ['$0.6(1+(3.1 La_\mathrm{t})^{-2}+(5.7 La_\mathrm{t})^{-4})\;$'...
    ' \textbf{[' str_rmse11 ']}'];
eq2_str = ['$0.64(1+0.098 La_\mathrm{t}^{-2})\;$'...
    ' \textbf{[' str_rmse12 ']}'];
lg = legend([p1, p2],eq_str,eq2_str);
lg.Location = 'NorthEast';
lg.Interpreter = 'Latex';
lg.Color = 'none';
lg.AutoUpdate = 'off';
ids2 = 1;
ids3 = 1:nWV;
for ids1=3:nBF
    xx = squeeze(xdat0(ids1,ids2,ids3));
    yy = squeeze(ydat0(ids1,ids2,ids3));
    xx(isnan(xx)) = [];
    yy(isnan(yy)) = [];
    [xx2, inds] = sort(xx);
    yy2 = yy(inds);
    px = plot(xx2(:),yy2(:),':','Clipping','off');
    px.Color = c_gray;
    px.LineWidth = 0.8;
end

% plot scatter
plotScatterAll(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims);
            
if (l_save_fig)
    figname = 'fig5b_sct_wps_La.pdf';
    postProcessFig(gcf, figname);
end


%% Figure 5c: wps vs LaSL^{-2} (ussl)
newFigure(l_save_fig);

xdat0 = lasl.^(-2);
ydat0 = wpsm(:,:,:,1)./utau.^2;
xlabel_str = '$La_\mathrm{SL}^{-2}=u^\mathrm{S}_\mathrm{SL}/{u^*}$';
ylabel_str = '$\langle \overline{w''^2}\rangle_{h_\mathrm{m}}/{u^*}^2$';
xlog = 0;
ylog = 0;
xlims = [0, 5];
ylims = [0, 6];

% --- ref line
hold on;
nx = 10000;
dx = (xlims(2)-xlims(1))/nx;
x = xlims(1):dx:xlims(2)*2;
y = 0.6.*(1+(1.5).^(-2).*x+(5.4).^(-4).*x.^2);
y2 = 0.398+0.48.*x.^(2/3);
y3 = 0.64+3.5.*exp(-2.69.*x.^(-0.5));
p1 = plot(x,y,'-','LineWidth',0.8,'Color',c_gray);
p2 = plot(x,y2,'--','LineWidth',0.8,'Color',c_gray);
p3 = plot(x,y3,'-.','LineWidth',0.8,'Color',c_gray);

xcl = squeeze(xdat0(1,:,:));
ycl = squeeze(ydat0(1,:,:));

ycl1 = 0.6.*(1+(1.5).^(-2).*xcl+(5.4).^(-4).*xcl.^2);
ycl2 = 0.398+0.48.*xcl.^(2/3);
ycl3 = 0.64+3.5.*exp(-2.69.*xcl.^(-0.5));

rmse21 = sqrt(mean((ycl1(:)-ycl(:)).^2));
rmse22 = sqrt(mean((ycl2(:)-ycl(:)).^2));
rmse23 = sqrt(mean((ycl3(:)-ycl(:)).^2));
str_rmse21 = sprintf('%4.2f',rmse21);
str_rmse22 = sprintf('%4.2f',rmse22);
str_rmse23 = sprintf('%4.2f',rmse23);
eq_str = ['$0.6 (1 + (1.5 La_\mathrm{SL})^{-2}+(5.4 La_\mathrm{SL})^{-4})\;$'...
    ' \textbf{[' str_rmse21 ']}'];
eq2_str = ['$0.398 + 0.480 La_\mathrm{SL}^{-4/3}\;$'...
    ' \textbf{[' str_rmse22 ']}'];
eq3_str = ['$0.640+3.50\exp(-2.69 La_\mathrm{SL})\;$'...
    ' \textbf{[' str_rmse23 ']}'];
lg = legend([p1, p2, p3],eq_str,eq2_str,eq3_str);
lg.Location = 'NorthEast';
lg.Interpreter = 'Latex';
lg.AutoUpdate = 'off';
ids2 = 1;
ids3 = 1:nWV;
for ids1=3:nBF
    xx = squeeze(xdat0(ids1,ids2,ids3));
    yy = squeeze(ydat0(ids1,ids2,ids3));
    xx(isnan(xx)) = [];
    yy(isnan(yy)) = [];
    [xx2, inds] = sort(xx);
    yy2 = yy(inds);
    px = plot(xx2(:),yy2(:),':','Clipping','off');
    px.Color = c_gray;
    px.LineWidth = 0.8;
end

% plot scatter
plotScatterAll(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims);

if (l_save_fig)
    figname = 'fig5c_sct_wps_LaSL.pdf';
    postProcessFig(gcf, figname);
end


%% Figure 6a: vke / tke, wstar
newFigure(l_save_fig);

% set data
xdat0 = (b0.*hb).^(2/3)./utau.^2;
ydat0 = 0.5.*wpsm(:,:,:,1)./tkem(:,:,:,1);

xlabel_str = '${w^*}^2/{u^*}^2$';
ylabel_str = 'VKE/TKE';
xlims = [0.1,30];
ylims = [0,0.5];
xlog = 1;
ylog = 0;
hold on;
rl = refline(0,1/3);
rl.Color = c_gray;
rl.LineWidth = 1.0;
rl.XData = xlims;

% plot scatter
plotScatterAll(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims);

% save figure
if (l_save_fig)
    figname = 'fig6a_sct_tke_wps_wstar.pdf';
    postProcessFig(gcf, figname);
end


%% Figure 6b: vke / tke, LaSL
newFigure(l_save_fig);

% set data
xdat0 = lasl.^(-2);
ydat0 = 0.5.*wpsm(:,:,:,1)./tkem(:,:,:,1);

xlabel_str = '$La_\mathrm{SL}^{-2}=u^\mathrm{S}_\mathrm{SL}/{u^*}$';
ylabel_str = 'VKE/TKE';
xlims = [0,5];
ylims = [0,0.5];
xlog = 0;
ylog = 0;
hold on;
rl = refline(0,1/3);
rl.Color = c_gray;
rl.LineWidth = 1.0;
rl.XData = xlims;

% plot scatter
plotScatterAll(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims);

% save figure
if (l_save_fig)
    figname = 'fig6b_sct_tke_wps_LaSL.pdf';
    postProcessFig(gcf, figname);
end


%% Figure 6c: vke / tke, La
newFigure(l_save_fig);

% set data
xdat0 = la.^(-2);
ydat0 = 0.5.*wpsm(:,:,:,1)./tkem(:,:,:,1);

xlabel_str = '$La_\mathrm{t}^{-2}=u^\mathrm{S}_0/{u^*}$';
ylabel_str = 'VKE/TKE';
xlims = [0,16];
ylims = [0,0.5];
xlog = 0;
ylog = 0;
hold on;
rl = refline(0,1/3);
rl.Color = c_gray;
rl.LineWidth = 1.0;
rl.XData = xlims;

% plot scatter
plotScatterAll(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims);

% save figure
if (l_save_fig)
    figname = 'fig6c_sct_tke_wps_La.pdf';
    postProcessFig(gcf, figname);
end


%% Figure 7a: entrainment vs surface buoyancy flux
newFigure(l_save_fig);

% set data
xdat0 = b0.*hb./utau.^3;
wgt = wpsm(:,:,:,1).^(3/2)./hb;
ydat0 = -stat_wb_mean./wgt;
ydat1 = -stat_wb_median./wgt;
err0_n = -(stat_wb_p75-stat_wb_median)./wgt;
err0_p = -(stat_wb_median-stat_wb_p25)./wgt;
xlabel_str = '${w^*}^3/{u^*}^3$';
ylabel_str = ['$-\overline{w''b''}_\mathrm{e} h_\mathrm{b}/'...
              '\langle \overline{w''^2}\rangle_{h_\mathrm{m}}^{3/2}$'];
xlog = 1;
ylog = 0;
ylims = [0.1,1.5];
xlims = [0.05,100];

% plot scatter
plotScatterAll(xdat0,ydat0,xlabel_str,ylabel_str,...
                xlog,ylog,xlims,ylims,...
                ydat1,err0_n,err0_p);

% save figure
if (l_save_fig)
    figname = 'fig7a_sct_wb2_wstar.pdf';
    postProcessFig(gcf, figname);
end


%% Figure 7b: entrainemnt vs LaSL^{-2}
newFigure(l_save_fig);

% set data
xdat0 = lasl.^(-2);
wgt = wpsm(:,:,:,1).^(3/2)./hb;
ydat0 = -(stat_wb_mean)./wgt;
ydat1 = -(stat_wb_median)./wgt;
err0_n = -(stat_wb_p75-stat_wb_median)./wgt;
err0_p = -(stat_wb_median-stat_wb_p25)./wgt;
xlabel_str = '$La_\mathrm{SL}^{-2}=u^\mathrm{S}_\mathrm{SL}/{u^*}$';
ylabel_str = ['$-\overline{w''b''}_\mathrm{e} h_\mathrm{b}/'...
              '\langle \overline{w''^2}\rangle_{h_\mathrm{m}}^{3/2}$'];
xlog = 0;
ylog = 0;
xlims = [0,5];
ylims = [0.1,1.5];

% plot scatter
plotScatterAll(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                ydat1,err0_n,err0_p);

% save figure
if (l_save_fig)
    figname = 'fig7b_sct_wb2_LaSL.pdf';
    postProcessFig(gcf, figname);
end


%% Figure 8a: entrainment vs surface buoyancy flux
newFigure(l_save_fig);

% set data
xdat0 = b0.*hb./utau.^3;
wgt = utau.^3./hb;
ydat0 = -stat_wb_mean./wgt;
ydat1 = -stat_wb_median./wgt;
err0_n = -(stat_wb_p75-stat_wb_median)./wgt;
err0_p = -(stat_wb_median-stat_wb_p25)./wgt;
xlabel_str = '${w^*}^3/{u^*}^3$';
ylabel_str = '$-\overline{w''b''}_\mathrm{e} h_\mathrm{b}/{u^*}^3$';
xlog = 1;
ylog = 1;
ylims = [0.1,20];
xlims = [0.05,100];

% --- ref line
hold on;
ids1 = 1:nBF;
ids3 = 1;
xdat = xdat0(ids1,:,ids3);
ydat = ydat0(ids1,:,ids3);
xdat(isnan(xdat)) = [];
ydat(isnan(ydat)) = [];
xx = xdat;
yy = ydat;
[ff,rr] = fit(xx(:),yy(:),'poly1');
ff1 = ff;
rr1 = rr;
nx = 10000;
dx = (xlims(2)-xlims(1))/nx;
x = 0:dx:xlims(2)*2;
y = ff.p1.*x+ff.p2;
y2 = fg.p10.*x+fg.p00;
plot(x,y,'-','LineWidth',1.0,'Color',c_gray);
plot(x,y2,'--','LineWidth',1.0,'Color',c_gray);
eq_str = sprintf('$y=%4.2g x + %4.2g$',ff.p1,ff.p2);
text(0.3, 3, eq_str,'Interpreter','Latex');
r2_str = sprintf('$r^2=%4.2f$',rr.rsquare);
text(0.3, 1.8, r2_str,'Interpreter','Latex');

% plot scatter
plotScatterAll(xdat0,ydat0,xlabel_str,ylabel_str,...
                xlog,ylog,xlims,ylims,...
                ydat1,err0_n,err0_p);

% save figure
if (l_save_fig)
    figname = 'fig8a_sct_wb_wstar.pdf';
    postProcessFig(gcf, figname);
end


%% Figure 8b: entrainemnt vs La^{-2}
newFigure(l_save_fig);

% set data
xdat0 = la.^-2;
wgt = utau.^3./hb;
soffset = 0.15;
% soffset = ff1.p2;
yoffset = soffset.*b0./wgt;
ydat0 = -(stat_wb_mean)./wgt-yoffset;
ydat1 = -(stat_wb_median)./wgt-yoffset;
err0_n = -(stat_wb_p75-stat_wb_median)./wgt;
err0_p = -(stat_wb_median-stat_wb_p25)./wgt;
xlabel_str = '$La_\mathrm{t}^{-2}=u^\mathrm{S}_0/{u^*}$';
ylabel_str = ['$(-\overline{w''b''}_\mathrm{e}-'...
              sprintf('%4.2g',soffset)...
              ' B_0) h_\mathrm{b}/{u^*}^3$'];

xlog = 0;
ylog = 0;
xlims = [0,16];
ylims = [-1.5,1.5];

% --- ref line
hold on;
ids = 1:nWV;
xdat = xdat0(1,:,ids);
ydat = ydat0(1,:,ids);
xdat(isnan(xdat)) = [];
ydat(isnan(ydat)) = [];
xx = xdat;
yy = ydat;
[ff,rr] = fit(xx(:),yy(:),'poly1');
ff2 = ff;
rr2 = rr;
nx = 1000;
dx = (xlims(2)-xlims(1))/nx;
x = xlims(1):dx:xlims(2)*2;
y = ff.p1.*x+ff.p2;
plot(x,y,'-k','LineWidth',1.0,'Color',c_gray);
eq_str = sprintf('$y = %4.2g x + %4.2g$',ff.p1,ff.p2);
text(10, -0.5, eq_str,'Interpreter','Latex');
r2_str = sprintf('$r^2=%4.2f$',rr.rsquare);
text(10, -1.0, r2_str,'Interpreter','Latex');

% plot scatter
plotScatterAll(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                ydat1,err0_n,err0_p);
            
% save figure
if (l_save_fig)
    figname = 'fig8b_sct_wb_La.pdf';
    postProcessFig(gcf, figname);
end


%% Figure 8c: entrainemnt vs LaSL^{-2}
newFigure(l_save_fig);

% set data
xdat0 = lasl.^(-2);
wgt = utau.^3./hb;
soffset = 0.15;
yoffset = soffset.*b0./wgt;
ydat0 = -(stat_wb_mean)./wgt-yoffset;
ydat1 = -(stat_wb_median)./wgt-yoffset;
err0_n = -(stat_wb_p75-stat_wb_median)./wgt;
err0_p = -(stat_wb_median-stat_wb_p25)./wgt;
xlabel_str = '$La_\mathrm{SL}^{-2}=u^\mathrm{S}_\mathrm{SL}/{u^*}$';
ylabel_str = ['$(-\overline{w''b''}_\mathrm{e}-'...
              sprintf('%4.2g',soffset)...
              ' B_0) h_\mathrm{b}/{u^*}^3$'];
xlog = 0;
ylog = 0;
xlims = [0,5];
ylims = [-1.5,1.5];

% --- ref line
hold on;
ids1 = 1;
ids3 = 1:nWV;
xdat = xdat0(ids1,:,ids3);
ydat = ydat0(ids1,:,ids3);
xdat(isnan(xdat)) = [];
ydat(isnan(ydat)) = [];
xx = xdat;
yy = ydat;
[ff,rr] = fit(xx(:),yy(:),'poly1');
ff3 = ff;
rr3 = rr;
nx = 1000;
dx = (xlims(2)-xlims(1))/nx;
x = xlims(1):dx:xlims(2)*2;
y = ff.p1.*x+ff.p2;
y2 = fg.p01.*x+fg.p00;
plot(x,y,'-','LineWidth',1.0,'Color',c_gray);
plot(x,y2,'--','LineWidth',1.0,'Color',c_gray);
eq_str = sprintf('$y=%4.2g x+%4.2g$',ff.p1,ff.p2);
text(3, -0.5, eq_str,'Interpreter','Latex');
r2_str = sprintf('$r^2=%4.2f$',rr.rsquare);
text(3, -1.0, r2_str,'Interpreter','Latex');

% plot scatter
plotScatterAll(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                ydat1,err0_n,err0_p);

% save figure
if (l_save_fig)
    figname = 'fig8c_sct_wb_LaSL.pdf';
    postProcessFig(gcf, figname);
end

%% Figure 9a: regime diagram for entrainment
newFigure(l_save_fig);

xlims = [0.1,20];   % Lasl
ylims = [0.01,200]; % -hb/L

c1 = 0.17;
c2 = 0.15;
c3 = 0.083;
nx = 500;
ny = 500;
xx = logspace(log10(xlims(1)),log10(xlims(2)),nx);
yy = logspace(log10(ylims(1)),log10(ylims(2)),ny);
zz = zeros(nx,ny);
zz1 = zz;
zz2 = zz;
zz3 = zz;
for i = 1:nx
for j = 1:ny
    zz1(j,i) = c1;
    zz2(j,i) = c3.*xx(i).^(-2);
    zz3(j,i) = c2.*yy(j);
end
end
zz = zz1 + zz2 + zz3;

% color countour (background)
hold on;
colormap('summer');
[C,p] = contourf(xx,yy,log10(zz));
p.ShowText = 'on';
p.LineColor = [0.3 0.3 0.3];
p.LabelSpacing = 500;
p.TextStep = 0.2;
clabel(C,p,'Color',[0.3 0.3 0.3]);

[~,p] = contour(xx,yy,zz1./zz);
p.LevelList = 0.9;
p.LineWidth = 1.5;
p.LineColor = 'k';
[~,p] = contour(xx,yy,zz2./zz);
p.LevelList = 0.9;
p.LineWidth = 1.5;
p.LineColor = 'k';
[~,p] = contour(xx,yy,zz3./zz);
p.LevelList = 0.9;
p.LineWidth = 1.5;
p.LineColor = 'k';

[~,p] = contour(xx,yy,zz1./zz);
p.LevelList = 0.6;
p.LineWidth = 1.5;
p.LineColor = 'k';
p.LineStyle = '--';
[~,p] = contour(xx,yy,zz2./zz);
p.LevelList = 0.6;
p.LineWidth = 1.5;
p.LineColor = 'k';
p.LineStyle = '--';
[~,p] = contour(xx,yy,zz3./zz);
p.LevelList = 0.6;
p.LineWidth = 1.5;
p.LineColor = 'k';
p.LineStyle = '--';

% white contour (parameter space in CESM)
hstdat = load('wstar3ustar3_LaSL.mat');
plot_dist_4p_l(hstdat.hst,hstdat.xi,hstdat.yi);

% scatter data points
ids3w = 2:5;
ids3s = 6:nWV;
xdat0 = lasl;
ydat0 = b0.*hb./utau.^3;
% group 1, wind waves
xdat = xdat0(:,1,ids3w);
ydat = ydat0(:,1,ids3w);
plot(xdat(:),ydat(:),'ow','LineWidth',1.0);
xdat = xdat0(:,2,ids3w);
ydat = ydat0(:,2,ids3w);
plot(xdat(:),ydat(:),'sw','LineWidth',1.0);
xdat = xdat0(:,3,ids3w);
ydat = ydat0(:,3,ids3w);
plot(xdat(:),ydat(:),'dw','LineWidth',1.0);
% group 2, swells
xdat = xdat0(:,1,ids3s);
ydat = ydat0(:,1,ids3s);
scatter(xdat(:),ydat(:),'ow','filled');
xdat = xdat0(:,2,ids3s);
ydat = ydat0(:,2,ids3s);
scatter(xdat(:),ydat(:),'sw','filled');
xdat = xdat0(:,3,ids3s);
ydat = ydat0(:,3,ids3s);
scatter(xdat(:),ydat(:),'dw','filled');
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('$La_\mathrm{SL}$','Interpreter','latex');
ylabel('$-h_\mathrm{b}/(\kappa L)$','Interpreter','latex');
xlim(xlims);
ylim(ylims);

tx = text(0.12, 0.06,'Langmuir','Interpreter','Latex');
tx.BackgroundColor = 'w';
tx = text(4, 0.06,'Wind','Interpreter','Latex');
tx.BackgroundColor = 'w';
tx = text(3, 20,'Convection','Interpreter','Latex');
tx.BackgroundColor = 'w';

if (l_save_fig)
    figname = 'fig9a_sct_wb_ParSpace.pdf';
    postProcessFig(gcf, figname, 5, 5);
end


%% functions

function h = plotScatterAll(xdat0,ydat0,...
                            xlabel_str, ylabel_str,...
                            xlog,ylog,xlims,ylims,...
                            varargin)
%
    % by default no error bar
    l_errbar = 0;
    % count arguments
    nArgs = length(varargin);
    if nArgs==3
        ydat1 = varargin{1};
        err0_n = varargin{2};
        err0_p = varargin{3};
        l_errbar = 1;
    elseif nArgs>0
        error(['h = plotScatterAll(xdat0,ydat0,'...
               'xlabel_str, ylabel_str,'...
               'xlog,ylog,xlims,ylims,'...
               '[ydat1],[err0_n],[err0_p])']);
    end

    hold on;
    % constants
    [nBF,nWD,nWV] = size(xdat0);
    id3 = 2:nWV;
    id3w = 2:5;
    id3s = 6:nWV;
    id1 = 2:nBF;
    marker_mat = ['o';'s';'d'];
    color_mat = ['r','b','m'];
    c_gray = [0.5 0.5 0.5];

    % plot
    % group 1: no wave, surface cooling > 5 W m^-2
    xdat = xdat0(id1,:,1);
    ydat = ydat0(id1,:,1);
    plot(xdat(:),ydat(:),'xk','LineWidth',1.0,'Clipping','off');
    if l_errbar
        err_neg = err0_n(id1,:,1);
        err_pos = err0_p(id1,:,1);
        ydat = ydat1(id1,:,1);
        errorbar(xdat(:),ydat(:),err_neg(:),err_pos(:),'+k',...
            'LineWidth',0.5,'Clipping','off');
    end
    % group 2: with wave, surface cooling > 5 W m^-2
    for i = 1:nWD
        xdat = xdat0(id1,i,id3);
        ydat = ydat0(id1,i,id3);
        imarker = marker_mat(i);
        plot(xdat(:),ydat(:),imarker,'LineWidth',1.0,'Color',c_gray,'Clipping','off');
        if l_errbar
            err_neg = err0_n(id1,i,id3);
            err_pos = err0_p(id1,i,id3);
            ydat = ydat1(id1,i,id3);
            imarker = '+';
            errorbar(xdat(:),ydat(:),err_neg(:),err_pos(:),...
                imarker,'LineWidth',0.5,'Color',c_gray,'Clipping','off');
        end
    end
    % group 3: no wave, surface cooling = 5 W m^-2
    for i = 1:nWD
        xdat = xdat0(1,i,1);
        ydat = ydat0(1,i,1);
        imarker = ['x' color_mat(i)];
        plot(xdat(:),ydat(:),imarker,'LineWidth',1.0);
        if l_errbar
            err_neg = err0_n(1,i,1);
            err_pos = err0_p(1,i,1);
            ydat = ydat1(1,i,1);
            imarker = ['+' color_mat(i)];
            errorbar(xdat(:),ydat(:),err_neg(:),err_pos(:),...
                imarker,'LineWidth',0.5);
        end
    end
    % group 4: wind waves, surface cooling = 5 W m^-2
    for i = 1:nWD
        xdat = xdat0(1,i,id3w);
        ydat = ydat0(1,i,id3w);
        imarker = [marker_mat(i) color_mat(i)];
        plot(xdat(:),ydat(:),imarker,'LineWidth',1.0);
        if l_errbar
            err_neg = err0_n(1,i,id3w);
            err_pos = err0_p(1,i,id3w);
            ydat = ydat1(1,i,id3w);
            imarker = ['+' color_mat(i)];
            errorbar(xdat(:),ydat(:),err_neg(:),err_pos(:),...
                imarker,'LineWidth',0.5);
        end
    end
    % group 5: swells, surface cooling = 5 W m^-2
    for i = 1:nWD
        xdat = xdat0(1,i,id3s);
        ydat = ydat0(1,i,id3s);
        imarker = [marker_mat(i) color_mat(i)];
        scatter(xdat(:),ydat(:),imarker,'filled','LineWidth',1.0);
        if l_errbar
            err_neg = err0_n(1,i,id3s);
            err_pos = err0_p(1,i,id3s);
            ydat = ydat1(1,i,id3s);
            imarker = ['+' color_mat(i)];
            errorbar(xdat(:),ydat(:),err_neg(:),err_pos(:),...
                imarker,'LineWidth',0.5);
        end
    end
    
    xlabel(xlabel_str,'Interpreter','Latex');
    ylabel(ylabel_str,'Interpreter','Latex');

    if xlog
        set(gca,'xscale','log');
    end
    if ylog
        set(gca,'yscale','log');
    end

    xlim(xlims);
    ylim(ylims);

    h = gcf;
end

function h = plot_dist_4p_l(hst,xi,yi)
% plot_dist_4p_l plot the highest 30%, 60%, 90% and 99%
%   centered distribution

    % find the isolines of pdf that enclose the area in which the
    % total probability is 30%, 60%, 90% and 99%
    hsum = sum(hst(:));
    hlist = sort(hst(:),'descend')./hsum;
    hcum = cumsum(hlist);
    vl = [0.3, 0.6, 0.9, 0.99];
    nv = numel(vl);
    vlev = zeros(1,nv);
    for i=1:nv
        [~,ind] = min(abs(hcum-vl(i)));
        vlev(i) = hlist(ind);
    end
    pdfData = hst./hsum;
    % plot log10(pdfData) to make use of the full colarbar
    pdfData(pdfData==0) = 1e-12;    % avoid -inf for log10(0)
    [~,h] = contour(xi,yi,pdfData');
    h.LevelListMode = 'manual';
    h.LevelList = vlev;
    h.LineColor = [0.9 0.9 0.9];
    h.LineWidth = 1.0;
    h.ShowText = 'on';
    h.TextList = vl;
end

function h = newFigure(l_save_fig)
% newFigure() creates new figure object. 
%   Visible if l_save_fig = 0, invisible otherwise.

    if (l_save_fig)
        figure('visible','off');
    else
        figure;
    end
    h = gcf;
end

function postProcessFig(gcf, filename, varargin)
% postProcessFig resizes a figure, convert it to eps format
%   postProcessFig(figname,[width, hight])

% counts arguments
% counts arguments
nArgs = length(varargin);
if nArgs==0
    wsize = 6;
    hsize = 4;
elseif nArgs==2
    if isnumeric(varargin{1}) && isnumeric(varargin{2})
        wsize = varargin{1};
        hsize = varargin{2};
    else
        error('Numeric arguments required...');
    end
else
    error('No arguments or pass in width and hight');
end

mleft = 0.25;
mdown = 0.25;
mright = wsize-0.25;
mup = hsize-0.25;

% Get x, y limits
xx = xlim;
yy = ylim;
% Change size
set(gcf,'Units','inches');
set(gcf,'PaperUnits','inches'); 
set(gcf,'PaperPositionMode', 'manual');
paperposition = [mleft mdown mright mup];
set(gcf,'PaperPosition',paperposition);
size = [wsize hsize];
set(gcf,'PaperSize',size);
position = [mleft mdown mright mup];
set(gcf,'Position',position);
% Change renderer, 'painters' works well for vector formats while
% 'opengl' works well for bitmap format
set(gcf, 'renderer', 'painters');
% Change font
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
        
% Export file
% check if contains inlet
h = findall(gcf,'Type','Axes');
if length(h)==2
    h(1).FontSize = 9;
end
xlim(xx);
ylim(yy);
exportgraphics(gcf, filename)
end