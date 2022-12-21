 close all; clear variables;
% This script plots scatters of a group of simulations

% Qing Li, 170130

% add path if necessary
if exist('qingAddPath.m','file') == 2
    qingAddPath
end
% add path for ncarleePflData
path(path,'../share');
% save figures or not (0 for testing)
l_save_fig = 0;
l_h = 2; % 1:h (buoyancy flux minimum),
         % 2:hm (N2 maximum),
         % 3:hb (dissipation threshold, 1e-9)
         
%% Set cases
% get profile group
% groupName = 'Langmuir';
% groupName = 'LangmuirIO';
% groupName = 'Res_Scale2';
% groupName = 'Convection';
groupName = 'LangmuirML'

%% run grp_cases.m
% This file should contain
%   caselist: a list of case name
%   lcolor: a list of line colors/types for profile plot
%   llabel: label for legend
%   lmarker: a list of markers for scatter plot
%   lmcolor: a list of marker colors
%   l_filled: a list of 1/0 for filled marker
scptName = [groupName '/grp_cases.m'];
if exist(scptName,'file')
    run(scptName);
else
    error('Group list %s does not exist. Stop.', scptName);
end

%% plot list
% group A
l_parSpace = 0;         % done, regime diagram (Belcher et al., 2012, BG12)
l_check_diss = 0;       % done, check dissipation versus BG12
% group B
l_wps_bflux = 0;        % done, mean w' squared versus w* squared
l_wps_la = 0;           % done, mean w' squared versus La_t^{-2}
l_wps_lasl = 1;         % done, mean w' squared versus La_{SL}^{-2}
l_wps_tke_bflux = 0;    % done, VKE/TKE, versus w* squared
l_wps_tke_lasl = 0;     % done, VKE/TKE, versus La_t^{-2}
l_wps_tke_la = 0;       % done, VKE/TKE, versus La_{SL}^{-2}
% group C
l_entr_bflux = 0;       % done, -min(w'b')/B0 versus u*^3/w*^3
l_entr_la = 0;          % done, -min(w'b') versus La_t^{-2}, with offset
l_entr_lasl = 1;        % done, -min(w'b') versus La_{SL}^{-2}, with offset
l_entr_check = 0;       % done, check the scaling for -min(w'b')
l_entr_parSpace = 0;    % done, regime diagram for entrainment.
% group D
l_wps_entr = 0;
% group E
l_tke_la = 0;           % done, mean TKE, versus La_t^{-2}
l_tke_lasl = 0;         % done, mean TKE, versus La_{SL}^{-2}
l_tke_bflux = 0;        % done, mean TKE, versus w* squared
% group F
l_diss_bflux = 0;       % to do
l_diss_la = 0;          % to do
l_diss_lasl = 0;        % to do
% group X
l_hbhm = 0;             % done, boundary layer depth versus mixed layer
l_hbhm_wstar = 0;       % done, hb/h, versus w*
l_hbkpp = 0;            % done, boundary layer depth, LES versus KPP
l_hbkpp_wstar = 0;      % done, boundary layer depth, KPP/LES versus w*

%% Get the data

% cases
nsize = size(caselist);
nc = nsize(1);

% gray
c_gray = [0.5 0.5 0.5];
% data dir
dir = '../AnalyzeSingleCase';
% variables
stat_wb_mean = zeros(nc,1); % -min(w'b'), mean
stat_wb_median = zeros(nc,1); % median, over an inertial period
stat_wb_p25 = zeros(nc,1); % 25th percentile
stat_wb_p75 = zeros(nc,1); % 75th percentile
utau = zeros(nc,1); % friction velocity
wpsm = zeros(nc,2); % h-averaged vertical tke * 2
tkem = zeros(nc,2); % h-averaged total tke, including SGS
h = zeros(nc,1); % boundary layer depth based on minimum wb
hm = zeros(nc,1); % mixed layer depth based on maximum N2
hb = zeros(nc,1); % boundary layer depth based on dissipation rate
hbkpp = zeros(nc,1); % boundary layer depth diagnosed in KPP
hbkpp0 = zeros(nc,1); % hb KPP
hbkpp1 = zeros(nc,1); % hb KPP, LW16 MA
hbkpp2 = zeros(nc,1); % hb KPP, RW16
hbkpp3 = zeros(nc,1); % hb KPP, LW16 EN
hbkpp4 = zeros(nc,1); % hb KPP, LF17
epsilonb = 1e-9;    % dissipation rate threshold to determine 
                    % the boundary layer
us = zeros(nc,1); % surface Stokes drift, at first level (-0.32 m)
b0 = zeros(nc,1); % surface buoyancy flux
eps = zeros(nc,2); % dissipation at z=h/2
Ds = zeros(nc,1); % stokes depth over h
lasl = zeros(nc,1); % surface layer averaged Langmuir number
laslp = zeros(nc,1); % surface layer averaged, projected Langmuir number
d0tketrans = zeros(nc,1); % depth at which the tke transport change sign
alphaL = zeros(nc,1); % angle between LC and wind, Lagrangian shear
alphaLOW = zeros(nc,1); % angle between LC and wind, law of the wall
for j=1:nc
    casename = caselist(j,:);
    fileName = [dir '/' casename '/MeanProfile.mat'];
    if exist(fileName,'file')
        fprintf('Case Name: %s; Loading...\n', casename);
        dat = load(fileName);
        utau(j,1) = dat.utau;
        h(j,1) = dat.h;
        hm(j,1) = dat.hm;
        zw = dat.z_w;
        zu = dat.z_u;
        stat_wb_mean(j,1) = dat.stat_wb.mean;
        stat_wb_median(j,1) = dat.stat_wb.median;
        stat_wb_p25(j,1) = dat.stat_wb.p25;
        stat_wb_p75(j,1) = dat.stat_wb.p75;
        wpsm(j,1) = dat.f.getHAvgPfl(dat.wps,zw,h(j,1));
        wpsm(j,2) = dat.f.getHAvgPfl(dat.wps,zw,hm(j,1));
        tkem(j,1) = dat.f.getHAvgPfl(dat.tke(1,:),zw,h(j,1));
        tkem(j,2) = dat.f.getHAvgPfl(dat.tke(1,:),zw,hm(j,1));
        us(j,1) = dat.stokes(1);
        b0(j,1) = dat.wtsfc.*dat.f.alpha.*dat.f.g;
        tmp = dat.dissip;
        eps(j,1) = interp1(zw,tmp,-dat.h/2,'linear');
        eps(j,2) = interp1(zw,tmp,-dat.hm/2,'linear');
        if isfield(dat,'hb')
            hb(j,1) = dat.hb;
        else
            hb(j,1) = -dat.f.getThresholdDepth(tmp,zw,epsilonb);
        end
        hbkpp(j,1) = dat.hb_kpp;
        % get KPP boundary layer depth
        stokes = dat.stokes;
        bbb = dat.bxym;
        uuu = dat.uxym;
        uul = dat.uxym+stokes;
        vvv = dat.vxym;
        w_s = zeros(size(zu));
        la = zeros(size(zu));
        us0 = dat.stokes(1);
        ustar = dat.utau;
        bf = -dat.wtsfc.*dat.f.alpha.*dat.f.g;
        surf_layer_ext = dat.f.kpp_get_constant('surf_layer_ext');
        for l=1:dat.f.nnz
            w_s(l) = dat.f.kpp_turbulent_scale(surf_layer_ext,...
                            -zu(l),bf,ustar,'s');
            stokessl = dat.f.getUsSL(stokes,zu,-zu(l));
            la(l) = sqrt(ustar./(stokessl-stokes(l)));
        end
        [hbkpp0(j,1), ~] = dat.f.kpp_boundary_layer_depth(bbb',...
                             uuu',vvv',zu',zw',w_s');
        [hbkpp1(j,1), ~] = dat.f.kpp_boundary_layer_depth(bbb',...
                             uuu',vvv',zu',zw',w_s','LW16',la');
        [hbkpp2(j,1), ~] = dat.f.kpp_boundary_layer_depth(bbb',...
                             uul',vvv',zu',zw',w_s','RW16',la');
        [hbkpp3(j,1), ~] = dat.f.kpp_boundary_layer_depth(bbb',...
                             uuu',vvv',zu',zw',w_s','LW16_En',la',us0);
        [hbkpp4(j,1), ~] = dat.f.kpp_boundary_layer_depth(bbb',...
                             uuu',vvv',zu',zw',w_s','LF17',la',-bf,ustar);
        % end KPP      
        sdepth = dat.f.getStokesDepth(dat.stokes,zu);
        Ds(j,1) = sdepth./dat.hm;
        stokessl = dat.f.getUsSL(dat.stokes,zu,dat.hm);
        [~,indhm] = min(abs(zw(:)+dat.hm));
        lasl(j,1) = sqrt(dat.utau./(stokessl-dat.stokes(indhm)));
        tmp = dat.tkeTrans;
        d0tketrans(j,1) = -dat.f.getPosDepth(tmp,zw);
        alphaL(j,1) = dat.f.getAlphaL(dat.uxym, dat.vxym,...
            dat.stokes.*cosd(thww{j}), dat.stokes.*sind(thww{j}),...
            dat.z_u, dat.h);
        alphaLOW(j,1) = dat.f.getAlphaLOW(ustar,us0,dat.h,thww{j},0.8);
        laslp(j,1) = sqrt(dat.utau.*cosd(alphaLOW(j,1))...
            ./((stokessl-dat.stokes(indhm)).*cosd(thww{j}-alphaLOW(j,1))));

    else
        fprintf('Case Name: %s; Not exist, skipping...\n', casename);
        utau(j,1) = NaN;
        h(j,1) = NaN;
        hm(j,1) = NaN;
        hb(j,1) = NaN;
        hbkpp(j,1) = NaN;
        hbkpp0(j,1) = NaN;
        hbkpp1(j,1) = NaN;
        hbkpp2(j,1) = NaN;
        hbkpp3(j,1) = NaN;
        hbkpp4(j,1) = NaN;
        stat_wb_mean(j,1) = NaN;
        stat_wb_median(j,1) = NaN;
        stat_wb_p25(j,1) = NaN;
        stat_wb_p75(j,1) = NaN;
        wpsm(j,:) = NaN;
        tkem(j,:) = NaN;
        us(j,1) = NaN;
        b0(j,1) = NaN;
        eps(j,:) = NaN;
        Ds(j,1) = NaN;
        lasl(j,1) = NaN;
        d0tketrans(j,1) = NaN;
        alphaL(j,1) = NaN;
        alphaLOW(j,1) = NaN;
        laslp(j,1) = NaN;
    end
end

% Choice of boundary layer depth for scaling
if l_h==1
    hx = h;
elseif l_h==2
    hx = hm;
elseif l_h==3
    hx = hb;
else
    error('Unsupported choice of l_h');
end

% Langmuir number
la = sqrt(utau./us);
% la(isinf(la))=1;

% h/LL
hL = b0.*hx./utau.^2./us;

fn_prefix = 'sct'; % figure name prefix
fn_suffix = ''; % figure name suffix


%% figure A1: parameter space
if l_parSpace
newFigure(l_save_fig);

% Reproduce Fig. 3 in Belcher et al., 2012
xlims = [0.1,10];   % La
ylims = [1e-3,1e3]; % h/LL

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

xdat0 = la;
ydat0 = hL;
xlog = 1;
ylog = 1;
xlims = [0.1 10];
ylims = [1e-3,1e3];
xlabel_str = '$La_\mathrm{t}$';
ylabel_str = '$h/L_\mathrm{L}$';

lmcolorw = lmcolor;
lmcolorw(:) = {'w'};

plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolorw,lmarker,l_filled);

tx = text(0.11, 2e-3,'Langmuir','Interpreter','Latex');
tx.BackgroundColor = 'w';
tx = text(3, 2e-2,'Wind','Interpreter','Latex');
tx.BackgroundColor = 'w';
tx = text(0.13, 1e2,'Convection','Interpreter','Latex');
tx.BackgroundColor = 'w';

if (l_save_fig)
    figname = 'sct_ParSpace.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,5,5);
end
end

%% figure A2: check dissipation scaling vs Belcher et al., 2012
if l_check_diss
newFigure(l_save_fig);

wgt = utau.^3./h;
xdat0 = 2.*(1-exp(-0.5.*la))+0.22.*la.^-2+0.3.*b0./wgt;
ydat0 = eps(:,1)./wgt;

xlabel_str = ['$A_\mathrm{s} +$ '...
              '$A_\mathrm{L} w_\mathrm{*L}^3/u_*^3+$ '...
              '$A_\mathrm{c} w_*^3/u_*^3$'];
ylabel_str = '$\epsilon|_{h/2} h/u_*^3$';
xlims = [1,50];
ylims = [1,50];
xlog = 1;
ylog = 1;
rl = refline(1,0);
rl.XData = xlims;
rl.YData = ylims;
rl.Color = c_gray;
daspect([1 1 1]);
% plot scatter
plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled);

if (l_save_fig)
    figname = 'sct_eps_check.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,5,5);
end
end


%% figure B1: wps vs surface buoyancy flux (w*^2)
if l_wps_bflux
newFigure(l_save_fig);

% set data
xdat0 = (b0.*hx).^(2/3)./utau.^2;
ydat0 = wpsm(:,1)./utau.^2;
xlabel_str = '$w_*^2/u_*^2$';
ylabel_str = '$\langle \overline{w''^2}\rangle_{h_\mathrm{b}}/u_*^2$';
xlog = 1;
ylog = 0;
xlims = [0.1, 30];
ylims = [0, 7];

% --- ref line
hold on;
nx = 10000;
dx = (xlims(2)-xlims(1))/nx;
x = xlims(1):dx:xlims(2)*2;
y = 0.3.*x;
y2 = x.*0+0.6;
plot(x,y,'-','LineWidth',1.0,'Color',c_gray);
plot(x,y2,'-','LineWidth',1.0,'Color',c_gray);
text(8, 6,'$y=0.3x$','Interpreter','Latex');
text(10, 0.9,'$y=0.6$','Interpreter','Latex');

% --- ref line

% plot scatter 
plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled);

% save figure
if (l_save_fig)
    figname = 'sct_wps_wstar.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end

%% figure B2: wps vs La (us0^2)
if l_wps_la
newFigure(l_save_fig);

xdat0 = la.^-2;
ydat0 = wpsm(:,1)./utau.^2;
xlabel_str = '$La_\mathrm{t}^{-2}=u^\mathrm{S}_0/u_*$';
ylabel_str = '$\langle \overline{w''^2}\rangle_{h_\mathrm{b}}/u_*^2$';
xlog = 0;
ylog = 0;
xlims = [0 20];
ylims = [0 7];

% --- ref line
hold on;
nx = 10000;
dx = (xlims(2)-xlims(1))/nx;
x = xlims(1):dx:xlims(2)*2;
y = 0.6.*(1+(3.1).^(-2).*x+(5.7).^(-4).*x.^2);
y2 = 0.64.*(1+0.098.*x);
p1 = plot(x,y,'-','LineWidth',0.8,'Color',c_gray);
p2 = plot(x,y2,'--','LineWidth',0.8,'Color',c_gray);
eq_str = '$0.6(1+(3.1 La_\mathrm{t})^{-2}+(5.7 La_\mathrm{t})^{-4})$';
eq2_str = '$0.64(1+0.098 La_\mathrm{t}^{-2})$';
lg = legend([p1, p2],eq_str,eq2_str);
lg.Location = 'NorthEast';
lg.Interpreter = 'Latex';
lg.Color = 'none';
% --- ref line

% plot scatter
plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled);

if (l_save_fig)
    figname = 'sct_wps_La.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end

%% figure B3: wps vs LaSL^{-2} (ussl^2)
if l_wps_lasl
newFigure(l_save_fig);

xdat0 = lasl.^(-2);
ydat0 = wpsm(:,1)./utau.^2;
xlabel_str = '$La_\mathrm{SL}^{-2}=u^\mathrm{S}_\mathrm{SL}/u_*$';
ylabel_str = '$\langle \overline{w''^2}\rangle_{h_\mathrm{b}}/u_*^2$';
xlog = 0;
ylog = 0;
xlims = [0 6];
ylims = [0 7];

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
eq_str = '$0.6 (1 + (1.5 La_\mathrm{SL})^{-2}+(5.4 La_\mathrm{SL})^{-4})$';
eq2_str = '$0.398 + 0.480 La_\mathrm{SL}^{-4/3}$';
eq3_str = '$0.640+3.50\exp(-2.69 La_\mathrm{SL})$';
lg = legend([p1, p2, p3],eq_str,eq2_str,eq3_str);
lg.Location = 'NorthEast';
lg.Interpreter = 'Latex';
% --- ref line

% plot scatter
plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled);

if (l_save_fig)
    figname = 'sct_wps_LaSL.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end

%% figure B4: wps vs tke, wstar
if l_wps_tke_lasl
newFigure(l_save_fig);

% set data
xdat0 = (b0.*hx).^(2/3)./utau.^2;
ydat0 = 0.5.*wpsm(:,2)./tkem(:,2);

xlabel_str = '$w_*^2/u_*^2$';
ylabel_str = 'VKE/TKE';
xlims = [0.1,30];
ylims = [0,0.5];
xlog = 1;
ylog = 0;

% plot scatter
plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled);

% save figure
if (l_save_fig)
    figname = 'sct_tke_wps_wstar.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end

%% figure B5: wps vs tke, lasl
if l_wps_tke_lasl
newFigure(l_save_fig);

% set data
xdat0 = lasl.^(-2);
ydat0 = 0.5.*wpsm(:,2)./tkem(:,2);

xlabel_str = '$La_\mathrm{SL}^{-2}=u^\mathrm{S}_\mathrm{SL}/u_*$';
ylabel_str = 'VKE/TKE';
xlims = [0,6];
ylims = [0,0.5];
xlog = 0;
ylog = 0;

% plot scatter
plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled);

% save figure
if (l_save_fig)
    figname = 'sct_tke_wps_LaSL.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end

%% figure B6: wps vs tke, la
if l_wps_tke_la
newFigure(l_save_fig);

% set data
xdat0 = la.^(-2);
ydat0 = 0.5.*wpsm(:,2)./tkem(:,2);

xlabel_str = '$La_\mathrm{t}^{-2}=u^\mathrm{S}_0/u_*$';
ylabel_str = 'VKE/TKE';
xlims = [0,20];
ylims = [0,0.5];
xlog = 0;
ylog = 0;

% plot scatter
plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled);

% save figure
if (l_save_fig)
    figname = 'sct_tke_wps_La.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end

%% figure C1: Entrainment vs surface buoyancy flux
if l_entr_bflux
newFigure(l_save_fig);

% set data
xdat0 = utau.^3./b0./hx;
wgt = b0;
ydat0 = -stat_wb_mean./wgt;
ydat1 = -stat_wb_median./wgt;
err0_n = -(stat_wb_p75-stat_wb_median)./wgt;
err0_p = -(stat_wb_median-stat_wb_p25)./wgt;
xlabel_str = '$u_*^3/w_*^3$';
ylabel_str = '$\overline{w''b''}_\mathrm{e} h_\mathrm{b}/w_*^3$';
xlog = 1;
ylog = 0;
ylims = [0,9];
xlims = [0.01,100];

% --- ref line
hold on;
par1 = 0.158;
par2 = 0.112;
nx = 10000;
dx = (xlims(2)-xlims(1))/nx;
x = 0:dx:xlims(2)*2;
y = par1.*x+par2;
y2 = x.*0+0.15;
plot(x,y,'-','LineWidth',1.0,'Color',c_gray);
plot(x,y2,'-','LineWidth',1.0,'Color',c_gray);
text(12,0.8,'$y=0.15$','Interpreter','Latex');
eq_str = sprintf('$y=%4.2f x + %4.2f$',par1,par2);
text(0.3, 3, eq_str,'Interpreter','Latex');
% --- ref line

% plot scatter
plotDataScatter(xdat0,ydat0,xlabel_str,ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled,...
                ydat1,err0_n,err0_p);

% add inlet
% create a new pair of axes inside current figure
axes('position',[.2 .6 .45 .3]);
xlims1 = [0.01,1];
ylims1 = [0,0.5];
% box on; % put box around new pair of axes
plotDataScatter(xdat0,ydat0,[],[],...
                xlog,ylog,xlims1,ylims1,...
                lmcolor,lmarker,l_filled,...
                ydat1,err0_n,err0_p);
plot(x,y,'-','LineWidth',1.0,'Color',c_gray);
plot(x,y2,'-','LineWidth',1.0,'Color',c_gray);

% save figure
if (l_save_fig)
    figname = 'sct_wb_wstar.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end

%% figure C2: entrainemnt vs la
if l_entr_la
newFigure(l_save_fig);

% set data
xdat0 = la.^-2;
wgt = utau.^3./hx;
soffset = 0.15;
% soffset = ff1.p2;
yoffset = soffset.*b0./wgt;
ydat0 = -(stat_wb_mean)./wgt-yoffset;
ydat1 = -(stat_wb_median)./wgt-yoffset;
err0_n = -(stat_wb_p75-stat_wb_median)./wgt;
err0_p = -(stat_wb_median-stat_wb_p25)./wgt;
xlabel_str = '$La_\mathrm{t}^{-2}=u^\mathrm{S}_0/u_*$';
ylabel_str = ['$(\overline{w''b''}_\mathrm{e}-'...
              sprintf('%4.2g',soffset)...
              ' B_0) h_\mathrm{b}/u_*^3$'];

xlog = 0;
ylog = 0;
xlims = [0,20];
ylims = [-2,5];

% --- ref line
hold on;
par1 = 0.0265;
par2 = 0.0130;
nx = 1000;
dx = (xlims(2)-xlims(1))/nx;
x = xlims(1):dx:xlims(2)*2;
y = par1.*x+par2;
plot(x,y,'-k','LineWidth',1.0,'Color',c_gray);
eq_str = sprintf('$y = %4.2g x + %4.2g$',par1,par2);
text(10, -1, eq_str,'Interpreter','Latex');
% --- ref line

% plot scatter
plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled,...
                ydat1,err0_n,err0_p);

% add inlet
% create a new pair of axes inside current figure
axes('position',[.52 .5 .35 .35]);
xlims1 = [0,16];
ylims1 = [0,0.6];
% box on; % put box around new pair of axes
plotDataScatter(xdat0,ydat0,[],[],...
                xlog,ylog,xlims1,ylims1,...
                lmcolor,lmarker,l_filled,...
                ydat1,err0_n,err0_p);
plot(x,y,'-','LineWidth',1.0,'Color',c_gray);
            
% save figure
if (l_save_fig)
    figname = 'sct_wb_La.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end

%% figure C3: entrainemnt vs laSL
if l_entr_lasl
newFigure(l_save_fig);

% set data
xdat0 = lasl.^(-2);
wgt = utau.^3./hx;
soffset = 0.15;
% soffset = ff1.p2;
yoffset = soffset.*b0./wgt;
ydat0 = -(stat_wb_mean)./wgt-yoffset;
ydat1 = -(stat_wb_median)./wgt-yoffset;
err0_n = -(stat_wb_p75-stat_wb_median)./wgt;
err0_p = -(stat_wb_median-stat_wb_p25)./wgt;
xlabel_str = '$La_\mathrm{SL}^{-2}=u^\mathrm{S}_\mathrm{SL}/u_*$';
ylabel_str = ['$(\overline{w''b''}_\mathrm{e}-'...
              sprintf('%4.2g',soffset)...
              ' B_0) h_\mathrm{b}/u_*^3$'];
xlog = 0;
ylog = 0;
xlims = [0,5];
ylims = [-1.5,1.5];

% --- ref line
hold on;
par1 = 0.083;
par2 = 0.17;
nx = 1000;
dx = (xlims(2)-xlims(1))/nx;
x = xlims(1):dx:xlims(2)*2;
y = par1.*x+par2;
plot(x,y,'-','LineWidth',1.0,'Color',c_gray);
eq_str = sprintf('$y=%4.2g x+%4.2g$',par1,par2);
text(3, -0.5, eq_str,'Interpreter','Latex');
% --- ref line

% plot scatter
plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled,...
                ydat1,err0_n,err0_p);

% add inlet
% create a new pair of axes inside current figure
% axes('position',[.42 .52 .45 .35]);
% xlims1 = [0,6];
% ylims1 = [0,0.6];
% % box on; % put box around new pair of axes
% plotDataScatter(xdat0,ydat0,[],[],...
%                 xlog,ylog,xlims1,ylims1,...
%                 lmcolor,lmarker,l_filled,...
%                 ydat1,err0_n,err0_p);
% plot(x,y,'-','LineWidth',1.0,'Color',c_gray);

% save figure
if (l_save_fig)
    figname = 'sct_wb_LaSL.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end

%% figure C4: check entrainemnt scaling
if l_entr_check
newFigure(l_save_fig);

% set data
c1 = 0.16;
c2 = 0.15;
c3 = 0.087;
% xdat0 = c1.*utau.^3./hx+c2.*b0+c3.*lasl.^(-2).*utau.^3./hx;
xdat0 = c1+c2.*b0.*hx./utau.^3+c3.*lasl.^(-2);
wgt = utau.^3./hx;
ydat0 = -(stat_wb_mean)./wgt;
ydat1 = -(stat_wb_median)./wgt;
err0_n = -(stat_wb_p75-stat_wb_median)./wgt;
err0_p = -(stat_wb_median-stat_wb_p25)./wgt;
% c1_str = sprintf('%4.2g',c1);
% c2_str = sprintf('%4.2g',c2);
% c3_str = sprintf('%4.2g',c3);
% xlabel_str = [c1_str '$u_*^3/h_\mathrm{b}+$ '...
%               c2_str '$w_*^3/h_\mathrm{b}+$ '...
%               c3_str '$u^\mathrm{S}_\mathrm{SL}u_*^2/h_\mathrm{b}$'];
xlabel_str = ['$c_\mathrm{b1}-$ '...
              '$c_\mathrm{b2} h_\mathrm{b}/L+$ '...
              '$c_\mathrm{b3} La_\mathrm{SL}^{-2}$'];
ylabel_str = '$\overline{w''b''}_\mathrm{e} h_\mathrm{b}/u_*^3$';
xlims = [1e-1,2];
ylims = [1e-1,2];
xlog = 1;
ylog = 1;
rl = refline(1,0);
rl.XData = xlims;
rl.YData = ylims;
rl.Color = c_gray;

% plot scatter
plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled,...
                ydat1,err0_n,err0_p);

daspect([1 1 1]);
% save figure
if (l_save_fig)
    figname = 'sct_wb_check.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,5,5);
end
end

%% figure C5: regime diagram for entrainment
if l_entr_parSpace
newFigure(l_save_fig);

xlims = [0.1,20];   % Lasl
ylims = [0.01,200]; % -h/L

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

xdat0 = la;
ydat0 = hL;
xlog = 1;
ylog = 1;
xlabel_str = '$La_\mathrm{SL}$';
ylabel_str = '$-h_\mathrm{b}/L$';

lmcolorw = lmcolor;
lmcolorw(:) = {'w'};

plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolorw,lmarker,l_filled);

tx = text(0.12, 0.06,'Langmuir','Interpreter','Latex');
tx.BackgroundColor = 'w';
tx = text(4, 0.06,'Wind','Interpreter','Latex');
tx.BackgroundColor = 'w';
tx = text(3, 20,'Convection','Interpreter','Latex');
tx.BackgroundColor = 'w';

if (l_save_fig)
    figname = 'sct_wb_ParSpace.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,5,5);
end
end

%% figure D1: wps vs entrainment
if l_wps_entr
newFigure(l_save_fig);

% set data
xdat0 = wpsm(:,2)./utau.^2;
wgt = utau.^3./hx;
ydat0 = -(stat_wb_mean)./wgt;
ydat1 = -(stat_wb_median)./wgt;
err0_n = -(stat_wb_p75-stat_wb_median)./wgt;
err0_p = -(stat_wb_median-stat_wb_p25)./wgt;
xlabel_str = '$\langle \overline{w''^2}\rangle_{h_\mathrm{b}}/u_*^2$';
ylabel_str = '$-\min(\overline{w''b''}) h_\mathrm{b}/u_*^3$';
xlims = [0,10];
ylims = [0,20];
xlog = 0;
ylog = 0;

% plot scatter
plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled,...
                ydat1,err0_n,err0_p);

% save figure
if (l_save_fig)
    figname = 'sct_wb_wps.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end

%% figure E1: tke vs La
if l_tke_la
newFigure(l_save_fig);

xdat0 = la.^(-2);
ydat0 = tkem(:,2)./utau.^2;
xlabel_str = '$La_\mathrm{t}^{-2}=u^\mathrm{S}_0/u_*$';
ylabel_str = '$0.5\langle \overline{u_i'' u_i''}\rangle_\mathrm{h}/u_*^2$';
xlog = 0;
ylog = 0;
xlims = [0, 20];
ylims = [0, 10];

% plot scatter
plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled);

if (l_save_fig)
    figname = 'sct_tke_La.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end

%% figure E2: tke vs LaSL
if l_tke_lasl
newFigure(l_save_fig);

xdat0 = lasl.^(-2);
ydat0 = tkem(:,2)./utau.^2;
xlabel_str = '$La_\mathrm{SL}^{-2}=u^\mathrm{S}_\mathrm{SL}/u_*$';
ylabel_str = '$0.5\langle \overline{u_i''u_i''}\rangle_\mathrm{h}/u_*^2$';
xlog = 0;
ylog = 0;
xlims = [0, 6];
ylims = [0, 10];

% plot scatter
plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled);

if (l_save_fig)
    figname = 'sct_tke_LaSL.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end

%% figure E3: tke vs bflux
if l_tke_bflux
newFigure(l_save_fig);

xdat0 = (b0.*hx).^(2/3)./utau.^2;
ydat0 = tkem(:,2)./utau.^2;
xlabel_str = '$w_*^2/u_*^2$';
ylabel_str = '$0.5\langle \overline{u_i''u_i''}\rangle_\mathrm{h}/u_*^2$';
xlog = 1;
ylog = 0;
xlims = [0.1, 30];
ylims = [0, 10];

% plot scatter
plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled);

if (l_save_fig)
    figname = 'sct_tke_wstar.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end

%% figure F1: Dissipation at h/2 vs surface buoyancy flux
if l_diss_bflux
newFigure(l_save_fig);

% set data
xdat0 = utau.^3./b0./hx;
wgt = b0;
ydat0 = eps(:,2)./wgt;
xlabel_str = '$u_*^3/w_*^3$';
ylabel_str = '$\epsilon|_{h_\mathrm{b}/2} h_\mathrm{b}/w_*^3$';
xlog = 1;
ylog = 0;
ylims = [0,50];
xlims = [0.01,100];

% --- ref line
% hold on;
% par1 = 0;
% par2 = 0;
nx = 10000;
dx = (xlims(2)-xlims(1))/nx;
x = 0:dx:xlims(2)*2;
% y = par1.*x+par2;
cc = 0.33;
cc_str = sprintf('%4.2g',cc);
y2 = x.*0+cc;
% plot(x,y,'-','LineWidth',1.0,'Color',c_gray);
plot(x,y2,'-','LineWidth',1.0,'Color',c_gray);
text(12,3,['$y=' cc_str '$'] ,'Interpreter','Latex');
% eq_str = sprintf('$y=%4.2f x + %4.2f$',par1,par2);
% text(0.3, 10, eq_str,'Interpreter','Latex');
% --- ref line

% plot scatter
plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled);

% add inlet
% create a new pair of axes inside current figure
axes('position',[.2 .6 .45 .3]);
xlims1 = [0.01,1];
ylims1 = [0,1];
% box on; % put box around new pair of axes
plotDataScatter(xdat0,ydat0,[],[],...
                xlog,ylog,xlims1,ylims1,...
                lmcolor,lmarker,l_filled);
% plot(x,y,'-','LineWidth',1.0,'Color',c_gray);
plot(x,y2,'-','LineWidth',1.0,'Color',c_gray);

% save figure
if (l_save_fig)
    figname = 'sct_eps_wstar.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end

%% figure F2: Dissipation at h/2 vs LaSL
if l_diss_bflux
newFigure(l_save_fig);

% set data
xdat0 = lasl.^(-2);
wgt = utau.^3./hx;
soffset = cc;
yoffset = soffset.*b0./wgt;
ydat0 = eps(:,2)./wgt-yoffset;
xlabel_str = '$La_\mathrm{SL}^{-2}$';
ylabel_str = ['$(\epsilon|_{h_\mathrm{b}/2}-'...
              sprintf('%4.2g',soffset) 'B_0) h_\mathrm{b}/u_*^3$'];
xlog = 0;
ylog = 0;
xlims = [0,6];
ylims = [0.5,3.5];

% --- ref line
% hold on;
% par1 = 0;
% par2 = 0;
% nx = 10000;
% dx = (xlims(2)-xlims(1))/nx;
% x = 0:dx:xlims(2)*2;
% y = par1.*x+par2;
% plot(x,y,'-','LineWidth',1.0,'Color',c_gray);
% eq_str = sprintf('$y=%4.2f x + %4.2f$',par1,par2);
% text(3, 1, eq_str,'Interpreter','Latex');
% --- ref line

% plot scatter
plotDataScatter(xdat0,ydat0,xlabel_str,ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled);

% save figure
if (l_save_fig)
    figname = 'sct_eps_LaSL.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end

%% figure F3: Dissipation at h/2 vs La
if l_diss_bflux
newFigure(l_save_fig);

% set data
xdat0 = la.^(-2);
wgt = utau.^3./hx;
soffset = cc;
yoffset = soffset.*b0./wgt;
ydat0 = eps(:,2)./wgt-yoffset;
xlabel_str = '$La_\mathrm{t}^{-2}$';
ylabel_str = ['$(\epsilon|_{h_\mathrm{b}/2}-'...
              sprintf('%4.2g',soffset) 'B_0) h_\mathrm{b}/u_*^3$'];
xlog = 0;
ylog = 0;
xlims = [0,20];
ylims = [0,5];

% --- ref line
% hold on;
% par1 = 0;
% par2 = 0;
% nx = 10000;
% dx = (xlims(2)-xlims(1))/nx;
% x = 0:dx:xlims(2)*2;
% y = par1.*x+par2;
% plot(x,y,'-','LineWidth',1.0,'Color',c_gray);
% eq_str = sprintf('$y=%4.2f x + %4.2f$',par1,par2);
% text(12, 1, eq_str,'Interpreter','Latex');
% --- ref line

% plot scatter
plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolor,lmarker,l_filled);

% save figure
if (l_save_fig)
    figname = 'sct_eps_La.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end

%% figure X1: boundary layer depth vs mixed layer depth
if l_hbhm
newFigure(l_save_fig);

xlabel_str = '$h$ (m)';
ylabel_str = '$h_\mathrm{b}$ (m)';
xlims = [40, 85];
ylims = [40, 85];
xlog = 0;
ylog = 0;

% plot scatter
hold on;
rl=refline(1,0);
rl.XData = xlims;
rl.YData = ylims;
rl.Color = c_gray;

xdat0 = h;
ydat0 = hm;
lmcolorw = lmcolor;
lmcolorw(:) = {'k'};
lmarkerw = lmarker;
lmarkerw(:) = {'o'};
p1 = plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolorw,lmarkerw,l_filled);
xdat0 = h;
ydat0 = hb;
lmcolorw(:) = {'r'};
lmarkerw(:) = {'+'};
p2 = plotDataScatter(xdat0,ydat0,xlabel_str, ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolorw,lmarkerw,l_filled);
            
lg = legend([p1,p2],'$z|\max(N^2)$','$z|\epsilon=10^{-9}$ m$^2$ s$^{-3}$');
lg.Location = 'NorthWest';
lg.Interpreter = 'Latex';

daspect([1 1 1]);
if (l_save_fig)
    figname = 'sct_HbHm.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,5,5);
end
end

%% figure X2: kpp boundary layer depth vs les boundary layer depth
if l_hbkpp
newFigure(l_save_fig);

xdat0 = hm;
ydat1 = -hbkpp0;
ydat2 = -hbkpp1;
ydat3 = -hbkpp2;
ydat4 = -hbkpp3;
ydat5 = -hbkpp4;
xlabel_str = '$h_\mathrm{b,LES}$ (m)';
ylabel_str = '$h_\mathrm{b,KPP}$ (m)';
xlims = [40, 85];
ylims = [40, 85];
xlog = 0;
ylog = 0;

hold on;
rl = refline(1,0);
rl.XData = xlims;
rl.YData = ylims;
rl.Color = c_gray;

lmcolorw = lmcolor;
lmarkerw = lmarker;
lmcolorw(:) = {'k'};
lmarkerw(:) = {'o'};
p1 = plotDataScatter(xdat0,ydat1,xlabel_str,ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolorw,lmarker,l_filled);
lmcolorw(:) = {'r'};
lmarkerw(:) = {'*'};
p2 = plotDataScatter(xdat0,ydat2,xlabel_str,ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolorw,lmarker,l_filled);
lmcolorw(:) = {'c'};
lmarkerw(:) = {'d'};
p3 = plotDataScatter(xdat0,ydat3,xlabel_str,ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolorw,lmarker,l_filled);
lmcolorw(:) = {'b'};
lmarkerw(:) = {'+'};
p4 = plotDataScatter(xdat0,ydat4,xlabel_str,ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolorw,lmarker,l_filled);
lmcolorw(:) = {'m'};
lmarkerw(:) = {'x'};
p5 = plotDataScatter(xdat0,ydat4,xlabel_str,ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolorw,lmarker,l_filled);

lg = legend([p1,p2,p4,p3,p5],...
            'LMD94','LW16-MA','LW16-EN','RW16','This study');
lg.Location = 'SouthEast';
lg.Interpreter = 'Latex';
daspect([1 1 1]);

if (l_save_fig)
    figname = 'sct_HbKPP.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,5,5);
end
end

%% figure X3: boundary layer depth / mixed layer depth vs -h/L
if l_hbhm_wstar
newFigure(l_save_fig);

ydat1 = hm./h;
ydat2 = hb./h;
ydat3 = -hbkpp0./h;
xdat0 = b0.*hm./utau.^3.;
xlims = [0.04, 200];
ylims = [0.9, 1.3];
xlog = 1;
ylog = 0;
xlabel_str = '$-h_\mathrm{b}/L$';
ylabel_str = '$h_\mathrm{b}/h$';

ynan = nan(size(xdat0));

hold on;
p1 = scatter(xdat0(:),ynan(:),'ko','filled'); 
p2 = scatter(xdat0(:),ynan(:),'bo','filled'); 
p3 = scatter(xdat0(:),ynan(:),'ro','filled'); 
rl = refline(0,1);
rl.XData = xlims;
rl.YData = [1, 1];
rl.Color = c_gray;

lmcolorw = lmcolor;
lmcolorw(:) = {'k'};
plotDataScatter(xdat0,ydat1,xlabel_str,ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolorw,lmarker,l_filled);
lmcolorw(:) = {'b'};
plotDataScatter(xdat0,ydat2,xlabel_str,ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolorw,lmarker,l_filled);
lmcolorw(:) = {'r'};
plotDataScatter(xdat0,ydat3,xlabel_str,ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolorw,lmarker,l_filled);

lg = legend([p1,p2,p3],...
    '$z|\max(N^2)$',...
    '$z|\epsilon=10^{-9}$ m$^2$ s$^{-3}$',...
    'KPP');
lg.Location = 'NorthEast';
lg.Interpreter = 'Latex';

if (l_save_fig)
    figname = 'sct_HbHm_wstar.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end


%% figure X4: kpp boundary layer depth / les boundary layer depth vs -h/L
if l_hbkpp_wstar
newFigure(l_save_fig);

ydat1 = -hbkpp0./hm;
ydat2 = -hbkpp2./hm;
ydat3 = -hbkpp3./hm;
ydat4 = -hbkpp4./hm;
xdat0 = b0.*hm./utau.^3.;
xlims = [0.04, 200];
ylims = [0.8, 1.6];
xlabel_str = '$-h_\mathrm{b}/L$';
ylabel_str = '$h_\mathrm{b,KPP}/h_\mathrm{b,LES}$';
xlog = 1;
ylog = 0;

hold on;
ynan = nan(size(xdat0));

p1 = scatter(xdat0(:),ynan(:),'ko','filled'); 
p3 = scatter(xdat0(:),ynan(:),'co','filled'); 
p4 = scatter(xdat0(:),ynan(:),'bo','filled'); 
p5 = scatter(xdat0(:),ynan(:),'ro','filled'); 
rl = refline(0,1);
rl.XData = xlims;
rl.YData = [1, 1];
rl.Color = c_gray;

lmcolorw = lmcolor;
lmcolorw(:) = {'k'};
plotDataScatter(xdat0,ydat1,xlabel_str,ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolorw,lmarker,l_filled);
lmcolorw(:) = {'c'};
plotDataScatter(xdat0,ydat2,xlabel_str,ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolorw,lmarker,l_filled);
lmcolorw(:) = {'b'};
plotDataScatter(xdat0,ydat3,xlabel_str,ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolorw,lmarker,l_filled);
lmcolorw(:) = {'r'};
plotDataScatter(xdat0,ydat4,xlabel_str,ylabel_str,...
                xlog,ylog,xlims,ylims,...
                lmcolorw,lmarker,l_filled);
            
lg = legend([p1,p4,p3,p5],'LMD94','LW16','RW16','This study');
lg.Location = 'NorthWest';
lg.Interpreter = 'Latex';

if (l_save_fig)
    figname = 'sct_HbKPP_wstar.fig';
    saveas(gcf,figname,'fig');
    postProcessFig(figname,6,4);
end
end

%% postprocess
% move figures to group dir         
if l_save_fig
    [s, r] = system(['mv ' fn_prefix '_* ' groupName]);
    if s~=0
        error('mv error');
    end
end
