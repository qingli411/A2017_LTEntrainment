close all; clear variables;
% load shared scripts
path(path,'../share');
path(path,'../share/cbrewer')

% color
global clightblue ccyan cmagenta cgray cgreen cyellow
clightblue = [0, 114, 189]./256; % light blue
ccyan = [0, 0.8, 0.8]; % dark cyan
cmagenta = [0.8, 0, 0.8]; % dark magenta
cgray = [0.5, 0.5, 0.5];    % gray
cgreen = [52, 168, 83]./256;    % dark green
cyellow = [251, 188, 5]./256;   % orange yellow

% case name
nameList = {'R8_BF05WD05WV00_ST01_ens03',...
            'R8_BF05WD05WV12_ST01_ens03',...
            'R8_BF05WD05WV00_fL00_ens03',...
            'R8_BF05WD05WV12_fL00_ens03'};
tStartList = [1,1,21601,1];
tEndList   = [32401,32401,25201,25201];
iziList = [64, 66, 69, 71];
% TKE budget xlims
tkex = [40, 40, 40, 40];
% Reynolds stress budget xlims
rsx11 = [15, 15, 30, 30];
rsx22 = [20, 20, 20, 20];
rsx33 = [15, 50, 15, 50];
rsx12 = [15,  5, 15,  5];
rsx13 = [40, 40, 40, 40];
rsx23 = [20, 50, 20, 50];
rsmult = [5, 5, 5, 5];
% Anisotropy budget xlims
aix11 = [3, 3, 8, 8];
aix22 = [3, 3, 3, 3];
aix33 = [3, 5, 5, 8];
aix12 = [2, 2, 1, 1];
aix13 = [5, 5, 8,12];
aix23 = [5, 5, 1, 1];
    
% select case
ic = 3;
casename = nameList{ic};
tStart = tStartList(ic);
tEnd = tEndList(ic);
izi = iziList(ic); 

% flags
l_save_fig = 0;     % save figure
l_tke = 1;          % plot TKE budget
l_Reynolds = 1;     % plot Reynolds stress budget
l_anisotropy = 0;   % plot anisotropy tensor budget
l_components = 1;   % save budget for each component of Reynolds stress
                    % and anisotropy tensor in individual figures
if ic == 1
    l_legend = 1;       % add legend
else
    l_legend = 0;
end

l_check_tke = 1;    % compare TKE directly calculated and that
                    % from Reynolds stress budget 

%% set parameters
get_dataRootDir;    % get dataRootDir and outRootDir

% set up profile data parameters
dataDir = [dataRootDir '/hist/'];
fPrefix = 'his.mp.vis';

% output directory
outDir = [outRootDir '/volume/' casename];
system(['mkdir -p ' outDir]);

% creating object ncarlesPflData
f = ncarlesPflData('caseName', casename,...
                   'dataDir',  dataDir,...
                   'fPrefix',  fPrefix,...
                   'tStart',   tStart,...
                   'tEnd',     tEnd);
par = f.getParameters;

% z levels
nnzcutoff = 128;
nnz = f.nnz;
nnzmax = min(nnzcutoff, nnz);
batag = f.alpha*f.g;
fcor = f.getVar('fcor');

% inertial period
% t_delta = 61094;
t_delta = 7200;
t_start = par.tend-t_delta;
t_end = par.tend;

% depth
zu = f.getVar('z_u')';
zw = f.getVar('z_w')';
dz = zu(2)-zu(1);

% get Stokes drift and stokes drift shear
tmp = f.getVar('stokes');
stokes = tmp(:, end)'; % Stokes drift at the last time step
clear tmp;
dir_x = f.getVar('dir_x');
dir_y = f.getVar('dir_y');
dstokesdzm = zeros(size(stokes));
dstokesdzm(1:end-1) = (stokes(1:end-1)-stokes(2:end))...
                    ./(zu(1:end-1)-zu(2:end));
dstokesdzm(end) = dstokesdzm(end-1);

% surface wind stress
tmp = f.getVar('utau');
utau = tmp(end);

% mean profile
if par.fcor > 0
    tmpv = f.getProfiles('inertial',[0,1]);
else
    tmpv = f.getProfiles('second',[par.tend-t_delta,par.tend]);
end
% boundary layer depth
hb = tmpv.mean_hb;
[~,ind_hb] = min(abs(zu+hb));

% weight
wgt = utau.^3./hb;
wgt2 = utau./hb;

%% calculate the budget profiles
pflind = 1:nnzmax;
pflind2 = pflind+1;
zpflu = zu(1:nnzmax)./hb;
zpflw = zu(1:nnzmax)./hb;
pflsize = size(zpflu);

% calculate the mean terms in Reynolds stress budget
% shear production (w level)
MSP11 = -2.*tmpv.uwle(pflind).*tmpv.dudz(pflind);
MSP22 = -2.*tmpv.vwle(pflind).*tmpv.dvdz(pflind);
MSP33 = zeros(pflsize);
MSP12 = -tmpv.uwle(pflind).*tmpv.dvdz(pflind)...
        -tmpv.vwle(pflind).*tmpv.dudz(pflind);
MSP13 = -tmpv.wps(pflind).*tmpv.dudz(pflind);
MSP23 = -tmpv.wps(pflind).*tmpv.dvdz(pflind);
% Stokes production (w level)
MStP11 = zeros(pflsize);
MStP22 = zeros(pflsize);
MStP33 = -2.*tmpv.uwle(pflind).*dstokesdzm(pflind).*dir_x...
         -2.*tmpv.vwle(pflind).*dstokesdzm(pflind).*dir_y;
MStP12 = zeros(pflsize);
tmp = tmpv.ups;
tmp2 = tmpv.uvle;
MStP13 = -0.5.*(tmp(pflind)+tmp(pflind2)).*dstokesdzm(pflind).*dir_x...
         -0.5.*(tmp2(pflind)+tmp2(pflind2)).*dstokesdzm(pflind).*dir_y;
tmp = tmpv.vps;
MStP23 = -0.5.*(tmp2(pflind)+tmp2(pflind2)).*dstokesdzm(pflind).*dir_x...
         -0.5.*(tmp(pflind)+tmp(pflind2)).*dstokesdzm(pflind).*dir_y;
% buoyancy production (w level)
MB11 = zeros(pflsize);
MB22 = zeros(pflsize);
MB33 = 2.*tmpv.wtle(pflind).*batag;
MB12 = zeros(pflsize);
tmp = tmpv.utle;
MB13 = 0.5.*(tmp(pflind)+tmp(pflind2)).*batag;
tmp = tmpv.vtle;
MB23 = 0.5.*(tmp(pflind)+tmp(pflind2)).*batag;
% TKE transport (u level)
tmp = zeros([1,nnz+1]);
tmp(2:end) = tmpv.uuwle;
tmp(1) = tmp(2);
MT11 = (tmp(pflind)-tmp(pflind2))./dz;
tmp(2:end) = tmpv.vvwle;
tmp(1) = tmp(2);
MT22 = (tmp(pflind)-tmp(pflind2))./dz;
tmp(2:end) = tmpv.wcube;
tmp(1) = tmp(2);
MT33 = (tmp(pflind)-tmp(pflind2))./dz;
tmp(2:end) = tmpv.uvwle;
tmp(1) = tmp(2);
MT12 = (tmp(pflind)-tmp(pflind2))./dz;
tmp(2:end) = tmpv.uwwle;
tmp(1) = tmp(2);
MT13 = (tmp(pflind)-tmp(pflind2))./dz;
tmp(2:end) = tmpv.vwwle;
tmp(1) = tmp(2);
MT23 = (tmp(pflind)-tmp(pflind2))./dz;
% pressure strain (w level)
MPS11 = -tmpv.udpdx(pflind)-tmpv.udpdx(pflind2);
MPS22 = -tmpv.vdpdy(pflind)-tmpv.vdpdy(pflind2);
MPS33 = -2.*tmpv.wdpdz(pflind);
MPS12 = -0.5.*(tmpv.udpdy(pflind)+tmpv.vdpdx(pflind)...
         +tmpv.udpdy(pflind2)+tmpv.vdpdx(pflind2));
MPS13 = -tmpv.udpdz(pflind)-0.5.*(tmpv.wdpdx(pflind)+tmpv.wdpdx(pflind2));
MPS23 = -tmpv.vdpdz(pflind)-0.5.*(tmpv.wdpdy(pflind)+tmpv.wdpdy(pflind2));
% Coriolis term (w level)
MC11 = fcor.*(tmpv.uvle(pflind)+tmpv.uvle(pflind2));
MC22 = -fcor.*(tmpv.uvle(pflind)+tmpv.uvle(pflind2));
MC33 = zeros(pflsize);
MC12 = 0.5.*fcor.*(tmpv.vps(pflind)+tmpv.vps(pflind2)...
         -tmpv.ups(pflind)-tmpv.ups(pflind2));
MC13 = fcor.*tmpv.vwle(pflind);
MC23 = fcor.*tmpv.uwle(pflind);
% transport of turbulent stress (w level)
tmp = tmpv.ttau11;
MTV11 =  -(tmp(pflind)-tmp(pflind2))./dz;
tmp = tmpv.ttau22;
MTV22 =  -(tmp(pflind)-tmp(pflind2))./dz;
tmp = tmpv.ttau33;
MTV33 =  -(tmp(pflind)-tmp(pflind2))./dz;
tmp = tmpv.ttau12;
MTV12 =  -(tmp(pflind)-tmp(pflind2))./dz;
tmp = tmpv.ttau13;
MTV13 =  -(tmp(pflind)-tmp(pflind2))./dz;
tmp = tmpv.ttau23;
MTV23 =  -(tmp(pflind)-tmp(pflind2))./dz;
% dissipation (u level)
MD11 = -tmpv.dsle11(pflind);
MD22 = -tmpv.dsle22(pflind);
MD33 = -tmpv.dsle33(pflind);
MD12 = -tmpv.dsle12(pflind);
MD13 = -tmpv.dsle13(pflind);
MD23 = -tmpv.dsle23(pflind);
% tendency
measure = 'second';
its = f.getTimeIndex(measure,t_start);
ite = f.getTimeIndex(measure,t_end);
tmp = f.getVar('ups');
MTD11 = (tmp(pflind,ite)-tmp(pflind,its))'./t_delta;
tmp = f.getVar('vps');
MTD22 = (tmp(pflind,ite)-tmp(pflind,its))'./t_delta;
tmp = f.getVar('wps');
MTD33 = (tmp(pflind,ite)-tmp(pflind,its))'./t_delta;
tmp = f.getVar('uvle');
MTD12 = (tmp(pflind,ite)-tmp(pflind,its))'./t_delta;
tmp = f.getVar('uwle');
MTD13 = (tmp(pflind,ite)-tmp(pflind,its))'./t_delta;
tmp = f.getVar('vwle');
MTD23 = (tmp(pflind,ite)-tmp(pflind,its))'./t_delta;

% TKE budget
MSP = 0.5.*(MSP11+MSP22+MSP33);
MStP = 0.5.*(MStP11+MStP22+MStP33);
MB = 0.5.*(MB11+MB22+MB33);
MT = 0.5.*(MT11+MT22+MT33);
MPS = 0.5.*(MPS11+MPS22+MPS33);
MTV = 0.5.*(MTV11+MTV22+MTV33);
MC = 0.5.*(MC11+MC22+MC33);
MD = 0.5.*(MD11+MD22+MD33);
MTD = 0.5*(MTD11+MTD22+MTD33);

% rearrange pressure stress and TKE transport
MT11 = MT11 + 2.*MPS./3;
MT22 = MT22 + 2.*MPS./3;
MT33 = MT33 + 2.*MPS./3;
MPS11 = MPS11 - 2.*MPS./3;
MPS22 = MPS22 - 2.*MPS./3;
MPS33 = MPS33 - 2.*MPS./3;
if l_check_tke
    % check consistency
    MSPb = tmpv.t_rprod;
    MStPb = tmpv.t_stokes;
    MBb = tmpv.wtle.*batag;
    MTb = tmpv.t_wq;
    MPSb = tmpv.t_wp;
    MTVb = (tmpv.t_tran-tmpv.t_wq-tmpv.t_wp);
    MCb = zeros(pflsize);
    MDb = -tmpv.t_dsle;
else
    MT = MT + MPS;
    MPS = MPS - MPS;
end
% a_{ij} budget
TKE = tmpv.englez(pflind);  % resolved TKE at w level
ASP11 = (MSP11-MSP./3)./2./TKE;
ASP22 = (MSP22-MSP./3)./2./TKE;
ASP33 = (MSP33-MSP./3)./2./TKE;
ASP12 = MSP12./2./TKE;
ASP13 = MSP13./2./TKE;
ASP23 = MSP23./2./TKE;

AStP11 = (MStP11-MStP./3)./2./TKE;
AStP22 = (MStP22-MStP./3)./2./TKE;
AStP33 = (MStP33-MStP./3)./2./TKE;
AStP12 = MStP12./2./TKE;
AStP13 = MStP13./2./TKE;
AStP23 = MStP23./2./TKE;

AB11 = (MB11-MB./3)./2./TKE;
AB22 = (MB22-MB./3)./2./TKE;
AB33 = (MB33-MB./3)./2./TKE;
AB12 = MB12./2./TKE;
AB13 = MB13./2./TKE;
AB23 = MB23./2./TKE;

AT11 = (MT11-MT./3)./2./TKE;
AT22 = (MT22-MT./3)./2./TKE;
AT33 = (MT33-MT./3)./2./TKE;
AT12 = MT12./2./TKE;
AT13 = MT13./2./TKE;
AT23 = MT23./2./TKE;

APS11 = (MPS11-MPS./3)./2./TKE;
APS22 = (MPS22-MPS./3)./2./TKE;
APS33 = (MPS33-MPS./3)./2./TKE;
APS12 = MPS12./2./TKE;
APS13 = MPS13./2./TKE;
APS23 = MPS23./2./TKE;

AC11 = (MC11-MC./3)./2./TKE;
AC22 = (MC22-MC./3)./2./TKE;
AC33 = (MC33-MC./3)./2./TKE;
AC12 = MC12./2./TKE;
AC13 = MC13./2./TKE;
AC23 = MC23./2./TKE;

ATV11 = (MTV11-MTV./3)./2./TKE;
ATV22 = (MTV22-MTV./3)./2./TKE;
ATV33 = (MTV33-MTV./3)./2./TKE;
ATV12 = MTV12./2./TKE;
ATV13 = MTV13./2./TKE;
ATV23 = MTV23./2./TKE;

AD11 = (MD11-MD./3)./2./TKE;
AD22 = (MD22-MD./3)./2./TKE;
AD33 = (MD33-MD./3)./2./TKE;
AD12 = MD12./2./TKE;
AD13 = MD13./2./TKE;
AD23 = MD23./2./TKE;

% anisotropy tensor
aij = f.getAnisotropicTensor(tmpv);
A11 = squeeze(aij(pflind,1,1))';
A12 = squeeze(aij(pflind,1,2))';
A13 = squeeze(aij(pflind,1,3))';
A22 = squeeze(aij(pflind,2,2))';
A23 = squeeze(aij(pflind,2,3))';
A33 = squeeze(aij(pflind,3,3))';

%% figure 1: mean TKE budget
if l_tke
    figure;
    hold on;
    mult = 10;
    indz = 2:ind_hb;
    zzu = zpflu(indz);
    zzw = zpflw(indz);
    p1 = pflplot2(MSP(indz)./wgt,  zzw, mult, '-b');
    p2 = pflplot2(MStP(indz)./wgt, zzw, mult, '-r');
    p3 = pflplot2(MB(indz)./wgt,   zzw, mult, '-k');
    p4 = pflplot2(MT(indz)./wgt,   zzw, mult, '-', 'Color', ccyan);
    p6 = pflplot2(MTV(indz)./wgt,  zzw, mult, '-', 'Color', cyellow);
    p7 = pflplot2(MD(indz)./wgt,   zzw, mult, '-', 'Color', cmagenta);
    if l_check_tke
        % check consistency
        p5 = pflplot2(MPS(indz)./wgt,  zzw, mult, '-', 'Color', cgreen);
        MTOT = MSP+MStP+MB+MT+MTV+MD+MPS+MTD;
        pflplot2(MTOT(indz)./wgt,  zzw, mult, '-', 'Color', cgray);
        pflplot2(MSPb(indz)./wgt,  zzw, mult, '--b');
        pflplot2(MStPb(indz)./wgt, zzw, mult, '--r');
        pflplot2(MBb(indz)./wgt,   zzw, mult, '--k');
        pflplot2(MTb(indz)./wgt,   zzu, mult, '--', 'Color', ccyan);
        pflplot2(MPSb(indz)./wgt,  zzw, mult, '--', 'Color', cgreen);
        pflplot2(MTVb(indz)./wgt,  zzw, mult, '--', 'Color', cyellow);
        pflplot2(MDb(indz)./wgt,   zzu, mult, '--', 'Color', cmagenta);
    end
    ylim([-1, 0]);
    xlim([-tkex(ic), tkex(ic)]);
    pflref2(zzu);
    pflrefentr(zw(izi)/hb);
    ylabel('$u_3/h_b$','Interpreter','latex');
    xlabel('TKE budget','Interpreter','latex');
    if l_legend
        lg = legend([p1, p2, p3, p4, p6, p7],...
            '${\cal P}$', '${\cal P}^S$',...
            '${\cal B}$', '${\cal T}$',...
            '${\cal T}^v$', '${\cal D}$');
        lg.Location = 'SouthWest';
        lg.Interpreter = 'latex';
    end

    if l_save_fig
        if ~l_legend
            figname = [outDir '/tke_budget_nolegend.fig'];
        else
            figname = [outDir '/tke_budget.fig'];
        end
        saveas(gcf, figname, 'fig');
        postProcessFig(figname, 4, 4);
    end
end

%% figure 2: mean Reynolds stress budget
if l_Reynolds
    mult = rsmult(ic);
    if l_components
        figure;
        xlims = [-rsx11(ic),rsx11(ic)];
        xlabel_str = '$\overline{u''_1 u''_1}$ budget';
        p = plot_reynolds_stress(MSP11, MStP11, MB11, MT11,...
            MPS11, MC11, MTV11, MD11,...
            zpflu, zpflw, ind_hb, izi, mult, wgt, xlims, xlabel_str);
        if l_legend
            lg = legend(p,...
                '${\cal P}_{ij}$', '${\cal P}^S_{ij}$',...
                '${\cal B}_{ij}$', '${\cal T}_{ij}$',...
                '$\Pi_{ij}$','${\cal C}_{ij}$',...
                '${\cal T}^v_{ij}$', '${\cal D}_{ij}$');
            lg.Location = 'SouthEast';
            lg.Interpreter = 'latex';
        end
        if l_save_fig
            figname = [outDir '/ReynoldsStress11_budget.fig'];
            saveas(gcf, figname, 'fig');
            postProcessFig(figname, 4, 4);
        end
        figure;
        xlims = [-rsx22(ic),rsx22(ic)];
        xlabel_str = '$\overline{u''_2 u''_2}$ budget';
        plot_reynolds_stress(MSP22, MStP22, MB22, MT22,...
            MPS22, MC22, MTV22, MD22,...
            zpflu, zpflw, ind_hb, izi, mult, wgt, xlims, xlabel_str);
        if l_save_fig
            figname = [outDir '/ReynoldsStress22_budget.fig'];
            saveas(gcf, figname, 'fig');
            postProcessFig(figname, 4, 4);
        end
        figure;
        xlims = [-rsx33(ic),rsx33(ic)];
        xlabel_str = '$\overline{u''_3 u''_3}$ budget';
        plot_reynolds_stress(MSP33, MStP33, MB33, MT33,...
            MPS33, MC33, MTV33, MD33,...
            zpflu, zpflw, ind_hb, izi, mult, wgt, xlims, xlabel_str);
        if l_save_fig
            figname = [outDir '/ReynoldsStress33_budget.fig'];
            saveas(gcf, figname, 'fig');
            postProcessFig(figname, 4, 4);
        end
        figure
        xlims = [-rsx12(ic),rsx12(ic)];
        xlabel_str = '$\overline{u''_1 u''_2}$ budget';
        plot_reynolds_stress(MSP12, MStP12, MB12, MT12,...
            MPS12, MC12, MTV12, MD12,...
            zpflu, zpflw, ind_hb, izi, mult, wgt, xlims, xlabel_str);
        if l_save_fig
            figname = [outDir '/ReynoldsStress12_budget.fig'];
            saveas(gcf, figname, 'fig');
            postProcessFig(figname, 4, 4);
        end
        figure;
        xlims = [-rsx13(ic),rsx13(ic)];
        xlabel_str = '$\overline{u''_1 u''_3}$ budget';
        plot_reynolds_stress(MSP13, MStP13, MB13, MT13,...
            MPS13, MC13, MTV13, MD13,...
            zpflu, zpflw, ind_hb, izi, mult, wgt, xlims, xlabel_str);
        if l_save_fig
            figname = [outDir '/ReynoldsStress13_budget.fig'];
            saveas(gcf, figname, 'fig');
            postProcessFig(figname, 4, 4);
        end
        figure;
        xlims = [-rsx23(ic),rsx23(ic)];
        xlabel_str = '$\overline{u''_2 u''_3}$ budget';
        plot_reynolds_stress(MSP23, MStP23, MB23, MT23,...
            MPS23, MC23, MTV23, MD23,...
            zpflu, zpflw, ind_hb, izi, mult, wgt, xlims, xlabel_str);
        if l_save_fig
            figname = [outDir '/ReynoldsStress23_budget.fig'];
            saveas(gcf, figname, 'fig');
            postProcessFig(figname, 4, 4);
        end
    else
        fig = figure;
        fig.Units = 'inches';
        fig.Position = [1 1 12 8];
        subplot(2, 3, 1)
        xlims = [-rsx11(ic),rsx11(ic)];
        xlabel_str = '$\overline{u''_1 u''_1}$ budget';
        plot_reynolds_stress(MSP11, MStP11, MB11, MT11,...
            MPS11, MC11, MTV11, MD11,...
            zpflu, zpflw, ind_hb, izi, mult, wgt, xlims, xlabel_str);
        subplot(2, 3, 2)
        xlims = [-rsx22(ic),rsx22(ic)];
        xlabel_str = '$\overline{u''_2 u''_2}$ budget';
        plot_reynolds_stress(MSP22, MStP22, MB22, MT22,...
            MPS22, MC22, MTV22, MD22,...
            zpflu, zpflw, ind_hb, izi, mult, wgt, xlims, xlabel_str);
        subplot(2, 3, 3)
        xlims = [-rsx33(ic),rsx33(ic)];
        xlabel_str = '$\overline{u''_3 u''_3}$ budget';
        plot_reynolds_stress(MSP33, MStP33, MB33, MT33,...
            MPS33, MC33, MTV33, MD33,...
            zpflu, zpflw, ind_hb, izi, mult, wgt, xlims, xlabel_str);
        subplot(2, 3, 4)
        xlims = [-rsx12(ic),rsx12(ic)];
        xlabel_str = '$\overline{u''_1 u''_2}$ budget';
        plot_reynolds_stress(MSP12, MStP12, MB12, MT12,...
            MPS12, MC12, MTV12, MD12,...
            zpflu, zpflw, ind_hb, izi, mult, wgt, xlims, xlabel_str);
        subplot(2, 3, 5)
        xlims = [-rsx13(ic),rsx13(ic)];
        xlabel_str = '$\overline{u''_1 u''_3}$ budget';
        p=plot_reynolds_stress(MSP13, MStP13, MB13, MT13,...
            MPS13, MC13, MTV13, MD13,...
            zpflu, zpflw, ind_hb, izi, mult, wgt, xlims, xlabel_str);
        lg = legend(p,...
            '${\cal P}_{ij}$', '${\cal P}^S_{ij}$',...
            '${\cal B}_{ij}$', '${\cal T}_{ij}$',...
            '$\Pi_{ij}$','${\cal C}_{ij}$',...
            '${\cal T}^v_{ij}$', '${\cal D}_{ij}$');
        lg.Location = 'SouthEast';
        lg.Interpreter = 'latex';
        subplot(2, 3, 6)
        xlims = [-rsx23(ic),rsx23(ic)];
        xlabel_str = '$\overline{u''_2 u''_3}$ budget';
        plot_reynolds_stress(MSP23, MStP23, MB23, MT23,...
            MPS23, MC23, MTV23, MD23,...
            zpflu, zpflw, ind_hb, izi, mult, wgt, xlims, xlabel_str);
    %     set(findall(gcf,'-property','FontSize'),'FontSize',12)
        if l_save_fig
            figname = [outDir '/ReynoldsStress_budget.fig'];
            saveas(gcf, figname, 'fig');
            [dir, name, ~] = fileparts(figname);
            print('-depsc2',[dir '/' name]);
        end
    end
end

%% figure 3: mean anisotropy tensor budget
if l_anisotropy
    if l_components
        figure;
        xlims = [-aix11(ic),aix11(ic)];
        xlabel_str = '$a_{11}$ budget';
        p = plot_anisotropy(ASP11, AStP11, AB11, AT11, APS11, AD11,...
            A11, TKE, MSP, MStP, MB, MT, MPS, MD,...
            zpflu, zpflw, ind_hb, izi, wgt2, xlims, xlabel_str);
        if l_legend
            lg = legend(p,...
                '${\cal P}_{ij}$', '${\cal P}^S_{ij}$',...
                '${\cal B}_{ij}$', '${\cal T}_{ij}$',...
                '$\Pi_{ij}$', '${\cal D}_{ij}$');
            lg.Location = 'East';
            lg.Interpreter = 'latex';
        end
        if l_save_fig
            figname = [outDir '/Anisotropy11_budget.fig'];
            saveas(gcf, figname, 'fig');
            postProcessFig(figname, 4, 4);
        end
        figure;
        xlims = [-aix22(ic),aix22(ic)];
        xlabel_str = '$a_{22}$ budget';
        plot_anisotropy(ASP22, AStP22, AB22, AT22, APS22, AD22,...
            A22, TKE, MSP, MStP, MB, MT, MPS, MD,...
            zpflu, zpflw, ind_hb, izi, wgt2, xlims, xlabel_str);
        if l_save_fig
            figname = [outDir '/Anisotropy22_budget.fig'];
            saveas(gcf, figname, 'fig');
            postProcessFig(figname, 4, 4);
        end
        figure;
        xlims = [-aix33(ic),aix33(ic)];
        xlabel_str = '$a_{33}$ budget';
        plot_anisotropy(ASP33, AStP33, AB33, AT33, APS33, AD33,...
            A33, TKE, MSP, MStP, MB, MT, MPS, MD,...
            zpflu, zpflw, ind_hb, izi, wgt2, xlims, xlabel_str);
        if l_save_fig
            figname = [outDir '/Anisotropy33_budget.fig'];
            saveas(gcf, figname, 'fig');
            postProcessFig(figname, 4, 4);
        end
        figure;
        xlims = [-aix12(ic),aix12(ic)];
        xlabel_str = '$a_{12}$ budget';
        plot_anisotropy(ASP12, AStP12, AB12, AT12, APS12, AD12,...
            A12, TKE, MSP, MStP, MB, MT, MPS, MD,...
            zpflu, zpflw, ind_hb, izi, wgt2, xlims, xlabel_str);
        if l_save_fig
            figname = [outDir '/Anisotropy12_budget.fig'];
            saveas(gcf, figname, 'fig');
            postProcessFig(figname, 4, 4);
        end
        figure;
        xlims = [-aix13(ic),aix13(ic)];
        xlabel_str = '$a_{13}$ budget';
        p=plot_anisotropy(ASP13, AStP13, AB13, AT13, APS13, AD13,...
            A13, TKE, MSP, MStP, MB, MT, MPS, MD,...
            zpflu, zpflw, ind_hb, izi, wgt2, xlims, xlabel_str);
        if l_save_fig
            figname = [outDir '/Anisotropy13_budget.fig'];
            saveas(gcf, figname, 'fig');
            postProcessFig(figname, 4, 4);
        end
        figure;
        xlims = [-aix23(ic),aix23(ic)];
        xlabel_str = '$a_{23}$ budget';
        plot_anisotropy(ASP23, AStP23, AB23, AT23, APS23, AD23,...
            A23, TKE, MSP, MStP, MB, MT, MPS, MD,...
            zpflu, zpflw, ind_hb, izi, wgt2, xlims, xlabel_str);
        if l_save_fig
            figname = [outDir '/Anisotropy23_budget.fig'];
            saveas(gcf, figname, 'fig');
            postProcessFig(figname, 4, 4);
        end
    else
        fig = figure;
        fig.Units = 'inches';
        fig.Position = [1 1 12 8];
        subplot(2, 3, 1)
        xlims = [-aix11(ic),aix11(ic)];
        xlabel_str = '$a_{11}$ budget';
        plot_anisotropy(ASP11, AStP11, AB11, AT11, APS11, AD11,...
            A11, TKE, MSP, MStP, MB, MT, MPS, MD,...
            zpflu, zpflw, ind_hb, izi, wgt2, xlims, xlabel_str);
        subplot(2, 3, 2)
        xlims = [-aix22(ic),aix22(ic)];
        xlabel_str = '$a_{22}$ budget';
        plot_anisotropy(ASP22, AStP22, AB22, AT22, APS22, AD22,...
            A22, TKE, MSP, MStP, MB, MT, MPS, MD,...
            zpflu, zpflw, ind_hb, izi, wgt2, xlims, xlabel_str);
        subplot(2, 3, 3)
        xlims = [-aix33(ic),aix33(ic)];
        xlabel_str = '$a_{33}$ budget';
        plot_anisotropy(ASP33, AStP33, AB33, AT33, APS33, AD33,...
            A33, TKE, MSP, MStP, MB, MT, MPS, MD,...
            zpflu, zpflw, ind_hb, izi, wgt2, xlims, xlabel_str);
        subplot(2, 3, 4)
        xlims = [-aix12(ic),aix12(ic)];
        xlabel_str = '$a_{12}$ budget';
        plot_anisotropy(ASP12, AStP12, AB12, AT12, APS12, AD12,...
            A12, TKE, MSP, MStP, MB, MT, MPS, MD,...
            zpflu, zpflw, ind_hb, izi, wgt2, xlims, xlabel_str);
        subplot(2, 3, 5)
        xlims = [-aix13(ic),aix13(ic)];
        xlabel_str = '$a_{13}$ budget';
        p=plot_anisotropy(ASP13, AStP13, AB13, AT13, APS13, AD13,...
            A13, TKE, MSP, MStP, MB, MT, MPS, MD,...
            zpflu, zpflw, ind_hb, izi, wgt2, xlims, xlabel_str);
        lg = legend(p,...
            '${\cal P}_{ij}$', '${\cal P}^S_{ij}$',...
            '${\cal B}_{ij}$', '${\cal T}_{ij}$',...
            '$\Pi_{ij}$', '${\cal D}_{ij}$');
        lg.Location = 'East';
        lg.Interpreter = 'latex';
        subplot(2, 3, 6)
        xlims = [-aix23(ic),aix23(ic)];
        xlabel_str = '$a_{23}$ budget';
        plot_anisotropy(ASP23, AStP23, AB23, AT23, APS23, AD23,...
            A23, TKE, MSP, MStP, MB, MT, MPS, MD,...
            zpflu, zpflw, ind_hb, izi, wgt2, xlims, xlabel_str);
    %     set(findall(gcf,'-property','FontSize'),'FontSize',12)
        if l_save_fig
            figname = [outDir '/Anisotropy_budget.fig'];
            saveas(gcf, figname, 'fig');
            [dir, name, ~] = fileparts(figname);
            print('-depsc2',[dir '/' name]);
        end
    end
end

%% functions

function p = pflplot3(xx, zz, varargin)
    nz = numel(zz);
    qnz = floor(nz/4);
    tqnz = floor(nz*3/4);
    indz1 = 1:qnz;
    indz2 = qnz:tqnz;
    indz3 = tqnz:nz;
    hold on;
    p = plot(xx(indz1),zz(indz1),varargin{:});
    p.LineWidth = 1.5;
    p2 = plot(xx(indz2).*10,zz(indz2),varargin{:});
    p2.LineWidth = 1.5;
    p3 = plot(xx(indz3).*25,zz(indz3),varargin{:});
    p3.LineWidth = 1.5;
end

function p = pflplot2(xx, zz, mult, varargin)
    nz = numel(zz);
    hnz = floor(nz/2);
    indz1 = 1:hnz;
    indz2 = hnz:nz;
    hold on;
    p = plot(xx(indz1),zz(indz1),varargin{:});
    p.LineWidth = 1.5;
    p2 = plot(xx(indz2).*mult,zz(indz2),varargin{:});
    p2.LineWidth = 1.5;
end

function rl = pflref3(zz)
    nz = numel(zz);
    qnz = floor(nz/4);
    tqnz = floor(nz*3/4);
    rl = refline(0, zz(qnz));
    rl.Color = 'k';
    rl.LineStyle = '-';
    rl.LineWidth = 0.75;
    rl2 = refline(0, zz(tqnz));
    rl2.Color = 'k';
    rl2.LineStyle = '-';
    rl2.LineWidth = 0.75;
end

function rl = pflref2(zz)
    nz = numel(zz);
    hnz = floor(nz/2);
    rl = refline(0, zz(hnz));
    rl.Color = 'k';
    rl.LineStyle = '-';
    rl.LineWidth = 0.75;
end

function rl = pflrefentr(zz)
    global cgray
    rl = refline(0, zz);
    rl.Color = cgray;
    rl.LineStyle = '--';
    rl.LineWidth = 0.75;
end

function p = plot_reynolds_stress(SP, StP, B, T, PS, C, TV, D,...
    zu, zw, ihb, izi, mult, wgt, xlims, xlabel_str)
    global ccyan cgreen cyellow cgray cmagenta
    inds = 2:ihb;
    p = gobjects(1,8);
    p(1)=pflplot2(SP(inds)./wgt, zw(inds), mult, '-b');
    p(2)=pflplot2(StP(inds)./wgt, zw(inds), mult, '-r');
    p(3)=pflplot2(B(inds)./wgt, zw(inds), mult, '-k');
    p(4)=pflplot2(T(inds)./wgt, zu(inds), mult, '-', 'Color', ccyan);
    p(5)=pflplot2(PS(inds)./wgt, zw(inds), mult, '-', 'Color', cgreen);
    p(6)=pflplot2(C(inds)./wgt, zw(inds), mult, '-', 'Color', cgray);
    p(7)=pflplot2(TV(inds)./wgt, zw(inds), mult, '-', 'Color', cyellow);
    p(8)=pflplot2(D(inds)./wgt, zu(inds), mult, '-', 'Color', cmagenta);
    ylim([-1, 0]);
    set(gca,'fontsize',12);
    ylabel('$x_3/h_b$','Interpreter','latex',...
        'FontSize', 14);
    xlim(xlims);
    xlabel(xlabel_str,'Interpreter','latex',...
        'FontSize', 14);
    pflref2(zu(inds));
    if izi > 0
        pflrefentr(zw(izi));
    end
end

function p = plot_anisotropy(SP, StP, B, T, PS, D, A, TKE,...
    MSP, MStP, MB, MT, MPS, MD, zu, zw, ihb, izi, wgt, xlims, xlabel_str)
    global ccyan cgreen cmagenta
    inds = 2:ihb;
    mult = 3;
    p = gobjects(1,6);
    hold on;
    p(1)=pflplot2(SP(inds)./wgt, zw(inds), mult, '-b');
    p(2)=pflplot2(StP(inds)./wgt, zw(inds), mult, '-r');
    p(3)=pflplot2(B(inds)./wgt, zw(inds), mult, '-k');
    p(4)=pflplot2(T(inds)./wgt, zu(inds), mult, '-', 'Color', ccyan);
    p(5)=pflplot2(PS(inds)./wgt, zw(inds), mult, '-', 'Color', cgreen);
    p(6)=pflplot2(D(inds)./wgt, zu(inds), mult, '-', 'Color', cmagenta);
    pflplot2(-A(inds)./TKE(inds).*MSP(inds)./wgt, zw(inds), mult, '--b');
    pflplot2(-A(inds)./TKE(inds).*MStP(inds)./wgt, zw(inds), mult, '--r');
    pflplot2(-A(inds)./TKE(inds).*MB(inds)./wgt, zw(inds), mult, '--k');
    pflplot2(-A(inds)./TKE(inds).*MT(inds)./wgt, zu(inds), mult,...
        '--', 'Color', ccyan);
    pflplot2(-A(inds)./TKE(inds).*MPS(inds)./wgt, zw(inds), mult,...
        '--', 'Color', cgreen);
    pflplot2(-A(inds)./TKE(inds).*MD(inds)./wgt, zu(inds), mult,...
        '--', 'Color', cmagenta);
    ylim([-1, 0]);
    set(gca,'fontsize',12);
    ylabel('$x_3/h_b$','Interpreter','latex',...
        'FontSize', 14);
    xlim(xlims);
    xlabel(xlabel_str,'Interpreter','latex',...
        'FontSize', 14);
    pflref2(zu(inds));
    if izi > 0
        pflrefentr(zw(izi));
    end
end
