 close all; clear variables;
% This script plots profiles of a group of simulations

% Qing Li, 170117

% add path if necessary
if exist('qingAddPath.m','file') == 2
    qingAddPath
end
% add path for ncarleePflData
path(path,'../share');
% save figures or not (0 for testing)
l_save_fig = 1;
l_h = 2; % 1:h (buoyancy flux minimum),
         % 2:hm (N2 maximum),
         % 3:hb (dissipation threshold, 1e-9)
         
%% Set cases
% get profile group
groupName = 'Langmuir';
% groupName = 'Langmuir2';
% groupName = 'Langmuir3';
% groupName = 'LangmuirIO';
% groupName = 'Res_Scale2';
% groupName = 'Convection';
% groupName = 'LowerBC';
% groupName = 'Coriolis';
% groupName = 'FCutoff';

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

%% Get data
nsize = size(caselist);
nc = nsize(1);
nnz = 256;  % maximum nz
epsilonb = 1e-9; % dissipation rate threshold to determine the hb
% gray
c_gray = [0.5 0.5 0.5];
% data dir
dir = '../AnalyzeSingleCase';
wps = zeros(nc,nnz);
wb = zeros(nc,nnz);
diss = zeros(nc,nnz);
tkeTrans = zeros(nc,nnz);
shearProd = zeros(nc,nnz);
stokesProd = zeros(nc,nnz);
stokes = zeros(nc,nnz);
tke = zeros(nc,nnz);
N2 = zeros(nc,nnz);
zw = zeros(nc,nnz);
zu = zeros(nc,nnz);
stokessl = zeros(nc,1);
h = zeros(nc,1);
hm = zeros(nc,1);
hb = zeros(nc,1);
utau = zeros(nc,1);
b0 = zeros(nc,1);
for j=1:nc
    casename = caselist(j,:);
    fileName = [dir '/' casename '/MeanProfile.mat'];
    if exist(fileName,'file')
        fprintf('Case Name: %s; Loading...\n', casename);
        dat = load(fileName);
        nzl = numel(dat.z_w);
        wps(j,1:nzl) = dat.wps;
        diss(j,1:nzl) = dat.dissip;
        tkeTrans(j,1:nzl) = dat.tkeTrans;
        shearProd(j,1:nzl) = dat.shearProd;
        stokesProd(j,1:nzl) = dat.stokesProd;
        stokes(j,1:nzl) = dat.stokes;
        tke(j,1:nzl) = dat.tke(1,:);
        wb(j,1:nzl) = dat.wb(1,:);
        N2(j,1:nzl) = dat.N2;
        zw(j,1:nzl) = dat.z_w;
        zu(j,1:nzl) = dat.z_u;
        utau(j,1) = dat.utau;
        b0(j,1) = dat.wtsfc.*dat.f.alpha.*dat.f.g;
        stokessl(j) = dat.f.getUsSL(dat.stokes,dat.z_u,dat.hm);
        h(j,1) = dat.h;
        hm(j,1) = dat.hm;
        if isfield(dat,'hb')
            hb(j,1) = dat.hb;
        else
            hb(j,1) = -getThresholdDepth(diss(j,:),zw(j,:),epsilonb);
        end
    else
        fprintf('Case Name: %s; Not exist, skipping...\n', casename);
        wps(j,:) = NaN;
        diss(j,:) = NaN;
        tkeTrans(j,:) = NaN;
        shearProd(j,:) = NaN;
        stokesProd(j,:) = NaN;
        stokes(j,:) = NaN;
        tke(j,:) = NaN;
        N2(j,:) = NaN;
        wb(j,:) = NaN;
        zw(j,:) = NaN;
        zu(j,:) = NaN;
        utau(j,1) = NaN;
        b0(j,1) = NaN;
        stokessl(j) = NaN;
        h(j,1) = NaN;
        hm(j,1) = NaN;
        hb(j,1) = NaN;
    end
end

%% plot figures
fn_prefix = 'pfl'; % figure name prefix
lh_error_str = 'l_h = 1 or 2 or 3';
ylab = '$z/h_\mathrm{b}$';
if l_h == 1
    hml = h;
    fn_suffix = '_h'; % figure name suffix
elseif l_h == 2
    hml = hm;
    fn_suffix = '_hm';
elseif l_h == 3
    hml = hb;
    fn_suffix = '_hb';
else
    error(lh_error_str);
end

ylims = [-1.2,0];
% wps
newFigure(l_save_fig);
datpfl = wps;
wgt = utau.^2;
xlab = '$\overline{w''^2}/{u^*}^2$';
zz = zw;
figname = [fn_prefix '_wps' fn_suffix];
xlog = 0;
xlims = [];
plotDataProfile(datpfl,wgt,xlab,ylab,xlims,ylims,xlog,...
                zz,hml,lcolor,[]);
saveFigure(l_save_fig, figname);
   
% wps scaled by w*^2
newFigure(l_save_fig);
datpfl = wps;
wgt = (b0.*hml).^(2/3);
xlab = '$\overline{w''^2}/{w^*}^2$';
zz = zw;
figname = [fn_prefix '_wps2' fn_suffix];
xlog = 0;
xlims = [];
plotDataProfile(datpfl,wgt,xlab,ylab,xlims,ylims,xlog,...
                zz,hml,lcolor,[]);
saveFigure(l_save_fig, figname);

% wb scaled by u*^3/hx
newFigure(l_save_fig);
datpfl = wb;
zz = zw;
xlab = '$\overline{w''b''}h_\mathrm{b}/{u^*}^3$';
wgt = utau.^3./hml;
figname = [fn_prefix '_wb' fn_suffix];
xlog = 0;
xlims = [];
if strcmp(groupName, 'Langmuir')
    xlims = [-0.8,0.8];
end
plotDataProfile(datpfl,wgt,xlab,ylab,xlims,ylims,xlog,...
                zz,hml,lcolor,[]);

if strcmp(groupName,'Convection')
    % add inlet
    % create a new pair of axes inside current figure
    axes('position',[.52 .2 .35 .35]);
    xlims1 = [-15,3];
    ylims1 = [-1.1,-0.7];
    % box on; % put box around new pair of axes
    plotDataProfile(datpfl,wgt,[],[],xlims1,ylims1,xlog,...
                zz,hml,lcolor,[]);
end
saveFigure(l_save_fig, figname);
   
% wb scaled by B0
newFigure(l_save_fig);
datpfl = wb;
zz = zw;
wgt = b0;
xlab = '$\overline{w''b''}/B_0$';
figname = [fn_prefix '_wb2' fn_suffix];
xlog = 0;
xlims = [];
plotDataProfile(datpfl,wgt,xlab,ylab,xlims,ylims,xlog,...
                zz,hml,lcolor,[]);
saveFigure(l_save_fig, figname);

% N2
newFigure(l_save_fig);
datpfl = N2;
zz = zu;
wgt = ones(size(datpfl));
xlab = '$N^2$';
figname = [fn_prefix '_N2' fn_suffix];
xlog = 0;
xlims = [];
plotDataProfile(datpfl,wgt,xlab,ylab,xlims,ylims,xlog,...
                zz,hml,lcolor,[]);
saveFigure(l_save_fig, figname);

% tke
newFigure(l_save_fig);
datpfl = tke;
zz = zw;
wgt = utau.^2;
xlab = 'Normalized TKE (${u^*}^2$) ';
figname = [fn_prefix '_tke' fn_suffix];
xlog = 0;
xlims = [];
plotDataProfile(datpfl,wgt,xlab,ylab,xlims,ylims,xlog,...
                zz,hml,lcolor,[]);
saveFigure(l_save_fig, figname);

% dissipation
newFigure(l_save_fig);
datpfl = diss;
zz = zw;
xlab = 'Normalized dissipation (${u^*}^3/h_\mathrm{b}$) ';
wgt = utau.^3./hml;
figname = [fn_prefix '_diss' fn_suffix];
xlog = 1;
xlims = [1e-2, 1e2];
plotDataProfile(datpfl,wgt,xlab,ylab,xlims,ylims,xlog,...
                zz,hml,lcolor,[]);
saveFigure(l_save_fig, figname);

% dissipation, scaled by B0
newFigure(l_save_fig);
datpfl = diss;
zz = zw;
xlab = 'Normalized dissipation ($B_0$) ';
wgt = b0;
figname = [fn_prefix '_diss2' fn_suffix];
xlog = 1;
xlims = [1e-2, 1e2];
plotDataProfile(datpfl,wgt,xlab,ylab,xlims,ylims,xlog,...
                zz,hml,lcolor,[]);
saveFigure(l_save_fig, figname);

% TKE transport
newFigure(l_save_fig);
datpfl = tkeTrans;
zz = zw;
xlab = 'Normalized TKE transport (${u^*}^3/h_\mathrm{b}$)';
wgt = utau.^3./hml;
figname = [fn_prefix '_tkeTrans' fn_suffix];
xlog = 1;
xlims = [1e-2, 1e2];
plotDataProfile(datpfl,wgt,xlab,ylab,xlims,ylims,xlog,...
                zz,hml,lcolor,[]);
saveFigure(l_save_fig, figname);

% shear production
newFigure(l_save_fig);
datpfl = shearProd;
zz = zw;
xlab = 'Normalized shear production (${u^*}^3/h_\mathrm{b}$)';
wgt = utau.^3./hml;
figname = [fn_prefix '_shearProd' fn_suffix];
xlog = 1;
xlims = [1e-2, 1e2];
plotDataProfile(datpfl,wgt,xlab,ylab,xlims,ylims,xlog,...
                zz,hml,lcolor,[]);
saveFigure(l_save_fig, figname);
            
% stokes + shear production or Lagrangian shear production
if strcmp(groupName, 'Langmuir')
    newFigure(l_save_fig);
    datpfl = shearProd+stokesProd;
    zz = zw;
    xlab = 'Normalized Lagrangian shear production (${u^*}^3/h_\mathrm{b}$)';
    wgt = utau.^3./hml;
    figname = [fn_prefix '_shearProdLagr' fn_suffix];
    xlog = 1;
    xlims = [1e-2, 1e2];
    plotDataProfile(datpfl,wgt,xlab,ylab,xlims,ylims,xlog,...
                    zz,hml,lcolor,[]);          
end
saveFigure(l_save_fig, figname);

% stokes drift
newFigure(l_save_fig);
datpfl = stokes;
zz = zu;
wgt = utau;
xlab = '$u^\mathrm{S}/{u^*}$';
figname = [fn_prefix '_stokes' fn_suffix];
xlog = 0;
xlims = [0, 12];
xx = zeros(nc,2);
yy = zeros(nc,2);
xx(:,1) = stokessl./wgt;
xx(:,2) = stokessl./wgt;
yy(:,1) = -ones(size(stokessl)).*0.64;
yy(:,2) = yy(:,1) - 0.2;
xx2(:,1) = stokes(:,1)./wgt;
xx2(:,2) = stokes(:,1)./wgt;
yy2(:,1) = -ones(size(stokessl)).*0.38;
yy2(:,2) = yy2(:,1) - 0.2;
plotDataProfile(datpfl,wgt,xlab,ylab,xlims,ylims,xlog,...
                zz,hml,lcolor,[],...
                xx,yy,xx2,yy2);          
saveFigure(l_save_fig, figname);



%% postprocess
% move figures to group dir         
if l_save_fig
    [s, r] = system(['mv ' fn_prefix '_* ' groupName]);
    if s~=0
        error('mv error');
    end
end

%% functions
function h = saveFigure(l_save_fig, figname)
% saveFigure() save and postprocess figures if l_save_fig = 1,
%   do nothing otherwise.

    if (l_save_fig)
        figName = [figname '.fig'];
        saveas(gcf,figName,'fig');
        postProcessFig(figName,5,5);
    end
end
