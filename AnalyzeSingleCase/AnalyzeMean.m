close all;
% clear all;
% This script analyze the mean profiles for a single case.
% It reads in the history files of NCAR LES in netCDF format
% produced by ncarles_his2nc.f90

% data are saved in the directory structure
% /dataDir/cases/fPrefix.tStartStr.tEndStr.nc
%
% Qing Li, 161108

% add path if necessary
if exist('qingAddPath.m','file') == 2
    qingAddPath
end
% add path for ncarleePflData
path(path,'../share');

% DEBUG
l_debug = 0;
% DEBUG

% save figures or not (0 for testing)
if l_debug
    l_save_fig = 0;
else
    l_save_fig = 1;
end

%% Input
%%% a list of cases
if l_debug
    % FOR TEST
    casename = 'R8_BF05WD05WV12_ST01_ens01';
    % FOR TEST
else
    casename = SET_CASENAME;
end
disp(casename);
% tStart = 1;
% tEnd = 30601;
% files, directories
get_dataRootDir;
dataDir = [dataRootDir '/hist'];
fPrefix = 'his.mp.vis';

[s0, r0] = system(['ls ' dataDir '/' casename ]);
if s0==0
    istart_str = r0(12:17);
    iend_str = r0(19:24);
    tStart = str2double(istart_str);
    tEnd = str2double(iend_str);
end

%% Initialization

% creating object ncarlesPflData
f = ncarlesPflData('caseName', casename,...
                   'dataDir',  dataDir,...
                   'fPrefix',  fPrefix,...
                   'tStart',   tStart,...
                   'tEnd',     tEnd);
nnz = f.nnz;
nscl = f.nscl;
time = f.time;

% get variable names
vars = f.vars;
% number of scalers
if isnan(nscl)
    nS = 1;
else
    nS = nscl;
end
% number of variables
nVar = length(vars);
% number of t related variables
nVS = f.nvscl;
% total number of variables
nVarTot = nVar+(nS-1).*nVS;

% set up for multiple scaler variables
% variable name
vname = cell(nVarTot,1);
for i = 1:nVar
    vname{i} = vars{i};
end
% flag for z: 1, z_w; 0, z_u
f_z = zeros(nVarTot,1);
f_z(1:nVar) = f.f_z;


%% Mean Profile
% get prameters
tmpp = f.getParameters;
% averaging window
if abs(tmpp.fcor) > 0
    l_avg_inertial = 1;
else
    l_avg_inertial = 0;
    te = tmpp.tend;
%     te = 88568;
    ts = te-61094;
end
% get profiles
if l_avg_inertial
    tmpv = f.getProfiles('inertial',[0,1]);
    [its, ite] = f.parsingAvgWindow('inertial',[0,1]);
    ts = f.time(its);
    te = f.time(ite);
    % entrainment statistics
    stat_wb = f.getStatEntrainment('inertial',[0,1]);
else
    [its, ite] = f.parsingAvgWindow('second',[ts, te]);
    tmpv = f.getProfiles('second',[ts, te]);
    % entrainment statistics
    stat_wb = f.getStatEntrainment('second',[ts, te]);
end
% read vertical grid
z_u = tmpv.z_u;
z_w = tmpv.z_w;
% get the mixed layer depth
tmp = tmpv.wtle+tmpv.wtsb;
[~,ind] = min(tmp);
hm = -z_w(ind);
% friction velocity
utau = tmpp.utau(end,:);
% boundary layer depth
hb_ts = f.getBLDepth('max_NSquared');
hb = -mean(hb_ts(its:ite));
% Coriolis parameter
fcor = tmpp.fcor;
% surface heat flux
wtsfc = tmpp.wtsfc;


%% Plot nondimensionalized profiles
disp('Plotting mean profiles...');
lColor = 'k';
% stokes drift
stokes = tmpv.stokes;
wgt = utau;
figname = 'stokes';
label = '$u_s/u_*$';
zz = z_u;
plotDataProfile(stokes,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

% buoyancy flux
wb = zeros(2,nnz);
wb(1,:) = (tmpv.wtle+tmpv.wtsb).*f.alpha.*f.g;
wb(2,:) = tmpv.wtsb.*f.alpha.*f.g;
wgt = [utau.^3./hb;utau.^3./hb];
figname = 'buoyflux';
label = '$\overline{w''b''}/(u_*^3/h_\mathrm{b})$';
zz = [z_w;z_w];
plotDataProfile(wb,wgt,label,zz,[hb;hb],{'-k';'--k'},[],l_save_fig,figname);

% TKE components
% total
tmp = tmpv.engsbz;
tmp(1:end-1) = 0.5.*(tmpv.engsbz(1:end-1)+tmpv.engsbz(2:end));
tke = zeros(2,nnz);
tke(1,:) = tmpv.englez+tmp;
tke(2,:) = tmp;
wgt = [utau.^2;utau.^2];
figname = 'tke';
label = 'TKE$/u_*^2$';
zz = [z_w;z_w];
plotDataProfile(tke,wgt,label,zz,[hb;hb],{'-k';'--k'},{'Total';'SGS'},l_save_fig,figname);

% ups
ups = tmpv.ups;
wgt = utau.^2;
figname = 'ups';
label = '$\overline{u''^2}/u_*^2$';
zz = z_u;
plotDataProfile(ups,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

% vps
vps = tmpv.vps;
wgt = utau.^2;
figname = 'vps';
label = '$\overline{v''^2}/u_*^2$';
zz = z_u;
plotDataProfile(vps,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

% wps
wps = tmpv.wps;
wgt = utau.^2;
figname = 'wps';
label = '$\overline{w''^2}/u_*^2$';
zz = z_w;
plotDataProfile(wps,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

% velocity covariances
uw = tmpv.uwle+tmpv.uwsb;
wgt = utau.^2;
figname = 'uw';
label = '$\overline{u''w''}/u_*^2$';
zz = z_w;
plotDataProfile(uw,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

vw = tmpv.vwle+tmpv.vwsb;
wgt = utau.^2;
figname = 'vw';
label = '$\overline{v''w''}/u_*^2$';
zz = z_w;
plotDataProfile(vw,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

uv = tmpv.uvle;
wgt = utau.^2;
figname = 'uv';
label = '$\overline{u''v''}/u_*^2$';
zz = z_u;
plotDataProfile(uv,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

% TKE budget
stokesProd = tmpv.t_stokes;
wgt = utau.^3./hb;
figname = 'stokesProd';
label = 'Stokes production (scaled)';
zz = z_w;
plotDataProfile(stokesProd,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

shearProd = tmpv.t_rprod+tmpv.t_sprod;
wgt = utau.^3./hb;
figname = 'shearProd';
label = 'Shear production (scaled)';
zz = z_w;
plotDataProfile(shearProd,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

tkeTrans = tmpv.t_tran;
wgt = utau.^3./hb;
figname = 'tkeTrans';
label = 'TKE transport (scaled)';
zz = z_w;
plotDataProfile(tkeTrans,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

dissip = tmpv.t_dsle;
wgt = utau.^3./hb;
figname = 'dissipation';
label = 'Dissipation (scaled)';
zz = z_w;
plotDataProfile(dissip,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

dat = tmpv.t_stokes+tmpv.t_rprod+tmpv.t_sprod;
wgt = utau.^3./hb;
figname = 'stokesShearProd';
label = 'Stokes + Shear production (scaled)';
zz = z_w;
plotDataProfile(dat,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

% mean u, v
uxym = tmpv.uxym;
wgt = utau;
figname = 'uxym';
label = '$\overline{u}/u_*$';
zz = z_u;
plotDataProfile(uxym,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

vxym = tmpv.vxym;
wgt = utau;
figname = 'vxym';
label = '$\overline{v}/u_*$';
zz = z_u;
plotDataProfile(vxym,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

dat = tmpv.uxym+tmpv.stokes;
wgt = utau;
figname = 'lagrangianUxym';
label = '$(\overline{u}+u_s)/u_*$';
zz = z_u;
plotDataProfile(dat,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

% tps
tps = tmpv.tps;
wgt = (utau.^2/hb./f.alpha./f.g).^2;
figname = 'tps';
label = '$\overline{T''^2}(gh_\mathrm{b}\alpha)^2/u_*^2$';
zz = z_u;
plotDataProfile(tps,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

% N2
N2 = tmpv.txym.*f.alpha.*f.g;
tmp = N2;
N2(1:end-1) = (tmp(1:end-1)-tmp(2:end))./(z_u(1:end-1)-z_u(2:end));
N2(end) = N2(end-1);
wgt = 1;
figname = 'N2';
label = '$N^2$ (s$^{-1}$)';
zz = z_w;
plotDataProfile(N2,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

% w skewness
wSkew = tmpv.wcube./tmpv.wps.^(1.5);
wgt = 1;
figname = 'wSkewness';
label = '$\overline{w''^3}/\overline{w''^2}^{3/2}$';
zz = z_w;
plotDataProfile(wSkew,wgt,label,zz,hb,{'-k'},[],l_save_fig,figname);

%% plot hodograph
newFigure(l_save_fig);
% from surface to hb
[~,ind_hb] = min(abs(z_w+hb));
datx = tmpv.uxym(1:ind_hb);
daty = tmpv.vxym(1:ind_hb);
wgt = utau;
p = plot(datx./wgt,daty./wgt,'-k');
p.LineWidth = 1.5;
xlabel('$\overline{u}/u_*$','Interpreter','latex');
ylabel('$\overline{v}/u_*$','Interpreter','latex');
% reference line
xx = xlim;
yy = ylim;
line('XData', [xx(1) xx(2)], 'YData', [0 0],'LineStyle','--',...
    'LineWidth', 1, 'Color', 'k');
if (xx(1)<0 && xx(2)>0)
    line('XData', [0 0], 'YData', [yy(1) yy(2)], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color', 'k');
end

if (l_save_fig)
    figName = 'hodograph.fig';
    saveas(gcf,figName,'fig');
    postProcessFig(figName,5,5);
end

%% Plot barycentric maps
disp('Plotting anisotropic barycentric maps...');
cb_label = '$z/h_\mathrm{b}$';
% get the barycentric coordinates
if l_avg_inertial
    [c,z] = f.getAnisotropicBarycentricCoord('inertial',[0,1]);
else
    [c,z] = f.getAnisotropicBarycentricCoord('second',[ts,te]);
end
% plot the anisotropic barycentric map
newFigure(l_save_fig);
plotAnisotropicBarycentricMap(c(1:ind_hb,:),z(1:ind_hb)./hb,0,cb_label);
caxis([-1.0 0]);
if (l_save_fig)
    figname = 'BarycentricMap';
    print('-depsc2', [figname '.eps']);
end
% plot the directions of eigen vectors
if l_avg_inertial
    [a,~] = f.getAnisotropicTensor('inertial',[0,1]);
else
    [a,~] = f.getAnisotropicTensor('second',[ts,te]);
end
cmin = zeros([ind_hb,3]);
cmax = zeros([ind_hb,3]);
lambda = zeros([ind_hb,3]);
for i=1:ind_hb
    % mean profile
    ai = squeeze(a(i,:,:));
    [lambda(i,:),cmax(i,:),cmin(i,:)] = eigMaxMin3(ai);
end
% principle axis corresponding to the biggest eigen vector
newFigure(l_save_fig);
plotEigenVectorDirectionMaxMin(lambda,cmax,cmin,z(1:ind_hb)./hb,1,cb_label);
caxis([-1.0 0]);
if (l_save_fig)
    figname = 'EigenVDir';
    print('-depsc2', [figname '.eps']);
end

%% anisotropic profile
% figure 1: ups, vps and wps normalized by tke
fig = newFigure(l_save_fig);
fig.Units = 'inches';
fig.Position = [1 1 10 6];
fsize = 12;
fsize2 = 14;
tke2 = 2.*tmpv.englez;
upsw = tmpv.ups;
upsw(1:end-1) = 0.5.*(tmpv.ups(1:end-1)+tmpv.ups(2:end));
vpsw = tmpv.vps;
vpsw(1:end-1) = 0.5.*(tmpv.vps(1:end-1)+tmpv.vps(2:end));
uvlew = tmpv.uvle;
uvlew(1:end-1) = 0.5.*(tmpv.uvle(1:end-1)+tmpv.uwle(2:end));
subplot('Position',[0.1 0.1 0.375 0.85]);
p1 = plot(upsw./tke2-1/3,z_w./hb,'-k','LineWidth',1.5);
ax = gca;
ax.FontSize = fsize;
hold on;
p2 = plot(vpsw./tke2-1/3,z_w./hb,'--k','LineWidth',1.5);
p3 = plot(tmpv.wps./tke2-1/3,z_w./hb,'-.k','LineWidth',1.5);
xlim([-0.5 0.5]);
ylim([-1.0 0]);
xlabel('$a_{ij}$','Interpreter','latex','FontSize',fsize2);
ylabel('$z/h_\mathrm{b}$','Interpreter','latex','FontSize',fsize2);
lg = legend([p1, p2, p3],'$11$',...
                         '$22$',...
                         '$33$');
lg.Interpreter = 'latex';
lg.FontSize = fsize;
lg.Location = 'southeast';

subplot('Position',[0.575 0.1 0.375 0.85]);
p4 = plot(uvlew./tke2,z_w./hb,'-k','LineWidth',1.5);
ax = gca;
ax.FontSize = fsize;
hold on;
p5 = plot(tmpv.uwle./tke2,z_w./hb,'--k','LineWidth',1.5);
p6 = plot(tmpv.vwle./tke2,z_w./hb,'-.k','LineWidth',1.5);
ylim([-1.0 0]);
xlim([-0.25 0.25]);
xlabel('$a_{ij}$','Interpreter','latex','FontSize',fsize2);
ylabel('$z/h_\mathrm{b}$','Interpreter','latex','FontSize',fsize2);
lg = legend([p4, p5, p6],'$12$',...
                         '$13$',...
                         '$23$');
lg.Interpreter = 'latex';
lg.FontSize = fsize;
lg.Location = 'southeast';

if (l_save_fig)
figname = 'AnisoPfl';
print('-depsc2', [figname '.eps']);
end

%% Save data
matfile = 'MeanProfile.mat';
save(matfile,'hm','hb','utau','stat_wb','f','fcor','wtsfc',...
             'wb','uw','vw','uv','z_u','z_w',...
             'wps','ups','vps','tps','uxym','vxym','N2','tke', 'stokes',...
             'wSkew','stokesProd','shearProd','tkeTrans','dissip');


%% Post processing
% save figures and make a web page
if (l_save_fig)
    disp('Making a web page...');
    % add path
    path1 = getenv('PATH');
    path1 = [path1 ':/usr/local/bin:/opt/local/bin'];
    setenv('PATH', path1);
    arcDir = outRootDir;
    fHTML = 'ncarles_mean_index.html';
    ts_str = sprintf('%g',ts);
    te_str = sprintf('%g',te);
    utau_str = sprintf('%4.2e',utau);
    wtsfc_str = sprintf('%4.2e',tmpp.wtsfc.*f.rho.*f.cp);
    fcor_str = sprintf('%4.2e',tmpp.fcor);
    system(['./doweb.sh ' casename ' ' fHTML ' ' arcDir ' ' matfile...
        ' ' ts_str ' ' te_str ' ' utau_str ' ' wtsfc_str ' ' fcor_str]);
    % Done
    disp('Done!');
end

