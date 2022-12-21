close all; clear variables;
% This script analyze the time series of profiles for a single case.
% It reads in the history files of NCAR LES in netCDF format
% produced by ncarles_his2nc.f90

% data are saved in the directory structure
% /dataDir/cases/fPrefix.tStartStr.tEndStr.nc
%
% Qing Li, 161117

% add path if necessary
if exist('qingAddPath.m','file') == 2
    qingAddPath
end
% add path for ncarleePflData
path(path,'../share');

% save figures or not (0 for testing)
l_save_fig = 0;

%% Input
%%% a list of cases
casename = 'R8_BF05WD05WV12_ST01_ens03';

% files, directories
get_dataRootDir;
% dataRootDir = '/Users/qingli/data_local/archive_les';
dataDir = [dataRootDir '/hist'];
fPrefix = 'his.mp.vis';

[s0, r0] = system(['ls ' dataDir '/' casename ]);
if s0==0
    istart_str = r0(12:17);
    iend_str = r0(19:25);
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
ntime = f.ntime;

% % get variable names
% vars = f.vars;
% % number of scalers
% if isnan(nscl)
%     nS = 1;
% else
%     nS = nscl;
% end
% % number of variables
% nVar = length(vars);
% % number of t related variables
% nVS = f.nvscl;
% % total number of variables
% nVarTot = nVar+(nS-1).*nVS;
%
% % set up for multiple scaler variables
% % variable name
% vname = cell(nVarTot,1);
% for i = 1:nVar
%     vname{i} = vars{i};
% end
% % flag for z: 1, z_w; 0, z_u
% f_z = zeros(nVarTot,1);
% f_z(1:nVar) = f.f_z;

% time of interest
[its, ite] = f.parsingAvgWindow('inertial',[0,2]);
% hour
thour = time./3600;
% get dt
dt = time;
dt(2:end-1) = 0.5.*(time(3:end)-time(1:end-2));
dt(1) = dt(2);
dt(end) = dt(end-1);
% z
z_w = f.getVar('z_w');
z_u = f.getVar('z_u');
% mixed layer depth
h = f.getBLDepth('max_Nsquared_smth');
mean_h = f.wgtMean(h(its:ite),dt(its:ite),2);
[v,h_inds] = min(abs(z_w-mean_h));
% for plot
c_gray = [0.5 0.5 0.5];
% colormap('jet');

%% Figure A1: Time series of entrainemnt

% buoyancy flux
tmp = f.getVar('wtle');
wtle = squeeze(tmp(:,1,:));
tmp = f.getVar('wtsb');
wtsb = squeeze(tmp(:,1,:));
utau = f.getVar('utau')';
wbpfl = (wtle+wtsb).*f.alpha.*f.g;
wbpfl2 = stretch_z(wbpfl, h, mean_h, z_w);
meandat = f.wgtMean(wbpfl(:,its:ite),dt(its:ite),2);
meandat2 = f.wgtMean(wbpfl2(:,its:ite),dt(its:ite),2);
[min_wb,hm_inds] = min(meandat2);
max_wb = max(meandat2);
inds1 = 2*hm_inds-h_inds;
inds2 = h_inds;
wb = wbpfl2(hm_inds,its:ite);
wb2 = mean(wbpfl2(inds1:inds2,its:ite),1);
wbp75 = prctile(wb,75);
wbp50 = prctile(wb,50);
wbp25 = prctile(wb,25);
ydepths = [floor((mean_h*1.4)/10)*10, 0];
mean_hm = z_w(hm_inds);

% create figure
lst = 0.08;
wid = 0.63;
ll = [lst, lst, lst, 0.8];
bb = [0.73,0.45 0.10, 0.45];
ww = [wid, wid, wid, 0.17];
hh = [0.22, 0.22, 0.25, 0.5];
cmin = min_wb.*1.1;
cmax = max_wb.*1.1;

fig=figure(1);
fig.Units = 'inches';
fig.Position = [1 1 12 8];
fig.Color = 'white';

% subfigure 1
subplot('position',[ll(1) bb(1) ww(1) hh(1)])
% map = qing_readNCLColormap('amwg256.rgb');
pcolor(thour(its:ite),z_w,wbpfl(:,its:ite)); shading flat;
% contourf(thour(its:ite),z_w,wbpfl(:,its:ite));
caxis([cmin, cmax]);
hold on;
p2 = plot(thour(its:ite),h(its:ite),'-w');
p2.LineWidth = 1.5;
ylabel('z (m)','Interpreter','latex');
ylim(ydepths);
text(0.03, 0.14,'(a)','Units','normalized','BackgroundColor','w',...
    'Interpreter','latex');

% subfigure 2
subplot('position',[ll(2) bb(2) ww(2) hh(2)]);
pcolor(thour(its:ite),z_w,wbpfl2(:,its:ite)); shading flat;
% contourf(thour(its:ite),z_w,wbpfl2(:,its:ite));
caxis([cmin, cmax]);
hold on;
hentr = ones(size(thour)).*mean_hm;
hmn = ones(size(thour)).*mean_h;
p2 = plot(thour(its:ite),hmn(its:ite),'-w');
p2.LineWidth = 1.5;
p3 = plot(thour(its:ite),hentr(its:ite),'--w');
p3.LineWidth = 1.5;
ylabel('z (m)','Interpreter','latex');
ylim(ydepths);
text(0.03, 0.14,'(b)','Units','normalized','BackgroundColor','w',...
    'Interpreter','latex');

% subfigure 3
subplot('position',[ll(3) bb(3) ww(3) hh(3)]);
box on;
hold on;
wbp = ones(size(wb)).*min_wb;
p2 = plot(thour(its:ite),wbp(:),'-','Color',c_gray,'LineWidth',0.8);
wbp = ones(size(wb)).*wbp25;
p3 = plot(thour(its:ite),wbp(:),'--','Color',c_gray,'LineWidth',0.8);
wbp = ones(size(wb)).*wbp50;
p4 = plot(thour(its:ite),wbp(:),'--','Color',c_gray,'LineWidth',0.8);
wbp = ones(size(wb)).*wbp75;
p5 = plot(thour(its:ite),wbp(:),'--','Color',c_gray,'LineWidth',0.8);
p1 = plot(thour(its:ite),wb,'-k');
p0 = plot(thour(its:ite),wb2,'-r');
p1.LineWidth = 1.5;
xlabel('Time (h)','Interpreter','latex');
ylabel('$\overline{w''b''}_\mathrm{e}$ (m$^2$ s$^{-3}$)','Interpreter','latex');
xlim([thour(its), thour(ite)]);
ylims = ylim;
text(0.03, 0.14,'(c)','Units','normalized','BackgroundColor','w',...
    'Interpreter','latex');

% axes('Position', [0.8 bb(3) 0.1 hh(3)], 'Visible', 'off');
% ylim(ylims)
% scatter(28, -1e-8,'o','Clipping','off')

% subfigure 4
subplot('position',[ll(4) bb(4) ww(4) hh(4)])
% newFigure(l_save_fig);
plot(meandat,z_w,'--k','LineWidth',1.5);
hold on;
plot(meandat2,z_w,'-k','LineWidth',1.5);
plot(min_wb, mean_hm,'o','Color',c_gray,'LineWidth',1.5);
neg = wbp50-wbp25;
pos = wbp75-wbp50;
errorbar(wbp50,mean_hm,neg,pos,'horizontal',...
    'Color',c_gray,'LineWidth',1.5);

ylimwb = [floor((mean_h*1.1)/5)*5, ceil((mean_hm*0.8)/5)*5];
ylim(ylimwb);
xlimwb = [cmin*1.1, -cmin*0.3];
xlim(xlimwb);
ylabel('z (m)','Interpreter','latex');
xlabel('$\overline{w''b''}_\mathrm{e}$','Interpreter','latex');
text(0.1, 0.1,'(d)','Units','normalized','BackgroundColor','w',...
    'Interpreter','latex');

set(findall(gcf,'-property','FontSize'),'FontSize',14)

%
% axes('Position', [ll(3) bb(3) 0.8 hh(3)], 'Visible', 'off');
% colorbar;
% caxis([cmin, cmax]);

if (l_save_fig)
    figName = 'ts_wb';
    saveas(gcf,[figName '.fig'],'fig');
%     print('-depsc2', [figName '.eps']);
    saveas(gcf,[figName '.eps'],'epsc');
%     export_fig 'ts_wb2.eps' -eps;

end

%% Figure A2: Time series of current
stokes = f.getVar('stokes');
uxym = f.getVar('uxym');
vxym = f.getVar('vxym');
dudz = f.getVar('dudz');
dvdz = f.getVar('dvdz');
M2 = dudz.^2+dvdz.^2;
M2_SZ = stretch_z(M2, h, mean_h, z_w);
M2hm = mean(M2_SZ(hm_inds:h_inds,its:ite),1);

% create figure
lst = 0.08;
wid = 0.63;
ll = [lst, lst, lst, 0.8];
bb = [0.73,0.45 0.10, 0.45];
ww = [wid, wid, wid, 0.17];
hh = [0.22, 0.22, 0.25, 0.5];
tmp = uxym(:,its:ite);
max_uxym = max(tmp(:));
min_uxym = min(tmp(:));
tmp = vxym(:,its:ite);
max_vxym = max(tmp(:));
min_vxym = min(tmp(:));
tmp = M2(floor(hm_inds/2):end,its:ite);
max_M2 = max(tmp(:));
min_M2 = min(tmp(:));

cmax = max([abs(max_uxym), abs(max_vxym), abs(min_uxym), abs(min_vxym)]);
cmin = -cmax;

cmax2 = max_M2;
cmin2 = min_M2;

fig=figure(2);
fig.Units = 'inches';
fig.Position = [1 1 12 8];
fig.Color = 'white';

% subfigure 1
subplot('position',[ll(1) bb(1) ww(1) hh(1)]);
pcolor(thour(its:ite),z_u,uxym(:,its:ite)+stokes(:,its:ite)); shading flat;
caxis([cmin, cmax]);
ylabel('z (m)','Interpreter','latex');
ylim(ydepths);
hold on;
p2 = plot(thour(its:ite),h(its:ite),'-w');
p2.LineWidth = 1.5;
text(0.03, 0.14,'(a)','Units','normalized','BackgroundColor','w',...
    'Interpreter','latex');

% subfigure 2
subplot('position',[ll(2) bb(2) ww(2) hh(2)]);
pcolor(thour(its:ite),z_u,vxym(:,its:ite)); shading flat;
caxis([cmin, cmax]);
ylabel('z (m)','Interpreter','latex');
ylim(ydepths);
hold on;
p2 = plot(thour(its:ite),h(its:ite),'-w');
p2.LineWidth = 1.5;
text(0.03, 0.14,'(b)','Units','normalized','BackgroundColor','w',...
    'Interpreter','latex');

% subfigure 3
subplot('position',[ll(3) bb(3) ww(3) hh(3)]);
pcolor(thour(its:ite),z_u,M2(:,its:ite)); shading flat;
caxis([cmin2, cmax2]);
ylabel('z (m)','Interpreter','latex');
ylim(ydepths);
hold on;
p2 = plot(thour(its:ite),h(its:ite),'-w');
p2.LineWidth = 1.5;
text(0.03, 0.14,'(c)','Units','normalized','BackgroundColor','w',...
    'Interpreter','latex');

set(findall(gcf,'-property','FontSize'),'FontSize',14)

% Figure A3: anisotropic barycentric plot
ups = f.getVar('ups');
vps = f.getVar('vps');
wps = f.getVar('wps');
uwle = f.getVar('uwle');
vwle = f.getVar('vwle');
uvle = f.getVar('uvle');

ups2 = stretch_z(ups, h, mean_h, z_w);
vps2 = stretch_z(vps, h, mean_h, z_w);
wps2 = stretch_z(wps, h, mean_h, z_w);
uvle2 = stretch_z(uvle, h, mean_h, z_w);
uwle2 = stretch_z(uwle, h, mean_h, z_w);
vwle2 = stretch_z(vwle, h, mean_h, z_w);

nt = ite-its+1;
c = zeros(nt, 3);
cmax = zeros(nt,3);
cmin = zeros(nt,3);
lambda = zeros(nt,3);
j = 1;
for i = its:ite
    a = f.anisotropy_tensor(ups2(hm_inds,i), vps2(hm_inds,i), wps2(hm_inds,i),...
                    uvle2(hm_inds,i),uwle2(hm_inds,i),vwle2(hm_inds,i));
    c(j,:) = f.barycentric_coord(a);
    [E,D] = eig(a);
    F = [D(1,1) D(2,2) D(3,3);E];
    FF = sortrows(F',-1);
    cmax(j,:) = FF(1,2:4)';
    cmin(j,:) = FF(3,2:4)';
    lambda(j,:) = FF(1:3,1);
    j = j+1;
end
figure(3)
plotAnisotropicBarycentricMap(c, M2hm, 1);

% principle axis corresponding to the biggest eigen vector
figure(4);
plotEigenVectorDirection(cmax.^2,cmin.^2,lambda,M2hm,0);

figure(5);
pcolor(thour(its:ite), z_w, M2_SZ(:,its:ite)); shading flat;
caxis([min_M2, max_M2]);
hold on;
p3 = plot(thour(its:ite),hentr(its:ite),'--w');
p2 = plot(thour(its:ite),hmn(its:ite),'-w');
ylim(ydepths);




%% functions
function out = stretch_z(data, h, mean_h, zz)

    ndims = size(data);
    nnz = ndims(1);
    ntime = ndims(2);
    zn_w = zz./mean_h;
    tmp = zeros(nnz+1,ntime);
    tmp(2:end,:) = data;
    zz_in = zeros(1,nnz+1);
    zz_in(2:end) = zz;
    zz_out = zn_w;
    out = data;
    for it = 1:ntime
        out(:,it) = interp1(zz_in./h(it),...
        tmp(:,it),zz_out,'linear',NaN);
    end
end
