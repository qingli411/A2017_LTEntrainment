close all; clear variables;
% load shared script
path(path,'../share');

%% data info
namelist = {'R8_BF05WD05WV00_ST01_ens03',...
            'R8_BF05WD05WV12_ST01_ens03',...
            'R8_BF1hWD00WV00_ST00_ens03',...
            'R8_BF05WD05WV00_fL00_ens03',...
            'R8_BF05WD05WV12_fL00_ens03'};
startlist = {'028800','028800','021600','021600','021600'};
endlist   = {'032400','032400','028800','025200','025200'};
tstart = startlist{1};
tend = endlist{1};
tstart3 = startlist{3};
tend3 = endlist{3};
tstart4 = startlist{4};
tend4 = endlist{4};
tstart5 = startlist{5};
tend5 = endlist{5};
casename = namelist{1};
casename2 = namelist{2};
casename3 = namelist{3};
casename4 = namelist{4};
casename5 = namelist{5};
h0 = 42;    % initial mixed layer depth
izi = [2, 6, 40, 63, 64, 65, 66, 67, 68];   % indices of saved x-y slices
% izi3 = [2, 6, 44, 77, 78, 79, 80, 81, 82];
% izi4 = [2, 6, 41, 68, 69, 70, 71, 72, 73];
idxentr = [5, 7, 7, 5, 7];
l_save_fig = 1;
c_gray = [0.5, 0.5, 0.5];
% filename
filename_xy = ['viz.vis.' tstart '.' tend '.xy.nc'];
filename3_xy = ['viz.vis.' tstart3 '.' tend3 '.xy.nc'];
filename4_xy = ['viz.vis.' tstart4 '.' tend4 '.xy.nc'];
filename5_xy = ['viz.vis.' tstart5 '.' tend5 '.xy.nc'];

% set up directories
get_dataRootDir;    % get dataRootDir and outRootDir
dataDir = [dataRootDir '/viz/'];
outDir = [outRootDir '/slice/' casename '/' tstart '-' tend];
system(['mkdir -p ' outDir]);

infile1 = [dataDir casename '/' filename_xy];
infile2 = [dataDir casename2 '/' filename_xy];
infile3 = [dataDir casename3 '/' filename3_xy];
infile4 = [dataDir casename4 '/' filename4_xy];
infile5 = [dataDir casename5 '/' filename5_xy];

[wpbp1, wpbpm1] = get_wb(infile1, idxentr(1));
[wpbp2, wpbpm2] = get_wb(infile2, idxentr(2));
[wpbp3, wpbpm3] = get_wb(infile3, idxentr(3));
[wpbp4, wpbpm4] = get_wb(infile4, idxentr(4));
[wpbp5, wpbpm5] = get_wb(infile5, idxentr(5));

nt = numel(wpbpm1);

% normalize the data
wpbp1_pct05 = prctile(wpbp1(:),  5);
wpbp1_pct95 = prctile(wpbp1(:), 95);
wpbp1_flt = wpbp1(wpbp1>wpbp1_pct05 & wpbp1<wpbp1_pct95);
wpbp2_pct05 = prctile(wpbp2(:),  5);
wpbp2_pct95 = prctile(wpbp2(:), 95);
wpbp2_flt = wpbp2(wpbp2>wpbp2_pct05 & wpbp2<wpbp2_pct95);
wpbp3_pct05 = prctile(wpbp3(:),  5);
wpbp3_pct95 = prctile(wpbp3(:), 95);
wpbp3_flt = wpbp3(wpbp3>wpbp3_pct05 & wpbp3<wpbp3_pct95);
wpbp4_pct05 = prctile(wpbp4(:),  5);
wpbp4_pct95 = prctile(wpbp4(:), 95);
wpbp4_flt = wpbp4(wpbp4>wpbp4_pct05 & wpbp4<wpbp4_pct95);
wpbp5_pct05 = prctile(wpbp5(:),  5);
wpbp5_pct95 = prctile(wpbp5(:), 95);
wpbp5_flt = wpbp5(wpbp5>wpbp5_pct05 & wpbp5<wpbp5_pct95);

wpbp1_dist = fitdist(wpbp1_flt(:), 'normal');
wpbp2_dist = fitdist(wpbp2_flt(:), 'normal');
wpbp3_dist = fitdist(wpbp3_flt(:), 'normal');
wpbp4_dist = fitdist(wpbp4_flt(:), 'normal');
wpbp5_dist = fitdist(wpbp5_flt(:), 'normal');

% wgt1 = wpbp1_dist.sigma;
% wgt2 = wpbp2_dist.sigma;
wgt1 = std(wpbp1(:));
wgt2 = std(wpbp2(:));
wgt3 = std(wpbp3(:));
wgt4 = std(wpbp4(:));
wgt5 = std(wpbp5(:));

wpbpnm1 = wpbp1./wgt1;
wpbpnm2 = wpbp2./wgt2;
wpbpnm3 = wpbp3./wgt3;
wpbpnm4 = wpbp4./wgt4;
wpbpnm5 = wpbp5./wgt5;

% pdf
[histy1, histedge1] = histcounts(wpbpnm1(:),'Normalization','pdf');
histx1 = 0.5.*(histedge1(1:end-1)+histedge1(2:end));
[histy2, histedge2] = histcounts(wpbpnm2(:),'Normalization','pdf');
histx2 = 0.5.*(histedge2(1:end-1)+histedge2(2:end));
[histy3, histedge3] = histcounts(wpbpnm3(:),'Normalization','pdf');
histx3 = 0.5.*(histedge3(1:end-1)+histedge3(2:end));
[histy4, histedge4] = histcounts(wpbpnm4(:),'Normalization','pdf');
histx4 = 0.5.*(histedge4(1:end-1)+histedge4(2:end));
[histy5, histedge5] = histcounts(wpbpnm5(:),'Normalization','pdf');
histx5 = 0.5.*(histedge5(1:end-1)+histedge5(2:end));

% fit a normal distribution
wpbpnm1_pct05 = prctile(wpbpnm1(:),  5);
wpbpnm1_pct95 = prctile(wpbpnm1(:), 95);
wpbpnm1_flt = wpbpnm1(wpbpnm1>wpbpnm1_pct05 & wpbpnm1<wpbpnm1_pct95);
wpbpnm2_pct05 = prctile(wpbpnm2(:),  5);
wpbpnm2_pct95 = prctile(wpbpnm2(:), 95);
wpbpnm2_flt = wpbpnm2(wpbpnm2>wpbpnm2_pct05 & wpbpnm2<wpbpnm2_pct95);
wpbpnm3_pct05 = prctile(wpbpnm3(:),  5);
wpbpnm3_pct95 = prctile(wpbpnm3(:), 95);
wpbpnm3_flt = wpbpnm3(wpbpnm3>wpbpnm3_pct05 & wpbpnm3<wpbpnm3_pct95);
wpbpnm4_pct05 = prctile(wpbpnm4(:),  5);
wpbpnm4_pct95 = prctile(wpbpnm4(:), 95);
wpbpnm4_flt = wpbpnm4(wpbpnm4>wpbpnm4_pct05 & wpbpnm4<wpbpnm4_pct95);
wpbpnm5_pct05 = prctile(wpbpnm5(:),  5);
wpbpnm5_pct95 = prctile(wpbpnm5(:), 95);
wpbpnm5_flt = wpbpnm5(wpbpnm5>wpbpnm5_pct05 & wpbpnm5<wpbpnm5_pct95);

wpbpnm1_dist = fitdist(wpbpnm1_flt(:), 'normal');
wpbpnm2_dist = fitdist(wpbpnm2_flt(:), 'normal');
wpbpnm3_dist = fitdist(wpbpnm3_flt(:), 'normal');
wpbpnm4_dist = fitdist(wpbpnm4_flt(:), 'normal');
wpbpnm5_dist = fitdist(wpbpnm5_flt(:), 'normal');

hold on;
tmp_dist = wpbpnm1_dist;
tmp_dist.mu = 0;
fity1 = tmp_dist.pdf(histx1);
plot(histx1,fity1,'--','LineWidth',1.5,'Color', c_gray);
fity2 = wpbpnm2_dist.pdf(histx2);
fity3 = wpbpnm3_dist.pdf(histx3);
fity4 = wpbpnm4_dist.pdf(histx4);
fity5 = wpbpnm5_dist.pdf(histx5);
% plot(histx2,fity2,'--r','LineWidth',1.5);

plot(histx3,histy3,'-','Color',c_gray,'LineWidth',1.5);
ylims = [1e-5,1e0];
rl = line([0,0],ylims);
rl.Color = 'k';
plot(histx1,histy1,'-b','LineWidth',1.5);
plot(histx2,histy2,'-r','LineWidth',1.5);
plot(histx4,histy4,'-c','LineWidth',1.5);
plot(histx5,histy5,'-m','LineWidth',1.5);

xlabel('Normalized buoyancy flux');
ylabel('PDF');
% xlim([-10e6, 10e6]);
set(gca, 'YScale', 'log');
ylim(ylims);
% output figure name
figname = [outDir '/wbPDF_wCT.fig'];
if l_save_fig
    saveas(gcf, figname, 'fig');
    postProcessFig(figname, 6, 4);
end

% flux Richardson number
[upwp1, upwpm1] = get_uw(infile1, idxentr(1));
[upwp2, upwpm2] = get_uw(infile2, idxentr(2));
[upwp3, upwpm3] = get_uw(infile3, idxentr(3));
[upwp4, upwpm4] = get_uw(infile4, idxentr(4));
[upwp5, upwpm5] = get_uw(infile5, idxentr(5));

Rf1 = wpbp1./upwp1;
Rf2 = wpbp2./upwp2;
Rf3 = wpbp3./upwp3;
Rf4 = wpbp4./upwp4;
Rf5 = wpbp5./upwp5;

% pdf
[histy1, histedge1] = histcounts(Rf1(:),'Normalization','pdf');
histx1 = 0.5.*(histedge1(1:end-1)+histedge1(2:end));
[histy2, histedge2] = histcounts(Rf2(:),'Normalization','pdf');
histx2 = 0.5.*(histedge2(1:end-1)+histedge2(2:end));
[histy3, histedge3] = histcounts(Rf3(:),'Normalization','pdf');
histx3 = 0.5.*(histedge3(1:end-1)+histedge3(2:end));
[histy4, histedge4] = histcounts(Rf4(:),'Normalization','pdf');
histx4 = 0.5.*(histedge4(1:end-1)+histedge4(2:end));
[histy5, histedge5] = histcounts(Rf5(:),'Normalization','pdf');
histx5 = 0.5.*(histedge5(1:end-1)+histedge5(2:end));

% plot
figure;
hold on;
% ylims = [1e-5,1e0];
% rl = line([0,0],ylims);
% rl.Color = 'k';
plot(histx1,histy1,'-b','LineWidth',1.5);
plot(histx2,histy2,'-r','LineWidth',1.5);
plot(histx4,histy4,'-c','LineWidth',1.5);
plot(histx5,histy5,'-m','LineWidth',1.5);

xlabel('$R_f$','Interpreter','latex');
ylabel('PDF');
% xlim([-10e6, 10e6]);
set(gca, 'YScale', 'log');
% ylim(ylims);

figname = [outDir '/RfPDF.fig'];
% if l_save_fig
%     saveas(gcf, figname, 'fig');
%     postProcessFig(figname, 6, 4);
% end

function [upwp, upwpm] = get_uw(infile, idx0)
    
    idx0p1 = idx0+1;
    tmp = ncread(infile,'w');
    wxy = squeeze(tmp(:,:,idx0,:));
    tmp = ncread(infile,'u');
    uxy0 = squeeze(tmp(:,:,idx0,:));
    uxy0p1 = squeeze(tmp(:,:,idx0p1,:));
    clear tmp;
    uxy = 0.5.*(uxy0+uxy0p1);
    uxym = mean(mean(uxy,1),2);
    wxym = mean(mean(wxy,1),2);
    upwp = (uxy-uxym).*(wxy-wxym);
    upwpm = squeeze(mean(mean(upwp,1),2));
end

function [wpbp, wpbpm] = get_wb(infile, idx0)
    alpha = 2.0000e-04;
    g = 9.81;
    % also read the buoyancy above an below the entrainment
    % layer for the flux limiter
    idx0m1 = idx0-1;
    idx0p1 = idx0+1;
    idx0p2 = idx0+2;
    tmp = double(ncread(infile,'w'));
    wxy0 = squeeze(tmp(:,:,idx0,:));
    tmp = double(ncread(infile,'t').*alpha.*g);
    bxy0 = squeeze(tmp(:,:,idx0,:));
    bxy0m1 = squeeze(tmp(:,:,idx0m1,:));
    bxy0p1 = squeeze(tmp(:,:,idx0p1,:));
    bxy0p2 = squeeze(tmp(:,:,idx0p2,:));
    % calculate the entrainment buoyancy fluc, apply flux limiter
    tmp = max(-wxy0, 0).*(bxy0+fluxLimiter(bxy0p1, bxy0, bxy0m1));
    tmp2 = min(-wxy0, 0).*(bxy0p1+fluxLimiter(bxy0, bxy0p1, bxy0p2));
    wpbp = -tmp-tmp2;
    wpbpm = squeeze(mean(mean(wpbp,1),2)); % mean entrainment flux
end