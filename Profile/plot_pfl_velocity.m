close all; clear variables;

l_save_fig = 1;

%% set parameters
get_dataRootDir;    % get dataRootDir and outRootDir

namelist = {'R8_BF05WD05WV00_ST01_ens03',...
            'R8_BF05WD05WV12_ST01_ens03',...
            'R8_BF05WD05WV00_fL00_ens03',...
            'R8_BF05WD05WV12_fL00_ens03'};
tStartList =[1, 1, 1, 1];
tEndList = [32401, 32401, 25201, 25201];
ize = [64, 66, 69, 71];
dt = 61094;
nnz = 256;
outdir = [outRootDir '/profiles/' namelist{1}];
c_gray = [0.5 0.5 0.5];
nc = numel(namelist);
u = zeros([nc,nnz]);
v = zeros([nc,nnz]);
wb = zeros([nc,nnz]);
hb = zeros([nc,1]);
izhhb = zeros([nc,1]);
wps = zeros([nc,nnz]);
izmwps = zeros([nc,1]);
wcube = zeros([nc,nnz]);
wbtmp = zeros([nc,nnz]);
zztmp = zeros([nc,nnz]);
for ic = 1:nc
    casename = namelist{ic};
    tStart = tStartList(ic);
    tEnd = tEndList(ic);
    f = initncarPflData(dataRootDir, casename, tStart, tEnd);
    pfl = f.getProfiles('second', [f.time(end)-dt, f.time(end)]);
    u(ic,:) = pfl.uxym+pfl.stokes;
    v(ic,:) = pfl.vxym;
    wb(ic,:) = (pfl.wtle+pfl.wtsb).*f.g.*f.alpha;
%     wb(ic,:) = (pfl.wtle).*f.g.*f.alpha;
    wps(ic,:) = pfl.wps;
    [~,izmwps(ic)] = max(wps(ic,:));
    wcube(ic,:) = pfl.wcube;
    hb(ic) = pfl.mean_hb;
    if ic == 1
        z_u = pfl.z_u;
        z_w = pfl.z_w;
        tmp = f.getVar('utau');
        utau = tmp(end);
    end
    if ic == 2
        ustokes = pfl.stokes;
    end
    [~,izhhb(ic)] = min(abs(hb(ic)./2+z_w));
end

% plot u
figure;
hold on;
linecolor = {'b', 'r', 'c', 'm'};
ylims = [-1.2, 0];
xlims = [-6, 24];
rl = line(xlims, [-1,-1]);
rl.Color = 'k';
rl = line([0,0], ylims);
rl.Color = 'k';
nanarray = ustokes*nan;
pstokes = plot(ustokes./utau, z_u./hb(2), '-.', 'Color', 'k', 'LineWidth', 1.5);
pu = plot(nanarray, z_u, '-', 'Color', 'k', 'LineWidth', 1.5);
pv = plot(nanarray, z_u, '--', 'Color', 'k', 'LineWidth', 1.5);
for ic = 1:nc
    plot(u(ic,:)./utau, z_u./hb(ic), '-',...
        'Color', linecolor{ic}, 'LineWidth', 1.5);
    plot(v(ic,:)./utau, z_u./hb(ic), '--',...
        'Color', linecolor{ic}, 'LineWidth', 1.5)
end
ylim(ylims);
xlim(xlims);
% xlabel('$\overline{u}_{i}^L/u^*$', 'Interpreter', 'latex');
% ylabel('$x_3/h_b$', 'Interpreter', 'latex');
% figname = [outdir '/FiguL.fig'];
xlabel('$\overline{u}^L/u^*$', 'Interpreter', 'latex');
ylabel('$z/h_b$', 'Interpreter', 'latex');
lg = legend([pu,pv,pstokes],'$\overline{u}$','$\overline{v}$','$u^S$');
lg.Location = 'SouthEast';
lg.Interpreter = 'latex';
figname = [outdir '/FiguL_xyz_s2.fig'];
if l_save_fig
    saveas(gcf, figname, 'fig');
    postProcessFig(figname, 4, 4);
end

% plot v
figure;
hold on;
ylims = [-1.2, 0];
xlims = [-5, 1];
rl = line(xlims, [-1,-1]);
rl.Color = 'k';
rl = line([0,0], ylims);
rl.Color = 'k';
for ic = 1:nc
    plot(v(ic,:)./utau, z_u./hb(ic), '-',...
        'Color', linecolor{ic}, 'LineWidth', 1.5)
end
ylim(ylims);
xlim(xlims);
% xlabel('$\overline{u}_2^L/u^*$', 'Interpreter', 'latex');
% ylabel('$x_3/h_b$', 'Interpreter', 'latex');
% figname = [outdir '/FigvL.fig'];
xlabel('$\overline{v}^L/u^*$', 'Interpreter', 'latex');
ylabel('$z/h_b$', 'Interpreter', 'latex');
figname = [outdir '/FigvL_xyz.fig'];
if l_save_fig
    saveas(gcf, figname, 'fig');
    postProcessFig(figname, 4, 4);
end

% plot wb
figure;
hold on;
ylims = [-1.2, 0];
xlims = [-2, 0.5];
rl = line(xlims, [-1,-1]);
rl.Color = 'k';
rl = line([0,0], ylims);
rl.Color = 'k';
for ic = 1:nc
    plot(wb(ic,:).*hb(ic)./utau.^3, z_w./hb(ic), '-',... 
        'Color', linecolor{ic}, 'LineWidth', 1.5);
    wbtmp(ic,:) = wb(ic,:).*hb(ic)./utau.^3;
    zztmp(ic,:) = z_w./hb(ic);
    plot(wb(ic,ize(ic)).*hb(ic)./utau.^3, z_w(ize(ic))./hb(ic), '+',...
        'Color', linecolor{ic}, 'LineWidth', 1.5);
end
tmp = interp1(zztmp(2,:),wbtmp(2,:),zztmp(1,:));
p1 = plot(tmp-wbtmp(1,:), zztmp(1,:), '--',...
    'Color', 'k', 'LineWidth', 1.5);
tmp = interp1(zztmp(4,:),wbtmp(4,:),zztmp(3,:));
p2 = plot(tmp-wbtmp(3,:), zztmp(3,:), '--',...
    'Color', c_gray, 'LineWidth', 1.5);
% lg = legend([p1,p2],'LT-R $-$ ST-R', 'LT-NR $-$ ST-NR');
% lg.Location = 'NorthWest';
% lg.Interpreter = 'latex';
ylim(ylims);
xlim(xlims);
xlabel('$\overline{u_3''b''} h_b/{u^*}^3$', 'Interpreter', 'latex');
ylabel('$x_3/h_b$', 'Interpreter', 'latex');
figname = [outdir '/Figwb.fig'];
% xlabel('$\overline{w''b''} h_b/{u^*}^3$', 'Interpreter', 'latex');
% ylabel('$z/h_b$', 'Interpreter', 'latex');
% figname = [outdir '/Figwb_xyz.fig'];
if l_save_fig
    saveas(gcf, figname, 'fig');
    postProcessFig(figname, 4, 4);
end

% plot wps
figure;
hold on;
ylims = [-1.2, 0];
xlims = [0, 4];
rl = line(xlims, [-1,-1]);
rl.Color = 'k';
% rl = line([0,0], ylims);
% rl.Color = 'k';
p = gobjects([nc,1]);
for ic = 1:nc
    p(ic) = plot(wps(ic,:)./utau.^2, z_w./hb(ic), '-',...
        'Color', linecolor{ic}, 'LineWidth', 1.5);
    plot(wps(ic,2)./utau.^2, z_w(2)./hb(ic), '+',...
        'Color', linecolor{ic}, 'LineWidth', 1.5);
    plot(wps(ic,izhhb(ic))./utau.^2, z_w(izhhb(ic))./hb(ic), '+',...
        'Color', linecolor{ic}, 'LineWidth', 1.5);
    plot(wps(ic,izmwps(ic))./utau.^2, z_w(izmwps(ic))./hb(ic), '+',...
        'Color', linecolor{ic}, 'LineWidth', 1.5);
    plot(wps(ic,ize(ic))./utau.^2, z_w(ize(ic))./hb(ic), '+',...
        'Color', linecolor{ic}, 'LineWidth', 1.5);
end
ylim(ylims);
xlim(xlims);
% xlabel('$\overline{{u_3''}^2} /{u^*}^2$', 'Interpreter', 'latex');
% ylabel('$x_3/h_b$', 'Interpreter', 'latex');
% figname = [outdir '/Figwps.fig'];
xlabel('$\overline{{w''}^2} /{u^*}^2$', 'Interpreter', 'latex');
ylabel('$z/h_b$', 'Interpreter', 'latex');
figname = [outdir '/Figwps_xyz.fig'];
lg = legend(p,'ST-R','LT-R','ST-NR','LT-NR');
lg.Location = 'SouthEast';
lg.Interpreter = 'latex';
if l_save_fig
    saveas(gcf, figname, 'fig');
    postProcessFig(figname, 4, 4);
end

% plot skewness
figure;
hold on;
ylims = [-1.2, 0];
xlims = [-1.2, 0.3];
rl = line(xlims, [-1,-1]);
rl.Color = 'k';
rl = line([0,0], ylims);
rl.Color = 'k';
for ic = 1:nc
    plot(wcube(ic,:)./wps(ic,:).^(3/2), z_w./hb(ic), '-',... 
        'Color', linecolor{ic}, 'LineWidth', 1.5)
    plot(wcube(ic,2)./wps(ic,2).^(3/2), z_w(2)./hb(ic), '+',...
        'Color', linecolor{ic}, 'LineWidth', 1.5);
    plot(wcube(ic,izhhb(ic))./wps(ic,izhhb(ic)).^(3/2), z_w(izhhb(ic))./hb(ic), '+',...
        'Color', linecolor{ic}, 'LineWidth', 1.5);
    plot(wcube(ic,izmwps(ic))./wps(ic,izmwps(ic)).^(3/2), z_w(izmwps(ic))./hb(ic), '+',...
        'Color', linecolor{ic}, 'LineWidth', 1.5);
    plot(wcube(ic,ize(ic))./wps(ic,ize(ic)).^(3/2), z_w(ize(ic))./hb(ic), '+',...
        'Color', linecolor{ic}, 'LineWidth', 1.5);
end
ylim(ylims);
xlim(xlims);
% xlabel('$\overline{{u_3''}^3}/(\overline{{u_3''}^2})^{3/2}$', 'Interpreter', 'latex');
% ylabel('$x_3/h_b$', 'Interpreter', 'latex');
% figname = [outdir '/FigSkewness.fig'];
xlabel('$\overline{{w''}^3}/(\overline{{w''}^2})^{3/2}$', 'Interpreter', 'latex');
ylabel('$z/h_b$', 'Interpreter', 'latex');
figname = [outdir '/FigSkewness_xyz.fig'];
if l_save_fig
    saveas(gcf, figname, 'fig');
    postProcessFig(figname, 4, 4);
end

%% functions
function f = initncarPflData(dataRootDir, casename, tStart, tEnd)

    % set up profile data parameters
    dataDir = [dataRootDir '/hist/'];
    fPrefix = 'his.mp.vis';
    % creating object ncarlesPflData
    f = ncarlesPflData('caseName', casename,...
                       'dataDir',  dataDir,...
                       'fPrefix',  fPrefix,...
                       'tStart',   tStart,...
                       'tEnd',     tEnd);
    disp(['Reading profile data ', f.filename]);
end

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