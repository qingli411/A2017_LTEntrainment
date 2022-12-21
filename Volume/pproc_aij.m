close all; clear variables;
casename = 'R8_BF05WD05WV00_ST01_ens01';
casename2 = 'R8_BF05WD05WV12_ST01_ens01';

path(path,'../share');
path(path,'../share/cbrewer')
get_dataRootDir;    % get dataRootDir
dataDir = [dataRootDir '/rest/'];
srcname = [dataDir casename '/les.F' ];
nnz = getValue(srcname, 'nnz');
Lz  = getValue(srcname, 'zl');   % z-domain size in meter
dz = Lz/nnz;
zw  = dz:dz:Lz;     % z for w
zu  = dz/2:dz:Lz;   % z for u, v, t

matfile = [casename '/data.mat'];
matfile2 = [casename2 '/data.mat'];

dd = load(matfile);
upsxym = dd.upsxym;
vpsxym = dd.vpsxym;
wpsuxym = dd.wpsuxym;
upvpxym = dd.upvpxym;
upwpxym = dd.upwpxym;
vpwpxym = dd.vpwpxym;
% zu = dd.zu;
% zw = dd.zw;
ind_hb = dd.ind_hb;
clear dd;
dd = load(matfile2);
upsxym2 = dd.upsxym;
vpsxym2 = dd.vpsxym;
wpsuxym2 = dd.wpsuxym;
upvpxym2 = dd.upvpxym;
upwpxym2 = dd.upwpxym;
vpwpxym2 = dd.vpwpxym;
% zu2 = dd.zu;
% zw2 = dd.zw;
ind_hb2 = dd.ind_hb;
clear dd;

% figure 2C: components of a_{ij}
fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 10 6];
fsize = 12;
fsize2 = 14;
tke2 = upsxym+vpsxym+wpsuxym;
a11 = upsxym./tke2-1/3;
a22 = vpsxym./tke2-1/3;
a33 = wpsuxym./tke2-1/3;
a12 = upvpxym./tke2;
a13 = upwpxym./tke2;
a23 = vpwpxym./tke2;
azz = -zu./zw(ind_hb);
tke2 = upsxym2+vpsxym2+wpsuxym2;
b11 = upsxym2./tke2-1/3;
b22 = vpsxym2./tke2-1/3;
b33 = wpsuxym2./tke2-1/3;
b12 = upvpxym2./tke2;
b13 = upwpxym2./tke2;
b23 = vpwpxym2./tke2;
bzz = -zu./zw(ind_hb2);
subplot('Position',[0.1 0.1 0.375 0.85]);
p1 = plot(a11(1:ind_hb),azz(1:ind_hb),'-k','LineWidth',1.5);
ax = gca;
ax.FontSize = fsize;
hold on;
p2 = plot(a22(1:ind_hb),azz(1:ind_hb),'--k','LineWidth',1.5);
p3 = plot(a33(1:ind_hb),azz(1:ind_hb),'-.k','LineWidth',1.5);
p1b = plot(b11(1:ind_hb2),bzz(1:ind_hb2),'-r','LineWidth',1.5);
p2b = plot(b22(1:ind_hb2),bzz(1:ind_hb2),'--r','LineWidth',1.5);
p3b = plot(b33(1:ind_hb2),bzz(1:ind_hb2),'-.r','LineWidth',1.5);
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
p4 = plot(a12(1:ind_hb),azz(1:ind_hb),'-k','LineWidth',1.5);
ax = gca;
ax.FontSize = fsize;
hold on;
p5 = plot(a13(1:ind_hb),azz(1:ind_hb),'--k','LineWidth',1.5);
p6 = plot(a23(1:ind_hb),azz(1:ind_hb),'-.k','LineWidth',1.5);
p4b = plot(b12(1:ind_hb2),bzz(1:ind_hb2),'-r','LineWidth',1.5);
p5b = plot(b13(1:ind_hb2),bzz(1:ind_hb2),'--r','LineWidth',1.5);
p6b = plot(b23(1:ind_hb2),bzz(1:ind_hb2),'-.r','LineWidth',1.5);
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