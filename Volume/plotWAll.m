close all; clear variables;
% load shared scripts
path(path,'../share');

%% flags
l_save_fig = 1;
c_gray = [0.5, 0.5, 0.5];

%% Parameters
get_dataRootDir;
dataDir = [dataRootDir '/rest/'];
namelist = {'R8_BF05WD05WV00_ST01_ens03',...
            'R8_BF05WD05WV12_ST01_ens03',...
            'R8_BF05WD05WV00_fL00_ens03',...
            'R8_BF05WD05WV12_fL00_ens03'};
restlist   = {'032401','032401','025201','025201'};
izilist = [64, 66, 69, 71];
% casename = 'R8_BF05WD05WV00_fL00_ens03';
% [s0, r0] = system(['ls ' dataDir casename ]);
% if ~s0
%     tttttt = r0(1:6);
% end
% tttttt = '025201';
ic = 4;
casename = namelist{ic};
tttttt = restlist{ic};
izi = izilist(ic);

% u file name
pattern = 'u\.[a-z]{2}\.[a-z]{3}[0-9]{3}';
dir = [dataDir casename '/' tttttt '/'];
[s0, r0] = system(['ls ' dir]);
if ~s0
    r1 = regexp(r0, pattern, 'match');
    fileIn = r1{1};
end
% output directory
outDir = [outRootDir '/volume/' casename '/' tttttt];
system(['mkdir -p ' outDir]);

% get grid
srcname = [dataDir casename '/les.F' ];
nnz = getValue(srcname, 'nnz');
nnx = getValue(srcname, 'nnx');
nny = getValue(srcname, 'nny');
Lx  = getValue(srcname, 'xl');   % x-domain size in meter
Ly  = getValue(srcname, 'yl');   % y-domain size in meter
Lz  = getValue(srcname, 'zl');   % z-domain size in meter
nscl = getValue(srcname, 'nscl'); % number of scalers
nvar = 4+nscl; % number of variables
dx = Lx/nnx;
x = dx/2:dx:Lx;     % center of the grid, x
dy = Ly/nny;
y = dy/2:dy:Ly;     % center of the grid, y
dz = Lz/nnz;
zw  = dz:dz:Lz;     % z for w
zu  = dz/2:dz:Lz;   % z for u, v, t
batag = 9.81/5000;

% read file
fid = fopen([dir fileIn]);
tmp = fread(fid, nvar*nnx*nny*nnz, 'double', 'l');
fclose(fid);
clear var;
var = reshape(tmp, nvar, nnx, nny, nnz);
clear tmp;

%% spectrum
% u, v, w at maximum wps depth
% get maximum wps depth
w = squeeze(var(3,:,:,:));
wps = zeros(1,nnz);
for i=1:nnz
    tmp = w(:,:,i);
    tmp2 = (tmp-mean(tmp(:))).^2;
    wps(i) = mean(tmp2(:));
end
[~,mi] = max(wps(:));

depth1 = zw(mi);
ind1 = mi;
fprintf('Depth: %g\n', depth1);
p1 = plotW(x, y, w(:,:,ind1)');

% u, v, w at hb/2
% get boundary layer depth
N2 = zeros(1,nnz);
txym = zeros(1,nnz);
for i=1:nnz
    tmp = squeeze(var(4,:,:,i));
    txym(i) = mean(tmp(:));
end
N2(1:end-1) = (txym(2:end)-txym(1:end-1))./dz;
N2(end) = N2(end-1);
[~,indhb] = max(N2(:));
hb = -zw(indhb);

mi = round(indhb/2);
depth2 = zw(mi);
ind2 = mi;
fprintf('Depth: %g\n', depth2);
p2 = plotW(x, y, w(:,:,ind2)');

% u, v, w near hb: z = 0.9*hb
% mi = round(indhb*0.9);
% depth3 = zw(mi);
% ind3 = mi;
depth3 = zw(izi);
ind3 = izi;
fprintf('Depth: %g\n', depth3);
p3 = plotW(x, y, w(:,:,ind3)');

% u, v, w near surface: z = 1.28 m
mi = 2;
depth4 = zw(mi);
ind4 = mi;
fprintf('Depth: %g\n', depth4);
p4 = plotW(x, y, w(:,:,ind4)');

% all figures together
fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 5 7];
hold on;
[xx,yy,zz] = meshgrid(x./Lx,y./Ly,zw(1:indhb)./hb);
zs = [zw(ind1), zw(ind2), zw(ind3)];
vol = zeros(nnx,nny,indhb);
for i=1:indhb
    tmp = w(:,:,i)';
    tmpstd = std(tmp(:));
    vol(:,:,i) = tmp./tmpstd;
end
vol_std = std(vol(:));
vol = vol./vol_std;
indsx = nnx/2+1:nnx;
indsy = 1:nny;

ps1 = slice(xx,yy,zz,vol,xx(indsx,indsy,ind4),yy(indsx,indsy,ind4),zz(indsx,indsy,ind4)); shading flat;
ps = slice(xx,yy,zz,vol,[],[],zs./hb); shading flat;
colorbar;
vmax = 3;
colormap(db2dr(-vmax,vmax));
xlabel('$x_1/L_1$', 'Interpreter', 'latex');
ylabel('$x_2/L_2$', 'Interpreter', 'latex');
zlabel('$x_3/h_b$',  'Interpreter', 'latex')
view(-25, 35);
xlim([0,1]);
ylim([0,1]);
zlim([-1,0]); 
daspect([1,1,0.6]);
colorbar('off');
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

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

if l_save_fig
    figname = [outDir '/w_all.fig'];
    saveas(gcf, figname, 'fig');
    [dir, name, ~] = fileparts(figname);
%     print('-depsc2',[dir '/' name]);
    print('-dpng', [dir '/' name],'-r600');
end



%% functions
function f = plotW(x, y, dat)
    figure;
    f = pcolor(x, y, dat);
    shading flat;
    colorbar;
    vmax = std(dat(:))*5;
    colormap(db2dr(-vmax,vmax));
    daspect([1, 1, 1]);
end

