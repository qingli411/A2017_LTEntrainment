% This script is used to test the calculation of root mean square 
% veritcal velocity and buoyancy flux profiles from 3d field data.
% The calculated profles are compared with the saved profile data.
%
% - Note that a flux limiter is used when calculating the buoyancy flux
% in NCAR LES, the save flux limiter is used here.
% - Miss matching the w level and t level causes significant error when
% calculating the buoyancy flux.

close all; clear variables;

path(path,'../share');
get_dataRootDir;
dataDir = [dataRootDir '/rest/'];
casename = 'R8_BF05WD05WV12_ST01_ens02';
tttttt = '028801';
% u file name
pattern = 'u\.[a-z]{2}\.[a-z]{3}[0-9]{3}';
dir = [dataDir casename '/' tttttt '/'];
[s0, r0] = system(['ls ' dir]);
if ~s0
    r1 = regexp(r0, pattern, 'match');
    fileIn = r1{1};
end

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

% set up profile data parameters
pflDir = [dataRootDir '/hist'];
fPrefix = 'his.mp.vis';
% creating object ncarlesPflData
f = ncarlesPflData('caseName', casename,...
                   'dataDir',  pflDir,...
                   'fPrefix',  fPrefix,...
                   'tStart',   1,...
                   'tEnd',     28801);

tmp = f.getVar('wtle');
wbpflts = squeeze(tmp(:,:,end)).*batag;

tmp = f.getVar('wps');
wpspflts = squeeze(tmp(:,end));

w = squeeze(var(3,:,:,:));
wxym = mean(mean(w,1),2);
wp = w-wxym;
wpsm = squeeze(mean(mean(wp.^2,1),2));

b = squeeze(var(4,:,:,:)).*batag;
bxym = mean(mean(b,1),2);
bp = b-bxym;
bpw = bp;
bpw(:,:,1:end-1) = 0.5.*(bp(:,:,1:end-1)+bp(:,:,2:end));
wpbp = wp.*bpw;
wpbpm = squeeze(mean(mean(wp.*bpw,1),2));
wpbpm3 = squeeze(mean(mean(wp.*bp,1),2));
idx = 2:nnz-2;
idxp1 = idx+1;
idxp2 = idx+2;
idxm1 = idx-1;
tmp = zeros(size(w));
tmp2 = zeros(size(w));
tmp(:,:,idx) = max(-w(:,:,idx), 0).*(b(:,:,idx)...
    +fluxLimiter(b(:,:,idxp1), b(:,:,idx), b(:,:,idxm1)));
tmp2(:,:,idx) = min(-w(:,:,idx), 0).*(b(:,:,idxp1)...
    +fluxLimiter(b(:,:,idx), b(:,:,idxp1), b(:,:,idxp2)));
tmp3 = -tmp-tmp2;
tmp3(:,:,1) = wpbp(:,:,1);
tmp3(:,:,end-2:end) = wpbp(:,:,end-2:end);
wpbpm2 = squeeze(mean(mean(tmp3,1),2));

figure;
plot(wpspflts, zw,'-k');
hold on;
plot(wpsm, zw, '+r');
xlabel('wps');
ylabel('depth');
ylim([-80, 0]);

figure;
p1 = plot(wbpflts, zw, '-k');
hold on;
p2 = plot(wpbpm, zw, 'xr');
p3 = plot(wpbpm2, zw, 'ob');
p4 = plot(wpbpm3, zw, '*c');
xlabel('wb');
ylabel('depth');
ylim([-80, 0]);
lg = legend([p1,p2,p3,p4], 'Pfl data', '3D field, no flux limiter',...
    '3D field, flux limiter', '3D field, mismatch z levels');
lg.Location = 'NorthWest';