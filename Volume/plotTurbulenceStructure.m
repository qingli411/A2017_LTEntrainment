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
casename = 'R8_BF05WD05WV12_ST01_ens03';
iwave = 1;  % 1 for no wave (WV00); 2 for Langmuir case (WV11)

% flags
l_save_fig = 1;     % 1: save figures (.fig and .eps)
l_update_RQ = 0;       % 1: reculculate all the R-Q data
l_update_pfl = 0;       % 1: reculculate all the profile data
l_lagrangian = 0;   % 1: account for Stokes drift in u

% flags for figures
l_wlevel = 0;       % w at different depth and spectrum
l_anisotropic = 1;  % anisotropic barycentric map
l_invariant = 0;    % isosurface of second invariant of velocity gradient
l_vorticity = 0;    % isosurface of vorticity modulus
l_RQmap = 0;        % RQ maps
l_tke = 0;          % TKE budget
l_Reynolds = 0;     % Reynolds stress budget

% check if mean profile data is needed
if l_lagrangian || l_anisotropic || l_Reynolds
    l_profile = 1;  % read profile data
else
    l_profile = 0;
end

% z levels
nnzcutoff = 128;

% inertial period
tio = 6.28e4;

%% set parameters
get_dataRootDir;    % get dataRootDir
dataDir = [dataRootDir '/rest/'];
% [s0, r0] = system(['ls ' dataDir casename ]);
% if ~s0
%     tttttt = r0(1:6);
% end
tttttt = '032401';
% u file name
inDir = [dataDir casename '/' tttttt];
[s0, r0] = system(['ls ' inDir '/u.mp.vis???']);
if ~s0
    [s1, r1] = system(['basename ' r0]);
    if ~s1
        fileIn = r1(1:end-1);
    end
end
% p file name
[s0, r0] = system(['ls ' inDir '/p.mp.vis???']);
if ~s0
    [s1, r1] = system(['basename ' r0]);
    if ~s1
        fileIn2 = r1(1:end-1);
    end
end
clear s0 r0 s1 r1;
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
rlat = getValue(srcname, 'rlat');   % latitude in degrees
fcor = 4.*pi.*sind(rlat)./86400;    % Coriolis parameter
nvar = 4+nscl; % number of variables
dx = Lx/nnx;
x = dx/2:dx:Lx;     % center of the grid, x
dy = Ly/nny;
y = dy/2:dy:Ly;     % center of the grid, y
dz = Lz/nnz;
zw  = dz:dz:Lz;     % z for w
zu  = dz/2:dz:Lz;   % z for u, v, t
batag = 9.81/5000;


% set up profile data parameters
dataDir = [dataRootDir '/hist/'];
fPrefix = 'his.mp.vis';
% [s0, r0] = system(['ls ' dataDir '/' casename ]);
% disp(['Reading profile data: ' r0(1:27)]);
% if s0==0
%     istart_str = r0(12:17);
%     iend_str = r0(19:24);
%     tStart = str2double(istart_str);
%     tEnd = str2double(iend_str);
% end
tStart = 1;
% tEnd = str2double(tttttt);
tEnd = 28801;
% creating object ncarlesPflData
f = ncarlesPflData('caseName', casename,...
                   'dataDir',  dataDir,...
                   'fPrefix',  fPrefix,...
                   'tStart',   tStart,...
                   'tEnd',     tEnd);
par = f.getParameters;

%% get data
if l_update_RQ || l_update_pfl     % update data
    % read uvwb
    fid = fopen([inDir '/' fileIn]);
    tmp = fread(fid, nvar*nnx*nny*nnz, 'double', 'l');
    fclose(fid);
    clear var fid;
    var = reshape(tmp, nvar, nnx, nny, nnz);
    clear tmp;

    % u, v, w, b, cut off the lower part of the domain
    nnzmax = min(nnzcutoff, nnz);
    u = squeeze(var(1,:,:,1:nnzmax));
    v = squeeze(var(2,:,:,1:nnzmax));
    w = squeeze(var(3,:,:,1:nnzmax));
    b = squeeze(var(4,:,:,1:nnzmax)).*batag;
    e = squeeze(var(5,:,:,1:nnzmax));
    clear var;
    % w and e on u level
    tmp = zeros([nnx, nny, nnzmax]);
    tmp(:,:,2:end) = w(:,:,1:end-1);
    w_zu = 0.5.*(tmp+w);
    tmp(:,:,2:end) = e(:,:,1:end-1);
    e_zu = 0.5.*(tmp+e);
    clear tmp;

    % get Stokes drift
    tmp = f.getVar('stokes');
    stokes = tmp(1:nnzmax, end)'; % Stokes drift at the last time step
    clear tmp;
    dir_x = f.getVar('dir_x');
    dir_y = f.getVar('dir_y');

    % read pressure
    fid = fopen([inDir '/' fileIn2]);
    tmp = fread(fid, nnx*nny*nnz, 'double', 'l');
    fclose(fid);
    clear var fid;
    var = reshape(tmp, nnx, nny, nnz);
    clear tmp;
    p = squeeze(var(:,:,1:nnzmax));
    clear var;

    % lagrangian u and v
    ul = u;
    vl = v;
    for i=1:nnx
        for j=1:nny
            ul(i,j,:) = squeeze(u(i,j,:))' + stokes(1,:).*dir_x;
            vl(i,j,:) = squeeze(v(i,j,:))' + stokes(1,:).*dir_y;
        end
    end
    % correct p (pstar -> p)
    p = p-e_zu./3-0.5.*(ul.^2+vl.^2+w_zu.^2);

    % get Lagrangian u if required
    if l_lagrangian
        u = ul;
        v = vl;
    end

end

%% calculate the vorticity and invariants of velocity gradient
if l_update_RQ
    % velocity gradients
    velgrad = zeros([3,3,nnx,nny,nnzmax]);
    velgrad(1,1,:,:,:) = get_ddx(u, dx);            % dudx
    velgrad(1,2,:,:,:) = get_ddy(u, dy);            % dudy
    velgrad(1,3,:,:,:) = get_ddz_u(u, dz, 'u');     % dudz
    velgrad(2,1,:,:,:) = get_ddx(v, dx);            % dvdx
    velgrad(2,2,:,:,:) = get_ddy(v, dy);            % dvdy
    velgrad(2,3,:,:,:) = get_ddz_u(v, dz, 'u');     % dvdz
    velgrad(3,1,:,:,:) = get_ddx(w_zu, dx);         % dwdx
    velgrad(3,2,:,:,:) = get_ddy(w_zu, dy);         % dwdy
    velgrad(3,3,:,:,:) = get_ddz_w(w, dz);          % dwdz

    % vorticity
    omega = zeros([3, nnx, nny, nnzmax]);
    omega(1,:,:,:) = velgrad(3,2,:,:,:)-velgrad(2,3,:,:,:);
    omega(2,:,:,:) = velgrad(1,3,:,:,:)-velgrad(3,1,:,:,:);
    omega(3,:,:,:) = velgrad(2,1,:,:,:)-velgrad(1,2,:,:,:);

    % second and third invariants of the velocity gradient tensor
    Q = zeros(nnx,nny,nnzmax);
    R = zeros(nnx,nny,nnzmax);
    for i=1:nnx
        for j=1:nny
            for k=1:nnzmax
                [~,Q(i,j,k),R(i,j,k)] = ...
                    invariant3(squeeze(velgrad(:,:,i,j,k)));
            end
        end
    end
    Q = -0.5.*Q;
    R = -1./3.*R;

    % symmetric & anti-symmetric component of the velocity gradiant tensor
    velgrad_s = zeros(size(velgrad));
    for i=1:3
        for j=1:3
            velgrad_s(i,j,:,:,:) = ...
                0.5.*(velgrad(i,j,:,:,:)+velgrad(j,i,:,:,:));
        end
    end
    % second and third invariants
    QS = zeros(nnx,nny,nnzmax);
    RS = zeros(nnx,nny,nnzmax);
    for i=1:nnx
        for j=1:nny
            for k=1:nnzmax
                [~,QS(i,j,k),RS(i,j,k)] = ...
                    invariant3(squeeze(velgrad_s(:,:,i,j,k)));
            end
        end
    end
    QS = -0.5.*QS;
    RS = -1./3.*RS;
    QA = Q-QS;

    % save mat file
    matfile = [inDir '/data.mat'];
    save(matfile, 'omega', 'Q', 'R', 'QS', 'RS', 'QA',...
         'zu', 'zw', 'x', 'y');
else
    matfile = [inDir '/data.mat'];
    load(matfile);
end

%% calculate mean profiles
if l_update_pfl
    % mean buoyancy profile
    bxym0 = mean(mean(b,1),2);
    N2 = squeeze(get_ddz_u(bxym0, dz, 'w'));
    bxym = squeeze(bxym0);
    clear bxym0;
    % get the index of boundary layer base
    [~, ind_hb] = max(N2);

    % boundary layer depth
    hb = -zw(ind_hb);
    % surface wind stress
    tmp = f.getVar('utau');
    utau = tmp(end);

    % mean horizontal velocity and shear
    um = mean(mean(u,1),2);
    vm = mean(mean(v,1),2);
    uxym = squeeze(um);
    vxym = squeeze(vm);
    dumdz = get_ddz_u(um, dz, 'u');
    dvmdz = get_ddz_u(vm, dz, 'u');
    dudzxym = squeeze(dumdz);
    dvdzxym = squeeze(dvmdz);

    % Stokes drift shear
    tmpus = zeros(size(um));
    tmpus(1,1,:) = stokes;
    dstokesdzm = squeeze(get_ddz_u(tmpus, dz, 'u'));
    clear tmpus;

    % mean Reynolds stress
    up = u-um;
    vp = v-vm;
    clear uxym0 vxym0;
    wp = w-mean(mean(w,1),2);
    wpu = w_zu-mean(mean(w_zu,1),2);
    upsxym = squeeze(mean(mean(up.*up,1),2));
    vpsxym = squeeze(mean(mean(vp.*vp,1),2));
    wpsxym = squeeze(mean(mean(wp.*wp,1),2));
    wpsuxym = squeeze(mean(mean(wpu.*wpu,1),2));    % u level
    upvpxym = squeeze(mean(mean(up.*vp,1),2));
    upwpxym = squeeze(mean(mean(up.*wpu,1),2));     % u level
    vpwpxym = squeeze(mean(mean(vp.*wpu,1),2));     % u level
    [~, ind_wps] = max(wpsxym);

    % Reynolds stress budget
    pflsize = size(wpsuxym);
    % buoyancy flux
    bp = b-mean(mean(b,1),2);
    B11 = zeros(pflsize);
    B22 = zeros(pflsize);
    B33 = 2.*squeeze(mean(mean(wpu.*bp,1),2));
    B12 = zeros(pflsize);
    B13 = squeeze(mean(mean(up.*bp,1),2));
    B23 = squeeze(mean(mean(vp.*bp,1),2));

    % shear production
    SP11 = -2.*upwpxym.*dudzxym;
    SP22 = -2.*vpwpxym.*dvdzxym;
    SP33 = zeros(pflsize);
    SP12 = -upwpxym.*dvdzxym-vpwpxym.*dudzxym;
    SP13 = -wpsuxym.*dudzxym;
    SP23 = -wpsuxym.*dvdzxym;

    % Stokes production
    StP11 = zeros(pflsize);
    StP22 = zeros(pflsize);
    StP33 = -2.*upwpxym.*dstokesdzm.*dir_x-2.*vpwpxym.*dstokesdzm.*dir_y;
    StP12 = zeros(pflsize);
    StP13 = -upsxym.*dstokesdzm.*dir_x-upvpxym.*dstokesdzm.*dir_y;
    StP23 = -upvpxym.*dstokesdzm.*dir_x-vpsxym.*dstokesdzm.*dir_y;

    % transport
    T11 = -squeeze(get_ddz_u(mean(mean(up.*up.*wpu,1),2),dz,'u'));
    T22 = -squeeze(get_ddz_u(mean(mean(vp.*vp.*wpu,1),2),dz,'u'));
    T33 = -squeeze(get_ddz_u(mean(mean(wpu.*wpu.*wpu,1),2),dz,'u'));
    T12 = -squeeze(get_ddz_u(mean(mean(up.*vp.*wpu,1),2),dz,'u'));
    T13 = -squeeze(get_ddz_u(mean(mean(up.*wpu.*wpu,1),2),dz,'u'));
    T23 = -squeeze(get_ddz_u(mean(mean(vp.*wpu.*wpu,1),2),dz,'u'));

    % Coriolis
    C11 = 2.*fcor.*upvpxym;
    C22 = -C11;
    C33 = zeros(pflsize);
    C12 = fcor.*(vpsxym-upsxym);
    C13 = fcor.*vpwpxym;
    C23 = -fcor.*upwpxym;

    % pressure strain
    pp = p-mean(mean(p,1),2);
    dpdx = get_ddx(pp,dx);
    dpdy = get_ddy(pp,dy);
    dpdz = get_ddz_u(pp,dz,'u');
    PS11 = -2.*squeeze(mean(mean(up.*dpdx,1),2));
    PS22 = -2.*squeeze(mean(mean(vp.*dpdy,1),2));
    PS33 = -2.*squeeze(mean(mean(wpu.*dpdz,1),2));
    PS12 = -squeeze(mean(mean(up.*dpdy+vp.*dpdx,1),2));
    PS13 = -squeeze(mean(mean(up.*dpdz+wpu.*dpdx,1),2));
    PS23 = -squeeze(mean(mean(vp.*dpdz+wpu.*dpdy,1),2));

    % correct pressure strain term and transport terms
%     dupdxm = squeeze(get_ddx(mean(mean(up.*pp,1),2),dx));
%     dvpdym = squeeze(get_ddy(mean(mean(vp.*pp,1),2),dy));
%     dwpdzm = squeeze(get_ddz_u(mean(mean(wpu.*pp,1),2),dz,'u'));
%     ptrans = 2/3.*(dupdxm+dvpdym+dwpdzm);
%     ptrans = -(PS11+PS22+PS33)./3;
%
%     PS11 = PS11+ptrans;
%     PS22 = PS22+ptrans;
%     PS33 = PS33+ptrans;
%     T11 = T11-ptrans;
%     T22 = T22-ptrans;
%     T33 = T33-ptrans;

    % dissipation to SGS
    dsl = (abs(2.25*dx*dy*dz))^(1/3);
    almin_c = 0.0001;
    stabmin = 1e-12;
    stab_c = 0.76;
    ck = 0.1;
    ceps = 0.93;
    vk = 0.4;
    csmag = sqrt(ck*sqrt(ck/ceps));
    dslk = zeros(size(e));
    for i=1:nnzmax
        dslk(:,:,i) = min(dsl,vk*abs(zw(i))/csmag);
    end
    almin = almin_c*dsl;
    dbdz = get_ddz_u(b,dz,'u');
    alk = dslk;
    inds = find(dbdz>stabmin);
    alk(inds) = stab_c.*sqrt(e(inds)./dbdz(inds));
    inds = find(dslk<alk);
    alk(inds) = dslk(inds);
    alk(alk<almin) = almin;
    vism = ck.*alk.*sqrt(e); % eddy viscosity
%     clear alk dbdz;
    dudx = get_ddx(u,dx);
    dudy = get_ddy(u,dy);
    dudz = get_ddz_u(u, dz, 'u');
    dvdx = get_ddx(v,dx);
    dvdy = get_ddy(v,dy);
    dvdz = get_ddz_u(v, dz, 'u');
    dwdx = get_ddx(w_zu,dx);
    dwdy = get_ddy(w_zu,dy);
    s11 = dudx;
    s22 = dvdy;
    s33 = get_ddz_w(w,dz);
    s12 = 0.5.*(dudy+dvdx);
    s13 = 0.5.*(dudz+get_ddx(w_zu,dx));
    s23 = 0.5.*(get_ddy(w_zu,dy)+dvdz);
    uzfluc = dudz - dumdz;
    vzfluc = dvdz - dvmdz;
    D11 = -4.*squeeze(mean(mean(vism.*(s11.^2+s12.*dudy+s13.*uzfluc),1),2));
    D22 = -4.*squeeze(mean(mean(vism.*(s12.*dvdx+s22.^2+s23.*vzfluc),1),2));
    D33 = -4.*squeeze(mean(mean(vism.*(s13.*dwdx+s23.*dwdy+s33.^2),1),2));
    D12 = -2.*squeeze(mean(mean(vism.*(s11.*dvdx+s12.*s22+s13.*vzfluc+...
        s12.*s11+s22.*dudy+s23.*uzfluc),1),2));
    D13 = -2.*squeeze(mean(mean(vism.*(s11.*dwdx+s12.*dwdy+s13.*s33+...
        s13.*s11+s23.*dudy+s33.*uzfluc),1),2));
    D23 = -2.*squeeze(mean(mean(vism.*(s12.*dwdx+s22.*dwdy+s23.*s33+...
        s13.*dvdx+s23.*s22+s33.*vzfluc),1),2));

    % stress transport
    TV11 = 4.*squeeze(get_ddz_u(mean(mean(vism.*s13.*up,1),2),'u'));
    TV22 = 4.*squeeze(get_ddz_u(mean(mean(vism.*s23.*vp,1),2),'u'));
    TV33 = 4.*squeeze(get_ddz_u(mean(mean(vism.*s33.*wpu,1),2),'u'));
    TV12 = 2.*squeeze(get_ddz_u(mean(mean(vism.*s13.*vp,1),2),'u'))+...
           2.*squeeze(get_ddz_u(mean(mean(vism.*s23.*up,1),2),'u'));
    TV13 = 2.*squeeze(get_ddz_u(mean(mean(vism.*s13.*wpu,1),2),'u'))+...
           2.*squeeze(get_ddz_u(mean(mean(vism.*s33.*up,1),2),'u'));
    TV23 = 2.*squeeze(get_ddz_u(mean(mean(vism.*s23.*wpu,1),2),'u'))+...
           2.*squeeze(get_ddz_u(mean(mean(vism.*s33.*vp,1),2),'u'));

    % save mat file
    matfile = [inDir '/pfldata.mat'];
    save(matfile,...
        'upsxym', 'vpsxym', 'wpsxym', 'upvpxym', 'upwpxym', 'vpwpxym',...
        'wpsuxym', 'bxym', 'uxym', 'vxym', 'dudzxym', 'dvdzxym', 'stokes',...
        'dstokesdzm', 'ind_hb', 'ind_wps', 'zu', 'zw', 'x', 'y',...
        'B11', 'B22', 'B33', 'B12', 'B13', 'B23',...
        'SP11', 'SP22', 'SP33', 'SP12', 'SP13', 'SP23',...
        'StP11', 'StP22', 'StP33', 'StP12', 'StP13', 'StP23',...
        'T11', 'T22', 'T33', 'T12', 'T13', 'T23',...
        'C11', 'C22', 'C33', 'C12', 'C13', 'C23',...
        'PS11', 'PS22', 'PS33', 'PS12', 'PS13', 'PS23',...
        'D11', 'D22', 'D33', 'D12', 'D13', 'D23',...
        'TV11', 'TV22', 'TV33', 'TV12', 'TV13', 'TV23',...
        'dir_x', 'dir_y');
else
    matfile = [inDir '/pfldata.mat'];
    load(matfile);
end

%% define subdomains
% horizontal
indsx1 = 1:nnx;
indsy1 = 1:nny;
indsx2 = 1:nnx;
indsy2 = 1:nny;
indsx3 = 1:nnx;
indsy3 = 1:nny;
indsx4 = 1:nnx;
indsy4 = 1:nny;

% vertical: upper 25%, middle 50% and lower 25% of the boundary layer
nnz_sub1 = floor(ind_hb/4);
nnz_sub2 = floor(ind_hb/4);
nnz_sub3 = floor(ind_hb/4);
nnz_half = floor(ind_hb/2);
indsz1 = 1:nnz_sub1;                    % upper 25%
indsz2 = ind_hb-nnz_sub2+1:ind_hb;      % lower 25%
indsz3 = nnz_sub1+1:nnz_half;           % middle upper 25%
indsz4 = nnz_half+1:ind_hb-nnz_sub2;    % middle lower 25%

if l_invariant || l_vorticity || l_RQmap
    % x, y, z mesh
    [xx1,yy1,zz1] = meshgrid(x(indsx1), y(indsy1), zu(indsz1));
    [xx2,yy2,zz2] = meshgrid(x(indsx2), y(indsy2), zu(indsz2));
    [xx3,yy3,zz3] = meshgrid(x(indsx3), y(indsy3), zu(indsz3));
    [xx4,yy4,zz4] = meshgrid(x(indsx4), y(indsy4), zu(indsz4));

    % vorticity, invariants of velocity gradient in the subdomains
    omega_sub1 = omega(:,indsx1, indsy1, indsz1);
    Q_sub1 = Q(indsx1, indsy1, indsz1);
    R_sub1 = R(indsx1, indsy1, indsz1);
    QS_sub1 = QS(indsx1, indsy1, indsz1);
    RS_sub1 = RS(indsx1, indsy1, indsz1);
    QA_sub1 = QA(indsx1, indsy1, indsz1);
    w_sub1 = w(indsx1, indsy1, indsz1);
    omega_sub2 = omega(:,indsx2, indsy2, indsz2);
    Q_sub2 = Q(indsx2, indsy2, indsz2);
    R_sub2 = R(indsx2, indsy2, indsz2);
    QS_sub2 = QS(indsx2, indsy2, indsz2);
    RS_sub2 = RS(indsx2, indsy2, indsz2);
    QA_sub2 = QA(indsx2, indsy2, indsz2);
    w_sub2 = w(indsx2, indsy2, indsz2);
    omega_sub3 = omega(:,indsx3, indsy3, indsz3);
    Q_sub3 = Q(indsx3, indsy3, indsz3);
    R_sub3 = R(indsx3, indsy3, indsz3);
    QS_sub3 = QS(indsx3, indsy3, indsz3);
    RS_sub3 = RS(indsx3, indsy3, indsz3);
    QA_sub3 = QA(indsx3, indsy3, indsz3);
    w_sub3 = w(indsx3, indsy3, indsz3);
    omega_sub4 = omega(:,indsx4, indsy4, indsz4);
    Q_sub4 = Q(indsx4, indsy4, indsz4);
    R_sub4 = R(indsx4, indsy4, indsz4);
    QS_sub4 = QS(indsx4, indsy4, indsz4);
    RS_sub4 = RS(indsx4, indsy4, indsz4);
    QA_sub4 = QA(indsx4, indsy4, indsz4);
    w_sub4 = w(indsx4, indsy4, indsz4);
end

%% figure 1: w at different levels and 1d spectrum
if l_wlevel
    [xx2d,yy2d] = meshgrid(x, y);
    figure;
    wdat = permute(w(:,:,1:ind_hb), [2, 1, 3]);
    cmax = max(abs(wdat(:)))*0.3;
    % w at the level of maximum wps
    subplot(2, 3, 1)
    ps1 = pcolor(xx2d,yy2d,wdat(:,:,ind_wps));
    shading flat;
    caxis([-cmax, cmax]);
    daspect([1 1 1]);
    xlabel('x (m)');
    ylabel('y (m)');
    title(['z = ' sprintf('%4.2f', zw(ind_wps)) ' m']);
    axis tight
    % w at the mid-level of the boundary layer
    subplot(2, 3, 2)
    ind_mhb = floor(ind_hb/2);
    ps2 = pcolor(xx2d,yy2d,wdat(:,:,ind_mhb));
    shading flat;
    % colorbar;
    caxis([-cmax, cmax]);
    daspect([1 1 1]);
    xlabel('x (m)');
    % ylabel('y (m)');
    title(['z = ' sprintf('%4.2f', zw(ind_mhb)) ' m']);
    axis tight
    % w at the level of maximum entrainment (z~0.9 hb)
    subplot(2, 3, 3)
    ind_wb = floor(ind_hb*0.9);
    ps3 = pcolor(xx2d,yy2d,wdat(:,:,ind_wb));
    shading flat;
    caxis([-cmax, cmax]);
    daspect([1 1 1]);
    xlabel('x (m)');
    % ylabel('y (m)');
    title(['z = ' sprintf('%4.2f', zw(ind_wb)) ' m']);
    axis tight
    % save figure
    figname = [outDir '/w_level.fig'];
    saveas(gcf, figname, 'fig');

    % save w in seperate files
    figure;
    pcolor(xx2d,yy2d,wdat(:,:,ind_wps));
    shading flat;
    caxis([-cmax, cmax]);
    daspect([1 1 1]);
    xlabel('x (m)');
    ylabel('y (m)');
    title(['z = ' sprintf('%4.2f', zw(ind_wps)) ' m']);
    axis tight
    figname = [outDir '/w_level_maxwps.fig'];
    saveas(gcf, figname, 'fig');
    postProcessFig(figname, 6, 4);

    figure;
    pcolor(xx2d,yy2d,wdat(:,:,ind_wb));
    shading flat;
    caxis([-cmax, cmax]);
    daspect([1 1 1]);
    xlabel('x (m)');
    % ylabel('y (m)');
    title(['z = ' sprintf('%4.2f', zw(ind_wb)) ' m']);
    axis tight
    figname = [outDir '/w_level_minwb.fig'];
    saveas(gcf, figname, 'fig');
    postProcessFig(figname, 6, 4);
end

%% figure 2: anisotropy triangle
if l_anisotropic
    % figure 2A: anisotropy triangle
    cb_label = '$x_3/h_b$';
    figure;
    figPropertySingle;
    c = zeros([ind_hb,3]);
    for i=1:ind_hb
        a = anisotropyTensor( upsxym(i), vpsxym(i), wpsuxym(i),...
                              upvpxym(i), upwpxym(i), vpwpxym(i));
        c(i,:) = barycentricCoord(a);
    end
    if par.fcor > 0
        [c2,~] = f.getAnisotropicBarycentricCoord('inertial',[0,1]);
    else
        [c2,~] = f.getAnisotropicBarycentricCoord('second',[par.tend-tio,par.tend]);
    end
    plotAnisotropicBarycentricMap2(c,c2(1:ind_hb,:),...
        -zu(1:ind_hb)./zu(ind_hb),1,cb_label);
    if l_save_fig
        % save figure
        figname = [outDir '/anisotropy.fig'];
        saveas(gcf, figname, 'fig');
        [inDir, name, ~] = fileparts(figname);
        print('-depsc2',[inDir '/' name]);
    end

    % figure 2B: direction of anisotropy versus direction of
    %   shear and velocity
    figure;
    figPropertySingle;
    cmin = zeros([ind_hb,3]);
    cmax = zeros([ind_hb,3]);
    lambda = zeros([ind_hb,3]);
    cmin2 = zeros([ind_hb,3]);
    cmax2 = zeros([ind_hb,3]);
    lambda2 = zeros([ind_hb,3]);
    if par.fcor > 0
        [apfl,~] = f.getAnisotropicTensor('inertial',[0,1]);
    else
        [apfl,~] = f.getAnisotropicTensor('second',[par.tend-tio,par.tend]);
    end
    for i=1:ind_hb
        % instantaneous profile
        a = anisotropyTensor( upsxym(i), vpsxym(i), wpsuxym(i),...
                      upvpxym(i), upwpxym(i), vpwpxym(i));
        [lambda(i,:),cmax(i,:),cmin(i,:)] = eigMaxMin3(a);
        % mean profile
        a2 = squeeze(apfl(i,:,:));
        [lambda2(i,:),cmax2(i,:),cmin2(i,:)] = eigMaxMin3(a2);
    end
    plotEigenVectorDirectionMaxMin2(lambda, cmax, cmin,...
        lambda2, cmax2, cmin2,...
        -zu(1:ind_hb)./zw(ind_hb),0,cb_label);
    % plot the mean direction of Lagrangian and Eulerian shear direction
    % near the surface (upper 25% of the boundary layer)
    % 1a. mean profile
    if par.fcor > 0
        tmpv = f.getProfiles('inertial',[0,1]);
    else
        tmpv = f.getProfiles('second',[par.tend-tio,par.tend]);
    end
    dudzmpfl = tmpv.dudz;
    dvdzmpfl = tmpv.dvdz;
    px = mean(dudzmpfl(indsz1));
    pxl = mean(dudzmpfl(indsz1))+mean(dstokesdzm(indsz1)).*dir_x;
    py = mean(dvdzmpfl(indsz1));
    pyl = mean(dvdzmpfl(indsz1))+mean(dstokesdzm(indsz1)).*dir_y;
    pr = sqrt(px.^2+py.^2);
    prl = sqrt(pxl.^2+pyl.^2);
    scatter(pxl./prl,pyl./prl,50,'s','filled',...
        'MarkerFaceColor',cgray,...
        'MarkerEdgeColor',cgray);
    scatter(px./pr,py./pr,50,'d','filled',...
        'MarkerFaceColor',cgray,...
        'MarkerEdgeColor',cgray);
    % 1b. instantaneous profile
    px = mean(dudzxym(indsz1));
    pxl = mean(dudzxym(indsz1))+mean(dstokesdzm(indsz1)).*dir_x;
    py = mean(dvdzxym(indsz1));
    pyl = mean(dvdzxym(indsz1))+mean(dstokesdzm(indsz1)).*dir_y;
    pr = sqrt(px.^2+py.^2);
    prl = sqrt(pxl.^2+pyl.^2);
    scatter(pxl./prl,pyl./prl,50,'s','filled',...
        'MarkerFaceColor','k',...
        'MarkerEdgeColor','k');
    scatter(px./pr,py./pr,50,'d','filled',...
        'MarkerFaceColor','k',...
        'MarkerEdgeColor','k');
    % plot the mean Lagrangian and Eulerian velocity direction
    % near the surface (upper 25% of the boundary layer)
    % 2a. mean profile
    uxympfl = tmpv.uxym;
    vxympfl = tmpv.vxym;
    px = mean(uxympfl(indsz1));
    pxl = mean(uxympfl(indsz1))+mean(stokes(indsz1)).*dir_x;
    py = mean(vxympfl(indsz1));
    pyl = mean(vxympfl(indsz1))+mean(stokes(indsz1)).*dir_y;
    pr = sqrt(px.^2+py.^2);
    prl = sqrt(pxl.^2+pyl.^2);
    scatter(pxl./prl,pyl./prl,50,'s','LineWidth',1.0,...
        'MarkerEdgeColor',cgray);
    scatter(px./pr,py./pr,50,'d','LineWidth',1.0,...
        'MarkerEdgeColor',cgray);
    % 2b. instantaneous profile
    px = mean(uxym(indsz1));
    pxl = mean(uxym(indsz1))+mean(stokes(indsz1)).*dir_x;
    py = mean(vxym(indsz1));
    pyl = mean(vxym(indsz1))+mean(stokes(indsz1)).*dir_y;
    pr = sqrt(px.^2+py.^2);
    prl = sqrt(pxl.^2+pyl.^2);
    scatter(pxl./prl,pyl./prl,50,'s','LineWidth',1.0,...
        'MarkerEdgeColor','k');
    scatter(px./pr,py./pr,50,'d','LineWidth',1.0,...
        'MarkerEdgeColor','k');
    if l_save_fig
        % save figure
        figname = [outDir '/anisoDirection.fig'];
        saveas(gcf, figname, 'fig');
        [inDir, name, ~] = fileparts(figname);
        print('-depsc2',[inDir '/' name]);
    end
end

%% figure 3: second invariant of velocity gradient
if l_invariant
    figure;
    dat = permute(Q_sub1, [2, 1, 3]);  % [x, y, z] -> [y, x, z]
    angle_view = [20, -10];
    angle_cam = [180, 75];
    isovalue = prctile(dat(:),90);
    plot_isosurface(xx1, yy1, zz1, dat, isovalue, angle_view, angle_cam);
    clear dat;
    % save figure
    figname = [outDir '/Q_dom1.fig'];
    saveas(gcf, figname, 'fig');

    figure;
    dat = permute(Q_sub2, [2, 1, 3]);  % [x, y, z] -> [y, x, z]
    angle_view = [20, 20];
    angle_cam = [180, -75];
    isovalue = prctile(dat(:),95);
    plot_isosurface(xx2, yy2, zz2, dat, isovalue, angle_view, angle_cam);
    clear dat;
    figname = [outDir '/Q_dom2.fig'];
    saveas(gcf, figname, 'fig');

    figure;
    dat = permute(Q_sub3, [2, 1, 3]);  % [x, y, z] -> [y, x, z]
    angle_view = [20, -20];
    angle_cam = [180, 75];
    isovalue = prctile(dat(:),95);
    plot_isosurface(xx3, yy3, zz3, dat, isovalue, angle_view, angle_cam);
    clear dat;
    figname = [outDir '/Q_dom3.fig'];
    saveas(gcf, figname, 'fig');

%     figure;
%     dat = permute(Q_sub4, [2, 1, 3]);  % [x, y, z] -> [y, x, z]
%     angle_view = [20, -10];
%     angle_cam = [180, 75];
%     isovalue = prctile(dat(:),90);
%     plot_isosurface(xx4, yy4, zz4, dat, isovalue, angle_view, angle_cam);
%     clear dat;
%     figname = [casename '/Q_dom4.fig'];
%     saveas(gcf, figname, 'fig');
end

%% figure 4: vorticity modulus
if l_vorticity
    figure;
    momega = squeeze(sqrt(omega_sub1(1,:,:,:).^2+...
        omega_sub1(2,:,:,:).^2+omega_sub1(3,:,:,:).^2));
    dat = permute(momega, [2, 1, 3]);   % [x, y, z] -> [y, x, z]
    angle_view = [20, -10];
    angle_cam = [180, 75];
    isovalue = prctile(dat(:),90);
    plot_isosurface(xx1, yy1, zz1, dat, isovalue, angle_view, angle_cam);
    clear dat;
    figname = [outDir '/VM_dom1.fig'];
    saveas(gcf, figname, 'fig');

    figure;
    momega = squeeze(sqrt(omega_sub2(1,:,:,:).^2+...
        omega_sub2(2,:,:,:).^2+omega_sub2(3,:,:,:).^2));
    dat = permute(momega, [2, 1, 3]);   % [x, y, z] -> [y, x, z]
    angle_view = [20, 20];
    angle_cam = [180, -75];
    isovalue = prctile(dat(:),95);
    plot_isosurface(xx2, yy2, zz2, dat, isovalue, angle_view, angle_cam);
    clear dat;
    figname = [outDir '/VM_dom2.fig'];
    saveas(gcf, figname, 'fig');

    figure;
    momega = squeeze(sqrt(omega_sub3(1,:,:,:).^2+...
        omega_sub3(2,:,:,:).^2+omega_sub3(3,:,:,:).^2));
    dat = permute(momega, [2, 1, 3]);   % [x, y, z] -> [y, x, z]
    angle_view = [20, -20];
    angle_cam = [180, 75];
    isovalue = prctile(dat(:),95);
    plot_isosurface(xx3, yy3, zz3, dat, isovalue, angle_view, angle_cam);
    clear dat;
    figname = [outDir '/VM_dom3.fig'];
    saveas(gcf, figname, 'fig');

%     figure;
%     momega = squeeze(sqrt(omega_sub4(1,:,:,:).^2+...
%         omega_sub4(2,:,:,:).^2+omega_sub4(3,:,:,:).^2));
%     dat = permute(momega, [2, 1, 3]);   % [x, y, z] -> [y, x, z]
%     angle_view = [20, -10];
%     angle_cam = [180, 75];
%     isovalue = prctile(dat(:),90);
%     plot_isosurface(xx4, yy4, zz4, dat, isovalue, angle_view, angle_cam);
%     clear dat;
%     figname = [casename '/VM_dom4.fig'];
%     saveas(gcf, figname, 'fig');
end

%% figure 5: R-Q map
if l_RQmap
    % dom 1
    SS_mean = -2.0.*mean(QS_sub1(:));   % S_{ij}S_{ij}
    inds_plt = find(w_sub1~=0);

    figure;
    xdata = R_sub1(inds_plt);
    ydata = Q_sub1(inds_plt);
    p = plot_RQmap(xdata, ydata, SS_mean);
    save_figure_RQ(p, [outDir '/RQ_dom1.fig'], l_save_fig);

    figure;
    xdata = RS_sub1(inds_plt);
    ydata = QS_sub1(inds_plt);
    p = plot_RSQSmap(xdata, ydata, SS_mean);
    save_figure_RQ(p, [outDir '/RSQS_dom1.fig'], l_save_fig);

    figure;
    xdata = QA_sub1(inds_plt);
    ydata = -QS_sub1(inds_plt);
    p = plot_QAQSmap(xdata, ydata, SS_mean);
    save_figure_RQ(p, [outDir '/QAQS_dom1.fig'], l_save_fig);

    % dom 2
    SS_mean = -2.0.*mean(QS_sub2(:));   % S_{ij}S_{ij}
    inds_plt = find(w_sub2~=0);

    figure;
    xdata = R_sub2(inds_plt);
    ydata = Q_sub2(inds_plt);
    p = plot_RQmap(xdata, ydata, SS_mean);
    save_figure_RQ(p, [outDir '/RQ_dom2.fig'], l_save_fig);

    figure;
    xdata = RS_sub2(inds_plt);
    ydata = QS_sub2(inds_plt);
    p = plot_RSQSmap(xdata, ydata, SS_mean);
    save_figure_RQ(p, [outDir '/RSQS_dom2.fig'], l_save_fig);

    figure;
    xdata = QA_sub2(inds_plt);
    ydata = -QS_sub2(inds_plt);
    p = plot_QAQSmap(xdata, ydata, SS_mean);
    save_figure_RQ(p, [outDir '/QAQS_dom2.fig'], l_save_fig);

    % dom 3
    SS_mean = -2.0.*mean(QS_sub3(:));   % S_{ij}S_{ij}
    inds_plt = find(w_sub3~=0);

    figure;
    xdata = R_sub3(inds_plt);
    ydata = Q_sub3(inds_plt);
    p = plot_RQmap(xdata, ydata, SS_mean);
    save_figure_RQ(p, [outDir '/RQ_dom3.fig'], l_save_fig);

    figure;
    xdata = RS_sub3(inds_plt);
    ydata = QS_sub3(inds_plt);
    p = plot_RSQSmap(xdata, ydata, SS_mean);
    save_figure_RQ(p, [outDir '/RSQS_dom3.fig'], l_save_fig);

    figure;
    xdata = QA_sub3(inds_plt);
    ydata = -QS_sub3(inds_plt);
    p = plot_QAQSmap(xdata, ydata, SS_mean);
    save_figure_RQ(p, [outDir '/QAQS_dom3.fig'], l_save_fig);
    
    % dom 4
    SS_mean = -2.0.*mean(QS_sub4(:));   % S_{ij}S_{ij}
    inds_plt = find(w_sub4~=0);

    figure;
    xdata = R_sub4(inds_plt);
    ydata = Q_sub4(inds_plt);
    p = plot_RQmap(xdata, ydata, SS_mean);
    save_figure_RQ(p, [outDir '/RQ_dom4.fig'], l_save_fig);

    figure;
    xdata = RS_sub4(inds_plt);
    ydata = QS_sub4(inds_plt);
    p = plot_RSQSmap(xdata, ydata, SS_mean);
    save_figure_RQ(p, [outDir '/RSQS_dom4.fig'], l_save_fig);

    figure;
    xdata = QA_sub4(inds_plt);
    ydata = -QS_sub4(inds_plt);
    p = plot_QAQSmap(xdata, ydata, SS_mean);
    save_figure_RQ(p, [outDir '/QAQS_dom4.fig'], l_save_fig);
end

%% figure 6: TKE budget
if l_tke
    wgt = utau.^3./hb;
    % snapshot
    tmp = f.getVar('t_rprod');
    MSP = tmp(1:nnzmax,end)./wgt;
    tmp = f.getVar('t_stokes');
    MStP = tmp(1:nnzmax,end)./wgt;
    tmp = f.getVar('wtle');
    MB = tmp(1:nnzmax,end).*batag./wgt;
    tmp = f.getVar('t_wq');
    MT = tmp(1:nnzmax,end)./wgt;
    tmp = f.getVar('t_wp');
    MPS = tmp(1:nnzmax,end)./wgt;
    tmp = f.getVar('t_tau');
    MTV = tmp(1:nnzmax,end)./wgt;
    MC = zeros(size(MSP));
    tmp = f.getVar('t_dsle');
    MD = -tmp(1:nnzmax,end)./wgt;

%     % mean profile
%     tmpv = f.getProfiles('inertial',[0,1]);
%     MSP = tmpv.t_rprod./wgt;
%     MStP = tmpv.t_stokes./wgt;
%     MB = tmpv.wtle.*batag./wgt;
%     MT = tmpv.t_wq./wgt;
%     MPS = tmpv.t_wp./wgt;
%     MTV = (tmpv.t_tran-tmpv.t_wq-tmpv.t_wp)./wgt;
%     MC = zeros(size(MSP));
%     MD = -tmpv.t_dsle./wgt;

%     % correct pressure strain and transport
%     MT = MT-ptrans./wgt;
%     MPS = MPS+ptrans./wgt;

    figure;
    hold on;
    indz = 2:ind_hb;
    zpfl = zu(indz)./hb;
    p1 = pflplot(MSP(indz),  zpfl, '-b');
    p2 = pflplot(MStP(indz), zpfl, '-r');
    p3 = pflplot(MB(indz),   zpfl, '-k');
    p4 = pflplot(MT(indz),   zpfl, '-', 'Color', ccyan);
    p5 = pflplot(MPS(indz),  zpfl, '-', 'Color', cgreen);
    p6 = pflplot(MTV(indz),  zpfl, '-', 'Color', cyellow);
    p7 = pflplot(MD(indz),   zpfl, '-', 'Color', cmagenta);
    ylim([-1, 0]);
    tkex = [50, 80];
    xlim([-tkex(iwave), tkex(iwave)]);
    pflref(zpfl);
    ylabel('$z/h_b$','Interpreter','latex');
    xlabel('TKE budget','Interpreter','latex');
    lg = legend([p1, p2, p3, p4, p5, p6, p7],...
        '${\cal P}$', '${\cal P}^S$',...
        '${\cal B}$', '${\cal T}$', '${\cal R}$',...
        '${\cal T}^v$', '${\cal D}$');
    lg.Location = 'SouthWest';
    lg.Interpreter = 'latex';

    if l_save_fig
        figname = [outDir '/tke_budget.fig'];
        saveas(gcf, figname, 'fig');
        postProcessFig(figname, 5, 5);
    end

%     % check if the calculation is consistent
%     SPTKE = 0.5.*(SP11+SP22+SP33)./wgt;
%     StPTKE = 0.5.*(StP11+StP22+StP33)./wgt;
%     BTKE = 0.5.*(B11+B22+B33)./wgt;
%     TTKE = 0.5.*(T11+T22+T33)./wgt;
%     PSTKE = 0.5.*(PS11+PS22+PS33)./wgt;
%     TVTKE = 0.5.*(TV11+TV22+TV33)./wgt;
%     CTKE = 0.5.*(C11+C22+C33)./wgt;
%     DTKE = 0.5.*(D11+D22+D33)./wgt;
%     pflplot(SPTKE(indz), zpfl, '--b');
%     pflplot(StPTKE(indz), zpfl, '--r');
%     pflplot(BTKE(indz), zpfl, '--k');
%     pflplot(TTKE(indz), zpfl, '--', 'Color', ccyan);
%     pflplot(PSTKE(indz), zpfl, '--', 'Color', cgreen);
%     pflplot(TVTKE(indz), zpfl, '--', 'Color', cyellow);
%     pflplot(CTKE(indz), zpfl, '--', 'Color', cgray);
%     pflplot(DTKE(indz), zpfl, '--', 'Color', cmagenta);
end

%% figure 7: Reynolds stress budget
if l_Reynolds
    wgt = utau.^3./hb;
    zpfl = zu(1:nnzmax)./hb;
    rsx11 = [30, 30];
    rsx22 = [40, 50];
    rsx33 = [50, 100];
    rsx12 = [30, 30];
    rsx13 = [80, 150];
    rsx23 = [80, 100];
    fig = figure;
    fig.Units = 'inches';
    fig.Position = [1 1 12 8];
    subplot(2, 3, 1)
    xlims = [-rsx11(iwave),rsx11(iwave)];
    xlabel_str = '$\overline{u''_1 u''_1}$ budget';
    plot_reynolds_stress(SP11, StP11, B11, T11, PS11, C11, TV11, D11,...
        zpfl, ind_hb, wgt, xlims, xlabel_str);
    subplot(2, 3, 2)
    xlims = [-rsx22(iwave),rsx22(iwave)];
    xlabel_str = '$\overline{u''_2 u''_2}$ budget';
    plot_reynolds_stress(SP22, StP22, B22, T22, PS22, C22, TV22, D22,...
        zpfl, ind_hb, wgt, xlims, xlabel_str);
    subplot(2, 3, 3)
    xlims = [-rsx33(iwave),rsx33(iwave)];
    xlabel_str = '$\overline{u''_3 u''_3}$ budget';
    plot_reynolds_stress(SP33, StP33, B33, T33, PS33, C33, TV33, D33,...
        zpfl, ind_hb, wgt, xlims, xlabel_str);
    subplot(2, 3, 4)
    xlims = [-rsx12(iwave),rsx12(iwave)];
    xlabel_str = '$\overline{u''_1 u''_2}$ budget';
    plot_reynolds_stress(SP12, StP12, B12, T12, PS12, C12, TV12, D12,...
        zpfl, ind_hb, wgt, xlims, xlabel_str);
    subplot(2, 3, 5)
    xlims = [-rsx13(iwave),rsx13(iwave)];
    xlabel_str = '$\overline{u''_1 u''_3}$ budget';
    p=plot_reynolds_stress(SP13, StP13, B13, T13, PS13, C13, TV13, D13,...
        zpfl, ind_hb, wgt, xlims, xlabel_str);
        lg = legend(p,...
        '${\cal P}_{ij}$', '${\cal P}^S_{ij}$',...
        '${\cal B}_{ij}$', '${\cal T}_{ij}$',...
        '${\cal R}_{ij}$','${\cal C}_{ij}$',...
        '${\cal T}^v_{ij}$', '${\cal D}_{ij}$');
    lg.Location = 'SouthEast';
    lg.Interpreter = 'latex';
    subplot(2, 3, 6)
    xlims = [-rsx23(iwave),rsx23(iwave)];
    xlabel_str = '$\overline{u''_2 u''_3}$ budget';
    plot_reynolds_stress(SP23, StP23, B23, T23, PS23, C23, TV23, D23,...
        zpfl, ind_hb, wgt, xlims, xlabel_str);

    set(findall(gcf,'-property','FontSize'),'FontSize',12)
    if l_save_fig
        figname = [outDir '/ReynoldsStress_budget.fig'];
        saveas(gcf, figname, 'fig');
        [dir, name, ~] = fileparts(figname);
        print('-depsc2',[dir '/' name]);
    end
end

%% functions

function p = pflplot(xx, zz, varargin)
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

function rl = pflref(zz)
    nz = numel(zz);
    qnz = floor(nz/4);
    tqnz = floor(nz*3/4);
    rl = refline(0, zz(qnz));
    rl.Color = 'k';
    rl.LineStyle = '--';
    rl.LineWidth = 1.0;
    rl2 = refline(0, zz(tqnz));
    rl2.Color = 'k';
    rl2.LineStyle = '--';
    rl2.LineWidth = 1.0;
end

function p = plot_reynolds_stress(SP, StP, B, T, PS, C, TV, D,...
    z, ihb, wgt, xlims, xlabel_str)
    global ccyan cgreen cyellow cgray cmagenta
    inds = 2:ihb;
    cgray = [0.5 0.5 0.5];
    p = gobjects(1,8);
    p(1)=pflplot(SP(inds)./wgt, z(inds), '-b');
    p(2)=pflplot(StP(inds)./wgt, z(inds), '-r');
    p(3)=pflplot(B(inds)./wgt, z(inds), '-k');
    p(4)=pflplot(T(inds)./wgt, z(inds), '-', 'Color', ccyan);
    p(5)=pflplot(PS(inds)./wgt, z(inds), '-', 'Color', cgreen);
    p(6)=pflplot(C(inds)./wgt, z(inds), '-', 'Color', cgray);
    p(7)=pflplot(TV(inds)./wgt, z(inds), '-', 'Color', cyellow);
    p(8)=pflplot(D(inds)./wgt, z(inds), '-', 'Color', cmagenta);
    ylim([-1, 0]);
    ylabel('$z/h_b$','Interpreter','latex');
    xlim(xlims);
    xlabel(xlabel_str,'Interpreter','latex');
    pflref(z(inds));
end


function [p, pc, h] = plot_isosurface(xx, yy, zz, dat, isovalue,...
            angle_view, angle_cam)
% plot isosurface of dat

    fcolor = [111, 206, 249]./256;
    p = patch(isosurface(xx,yy,zz,dat,isovalue));
    p.FaceColor = fcolor;
    p.EdgeColor = 'none';
    p.AmbientStrength = 0.4;
    isonormals(xx,yy,zz,dat,p);
    pc = patch(isocaps(xx,yy,zz,dat,isovalue));
    pc.FaceColor = 'interp';
    pc.EdgeColor = 'none';
    axis tight
    daspect([1 1 1]);
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');
    view(angle_view(1), angle_view(2));
    h = camlight(angle_cam(1), angle_cam(2));
    lighting gouraud
end

function p = plot_RQmap(xdata, ydata, wgt)
% plot (R, Q) map

    xdata = xdata./wgt.^1.5; % normalization
    ydata = ydata./wgt;
    xpltmin = -1.5;
    xpltmax = 1.5;
    ypltmin = -3;
    ypltmax = 3;
    xmin = min(prctile(xdata(:), 0.1), xpltmin);
    xmax = max(prctile(xdata(:), 99.9), xpltmax);
    ymin = min(prctile(ydata(:), 0.1), ypltmin);
    ymax = max(prctile(ydata(:), 99.9), ypltmax);
    [hst, xi, yi] = jointDist(xdata, xmin, xmax, ydata, ymin, ymax);
    plot_dist_4p(hst,xi,yi);
    hold on;
    xlims = [xpltmin xpltmax];
    ylims = [ypltmin ypltmax];
    dyref = (0-ylims(1))/100;
    yref = ylims(1):dyref:0;
    xref = 2.*sqrt(3)./9.*(-yref).^(3/2);
    plot(xref, yref, '--k', 'LineWidth', 1.5);
    plot(-xref, yref, '--k', 'LineWidth', 1.5);
    plot([0, 0], ylims, '--k', 'LineWidth', 1.5);
    pbaspect([1 1 1]);
    xlim(xlims);
    ylim(ylims);
    xlabel_str = '$R/\langle S_{ij}S_{ij}\rangle^{3/2}$';
    ylabel_str = '$Q/\langle S_{ij}S_{ij}\rangle$';
    xlabel(xlabel_str, 'Interpreter', 'latex');
    ylabel(ylabel_str, 'Interpreter', 'latex');
    p = gcf;
end

function p = plot_RSQSmap(xdata, ydata, wgt)
% plot (R_S, Q_S) map

    xdata = xdata./wgt.^1.5; % normalization
    ydata = ydata./wgt;
    xpltmin = -1;
    xpltmax = 1;
    ypltmin = -2;
    ypltmax = 0;
    xmin = min(prctile(xdata(:), 0.1), xpltmin);
    xmax = max(prctile(xdata(:), 99.9), xpltmax);
    ymin = min(prctile(ydata(:), 0.1), ypltmin);
    ymax = ypltmax;
    [hst, xi, yi] = jointDist(xdata, xmin, xmax, ydata, ymin, ymax);
    plot_dist_4p(hst,xi,yi);
    hold on;
    xlims = [xpltmin xpltmax];
    ylims = [ypltmin ypltmax];
    dyref = (0-ylims(1))/100;
    yref = ylims(1):dyref:0;
    r = [-0.5, 0, 0.5, 1];
    for i=1:numel(r)
        a = r(i);
        xref = (-yref).^(3/2).*a.*(1+a).*(1+a+a.^2).^(-3/2);
        plot(xref, yref, '--k', 'LineWidth', 1.5);
    end
    pbaspect([1 1 1]);
    xlim(xlims);
    ylim(ylims);
    xlabel_str = '$R_S/\langle S_{ij}S_{ij}\rangle^{3/2}$';
    ylabel_str = '$Q_S/\langle S_{ij}S_{ij}\rangle$';
    xlabel(xlabel_str, 'Interpreter', 'latex');
    ylabel(ylabel_str, 'Interpreter', 'latex');
    p = gcf;
end

function p = plot_QAQSmap(xdata, ydata, wgt)
% plot (Q_A, -Q_S) map

    xdata = xdata./wgt; % normalization
    ydata = ydata./wgt;
    xpltmin = 0;
    xpltmax = 3;
    ypltmin = 0;
    ypltmax = 3;
    xmin = xpltmin;
    xmax = max(prctile(xdata(:), 99.9),xpltmax);
    ymin = ypltmin;
    ymax = max(prctile(ydata(:), 99.9),ypltmax);
    [hst, xi, yi] = jointDist(xdata, xmin, xmax, ydata, ymin, ymax);
    plot_dist_4p(hst,xi,yi);
    hold on;
    xlims = [xpltmin xpltmax];
    ylims = [ypltmin ypltmax];
    xmax = max(xlims(2),ylims(2));
    plot([0, xmax], [0, xmax], '--k', 'LineWidth', 1.5);
    pbaspect([1 1 1]);
    xlim([0, xmax]);
    ylim([0, xmax]);
    xlabel_str = '$Q_A/\langle S_{ij}S_{ij}\rangle$';
    ylabel_str = '$-Q_S/\langle S_{ij}S_{ij}\rangle$';
    xlabel(xlabel_str, 'Interpreter', 'latex');
    ylabel(ylabel_str, 'Interpreter', 'latex');
    p = gcf;
end

function h = plot_dist_4p(hst,xi,yi)
% plot_dist_4p plot the highest 50%, 75%, 90% and 95%
%   centered distribution

    % find the isolines of pdf that enclose the area in which the
    % total probability is 50%, 75%, 90% and 95%
    hsum = sum(hst(:));
    hlist = sort(hst(:),'descend')./hsum;
    hcum = cumsum(hlist);
    vl = [0.5, 0.75, 0.9, 0.95];
    nv = numel(vl);
    vlev = zeros(1,nv);
    for i=1:nv
        [~,ind] = min(abs(hcum-vl(i)));
        vlev(i) = hlist(ind);
    end
    pdfData = hst./hsum;
    % plot log10(pdfData) to make use of the full colarbar
    pdfData(pdfData==0) = 1e-12;    % avoid -inf for log10(0)
    [~,h] = contourf(xi,yi,log10(pdfData'));
    h.LevelListMode = 'manual';
    h.LevelList = log10(vlev);
    h.ShowText = 'on';
    h.TextList = vl;

    caxis([log10(vlev(end)) log10(vlev(1))]);
    tmp = colormap;
    inds = [2, 2, 15, 15, 15, 15, 15, 15, 35, 35, 35, 35, 35, 35, 60, 60];
    my_colorm = tmp(inds,:);
    colormap(my_colorm);
end

function save_figure_RQ(p, figname, l_save_fig)
% save figure

    if l_save_fig
        saveas(p, figname, 'fig');
        postProcessFig(figname, 6, 4);
    end
end

function res = get_ddx(dat, dx)
% get_ddx returns the derivative of dat with respect to x using
% central differences
% dat(nnx, nny, nnz) is a three dimensional array
% periodical boundary for x and y

    ndat = size(dat);
    datp1 = zeros(ndat);
    datm1 = zeros(ndat);
    datp1(1:end-1,:,:) = dat(2:end,:,:);
    datp1(end,:,:) = dat(1,:,:);
    datm1(2:end,:,:) = dat(1:end-1,:,:);
    datm1(1,:,:) = dat(end,:,:);
    res = 0.5./dx.*(datp1-datm1);
end

function res = get_ddy(dat, dy)
% get_ddy returns the derivative of dat with respect to y using
% central differences
% dat(nnx, nny, nnz) is a three dimensional array
% periodical boundary for x and y

    ndat = size(dat);
    datp1 = zeros(ndat);
    datm1 = zeros(ndat);
    datp1(:,1:end-1,:) = dat(:,2:end,:);
    datp1(:,end,:) = dat(:,1,:);
    datm1(:,2:end,:) = dat(:,1:end-1,:);
    datm1(:,1,:) = dat(:,end,:);
    res = 0.5./dy.*(datp1-datm1);
end

function res = get_ddz_w(dat, dz)
% get_ddz returns the derivative of dat with respect to z on u level
% dat(nnx, nny, nnz) is a three dimensional array
% dat on w level

    ndat = size(dat);
    datp1 = dat;
    datm1 = zeros(ndat);
    datm1(:,:,2:end) = dat(:,:,1:end-1);  % zero at the surface
    res = (datp1-datm1)./dz;
end

function res = get_ddz_u(dat, dz, varargin)
% get_ddz returns the derivative of dat with respect to z on u level
% dat(nnx, nny, nnz) is a three dimensional array
% dat on u level

    nVar = numel(varargin);
    if nVar==0
        l_level = 1;    % 1: w-level, 0: u-level, w-level by default
    elseif nVar==1
        if strcmp(varargin{1},'w')
            l_level = 1;
        elseif strcmp(varargin{1},'u')
            l_level = 0;
        else
            error('Input "w" for w-level, "u" for u-level.');
        end
    else
        error('Require 2 or 3 arguments, got %i.',nVar+2);
    end
    ndat = size(dat);
    datp1 = zeros(ndat);
    datp1(:,:,1:end-1) = dat(:,:,2:end);
    datp1(:,:,end) = dat(:,:,end);
    tmp = (datp1-dat)./dz;  % ddz on w level
    if l_level
        res = tmp;
    else
        tmp2 = zeros(ndat);
        tmp2(:,:,2:end) = tmp(:,:,1:end-1);
        res = 0.5.*(tmp+tmp2); % ddz on u level
    end
end
