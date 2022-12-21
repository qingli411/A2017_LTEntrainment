if l_update_data
    disp('Reading data...');
    % read data
    w_xy = ncread([dataDir casename '/' filename_xy],'w');
    w_xz = ncread([dataDir casename '/' filename_xz],'w');
    w_yz = ncread([dataDir casename '/' filename_yz],'w');
    u_xy = ncread([dataDir casename '/' filename_xy],'u');
    v_xy = ncread([dataDir casename '/' filename_xy],'v');
    t_xy = ncread([dataDir casename '/' filename_xy],'t');
    t_xz = ncread([dataDir casename '/' filename_xz],'t');
    t_yz = ncread([dataDir casename '/' filename_yz],'t');
    x = ncread([dataDir casename '/' filename_xy],'x');
    y = ncread([dataDir casename '/' filename_xy],'y');
    z = ncread([dataDir casename '/' filename_xz],'zw');
    time = ncread([dataDir casename '/' filename_xz],'time');
    time_mm = time./60;

    nsize = size(w_xy);
    nx = nsize(1);
    ny = nsize(2);
    nl = nsize(3);
    nt = nsize(4);
    nsize = size(w_xz);
    nz = nsize(2);

    % get ups, vps, wps, upvp, upwp, vpwp
    um_xy = squeeze(mean(mean(u_xy,1),2));
    vm_xy = squeeze(mean(mean(v_xy,1),2));
    wm_xy = squeeze(mean(mean(w_xy,1),2));

    up = u_xy;
    vp = v_xy;
    wp = w_xy;

    for i=1:nl
        for j=1:nt
            up(:,:,i,j) = u_xy(:,:,i,j)-um_xy(i,j);
            vp(:,:,i,j) = v_xy(:,:,i,j)-vm_xy(i,j);
            wp(:,:,i,j) = w_xy(:,:,i,j)-wm_xy(i,j);
        end
    end

    ups  = squeeze(mean(mean(up.*up,1),2));
    vps  = squeeze(mean(mean(vp.*vp,1),2));
    wps  = squeeze(mean(mean(wp.*wp,1),2));
    upvp = squeeze(mean(mean(up.*vp,1),2));
    upwp = squeeze(mean(mean(up.*wp,1),2));
    vpwp = squeeze(mean(mean(vp.*wp,1),2));
    clear up vp wp um_xy vm_xy u_xy v_xy;

    % swap data blocks to move xz and yz slices from the center to the side
    if l_swap
        tmp = w_xy;
        indx1 = 1:nx/2;
        indx2 = nx/2+1:nx;
        indy1 = 1:ny/2;
        indy2 = ny/2+1:ny;
        w_xy(indx1,indy1,:,:) = tmp(indx2,indy2,:,:);
        w_xy(indx2,indy1,:,:) = tmp(indx1,indy2,:,:);
        w_xy(indx1,indy2,:,:) = tmp(indx2,indy1,:,:);
        w_xy(indx2,indy2,:,:) = tmp(indx1,indy1,:,:);
        clear tmp;
        tmp = w_xz;
        w_xz(indx1,:,:,:) = tmp(indx2,:,:,:);
        w_xz(indx2,:,:,:) = tmp(indx1,:,:,:);
        clear tmp;
        tmp = w_yz;
        w_yz(indy1,:,:,:) = tmp(indy2,:,:,:);
        w_yz(indy2,:,:,:) = tmp(indy1,:,:,:);
        clear tmp;
    end

    save([dataDir casename '/slice.mat']);
else
    load([dataDir casename '/slice.mat']);
    disp('Loading workspace...');
end
