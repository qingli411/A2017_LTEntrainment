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
nc = numel(namelist);

h0 = 42;    % initial mixed layer depth
izi = [2, 6, 40, 63, 64, 65, 66, 67, 68];   % indices of saved x-y slices
% izi3 = [2, 6, 44, 77, 78, 79, 80, 81, 82];
% izi4 = [2, 6, 41, 68, 69, 70, 71, 72, 73];
idxentr = [5, 7, 7, 5, 7];
l_save_fig = 1;
c_gray = [0.5, 0.5, 0.5];

% set up directories
get_dataRootDir;    % get dataRootDir and outRootDir
dataDir = [dataRootDir '/viz/'];

for k = 1:nc
    % filename
    casename = namelist{k};
    fprintf([casename '\n']);
    filename_xy = ['viz.vis.' startlist{k} '.' endlist{k} '.xy.nc'];
    outDir = [outRootDir '/slice/' casename '/' startlist{k} '-' endlist{k}];
    system(['mkdir -p ' outDir]);
    % input file name
    infile = [dataDir casename '/' filename_xy];
    % get w prime and b prime
    [wp, bp] = get_wp_bp(infile, idxentr(k));
    % weight
    wp_wgt = std(wp(:));
    bp_wgt = std(bp(:));
    
    if k == 1
        wp0_wgt = wp_wgt;
        bp0_wgt = bp_wgt;
    end
    wp_wgt/wp0_wgt
    bp_wgt/bp0_wgt
    
    % plot
    figure;
    xlims = [-5, 5];
    ylims = [-5, 5];
    [hst, xi, yi] = jointDist(wp(:)./wp_wgt, xlims(1), xlims(2), bp(:)./bp_wgt, ylims(1), ylims(2));
    plot_dist_4p(hst,xi,yi);
    xlabel('$u_3^\prime/\sigma(u_3^\prime)$', 'Interpreter', 'latex');
    ylabel('$b^\prime/\sigma(b^\prime)$', 'Interpreter', 'latex');
    xlim(xlims);
    ylim(ylims);
    rl = line(xlims,[0,0]);
    rl.Color = 'k';
    rl = line([0,0],ylims);
    rl.Color = 'k';
    pbaspect([1 1 1]);

    % output figure name
    figname = [outDir '/wbJointPDF.fig'];
    if l_save_fig
        saveas(gcf, figname, 'fig');
        postProcessFig(figname, 5, 5);
    end
end

function [wp, bp] = get_wp_bp(infile, idx0)
    alpha = 2.0000e-04;
    g = 9.81;
    % also read the buoyancy above an below the entrainment
    % layer for the flux limiter
    idx0m1 = idx0-1;
    idx0p1 = idx0+1;
    idx0p2 = idx0+2;
    tmp_w = double(ncread(infile,'w'));
    tmp_b = double(ncread(infile,'t').*alpha.*g);
    [nx, ny, ~, nt] = size(tmp_w);
    wp = zeros(nx, ny, nt);
    bp = zeros(nx, ny, nt);
    for i = 1:nt
        wxy0 = squeeze(tmp_w(:,:,idx0,i));
        bxy0 = squeeze(tmp_b(:,:,idx0,i));
        bxy0m1 = squeeze(tmp_b(:,:,idx0m1,i));
        bxy0p1 = squeeze(tmp_b(:,:,idx0p1,i));
        bxy0p2 = squeeze(tmp_b(:,:,idx0p2,i));
        % calculate the entrainment buoyancy fluc, apply flux limiter
        tmp = max(-wxy0, 0).*(bxy0+fluxLimiter(bxy0p1, bxy0, bxy0m1));
        tmp2 = min(-wxy0, 0).*(bxy0p1+fluxLimiter(bxy0, bxy0p1, bxy0p2));
        wpbp = -tmp-tmp2;
        wp(:,:,i) = wxy0;
        % implied buoyancy perturbation
        tmp = wpbp./wxy0;
        bp(:,:,i) = tmp-mean(tmp(:));
    end
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
    cmap = colormap;
    inds = [2, 2, 15, 15, 15, 15, 15, 15, 35, 35, 35, 35, 35, 35, 60, 60];
    my_colorm = cmap(inds,:);
    colormap(my_colorm);
    
    % percentage of four quadrants
    [nx,ny] = size(hst);
    tmp = hst(nx/2+1:nx,ny/2+1:ny);
    pqd1 = sum(tmp(:))/hsum;
    tmp = hst(1:nx/2,ny/2+1:ny);
    pqd2 = sum(tmp(:))/hsum;
    tmp = hst(1:nx/2,1:ny/2);
    pqd3 = sum(tmp(:))/hsum;
    tmp = hst(nx/2+1:nx,1:ny/2);
    pqd4 = sum(tmp(:))/hsum;
    text(-4, -4, sprintf('(III) %6.2f %%',pqd3*100));
    text(1, 4, sprintf('(I) %6.2f %%',pqd1*100));
    text(-4, 4, sprintf('(II) %6.2f %%',pqd2*100));
    text(1, -4, sprintf('(IV) %6.2f %%',pqd4*100));
end