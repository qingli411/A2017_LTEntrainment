close all; clear variables;
% This script plots the parameter space in CESM estimated from
% the daily averaged output of a typical simulated year in a
% fully coupled CESM-WW3 simulation (PI-WW3 in Li et al., 2017)

% load data
dat = load('LF17_Data_ParSpaceCESM.mat');
laturb = dat.laturb;
lasl = dat.lasl;
bflux = dat.bflux;
hb = dat.hb;
ustar = dat.ustar;

%% Figure 1: w*^3/u*^2uS - LaTurb regime diagram
xdata  = laturb;
ydata = bflux.*hb./ustar.^3.*laturb.^2;
figure;
tmp = isnan(ydata);
nnan = sum(tmp(:)); % number of NaNs
inds = find(ydata<0);
tmp = ydata<0;
nneg = sum(tmp(:)); % number of negtive values
ntot = numel(ydata); % total number
ydata(inds) = NaN;
xdata(inds) = NaN;
xdata = log10(xdata);
ydata = log10(ydata);
rpos = (ntot-nnan-nneg)./(ntot-nnan);
fprintf('The ratio of used points: %f\n', rpos);

xpltmin = -1;
xpltmax = 1;
ypltmin = -3;
ypltmax = 3;
xmin = min(prctile(xdata(:), 0.1), xpltmin);
xmax = max(prctile(xdata(:), 99.9), xpltmax);
ymin = min(prctile(ydata(:), 0.1), ypltmin);
ymax = max(prctile(ydata(:), 99.9), ypltmax);
[hst, xi, yi] = jointDist(xdata, xmin, xmax, ydata, ymin, ymax);
xi = 10.^xi;
yi = 10.^yi;
plot_dist_4p(hst,xi,yi);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
pbaspect([1, 1, 1]);
ystr = '$h_\mathrm{b}/L_\mathrm{L}$';
xstr = '$La_\mathrm{t}$';
xlabel(xstr,'Interpreter','latex');
ylabel(ystr,'Interpreter','latex');
xlim([1e-1, 1e1]);
ylim([1e-3, 1e3]);

save('wstar3ustar2uS_LaTurb.mat','hst','xi','yi');


%% Figure 2: w*^3/u*^3 - LaSL regime diagram
xdata = lasl;
ydata = bflux.*hb./ustar.^3;
figure;
tmp = isnan(ydata);
nnan = sum(tmp(:)); % number of NaNs
inds = find(ydata<0);
tmp = ydata<0;
nneg = sum(tmp(:)); % number of negtive values
ntot = numel(ydata); % total number
ydata(inds) = NaN;
xdata(inds) = NaN;
xdata = log10(xdata);
ydata = log10(ydata);
rpos = (ntot-nnan-nneg)./(ntot-nnan);
fprintf('The ratio of used points: %f\n', rpos);

xpltmin = -1;
xpltmax = 1;
ypltmin = -2;
ypltmax = 2;
xmin = min(prctile(xdata(:), 0.1), xpltmin);
xmax = max(prctile(xdata(:), 99.9), xpltmax);
ymin = min(prctile(ydata(:), 0.1), ypltmin);
ymax = max(prctile(ydata(:), 99.9), ypltmax);
[hst, xi, yi] = jointDist(xdata, xmin, xmax, ydata, ymin, ymax);
xi = 10.^xi;
yi = 10.^yi;
plot_dist_4p(hst,xi,yi);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
pbaspect([1, 1, 1]);
ystr = '$-h_\mathrm{b}/(\kappa L)$';
xstr = '$La_\mathrm{SL}$';
xlabel(xstr,'Interpreter','latex');
ylabel(ystr,'Interpreter','latex');
xlim([1e-1, 2e1]);
ylim([1e-2, 1e2]);

save('wstar3ustar3_LaSL.mat','hst','xi','yi');


%% functions
function [hst, xi, yi] = jointDist(dat1,vmin1,vmax1,dat2,vmin2,vmax2)
% jointDist joint distribution of two variabels
%   [hst, xi, yi] = jointDist(dat1, vmin1, vmax1, dat2, vmin2, vmax2 )
%   returns the joint distribution (number of points in each bins) of
%   dat1 and dat2. 100 bins for each are used, with the range set by 
%   [vmin1, vmax1] and [vmin2, vmax2], respectively. The edges of each
%   bin are returned in xi and yi.
%
%   See also hist3

    newsize = numel(dat1);
    data = zeros(2,newsize);
    data(1,:) = dat1(:);
    data(2,:) = dat2(:);
    xi = linspace(vmin1,vmax1,100);
    yi = linspace(vmin2,vmax2,100);
    hst = hist3(data','Edge',{xi',yi'});
end

function h = plot_dist_4p(hst,xi,yi)
% plot_dist_4p plot the highest 30%, 60%, 90% and 99%
%   centered distribution

    % find the isolines of pdf that enclose the area in which the
    % total probability is 50%, 75%, 90% and 95%
    hsum = sum(hst(:));
    hlist = sort(hst(:),'descend')./hsum;
    hcum = cumsum(hlist);
%     vl = [0.5, 0.75, 0.9, 0.95];
    vl = [0.3, 0.6, 0.9, 0.99];
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
    
    clim([log10(vlev(end)) log10(vlev(1))]);
    tmp = colormap(parula(64));
%     % [0.5, 0.75, 0.9, 0.95]
%     inds = [2, 2, 15, 15, 15, 15, 15, 15, 35, 35, 35, 35, 35, 35, 60, 60];
    % [0.3, 0.6, 0.9, 0.95]
    inds = [2, 2, 15, 15, 15, 15, 15, 15, 15, 35, 35, 35, 35, 35, 60, 60];
    my_colorm = tmp(inds,:);
    colormap(my_colorm);
end
