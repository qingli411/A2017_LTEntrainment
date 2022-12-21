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