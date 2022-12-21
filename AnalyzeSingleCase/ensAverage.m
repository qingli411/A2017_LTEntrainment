close all; clear all;
% This script do ensemble average on the profile data
% Qing Li, 170123

% add path for ncarleePflData
path(path,'../share');

% set case name
casename = 'R7_BF05WD05WV11_ST00';

% get ensemble number
[s, r] = system(['for f in ' casename '_ens[0-9][0-9]; do echo $f; done']);
if s ~= 0
    error('ls error');
end
rs = strsplit(r);
% remove empty cell
for i=1:numel(rs)
    if isempty(rs{i})
        rs(i) = [];
    end
end
nr = numel(rs);

% sum over all cases
dn = 'MeanProfile.mat';
navg = nr;
for i=1:nr
    ensnum = sprintf('%02d',i);
    cn = rs{i};
    fileName = [cn '/' dn];
    if exist(fileName,'file')
        fprintf('Case Name: %s; Loading...\n', cn);
        dat = load(fileName);
        % initialization
        if i==1
            fldnames = fieldnames(dat);
            nf = numel(fldnames);
            datnew = dat;
        else
            for j=1:nf
                if strcmp(fldnames{j},'stat_wb')
                    datnew.stat_wb.mean = datnew.stat_wb.mean...
                                        + dat.stat_wb.mean;
                    datnew.stat_wb.median = datnew.stat_wb.median...
                                        + dat.stat_wb.median;
                    datnew.stat_wb.p25 = datnew.stat_wb.p25...
                                        + dat.stat_wb.p25;
                    datnew.stat_wb.p75 = datnew.stat_wb.p75...
                                        + dat.stat_wb.p75;
                elseif ~strcmp(fldnames{j},'f')
                    datnew.(fldnames{j}) = datnew.(fldnames{j})...
                                         + dat.(fldnames{j});
                end
            end
        end
    else
        fprintf('Case Name: %s; data not exist, skipping...\n', cn);
        navg = navg - 1;
    end
end
% do average
for j=1:nf
    if strcmp(fldnames{j},'stat_wb')
        datnew.stat_wb.mean = datnew.stat_wb.mean./navg;
        datnew.stat_wb.median = datnew.stat_wb.median./navg;
        datnew.stat_wb.p25 = datnew.stat_wb.p25./navg;
        datnew.stat_wb.p75 = datnew.stat_wb.p75./navg;
    elseif ~strcmp(fldnames{j},'f')
        datnew.(fldnames{j}) = datnew.(fldnames{j})./navg;
    end
end
% save data
dirout = [casename '_ensmn'];
if ~exist(dirout,'dir')
    [s, r] = system(['mkdir ' dirout]);
    if s ~= 0
        error('mkdir error');
    end
end
save([dirout '/' dn], '-struct', 'datnew');