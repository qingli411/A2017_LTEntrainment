close all; clear variables;
% This script perform AnalyzeMean.m on all cases that satisfy the search
% condition
% Qing Li, 161126

[s, r] = system('cat caselist3.txt');
if s ~= 0
    error('ls error');
end

rs = strsplit(r);

nr = numel(rs);

for i=1:nr
    if ~isempty(rs{i})
        system(['sed "s#SET_CASENAME#''' char(rs{i})...
        '''#g" ./AnalyzeMean.m > AnalyzeMeanTMP.m']);
        AnalyzeMeanTMP;
        system('rm AnalyzeMeanTMP.m');
    end
end