[s,r] = system('hostname | tr -d ''\n''');
if s==0
    if strcmp(r,'gc134-chace.geo.brown.edu')
%         dataRootDir = '/Users/qingli/data_local/archive_les';
        dataRootDir = '/Volumes/Qing_Work/data/NCARLES/archive_les';
        outRootDir = '/Users/qingli/work/ncarles/results';
    elseif strcmp(r,'Turen.local') || strcmp(r,'Turen.home')...
            || strcmp(r,'turen.devices.brown.edu')
        dataRootDir = '/Volumes/scratch/data/NCARLES/archive_les';
        outRootDir = '/Users/qingli/work/ncarles/results';
    else
        error('Hostname %s not recognized...',r);
    end
else
    error('Hostname error...')
end
