function f = initncarPflData(dataRootDir, casename)

    % set up profile data parameters
    dataDir = [dataRootDir '/hist/'];
    fPrefix = 'his.mp.vis';
    [s0, r0] = system(['ls ' dataDir '/' casename ]);
    disp(['Reading profile data: ' r0(1:27)]);
    if s0==0
        istart_str = r0(12:17);
        iend_str = r0(19:24);
        tStart = str2double(istart_str);
        tEnd = str2double(iend_str);
    end
    % creating object ncarlesPflData
    f = ncarlesPflData('caseName', casename,...
                       'dataDir',  dataDir,...
                       'fPrefix',  fPrefix,...
                       'tStart',   tStart,...
                       'tEnd',     tEnd);
end