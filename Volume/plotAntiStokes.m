close all; clear variables;
% load shared scripts
path(path,'../share');
path(path,'../share/cbrewer')

casename = 'R8_BF05WD05WV00_fL00_ens01';
casename2 = 'R8_BF05WD05WV12_fL00_ini01';

get_dataRootDir;    % get dataRootDir

outDir = [ outRootDir '/profiles'];

f1 = init_ncarlesPflData(dataRootDir, casename);
f2 = init_ncarlesPflData(dataRootDir, casename2);

% pfl1 = f1.getProfiles('inertial',[0,1]);
% pfl2 = f2.getProfiles('inertial',[0,1]);

tint = 61095;
pfl1 = f1.getProfiles('second',[f1.time(end)-tint, f1.time(end)]);
pfl2 = f2.getProfiles('second',[f2.time(end)-tint, f2.time(end)]);

tmp = f1.getVar('utau');
utau1 = tmp(end);
tmp = f2.getVar('utau');
utau2 = tmp(end);

figure;
p1=plot(pfl1.uxym./utau1, pfl1.z_u./pfl1.mean_hb,...
    '-k', 'LineWidth', 1.5);
hold on;
p2=plot(pfl2.uxym./utau2, pfl2.z_u./pfl2.mean_hb,...
    '-r', 'LineWidth', 1.5);
p3=plot(pfl2.stokes./utau2, pfl2.z_u./pfl2.mean_hb,...
    '--r', 'LineWidth', 1.5);
p4=plot((pfl2.uxym+pfl2.stokes)./utau2, pfl2.z_u./pfl2.mean_hb,...
    '-.r', 'LineWidth', 1.5);
ylim([-1,0]);
xlabel('$\overline{u}/u^*$', 'Interpreter', 'latex');
ylabel('$z/h_\mathrm{b}$', 'Interpreter', 'latex');
lg = legend([p1,p2,p3,p4],...
    'ST-Eulerian','LT-Eulerian','LT-StokesDrift','LT-Lagrangian');
lg.Location = 'SouthEast';
lg.Interpreter = 'latex';

figname = [outDir '/uxymPfl_STLTfL00IC.fig'];
saveas(gcf, figname, 'fig');
postProcessFig(figname, 4, 4);

figure;
pb1 = plot(pfl1.vxym./utau1, pfl1.z_u./pfl1.mean_hb,...
    '-k', 'LineWidth', 1.5);
hold on;
pb2 = plot(pfl2.vxym./utau2, pfl2.z_u./pfl2.mean_hb,...
    '-r', 'LineWidth', 1.5);
ylim([-1,0]);
xlabel('$\overline{v}/u^*$', 'Interpreter', 'latex');
ylabel('$z/h_\mathrm{b}$', 'Interpreter', 'latex');
lg = legend([pb1,pb2],...
    'ST-Eulerian','LT-Eulerian');
lg.Location = 'SouthWest';
lg.Interpreter = 'latex';

figname = [outDir '/vxymPfl_STLTfL00IC.fig'];
saveas(gcf, figname, 'fig');
postProcessFig(figname, 4, 4);

function f = init_ncarlesPflData(dataRootDir, casename)

    % set up profile data parameters
    dataDir = [dataRootDir '/hist'];
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
