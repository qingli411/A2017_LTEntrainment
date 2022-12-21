close all; clear variables;

l_save_fig = 1;

path(path,'../share');
get_dataRootDir;

casename1 = 'R8_BF05WD10WV00_ST01_ens02';
casename2 = 'R8_BF05WD10WV11_ST01_ens02';
tttttt = '079201';
outDir1 = [outRootDir '/volume/' casename1 '/' tttttt];
outDir2 = [outRootDir '/volume/' casename2 '/' tttttt];

dat1 = load([outDir1 '/wpsdata.mat']);
dat2 = load([outDir2 '/wpsdata.mat']);

f1 = initncarPflData(dataRootDir, casename1);
f2 = initncarPflData(dataRootDir, casename2);

tmp = f1.getVar('utau');
utau1 = tmp(end);
tmp = f2.getVar('utau');
utau2 = tmp(end);

fig=figure;
fig.Units = 'inches';
fig.Position = [1 1 5 7.5];
hold on;
dat = [dat1, dat2];
wgt = [utau1.^2, utau2.^2];
linestyl = {'-', '--'};
markerstyl = {'*','o'};
for i=1:2
    wps = dat(i).wps;
    ups = dat(i).ups;
    vps = dat(i).vps;
    zw = dat(i).zw;
    mi0 = dat(i).indhb;
    ind1 = dat(i).ind1;
    ind2 = dat(i).ind2;
    ind3 = dat(i).ind3;
    ind4 = dat(i).ind4;
    plot(wps./wgt(i), -zw./zw(mi0), [linestyl{i} 'r'], 'LineWidth', 1.5);
    plot(wps([ind1, ind2, ind3, ind4])./wgt(i),...
       -[zw(ind1), zw(ind2), zw(ind3), zw(ind4)]./zw(mi0),...
       [markerstyl{i} 'r'], 'LineWidth', 1.5);
    plot(ups./wgt(i), -zw./zw(mi0), [linestyl{i} 'k'], 'LineWidth', 1.5);
    plot(ups([ind1, ind2, ind3, ind4])./wgt(i),...
       -[zw(ind1), zw(ind2), zw(ind3), zw(ind4)]./zw(mi0),...
       [markerstyl{i} 'k'], 'LineWidth', 1.5);
    plot(vps./wgt(i), -zw./zw(mi0), [linestyl{i} 'b'], 'LineWidth', 1.5);
    plot(vps([ind1, ind2, ind3, ind4])./wgt(i),...
       -[zw(ind1), zw(ind2), zw(ind3), zw(ind4)]./zw(mi0),...
       [markerstyl{i} 'b'], 'LineWidth', 1.5);
end
xlabel('$\overline{u_i''^2}/u_*^2$','Interpreter', 'latex');
ylabel('$z/h_\mathrm{b}$', 'Interpreter', 'latex');
xlim([0, 5]);
ylim([-1.0, 0]);
pbaspect([1, 2, 1]);

a = findall(gcf,'-property','FontSize');
for i=1:numel(a)
    if strcmp(a(i).Type,'axes')
        a(i).FontSize = 12;
    elseif strcmp(a(i).Type,'text')
        a(i).FontSize = 14;
    else
        a(i).FontSize = 12;
    end
end

if l_save_fig
    figname = [outDir1 '/profile_wps_all.fig'];
    saveas(gcf, figname, 'fig');
    [dir, name, ~] = fileparts(figname);
    print('-depsc2',[dir '/' name]);
end
