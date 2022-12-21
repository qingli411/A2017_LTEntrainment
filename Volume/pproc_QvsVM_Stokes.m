close all; clear variables;
% case name
casename = 'R8_BF05WD10WV11_ST00_ens01';
figs = {'Q_dom4',...
        'Q_dom4b',...
        'VM_dom4',...
        'VM_dom4b'};
nf = numel(figs);

wsize = 6;
hsize = 2.5;
mleft = 0.15;
mdown = 0.15;
mright = wsize-0.15;
mup = hsize-0.15;

for i=1:nf
    filename = [casename '/' figs{i} '.fig'];
    % load figure
    [dir, name, ~] = fileparts(filename);
    open(filename);

    % Change size
    set(gcf,'Units','inches');
    set(gcf,'PaperUnits','inches'); 
    set(gcf,'PaperPositionMode', 'manual');
    paperposition = [mleft mdown mright mup];
    set(gcf,'PaperPosition',paperposition);
    size = [wsize hsize];
    set(gcf,'PaperSize',size);
    position = [mleft mdown mright mup];
    set(gcf,'Position',position);
    set(findall(gcf,'-property','FontSize'),'FontSize',10)
    print('-depsc2',[dir '/' name]);
end