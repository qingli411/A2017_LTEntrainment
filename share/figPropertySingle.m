% set figure properties
wsize = 6;
hsize = 4;
mleft = 0.25;
mdown = 0.25;
mright = wsize-0.25;
mup = hsize-0.25;
set(gcf,'Units','inches');
set(gcf,'PaperUnits','inches'); 
set(gcf,'PaperPositionMode', 'manual');
paperposition = [mleft mdown mright mup];
set(gcf,'PaperPosition',paperposition);
size = [wsize hsize];
set(gcf,'PaperSize',size);
position = [mleft mdown mright mup];
set(gcf,'Position',position);
% Change renderer, 'painters' works well for vector formats while
% 'opengl' works well for bitmap format
set(gcf, 'renderer', 'painters');
