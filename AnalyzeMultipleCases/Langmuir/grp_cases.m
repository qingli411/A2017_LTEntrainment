% case list for group Langmuir
% this group contains simulations under weak wind U_{10} = 5 m s^-1
% and surface cooling Q_0 = -5 W m^-2, Stokes drift calculated from 
% DHH wind wave spectrum and also monochromatic wave

caselist = [...
    'R7_BF05WD05WV00_ST00_ens01';...
    'R7_BF05WD05WV01_ST00_ens01';...
    'R7_BF05WD05WV02_ST00_ens01';...
    'R7_BF05WD05WV03_ST00_ens01';...
    'R7_BF05WD05WV04_ST00_ens01';...
    'R7_BF05WD05WV12_ST00_ens01';...
    'R7_BF05WD05WV13_ST00_ens01';...
    'R7_BF05WD05WV14_ST00_ens01';...
            ];

% plot
lcolor = {'-k','-b','--b','-.b',':b','--c','-.c',':c'};
llabel = {'No Wave',...
          'DHH, $C_\mathrm{p}/U_{10} = 1.2$',...
          'DHH, $C_\mathrm{p}/U_{10} = 1.0$',...
          'DHH, $C_\mathrm{p}/U_{10} = 0.8$',...
          'DHH, $C_\mathrm{p}/U_{10} = 0.6$',...
          'MC, $k = 0.10$ m$^{-1}$, $A = 0.80$ m',...
          'MC, $k = 0.10$ m$^{-1}$, $A = 0.53$ m',...
          'MC, $k = 0.10$ m$^{-1}$, $A = 0.40$ m',...
          };
      
lmarker = {'o','s','d','^','v','o','o','o'};
lmcolor = {'k','b','b','b','b','c','c','c'};
l_filled = {0,0,0,0,0,1,1,1};
