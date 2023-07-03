function  createFigDLvsConventionalFitting(FFmaps, dlMaps, errormaps, dlErrormaps)
% function  createFigDLvsConventionalFitting(FFmaps, dlMaps, errormaps, dlErrormaps)

%Helper function for visualising performance of DL methods over a range of
%FF values, compared to conventional fitting

%Written by Dr Tim Bray
%t.bray@ucl.ac.uk
%November 2022

%Inputs
% FFmaps is structure containing FFmaps from conventional fittings
% dlMaps is structure containing FFmaps from DL fittings
% Errormaps is structure containing errormaps from conventional fittings
% dlErrormaps is structure containing errormaps from DL fittings

figure('Name', 'MAGORINO parameter error')
subplot(2,3,1)
image(FFmaps.standard,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 1];
FigLabels;
title('Gaussian magnitude FF maps')
colorbar

subplot(2,3,2)
image(FFmaps.Rician,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 1];
FigLabels;
title('Rician magnitude FF maps')
colorbar

subplot(2,3,3)
image(dlMaps.ff,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 1];
FigLabels;
title('DL FF maps')
colorbar

subplot(2,3,4)
image(errormaps.FFstandard,'CDataMapping','scaled')
ax=gca;
ax.CLim=[-.2 .2];
FigLabels;
title('Gaussian magnitude FF error')
colorbar

subplot(2,3,5)
image(errormaps.FFRician,'CDataMapping','scaled')
ax=gca;
ax.CLim=[-.2 .2];
FigLabels;
title('Rician magnitude FF error')
colorbar

subplot(2,3,6)
image(dlErrormaps.ff,'CDataMapping','scaled')
ax=gca;
ax.CLim=[-.2 .2];
FigLabels;
title('DL error maps')
colorbar

%1.5 Print error
sum(abs(errormaps.FFstandard),'all')
sum(abs(errormaps.FFRician),'all')
sum(abs(dlErrormaps.ff),'all')