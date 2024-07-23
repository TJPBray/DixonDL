function  createFigDLvsConventionalFitting(FFmaps, dlMaps, errormaps, dlErrormaps, dlSdMaps, sdMaps)
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

figure('Name', 'Parameter error for DL vs conventional fitting')
subplot(2,3,1)
image(FFmaps.standard,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 1];
FigLabels;
title('Gaussian FF')
h=colorbar
h.Label.String = "FF";
h.Label.FontSize = 12;

subplot(2,3,2)
image(FFmaps.Rician,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 1];
FigLabels;
title('Rician FF')
h=colorbar
h.Label.String = "FF";
h.Label.FontSize = 12;

subplot(2,3,3)
image(dlMaps.ff,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 1];
FigLabels;
title('DL FF')
h=colorbar
h.Label.String = "FF";
h.Label.FontSize = 12;

subplot(2,3,4)
image(errormaps.FFstandard,'CDataMapping','scaled')
ax=gca;
ax.CLim=[-.61 .61];
FigLabels;
title('Gaussian FF error')
h=colorbar
h.Label.String = "FF error";
h.Label.FontSize = 12;

subplot(2,3,5)
image(errormaps.FFRician,'CDataMapping','scaled')
ax=gca;
ax.CLim=[-.61 .61];
FigLabels;
title('Rician FF error')
h=colorbar
h.Label.String = "FF error";
h.Label.FontSize = 12;

subplot(2,3,6)
image(dlErrormaps.ff,'CDataMapping','scaled')
ax=gca;
ax.CLim=[-.61 .61];
FigLabels;
title(['DL PDFF error'])
h=colorbar
h.Label.String = "FF error";
h.Label.FontSize = 12;

%1.5 Print mean absolute error
sum(abs(errormaps.FFstandard),'all')
sum(abs(errormaps.FFRician),'all')
sum(abs(dlErrormaps.ff),'all')

%1.6 Compare SD

figure('Name', 'SD comparison')
subplot(2,2,1)
image(sdMaps.FFstandard,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 1];
FigLabels;
title('Conventional fitting FF SD')
h=colorbar
h.Label.String = "FF SD";
h.Label.FontSize = 12;

subplot(2,2,2)
image(dlSdMaps.ff,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 1];
FigLabels;
title('DL FF SD')
h=colorbar
h.Label.String = "FF SD";
h.Label.FontSize = 12;

subplot(2,2,3)
image(sdMaps.R2standard,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 .2];
FigLabels;
title('Conventional fitting R2* SD')
h=colorbar
h.Label.String = "R2* SD (s^-^1)";
h.Label.FontSize = 12;

subplot(2,2,4)
image(dlSdMaps.r2,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 .2];
FigLabels;
title('DL R2* SD')
h=colorbar
h.Label.String = "R2* SD (s^-^1)";
h.Label.FontSize = 12;

%1.5 Print mean absolute error
sum(abs(sdMaps.FFstandard),'all')
sum(abs(dlSdMaps.ff),'all')
sum(abs(sdMaps.R2standard),'all')
sum(abs(dlSdMaps.r2),'all')