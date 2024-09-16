function  createFigDLvsConventionalFitting(FFmaps, R2maps, dlMaps, errormaps, dlErrormaps, dlSdMaps, sdMaps)
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


%% PDFF parameter error 

figure('Name', 'Parameter error for PDFF: DL vs conventional fitting')
subplot(3,3,1)
image(FFmaps.standard,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 1];
FigLabels;
title('Gaussian FF')
h=colorbar
h.Label.String = "FF";
h.Label.FontSize = 12;

subplot(3,3,2)
image(FFmaps.Rician,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 1];
FigLabels;
title('Rician FF')
h=colorbar
h.Label.String = "FF";
h.Label.FontSize = 12;

subplot(3,3,3)
image(dlMaps.ff,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 1];
FigLabels;
title('DL FF')
h=colorbar
h.Label.String = "FF";
h.Label.FontSize = 12;

subplot(3,3,4)
image(errormaps.FFstandard,'CDataMapping','scaled')
ax=gca;
ax.CLim=[-1 1];
FigLabels;
title('Gaussian FF error')
h=colorbar
h.Label.String = "FF error";
h.Label.FontSize = 12;

subplot(3,3,5)
image(errormaps.FFRician,'CDataMapping','scaled')
ax=gca;
ax.CLim=[-1 1];
FigLabels;
title('Rician FF error')
h=colorbar
h.Label.String = "FF error";
h.Label.FontSize = 12;

subplot(3,3,6)
image(dlErrormaps.ff,'CDataMapping','scaled')
ax=gca;
ax.CLim=[-1 1];
FigLabels;
title(['DL FF error'])
h=colorbar
h.Label.String = "FF error";
h.Label.FontSize = 12;

subplot(3,3,7)
image(sdMaps.FFstandard,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 1];
FigLabels;
title('Gaussian FF SD')
h=colorbar
h.Label.String = "FF SD";
h.Label.FontSize = 12;

subplot(3,3,8)
image(sdMaps.FFRician,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 1];
FigLabels;
title('Rician FF SD')
h=colorbar
h.Label.String = "FF SD";
h.Label.FontSize = 12;

subplot(3,3,9)
image(dlSdMaps.ff,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 1];
FigLabels;
title('DL FF SD')
h=colorbar
h.Label.String = "FF SD";
h.Label.FontSize = 12;

%Print mean absolute error
sum(abs(errormaps.FFstandard),'all')
sum(abs(errormaps.FFRician),'all')
sum(abs(dlErrormaps.ff),'all')

%% R2* parameter error

figure('Name', 'Parameter error for R_2^*: DL vs conventional fitting')

subplot(3,3,1)
image(1000*R2maps.standard,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 500];
FigLabels;
title('Gaussian R_2^*')
h=colorbar
h.Label.String = "R_2^* (s^-^1)";
h.Label.FontSize = 12;

subplot(3,3,2)
image(1000*R2maps.Rician,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 500];
FigLabels;
title('Rician R_2^*')
h=colorbar
h.Label.String = "R_2^* (s^-^1)";
h.Label.FontSize = 12;

subplot(3,3,3)
image(1000*dlMaps.r2,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 500];
FigLabels;
title('DL R_2^*')
h=colorbar
h.Label.String = "R_2^* (s^-^1)";
h.Label.FontSize = 12;

subplot(3,3,4)
image(1000*errormaps.R2standard,'CDataMapping','scaled')
ax=gca;
ax.CLim=[-100 100];
FigLabels;
title('Gaussian R_2^* error')
h=colorbar
h.Label.String = "R_2^* error (s^-^1)";
h.Label.FontSize = 12;

subplot(3,3,5)
image(1000*errormaps.R2Rician,'CDataMapping','scaled')
ax=gca;
ax.CLim=[-100 100];
FigLabels;
title('Rician R_2^* error')
h=colorbar
h.Label.String = "R_2^* error (s^-^1)";
h.Label.FontSize = 12;

subplot(3,3,6)
image(1000*dlErrormaps.R2,'CDataMapping','scaled')
ax=gca;
ax.CLim=[-100 100];
FigLabels;
title(['DL R_2^* error'])
h=colorbar
h.Label.String = "R_2^* error (s^-^1)";
h.Label.FontSize = 12;


subplot(3,3,7)
image(1000*sdMaps.R2standard,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 200];
FigLabels;
title('Gaussian R_2^* SD')
h=colorbar
h.Label.String = "R2* SD (s^-^1)";
h.Label.FontSize = 12;

subplot(3,3,8)
image(1000*sdMaps.R2Rician,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 200];
FigLabels;
title('Rician R_2^* SD')
h=colorbar
h.Label.String = "R2* SD (s^-^1)";
h.Label.FontSize = 12;

subplot(3,3,9)
image(1000*dlSdMaps.r2,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 200];
FigLabels;
title('DL R_2^* SD')
h=colorbar
h.Label.String = "R2* SD (s^-^1)";
h.Label.FontSize = 12;

%1.5 Print mean absolute error
sum(abs(sdMaps.FFstandard),'all')
sum(abs(dlSdMaps.ff),'all')
sum(abs(sdMaps.R2standard),'all')
sum(abs(dlSdMaps.r2),'all')


