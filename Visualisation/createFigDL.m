function  [maps, errormaps, sdMaps] = createFigDL(predictionVec,sdVec,yTest,FFvals,R2vals,figTitle)
% function  createFigDL(predictionVec,ytest)

%Helper function for visualising performance of DL methods over a range of
%FF values

%Written by Dr Tim Bray
%t.bray@ucl.ac.uk
%November 2022

%Inputs
% predictionVec is an nx2 prediction of FF and R2* values
% sdVec is an nx2 matrix of SD FF and R2* values
% yTest is an nx2 matrix of ground truth FF and R2* values
% FFvals is an nx1 vector of FF values used for testing
% R2vals is an nx1 vector of R2* values used for testing
% title is a string specifying the title of the figure

%Split prediction vector into FF and R2*
predictionVecFF=predictionVec(:,1);
predictionVecR2=predictionVec(:,2);

%FF 
%Get errorgrids
ffErrorVec=predictionVecFF-yTest(:,1);
r2ErrorVec=predictionVecR2-yTest(:,2);

% Reshape prediction data for plotting
ffPredictions = reshape(predictionVecFF,[numel(FFvals) numel(R2vals)]);
r2Predictions = reshape(predictionVecR2,[numel(FFvals) numel(R2vals)]);
ffError = reshape(ffErrorVec,[numel(FFvals) numel(R2vals)]);
r2Error = reshape(r2ErrorVec,[numel(FFvals) numel(R2vals)]);

%SD
% Split prediction vector into FF and R2*
sdVecFF=sdVec(:,1);
sdVecR2=sdVec(:,2);

%Get errorgrids
ffSdVec=sdVecFF-yTest(:,1);
r2SdVec=sdVecR2-yTest(:,2);

% Reshape prediction data for plotting
ffSd = reshape(sdVecFF,[numel(FFvals) numel(R2vals)]);
r2Sd = reshape(sdVecR2,[numel(FFvals) numel(R2vals)]);


%Plot 

figure('Name', figTitle)

subplot(3,2,1)
image(ffPredictions,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 1];
FigLabels;
colorbar
xticks([1:2:21]);
xticklabels({'0','.50', '100', '150', '200', '250', '300', '350', '400', '450','500'});
xlabel('R2* (s^-^1)','FontSize',12)
yticks([0 10 20 30 40 50 60 70 80 90 100]);
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});
ylabel('Fat fraction','FontSize',12)
title('FF values')
h=colorbar;
h.Label.String = "FF";
h.Label.FontSize = 12;
view(0,90)


subplot(3,2,2)
image(1000*r2Predictions,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 500];
FigLabels;
xticks([1:2:21]);
xticklabels({'0','50', '100', '150', '200', '250', '300', '350', '400', '450','500'});
xlabel('R2* (s^-^1)','FontSize',12)
yticks([0 10 20 30 40 50 60 70 80 90 100]);
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});
ylabel('Fat fraction','FontSize',12)
title('R2* values')
h=colorbar;
h.Label.String = "R2* (s^-^1)";
h.Label.FontSize = 12;
view(0,90)

subplot(3,2,3)
image(ffError,'CDataMapping','scaled')
ax=gca;
ax.CLim=[-.2 .2];
FigLabels;
colorbar
xticks([1:2:21]);
xticklabels({'0','.50', '100', '150', '200', '250', '300', '350', '400', '450','500'});
xlabel('R2* (s^-^1)','FontSize',12)
yticks([0 10 20 30 40 50 60 70 80 90 100]);
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});
ylabel('Fat fraction','FontSize',12)
title('FF error')
h=colorbar;
h.Label.String = "FF error";
h.Label.FontSize = 12;
view(0,90)

subplot(3,2,5)
image(ffSd,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 .3];
FigLabels;
colorbar
xticks([1:2:21]);
xticklabels({'0','.50', '100', '150', '200', '250', '300', '350', '400', '450','500'});
xlabel('R2* (s^-^1)','FontSize',12)
yticks([0 10 20 30 40 50 60 70 80 90 100]);
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});
ylabel('Fat fraction','FontSize',12)
title('FF SD')
h=colorbar;
h.Label.String = "FF";
h.Label.FontSize = 12;
view(0,90)

subplot(3,2,4)
image(1000*r2Error,'CDataMapping','scaled')
ax=gca;
ax.CLim=[-100 100];
FigLabels;
colorbar
xticks([1:2:21]);
xticklabels({'0','.50', '100', '150', '200', '250', '300', '350', '400', '450','500'});
xlabel('R2* (s^-^1)','FontSize',12)
yticks([0 10 20 30 40 50 60 70 80 90 100]);
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});
ylabel('Fat fraction','FontSize',12)
title('R2* error')
h=colorbar;
h.Label.String = "R2* error (s^-^1)";
h.Label.FontSize = 12;
view(0,90)

subplot(3,2,6)
image(1000*r2Sd,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 100];
FigLabels;
xticks([1:2:21]);
xticklabels({'0','50', '100', '150', '200', '250', '300', '350', '400', '450','500'});
xlabel('R2* (s^-^1)','FontSize',12)
yticks([0 10 20 30 40 50 60 70 80 90 100]);
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});
ylabel('Fat fraction','FontSize',12)
title('R2* SD')
h=colorbar;
h.Label.String = "R2* (s^-^1)";
h.Label.FontSize = 12;
view(0,90)

%% Export 
maps.ff = ffPredictions;
maps.r2 = r2Predictions;

errormaps.ff = ffError;
errormaps.R2 = r2Error;

sdMaps.ff = ffSd;
sdMaps.r2 = r2Sd;

end

