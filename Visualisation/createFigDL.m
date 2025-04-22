function  [maps, errormaps, sdMaps] = createFigDL(predictionVec,sdVec,yTest,FFvals,R2vals,S0val,figTitle)
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
% S0val is a scalar S0 value used for testing
% title is a string specifying the title of the figure

%Split prediction vector into FF and R2*
predictionVecFF=predictionVec(:,1);
predictionVecR2=predictionVec(:,2);
try
    predictionVecS0=predictionVec(:,3); %Do this if present (only for combined predictions in current implementation)
catch
end



%FF
%Get errorgrids
ffErrorVec=predictionVecFF-yTest(:,1);
r2ErrorVec=predictionVecR2-yTest(:,2);
try
    s0ErrorVec=predictionVecS0-S0val;
catch
end

% Reshape prediction data for plotting
ffPredictions = reshape(predictionVecFF,[numel(FFvals) numel(R2vals)]);
r2Predictions = reshape(predictionVecR2,[numel(FFvals) numel(R2vals)]);
try
    s0Predictions = reshape(predictionVecS0,[numel(FFvals) numel(R2vals)]);
catch
end

ffError = reshape(ffErrorVec,[numel(FFvals) numel(R2vals)]);
r2Error = reshape(r2ErrorVec,[numel(FFvals) numel(R2vals)]);
try
    s0Error = reshape(s0ErrorVec,[numel(FFvals) numel(R2vals)]);
catch
end

%SD
% Split prediction vector into FF and R2*
sdVecFF=sdVec(:,1);
sdVecR2=sdVec(:,2);
try
    sdVecS0=sdVec(:,3);
catch
end

%Get errorgrids
% ffSdVec=sdVecFF-yTest(:,1);
% r2SdVec=sdVecR2-yTest(:,2);
% s0SdVec=sdVecS0-S0val;

% Reshape prediction data for plotting
ffSd = reshape(sdVecFF,[numel(FFvals) numel(R2vals)]);
r2Sd = reshape(sdVecR2,[numel(FFvals) numel(R2vals)]);
try
    s0Sd = reshape(sdVecS0,[numel(FFvals) numel(R2vals)]);
catch
end

%Plot FF and R2*

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


%Plot S0 if available
if exist('s0Predictions','var') == 1

    figure('Name', figTitle)

    subplot(3,1,1)
    image(s0Predictions,'CDataMapping','scaled')
    ax=gca;
    ax.CLim=[0.8*S0val 1.2*S0val];
    FigLabels;
    colorbar
    xticks([1:2:21]);
    xticklabels({'0','.50', '100', '150', '200', '250', '300', '350', '400', '450','500'});
    xlabel('R2* (s^-^1)','FontSize',12)
    yticks([0 10 20 30 40 50 60 70 80 90 100]);
    yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});
    ylabel('Fat fraction','FontSize',12)
    title('S0 values')
    h=colorbar;
    h.Label.String = "S_0";
    h.Label.FontSize = 12;
    view(0,90)

    subplot(3,1,2)
    image(s0Error,'CDataMapping','scaled')
    ax=gca;
    ax.CLim=[-0.1*S0val 0.1*S0val];
    FigLabels;
    colorbar
    xticks([1:2:21]);
    xticklabels({'0','.50', '100', '150', '200', '250', '300', '350', '400', '450','500'});
    xlabel('R2* (s^-^1)','FontSize',12)
    yticks([0 10 20 30 40 50 60 70 80 90 100]);
    yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});
    ylabel('Fat fraction','FontSize',12)
    title('S0 error')
    h=colorbar;
    h.Label.String = "S_0 error";
    h.Label.FontSize = 12;
    view(0,90)

    subplot(3,1,3)
    image(s0Sd,'CDataMapping','scaled')
    ax=gca;
    ax.CLim=[0 0.1*S0val];
    FigLabels;
    colorbar
    xticks([1:2:21]);
    xticklabels({'0','.50', '100', '150', '200', '250', '300', '350', '400', '450','500'});
    xlabel('R2* (s^-^1)','FontSize',12)
    yticks([0 10 20 30 40 50 60 70 80 90 100]);
    yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});
    ylabel('Fat fraction','FontSize',12)
    title('S0 SD')
    h=colorbar;
    h.Label.String = "S_0 SD";
    h.Label.FontSize = 12;
    view(0,90)

else ;
end

%% Export
maps.ff = ffPredictions;
errormaps.ff = ffError;
sdMaps.ff = ffSd;


maps.r2 = r2Predictions;
errormaps.R2 = r2Error;
sdMaps.r2 = r2Sd;


try
    maps.s0 = s0Predictions;
    errormaps.S0 = s0Error;
    sdMaps.s0 = s0Sd;
catch
end


end

