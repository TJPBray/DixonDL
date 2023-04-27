function [tableVar] = tabulateCoeffs(ff,regressionModels,site, protocol, ReferenceValues,saveFolderInfo)
% function tableVar = tabulateCoeffs[regressionModels,saveFolderInfo]

%% Create subplots with reference values to plot on later
figure
% Gaussian subplots
g(1)=subplot(2,2,1)
plot(ReferenceValues.FF,ReferenceValues.FF,'k-','LineWidth',1)
title('1.5T protocol 1')

g(2)=subplot(2,2,2)
plot(ReferenceValues.FF,ReferenceValues.FF,'k-','LineWidth',1)
title('1.5T protocol 2')

g(3)=subplot(2,2,3)
plot(ReferenceValues.FF,ReferenceValues.FF,'k-','LineWidth',1)
title('3T protocol 1')

g(4)=subplot(2,2,4)
plot(ReferenceValues.FF,ReferenceValues.FF,'k-','LineWidth',1)
title('3T protocol 2')

% % Rician subplots
% r(1)=subplot(4,2,2)
% plot(ReferenceValues.FF,ReferenceValues.FF,'k-','LineWidth',1)
% 
% r(2)=subplot(4,2,4)
% plot(ReferenceValues.FF,ReferenceValues.FF,'k-','LineWidth',1)
% 
% r(3)=subplot(4,2,6)
% plot(ReferenceValues.FF,ReferenceValues.FF,'k-','LineWidth',1)
% 
% r(4)=subplot(4,2,8)
% plot(ReferenceValues.FF,ReferenceValues.FF,'k-','LineWidth',1)


%% Create colours and symbols for display in loops
colours = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0; 0 0.4 0.7]; %red green blue cyan magenta yellow otherblue (semitransparent)
symbols = {'x','+','>','d','o','<','x'};
sizes = [8 6 5 5 4 3 2];

%% Loop over protocols/sites
for n = 1:numel(regressionModels)

%% 1. For DL

%Rsquared
Rsquared(n,1)=regressionModels{n}.Rsquared.Ordinary;

%Coefficients
slope(n,1)=table2array(regressionModels{n}.Coefficients(2,1));
intercept(n,1)=table2array(regressionModels{n}.Coefficients(1,1));
    
%95% CIs
cis = coefCI(regressionModels{n});
slopeLowerCI(n,1)=cis(2,1);
slopeUpperCI(n,1)=cis(2,2);
interceptLowerCI(n,1)=cis(1,1);
interceptUpperCI(n,1)=cis(1,2);

axes(g(protocol(n))) %Select axis for relevant site
hold on 
plot(ReferenceValues.FF,ff{n}.DL.median,'--','LineWidth',1,'Color',colours(site(n),:),'Marker',symbols{site(n)},'MarkerSize',sizes(site(n)))
xlabel('Reference FF','FontSize',12)
ylabel('DL PDFF','FontSize',12)
xlim([0 1])
ylim([0 1])
hold off


%% 2. Protocol name
protocolInfo = saveFolderInfo(n+2,1);
protName{n}=protocolInfo.name;

end


%% 3. Create table with regression values for export

import mlreportgen.dom.*
setDefaultNumberFormat("%0.3f");

varNames= {'Protocol name','Protocol number','Site','R^2', 'Slope', 'Slope 95% lower', 'Slope 95% upper', 'Intercept', 'Intercept 95% lower', 'Intercept 95% upper'}
tableVar = table(protName',protocol',site',Rsquared,slope,slopeLowerCI,slopeUpperCI,intercept,interceptLowerCI,interceptUpperCI,'VariableNames',varNames)

%Sort by protocol number before export
tableVar = sortrows(tableVar,2,'ascend');


end

