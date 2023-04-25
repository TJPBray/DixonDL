

function testOnSimulatedData(nets,settings)
%function testOnSimulatedData(nets,settings)

%% 6.  Create simulated test dataset with constant spacing of FF values and R2* values (to aid visualisation)

%Select values for test dataset
r2max=0.5;

%Extract echotimes
echotimes = settings.echotimes;
noiseSD = settings.noiseSD;
tesla = settings.fieldStrength;

S0=1;
FFvals = (0:0.01:1)';
R2vals = (0:0.05*r2max:r2max);

%Select number of repetitions

%Call helper function (sVecFixedSpacing) to generate vectors of values with
%fixed spacing
[paramVec, sVecNoiseFree] = sVecFixedSpacing(S0,FFvals,R2vals,echotimes);

%Create yTest parameter vector
yTest = paramVec;

%Loop over noise instantiations 

reps = 100;

for r = 1:reps

%Create noise
realnoise=(noiseSD)*randn(size(sVecNoiseFree,1),numel(echotimes));
imagnoise=1i*(noiseSD)*randn(size(sVecNoiseFree,1),numel(echotimes));

%Add noise to signal 
sVec = sVecNoiseFree + realnoise + imagnoise;

%Magnitude
sVec = abs(sVec);

xTestNoNorm(:,:,r) = sVec;

%% Normalisation
sVecNorm = normaliseSignals(sVec,settings);

xTest(:,:,r) = sVecNorm;

% Complex
% xTest(:,:,r) = horzcat(real(sVec),imag(sVec));

end

%% 8.  Generate prediction from either net1 (water-dominant solution training)
%or net2 (fat-dominant solution training) depending on SSE

%8.1 Get predictions
likVec2(1,1:2121,1:100) = zeros;
likVec2(1,1:2121,1:100) = zeros;

choiceVecRic(1:2121,1,1:100) = zeros;
choiceVecSSE(1:2121,1,1:100) = zeros;
choiceVec = choiceVecSSE;

%% Loop over 'voxels' / noise instantiations (conventional for loop)
for r = 1:reps

%Net1
predictionVec1(:,:,r)=nets.net1.predict(xTest(:,:,r));
[likVec1(:,:,r),sseVec1(:,:,r)] = sseVecCalc (echotimes, tesla, predictionVec1(:,:,r), xTestNoNorm(:,:,r), noiseSD);

%Net2
predictionVec2(:,:,r)=nets.net2.predict(xTest(:,:,r));
[likVec2(:,:,r),sseVec2(:,:,r)] = sseVecCalc (echotimes, tesla, predictionVec2(:,:,r), xTestNoNorm(:,:,r), noiseSD);

%8.2 Create binary vector to choose between values
choiceVecRic(:,:,r)=(likVec1(:,:,r)>likVec2(:,:,r))';
choiceVecSSE(:,:,r)=(sseVec1(:,:,r)<sseVec2(:,:,r))';

choiceVec = choiceVecSSE;

%8.3 Create predictionVec with best likelihood estimates:
predictionVec3(:,:,r) = choiceVec(:,:,r).*predictionVec1(:,:,r) + (1-choiceVec(:,:,r)).*predictionVec2(:,:,r);

%8.4 Create predictionVec by using the third network 
% 
% xTest2(:,:,r) = horzcat(sseVec1(:,:,r)', sseVec2(:,:,r)', predictionVec1(:,:,r), predictionVec2(:,:,r));
% 
% predictionVec4(:,:,r) = net3.predict(xTest2(:,:,r));

end

% %% Loop over 'voxels' / noise instantiations (parfor)
% parfor r = 1:reps
% 
% %Net1
% predictionVec1=nets.net1.predict(xTest(:,:,r));
% [likVec1,sseVec1] = sseVecCalc (echotimes, tesla, predictionVec1, xTest(:,:,r), noiseSD);
% 
% %Net2
% predictionVec2=nets.net2.predict(xTest(:,:,r));
% [likVec2,sseVec2] = sseVecCalc (echotimes, tesla, predictionVec2, xTest(:,:,r), noiseSD);
% 
% %8.2 Create binary vector to choose between values
% choiceVecRic=(likVec1>likVec2)';
% choiceVecSSE=(sseVec1<sseVec2)';
% 
% choiceVec = choiceVecSSE;
% 
% %8.3 Create predictionVec with best likelihood estimates:
% predictionVec3(:,:,r) = choiceVec.*predictionVec1 + (1-choiceVec).*predictionVec2;
% 
% %8.4 Create predictionVec by using the third network 
% % 
% % xTest2(:,:,r) = horzcat(sseVec1(:,:,r)', sseVec2(:,:,r)', predictionVec1(:,:,r), predictionVec2(:,:,r));
% % 
% % predictionVec4(:,:,r) = net3.predict(xTest2(:,:,r));
% 
% end


%% 9. Visualise predicted values vs ground truth (use all data to ease visualisation)
dispSl = 1;

createFigDL(predictionVec1(:,:,dispSl), yTest,FFvals,R2vals, 'Low FF net, mean values');
createFigDL(predictionVec2(:,:,dispSl), yTest,FFvals,R2vals, 'High FF net, mean values');
createFigDL(predictionVec3(:,:,dispSl), yTest,FFvals,R2vals,'Likelihood combined nets, mean values');
% createFigDL(predictionVec4(:,:,dispSl), yTest,FFvals,R2vals,'Third network, mean values');

createFigDL(mean(predictionVec1,3), yTest,FFvals,R2vals, 'Low FF net, mean values');
createFigDL(mean(predictionVec2,3), yTest,FFvals,R2vals, 'High FF net, mean values');
createFigDL(mean(predictionVec3,3), yTest,FFvals,R2vals,'Likelihood combined nets, mean values');
% createFigDL(mean(predictionVec4,3), yTest,FFvals,R2vals,'Third network, mean values');

createFigDL(std(predictionVec1,0,3), yTest,FFvals,R2vals, 'Low FF net, std');
createFigDL(std(predictionVec2,0,3), yTest,FFvals,R2vals, 'High FF net, std');
createFigDL(std(predictionVec3,0,3), yTest,FFvals,R2vals,'Likelihood combined nets, std');
% createFigDL(std(predictionVec4,0,3), yTest,FFvals,R2vals,'Third network, std');

% Visualise comparison against conventional fitting
load('/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO/Experiments/experimentResults/Rev2_Simulation_Results_SNR60_R1000.mat')

figure('Name', 'MAGORINO parameter error')
subplot(2,2,1)
image(FFmaps.standard(:,1:6),'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 1];
FigLabels;
title('Gaussian magnitude FF maps')
colorbar

subplot(2,2,2)
image(FFmaps.Rician(:,1:6),'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 1];
FigLabels;
title('Rician magnitude FF maps')
colorbar

subplot(2,2,3)
image(errormaps.FFstandard(:,1:6),'CDataMapping','scaled')
ax=gca;
ax.CLim=[-1 1];
FigLabels;
title('Gaussian magnitude FF error')
colorbar

subplot(2,2,4)
image(errormaps.FFRician(:,1:6),'CDataMapping','scaled')
ax=gca;
ax.CLim=[-1 1];
FigLabels;
title('Rician magnitude FF error')
colorbar


