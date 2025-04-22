

function [dlMaps, dlErrormaps, sdMaps] = testOnSimulatedData(nets,settings)
%function testOnSimulatedData(nets,settings)

%% 6.  Create simulated test dataset with constant spacing of FF values and R2* values (to aid visualisation)

%Select values for test dataset
r2max=0.5;

%Extract echotimes
echotimes = settings.echotimes;
tesla = settings.fieldStrength;

%Set parameters for simulation
S0=10;
FFvals = (0:0.01:1)';
R2vals = (0:0.05*r2max:r2max);

%Set noise SD based on SNR
noiseSD = S0/settings.SNR;
settings.noiseSD = settings;

%Select number of repetitions

%Call helper function (sVecFixedSpacing) to generate vectors of values with
%fixed spacing
[paramVec, sVecNoiseFree] = sVecFixedSpacing(S0,FFvals,R2vals,settings);

%Create yTest parameter vector
yTest = paramVec;


%If noise free simulation, set noiseSD = 0
if settings.noisyTesting == 0

    settings.sigmaEst = 0
elseif settings.noisyTesting == 1
    % Otherwise, estimate noiseSD using test voxels
    settings.sigmaEst = estimateSigmaForSimulation(S0, settings, 100);
else ;
end


%Loop over noise instantiations
reps = 100;

parfor r = 1:reps

    %Create noise
    realnoise=(noiseSD)*randn(size(sVecNoiseFree,1),numel(echotimes));
    imagnoise=1i*(noiseSD)*randn(size(sVecNoiseFree,1),numel(echotimes));

    %Add noise to signal
    if settings.noisyTesting == 0
        sVec = sVecNoiseFree;
    elseif settings.noisyTesting == 1
        sVec = sVecNoiseFree + realnoise + imagnoise;
    end

    %Magnitude
    sVec = abs(sVec);

    xTestNoNorm(:,:,r) = sVec;

    % Normalisation
    sVecNorm = normaliseSignals(sVec,settings);

    xTest(:,:,r) = sVecNorm;

    % Complex
    % xTest(:,:,r) = horzcat(real(sVec),imag(sVec));

end

%% Choose non-normalised data (for experimentation only)
% xTest = xTestNoNorm;

%% 8.  Generate prediction from either net1 (water-dominant solution training)
%or net2 (fat-dominant solution training) depending on SSE

%8.1 Get predictions
likVec2(1,1:2121,1:100) = zeros;
likVec2(1,1:2121,1:100) = zeros;

choiceVecRic(1:2121,1,1:100) = zeros;
choiceVecSSE(1:2121,1,1:100) = zeros;
choiceVec = choiceVecSSE;


%% Loop over 'voxels' / noise instantiations (conventional for loop)
tic
for r = 1:reps

    %Get net1 predictions
    predictionVec1(:,:,r)= nets.net1.predict(xTest(:,:,r));

    %Get net2 predictions
    predictionVec2(:,:,r)= nets.net2.predict(xTest(:,:,r));

    %Combine predictions
    predictionVec3(:,:,r) = combinePredictions(predictionVec1(:,:,r), predictionVec2(:,:,r), settings, xTestNoNorm(:,:,r), xTest(:,:,r));

end
toc

%Get predictions for all noise instantiations together (doesn't work yet)

% %First reorganise data
% xTestPermuted = permute(xTest,[1 3 2]);
% xTestReshaped = reshape(xTestPermuted,[size(xTest,1)*size(xTest,3),size(xTest,2),1]);
%
% %Get predictions
% prediction1 =nets.net1.predict(xTestReshaped);
% prediction2 =nets.net2.predict(xTestReshaped);
%
% %With image-based normalisation for likelihood calc
% prediction3 = combinePredictions(prediction1, prediction2, settings, xTestReshaped, xTestReshaped);
%
% prediction1 = reshape(prediction1,size(xTest,1),reps,2);
% prediction1 = permute(prediction1,[1 3 2]);
% prediction2 = reshape(prediction2,size(xTest,1),reps,2);
% prediction2 = permute(prediction2,[1 3 2]);
% prediction3 = reshape(prediction3,size(xTest,1),reps,3);
% prediction3 = permute(prediction3,[1 3 2]);


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


%Show summaries for single instantion, mean and SD over instantions

%dispSl specifies the instantiation
dispSl = 1;

createFigDL(predictionVec1(:,:,dispSl), std(predictionVec1,0,3), yTest,FFvals,R2vals,S0, 'Low FF net, mean values');
createFigDL(predictionVec2(:,:,dispSl), std(predictionVec2,0,3), yTest,FFvals,R2vals,S0, 'High FF net, mean values');
createFigDL(predictionVec3(:,:,dispSl), std(predictionVec3,0,3), yTest,FFvals,R2vals,S0,'Likelihood combined nets, mean values');
% createFigDL(predictionVec4(:,:,dispSl), yTest,FFvals,R2vals,'Third network, mean values');

createFigDL(mean(predictionVec1,3), std(predictionVec1,0,3), yTest,FFvals,R2vals,S0, 'Low FF net, mean values');
createFigDL(mean(predictionVec2,3), std(predictionVec2,0,3), yTest,FFvals,R2vals,S0, 'High FF net, mean values');
[dlMaps,dlErrormaps,sdMaps] = createFigDL(mean(predictionVec3,3), std(predictionVec3,0,3), yTest,FFvals,R2vals,S0,'Likelihood combined nets, mean values');
% createFigDL(mean(predictionVec4,3), yTest,FFvals,R2vals,'Third network, mean values');

% createSDFigDL(std(predictionVec1,0,3), yTest,FFvals,R2vals, 'Low FF net, std');
% createSDFigDL(std(predictionVec2,0,3), yTest,FFvals,R2vals, 'High FF net, std');
% [sdMaps] = createSDFigDL(std(predictionVec3,0,3), yTest,FFvals,R2vals,'Likelihood combined nets, std');
% createFigDL(std(predictionVec4,0,3), yTest,FFvals,R2vals,'Third network, std');

%% Show a histogram of values for chosen ground truth FF value
% gtFF = 0.94;
% gtR2 = 0.0;
%
% index = find(yTest(:,1)==gtFF & yTest(:,2)==gtR2);
%
% %Get values for chosen ground truth
% ffValues = predictionVec3(index,1,:);
% r2Values = predictionVec3(index,2,:);
%
% %Reshape
% ffValues = reshape(ffValues,[reps 1]);
% r2Values = reshape(r2Values,[reps 1]);
%
% %Plot
% figure
% s1=scatter(ffValues,r2Values);
% xticks([0 0.2 0.4 0.6 0.8 1.0]);
% xticklabels({'0','0.2','0.4', '0.6', '0.8','1.0'});
% yticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]);
% yticklabels({'0','100', '200', '300', '400', '500', '600', '700', '800', '900','1000'});
% % set(s1,'YGrid','on','GridAlpha',0.5)
% title('FF / R2* scatterplot (Gaussian)')
% ylabel('R2* (s^-1)','FontSize',12)
% xlabel('PDFF estimate','FontSize',12)
% xlim([0 1])
% ylim([0 1]); %get ylim to enable setting for next plot
% hold on
% plot(gtFF,gtR2,'rd','MarkerSize',8,'MarkerFaceColor','red') %..add ground truth as point
% hold off
% legend('Estimates','Ground truth')


end






