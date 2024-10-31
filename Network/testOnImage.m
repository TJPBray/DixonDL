%% Uses training neural networks to fit an image 

function predictions = testOnImage(image,nets,settings)
%function predictions = testOnImage(image,nets,settings)

%% Get S0 estimate from whole image 

%First pick slice 
sliceIntensity = mean(image,[1 2]);
[maxSlice,slice] = max(sliceIntensity);

%Filter image 
% filteredImage = imgaussfilt(image(:,:,slice))

%Take 99th percentile 
prc = prctile(image(:,:,2),99,'all');
maxImage = max(image(:,:,:),[],'all');

%% Loop over voxels 

%Specify indent
indent = 0;

%Prefill arrays
rawPrediction1 = zeros(size(image,1),size(image,2),2);
rawPrediction2 = rawPrediction1;

%Vectorised implementation of prediction
signalsMat = reshape(image,[size(image,1)*size(image,2),size(image,3)]);

%Normalise signals
normalisedSignalsMat = normaliseSignals(signalsMat,settings);

%Get predictions
prediction1 =nets.net1.predict(normalisedSignalsMat);
prediction2 =nets.net2.predict(normalisedSignalsMat);

%With image-based normalisation for likelihood calc
prediction3 = combinePredictions(prediction1, prediction2, settings, signalsMat, normalisedSignalsMat);

%Reshape
prediction1 = reshape(prediction1, [size(image,1) size(image,2) 2]);
prediction2 = reshape(prediction2, [size(image,1) size(image,2) 2]);
prediction3 = reshape(prediction3, [size(image,1) size(image,2) 3]);


% %Start loop
% for y = indent+1:(size(image,1)-indent)

% parfor x = indent+1:(size(image,2)-indent)
% 
% signals = image(y,x,:);
% signals = reshape(signals,[1 6]);
% 
% normalisedSignals = normaliseSignals(signals,settings);
% 
% %Get predictions
% prediction1(y,x,:)=nets.net1.predict(normalisedSignals);
% prediction2(y,x,:)=nets.net2.predict(normalisedSignals);
% 
% %Combine predictions
% %With image-based normalisation for likelihood calc
% prediction3(y,x,:) = combinePredictions(prediction1(y,x,:), prediction2(y,x,:), settings, signals/maxS);

%With voxel-based normalisation for likelihood calc
% prediction3(y,x,:) = combinePredictions(prediction1(y,x,:), prediction2(y,x,:), settings, normalisedSignals);
% % 
% 
% %Store 'raw' predictions with no clipping
% rawPred1 = pred1;
% rawPred2 = pred2;
% rawPrediction1(y,x,:) = rawPred1(:);
% rawPrediction2(y,x,:) = rawPred2(:);
% 
% % Eliminate extreme values of parameters before calculating sse / likelihood
% pred1(pred1>1) = 1; 
% pred1(pred1<0) = 0;
% prediction1(y,x,:) = pred1;
% 
% pred2(pred2>1) = 1; 
% pred2(pred2<0) = 0;
% prediction2(y,x,:)= pred2;
% 
% %Calculate sse/likelihood
% % [lik1,sse1] = sseVecCalc (settings.echotimes, settings.fieldStrength, pred1, signals, settings.noiseSD);
% % [lik2,sse2] = sseVecCalc (settings.echotimes, settings.fieldStrength, pred2, signals, settings.noiseSD);
% 
% [lik1,sse1] = sseVecCalc (settings.echotimes, settings.fieldStrength, pred1, normalisedSignals, settings.noiseSD);
% [lik2,sse2] = sseVecCalc (settings.echotimes, settings.fieldStrength, pred2, normalisedSignals, settings.noiseSD);
% 
% %Create binary choice value to choose between values
% choiceRic=(lik1>lik2);
% choiceSSE=(sse1<sse2);
% 
% choiceImage(y,x) = choiceSSE;
% 
% %Create predictionVec with best likelihood estimates:
% prediction3(y,x,:) = choiceSSE.*prediction1(y,x,:) + (1-choiceSSE).*prediction2(y,x,:);
% 
% % If extreme R2*, assume incorrect and make a choice on this basis 
% if rawPred2(2) < 0
%     
% prediction3(y,x,:) = pred1;
% 
% % If extreme PDFF, assume incorrect and make a choice on this basis 
% elseif rawPred2(1) > 1 
% 
% prediction3(y,x,:) = pred1;
% 
% %Otherwise, choose between predictions based on SSE / likelihood
% else ;
% end

% end
% 
% y
% end

%% Output predictions
%With parameter clipping prior to choice
predictions.prediction1 = prediction1;
predictions.prediction2 = prediction2;
predictions.prediction3 = prediction3;

%Unclipped
predictions.rawPrediction1 = rawPrediction1;
predictions.rawPrediction2 = rawPrediction2;

end

