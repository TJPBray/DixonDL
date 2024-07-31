function combinedPredictions = combinePredictions(pred1, pred2, settings, signals)
% function combinedPredictions = combinePredictions(pred1, pred2)

%Takes predictions from the two networks (water network and fat network)
%and chooses the best one based on likelihood

% Input:
% pred1 is m x n matrix of predictions from first network
% pred2 is m x n matrix of predictions from second network
% settings is a structure containing the field strenght, echotimes and noise sigma
% signals is a m x t matrix containing the signals (note that if the
% predictions were generated using normalised signals, this signal should
% also be normalised)

% Output:
% combinedPredictions is the chosen predictions from pred1 and pred2

%% Get settings from settings structure
noiseSD = settings.noiseSD;
echotimes = settings.echotimes;
tesla = settings.fieldStrength;

%Store 'raw' predictions with no clipping
rawPred1 = pred1;
rawPred2 = pred2;

% % Eliminate extreme values of parameters before calculating sse / likelihood
pred1(pred1>1) = 1; 
pred1(pred1<0) = 0;

pred2(pred2>1) = 1; 
pred2(pred2<0) = 0;

%% Get likelihood for different predictions

%Net1
[lik1,sse1] = sseVecCalc (echotimes, tesla, pred1, signals, noiseSD);

%Net2
[lik2,sse2] = sseVecCalc (echotimes, tesla, pred2, signals, noiseSD);

%Create binary vector to choose between values
choiceVecRic =(lik1>lik2)';
choiceVecSSE=(sse1<sse2)';

choiceVec = choiceVecSSE;


%% Create predictionVec with best likelihood estimates:

%Get combined predictions
combinedPredictions = choiceVec.*pred1 + (1-choiceVec).*pred2;

%% Use information about extreme values
% If extreme R2*, assume incorrect and make a choice on this basis 
% 

for k = 1:numel(size(pred1,1))

if rawPred2(k,2) < 0
combinedPredictions(k,:) = pred1(k,:);

% % If extreme high PDFF from fat network, assume incorrect and make a choice on this basis 
% elseif rawPred2(k,1) > 1 
% 
% combinedPredictions(k,:) = pred1(k,:);

%Otherwise, choose between predictions based on SSE / likelihood
else ;
end

% If extreme low PDFF from water network, assume incorrect and make a choice on this basis 

%Only do this if the relevant option is turned on 
% if settings.tuning.useExtremeLowPdff == 1;
% 
%     if rawPred1(k,1) < -0.05
% 
%         combinedPredictions(k,:) = pred2(k,:);
%     else ;
%     end
% 
% else ;
% end

end

end





