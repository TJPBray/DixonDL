function combinedPredictions = combinePredictions(pred1, pred2, settings, signals, normsignals)
% function combinedPredictions = combinePredictions(pred1, pred2)

%Takes predictions from the two networks (water network and fat network)
%and chooses the best one based on likelihood

% Input:
% pred1 is m x n matrix of predictions from first network
% pred2 is m x n matrix of predictions from second network
% settings is a structure containing the field strenght, echotimes and noise sigma
% signals is a m x t matrix containing the signals 
% normsignals is a m x t matrix containing the normalised signals (both
% pieces of information are available)

% Output:
% combinedPredictions is the chosen predictions from pred1 and pred2

%% Get settings from settings structure
echotimes = settings.echotimes;
tesla = settings.fieldStrength;
sigmaEst = settings.sigmaEst;

%Store 'raw' predictions with no clipping
rawPred1 = pred1;
rawPred2 = pred2;

% % Eliminate extreme values of parameters before calculating sse / likelihood
% pred1(pred1>1) = 1; 
% pred1(pred1<0) = 0;
% 
% pred2(pred2>1) = 1; 
% pred2(pred2<0) = 0;

%% Calculate s0 estimates ahead of likelihood calculation
%Use the fact that S = A*S0 where S0 is a scalar factor and A is the
%remainder of the signal model 

% Needs vectorising 

%Net1
Anet1 = abs(MultiPeakFatSingleR2(echotimes,tesla,pred1(:,1),1-pred1(:,1),pred1(:,2),0));

%Net2
Anet2 = abs(MultiPeakFatSingleR2(echotimes,tesla,pred2(:,1),1-pred2(:,1),pred2(:,2),0));

% Loop implementation
parfor k = 1:size(Anet1,1)

    a1 = Anet1(k,:)';
    a2 = Anet2(k,:)';

    svox = signals(k,:)';

    s0est1(k,1) = pinv(a1)*svox;
    s0est2(k,1) = pinv(a2)*svox;

end

% Add S0 estimates to predictions 
pred1(:,3) = s0est1;
pred2(:,3) = s0est2;


%% Get likelihood for different predictions

%Net1
[sse1,lik1] = sseVecCalc (echotimes, tesla, pred1, signals, sigmaEst);

%Net2
[sse2,lik2] = sseVecCalc (echotimes, tesla, pred2, signals, sigmaEst);

%Create binary vector to choose between values
choiceVecRic =(lik1>lik2)';
choiceVecSSE=(sse1<sse2)';

if sigmaEst == 0
choiceVec = choiceVecSSE;  % Rician distribution becomes Gaussian for sigma = 0, so use Gaussian to avoid NaN 
else 
choiceVec = choiceVecRic;
end

%% Create predictionVec with best likelihood estimates:

%Get combined predictions
combinedPredictions = choiceVec.*pred1 + (1-choiceVec).*pred2;

%% Use information about implausible values


if settings.tuning.useImplausibleValues == 1;

% If extreme R2*, assume incorrect and make a choice on this basis 
parfor k = 1:size(pred1,1)

% If negative R2* from fat network, assume incorrect and choose water
% network
if rawPred2(k,2) < 0
combinedPredictions(k,:) = pred1(k,:);

% If negative R2* from water network, assume incorrect and choose fat network
elseif rawPred1(k,2) < 0
combinedPredictions(k,:) = pred2(k,:);    
    
else; 
end


% If extreme high PDFF from fat network, assume incorrect and make a choice on this basis 
% if rawPred2(k,1) > 1 
% combinedPredictions(k,:) = pred1(k,:);
% else ;
% end

% If extreme low PDFF from water network, assume incorrect and make a choice on this basis 

%Only do this if the relevant option is turned on 

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





