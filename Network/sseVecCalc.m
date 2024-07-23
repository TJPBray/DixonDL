function [likVec,sseVec] = sseVecCalc (echotimes, tesla, predictionVec, S, sigma) 
% function likVec,sseVec = sseVecCalc (echotimes, tesla, predictionVec, S) 

%Calculates a vector of SSE values from the predictions of neural network

% Input:
% 

% Output:
% likVec is an output of likelihood values based on the above

% 1. Get n
n=size(predictionVec,1);

%2. Loop over n values
for k=1:n
    
    %2.1 First specify components of p
    p(1)=predictionVec(k,1);
    p(2)=1-p(1);
    p(3)=predictionVec(k,2);
    p(4)=0;

Smeasured=S(k,:);

    %2.2 Calculate likelihood (use Gaussian initially)
    [gaussianLikVec(k),sseVec(k)] = R2Obj(p,echotimes,tesla,Smeasured,sigma); %Use sigma of 0 - note that likelilihood values with be NaN
   
    %2.2 Calculate likelihood (Rician)
     %2.2 Calculate likelihood (use Gaussian initially)
    [ricianLikVec(k),sseVec(k)] = R2Obj(p,echotimes,tesla,Smeasured,sigma); %Use sigma of 0 - note that likelilihood values with be NaN

end

likVec = ricianLikVec;

end

