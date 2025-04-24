function [sseVec,likVec] = sseVecCalc (echotimes, tesla, predictionVec, S, sigma)
% function sseVec, likVec = sseVecCalc (echotimes, tesla, predictionVec, S)

%Calculates a vector of SSE values from the predictions of neural network

% Input:
%

% Output:
% likVec is an output of likelihood values based on the above

%% 1. Get n
n=size(predictionVec,1);

% Prefill arrays
p = zeros(4,1);


%% 2. Loop implementation over n values (slower)

% parfor k=1:n
% 
%     %2.1 First specify components of p
%     f=predictionVec(k,1)*predictionVec(k,3);
%     w=(1-predictionVec(k,1))*predictionVec(k,3);
%     v=predictionVec(k,2);
%     fB=0;
% 
%     p = [f; w; v; fB];
%     p = double(p);
% 
%     Smeasured=S(k,:);
% 
%     %2.2 Calculate likelihood (use Gaussian initially)
%     [gaussianLikVec(k),sseVec(k)] = R2Obj(p,echotimes,tesla,Smeasured,sigma); %Use sigma of 0 - note that likelilihood values with be NaN
% 
%     %2.2 Calculate likelihood (Rician)
%     [ricianLikVec(k)] = R2RicianObj(p,echotimes,tesla,Smeasured,sigma); %Use sigma of 0 - note that likelilihood values with be NaN
% 
% end

%% 2. Vectorised implementation (faster)

    %2.1 First specify components of p
    f=predictionVec(:,1).*predictionVec(:,3);
    w=(1-predictionVec(:,1)).*predictionVec(:,3);
    v=predictionVec(:,2);
    fB=zeros(size(v));

    p = [f w v fB]';
    p = double(p);
  
    Smeasured=S;

    %2.2 Calculate likelihood (use Gaussian initially)
    [gaussianLikVec,sseVec] = R2Obj(p,echotimes,tesla,Smeasured,sigma); %Use sigma of 0 - note that likelilihood values with be NaN

    gaussianLikVec = gaussianLikVec';
    sseVec = sseVec';

    %2.2 Calculate likelihood (Rician)
    [ricianLikVec] = R2RicianObj(p,echotimes,tesla,Smeasured,sigma); %Use sigma of 0 - note that likelilihood values with be NaN

    ricianLikVec = ricianLikVec';

%%    

likVec = ricianLikVec;

end

