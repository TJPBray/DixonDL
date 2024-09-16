function [normSignals, s0estimates] = normaliseSignals(signals,settings)
%function [normSignals , s0estimates] = normaliseSignals(signals)

%Approximates S0 and normalises by that (a surrogate for a 'first guess'
%with deep learning fitting)

%Inputs:

%signals is m x t matrix where m is the number of instantiations and t is
%the number of echo times

%settings contains settings.echotimes, which is a t X 1 column vector containing
%the echotimes, and settings.fieldStrength

%Outputs:
%normSignals is m x t matrix of normalised signals where m is the number of instantiations and t is
%the number of echo times

%s0estimates is m x 1 vector of s0 estimates

%t.bray@ucl.ac.uk


%% 1. Find echotimes closest to inphase

echotimes = settings.echotimes;

% If 3T
if round(settings.fieldStrength) == 3

ip1Distance = echotimes - 2.3;
[min1,ind1] = min(abs(ip1Distance));

ip2Distance = echotimes - 4.6;
[min2,ind2] = min(abs(ip2Distance));

%If 1.5T
elseif round(2*settings.fieldStrength) == 3 %Rounding accounts for scanners close to but not exactly 1.5 or 3T)

ip1Distance = echotimes - 4.6;
[min1,ind1] = min(abs(ip1Distance));

ip2Distance = echotimes - 9.2;
[min2,ind2] = min(abs(ip2Distance));

else ; 
end

%% Get dimensions of signals (may be matrix with multiple different signals) and loop through each set of signals

%Set number of loops
dim = size(signals,1);

%Prefill arrays
signalRow = zeros(1,size(signals,2));
s0estimates = zeros(dim,1);

parfor k = 1:dim

%% 2. Find s0 estimate for each set of signals

%2.1 Get signals for each row

signalRow = signals(k,:);

%2.2. Approximate S0 by

% (i) extrapolating back to S0 from ind1 and ind2

r2estimate(k,:) = -log(signalRow(ind1)/signalRow(ind2))...
    /(echotimes(ind1) - echotimes(ind2));

s0estimate1 = signalRow(ind1)*exp(r2estimate(k,:)*echotimes(ind1));

% or

% (ii) taking the max of the signal vectors

s0estimate2 = max(signalRow);

%Choose the higher estimate of these two options
s0estimate3 = max(s0estimate1, s0estimate2);

s0estimates(k) = s0estimate3;

%2.3 Incorporate inaccurate S0 normalisation if this has been specified
% if isfield(settings,'normInaccuracyConstant') == 1
% s0estimates(k,1) = s0estimates(k,:)*settings.normInaccuracyConstant;
% else;
% end

%% 3. Strategy 1: Normalise by different s0 estimate for each voxel
normSignals(k,:) = signalRow ./ s0estimate3; 

end

%% 4. Strategy 2: Normalise by mean S0 (outside loop)
% Specify correction factor for particular protocol (for now just leave as
% % 1)
% correctionFactor = 1;

% %Normalise
% normSignals = signals / max(signals,[],'all'); 

end



