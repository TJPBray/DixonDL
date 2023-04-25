function normSignals = normaliseSignals(signals,settings)
%function normSignals = normaliseSignals(signals)

%Approximates S0 and normalises by that (a surrogate for a 'first guess'
%with deep learning fitting)

%Inputs:

%signals is n x t matrix where n is the number of instantiations and t is
%the number of echo times

%settings contains settings.echotimes, which is a t X 1 column vector containing
%the echotimes, and settings.fieldStrength

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

%% 2. Get dimensions of signals (may be matrix with multiple different signals) and loop through each set of signals

dim = size(signals,1);

for k = 1:dim

%% 2. Find s0 estimate for each set of signals

%2.1 Get signals for each row

signalRow = signals(k,:);

%2.2. Approximate S0 by interpolating between ind1 and ind2

r2estimate = -log(signalRow(ind1)/signalRow(ind2))...
    /(echotimes(ind1) - echotimes(ind2));

s0estimate = signalRow(ind1)*exp(r2estimate*echotimes(ind1));

%% 3. Normalise by s0 estimate

normSignals(k,:) = signalRow ./ s0estimate; 

end



