function sigmaEst = estimateSigmaForSimulation(s0, settings, reps)
% function sigmaEst = estimateSigmaForSimulation(settings, reps)

% Generate simulations for a small number of voxels (e.g. 100) to get sigma estimate
% Finds SD across voxels (homogenous region)


%Specify ground truth parameter values
F=0;
W=s0;
v=0;

%Specify echotimes
echotimes=settings.echotimes;

%Specify true noise SD
trueNoiseSD = s0/settings.SNR;

% Initalise GT
GT=struct();
GT.p = [0 s0 0 0];

%Simulate noise-free signal
Snoisefree=MultiPeakFatSingleR2(echotimes,settings.fieldStrength,F,W,v,0);

% Specify ground truth signal
GT.S = Snoisefree;

parfor n=1:reps

%Create noise
noise=trueNoiseSD*randn(1,numel(echotimes));

%Add noise
Snoisy(n,:)=Snoisefree+noise;

end

sigmaEst = std(Snoisy(:,1));

end
