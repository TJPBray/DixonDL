%% Fat-water deep neural network (DNN) w/ Matlab Deep Learning Toolbox

function nets = trainNetworks(settings)
%function nets = trainNetworks(settings)

%Input:
%Settings is a structure containing TE values an field strength

%Output:
%nets is a structure containing two neural networks ('fat net' and 'water net')

%% 0.0 Specify field strength and echotimes
%Specify field strength 
tesla = settings.fieldStrength;
echotimes=settings.echotimes;

%% 1.0 Synthesise training/validation data using multipeak fat model 

%Specify random seed
rng(5)

%Specify dataset size
sz = 1000;

%Specify curtail factor to restrict R2* values (avoids ambiguity at higher
%R2* due to increased peak width) 
curtail = 1; %1 = no restriction of training range

%Specify dataset size for reduced dataset 
sz(2) = sz*curtail; 

% Create loop to enable training twice (with two different ffRanges)
for k = 1:2

% 1.1 First try uniformly spaced samples over the parameter space (vary
% both FF and R2*)

% Specify S0 (au)
S0 = 1;

%Specify SNR
SNR=60;

% 1.2 Specify ff Range

    %If using two sets of values, specify switch point
    switchPoint = 0.63;
    ffRanges = [0 switchPoint; switchPoint 1];

    %Select range depending on value of k 
    ffRange = ffRanges(k,:);
    
    %Create vector of ff values
    FFvec=ffRange(1) + (ffRange(2)-ffRange(1))*rand(sz(k),1);


% 1.3 Specify R2* range
r2max=0.5;

    %If using two sets of values, specify switch point
    r2Ranges = [0 r2max; 0 curtail*r2max]; %Restrict high FF R2* values to plausible range

    %Select range depending on value of k 
    r2Range = r2Ranges(k,:);

    %Create vector of R2* values
    R2starvec=r2Range(2)*rand(sz(k),1);

%Specify F, W and R2* values
Fvec=S0*FFvec;
Wvec=S0-Fvec;

%Concatenate vectors chosen for training
trainingParams=horzcat(FFvec,R2starvec);

%Define fB
fB = 0;

% (normalised) signal samples
sNoiseFree = MultiPeakFatSingleR2(echotimes,tesla,Fvec,Wvec,R2starvec,fB);

% Generate noise-free training data
Sreal = real(sNoiseFree);
Simag = imag(sNoiseFree);

sCompNoiseFree = horzcat(Sreal,Simag);
sMagNoiseFree = abs(sNoiseFree);

%Create noise
noiseSD=1/SNR;

realnoise=(noiseSD)*randn(sz(k),numel(echotimes));
imagnoise=1i*(noiseSD)*randn(sz(k),numel(echotimes));

% Add noise to signal to create noisy signal
sCompNoisy = sNoiseFree + realnoise + imagnoise; 

%Get noise magnitude data
sMagNoisy=abs(sCompNoisy);

%Reformat complex data for DL
sCompNoisy = horzcat(real(sCompNoisy),imag(sCompNoisy));

%Choose which data to use for training
S = sMagNoiseFree;

%% Normalise by sum of signals over echo times

Snorm = normaliseSignals(S,settings);

S = Snorm;

%% 2.0 Split synthesised data into the training and validation set
%
% This is now done with matlab's built-in cvpartition tool set.  Initially,
% used setdiff, which orders the data by default.  This ends up corrupting
% the association between the input and the output set.  This could be
% fixed with using the 'stable' option of setdiff.

%2.1 Create randomly spaced training and validation datasets

% percentage of the data to be held out as validation
hPercentage = 0.2;

% use matlab's built-in cvpartition
hPartition = cvpartition(sz(k), 'Holdout', hPercentage);

% get indices of the training and validation set
idxTrain = training(hPartition);
idxValidation = test(hPartition);

% extract the training set
xTrain = S(idxTrain,:);
yTrain = trainingParams(idxTrain,:);

% extract the validation set
xValidation = S(idxValidation,:);
yValidation = trainingParams(idxValidation,:);

% create a separate test set with values on a grid 

%% 3.0 Build a minimal DNN

% number of features
numOfFeatures = size(S,2);

% name of the input
inputName = 'Signal';

% number of output
numOfOutput = 2;

% name of the output
outputName = 'FF R2*';

% create the layers, including relu layer
layers = [
    featureInputLayer(numOfFeatures, 'Name', inputName);
    fullyConnectedLayer(numOfFeatures, 'Name', 'fc1');
    fullyConnectedLayer(numOfFeatures, 'Name', 'fc2');
    fullyConnectedLayer(numOfFeatures, 'Name', 'fc3');
    fullyConnectedLayer(numOfFeatures, 'Name', 'fc4');
    fullyConnectedLayer(numOfOutput, 'Name', 'fc5');
    regressionLayer('Name', outputName);
    ];

% layers = [
%     featureInputLayer(numOfFeatures, 'Name', inputName);
%     fullyConnectedLayer(numOfOutput, 'Name', 'fc1');
%     regressionLayer('Name', outputName);
%     ];

% number of layers
numOfLayers = size(layers, 1);

% % visualise the layers
% analyzeNetwork(layers)

%% 4.0 Set up the training options

% set up the training options with Stochastic Gradient Descent
%
% mini-batch size changed from default (128) to 64
%
% Note that Matlab implementation appears to discard the last few training
% samples that do not completely fill up a mini-batch.
%
options = trainingOptions('sgdm', ...
    'MaxEpochs',2000, ...
    'InitialLearnRate',1e-2, ...
    'MiniBatchSize', 100, ...
    'Verbose',false, ...
    'Plots','training-progress'); %No regularisation as low FF values should not be preferred

% include the validation data
options.ValidationData = {xValidation, yValidation};

%% 5.0 Training

% Run the training

% Name networks according to the loop value (net1 for low FF training, net2
% for high FF training)
if k == 1
    nets.net1 = trainNetwork(xTrain, yTrain, layers, options);

elseif k == 2
    nets.net2 = trainNetwork(xTrain, yTrain, layers, options);

else ;
end


end %End loop over FF training range values

end
