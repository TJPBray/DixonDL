%% Fat-water deep neural network (DNN) w/ Matlab Deep Learning Toolbox

function nets = trainSingleNetwork(settings)
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
sz = 20000;

%Specify curtail factor to restrict R2* values (avoids ambiguity at higher
%R2* due to increased peak width) 
curtail = 1; %1 = no restriction of training range

% 1.1 First try uniformly spaced samples over the parameter space (vary
% both FF and R2*)

% Specify S0 (au)
S0 = 1;

%Specify SNR
noiseSD = settings.noiseSD;

% 1.2 Specify ff Range

    %Select range depending on value of k 
    ffRange = [0 1];
    
    %Create vector of ff values
    FFvec=ffRange(1) + (ffRange(2)-ffRange(1))*rand(sz,1);

% 1.3 Specify R2* range
    r2max=0.5;
    r2Range = [0 r2max]; %Restrict high FF R2* values to plausible range

    %Create vector of R2* values
    R2starvec=r2Range(2)*rand(sz(1),1);

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
realnoise=(noiseSD)*randn(sz(1),numel(echotimes));
imagnoise=(noiseSD)*randn(sz(1),numel(echotimes));

noise = realnoise + 1i*imagnoise;

% %Visualise noise
% figure
% subplot(1,2,1)
% scatter(real(noise),imag(noise))
% title('Complex noise')
% 
% subplot(1,2,2)
% hist(abs(noise(:,1)),20)
% title('Magnitude of noise')

% Add noise to signal to create noisy signal
sCompNoisy = sNoiseFree + noise; 

%Get noise magnitude data
sMagNoisy=abs(sCompNoisy);

%Reformat complex data for DL
sCompNoisy = horzcat(real(sCompNoisy),imag(sCompNoisy));

%Choose which data to use for training
S = sMagNoiseFree;

%% Normalise signals (divide by estimated S0)

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
hPartition = cvpartition(sz(1), 'Holdout', hPercentage);

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

%% 3.0 Build a DNN

% number of features
numOfFeatures = size(S,2);

% name of the input
inputName = 'Signal';

% number of output
numOfOutput = 2;

% name of the output
outputName = 'FF R2*';

% create the layers, including relu layer
% layers = [
%     featureInputLayer(numOfFeatures, 'Name', inputName);
%     fullyConnectedLayer(numOfFeatures, 'Name', 'fc1');
%     fullyConnectedLayer(numOfFeatures, 'Name', 'fc2');
%     fullyConnectedLayer(numOfFeatures, 'Name', 'fc3');
%     fullyConnectedLayer(numOfFeatures, 'Name', 'fc4');
%     fullyConnectedLayer(numOfOutput, 'Name', 'fc5');
%     regressionLayer('Name', outputName);
%     ];

% Create the layers, including elu layer
layers = [
    featureInputLayer(numOfFeatures, 'Name', inputName);
    fullyConnectedLayer(numOfFeatures, 'Name', 'fc1');
    eluLayer;
    fullyConnectedLayer(numOfFeatures, 'Name', 'fc2');
    eluLayer;
    fullyConnectedLayer(numOfFeatures, 'Name', 'fc3');
    eluLayer;
    fullyConnectedLayer(numOfFeatures, 'Name', 'fc4');
    eluLayer;
    fullyConnectedLayer(numOfOutput, 'Name', 'fc5');
    regressionLayer('Name', outputName);
    ];

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
options = trainingOptions('adam', ...
    'MaxEpochs',1000, ...
    'InitialLearnRate',1e-3, ...
    'MiniBatchSize', 32, ...
    'ValidationPatience', 20, ....
    'L2Regularization',0,...
    'Verbose',false, ...
    'Plots','training-progress'); %No regularisation as low FF values should not be preferred

% include the validation data
options.ValidationData = {xValidation, yValidation};

%% 5.0 Training

% Run the training
    nets.net1 = trainNetwork(xTrain, yTrain, layers, options);
    
% Export duplicate of first network (simplifies code as avoids need to create a separate script for only one network)         
    nets.net2 = nets.net1; 

end
