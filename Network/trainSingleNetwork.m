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
sz = 100000;

%Specify curtail factor to restrict R2* values (avoids ambiguity at higher
%R2* due to increased peak width) 
curtail = 1; %1 = no restriction of training range

% 1.1 First try uniformly spaced samples over the parameter space (vary
% both FF and R2*)

% Specify S0 (au)
S0 = 50;

%Specify SNR
% noiseSD = S0/settings.SNR;

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

%Create noise (uniform SNR)
% realnoise=(noiseSD)*randn(sz(k),numel(echotimes));
% imagnoise=1i*(noiseSD)*randn(sz(k),numel(echotimes));

%Varied SNR
%Set up matrices to allow graded SNR over the full range
snrHigh = 120;
snrLow = 20; 
snrRange = snrHigh - snrLow;
snrVec = snrLow + snrRange*rand(sz,1);
noiseSdVec = S0./snrVec; 
noiseSdMat = repmat(noiseSdVec,1,numel(echotimes));

realnoise=noiseSdMat.*randn(sz,numel(echotimes));
imagnoise=1i*noiseSdMat.*randn(sz,numel(echotimes));

noise = realnoise + imagnoise; 

%Visualise noise
figure
subplot(1,2,1)
scatter(real(noise),imag(noise))
title('Complex noise')

subplot(1,2,2)
hist(abs(noise(:,1)),20)
title('Magnitude of noise')

% Add noise to signal to create noisy signal
sCompNoisy = sNoiseFree + noise; 

%Get noise magnitude data
sMagNoisy=abs(sCompNoisy);

%Reformat complex data for DL
sCompNoisy = horzcat(real(sCompNoisy),imag(sCompNoisy));

%Choose which data to use for training
if settings.noisyTraining == 0; 
S = sMagNoiseFree;
elseif settings.noisyTraining == 1; 
S = sMagNoisy;    
else ;
end

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

% create the layers, including elu layers (consider 9 layers as per Gyori to
% optimise fitting as far as possible)

% layers = [
%     featureInputLayer(numOfFeatures, 'Name', inputName);
%     fullyConnectedLayer(numOfFeatures, 'Name', 'fc1');
%     eluLayer;
%     fullyConnectedLayer(numOfFeatures, 'Name', 'fc2');
%     eluLayer;
%     fullyConnectedLayer(numOfOutput, 'Name', 'fc3');
%     regressionLayer('Name', outputName);
%     ];
% 
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

% layers = [
%     featureInputLayer(numOfFeatures, 'Name', inputName);
%     fullyConnectedLayer(numOfFeatures, 'Name', 'fc1');
%     eluLayer;
%     fullyConnectedLayer(numOfFeatures, 'Name', 'fc2');
%     eluLayer;
%     fullyConnectedLayer(numOfFeatures, 'Name', 'fc3');
%     eluLayer;
%     fullyConnectedLayer(numOfFeatures, 'Name', 'fc4');
%     eluLayer;
%     fullyConnectedLayer(numOfFeatures, 'Name', 'fc5');
%     eluLayer;
%     fullyConnectedLayer(numOfFeatures, 'Name', 'fc6');
%     eluLayer;
%     fullyConnectedLayer(numOfOutput, 'Name', 'fc7');
%     regressionLayer('Name', outputName);
%     ];


%Wide and deep
% layers = [
%     featureInputLayer(numOfFeatures, 'Name', inputName);
%     fullyConnectedLayer(2*numOfFeatures, 'Name', 'fc1');
%     eluLayer;
%     fullyConnectedLayer(2*numOfFeatures, 'Name', 'fc2');
%     eluLayer;
%     fullyConnectedLayer(2*numOfFeatures, 'Name', 'fc3');
%     eluLayer;
%     fullyConnectedLayer(2*numOfFeatures, 'Name', 'fc4');
%     eluLayer;
%     fullyConnectedLayer(2*numOfFeatures, 'Name', 'fc5');
%     eluLayer;
%     fullyConnectedLayer(2*numOfFeatures, 'Name', 'fc6');
%     eluLayer;
%     fullyConnectedLayer(2*numOfFeatures, 'Name', 'fc7');
%     eluLayer;
%     fullyConnectedLayer(2*numOfFeatures, 'Name', 'fc8');
%     eluLayer;
%     fullyConnectedLayer(numOfOutput, 'Name', 'fc9');
%     regressionLayer('Name', outputName);
%     ];

%Superwide (3x) and superdeep (12 fully connected layers)
% layers = [
%     featureInputLayer(numOfFeatures, 'Name', inputName);
%     fullyConnectedLayer(3*numOfFeatures, 'Name', 'fc1');
%     eluLayer;
%     fullyConnectedLayer(3*numOfFeatures, 'Name', 'fc2');
%     eluLayer;
%     fullyConnectedLayer(3*numOfFeatures, 'Name', 'fc3');
%     eluLayer;
%     fullyConnectedLayer(3*numOfFeatures, 'Name', 'fc4');
%     eluLayer;
%     fullyConnectedLayer(3*numOfFeatures, 'Name', 'fc5');
%     eluLayer;
%     fullyConnectedLayer(3*numOfFeatures, 'Name', 'fc6');
%     eluLayer;
%     fullyConnectedLayer(3*numOfFeatures, 'Name', 'fc7');
%     eluLayer;
%     fullyConnectedLayer(3*numOfFeatures, 'Name', 'fc8');
%     eluLayer;
%     fullyConnectedLayer(3*numOfFeatures, 'Name', 'fc9');
%     eluLayer;
%     fullyConnectedLayer(3*numOfFeatures, 'Name', 'fc10');
%     eluLayer;
%     fullyConnectedLayer(3*numOfFeatures, 'Name', 'fc11');
%     eluLayer;
%     fullyConnectedLayer(numOfOutput, 'Name', 'fc12');
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
options = trainingOptions('adam', ...
    'MaxEpochs',50, ...
    'OutputNetwork','best-validation-loss',...
    'ValidationData', {xValidation, yValidation},...
    'InitialLearnRate',1e-3, ...
    'MiniBatchSize', 32, ...
    'L2Regularization',0,... %No regularisation as low FF values should not be preferred
    'Verbose',false, ...
    'Plots','training-progress'); 
       
%     'ValidationPatience', 50, ....

%% 5.0 Training

% Run the training
    [nets.net1,nets.info1] = trainNetwork(xTrain, yTrain, layers, options);
    
% Export duplicate of first network (simplifies code as avoids need to create a separate script for only one network)         
    nets.net2 = nets.net1; 

end
