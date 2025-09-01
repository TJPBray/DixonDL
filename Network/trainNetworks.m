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

%Create vector of ff values for full training distribution
FFvecAll=ffRange(1) + (ffRange(2)-ffRange(1))*rand(sz,1);

% 1.3 Specify R2* range
r2max=0.5;
r2Range = [0 r2max]; %Restrict high FF R2* values to plausible range

%Create vector of R2* values
R2starvecAll=r2Range(2)*rand(sz(1),1);

%Specify F, W and R2* values
Fvec=S0*FFvecAll;
Wvec=S0-Fvec;

%Define fB
fB = 0;

% (normalised) signal samples
sNoiseFree = MultiPeakFatSingleR2(echotimes,tesla,Fvec,Wvec,R2starvecAll,fB);

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

sAll = S;

%% Create loop to enable training twice (with two different ffRanges)
for k = 1:2

    % 1.1 Specify switch point
    switchPoint = 0.58;

    % 1.2 Create vector of ff values for relevant region of paramter space

    if k == 1
        FFvec = FFvecAll(FFvecAll<=switchPoint);
        R2starvec = R2starvecAll(FFvecAll<=switchPoint);
        S = sAll((FFvecAll<=switchPoint),:);

    elseif k == 2

        FFvec = FFvecAll(FFvecAll>switchPoint);
        R2starvec = R2starvecAll(FFvecAll>switchPoint);
        S = sAll((FFvecAll>switchPoint),:);

    else ;
    end

    figure, scatter(FFvec,R2starvec)
    xlim([0 1])
    ylim([0 .5])

    %1.3 Specify F, W and R2* values
    Fvec=S0*FFvec;
    Wvec=S0-Fvec;

    %Concatenate vectors chosen for training
    trainingParams=horzcat(FFvec,R2starvec);


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
    hPartition = cvpartition(size(FFvec,1), 'Holdout', hPercentage);

    % get indices of the training and validation set
    idxTrain = training(hPartition);
    idxValidation = test(hPartition);

    % extract the training set
    xTrain = S(idxTrain,:);
    yTrain = trainingParams(idxTrain,:);

    % extract the validation set
    xValidation = S(idxValidation,:);
    yValidation = trainingParams(idxValidation,:);

    %% 3.0 Build a minimal DNN

    % number of features
    numOfFeatures = size(S,2);

    % name of the input
    inputName = 'Signal';

    % number of output
    numOfOutput = 2;

    % name of the output
    outputName = 'FF R2*';

    % create the layers, including elu layers
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
        ];

    % number of layers
    numOfLayers = size(layers, 1);

    % % visualise the layers
    % analyzeNetwork(layers)

    % create a network from the specified layers and initialise any unset learnable and state parameters
    if k == 1
    nets.net1 = dlnetwork(layers)
    elseif k==2
    nets.net2 = dlnetwork(layers)
    else ; 
    end

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
        'ValidationPatience', 50, ....
        'OutputNetwork','best-validation-loss',...
        'ValidationData', {xValidation, yValidation},...
        'InitialLearnRate',1e-3, ...
        'MiniBatchSize', 32, ... 
        'L2Regularization',0,...
        'Verbose',false, ...
        'Plots','training-progress');   %No regularisation as low FF values should not be preferred

    %% 5.0 Training

    % Run the training

    % Name networks according to the loop value (net1 for low FF training, net2
    % for high FF training)
    if k == 1
        [nets.net1,nets.info1] = trainnet(xTrain, yTrain, nets.net1, 'mse', options);

    elseif k == 2
        [nets.net2,nets.info2] = trainnet(xTrain, yTrain, nets.net2, 'mse', options);

    else ;
    end


end %End loop over FF training range values

end
