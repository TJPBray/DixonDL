% Master script
% Runs the various simulations and generates figures
% Author: Tim Bray (t.bray@ucl.ac.uk)

%% Set user path
% userpath('/Users/tjb57/Dropbox/MATLAB/');
up = userpath;

%% 1. Simulations

% Specify settings
% settings.echotimes = [1.1:1.1:13.2]';
load(strcat(up,'/Fat-water MAGORINO/Data/Hernando data/site1_begin_3T_protocol1.mat'))
settings.echotimes = (1000*imDataAll.TE);

settings.fieldStrength = 3;

%Specify tuning parameters for inference (choice of network output)
settings.tuning.useImplausibleValues = 0;
settings.tuning.clipOutputs = 0;

%% 1.1 First single network only
% Note: trainSingleNetworks actually returns two networks but one is a duplicate of the first - simplifies code by avoiding need for separate single-network script)

%Train in the absence of noise
settings.noisyTraining = 0;
net = trainSingleNetwork(settings);
trainingCurves(net,'Noise free training')

%Train in the presence of noise
settings.noisyTraining = 1;
net = trainSingleNetwork(settings);
trainingCurves(net,'Noisy training')

%Test on simulation data (noise-free)
settings.SNR = inf;
settings.noisyTesting = 0;
[dlMaps,dlErrormaps,dlSDMaps] = testOnSimulatedData(net,settings)
meanFfError = mean(abs(dlErrormaps.ff),'all')
meanR2Error = mean(abs(dlErrormaps.R2),'all')

%Test on simulation data (with noise)
settings.SNR = 60;
settings.noisyTesting = 1;
[dlMaps,dlErrormaps,dlSDMaps] = testOnSimulatedData(net,settings)
meanFfError = mean(abs(dlErrormaps.ff),'all')
meanR2Error = mean(abs(dlErrormaps.R2),'all')

%Show comparison vs conventional fitting
%Load conventional fitting data for comparison
load(strcat(up,'/Fat-water MAGORINO DL/DixonDL/Results/SimulationResults/6 echoes/MAGORINO_Simulation_Results_SNR60_R100_fineGrain_6echoes.mat'))

%Visualise comparison against conventional fitting
createFigDLvsConventionalFitting(FFmaps, R2maps, dlMaps, errormaps, dlErrormaps, dlSDMaps, sdmaps)

%% 1.2 Dual networks

%Train in the absence of noise
settings.noisyTraining = 0;
nets = trainNetworks(settings);
trainingCurves(nets)

%Train in the presence of noise
settings.noisyTraining = 1;
nets = trainNetworks(settings);
trainingCurves(nets)

%Test on simulation data (noise-free)
settings.SNR = inf;
settings.noisyTesting = 0;
[dlMaps,dlErrormaps,dlSDMaps] = testOnSimulatedData(nets,settings)
meanFfError = mean(abs(dlErrormaps.ff),'all')
meanR2Error = mean(abs(dlErrormaps.R2),'all')

%Test on simulation data (with noise)
settings.SNR = 60;
settings.noisyTesting = 1;
[dlMaps,dlErrormaps,dlSDMaps] = testOnSimulatedData(nets,settings)
meanFfError = mean(abs(dlErrormaps.ff),'all')
meanR2Error = mean(abs(dlErrormaps.R2),'all')

%Show comparison vs conventional fitting
%Load conventional fitting data for comparison
load(strcat(up,'Fat-water MAGORINO DL/DixonDL/Results/SimulationResults/6 echoes/MAGORINO_Simulation_Results_SNR60_R100_fineGrain_12echoes_withComplex.mat'))

%Visualise comparison against conventional fitting
createFigDLvsConventionalFitting(FFmaps, R2maps, dlMaps, errormaps, dlErrormaps, dlSDMaps,sdmaps)

% %% 1.2 Dual networks, investigation of normalisation
%
% %Specify inaccurate normalisation (normInaccuracy is multipled by the estimated S0 prior to normalisation: 0.9 means S0 estimate is too low , 1.1
% %means S0 estimate is too high)
% settings.normInaccuracyConstant = 1.1;
%
% %Test on simulation data (noise-free)
% [dlMaps,dlErrormaps,dlSDMaps] = testOnSimulatedData(nets,settings)
%
% %Test on simulation data (with noise)
% settings.SNR = 60;
% settings.noiseSD = 1/settings.SNR;
%
% [dlMaps,dlErrormaps,dlSDMaps] = testOnSimulatedData(nets,settings)
%
% %Show comparison vs conventional fitting
% %Load conventional fitting data for comparison
% load('/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO DL/DixonDL/Results/SimulationResults/MAGORINO_Simulation_Results_SNR60_R100_fineGrain.mat')
%
% %Visualise comparison against conventional fitting
% createFigDLvsConventionalFitting(FFmaps, dlMaps, errormaps, dlErrormaps, dlSDMaps,sdmaps)



% %% 1.3 Test accuracy of normalisation
%
% % Load chosen dataset to get settings for that dataset
% load('site1_begin_3T_protocol2.mat')
%
% % Get settings
% settings.echotimes = 1000*imDataAll.TE;
% settings.fieldStrength = imDataAll.FieldStrength;
%
% % Provide guess for SNR and noiseSD
% SNRguess = 50;
% % settings.SNR = SNRguess;
% % settings.noiseSD = 1/SNRguess;
% settings.noiseSD = 0;
%
% %Test normalisation
% testNormalisation(settings)


%% 2. Train and implement network with echotimes corresponding to specific SUBJECT dataset
clear all

%2.1 Load chosen dataset

% Example cases:
load(strcat(up,'/Fat-water MAGORINO/Data/Subjects/FW101_bipolar.mat'));
load(strcat(up,'/Fat-water MAGORINO/Data/Subjects/FW101_bipolar_sigmaRoi.mat'));
load(strcat(up,'/Fat-water MAGORINO/Experiments/experimentResults/Rev2_FW101_results')); prc = 60; yDisp = (20:260);
slice = 20; sigmaEst = 33.6073;
%
load(strcat(up,'/Fat-water MAGORINO/Data/Subjects/LegsSwapData.mat'));
load(strcat(up,'/Fat-water MAGORINO/Data/Subjects/LegsSwapData_sigmaRoi.mat'));
load(strcat(up,'/Fat-water MAGORINO/Experiments/experimentResults/Rev2_LegsSwapData_results.mat')); prc = 86; yDisp = (1:170);
slice = 20; sigmaEst = 5.3559;

load(strcat(up,'/Fat-water MAGORINO/Data/Subjects/wbMriAbdomen.mat'));
load(strcat(up,'/Fat-water MAGORINO/Data/Subjects/wbMriAbdomen_sigmaRoi.mat')); 
load(strcat(up,'/Fat-water MAGORINO/Experiments/experimentResults/wbMRI_abdoPelvis.mat')); prc = 86; yDisp = (1:170);prc = 65; yDisp = (1:192);
slice = 40; sigmaEst = 3.4585;
%
load(strcat(up,'/Fat-water MAGORINO/Data/Subjects/wbMriThorax.mat'));
load(strcat(up,'/Fat-water MAGORINO/Data/Subjects/wbMriThorax_sigmaRoi.mat'));
load(strcat(up,'/Fat-water MAGORINO/Experiments/experimentResults/wbmri_thorax.mat')); prc = 68; yDisp = (1:192);
slice = 40; sigmaEst = 5.5703;

% Others:
% load('/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO/Data/Subjects/FW111_monopolar.mat')

%Convert to imData if not already
try
    imDataAll = imDataParams;
catch
end

%Mask out background to avoid unwanted 'noise' in plot
mask = zeros (size(imDataAll.images(:,:,slice,1,1)));
thresh = prctile(abs(imDataAll.images(:,:,:,1,1)),prc,'all');
mask(abs(imDataAll.images(:,:,slice,1,1))>thresh) = 1;
mask = imfill(mask,"holes");
figure, imshow(mask)

%Convert ROI to single structure
roi.mask=BW{1};
roi.slice=slice;
roi.size=sz_vol;

%Get sigma estimate
% imDataAll.fittingIndent=0; %Specify no indent
% [mapsWithSigma,sigmaFromRoi,sigmaFromFit] = FitImageSigma(imDataAll,roi);

%2.2 Get settings
settings.echotimes = 1000*imDataAll.TE';
settings.fieldStrength = imDataAll.FieldStrength;

%2.3 Provide estimate of noiseSD
settings.sigmaEst = sigmaEst;

%2.4 Train networks with settings (echotimes, fieldstrength) for chosen
%dataset
settings.noisyTraining = 1;
nets = trainNetworks(settings);
net = trainSingleNetwork(settings);

%Specify tuning parameters for inference (choice of network output)
settings.tuning.useImplausibleValues = 1;
settings.tuning.clipOutputs = 0;

%2.5 Check network performance on simulation data
% testOnSimulatedData(nets,settings)

%2.6 Use networks to get predicts for the chosen dataset

%Get the signals for the chosen slice
for echoN=1:numel(settings.echotimes)
    complexSlice(:,:,echoN)=imDataAll.images(:,:,slice,1,echoN);
    magnitudeSlice(:,:,echoN)=abs(complexSlice(:,:,echoN));

end

% %Remove extreme outlying values
% threshold = prctile(magnitudeSlice,95,'all');
% for y = 1:(size(magnitudeSlice))
%     for x = 1:(size(magnitudeSlice))
%
%         values = magnitudeSlice(y,x,:);
%         if max(values,[],'all') > threshold
%             magnitudeSlice(y,x,:) = values .* (threshold/max(values,[],'all'));
%         else ;
%         end
%     end
% end

% Concatenate for time calculation (does not necessarily need to time
% saving compared to a loop with a CPU)
sliceN = 1;
repSlice = repmat(magnitudeSlice,[sliceN 1 1]);

%Call the function
tic
predictions = testOnImage(magnitudeSlice,nets,settings);
toc

%2.8 Get maps
ffLow = predictions.prediction1(:,:,1);
ffHigh = predictions.prediction2(:,:,1);
ffMap = predictions.prediction3(:,:,1);

r2Low = predictions.prediction1(:,:,2);
r2High = predictions.prediction2(:,:,2);
r2Map = predictions.prediction3(:,:,2);

%Remove NaNs
ffMap(isnan(ffMap)) = 0;
r2Map(isnan(r2Map)) = 0;

%2.7 Display maps

%Mask out background to avoid unwanted 'noise' in plot
mask = zeros (size(imDataAll.images(:,:,slice,1,1)));
thresh = prctile(abs(imDataAll.images(:,:,:,1,1)),prc,'all');
mask(abs(imDataAll.images(:,:,slice,1,1))>thresh) = 1;
mask = imfill(mask,"holes");
figure, imshow(mask)

%Create plots
figure
subplot(2,4,1)
s1=imshow(ffLow(yDisp,:).*mask(yDisp,:),[0 1])
set(s1,'Alphadata',mask(yDisp,:))
title('PDFF - Water network output')
colormap(gca,'parula')
colorbar

subplot(2,4,2)
s2=imshow(ffHigh(yDisp,:).*mask(yDisp,:),[0 1])
set(s2,'Alphadata',mask(yDisp,:))
title('PDFF - Fat network output')
colormap(gca,'parula')
colorbar

subplot(2,4,3)
s3=imshow(ffMap(yDisp,:).*mask(yDisp,:),[0 1])
set(s3,'Alphadata',mask(yDisp,:))
title('PDFF - Likelihood-chosen output')
colormap(gca,'parula')
colorbar

subplot(2,4,5)
s5=imshow(1000*r2Low(yDisp,:).*mask(yDisp,:),[0 200])
set(s5,'Alphadata',mask(yDisp,:))
title('R2* - Water network output')
colormap(gca,'copper')
colorbar

subplot(2,4,6)
s6=imshow(1000*r2High(yDisp,:).*mask(yDisp,:),[0 200])
set(s6,'Alphadata',mask(yDisp,:))
title('R2* - Fat network output')
colormap(gca,'copper')
colorbar

subplot(2,4,7)
s7=imshow(1000*r2Map(yDisp,:).*mask(yDisp,:),[0 200])
set(s7,'Alphadata',mask(yDisp,:))
title('R2* - Likelihood-chosen output')
colormap(gca,'copper')
colorbar

try
    subplot(2,4,4)
    s4 = imshow(maps.FFrician(yDisp,:).*mask(yDisp,:),[0 1],'InitialMagnification', 400)
    set(s4,'Alphadata',mask(yDisp,:))
    colormap(gca, 'parula')
    title('PDFF - Conventional fitting')
    h=colorbar
    h.Label.String = "PDFF";
    h.FontSize=12;

    subplot(2,4,8)
    s8 = imshow(1000*maps.R2rician(yDisp,:).*mask(yDisp,:),[0 200])
    set(s8,'Alphadata',mask(yDisp,:))
    colormap(gca, 'copper')
    title('R2* - Conventional fitting')
    h=colorbar
    h.Label.String = "R2* (s^-^1)";
    h.FontSize=12;

catch
end



% %Interrogate chosen voxel
% xcoord = 100; ycoord = 160;
%
% voxelSignal = magnitudeSlice(ycoord,xcoord,:);
%
% figure, plot(settings.echotimes,squeeze(voxelSignal))
% ylabel('Signal')
% xlabel('TE')
%
% % Get prediction for voxel
% predictions = testOnImage(voxelSignal,nets,settings)
%

% 2.7 Run conventional fitting for comparison

[FFmaps,errormaps,sdmaps,residuals] = Simulate_Values(50, 0, 20)



%% Display comparison with conventional fitting


% mask = ones(320,320);

%Show images comparing the two
figure
subplot(3,2,1)
s1 = imshow(maps.FFrician(yDisp,:).*mask(yDisp,:),[0 1],'InitialMagnification', 400)
set(s1,'Alphadata',mask(yDisp,:))
colormap(gca, 'parula')
title('PDFF - Conventional fitting')
h=colorbar
h.Label.String = "PDFF";
h.FontSize=12;

subplot(3,2,3)
s2 = imshow(ffMap(yDisp,:).*mask(yDisp,:),[0 1])
set(s2,'Alphadata',mask(yDisp,:))
colormap(gca, 'parula')
title('PDFF - DL')
h=colorbar
h.Label.String = "PDFF";
h.FontSize=12;

subplot(3,2,5)
s3 = imshow(ffMap(yDisp,:).*mask(yDisp,:) - maps.FFrician(yDisp,:).*mask(yDisp,:),[-0.2 0.2])
set(s3,'Alphadata',mask(yDisp,:))
colormap(gca, 'parula')
title('PDFF difference (DL - conventional)')
h=colorbar
h.Label.String = "PDFF difference";
h.FontSize=12;

% subplot(3,2,7)
% scatter(reshape(maps.FFrician.*mask,1,[]),reshape(ffMap.*mask,1,[]),'filled','MarkerFaceAlpha',0.1)
% title('PDFF difference (DL - conventional)')
% xlim([0 1])
% ylim([0 1])
% xlabel('MAGORINO FF')
% ylabel('RAIDER FF')

% Do the same for R2*
subplot(3,2,2)
s4 = imshow(1000*maps.R2rician(yDisp,:).*mask(yDisp,:),[0 500])
set(s4,'Alphadata',mask(yDisp,:))
colormap(gca,'copper')
title('R2* - Conventional fitting')
h=colorbar
h.Label.String = "R2* (s^-^1)";
h.FontSize=12;

subplot(3,2,4)
s5 = imshow(1000*r2Map(yDisp,:).*mask(yDisp,:),[0 500])
set(s5,'Alphadata',mask(yDisp,:))
colormap(gca, 'copper')
title('R2* - DL')
h=colorbar
h.Label.String = "R2* (s^-^1)";
h.FontSize=12;

subplot(3,2,6)
s6 = imshow(1000*r2Map(yDisp,:).*mask(yDisp,:) - 1000*maps.R2rician(yDisp,:).*mask(yDisp,:),[-50 50])
set(s6,'Alphadata',mask(yDisp,:))
colormap(gca, 'copper')
title('R2* difference (DL - conventional)')
h=colorbar
h.Label.String = "R2* difference (s^-^1)";
h.FontSize=12;

% Extract values for scatterplots

%Creat or load ROIs
% newanal2(ffMap)
load('/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO DL/DixonDL/Results/SubjectResults/wbMriThorax_ROIS.mat')
load('/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO DL/DixonDL/Results/SubjectResults/legsData_ROIs.mat')

%Define reference maps
MAGORINOslice = maps.FFrician;
MAGORINOR2slice = maps.R2rician;

%Extract ROI values for different organs

%% For thorax
%Liver
roiMask1 = BW{1};
ffLiver = reshape(ffMap(roiMask1==1),1,[]);
ffLiverRef = reshape(MAGORINOslice(roiMask1==1),1,[]);
r2Liver = reshape(r2Map(roiMask1==1),1,[]);
r2LiverRef = reshape(MAGORINOR2slice(roiMask1==1),1,[]);


%Kidney
roiMask2 = BW{2} + BW{3};
ffKidney = reshape(ffMap(roiMask2==1),1,[]);
ffKidneyRef = reshape(MAGORINOslice(roiMask2==1),1,[]);
r2Kidney = reshape(r2Map(roiMask2==1),1,[]);
r2KidneyRef = reshape(MAGORINOR2slice(roiMask2==1),1,[]);

%Bone
roiMask3 = BW{4} + BW{5} + BW{6};
ffBone = reshape(ffMap(roiMask3==1),1,[]);
ffBoneRef = reshape(MAGORINOslice(roiMask3==1),1,[]);
r2Bone = reshape(r2Map(roiMask3==1),1,[]);
r2BoneRef = reshape(MAGORINOR2slice(roiMask3==1),1,[]);

%Visceral fat
roiMask4 = BW{7};
ffFat = reshape(ffMap(roiMask4==1),1,[]);
ffFatRef = reshape(MAGORINOslice(roiMask4==1),1,[]);
r2Fat = reshape(r2Map(roiMask4==1),1,[]);
r2FatRef = reshape(MAGORINOR2slice(roiMask4==1),1,[]);


%Scatterplots of FF and R2* (RAIDER vs MAGORINO)
figure
scatter(ffFatRef,ffFat,60,'filled','MarkerFaceAlpha',.3)
hold on
scatter(ffKidneyRef,ffKidney,60,'filled','MarkerFaceAlpha',.3)
scatter(ffBoneRef,ffBone,60,'filled','MarkerFaceAlpha',.3)
scatter(ffLiverRef,ffLiver,60,'filled','MarkerFaceAlpha',.3)
plot([-.1 1.1], [-.1 1.1],'k')
hold off
xlim([0 1])
ylim([0 1])
legend('Visceral fat','Kidney','Bone','Liver')
xlabel('MAGORINO PDFF','FontSize',12)
ylabel('RAIDER PDFF','FontSize',12)


figure
scatter(r2FatRef,r2Fat,60,'filled','MarkerFaceAlpha',.3)
hold on
scatter(r2KidneyRef,r2Kidney,60,'filled','MarkerFaceAlpha',.3)
scatter(r2BoneRef,r2Bone,60,'filled','MarkerFaceAlpha',.3)
scatter(r2LiverRef,r2Liver,60,'filled','MarkerFaceAlpha',.3)
plot([-.1 1.1], [-.1 1.1],'k')
hold off
xlim([0 0.5])
ylim([0 0.5])
xticks([0 .1 .2 .3 .4 .5])
yticks([0 .1 .2 .3 .4 .5])
xticklabels([0 100 200 300 400 500])
yticklabels([0 100 200 300 400 500])
legend('Visceral fat','Kidney','Bone','Liver')
xlabel('MAGORINO R_2* (s^-^1)','FontSize',12)
ylabel('RAIDER R_2* (s^-^1)','FontSize',12)

% Combine values for regression

ffRefAllTissues = horzcat(ffLiverRef,ffKidneyRef, ffBoneRef,ffFatRef);
ffAllTissues = horzcat(ffLiver,ffKidney, ffBone,ffFat);
mdl = fitlm(ffRefAllTissues,ffAllTissues)
figure, scatter(ffRefAllTissues,ffAllTissues)
hold on
plot([0 1],[0 1])
hold off
%
figure

ffVals = reshape(ffMap(mask==1),1,[]);
MAGORINOslice = maps.FFrician;
ffRefVals = reshape(MAGORINOslice(mask==1),1,[]);
MAGORINOR2slice = maps.R2rician;
r2RefVals = reshape(MAGORINOR2slice(mask==1),1,[]);
ffDifference = ffVals - ffRefVals;

% Make scatterplots
figure
subplot(1,2,1)
scatter(ffRefVals,ffVals,'filled','MarkerFaceAlpha',0.01)
hold on
plot([-.1 1.1], [-.1 1.1],'k')
xlim([-.1 1.1])
ylim([-.1 1.1])

subplot(1,2,2)
scatter3(r2RefVals,ffRefVals,ffDifference,15,ffDifference,'o','filled','MarkerFaceAlpha',0.5)
ylim([0 1])
xlim([0 0.5])
ylabel('FF')
xlabel('R2*')
zlabel('FF difference')
view(0,90)
ax=gca;
ax.CLim=[-1 1];
set(gca, 'YDir','reverse')
colorbar

%% For legs
%Liver
roiMask1 = BW{1};
ffMuscle = reshape(ffMap(roiMask1==1),1,[]);
ffMuscleRef = reshape(MAGORINOslice(roiMask1==1),1,[]);
r2Muscle = reshape(r2Map(roiMask1==1),1,[]);
r2MuscleRef = reshape(MAGORINOR2slice(roiMask1==1),1,[]);

%Bone
roiMask2 = BW{2};
ffBone = reshape(ffMap(roiMask2==1),1,[]);
ffBoneRef = reshape(MAGORINOslice(roiMask2==1),1,[]);
r2Bone = reshape(r2Map(roiMask2==1),1,[]);
r2BoneRef = reshape(MAGORINOR2slice(roiMask2==1),1,[]);

%Subcutaneous fat
roiMask3 = BW{3};
ffFat = reshape(ffMap(roiMask3==1),1,[]);
ffFatRef = reshape(MAGORINOslice(roiMask3==1),1,[]);
r2Fat = reshape(r2Map(roiMask3==1),1,[]);
r2FatRef = reshape(MAGORINOR2slice(roiMask3==1),1,[]);


%Scatterplots of FF and R2* (RAIDER vs MAGORINO)
figure
scatter(ffMuscleRef,ffMuscle,60,'filled','MarkerFaceAlpha',.3)
hold on
scatter(ffBoneRef,ffBone,60,'filled','MarkerFaceAlpha',.3)
scatter(ffFatRef,ffFat,60,'filled','MarkerFaceAlpha',.3)
hold off
xlim([0 1])
ylim([0 1])
legend('Muscle','Bone','Subcutaneous fat')
xlabel('MAGORINO PDFF','FontSize',12)
ylabel('RAIDER PDFF','FontSize',12)


figure
scatter(r2MuscleRef,r2Muscle,60,'filled','MarkerFaceAlpha',.3)
hold on
scatter(r2BoneRef,r2Bone,60,'filled','MarkerFaceAlpha',.3)
scatter(r2FatRef,r2Fat,60,'filled','MarkerFaceAlpha',.3)
hold off
xlim([0 0.5])
ylim([0 0.5])
xticks([0 .1 .2 .3 .4 .5])
yticks([0 .1 .2 .3 .4 .5])
xticklabels([0 100 200 300 400 500])
yticklabels([0 100 200 300 400 500])
legend('Muscle','Bone','Subcutaneous fat')
xlabel('MAGORINO R_2* (s^-^1)','FontSize',12)
ylabel('RAIDER R_2* (s^-^1)','FontSize',12)


%% Display comparison with vendor

load('/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO DL/DixonDL/Results/SubjectResults/Old/wbMriThorax.mat', 'ffMDQ')
load('/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO DL/DixonDL/Results/SubjectResults/Old/wbMriThorax.mat', 'r2MDQ')

%Mask out background
mask = zeros (size(magnitudeSlice(:,:,1)));
thresh = prctile(abs(magnitudeSlice),55,'all');
mask(magnitudeSlice(:,:,1)>thresh) = 1;
mask = imfill(mask,"holes");
figure, imshow(mask)

% Plot
figure
p1 = subplot(2,3,1)
imshow(ffMap(10:182,10:296),[0 1])
title('PDFF - DL')
h=colorbar
colormap(p1,'parula')
h.Label.String = "PDFF";
h.FontSize=12;

p2 = subplot(2,3,4)
imshow(1000*r2Map(10:182,10:296),[0 200])
title('R2* DL')
h=colorbar
h.Label.String = "R2* (s^-^1)";
h.FontSize=12;

p3 = subplot(2,3,2)
imshow(ffMDQ(10:182,10:296,slice)*0.001,[0 1])
title('PDFF - vendor')
h=colorbar
colormap(p3, 'parula')
h.Label.String = "PDFF";
h.FontSize=12;

p4 = subplot(2,3,5)
imshow(r2MDQ(10:182,10:296,slice),[0 200])
title('R2* - vendor')
h=colorbar
h.Label.String = "R2* (s^-^1)";
h.FontSize=12;

p5 = subplot(2,3,3)
imshow(ffMap(10:182,10:296) - ffMDQ(10:182,10:296,slice)*0.001,[-1 1])
title('PDFF difference (DL - vendor)')
h=colorbar
h.Label.String = "PDFF difference";
h.FontSize=12;


p6 = subplot(2,3,6)
imshow(1000*r2Map(10:182,10:296) - r2MDQ(10:182,10:296,slice),[-200 200])
title('R2* difference (DL - vendor)')
h=colorbar
h.Label.String = "R2* difference (s^-^1)";
h.FontSize=12;

% Show scatterplot
ffVals = reshape(ffMap(mask==1),1,[]);
MDQFFslice = ffMDQ(:,:,slice);
ffRefVals = 0.001*reshape(MDQFFslice(mask==1),1,[]);
MDQR2slice = r2MDQ(:,:,slice);
r2RefVals = 0.001*reshape(MDQR2slice(mask==1),1,[]);
ffDifference = ffVals - ffRefVals;


figure
subplot(1,2,1)
scatter(ffRefVals,ffVals,'filled','MarkerFaceAlpha',0.01)
hold on
plot([-.1 1.1], [-.1 1.1],'k')
xlim([-.1 1.1])
ylim([-.1 1.1])

subplot(1,2,2)
scatter3(r2RefVals,ffRefVals,ffDifference,15,ffDifference,'o','filled')
ylim([0 1])
xlim([0 0.5])
ylabel('FF')
xlabel('R2*')
zlabel('FF difference')
view(0,90)
ax=gca;
ax.CLim=[-1 1];
set(gca, 'YDir','reverse')
colorbar

%% 3A. TRAIN network with echotimes corresponding to PHANTOM dataset

%3.1 Define folders for import of multiecho data and ROIs
imageFolder='/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO/Data/Hernando data';
roiFolder='/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO/Data/Hernando ROIs';

%Folder for saving
saveFolder='/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO DL/DixonDL/Results/HernandoPhantomNetworks';

load('/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO DL/DixonDL/Results/sigmaEstimates.mat')

%3.2 Get image and ROI folder info
imageFolderInfo=dir(imageFolder);
roiFolderInfo=dir(roiFolder);
saveFolderInfo=dir(saveFolder);

%3.3. Specify reference values in phantoms
% Create grid of reference values based on phantom structure
ReferenceValues.FF = ([0; 0.026; 0.053; 0.079; 0.105; 0.157; 0.209; 0.312; 0.413; 0.514; 1]);
% ReferenceValues.R2 = repelem([0; 0; 0; 0],1,5);

%3.4 Specify slice (use 1 if only 1 slice)
sl=2;

%3.5 Loop over datasets to fit each dataset
for n=1:28

    %Clear unwanted variables
    clear complexSlice magnitudeSlice

    %3.5 Define filenames for multiecho data and ROIs
    dataFileName = imageFolderInfo(n+2).name
    roiFileName = roiFolderInfo(n+2).name;

    %3.6 Import data
    struct=load(fullfile(imageFolder,dataFileName));
    imData=struct.imDataAll;
    fwmc_ff=struct.fwmc_ff;
    fwmc_r2star=struct.fwmc_r2star;

    %3.7 Import ROIs
    mask = niftiread(fullfile(roiFolder,roiFileName));
    roi.mask=mask(:,:,sl);
    roi.slice=sl;

    % phantomImage = niftiread(fullfile(folder,dataFileName));

    % % Show check image for alignment
    % figure
    % subplot(2,3,1)
    % imshow(phantomROIs(:,:,1),[])
    % title('ROIs')
    %
    % subplot(2,3,2)
    % imshow(abs(imData.images(:,:,sl,1,1)),[])
    % title('Magnitude image of phantom (first echo)')
    %
    % subplot(2,3,3)
    % imshow(fwmc_ff(:,:,2),[0 100])
    % title('PDFF map of phantom')
    %
    % subplot(2,3,4)
    % imshow(phantomROIs(:,:,sl),[0 12])
    % title('Phantom ROIs')

    % Perform fitting of multiecho data, using sigma estimate from phantom ROIs

    % %3.8 Specify indent (to avoid fitting dead space at edge of phantom) and
    % %filtersize
    % imData.fittingIndent=30;

    %3.9 Get settings
    settings.echotimes = 1000*imData.TE;
    settings.fieldStrength = imData.FieldStrength;

    %Specify echotimes shape to ensure consistently correct
    settings.echotimes = reshape(settings.echotimes,[6 1]);

    %Use existing sigma estimates
    settings.sigmaEst = sigmaEstimates(n);

    %3.11 Train networks with settings (echotimes, fieldstrength) for chosen
    %dataset
    settings.noisyTraining = 1;
    nets = trainNetworks(settings);

    %3.12 Check performance on simulation data
    % testOnSimulatedData(nets,settings)

    %3.13 Get the signals for the chosen slice
    for echoN=1:numel(settings.echotimes)
        complexSlice(:,:,echoN)=imData.images(:,:,sl,1,echoN);
        magnitudeSlice(:,:,echoN)=abs(complexSlice(:,:,echoN));
    end

    % %3.17 Save variables (mapsWithSigma,filteredSigma,maps)
    saveFileName = strcat('MAPS_', dataFileName);
    save(fullfile(saveFolder,saveFileName), 'nets');
end


%% 3B. IMPLEMENT networks with echotimes corresponding to PHANTOM datasets
%3.1 Define folders for import of multiecho data and ROIs
imageFolder='/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO/Data/Hernando data';
roiFolder='/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO/Data/Hernando ROIs';

%Network folder for loading from
loadFolder='/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO DL/DixonDL/Results/HernandoPhantomNetworks';
% loadFolder='/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO DL/DixonDL/Results/HernandoPhantomNetworks';

%Folder for saving results to
saveFolder='/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO DL/DixonDL/Results/HernandoPhantomResults';

% Load sigma estimates
load('sigmaEstimates.mat')

%3.2 Get image and ROI folder info
imageFolderInfo=dir(imageFolder);
roiFolderInfo=dir(roiFolder);
loadFolderInfo=dir(loadFolder);
saveFolderInfo=dir(saveFolder);

%3.3. Specify reference values in phantoms
% Create grid of reference values based on phantom structure
ReferenceValues.FF = ([0; 0.026; 0.053; 0.079; 0.105; 0.157; 0.209; 0.312; 0.413; 0.514; 1]);
% ReferenceValues.R2 = repelem([0; 0; 0; 0],1,5);

%3.4 Specify slice (use 1 if only 1 slice)
sl=2;

tic
%3.5 Loop over datasets to fit each dataset
for n=1:28

    %Clear unwanted variables
    clear complexSlice magnitudeSlice nets predictions

    %3.5 Define filenames for multiecho data and ROIs
    dataFileName = imageFolderInfo(n+2).name
    roiFileName = roiFolderInfo(n+2).name;

    %3.6 Import data
    struct=load(fullfile(imageFolder,dataFileName));
    imData=struct.imDataAll;
    fwmc_ff=struct.fwmc_ff;
    fwmc_r2star=struct.fwmc_r2star;

    %3.7 Import ROIs
    mask = niftiread(fullfile(roiFolder,roiFileName));
    roi.mask=mask(:,:,sl);
    roi.slice=sl;

    % phantomImage = niftiread(fullfile(folder,dataFileName));

    % % Show check image for alignment
    % figure
    % subplot(2,3,1)
    % imshow(phantomROIs(:,:,1),[])
    % title('ROIs')
    %
    % subplot(2,3,2)
    % imshow(abs(imData.images(:,:,sl,1,1)),[])
    % title('Magnitude image of phantom (first echo)')
    %
    % subplot(2,3,3)
    % imshow(fwmc_ff(:,:,2),[0 100])
    % title('PDFF map of phantom')
    %
    % subplot(2,3,4)
    % imshow(phantomROIs(:,:,sl),[0 12])
    % title('Phantom ROIs')

    % Perform fitting of multiecho data, using sigma estimate from phantom ROIs

    % %3.8 Specify indent (to avoid fitting dead space at edge of phantom) and
    % %filtersize
    % imData.fittingIndent=30;

    %3.9 Get settings
    imData.echotimes = 1000*imData.TE;
    imData.echotimes = reshape(imData.echotimes,[6 1]); %Specify echotimes shape to ensure consistently correct
    imData.fieldStrength = imData.FieldStrength;
    settings.echotimes=imData.echotimes;
    settings.fieldStrength = imData.fieldStrength;


    %3.10 Get estimate of sigma

    % %Get sigma estimate
    % [mapsWithSigma,sigmaFromRoi,sigmaFromFit] = FitImageSigma(imData,roi); %NB performs sigma where  mask elements == 1 (i.e. PDFF 0)
    %
    % %Add to vector of sigma estimates for all phantoms
    % sigmaEstimates(n) = sigmaFromFit;

    %Use existing sigma estimates
    settings.sigmaEst = sigmaEstimates(n);

    %3.11 Load networks corresponding to chosen dataset
    loadFileName = strcat('MAPS_', dataFileName);
    load(fullfile(loadFolder,loadFileName));

    %Specify tuning parameters for inference (choice of network output)
    settings.tuning.useImplausibleValues = 1;
    settings.tuning.clipOutputs = 1;

    %3.12 Check performance on simulation data
    % testOnSimulatedData(nets,settings)

    %3.13 Get the signals for the chosen slice
    for echoN=1:numel(imData.echotimes)
        complexSlice(:,:,echoN)=imData.images(:,:,sl,1,echoN);
        magnitudeSlice(:,:,echoN)=abs(complexSlice(:,:,echoN));
    end

    %3.14 Get the predictions for the chosen slice
    predictions = testOnImage(magnitudeSlice,nets,settings);

    %3.15 Generate maps
    ffLow = predictions.prediction1(:,:,1);
    ffHigh = predictions.prediction2(:,:,1);
    ffMap = predictions.prediction3(:,:,1);

    r2Low = predictions.prediction1(:,:,2);
    r2High = predictions.prediction2(:,:,2);
    r2Map = predictions.prediction3(:,:,2);

    %3.16 Display maps

    figure
    subplot(2,3,1)
    imshow(ffLow,[0 1])
    title('PDFF - Water network output')
    colormap('parula')
    colorbar

    subplot(2,3,2)
    imshow(ffHigh,[0 1])
    title('PDFF - Fat network output')
    colormap('parula')
    colorbar

    subplot(2,3,3)
    imshow(ffMap,[0 1])
    title('PDFF - Likelihood-chosen output')
    colormap('parula')
    colorbar

    subplot(2,3,4)
    imshow(r2Low,[0 0.2])
    title('R2* - Water network output')
    colorbar

    subplot(2,3,5)
    imshow(r2High,[0 0.2])
    title('R2* - Fat network output')
    colorbar

    subplot(2,3,6)
    imshow(r2Map,[0 0.2])
    title('R2* - Likelihood-chosen output')
    colorbar

    %3.17 Save variables (mapsWithSigma,filteredSigma,maps)

    saveFileName = strcat('MAPS_', dataFileName);
    save(fullfile(saveFolder,saveFileName), 'predictions','nets');

end
toc

%% 4. Hernando phantom data ROI analysis (split off to enable rapid modification of figures)

%Folder for analysis
saveFolder='/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO DL/DixonDL/Results/HernandoPhantomResults';
saveFolderInfo=dir(saveFolder);

% for n=1:(numel(saveFolderInfo)-2)
for n=1:28

    %3.1 Define filenames for multiecho data and ROIs
    dataFileName = imageFolderInfo(n+2).name
    roiFileName = roiFolderInfo(n+2).name;

    %3.2 Import data
    struct=load(fullfile(imageFolder,dataFileName));
    imData=struct.imDataAll;
    fwmc_ff=struct.fwmc_ff;
    fwmc_r2star=struct.fwmc_r2star;

    %3.3 Import ROIs
    phantomROIs = niftiread(fullfile(roiFolder,roiFileName));

    %3.4 Load MAGORINO maps (created earlier)
    load(fullfile(saveFolder,strcat('MAPS_', dataFileName)));

    ffMap = predictions.prediction3(:,:,1);
    r2Map = predictions.prediction3(:,:,2);


    %3.5 Generate example figure for site 1, 3T, protocol 2
    if n==4

        % Generate mask
        im= abs(imData.images(:,:,2,1,2));
        thresh = prctile(reshape(im,1,[]),94);
        mask = im>thresh;
        figure, imshow(mask,[])

        %Use MAGORINO as a comparator
        load(fullfile('/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO/Data/Hernando results Revision2','MAPS_site1_begin_3T_protocol2.mat'));

        figure

        subplot(1,3,1)
        s1 = imshow(maps.FFrician(80:160,50:200),[0 1])
        set(s1,'Alphadata',mask(80:160,50:200))
        a=colorbar
        colormap('parula')
        ylabel(a,'MAGORINO PDFF','FontSize',12)

        subplot(1,3,2)
        s2 = imshow(ffMap(80:160,50:200),[0 1])
        set(s2,'Alphadata',mask(80:160,50:200))
        a=colorbar
        colormap('parula')
        ylabel(a,'RAIDER PDFF','FontSize',12)

        subplot(1,3,3)
        s3 = imshow(ffMap(80:160,50:200) - maps.FFrician(80:160,50:200),[-0.1 0.1])
        set(s3,'Alphadata',mask(80:160,50:200))
        a=colorbar
        colormap('parula')
        ylabel(a,'RAIDER PDFF - MAGORINO PDFF','FontSize',12)


    else ;
    end

    %3.5 Phantom ROI analysis
    [ff{n},regressionModels{n}] = PhantomRoiAnalysisDL(predictions,phantomROIs(:,:,sl),ReferenceValues,fwmc_ff,fwmc_r2star);
end


%3.7 Tabulate coefficients and get figures

site = [1 1 1 1 7 7 7 7 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6];
protocol = [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4];

[tableVar] = tabulateCoeffsDL(ff,regressionModels,site, protocol,ReferenceValues,saveFolderInfo);

legend('Unity', 'Site 1', 'Site 2', 'Site 3', 'Site 4', 'Site 5', 'Site 6', 'Site 7')



