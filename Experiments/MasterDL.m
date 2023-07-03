% Master script
% Runs the various simulations and generates figures
% Author: Tim Bray (t.bray@ucl.ac.uk)

%% 1. Simulations

%% First single network only 
% Note: trainSignleNetworks actually returns two networks but one is a duplicate of the first - simplifies code by avoiding need for separate single-network script) 

%1.1 Specify settings
settings.echotimes = [1.1:1.1:13.2]';
settings.fieldStrength = 3;

% settings.SNR = 60;
% settings.noiseSD = 1/settings.SNR;

%For noise-free simulations, set settings.noiseSD = 0; 
settings.noiseSD = 0;

%1.2 Train single network with uniform training distribution
net = trainSingleNetwork(settings);

%1.3 Test on simulation data
[dlMaps,dlErrormaps] = testOnSimulatedData(net,settings)

%1.4 Show comparison vs conventional fitting

%Load conventional fitting data for comparison
load('/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO DL/DixonDL/Results/SimulationResults/MAGORINO_Simulation_Results_SNR60_R100_fineGrain.mat')

%Visualise comparison against conventional fitting
createFigDLvsConventionalFitting(FFmaps, dlMaps, errormaps, dlErrormaps)


%% Dual networks 

%1.4 Train networks with settings (echotimes, fieldstrength) for chosen
%dataset
nets = trainNetworks(settings);

%1.5 Test on simulation data
[dlMaps,dlErrormaps] = testOnSimulatedData(nets,settings)

%1.6 Show comparison vs conventional fitting
%Load conventional fitting data for comparison
load('/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO DL/DixonDL/Results/MAGORINO_Simulation_Results_SNR60_R100_fineGrain.mat')

% Visualise comparison against conventional fitting
createFigDLvsConventionalFitting(FFmaps, dlMaps, errormaps, dlErrormaps)

%% 2. Train and implement network with echotimes corresponding to specific SUBJECT dataset
clear all

%2.1 Load chosen dataset
% load('/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO/Data/Subjects/FW101_bipolar.mat')
% load('/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO/Data/Subjects/FW111_monopolar.mat')
load('/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO/Data/Subjects/LegsSwapData.mat')

%Convert to imData if not already
try
imDataAll = imDataParams;
catch
end

%2.2 Get settings
settings.echotimes = 1000*imDataAll.TE';
settings.fieldStrength = imDataAll.FieldStrength;

%2.3 Provide guess for SNR and noiseSD
SNRguess = 50;
settings.SNR = SNRguess;
settings.noiseSD = 1/SNRguess;

%2.4 Train networks with settings (echotimes, fieldstrength) for chosen
%dataset
nets = trainNetworks(settings);

%2.5 Check working ok on simulation data
testOnSimulatedData(nets,settings)

%2.6 Use networks to get predicts for the chosen dataset

%Specify the slice 
slice = 20;

%Get the signals for the chosen slice
for echoN=1:numel(settings.echotimes)  
complexSlice(:,:,echoN)=imDataAll.images(:,:,slice,1,echoN);
magnitudeSlice(:,:,echoN)=abs(complexSlice(:,:,echoN));

%Include normalisation by the first echo
% magnitudeSlice(:,:,:) = magnitudeSlice(:,:,:) / magnitudeSlice(:,:,1);

end

%Call the function
predictions = testOnImage(magnitudeSlice,nets,settings);

%2.8 Get maps
ffLow = predictions.prediction1(:,:,1);
ffHigh = predictions.prediction2(:,:,1);
ffMap = predictions.prediction3(:,:,1);

r2Low = predictions.rawPrediction1(:,:,2);
r2High = predictions.rawPrediction2(:,:,2);
r2Map = predictions.prediction3(:,:,2);

%2.7 Display maps

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
imshow(1000*r2Low,[0 200])
title('R2* - Water network output')
colorbar

subplot(2,3,5)
imshow(1000*r2High,[0 200])
title('R2* - Fat network output')
colorbar

subplot(2,3,6)
imshow(1000*r2Map,[0 200])
title('R2* - Likelihood-chosen output')
colorbar

%2.8 Display comparison with conventional fitting

load('/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO/Experiments/experimentResults/LegsSwapDataResults.mat')

yDisp = (65:175);

figure
subplot(3,2,1)
imshow(maps.FFstandard(yDisp,:),[0 1])
colormap('parula')
title('PDFF - Conventional fitting')
h=colorbar
h.Label.String = "PDFF";
h.FontSize=12; 

subplot(3,2,3)
imshow(ffMap(yDisp,:),[0 1])
colormap('parula')
title('PDFF - DL')
h=colorbar
h.Label.String = "PDFF";
h.FontSize=12; 

subplot(3,2,5)
imshow(ffMap(yDisp,:) - maps.FFstandard(yDisp,:),[-1 1])
title('PDFF difference (DL - conventional)')
h=colorbar
h.Label.String = "PDFF difference";
h.FontSize=12; 

subplot(3,2,2)
imshow(1000*maps.R2standard(yDisp,:),[0 500])
title('R2* - Conventional fitting')
h=colorbar
h.Label.String = "R2* (s^-^1)";
h.FontSize=12; 

subplot(3,2,4)
imshow(1000*r2Map(yDisp,:),[0 500])
title('R2* - DL')
h=colorbar
h.Label.String = "R2* (s^-^1)";
h.FontSize=12; 

subplot(3,2,6)
imshow(1000*r2Map(yDisp,:) - 1000*maps.R2standard(yDisp,:),[-500 500])
title('R2* difference (DL - conventional)')
h=colorbar
h.Label.String = "R2* difference (s^-^1)";
h.FontSize=12; 

%% 3. Train and implement network with echotimes corresponding to PHANTOM dataset

%3.1 Define folders for import of multiecho data and ROIs
imageFolder='/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO/Data/Hernando data';
roiFolder='/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO/Data/Hernando ROIs';

%Folder for saving
saveFolder='/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO DL/Results';

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
for n=23:28

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

%3.10 Provide guess for SNR and noiseSD
if round(settings.fieldStrength) == 3 %If field strength is 3
    SNRguess = 50;

elseif round(2*settings.fieldStrength) == 3 %If field strength is 1.5
    SNRguess = 30;
else ;
end

settings.noiseSD = 1/SNRguess;

%3.11 Train networks with settings (echotimes, fieldstrength) for chosen
%dataset
nets = trainNetworks(settings);

%3.12 Check performance on simulation data
% testOnSimulatedData(nets,settings)

%3.13 Get the signals for the chosen slice
for echoN=1:numel(settings.echotimes)  
complexSlice(:,:,echoN)=imData.images(:,:,sl,1,echoN);
magnitudeSlice(:,:,echoN)=abs(complexSlice(:,:,echoN));
end

%3.14 Get the predictions for the chosen slice
predictions = testOnImage(magnitudeSlice,nets,settings);

%3.15 Generate maps
ffLow = predictions.prediction1(:,:,1);
ffHigh = predictions.prediction2(:,:,1);
ffMap = predictions.prediction3(:,:,1);

r2Low = predictions.rawPrediction1(:,:,2);
r2High = predictions.rawPrediction2(:,:,2);
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

% %3.17 Save variables (mapsWithSigma,filteredSigma,maps)

saveFileName = strcat('MAPS_', dataFileName);
save(fullfile(saveFolder,saveFileName), 'predictions');

end


%% 4. Hernando phantom data ROI analysis (split off to enable rapid modification of figures)

%Folder for analysis
saveFolder='/Users/tjb57/Dropbox/MATLAB/Fat-water MAGORINO DL/DixonDL/Results/Hernando phantom';
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

figure

subplot(1,2,1)
imshow(ffMap,[0 1])
a=colorbar
colormap('parula')
ylabel(a,'PDFF','FontSize',12)

subplot(1,2,2)
imshow(1000*r2Map,[0 1000])
a=colorbar
colormap('parula')
ylabel(a,'R2*','FontSize',12)

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



