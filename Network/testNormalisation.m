
function testNormalisation(settings)
%Evaluates the quality of normalisation for a given dataset

%Select values for test dataset
r2max=0.5;

%Extract echotimes
echotimes = settings.echotimes;
noiseSD = settings.noiseSD;
tesla = settings.fieldStrength;

S0=1;
FFvals = (0:0.01:1)';
R2vals = (0:0.05*r2max:r2max);

%Select number of repetitions

%Call helper function (sVecFixedSpacing) to generate vectors of values with
%fixed spacing
[paramVec, sVecNoiseFree] = sVecFixedSpacing(S0,FFvals,R2vals,settings);

%Create yTest parameter vector
yTest = paramVec;

%Loop over noise instantiations 

% reps = 100;
% 
% for r = 1:reps

%Create noise
realnoise=(noiseSD)*randn(size(sVecNoiseFree,1),numel(echotimes));
imagnoise=1i*(noiseSD)*randn(size(sVecNoiseFree,1),numel(echotimes));

%Add noise to signal 
sVec = sVecNoiseFree + realnoise + imagnoise;

%Magnitude
sVec = abs(sVec);

% Normalisation
[sVecNorm,s0estimates] = normaliseSignals(sVec,settings);

% Reshape prediction data for plotting
s0estimatesGrid = reshape(s0estimates,[numel(FFvals) numel(R2vals)]);

%Visualise s0 estimates
figure
image(s0estimatesGrid,'CDataMapping','scaled')
ax=gca;
ax.CLim=[0 2];
FigLabels;
colorbar
xticks([1:2:21]);
xticklabels({'0','.50', '100', '150', '200', '250', '300', '350', '400', '450','500'});
xlabel('R2* (s^-^1)','FontSize',12)
yticks([0 10 20 30 40 50 60 70 80 90 100]);
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});
ylabel('Fat fraction','FontSize',12)
title('s0 estimates over parameter space')
h=colorbar
h.Label.String = "s0estimates";
h.Label.FontSize = 12;
view(0,90)

%% Plot signals for chosen FF, R2*
ffDisp = 0.1; 
r2Disp = 0.2;

ind = find(paramVec(:,1)==ffDisp & paramVec(:,2)==r2Disp);

figure
plot(echotimes, sVec(ind,:))
xlabel('Echotime (ms)')
ylabel('Signal')

s0estimates(ind)

end

