function trainingCurves(net)

%Get RMSE
trainingRMSE = net.info1.TrainingRMSE;
validationRMSE = net.info1.ValidationRMSE;

% Define scale factor for max
scaling = [0.1 0.45];


%% 
figure
subplot(1,2,1)
plot(smooth(trainingRMSE,100))
ylim(scaling*max(trainingRMSE))

subplot(1,2,2)
plot(validationRMSE(validationRMSE>0),'k')
ylim(scaling*max(trainingRMSE))

end


