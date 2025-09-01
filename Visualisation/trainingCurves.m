function trainingCurves(net,titleText)

%Get RMSE
trainingRMSE = net.info1.TrainingRMSE;
validationRMSE = net.info1.ValidationRMSE;

% Define scale factor for max
scaling = [0.2 0.3];

%% 
figure('Name',titleText)
subplot(1,2,1)
plot(smooth(trainingRMSE,5000),'LineWidth',2,'MarkerSize',6,'Marker','.')
ylim(scaling*max(trainingRMSE))
ylim([0 .5])
xlim([0 numel(trainingRMSE)])
xticks([0 0.2 0.4 0.6 0.8 1]*numel(trainingRMSE))
xticklabels([0 10 20 30 40 50])
ylabel('Training RMSE','FontSize',14)
xlabel('Epochs','FontSize',14)

subplot(1,2,2)
plot(smooth(validationRMSE(validationRMSE>0),100),'k', 'LineWidth',2,'MarkerSize',6,'Marker','.')
ylim(scaling*max(trainingRMSE))
ylim([0 .5])
ylabel('Validation RMSE','FontSize',16)
xticks([0 0.2 0.4 0.6 0.8 1]*numel(validationRMSE(validationRMSE>0)))
xticklabels([0 10 20 30 40 50])
xlim([0 numel(validationRMSE(validationRMSE>0))])
xlabel('Epochs','FontSize',14)
end


