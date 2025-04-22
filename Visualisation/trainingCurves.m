function trainingCurves(net)

%Get RMSE
trainingRMSE = net.info1.TrainingRMSE;
validationRMSE = net.info1.ValidationRMSE;

% Define scale factor for max
scaling = [0.2 0.3];


%% 
figure
subplot(1,2,1)
plot(smooth(trainingRMSE,5000),'LineWidth',2,'MarkerSize',6,'Marker','.')
ylim(scaling*max(trainingRMSE))
ylabel('Training RMSE','FontSize',16)

subplot(1,2,2)
plot(smooth(validationRMSE(validationRMSE>0),100),'k', 'LineWidth',2,'MarkerSize',6,'Marker','.')
ylim(scaling*max(trainingRMSE))
ylabel('Validation RMSE','FontSize',16)
end


