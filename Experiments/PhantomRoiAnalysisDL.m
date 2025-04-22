%% Extract values from PDFF maps and phantom ROIs and generate agreement plots against reference standard
% Author:
% Tim Bray, t.bray@ucl.ac.uk
% January 2022

function [ff,regressionModels] = PhantomRoiAnalysis(maps,phantomROIs,ReferenceValues,ffVendor,r2starVendor)


%% 1. Check if cropping required and if so do this

if numel(maps.prediction3(:,:,1)) < numel(phantomROIs)
    sz = size(maps.prediction3) ;
    phantomROIs = phantomROIs(1:sz(1),1:sz(2),:);
else ;
end

%% 2. Get measurements for each ROI

%2.1 For FF measurements
[ff.DL.mean, ff.DL.median, ff.DL.sd, ff.DL.n] = RoiStats(maps.prediction3(:,:,1),phantomROIs);
[ff.vendor.mean, ff.vendor.median, ff.vendor.sd, ff.vendor.n] = RoiStats(ffVendor,phantomROIs);

%% 3. Reshape for visualisation

%3.1 For FF measurements
% ff.standard.meangrid = reshape(ff.standard.mean,[4,5]);


%% 4. Perform linear regression

% 5.1 For FF
regressionModels = fitlm(ReferenceValues.FF,ff.DL.median); %RM denotes regression model


%% Choose whether to display data
disp = 0;

if disp ==1


    %% 4. Display phantom images

    figure
    subplot(1,2,1)
    imshow(maps.prediction3(:,:,1),[0 1])
    title('Fat fraction, DL')
    colormap('parula')
    colorbar

    subplot(1,2,2)
    imshow(0.01*ffVendor(:,:,1),[0 1])
    title('Fat fraction, Vendor')
    colormap('parula')
    colorbar

    figure
    subplot(1,2,1)
    imshow(1000*maps.prediction3(:,:,2),[0 100])
    title('R2*, DL')
    colormap('parula')
    colorbar

    subplot(1,2,2)
    imshow(r2starVendor(:,:,1),[0 100])
    title('R2*, Vendor')
    colormap('parula')
    colorbar


    %% 6. Display agreement plots


    % 6.1 For FF
    figure('Name', 'FF')
    subplot(2,4,1)
    plot(ReferenceValues.FF,ff.DL.median,'x-','LineWidth',1)
    hold on
    plot(ReferenceValues.FF,0.01*ff.vendor.median,'x--','LineWidth',1)
    xlabel('Reference FF')
    ylabel('PDFF')
    legend('DL','Hernando complex fitting')
    xlim([0 1])
    ylim([0 1])

    %% 7. Display values as grids

    % % 5.1 For FF
    % figure('Name', 'FF')
    %
    % s1=subplot(2,4,1)
    % image(ReferenceValues.FF,'CDataMapping','scaled')
    % ax=gca;
    % ax.CLim=[0 1];
    % xticks([1 2 3 4 5]);
    % xticklabels({'0','20', '40', '50', '60'});
    % xlabel('Fat fraction (%)','FontSize',12)
    % yticks([1 2 3 4]);
    % yticklabels({'150','100','50','0'});
    % ylabel('Bone mineral density','FontSize',12)
    % title('Reference FF')
    % colorbar
    %
    % s1=subplot(2,4,2)
    % image(ff.standard.meangrid,'CDataMapping','scaled')
    % ax=gca;
    % ax.CLim=[0 1];
    % xticks([1 2 3 4 5]);
    % xticklabels({'0','20', '40', '50', '60'});
    % xlabel('Fat fraction (%)','FontSize',12)
    % yticks([1 2 3 4]);
    % yticklabels({'150','100','50','0'});
    % ylabel('Bone mineral density','FontSize',12)
    % title('Gaussian FF')
    % colorbar
    %
    % s1=subplot(2,4,3)
    % image(ff.rician.meangrid,'CDataMapping','scaled')
    % ax=gca;
    % ax.CLim=[0 1];
    % xticks([1 2 3 4 5]);
    % xticklabels({'0','20', '40', '50', '60'});
    % xlabel('Fat fraction (%)','FontSize',12)
    % yticks([1 2 3 4]);
    % yticklabels({'150','100','50','0'});
    % ylabel('Bone mineral density','FontSize',12)
    % title('Rician FF')
    % colorbar
    %
    % s1=subplot(2,4,4)
    % image(ff.complex.meangrid,'CDataMapping','scaled')
    % ax=gca;
    % ax.CLim=[0 1];
    % xticks([1 2 3 4 5 6 7 8 9 10 11]);
    % xticklabels({'0','.1', '.2', '.3', '.4', '.5', '.6', '.7', '.8', '.9','1.0'});
    % xlabel('R2* (ms^-^1)','FontSize',12)
    % yticks([1 6 11 16 21 26 31 36 41 46 51]);
    % yticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    % ylabel('Fat fraction (%)','FontSize',12)
    % title('Complex FF')
    % colorbar
    %
    % s1=subplot(2,4,6)
    % image(ff.standard.meangrid-ReferenceValues.FF,'CDataMapping','scaled')
    % ax=gca;
    % ax.CLim=[-1 1];
    % xticks([1 2 3 4 5]);
    % xticklabels({'0','20', '40', '50', '60'});
    % xlabel('Fat fraction (%)','FontSize',12)
    % yticks([1 2 3 4]);
    % yticklabels({'150','100','50','0'});
    % ylabel('Bone mineral density','FontSize',12)
    % title('Gaussian FF error')
    % colorbar
    %
    % s1=subplot(2,4,7)
    % image(ff.rician.meangrid-ReferenceValues.FF,'CDataMapping','scaled')
    % ax=gca;
    % ax.CLim=[-1 1];
    % xticks([1 2 3 4 5 6 7 8 9 10 11]);
    % xticklabels({'0','.1', '.2', '.3', '.4', '.5', '.6', '.7', '.8', '.9','1.0'});
    % xlabel('R2* (ms^-^1)','FontSize',12)
    % yticks([1 6 11 16 21 26 31 36 41 46 51]);
    % yticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    % ylabel('Fat fraction (%)','FontSize',12)
    % title('Rician FF error')
    % colorbar
    %
    % s1=subplot(2,4,8)
    % image(ff.complex.meangrid-ReferenceValues.FF,'CDataMapping','scaled')
    % ax=gca;
    % ax.CLim=[-1 1];
    % xticks([1 2 3 4 5 6 7 8 9 10 11]);
    % xticklabels({'0','.1', '.2', '.3', '.4', '.5', '.6', '.7', '.8', '.9','1.0'});
    % xlabel('R2* (ms^-^1)','FontSize',12)
    % yticks([1 6 11 16 21 26 31 36 41 46 51]);
    % yticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    % ylabel('Fat fraction (%)','FontSize',12)
    % title('Complex FF error')
    % colorbar
    %
    % % 5.1 For R2star
    %
    % figure('Name', 'R2*')
    % s1=subplot(2,3,1)
    % image(r2.standard.meangrid,'CDataMapping','scaled')
    % ax=gca;
    % ax.CLim=[0 0.3];
    % xticks([1 2 3 4 5]);
    % xticklabels({'0','20', '40', '50', '60'});
    % xlabel('Fat fraction (%)','FontSize',12)
    % yticks([1 2 3 4]);
    % yticklabels({'150','100','50','0'});
    % ylabel('Bone mineral density','FontSize',12)
    % title('R2* Gaussian/MAGO')
    % colorbar
    %
    % s1=subplot(2,3,2)
    % image(r2.rician.meangrid,'CDataMapping','scaled')
    % ax=gca;
    % ax.CLim=[0 0.3];
    % xticks([1 2 3 4 5]);
    % xticklabels({'0','20', '40', '50', '60'});
    % xlabel('Fat fraction (%)','FontSize',12)
    % yticks([1 2 3 4]);
    % yticklabels({'150','100','50','0'});
    % ylabel('Bone mineral density','FontSize',12)
    % title('R2* Rician/MAGORINO')
    % colorbar
    %
    % s1=subplot(2,3,3)
    % image(r2.complex.meangrid,'CDataMapping','scaled')
    % ax=gca;
    % ax.CLim=[0 0.3];
    % xticks([1 2 3 4 5]);
    % xticklabels({'0','20', '40', '50', '60'});
    % xlabel('Fat fraction (%)','FontSize',12)
    % yticks([1 2 3 4]);
    % yticklabels({'150','100','50','0'});
    % ylabel('Bone mineral density','FontSize',12)
    % title('R2* complex MAGO-equivalent')
    % colorbar
    %
    % s1=subplot(2,3,5)
    % image(r2.rician.meangrid-r2.standard.meangrid,'CDataMapping','scaled')
    % ax=gca;
    % ax.CLim=[-0.03 0.03];
    % xticks([1 2 3 4 5]);
    % xticklabels({'0','20', '40', '50', '60'});
    % xlabel('Fat fraction (%)','FontSize',12)
    % yticks([1 2 3 4]);
    % yticklabels({'150','100','50','0'});
    % ylabel('Bone mineral density','FontSize',12)
    % title('R2* Rician/MAGORINO - Gaussian/MAGO')
    % colorbar
    %
    % s1=subplot(2,3,6)
    % image(r2.complex.meangrid-r2.standard.meangrid,'CDataMapping','scaled')
    % ax=gca;
    % ax.CLim=[-0.03 0.03];
    % xticks([1 2 3 4 5]);
    % xticklabels({'0','20', '40', '50', '60'});
    % xlabel('Fat fraction (%)','FontSize',12)
    % yticks([1 2 3 4]);
    % yticklabels({'150','100','50','0'});
    % ylabel('Bone mineral density','FontSize',12)
    % title('R2* complex MAGO-equivalent - Gaussian/MAGO')
    % colorbar
    %
    % % s1=subplot(1,4,3)
    % % image(r2.vendor.meangrid,'CDataMapping','scaled')
    % % ax=gca;
    % % ax.CLim=[0 0.3];
    % % xticks([1 2 3 4 5]);
    % % xticklabels({'0','20', '40', '50', '60'});
    % % xlabel('Fat fraction (%)','FontSize',12)
    % % yticks([1 2 3 4]);
    % % yticklabels({'150','100','50','0'});
    % % ylabel('Bone mineral density','FontSize',12)
    % % title('R2* vendor')
    % % colorbar
    %
    % % figure('Name', 'R2star')
    % % s1=subplot(1,3,1)
    % % image(ReferenceValues.R2,'CDataMapping','scaled')
    % % ax=gca;
    % % ax.CLim=[0 150];
    % % xticks([1 2 3 4 5]);
    % % xticklabels({'0','20', '40', '50', '60'});
    % % xlabel('Fat fraction (%)','FontSize',12)
    % % yticks([1 2 3 4]);
    % % yticklabels({'150','100','50','0'});
    % % ylabel('Bone mineral density','FontSize',12)
    % % title('BMD')
    % % colorbar

else ;
end



end


