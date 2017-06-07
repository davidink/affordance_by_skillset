clear;
close all;

load('gp_cv_nogroup_result_bettersample.mat');
gp_mean_rmse = mean(rmse);
gp_std_rmse = std(rmse);

load('gp_cv_nogroup_result_lessfeat.mat');

load('gp_cv_group_gmm4.mat');
%load('gp_cv_group_gmm5.mat');
mean_rmse = [free_rmse_mean;cont_rmse_mean];
std_rmse = [free_rmse_std;cont_rmse_std];

%load('gp_cv_group_gmm5_lessfeat.mat');
load('gp_cv_group_gmm4_lessfeat.mat');

mean_rmse = [mean_rmse;free_rmse_mean;cont_rmse_mean;gp_mean_rmse;nogroup_rmse_mean];
std_rmse = [std_rmse;free_rmse_std;cont_rmse_std;gp_std_rmse;nogroup_rmse_std];

theta_scale = 10;
mean_rmse(:,3) = mean_rmse(:,3)/theta_scale;
std_rmse(:,3) = std_rmse(:,3)/theta_scale; 

figure
hold on
hb = bar(1:3,mean_rmse(:,1:3)');
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(xData,mean_rmse(ib,1:3),std_rmse(ib,1:3),'k.')
end

hb(1).FaceColor = [1 0 0];
hb(2).FaceColor = [0.0 1 0];
hb(3).FaceColor = [0.0 0 1];
hb(4).FaceColor = [0.0 1 1];
hb(5).FaceColor = [1 0 1];
hb(6).FaceColor = [1 1 0];

Labels = {'y', 'z', 'theta'};
set(gca, 'XTick', 1:3, 'XTickLabel', Labels);

legend('Free-space push with all features','Constraint push with all features','Free-space push without contextual features','Constraint push without contextual features','No clustering','No clustering without contextual features','Location','northwest');
title('RMSE of Affordance Prediction Model');

y_lim = get(gca,'ylim');
ylabel('meter');
box off
% Create second Y axes on the right.
a2 = axes('YAxisLocation', 'Right');
% Hide second plot.
set(a2, 'color', 'none');
set(a2, 'XTick', []);
% Set scala for second Y.
set(a2, 'YLim', [y_lim(1)*theta_scale y_lim(2)*theta_scale]);
ylabel('radian');