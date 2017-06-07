clear;
close all;

load('gp_cv_nogroup_result_bettersample.mat');
gp_mean_rmse = mean(rmse);
gp_std_rmse = std(rmse);

% load('gp_cv_nogroup_result_half.mat');
% gp_mean_rmse = nogroup_rmse_mean;
% gp_std_rmse = nogroup_rmse_std;

%load('gp_cv_group.mat');
load('gp_cv_group_gmm4.mat');
free_rmse_mean_full_feat = free_rmse_mean;
free_rmse_std_full_feat = free_rmse_std;
cont_rmse_mean_full_feat = cont_rmse_mean;
cont_rmse_std_full_feat = cont_rmse_std;

load('gp_cv_group_gmm4_lessfeat.mat');
free_rmse_mean_less_feat = free_rmse_mean;
free_rmse_std_less_feat = free_rmse_std;
cont_rmse_mean_less_feat = cont_rmse_mean;
cont_rmse_std_less_feat = cont_rmse_std;

mean_rmse_group_fs = mean([free_rmse_mean_less_feat;cont_rmse_mean_full_feat]);
std_rmse_group_fs = mean([free_rmse_std_less_feat;cont_rmse_std_full_feat]);

mean_rmse_group = mean([free_rmse_mean_full_feat;cont_rmse_mean_full_feat]);
std_rmse_group = mean([free_rmse_std_full_feat;cont_rmse_std_full_feat]);

% mean_rmse = [free_rmse_mean;cont_rmse_mean;gp_mean_rmse];
% std_rmse = [free_rmse_std;cont_rmse_std;gp_std_rmse];

mean_rmse = [gp_mean_rmse;mean_rmse_group;mean_rmse_group_fs];
std_rmse = [gp_std_rmse;std_rmse_group;std_rmse_group_fs];


theta_scale = 10;
mean_rmse(:,3) = mean_rmse(:,3)/theta_scale;
std_rmse(:,3) = std_rmse(:,3)/theta_scale; 

% load('gp_cv_group_gmm5.mat');
% mean_rmse = [free_rmse_mean;cont_rmse_mean;gp_mean_rmse];
% std_rmse = [free_rmse_std;cont_rmse_std;gp_std_rmse];

figure
hold on
hb = bar(1:3,mean_rmse');
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(xData,mean_rmse(ib,:),std_rmse(ib,:),'k.')
end

hb(1).FaceColor = [1 0 0];
hb(2).FaceColor = [0.0 1 0];
hb(3).FaceColor = [0.0 0 1];

Labels = {'y_p', 'z_p','theta'};
set(gca, 'XTick', 1:3, 'XTickLabel', Labels);

legend('Without hierarchical','Hierarchical modeling','Hierarchical with selected features','Location','northwest');
title('RMSE of Forward Affordance Prediction Models');

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