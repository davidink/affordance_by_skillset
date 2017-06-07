clear;
close all;

load('gp_cv_nogroup_result_bettersample.mat');
gp_mean_rmse = mean(rmse);
gp_std_rmse = std(rmse);

load('blr_cv_result.mat');
blr_mean_rmse = mean(rmse);
blr_std_rmse = std(rmse);
% 
% bar_data_mean = [blr_mean_rmse' gp_mean_rmse'];
% bar_data_std = [blr_std_rmse' gp_std_rmse'];
% b = bar(bar_data_mean);
% errorbar(1:2,bar_data_mean,bar_data_std,'.');
% 
% Labels = {'y', 'z', 'theta'};
% set(gca, 'XTick', 1:3, 'XTickLabel', Labels);

theta_scale = 10;
mean_rmse = [blr_mean_rmse;gp_mean_rmse];
std_rmse = [blr_std_rmse; gp_std_rmse];
mean_rmse(:,3) = mean_rmse(:,3)/theta_scale;
std_rmse(:,3) = std_rmse(:,3)/theta_scale; 
fig_bar = figure;
hold on;
hb = bar(1:3,mean_rmse(:,1:3)');
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData + hb(ib).XOffset;
    errorbar(xData,mean_rmse(ib,1:3),std_rmse(ib,1:3),'k.')
end
hb(1).FaceColor = [1 0 0 ];
hb(2).FaceColor = [0.0 1 0 ];
Labels = {'y_p', 'z_p','theta'};
set(gca, 'XTick', 1:3, 'XTickLabel', Labels);
legend('Bayesian Linear Regression','Gaussian Process','Location','northwest');
title('RMSE of Forward Affordance Prediction Model');

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

% subplot(1,2,2);
% hold on;
% hb2 = bar(1:2',mean_rmse(:,3)',0.05);
% %hb = bar(1:3,mean_rmse');
% % For each set of bars, find the centers of the bars, and write error bars
% pause(0.1); %pause allows the figure to be created
% % for ib = 1:numel(hb2)
% %     %XData property is the tick labels/group centers; XOffset is the offset
% %     %of each distinct group
% %     xData = hb2(ib).XData + hb2(ib).XOffset;
% %     errorbar(xData,mean_rmse(ib,3),std_rmse(ib,3),'k.')
% % end
% hb2(1).FaceColor = [1 0 0 ];
% hb2(2).FaceColor = [0.0 1 0 ];
% Labels = {'theta'};
% set(gca, 'XTick', 1:1, 'XTickLabel', Labels);


