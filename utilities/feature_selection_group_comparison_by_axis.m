clear;
close all;

load('feature_selection_group1.mat');
mean_rmse = [rmse_mean_feat_obj;rmse_mean_feat_objshp;rmse_mean_feat_all];
std_rmse = [rmse_std_feat_obj;rmse_std_feat_objshp;rmse_std_feat_all];

load('feature_selection_group2.mat');

mean_group_2_obj = mean([rmse_mean_feat_a_obj;rmse_mean_feat_b_obj]);
mean_group_2_objshp = mean([rmse_mean_feat_a_objshp;rmse_mean_feat_b_objshp]);
mean_group_2_all = mean([rmse_mean_feat_a_all;rmse_mean_feat_b_all]);

std_group_2_obj = mean([rmse_std_feat_a_obj;rmse_std_feat_b_obj]);
std_group_2_objshp = mean([rmse_std_feat_a_objshp;rmse_std_feat_b_objshp]);
std_group_2_all = mean([rmse_std_feat_a_all;rmse_std_feat_b_all]);

mean_rmse = [mean_rmse;mean_group_2_obj;mean_group_2_objshp;mean_group_2_all];
std_rmse = [std_rmse;std_group_2_obj;std_group_2_objshp;std_group_2_all];

load('feature_selection_group3.mat');
group_3_data_obj=[];
group_3_data_objshp=[];
group_3_data_all=[];
group_3_std_obj=[];
group_3_std_objshp=[];
group_3_std_all=[];
for i=1:3
    group_3_data_obj = [group_3_data_obj;rmse_mean_feat_obj{i}];
    group_3_data_objshp = [group_3_data_objshp;rmse_mean_feat_objshp{i}];
    group_3_data_all = [group_3_data_all;rmse_mean_feat_all{i}];
    group_3_std_obj = [group_3_std_obj;rmse_std_feat_obj{i}];
    group_3_std_objshp = [group_3_std_obj;rmse_std_feat_objshp{i}];
    group_3_std_all = [group_3_std_all;rmse_std_feat_all{i}];
end

mean_rmse = [mean_rmse;mean(group_3_data_obj);mean(group_3_data_objshp);mean(group_3_data_all)];
std_rmse = [std_rmse;mean(group_3_std_obj);mean(group_3_std_objshp);mean(group_3_std_all)];

load('feature_selection_group4.mat');
group_4_data_obj=[];
group_4_data_objshp=[];
group_4_data_all=[];
group_4_std_obj=[];
group_4_std_objshp=[];
group_4_std_all=[];
for i=1:4
    group_4_data_obj = [group_4_data_obj;rmse_mean_feat_obj{i}];
    group_4_data_objshp = [group_4_data_objshp;rmse_mean_feat_objshp{i}];
    group_4_data_all = [group_4_data_all;rmse_mean_feat_all{i}];
    group_4_std_obj = [group_4_std_obj;rmse_std_feat_obj{i}];
    group_4_std_objshp = [group_4_std_obj;rmse_std_feat_objshp{i}];
    group_4_std_all = [group_4_std_all;rmse_std_feat_all{i}];
end

mean_rmse = [mean_rmse;mean(group_4_data_obj);mean(group_4_data_objshp);mean(group_4_data_all)];
std_rmse = [std_rmse;mean(group_4_std_obj);mean(group_4_std_objshp);mean(group_4_std_all)];

% load('gp_cv_nogroup_result_bettersample.mat');
% gp_mean_rmse = mean(rmse);
% gp_std_rmse = std(rmse);
% 
% % load('gp_cv_nogroup_result_half.mat');
% % gp_mean_rmse = nogroup_rmse_mean;
% % gp_std_rmse = nogroup_rmse_std;
% 
% %load('gp_cv_group.mat');
% load('gp_cv_group_gmm4.mat');
% free_rmse_mean_full_feat = free_rmse_mean;
% free_rmse_std_full_feat = free_rmse_std;
% cont_rmse_mean_full_feat = cont_rmse_mean;
% cont_rmse_std_full_feat = cont_rmse_std;
% 
% load('gp_cv_group_gmm4_lessfeat.mat');
% free_rmse_mean_less_feat = free_rmse_mean;
% free_rmse_std_less_feat = free_rmse_std;
% cont_rmse_mean_less_feat = cont_rmse_mean;
% cont_rmse_std_less_feat = cont_rmse_std;
% 
% mean_rmse_group_fs = mean([free_rmse_mean_less_feat;cont_rmse_mean_full_feat]);
% std_rmse_group_fs = mean([free_rmse_std_less_feat;cont_rmse_std_full_feat]);
% 
% mean_rmse_group = mean([free_rmse_mean_full_feat;cont_rmse_mean_full_feat]);
% std_rmse_group = mean([free_rmse_std_full_feat;cont_rmse_std_full_feat]);

% mean_rmse = [free_rmse_mean;cont_rmse_mean;gp_mean_rmse];
% std_rmse = [free_rmse_std;cont_rmse_std;gp_std_rmse];

% mean_rmse = [gp_mean_rmse;mean_rmse_group;mean_rmse_group_fs];
% std_rmse = [gp_std_rmse;std_rmse_group;std_rmse_group_fs];

y_p_mean = [mean_rmse(1:3,1) mean_rmse(4:6,1) mean_rmse(7:9,1) mean_rmse(10:12,1)];
y_p_std = [std_rmse(1:3,1) std_rmse(4:6,1) std_rmse(7:9,1) std_rmse(10:12,1)];

z_p_mean = [mean_rmse(1:3,2) mean_rmse(4:6,2) mean_rmse(7:9,2) mean_rmse(10:12,2)];
z_p_std = [std_rmse(1:3,2) std_rmse(4:6,2) std_rmse(7:9,2) std_rmse(10:12,2)];

theta_p_mean = [mean_rmse(1:3,3) mean_rmse(4:6,3) mean_rmse(7:9,3) mean_rmse(10:12,3)];
theta_p_std = [std_rmse(1:3,3) std_rmse(4:6,3) std_rmse(7:9,3) std_rmse(10:12,3)];

% theta_scale = 10;
% mean_rmse(:,3) = mean_rmse(:,3)/theta_scale;
% std_rmse(:,3) = std_rmse(:,3)/theta_scale; 

% load('gp_cv_group_gmm5.mat');
% mean_rmse = [free_rmse_mean;cont_rmse_mean;gp_mean_rmse];
% std_rmse = [free_rmse_std;cont_rmse_std;gp_std_rmse];

FigHandle = figure('Position', [100, 100, 500, 500]);
hold on
hb = bar(1:4,y_p_mean');
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(xData,y_p_mean(ib,:),y_p_std(ib,:),'k.')
end

hb(1).FaceColor = [1 0 0];
hb(2).FaceColor = [0.0 1 0];
hb(3).FaceColor = [0.0 0 1];
% hb(4).FaceColor = [1 0 0];
% hb(5).FaceColor = [0.0 1 0];
% hb(6).FaceColor = [0.0 0 1];
% hb(7).FaceColor = [1 0 0];
% hb(8).FaceColor = [0.0 1 0];
% hb(9).FaceColor = [0.0 0 1];
% hb(10).FaceColor = [1 0 0];
% hb(11).FaceColor = [0.0 1 0];
% hb(12).FaceColor = [0.0 0 1];

Labels = {'1', '2','3','4'};
set(gca, 'XTick', 1:4, 'XTickLabel', Labels);

legend('Object features','Object+Contact features','Object+Contact+Context features','Location','northeast');
title('RMSE of Prediction Models for \Deltay_p');
set(gca,'ylim',[0 0.05]);
ylabel('meter');