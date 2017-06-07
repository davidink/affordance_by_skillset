clear all
close all

pred_err = [];
load 'rslts/sub_1_feat40feat40_cv1.mat'
pred_err = [pred_err;pred_test_err];
load 'rslts/sub_1_feat40feat40_cv2.mat'
pred_err = [pred_err;pred_test_err];
load 'rslts/sub_1_feat40feat40_cv3.mat'
pred_err = [pred_err;pred_test_err];
pred_err_sub{1} = pred_err;

pred_err = [];
load 'rslts/sub_2_feat40feat40_cv1.mat'
pred_err = [pred_err;pred_test_err];
load 'rslts/sub_2_feat40feat40_cv2.mat'
pred_err = [pred_err;pred_test_err];
load 'rslts/sub_2_feat40feat40_cv3.mat'
pred_err = [pred_err;pred_test_err];
pred_err_sub{2} = pred_err;

pred_err = [];
load 'rslts/sub_3_feat40feat40_cv1.mat'
pred_err = [pred_err;pred_test_err];
load 'rslts/sub_3_feat40feat40_cv2.mat'
pred_err = [pred_err;pred_test_err];
load 'rslts/sub_3_feat40feat40_cv3.mat'
pred_err = [pred_err;pred_test_err];
pred_err_sub{3} = pred_err;

pred_err = [];
load 'rslts/sub_4_feat40feat40_cv1.mat'
pred_err = [pred_err;pred_test_err];
load 'rslts/sub_4_feat40feat40_cv2.mat'
pred_err = [pred_err;pred_test_err];
load 'rslts/sub_4_feat40feat40_cv3.mat'
pred_err = [pred_err;pred_test_err];
pred_err_sub{4} = pred_err;

for i=1:4
    avg_pred_err(i,:) = mean(pred_err_sub{i});
    std_pred_err(i,:) = std(pred_err_sub{i});
end

theta_scale = 10;
avg_pred_err(:,3) = avg_pred_err(:,3)/theta_scale;
std_pred_err(:,3) = std_pred_err(:,3)/theta_scale; 
avg_pred_err(:,6) = avg_pred_err(:,6)/theta_scale;
std_pred_err(:,6) = std_pred_err(:,6)/theta_scale; 

figure
hold on
hb = bar(1:6,avg_pred_err');
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(xData,avg_pred_err(ib,:),std_pred_err(ib,:),'k.')
end

% hb(1).FaceColor = [1 0 0];
% hb(2).FaceColor = [0.0 1 0];
% hb(3).FaceColor = [0.0 0 1];
% hb(4).FaceColor = [0.5 0.5 0];
% hb(5).FaceColor = [0.5 0 0.5];
% hb(6).FaceColor = [0.0 0.5 0.5];

Labels = {'y_p', 'z_p','theta','y_p', 'z_p','theta'};
set(gca, 'XTick', 1:6, 'XTickLabel', Labels);

legend('No submodeling','Sub2','Sub3','Sub4','Sub5','Sub6','Location','northwest');
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

