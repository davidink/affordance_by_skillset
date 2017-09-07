close all
clear all
%dataFolder = 'data/trajectories/push_1_object_4/';
testdataFolder = 'data/test/';

% monolithic
num_test_scene = size(dir([testdataFolder 'test*']),1);


for i=1:num_test_scene
    obj1_pred = [];
    obj2_pred = [];
    obj1_gt = [];
    obj2_gt = [];
    for k=1:5
        load([testdataFolder 'mono_prediction_result_' num2str(i) '_' num2str(k) '.mat']);
        obj1_pred = [obj1_pred;prediction(1,:)];
        obj2_pred = [obj2_pred;prediction(2,:)];
        obj1_gt = [obj1_gt;ground_truth(1,:)];
        obj2_gt = [obj2_gt;ground_truth(2,:)];
    end
    obj1_diff = obj1_pred-obj1_gt;
    for k=1:5
        for j=4:6
            if obj1_diff(k,j) > pi
                obj1_diff(k,j) = obj1_diff(k,j) - 2*pi;
            else if obj1_diff(k,j) < -pi
                obj1_diff(k,j) = obj1_diff(k,j) + 2*pi;
                end
            end
        end
    end
    mono_obj1_pred_rms(i,:) = rms(obj1_diff);
    obj2_diff = obj2_pred-obj2_gt;
    for k=1:5
        for j=4:6
            if obj2_diff(k,j) > pi
                obj2_diff(k,j) = obj2_diff(k,j) - 2*pi;
            else if obj2_diff(k,j) < -pi
                obj2_diff(k,j) = obj2_diff(k,j) + 2*pi;
                end
            end
        end
    end
    mono_obj1_pred_rms(i,:) = rms(obj1_diff);
    mono_obj2_pred_rms(i,:) = rms(obj2_diff);
end

mono_obj1_pred_err = mean(mono_obj1_pred_rms);
mono_obj1_pred_std = std(mono_obj1_pred_rms);
mono_obj2_pred_err = rms(mono_obj2_pred_rms);
mono_obj2_pred_std = std(mono_obj2_pred_rms);

% obj1_pred = [];
% obj2_pred = [];
% obj1_gt = [];
% obj2_gt = [];
for i=1:num_test_scene
    load([testdataFolder 'modular_prediction_result_' num2str(i) '.mat']);
%     obj1_pred_err(i,:) = mean(obj_trajectory{1}-obj_trajectory_gt{1}(1:size(obj_trajectory{1},1),:));
%     obj2_pred_err(i,:) = mean(obj_trajectory{2}-obj_trajectory_gt{2}(1:size(obj_trajectory{1},1),:));
    obj1_pred_err(i,:) = (obj_trajectory{1}(end,:)-obj_trajectory_gt{1}(size(obj_trajectory{1},1),:));
    obj2_pred_err(i,:) = (obj_trajectory{2}(end,:)-obj_trajectory_gt{2}(size(obj_trajectory{2},1),:));
    for j=4:6
        if obj1_pred_err(i,j) > pi
            obj1_pred_err(i,j) = obj1_pred_err(i,j) - 2*pi;
        else if obj1_pred_err(i,j) < -pi
            obj1_pred_err(i,j) = obj1_pred_err(i,j) + 2*pi;
            end
        end
    end
end

modu_obj1_pred_err = rms(obj1_pred_err);
modu_obj1_pred_std = std(obj1_pred_err);
modu_obj2_pred_err = rms(obj2_pred_err);
modu_obj2_pred_std = std(obj2_pred_err);

avg_pred_err = [mono_obj1_pred_err;modu_obj1_pred_err;mono_obj2_pred_err;modu_obj2_pred_err];
avg_std_err = [mono_obj1_pred_std;modu_obj1_pred_std;mono_obj2_pred_std;modu_obj2_pred_std];

theta_scale = 10;
avg_pred_err(:,4:6) = avg_pred_err(:,4:6)/theta_scale;
avg_std_err(:,4:6) = avg_std_err(:,4:6)/theta_scale; 

figure;
hold on
hb = bar(1:6,avg_pred_err');
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(xData,avg_pred_err(ib,:),avg_std_err(ib,:),'k.')
end

% hb(1).FaceColor = [1 0 0];
% hb(2).FaceColor = [0.0 1 0];
% hb(3).FaceColor = [0.0 0 1];
% hb(4).FaceColor = [0.5 0.5 0];
% hb(5).FaceColor = [0.5 0 0.5];
% hb(6).FaceColor = [0.0 0.5 0.5];

Labels = {'x', 'y','z','rx', 'ry','rz'};
set(gca, 'XTick', 1:6, 'XTickLabel', Labels);

legend('mono obj1','modular obj1','mono obj2','modular obj2','Location','northwest');


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

