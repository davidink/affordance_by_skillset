%function [feat rslt] = load_test_dataset()
clear all
close all

addpath('utilities');
dataFolder='~/catkin_ws/data/';
dataSaveFolder=[dataFolder 'data_push_result/'];
load([dataFolder 'models/ARmodels/pair_map.mat']);

current_exp_num = size(dir(dataSaveFolder),1)-2;

fig_table = figure;
hold on;

for i=1:current_exp_num
    disp(['loading ' num2str(i) 'th experiment..']);
    load([dataSaveFolder num2str(i) '/push_command.csv']);
    load([dataSaveFolder num2str(i) '/push_result.csv']);
    num_obj = size(push_result,1)/2;
    obj_poses_wID = push_result(1:num_obj,:);
    obj_poses_wID_ap = push_result(num_obj+1:end,:);
    
    %align ids
    if obj_poses_wID_ap(1,1) ~= obj_poses_wID(1,1)
        %swap lines
        temp = obj_poses_wID_ap(1,:);
        obj_poses_wID_ap(1,:) = obj_poses_wID_ap(2,:);
        obj_poses_wID_ap(2,:) = temp;
    end
    
    % obtain pointcloud from model with pose
    [pointcloud_bp pointcloud_norm_bp] = load_pointcloud(obj_poses_wID); %before push
    %[pointcloud_ap pointcloud_norm_ap] = load_pointcloud(obj_poses_wID_ap); %after push

    clf(fig_table);
    visualize_objs(fig_table, obj_poses_wID, [0 1 0]); %before push
    %visualize_objs(fig_table, obj_poses_wID_ap, [0 0 1]); %after push
    view(0,90)  % XY
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal;

    % load the gripper
    load('PR2_gripper.mat');
    gripper_points= modelpoints;
    gripper_norms = normpoints;
    
    start_gripper_pose = [qGetR(push_command(1,4:7)) push_command(1,1:3)'; 0 0 0 1];
    
    modelpoints=start_gripper_pose*[gripper_points'; ones(1,size(gripper_points,1))];
    %normalpoints = crnt_obj_pose*[normalPoints'; ones(1,size(normalPoints,1))];
    normalpoints = [start_gripper_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[gripper_norms'; ones(1,size(gripper_norms,1))];
    gripper_points=modelpoints(1:3,:)';
    gripper_norms=normalpoints(1:3,:)';
    
    plot3(gripper_points(:,1),gripper_points(:,2),gripper_points(:,3),'Color',[1 0 0],'Marker','.','Linestyle','none');
    hold on;
    plotCoord(start_gripper_pose(1:3,4)',start_gripper_pose(1:3,1:3),0.025);
    
    if sum((obj_poses_wID(1,2:4)-push_command(1,1:3)).^2) < sum((obj_poses_wID(2,2:4)-push_command(1,1:3)).^2)
        clsr_obj_label = 1;
    else
        clsr_obj_label = 2;
    end
    
    %label random points whether they are in contact or not
    features = [];
    labels = [];
    point_idx = randi(size(pointcloud_bp{clsr_obj_label},1),10,1);
    %plot3(pointcloud_bp{clsr_obj_label}(point_idx,1),pointcloud_bp{clsr_obj_label}(point_idx,2),pointcloud_bp{clsr_obj_label}(point_idx,3),'color',[1 0 1],'Marker','x','Linestyle','none');
    for j=1:size(point_idx,1)
        feature = compute_contact_feat_lr(pointcloud_bp{clsr_obj_label}(point_idx(j),:), pointcloud_norm_bp{clsr_obj_label}(point_idx(j),:),gripper_points,gripper_norms);
        plot3(pointcloud_bp{clsr_obj_label}(point_idx(j),1),pointcloud_bp{clsr_obj_label}(point_idx(j),2),pointcloud_bp{clsr_obj_label}(point_idx(j),3),'color',[1 0 1],'Marker','x','Linestyle','none','MarkerSize',10);
        prompt = 'is this point in contact? 0:no 1:yes';
        contact_label = input(prompt);
        features = [features;feature];
        labels = [labels;contact_label];
    end
    
    %label probable points of contacts
    push_dir = push_command(2,1:3) - push_command(1,1:3);
    push_dir = push_dir./sqrt(sum(push_dir.^2));
    
    eps = 0.1;
    norm_idx = [];
    for j=1:size(pointcloud_bp{clsr_obj_label},1)
        diff = abs(push_dir-pointcloud_norm_bp{clsr_obj_label}(j,:));
        cross_val = cross(push_dir,pointcloud_norm_bp{clsr_obj_label}(j,:));
        if push_dir * pointcloud_norm_bp{clsr_obj_label}(j,:)' < -0.45
            norm_idx = [norm_idx;j];
        end
    end
%     plot3(pointcloud_bp{clsr_obj_label}(norm_idx,1),pointcloud_bp{clsr_obj_label}(norm_idx,2),pointcloud_bp{clsr_obj_label}(norm_idx,3),'color',[0 1 1],'Marker','o','Linestyle','none','MarkerSize',10);
    
    point_idx = randi(size(norm_idx,1),10,1);
    for j=1:size(point_idx,1)
        feature = compute_contact_feat_lr(pointcloud_bp{clsr_obj_label}(norm_idx(point_idx(j)),:), pointcloud_norm_bp{clsr_obj_label}(norm_idx(point_idx(j)),:),gripper_points,gripper_norms);
        plot3(pointcloud_bp{clsr_obj_label}(norm_idx(point_idx(j)),1),pointcloud_bp{clsr_obj_label}(norm_idx(point_idx(j)),2),pointcloud_bp{clsr_obj_label}(norm_idx(point_idx(j)),3),'color',[1 0 1],'Marker','*','Linestyle','none','MarkerSize',15);
        prompt = 'is this point in contact? 0:no 1:yes';
        contact_label = input(prompt);
        features = [features;feature];
        labels = [labels;contact_label];
    end
    
    
    
%     plotCoord(crnt_obj_pose(1:3,4)',crnt_obj_pose(1:3,1:3),0.025);
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     axis equal
%     
    save(['feat_label_contact_' num2str(i) '.mat'],'features','labels');
    pause(0.1);
    
end

% save features