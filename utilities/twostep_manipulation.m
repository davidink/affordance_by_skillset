clear all;
close all;

%Start
dataFolder='~/catkin_ws/data/';
% if getenv('OS')=='Windows_NT'
%     dataFolder='c:/Users/David/Dropbox/Research/Code/data/';    
% end
dataSaveFolder=[dataFolder 'data_push_result/'];

%addpath(dataFolder);
%load 'PR2_boundary.mat'
load([dataFolder 'models/ARmodels/pair_map.mat']);
load('WLRWGP_mdl47.mat');
set_num = 7;
num_submodel = 4;
if set_num == 1
    feat_WLR_idx = [1:8]; % object features 
    feat_WGP_idx = [1:8]; % object features
end
if set_num == 2
    feat_WLR_idx = [1:21]; % obj+contact 
    feat_WGP_idx = [1:21]; 
end
if set_num == 3
    feat_WLR_idx = [22:101]; % context only
    feat_WGP_idx = [22:101]; 
end
if set_num == 4
    feat_WLR_idx = [1:101]; % obj+contact+context 
    feat_WGP_idx = [1:101]; 
end
if set_num == 5
    feat_WLR_idx = [1:21]; % obj+contact 
    feat_WGP_idx = [22:101]; %ctxt only
end
if set_num == 6
    feat_WLR_idx = [1:21]; % obj+contact 
    feat_WGP_idx = [1:101]; %all 
end
if set_num == 7
    feat_WLR_idx = [22:101]; % obj+contact 
    feat_WGP_idx = [1:21]; %obj only
end
if set_num == 8
    feat_WLR_idx = [22:101]; % obj+contact 
    feat_WGP_idx = [1:101]; % all
end

%Start
%load 'PR2_boundary.mat'

% get the starting state
% got_scene = false;
% disp('Acquiring the initial scene..');
% % [status_recognition result_recognition] = system(['rosrun ar_track_alvar request_obj_ids']);
% ar_ids = load([dataFolder 'current_obj_ids.csv']);
% 
% %convert into global frames
% obj_poses_wID=[];
% for i=1:size(ar_ids,1)
%     %crnt_pose = convert_ar_to_obj_pose(ar_ids(i,1),ar_ids(i,2:8));
%     obj_poses_wID(i,:) = [ar_ids(i,1) convert_ar_to_obj_pose(ar_ids(i,1),ar_ids(i,2:8))];
% end

load('exp/exp_push_5.mat');
obj_poses_wID = push_predict_dis{min_idx(1),min_idx(2)};

% load([dataSaveFolder '70/push_result.csv']);
% obj_poses_wID = push_result(1:2,:);
% obj_poses_wID_ap = push_result(3:4,:);
%load([dataSaveFolder '35/push_command.csv']);



fig_table = figure;
hold on;
visualize_objs(fig_table, obj_poses_wID, [0 1 0]);
%visualize_objs(fig_table, obj_poses_wID_ap, [0 0 1]);
camroll(-90);
axis([0.3 0.9 -0.4 0.4]);

%set targets
tar_id = 1;
sec_id = 2;

% sample candidate pushes
cand_cnt1 = 20;
cand_cnt2 = 5;
push_cand = [];
push_r = 0.25;
drawArrow3 = @(x,y,varargin) quiver3( x(1),x(2),x(3),y(1)-x(1),y(2)-x(2),y(3)-x(3),varargin{:} );

for i=1:cand_cnt1
    % sample random starting position
    smp_r(i) = rand/10 + 0.15; % r = 0.13 ~ 0.23m
    smp_theta1(i) = rand * pi * 2; % theta = 0 ~ 360 degree;
    cand_push_start(i,1) = obj_poses_wID(tar_id,2) +smp_r(i) * cos(smp_theta1(i)); % x
    cand_push_start(i,2) = obj_poses_wID(tar_id,3) +smp_r(i) * sin(smp_theta1(i)); % y
    cand_push_start(i,3) = obj_poses_wID(tar_id,4); % z
    for j=1:cand_cnt2
        smp_theta2(j) = rand*pi/2 + 0.75*pi + smp_theta1(i);
        cand_push_end(j,1) = cand_push_start(i,1) + push_r*cos(smp_theta2(j));
        cand_push_end(j,2) = cand_push_start(i,2) + push_r*sin(smp_theta2(j));
        cand_push_end(j,3) = cand_push_start(i,3);
        smp_theta_rot = smp_theta2(j)+ pi;
        m = [0 -sin(smp_theta_rot) cos(smp_theta_rot); 0 cos(smp_theta_rot) sin(smp_theta_rot); -1 0 0];
        qw = sqrt(1 + m(1,1)+m(2,2)+m(3,3))/2;
        cand_push_end(j,5) = (m(3,2) - m(2,3))/(4*qw);
        cand_push_end(j,6) = (m(1,3) - m(3,1))/(4*qw);
        cand_push_end(j,7) = (m(2,1) - m(1,2))/(4*qw);
        cand_push_end(j,4) = qw;
        push_command{i,j} = [[cand_push_start(i,1:3) cand_push_end(j,4:7)]; cand_push_end(j,1:7)]; 

        % Draw candidate push
        figure(fig_table);
        drawArrow3([cand_push_start(i,1:3)],[cand_push_end(j,1:3)],'LineWidth',2);
        
    end
end

[pointcloud_bp pointcloud_norm_bp] = load_pointcloud(obj_poses_wID); %before push
for i=1:cand_cnt1
    for j=1:cand_cnt2
        % obtain pointcloud from model with pose
%         fig_test = figure;
%         hold on;
%         visualize_objs(fig_test, obj_poses_wID, [0 1 0]);
%         drawArrow3(push_command{i,j}(1,1:3),push_command{i,j}(2,1:3),'LineWidth',2);        
        for k=1:size(pointcloud_bp,2)
            % rotate pointcloud to pushing frame
            [pf_transform pointcloud_bp_pf{k}] = get_pf_transform(pointcloud_bp{k}, push_command{i,j});
            normalpoints = [pf_transform(1:3,1:3) [0 0 0]';0 0 0 1]*[pointcloud_norm_bp{k}'; ones(1,size(pointcloud_norm_bp{k},1))];
            pointcloud_norm_bp_pf{k}=normalpoints(1:3,:)';
            obj_poses_pf_bp{k} = pf_transform * [qGetR(obj_poses_wID(k,5:8)) obj_poses_wID(k,2:4)'; 0 0 0 1];            
        end
%         if obj_poses_pf_bp{1}(3,4) < obj_poses_pf_bp{2}(3,4) 
%             if obj_poses_pf_bp{1}(3,4) < 0
%                 tar_id = 2;
%                 sec_id =1;
%             else
%                 tar_id = 1;
%                 sec_id = 2;
%             end
%         else if obj_poses_pf_bp{2}(3,4) < 0 % this object is behind the push
%                 tar_id = 1;
%                 sec_id = 2;
%             else
%                 tar_id = 2;
%                 sec_id = 1;        
%             end
%         end
        % extract features
        features = extract_feature(pointcloud_bp_pf, pointcloud_norm_bp_pf, obj_poses_pf_bp, push_command);
        WLR_pred = mdl_WLR.computeLikelihood(features(:,feat_WLR_idx),[]);
        for a=1:6
            for m=1:num_submodel
                [WGP_pred{a,m}.pred_mu WGP_pred{a,m}.pred_ss] = wgp_pred(mdl_WGP{a,m}.hyp, mdl_WGP{a,m}.weight, mdl_WGP{a,m}.input, mdl_WGP{a,m}.target, features(:,feat_WGP_idx));
            end
        end
        pred_label = [];
        for m=1:num_submodel
            if WLR_pred(1,m) == max(WLR_pred(1,:))
                pred_label = m;
            end
        end
        pred_mu = zeros(1,6);
        for k=1:6
            pred_mu(1,k) = WGP_pred{k,pred_label}.pred_mu;
            pred_ss(1,k) = WGP_pred{k,pred_label}.pred_ss;
            if isnan(pred_mu(1,k))
                halt=0;
            end
            if k==5 & pred_mu(1,k) < 0
                pred_mu(1,k) = 0;
            end
        end
        % check if thing is behind the push, it shouldn't move
        for k=1:size(pointcloud_bp,2)
            if obj_poses_pf_bp{k}(3,4) < 0 % behind the push
                pred_mu(1,4) = 0;
                pred_mu(1,5) = 0;
                pred_mu(1,6) = 0;
            end
        end        
        % project expected place of pushes
        for k=1:size(pointcloud_bp,2)
            if k== tar_id
                obj_poses_pf_ap_pred{k}(1,4) = obj_poses_pf_bp{tar_id}(1,4); % x1        
                obj_poses_pf_ap_pred{k}(2,4) = obj_poses_pf_bp{tar_id}(2,4) + pred_mu(1); % y1
                obj_poses_pf_ap_pred{k}(3,4) = obj_poses_pf_bp{tar_id}(3,4) + pred_mu(2); % z1
                R = [1 0 0; 0 cos(-pred_mu(3)) -sin(-pred_mu(3));0 sin(-pred_mu(3)) cos(-pred_mu(3))];
                obj_poses_pf_ap_pred{k}(1:3,1:3) = R*obj_poses_pf_bp{k}(1:3,1:3);
                new_rot(k) = atan2(obj_poses_pf_ap_pred{k}(2,1),obj_poses_pf_ap_pred{k}(3,1));
                obj_poses_pf_ap_pred{k}(4,:) = [0 0 0 1];
            else
                obj_poses_pf_ap_pred{k}(1,4) = obj_poses_pf_bp{sec_id}(1,4); % x1        
                obj_poses_pf_ap_pred{k}(2,4) = obj_poses_pf_bp{sec_id}(2,4) + pred_mu(4); % y1
                obj_poses_pf_ap_pred{k}(3,4) = obj_poses_pf_bp{sec_id}(3,4) + pred_mu(5); % z1
                R = [1 0 0; 0 cos(-pred_mu(6)) -sin(-pred_mu(6));0 sin(-pred_mu(6)) cos(-pred_mu(6))];
                obj_poses_pf_ap_pred{k}(1:3,1:3) = R*obj_poses_pf_bp{k}(1:3,1:3);
                new_rot(k) = atan2(obj_poses_pf_ap_pred{k}(2,1),obj_poses_pf_ap_pred{k}(3,1));
                obj_poses_pf_ap_pred{k}(4,:) = [0 0 0 1];
            end
            RT = pf_transform^-1 * obj_poses_pf_ap_pred{k};
            obj_poses_wID_ap_pred(k,:) = [obj_poses_wID(k,1) RT(1:3,4)' qGetQ(RT(1:3,1:3))'];            
        end        
%         visualize_objs(fig_test, obj_poses_wID_ap_pred, [0 0 1]);
%         camroll(-90);
%         axis([0.3 0.9 -0.4 0.4]);
        push_ang = atan2(push_command{i,j}(2,2)-push_command{i,j}(2,1),push_command{i,j}(1,2)-push_command{i,j}(1,1));
        push_predict_rot(i,j) = new_rot(sec_id)-push_ang;
        push_predict_dis{i,j} = obj_poses_wID_ap_pred;
    end
end

% search for second object to be at certain degree as possible
% min_x = 100;
% for i=1:cand_cnt1
%     for j=1:cand_cnt2
%         dist = sqrt((push_predict_dis{i,j}(tar_id,2)-push_predict_dis{i,j}(sec_id,2))^2+(push_predict_dis{i,j}(tar_id,3)-push_predict_dis{i,j}(sec_id,3))^2);
%         if min_x > dist & dist > 0.18
%             min_x = dist;
%             min_idx = [i;j];
%         end
%     end
% end

% max_y = -100;
% for i=1:cand_cnt1
%     for j=1:cand_cnt2
%         if push_command{i,j}(1,1) < 0.74 & push_command{i,j}(2,1) < 0.74
%         %dist = sqrt((push_predict_dis{i,j}(tar_id,2)-push_predict_dis{i,j}(sec_id,2))^2+(push_predict_dis{i,j}(tar_id,3)-push_predict_dis{i,j}(sec_id,3))^2);
%             if max_y < push_predict_dis{i,j}(sec_id,3)
%                 max_y = push_predict_dis{i,j}(sec_id,3);
%                 min_idx = [i;j];
%             end
%         end
%     end
% end

min_x = 100;
for i=1:cand_cnt1
    for j=1:cand_cnt2
        if push_command{i,j}(1,1) < 0.74 & push_command{i,j}(2,1) < 0.74
        %dist = sqrt((push_predict_dis{i,j}(tar_id,2)-push_predict_dis{i,j}(sec_id,2))^2+(push_predict_dis{i,j}(tar_id,3)-push_predict_dis{i,j}(sec_id,3))^2);
            if obj_poses_wID(sec_id,2) > push_predict_dis{i,j}(sec_id,2)
                min_x = push_predict_dis{i,j}(sec_id,2);
                min_idx = [i;j];
            end
        end
    end
end

% desired_rot = pi/2;
% min_err = 100;
% for i=1:cand_cnt1
%     for j=1:cand_cnt2
%         if min_err > abs(push_predict_rot(i,j)-desired_rot)
%             min_err = abs(push_predict_rot(i,j)-desired_rot);
%             min_idx = [i;j];
%         end
%     end
% end

fig_best = figure;
hold on;
visualize_objs(fig_best, obj_poses_wID, [0 1 0]);
drawArrow3(push_command{min_idx(1),min_idx(2)}(1,1:3),push_command{min_idx(1),min_idx(2)}(2,1:3),'LineWidth',2); 
visualize_objs(fig_best, push_predict_dis{min_idx(1),min_idx(2)}, [0 0 1]);
camroll(-90);
axis([0.3 0.9 -0.4 0.4]);

min_idx=[13;5];
fig_best = figure;
hold on;
visualize_objs(fig_best, obj_poses_wID, [0 1 0]);
drawArrow3(push_command{min_idx(1),min_idx(2)}(1,1:3),push_command{min_idx(1),min_idx(2)}(2,1:3),'LineWidth',2); 
visualize_objs(fig_best, push_predict_dis{min_idx(1),min_idx(2)}, [0 0 1]);
camroll(-90);
axis([0.3 0.9 -0.4 0.4]);



push_start_pose = [0.6 -0.3 0.75 -0.5 0.5 0.5 0.5]; %x y z qx qy qz qw
push_length = 0.25;
push_end_pose = push_start_pose;
push_end_pose(2) = push_end_pose(2) + push_length;

% execute the push
current_push = [push_start_pose push_end_pose];
csvwrite([dataFolder,'current_push_command.txt'],current_push);
disp('Executing the push..');
[push_status, push_result] = system('rosrun affordance_prediction execute_simple_push');

if size(strfind(push_result,'successful'),1) == 0
    disp('Pushing was not successful..');
else
    disp('Pushing successful!');
    %update the state
    disp('Acquiring the scene..');
    [status_recognition result_recognition] = system(['rosrun ar_track_alvar request_obj_ids']);
    ar_ids_ap = load([dataFolder 'current_obj_ids.csv']);

    %convert into global frames
    obj_poses_wID_ap=[];
    for i=1:size(ar_ids,1)
        %crnt_pose = convert_ar_to_obj_pose(ar_ids(i,1),ar_ids(i,2:8));
        obj_poses_wID_ap(i,:) = [ar_ids_ap(i,1) convert_ar_to_obj_pose(ar_ids_ap(i,1),ar_ids_ap(i,2:8))];
    end

    visualize_objs(fig_table, obj_poses_wID_ap, [0 0 1]);
    axis equal
    camroll(-90);
    
    csvwrite([current_save_folder,'push_command.csv'],[push_start_pose;push_end_pose]);
    csvwrite([current_save_folder,'push_result.csv'],[obj_poses_wID;obj_poses_wID_ap]);
end




% set the goal state
%goal_state = [0.45 -0.48 0.74]; % edge of the table
% load 'SVM_model_grasp_global_feat.mat';
% SVMmodle_grasp = SVMModel;
goal_pose = gen_grasp_pose(table_dim);
goal_state = [table_dim(1)+(table_dim(2)-table_dim(1))/4 table_dim(3)+0.05 table_dim(6)];
%goal_state = goal_pose(1:3);

manipulation_done = false;
crnt_obj_pose = obj_start_pose;
fig3d = figure;
replan_cnt = 0;

while ~manipulation_done
    
    % compute the plan
    disp('Computing a plan..');
    plan_set = compute_plan(crnt_obj_pose,obs_pose,goal_state,table_dim);
    %load 'plan_set.mat'

    close(figure(fig3d));
    fig3d = figure;
    plot3([table_dim(1) table_dim(2) table_dim(2) table_dim(1) table_dim(1)],[table_dim(3) table_dim(3) table_dim(4) table_dim(4) table_dim(3)],[table_dim(6) table_dim(6) table_dim(6) table_dim(6) table_dim(6)],'LineWidth',3);
    axis equal;
    xlabel('x(m)');
    ylabel('y(m)');
    axis([table_dim(1)-0.1 table_dim(2)+0.1 table_dim(3)-0.1 table_dim(4)+0.1]);
    camroll(90);
    hold on;
    visualize_model(obj_name,crnt_obj_pose,[0 1 0]);
    visualize_model(obs_name,obs_pose,[1 0 0]);
    %patch(fv_upright_right,'FaceColor','g','FaceAlpha',0.2);
    %patch(fv_upright_left,'FaceColor','g','FaceAlpha',0.2);
    %patch(structure_sidegrasp,'FaceColor','b','FaceAlpha',0.2);
    %figure(fig3d);
    for i=1:size(plan_set,1)
        visualize_model(obj_name, plan_set(i,5:11), [0 1-0.05*i 0]);
    end

    % Execute step by step and check 
  for i=2:size(plan_set,1)
        %if i ==1 % start state
        %else
        current_push = [plan_set(i,12:18) plan_set(i-1,5:7) plan_set(i,15:18)];
        csvwrite([dataFolder,'current_push_command.txt'],current_push);

        % execute push
        disp('Executing the push..');
        push_success = 0;
        [push_status, push_result] = system('rosrun affordance_prediction execute_simple_push');
        if size(strfind(push_result,'successful'),1) == 0
            disp('Pushing was not successful, re-planning..');
            push_success = 0;
        else push_success =1;
        end
        
        %update the state
        disp('Updating the scene..');
        updated = false;
        while ~updated
            [status,result] = system('rosrun affordance_prediction save_tabletop_pointcloud');
            err = strfind(result,'error');
            if size(err,1) ==0
                updated = true;
                disp('Scene updated');
            end
        end
        [got_obj rslt_obj_pose(i-1,:)] = extract_obj_pose();
        %system('rosrun affordance_prediction acquire_obj_obs_pose');
        %obj_obs_pose = csvread([dataFolder,'current_obj_obs_pose.txt']);
        %rslt_obj_pose(i-1,:) = obj_obs_pose(1,:);
        %obs_crnt_pose = obj_obs_pose(2,:);

        figure(fig3d);
        visualize_model(obj_name, rslt_obj_pose(i-1,:), [0 0 1]);
        crnt_err(i-1) = rms(rslt_obj_pose(i-1,1:3)-plan_set(i,5:7))
        
        if push_success == 0
            crnt_obj_pose = rslt_obj_pose(i-1,:);
            disp('Push failed, re-planning...');
            save(['plan_execution_result_' num2str(replan_cnt) '.mat'],'plan_set','crnt_err','rslt_obj_pose');
            replan_cnt = replan_cnt +1;
            break;
        end
        
        if rms(rslt_obj_pose(i-1,1:3)-goal_state) < 0.05 
            manipulation_done = true;
            disp('Reached the goal position!!');
            save(['plan_execution_result_' num2str(replan_cnt) '.mat'],'plan_set','crnt_err','rslt_obj_pose');
            crnt_obj_pose = rslt_obj_pose(i-1,:);
            break;
        end
        
        if crnt_err(i-1) > 0.05
            crnt_obj_pose = rslt_obj_pose(i-1,:);
            disp('Error over threshold, re-planning...');
            save(['plan_execution_result_' num2str(replan_cnt) '.mat'],'plan_set','crnt_err','rslt_obj_pose');
            replan_cnt = replan_cnt +1;
            break;
        end
        pause(0.05);
    end    
    
    % check whether the object is in the region to grasp
    % execute the grasping planner
    
end

% Grasp!!
%update the state
disp('Now grasping the object..');
disp('Updating the scene..');
% updated = false;
% while ~updated
%     [status,result] = system('rosrun affordance_prediction save_tabletop_pointcloud');
%     err = strfind(result,'error');
%     if size(err,1) ==0
%         updated = true;
%         disp('Scene updated');
%     end
% end
final_obj_pose = crnt_obj_pose;
gripper_pose_above = [final_obj_pose(1) table_dim(3)-0.08 final_obj_pose(3)+0.05 0.7071 0 -0.7071 0];
gripper_pose_grasp = [final_obj_pose(1) table_dim(3)-0.03 final_obj_pose(3)+0.05 0.5 0.5 -0.5 -0.5];
in_above = in_polyhedron(fv_upright_right,gripper_pose_above(1:3));
in_grasp = in_polyhedron(fv_sidegrasp_right,gripper_pose_grasp(1:3));

if in_above&in_grasp == 1
    %run the grasp
    csvwrite([dataFolder,'current_grasp_command.txt'],[gripper_pose_above gripper_pose_grasp]);
    system('rosrun affordance_prediction execute_simple_grasp');
end
