clear all;
close all;

%Start
dataFolder='~/catkin_ws/data/';
%addpath(dataFolder);
load 'PR2_boundary.mat'

% get the starting state
system('rosrun affordance_prediction acquire_obj_obs_pose');
% pause(3);

obj_obs_pose = csvread([dataFolder,'start_obj_obs_pose.txt']);
obj_start_pose = obj_obs_pose(1,:);
obs_pose = obj_obs_pose(2,:);

obj_name = 'Book1';
obs_name = 'Battery1';

table_dim = csvread([dataFolder,'start_table_dimension.txt']);

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