clear all;
close all;

%Start
dataFolder='~/catkin_ws/data/';
graspFolder =[dataFolder 'data_grasp/'];
crnt_trial_num = size(dir(graspFolder),1) - 1;
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

% Grasp!!
%update the state
disp('Now grasping the object..');

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

[got_obj crnt_obj_pose] = extract_obj_pose();

fig3d = figure;
plot3([table_dim(1) table_dim(2) table_dim(2) table_dim(1) table_dim(1)],[table_dim(3) table_dim(3) table_dim(4) table_dim(4) table_dim(3)],[table_dim(6) table_dim(6) table_dim(6) table_dim(6) table_dim(6)],'LineWidth',3);
axis equal;
xlabel('x');
ylabel('y');
axis([table_dim(1)-0.1 table_dim(2)+0.1 table_dim(3)-0.1 table_dim(4)+0.1]);
camroll(90);
hold on;
visualize_model(obj_name,crnt_obj_pose,[0 1 0]);
%visualize_model(obj_name,obj_start_pose,[0 1 0]);
visualize_model(obs_name,obs_pose,[1 0 0]);
patch(fv_sidegrasp_right,'FaceColor','b','FaceAlpha',0.2);
patch(fv_sidegrasp_left,'FaceColor','b','FaceAlpha',0.2);
%patch(structure_sidegrasp,'FaceColor','b','FaceAlpha',0.2);

prom_global_grasp = 'Is the object graspable? (0: no, 1: yes)';
flg_global_grasp = input(prom_global_grasp);

% generate sample points of next state
% sample a pushing position
cand_cnt = 10;

smp_r = 0.1;

push_cand = [];
for i=1:cand_cnt
    smp_theta(i) = i*((3.141592 * 2)/cand_cnt);
    push_cand(i,1) = crnt_obj_pose(1) + smp_r * cos(smp_theta(i));
    push_cand(i,2) = crnt_obj_pose(2) + smp_r * sin(smp_theta(i));
    push_cand(i,3) = crnt_obj_pose(3);

    smp_theta_rot = smp_theta(i)+ 3.141592;
    m = [0 -sin(smp_theta_rot) cos(smp_theta_rot); 0 cos(smp_theta_rot) sin(smp_theta_rot); -1 0 0];
    qw = sqrt(1 + m(1,1)+m(2,2)+m(3,3))/2;
    push_cand(i,4) = (m(3,2) - m(2,3))/(4*qw);
    push_cand(i,5) = (m(1,3) - m(3,1))/(4*qw);
    push_cand(i,6) = (m(2,1) - m(1,2))/(4*qw);
    push_cand(i,7) = qw;
    
    plot3(push_cand(i,1),push_cand(i,2),push_cand(i,3),'xr');
end

push_cand = [crnt_obj_pose;obs_pose;push_cand];
csvwrite([dataFolder,'sample_push_candidates.txt'],push_cand);

% request features of sample points
system('rosrun affordance_prediction compute_feature_from_cand');

% compute the output from model
% load the features
features_cand = csvread([dataFolder,'computed_features.txt']);

grasp_data = [];

if flg_global_grasp ==0 % no graspable points
    for i=1:cand_cnt
        flg_local_grasp(i) = 0;
    end
else
    fig_push_frame = figure;
    for i=1:cand_cnt
        visualize_pointcloud(['candidate_pointcloud_' num2str(i) '.pcd'],fig_push_frame);
        prom_grasp = 'Is the object graspable? (0: no, 1: yes)';
        flg_local_grasp(i) = input(prom_grasp);
    end
end

%grasp_data = [globally graspable? locally graspable? global_obj_coord table_dim features];
for i=1:cand_cnt
    grasp_data = [grasp_data;flg_global_grasp flg_local_grasp(i) crnt_obj_pose table_dim features_cand(i,:)];
end

%cd(graspFolder);
mkdir([graspFolder num2str(crnt_trial_num)]);
%cd(num2str(crnt_trial_num));
csvwrite([graspFolder num2str(crnt_trial_num) '/grasp_data.txt'],grasp_data);

disp('Done!');



% final_obj_pose = crnt_obj_pose;
% gripper_pose_above = [final_obj_pose(1) table_dim(3)-0.08 final_obj_pose(3)+0.05 0.7071 0 -0.7071 0];
% gripper_pose_grasp = [final_obj_pose(1) table_dim(3)-0.03 final_obj_pose(3)+0.05 0.5 0.5 -0.5 -0.5];
% in_above = in_polyhedron(fv_upright_right,gripper_pose_above(1:3));
% in_grasp = in_polyhedron(fv_sidegrasp_right,gripper_pose_grasp(1:3));
% 
% if in_above&in_grasp == 1
%     %run the grasp
%     csvwrite([dataFolder,'current_grasp_command.txt'],[gripper_pose_above gripper_pose_grasp]);
%     system('rosrun affordance_prediction execute_simple_grasp');
% end