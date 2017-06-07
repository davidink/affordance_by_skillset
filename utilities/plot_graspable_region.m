clear all;
close all;

%Start
dataFolder='~/catkin_ws/data/';
graspFolder =[dataFolder 'data_grasp/'];
crnt_trial_num = size(dir(graspFolder),1) - 1;
%addpath(dataFolder);
load 'PR2_boundary.mat'

obj_obs_pose = csvread([dataFolder,'start_obj_obs_pose.txt']);
obj_start_pose = obj_obs_pose(1,:);
obs_pose = obj_obs_pose(2,:);

obj_name = 'Book1';
obs_name = 'Battery1';

table_dim = csvread([dataFolder,'start_table_dimension.txt']);

fig3d = figure;
plot3([table_dim(1) table_dim(2) table_dim(2) table_dim(1) table_dim(1)],[table_dim(3) table_dim(3) table_dim(4) table_dim(4) table_dim(3)],[table_dim(6) table_dim(6) table_dim(6) table_dim(6) table_dim(6)],'LineWidth',3);
axis equal;
xlabel('x');
ylabel('y');
axis([table_dim(1)-0.1 table_dim(2)+0.1 table_dim(3)-0.1 table_dim(4)+0.1]);
camroll(90);
hold on;
%visualize_model(obj_name,crnt_obj_pose,[0 1 0]);
%visualize_model(obj_name,obj_start_pose,[0 1 0]);
%visualize_model(obs_name,obs_pose,[1 0 0]);
patch(fv_sidegrasp_right,'FaceColor','b','FaceAlpha',0.2);
patch(fv_sidegrasp_left,'FaceColor','b','FaceAlpha',0.2);
%patch(structure_sidegrasp,'FaceColor','b','FaceAlpha',0.2);
