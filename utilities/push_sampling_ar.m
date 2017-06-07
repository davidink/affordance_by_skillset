clear all;
close all;

%Start
dataFolder='~/catkin_ws/data/';
dataSaveFolder=[dataFolder 'data_push_result/'];
%addpath(dataFolder);
%load 'PR2_boundary.mat'
load([dataFolder 'models/ARmodels/pair_map.mat']);

current_exp_num = size(dir(dataSaveFolder),1)-2 +1;
current_save_folder = [dataSaveFolder num2str(current_exp_num) '/'];
mkdir(current_save_folder);

% get the starting state
got_scene = false;
disp('Acquiring the initial scene..');
[status_recognition result_recognition] = system(['rosrun ar_track_alvar request_obj_ids']);
ar_ids = load([dataFolder 'current_obj_ids.csv']);

%convert into global frames
obj_poses_wID=[];
for i=1:size(ar_ids,1)
    %crnt_pose = convert_ar_to_obj_pose(ar_ids(i,1),ar_ids(i,2:8));
    obj_poses_wID(i,:) = [ar_ids(i,1) convert_ar_to_obj_pose(ar_ids(i,1),ar_ids(i,2:8))];
end

fig_table = figure;
hold on;

visualize_objs(fig_table, obj_poses_wID, [0 1 0]);

push_start_pose = [0.6 -0.3 0.75 -0.5 0.5 0.5 0.5]; %x y z qw qx qy qz
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