clear all;
close all;

%Start
dataFolder='~/catkin_hydro_ws/data/';
% if getenv('OS')=='Windows_NT'
%     dataFolder='c:/Users/David/Dropbox/Research/Code/data/';    
% end
dataSaveFolder=[dataFolder 'data_push_result/'];
drawArrow3 = @(x,y,varargin) quiver3( x(1),x(2),x(3),y(1)-x(1),y(2)-x(2),y(3)-x(3),varargin{:} );

%addpath(dataFolder);
%load 'PR2_boundary.mat'
load([dataFolder 'models/ARmodels/pair_map.mat']);

load('exp/exp_push_4-2.mat');
%load('exp/exp_push2-1.mat');

%obj_poses_wID = push_predict_dis{min_idx(1),min_idx(2)};

fig_table = figure;
hold on;
visualize_objs(fig_table, obj_poses_wID, [0 1 0]);
%visualize_objs(fig_table, obj_poses_wID_ap, [0 0 1]);


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
visualize_objs(fig_table, obj_poses_wID, [0 0 1]);

camroll(-90);
axis([0.4 1.1 -0.4 0.4]);
drawArrow3(push_command{min_idx(1),min_idx(2)}(1,1:3),push_command{min_idx(1),min_idx(2)}(2,1:3),'LineWidth',2); 
%visualize_objs(fig_table, obj_poses_wID_ap, [0 0 1]);

%execute push
push_ang = atan2(push_command{min_idx(1),min_idx(2)}(2,2)-push_command{min_idx(1),min_idx(2)}(1,2),push_command{min_idx(1),min_idx(2)}(2,1)-push_command{min_idx(1),min_idx(2)}(1,1));
R = [1 0 0;0 cos(push_ang+pi/2) -sin(push_ang+pi/2);0 sin(push_ang+pi/2) cos(push_ang+pi/2)];
r_bq = qGetR([-0.5 0.5 0.5 0.5]);
push_q = qGetQ(r_bq*R);
push_start_pose = [push_command{min_idx(1),min_idx(2)}(1,1:3) push_q'];

push_end_pose = [push_command{min_idx(1),min_idx(2)}(2,1:3) push_q'];

% execute the push
current_push = [push_start_pose push_end_pose];
csvwrite([dataFolder,'current_push_command.txt'],current_push);
disp('Executing the push..');
[push_status, push_result] = system('rosrun affordance_prediction execute_simple_push');
visualize_objs(fig_table, push_predict_dis{min_idx(1),min_idx(2)}, [1 0 0]);

got_scene = false;
disp('Acquiring the initial scene..');
[status_recognition result_recognition] = system(['rosrun ar_track_alvar request_obj_ids']);
ar_ids = load([dataFolder 'current_obj_ids.csv']);

%convert into global frames
obj_poses_wID_ap=[];
for i=1:size(ar_ids,1)
    %crnt_pose = convert_ar_to_obj_pose(ar_ids(i,1),ar_ids(i,2:8));
    obj_poses_wID_ap(i,:) = [ar_ids(i,1) convert_ar_to_obj_pose(ar_ids(i,1),ar_ids(i,2:8))];
end
visualize_objs(fig_table, obj_poses_wID_ap, [0 1 1]);

