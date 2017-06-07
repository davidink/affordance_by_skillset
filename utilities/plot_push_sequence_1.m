clear all;
close all;

%Start
dataFolder='~/catkin_ws/data/';
% if getenv('OS')=='Windows_NT'
%     dataFolder='c:/Users/David/Dropbox/Research/Code/data/';    
% end
dataSaveFolder=[dataFolder 'data_push_result/'];
drawArrow3 = @(x,y,varargin) quiver3( x(1),x(2),x(3),y(1)-x(1),y(2)-x(2),y(3)-x(3),varargin{:} );

load('exp/exp_push_1.mat');
%obj_poses_wID = push_predict_dis{min_idx(1),min_idx(2)};

fig1 = figure;
hold on;
visualize_objs(fig1, obj_poses_wID, [0 1 0]);
%visualize_objs(fig_table, obj_poses_wID_ap, [0 0 1]);
drawArrow3([push_command{min_idx(1),min_idx(2)}(1,1:2) 0.8],[push_command{min_idx(1),min_idx(2)}(2,1:2) 0.8],'LineWidth',2,'Color',[1 0 0]); 
visualize_objs(fig1, push_predict_dis{min_idx(1),min_idx(2)}, [0 0 1]);
camroll(-90);
axis([0.40 1.05 -0.5 0.5]);

load('exp/exp_push_2.mat');
fig2 = figure;
hold on;
visualize_objs(fig2, obj_poses_wID, [0 1 0]);
%visualize_objs(fig_table, obj_poses_wID_ap, [0 0 1]);
drawArrow3([push_command{min_idx(1),min_idx(2)}(1,1:2) 0.8],[push_command{min_idx(1),min_idx(2)}(2,1:2) 0.8],'LineWidth',2,'Color',[1 0 0]); 
visualize_objs(fig2, push_predict_dis{min_idx(1),min_idx(2)}, [0 0 1]);
camroll(-90);
axis([0.40 1.05 -0.5 0.5]);