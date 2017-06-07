addpath('~/catkin_ws/src/affordance_learning/simple_robot_vision/models/ObjectModeling/Utilities/');
file_path = '~/catkin_ws/src/affordance_learning/simple_robot_vision/models/';    
[obj_pcl obj_color]=loadCloudColor([file_path,'Book1.pcd']);
[obs_pcl obs_color]=loadCloudColor([file_path,'Battery1.pcd']);

model_rot = [0 0 -1;0 1 0; 1 0 0];
obj_pcl = obj_pcl*(model_rot);;
obs_pcl = obs_pcl*(model_rot);

%     figure;
%     hold on;
prj_obj_pcl = [obj_pcl(:,1) obj_pcl(:,2)];
prj_obs_pcl = [obs_pcl(:,1) obs_pcl(:,2)];
bound_obj = boundary(prj_obj_pcl);
bound_obs = boundary(prj_obs_pcl);

obj_bound_pcl = [obj_pcl(bound_obj,1) obj_pcl(bound_obj,2)];
obs_bound_pcl = [obs_pcl(bound_obs,1) obs_pcl(bound_obs,2)];

figure;
hold on;
plot(obj_pcl(bound_obj,1),obj_pcl(bound_obj,2),'.r');
plot(obs_pcl(bound_obs,1),obs_pcl(bound_obs,2),'.b');

