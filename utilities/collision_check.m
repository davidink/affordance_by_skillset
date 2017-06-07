function result = collision_check(obj_name,exp_obj_pose,obs_name,obs_pose)
%     addpath('~/catkin_ws/src/affordance_learning/simple_robot_vision/models/ObjectModeling/Utilities/');
%     file_path = '~/catkin_ws/src/affordance_learning/simple_robot_vision/models/';    
%     [obj_pcl obj_color]=loadCloudColor([file_path,obj_name,'.pcd']);
%     [obs_pcl obs_color]=loadCloudColor([file_path,obs_name,'.pcd']);
%     
%     model_rot = [0 0 -1;0 1 0; 1 0 0];
%     rotm = quat2rotm(exp_obj_pose(4:7));
%     obj_pcl = obj_pcl*(rotm*model_rot)';
%     for i=1:size(obj_pcl,1)
%         obj_pcl(i,1) = obj_pcl(i,1) + exp_obj_pose(1); %x
%         obj_pcl(i,2) = obj_pcl(i,2) + exp_obj_pose(2); %y
%         obj_pcl(i,3) = obj_pcl(i,3) + exp_obj_pose(3); %z
%     end
%     
%     rotm = quat2rotm(obs_pose(4:7));
%     obs_pcl = obs_pcl*(rotm*model_rot)';
%         for i=1:size(obs_pcl,1)
%         obs_pcl(i,1) = obs_pcl(i,1) + obs_pose(1); %x
%         obs_pcl(i,2) = obs_pcl(i,2) + obs_pose(2); %y
%         obs_pcl(i,3) = obs_pcl(i,3) + obs_pose(3); %z
%     end
%     
% %     figure;
% %     hold on;
%     prj_obj_pcl = [obj_pcl(:,1) obj_pcl(:,2)];
%     prj_obs_pcl = [obs_pcl(:,1) obs_pcl(:,2)];
%     bound_obj = boundary(prj_obj_pcl);
%     bound_obs = boundary(prj_obs_pcl);
%     
%     obj_bound_pcl = [obj_pcl(bound_obj,1) obj_pcl(bound_obj,2)];
%     obs_bound_pcl = [obs_pcl(bound_obs,1) obs_pcl(bound_obs,2)];
    
    dist_cent = sqrt((exp_obj_pose(1)-obs_pose(1))^2+(exp_obj_pose(2)-obs_pose(2))^2);
    
%     obj_dim_max = -1;
%     for i=1:size(obj_bound_pcl,1)
%         obj_dim = sqrt((exp_obj_pose(1)-obj_bound_pcl(i,1))^2+(exp_obj_pose(2)-obj_bound_pcl(i,2))^2);
%         if obj_dim > obj_dim_max
%             obj_dim_max = obj_dim;
%         end
%     end
%     obs_dim_max = -1;   
%     for i=1:size(obs_bound_pcl,1)
%         obs_dim = sqrt((obs_pose(1)-obs_bound_pcl(i,1))^2+(obs_pose(2)-obs_bound_pcl(i,2))^2);
%         if obs_dim > obs_dim_max
%             obs_dim_max = obs_dim;
%         end
%     end
    obj_dim_max = 0.1527;
    obs_dim_max = 0.1218;
    
    if dist_cent < obj_dim_max + obs_dim_max
        result = true;
    else
        result = false;
    end
    
%     min_dist = min(min(pdist2(obj_bound_pcl,obs_bound_pcl)));
%     
%     if min_dist < 0.01
%         result = true;
%         figure;
%         hold on;
%         plot(obj_pcl(bound_obj,1),obj_pcl(bound_obj,2),'.r');
%         plot(obs_pcl(bound_obs,1),obs_pcl(bound_obs,2),'.b');
%         close;
%         
%     else result = false;
%     end
    
%     plot3(obj_pcl(:,1),obj_pcl(:,2),obj_pcl(:,3),'.r')
%     plot3(obs_pcl(:,1),obs_pcl(:,2),obs_pcl(:,3),'.b')
        

end