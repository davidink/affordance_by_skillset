function [got_obj obj_pose] = extract_obj_pose()
    addpath('~/catkin_ws/src/affordance_learning/simple_robot_vision/models/ObjectModeling/Utilities/');
    file_path = '~/catkin_ws/src/affordance_learning/simple_robot_vision/models/';
   
    obs_name = 'Battery1';
    [modelpoints modelcolors]=loadCloudColor([file_path,obs_name,'.pcd']);

    dataFolder='~/catkin_ws/data/';
    obj_obs_pose = csvread([dataFolder,'start_obj_obs_pose.txt']);
    obj_start_pose = obj_obs_pose(1,:);
    obs_pose = obj_obs_pose(2,:);
        
    model_rot = [0 0 -1;0 1 0; 1 0 0];
    rotm = quat2rotm(obs_pose(4:7));
    modelpoints= modelpoints*(rotm*model_rot)';
    
    for i=1:size(modelpoints,1)
        modelpoints(i,1) = modelpoints(i,1) + obs_pose(1); %x
        modelpoints(i,2) = modelpoints(i,2) + obs_pose(2); %y
        modelpoints(i,3) = modelpoints(i,3) + obs_pose(3); %z
    end
    
    obs_points = modelpoints;
    
    [modelpoints modelcolors]=loadCloudColor([dataFolder,'current_tabletop_pointcloud.pcd']);
    tabletop_points = modelpoints;
    
%     figure;
%     hold on;
%     plot3(obs_points(:,1),obs_points(:,2),obs_points(:,3),'Color',[1 0 0],'Marker','.','Linestyle','none');
%     plot3(tabletop_points(:,1),tabletop_points(:,2),tabletop_points(:,3),'Color',[0 1 0],'Marker','.','Linestyle','none');
    
    table_prj_points = [tabletop_points(:,1) tabletop_points(:,2)];
    obs_prj_points = [obs_points(:,1) obs_points(:,2)];
    
    obs_bound = boundary(obs_prj_points);
    obs_bound_pt = obs_prj_points(obs_bound,:);
    
%     figure;
%     hold on;
%     plot(obs_bound_pt(:,1),obs_bound_pt(:,2),'.b');
%     plot(table_prj_points(:,1),table_prj_points(:,2),'.g');
    in = inpolygon(table_prj_points(:,1),table_prj_points(:,2),obs_bound_pt(:,1),obs_bound_pt(:,2));
    
    obj_cloud_size = 0;
    obj_points=[];
    for i=1:size(tabletop_points,1)
        if in(i) ==0
            obj_cloud_size = obj_cloud_size +1;
            obj_points(obj_cloud_size,:) = tabletop_points(i,:);
        end
    end
%     figure;
%     plot3(obj_points(:,1),obj_points(:,2),obj_points(:,3),'Color',[0 0 1],'Marker','.','Linestyle','none');

    if obj_cloud_size < 100
        got_obj = false;
        obj_pose = [0 0 0 0 0 0 0];
    else
        obj_cent = mean(obj_points);
        xy_cov = cov(obj_points(:,1),obj_points(:,2));
        [Vec Ev] = eig(xy_cov);
        a= xy_cov(1,1); b= xy_cov(1,2); c= xy_cov(2,1); d= xy_cov(2,2);
        D = a*d-b*c;
        T= a+d;
        L1 = T/2 + sqrt(T^2/4-D);
        L2 = T/2 - sqrt(T^2/4-D);
        e1(1) = L1 -d;
        e1(2) = c;
        e2(1) = L2 -d;
        e2(2) = c;
        e1 = e1/norm(e1);
        e2 = e2/norm(e2);

        if e1(1) <0 
            e1 = -e1;
        end
        if e2(2) < 0
            e2 = -e2;
        end

        rotm = [0 -e2(1) -e1(1); 0 -e2(2) -e1(2); -1 0 0];
        obj_quat = rotm2quat(rotm);

        obj_pose = [obj_cent obj_quat];
        got_obj = true;
    end
end