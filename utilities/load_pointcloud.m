function [pointcloud pointcloud_norm] = load_pointcloud(obj_poses_wID)
    %Start
    addpath('~/catkin_ws/src/affordance_learning/simple_robot_vision/models/ObjectModeling/Utilities/');
    dataFolder='~/catkin_ws/data/';
%     if getenv('OS')=='Windows_NT'
%         addpath('c:/Users/David/Dropbox/Research/Code/simple_robot_vision/models/ObjectModeling/Utilities/');
%         dataFolder='c:/Users/David/Dropbox/Research/Code/data/';    
%     end
    modelFolder=[dataFolder 'models/'];
    load([dataFolder 'models/ARmodels/pair_map.mat']);
    
    for i=1:size(obj_poses_wID,1)
        crnt_ar_id = obj_poses_wID(i,1);
        
        %load model with ARtagID
        for j=1:size(pair_map,1)
            if pair_map{j,2} == num2str(crnt_ar_id)
                crnt_modelName = pair_map{j,1};
            end
        end
        ARtagID = crnt_ar_id;
        
        crnt_obj_pose = [qGetR(obj_poses_wID(i,5:8)) obj_poses_wID(i,2:4)'; 0 0 0 1];
        
        %Load model and transform
        [modelPointsNnormals, modelColorNzeros]=loadCloudColor([modelFolder,crnt_modelName,'.pcd']);
        modelPoints=modelPointsNnormals(1:2:end-1,:);
        normalPoints = modelPointsNnormals(2:2:end,:);
        
        modelpoints=crnt_obj_pose*[modelPoints'; ones(1,size(modelPoints,1))];
        %normalpoints = crnt_obj_pose*[normalPoints'; ones(1,size(normalPoints,1))];
        normalpoints = [crnt_obj_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[normalPoints'; ones(1,size(normalPoints,1))];
        modelpoints=modelpoints(1:3,:)';
        normalpoints=normalpoints(1:3,:)';
        
        pointcloud{i} = modelpoints;
        pointcloud_norm{i} = normalpoints;
           
    end

end