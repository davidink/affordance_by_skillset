function obj_pose = convert_ar_to_obj_pose(ar_id, ar_pose)
    
    %Start
    addpath('~/catkin_hydro_ws/src/simple_robot_vision/models/ObjectModeling/Utilities/');
    dataFolder='~/catkin_hydro_ws/data/';
    modelFolder=[dataFolder 'models/'];
    
    %load model with ARtagID
    load([dataFolder 'models/ARmodels/pair_map.mat']);
    for i=1:size(pair_map,1)
        if pair_map{i,2} == num2str(ar_id)
            modelName = pair_map{i,1};
        end
    end
    %modelname = target_obj;
    %modelName = ['DKCheezitBox1_wNorm'];
    ARtagID = ar_id;
    
    load([modelFolder 'ARmodels/Transform_',modelName,'_',num2str(ARtagID),'.mat']);     
    
    arTagPose=[qGetR(ar_pose(4:7)) ar_pose(1:3)'; 0 0 0 1];
    obj_poseR=arTagPose*arTagH;

    obj_pose = [obj_poseR(1:3,4)' qGetQ(obj_poseR(1:3,1:3))'];
end