function pose = request_obj_pose(target_obj)
    %Start
    addpath('~/catkin_ws/src/affordance_learning/simple_robot_vision/models/ObjectModeling/Utilities/');
    dataFolder='~/catkin_ws/data/';
    modelFolder=[dataFolder 'models/'];
    
    %load model with ARtagID
    %modelname = target_obj;
    modelName = ['DKCheezitBox1_wNorm'];
    ARtagID = 0;
    
    %Get new scene
    %system('timeout 10s roslaunch simple_robot_vision kinect_dump.launch');
    %pause(1);
    
    %load transform of ARtag
    load([modelFolder 'ARmodels/Transform_',modelName,'_',num2str(ARtagID),'.mat']);
    %load ARtag position
    system(['rosrun ar_track_alvar request_obj_poses _objID:=',num2str(ARtagID)]);
    pose = load([dataFolder 'obj_pose.csv']);
  
    arTagPose=[qGetR(pose(5:8)) pose(2:4)'; 0 0 0 1];
    %arTagPose=[quat2rotm(pose(5:8)) pose(2:4)'; 0 0 0 1];
    newobjPose=arTagPose*arTagH;

    %Load model and transform
    [modelPointsNnormals, modelColorNzeros]=loadCloudColor([modelFolder,modelName,'.pcd']);
    modelPoints=modelPointsNnormals(1:2:end-1,:);
    %modelpoints=bsxfun(@plus, modelPoints*newobjPose(1:3,1:3), newobjPose(1:3,4)');
    modelpoints=newobjPose*[modelPoints'; ones(1,size(modelPoints,1))];
    modelpoints=modelpoints(1:3,:)';
    
    %Load scene points
    scenePC=loadpcd('~/catkin_ws/data/Kinect_Scene_1.pcd');
    scenePC=scenePC(1:3,:)';

    %PLOT
    figure
    hold on
    plot3(scenePC(:,1),scenePC(:,2),scenePC(:,3),'.');
    %plotCoord(pose(2:4),quat2rotm(pose(5:end)),0.025);
    plot3(modelpoints(:,1),modelpoints(:,2),modelpoints(:,3),'.g')
    axis equal

    
%     model_rot = [0 0 -1;0 1 0; 1 0 0];
%     rotm = quat2rotm(obj_pose(4:7));
%     modelpoints= modelpoints*(rotm*model_rot)';
%     
%     for i=1:size(modelpoints,1)
%         modelpoints(i,1) = modelpoints(i,1) + obj_pose(1); %x
%         modelpoints(i,2) = modelpoints(i,2) + obj_pose(2); %y
%         modelpoints(i,3) = modelpoints(i,3) + obj_pose(3); %z
%     end

%     figure;
%     plot3(modelPoints(:,1),modelPoints(:,2),modelPoints(:,3),'Color',[1 0 0],'Marker','.','Linestyle','none');

    
    
    %load transform bew ARtag and the model
    
    
    
    

end