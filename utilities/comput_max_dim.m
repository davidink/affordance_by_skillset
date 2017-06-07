function output = comput_max_dim(obj_name)
    addpath('~/catkin_ws/src/affordance_learning/simple_robot_vision/models/ObjectModeling/Utilities/');
    file_path = '~/catkin_ws/src/affordance_learning/simple_robot_vision/models/';
    
    [modelpoints modelcolors]=loadCloudColor([file_path,obj_name,'.pcd']);
    model_rot = [0 0 -1;0 1 0; 1 0 0];
    %rotm = quat2rotm(obj_pose(4:7));
    modelpoints= modelpoints*(model_rot)';
    
%     plot3(modelpoints(:,1),modelpoints(:,2),modelpoints(:,3),'.r');
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
    xypoints = [modelpoints(:,2) modelpoints(:,3)];
    
    r= sqrt(xypoints(:,1).^2 + xypoints(:,2).^2);
    output = max(r);  

end