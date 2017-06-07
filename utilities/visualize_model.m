function output = visualize_model(obj_name, obj_pose, color)
    addpath('~/catkin_ws/src/affordance_learning/simple_robot_vision/models/ObjectModeling/Utilities/');

    file_path = '~/catkin_ws/src/affordance_learning/simple_robot_vision/models/';
    
    [modelpoints modelcolors]=loadCloudColor([file_path,obj_name,'.pcd']);
    
    %obj_theta = 1;
    
    model_rot = [0 0 -1;0 1 0; 1 0 0];
        
%     Ang1=0*pi/180; %Rotate around table normal
%     Ang2=-90*pi/180; %Rotate around local y-axis
%     Ang3=90*pi/180; %Rotate around 
%     
%     Ra=[cos(Ang1) -sin(Ang1) 0; sin(Ang1) cos(Ang1) 0; 0 0 1];
%     Rb=[cos(Ang2) 0 sin(Ang2); 0 1 0; -sin(Ang2) 0 cos(Ang2)];
%     Rc=[cos(Ang3) -sin(Ang3) 0; sin(Ang3) cos(Ang3) 0; 0 0 1];
% 
     rotm = quat2rotm(obj_pose(4:7));
%     R=Ra*Rb*Rc*rotm;
     modelpoints= modelpoints*(rotm*model_rot)';
    
    for i=1:size(modelpoints,1)
        modelpoints(i,1) = modelpoints(i,1) + obj_pose(1); %x
        modelpoints(i,2) = modelpoints(i,2) + obj_pose(2); %y
        modelpoints(i,3) = modelpoints(i,3) + obj_pose(3); %z
    end
    
    output = modelpoints;
    
    plot3(modelpoints(:,1),modelpoints(:,2),modelpoints(:,3),'Color',color,'Marker','.','Linestyle','none');
%     modelpoints=modelpoints*R';
%     plot3(modelpoints(:,1),modelpoints(:,2),modelpoints(:,3),'.r')

%     view(90,90);
%     xlabel('x');
%     ylabel('x');
%     zlabel('x');
% 
%     axis equal
% 
%     qw=sqrt(1+R(1,1)+R(2,2)+R(3,3))/2
%     qx=(R(3,2)-R(2,3))/(4*qw)
%     qy=(R(1,3)-R(3,1))/(4*qw)
%     qz=(R(2,1)-R(1,2))/(4*qw)
% 
%     q=[qw qx qy qz]
end