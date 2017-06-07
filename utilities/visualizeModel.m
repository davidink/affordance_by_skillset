function visualize_model(obj_name, obj_pose)
    addpath('~/catkin_ws/src/affordance_learning/simple_robot_vision/models/ObjectModeling/Utilities/')

    file_path = '~/catkin_ws/src/affordance_learning/simple_robot_vision/models/'
    
    [modelpoints modelcolors]=loadCloudColor([file_path,obj_name,'.pcd');
    
    Ang1=0*pi/180; %Rotate around table normal
    Ang2=-90*pi/180; %Rotate around local y-axis
    Ang3=90*pi/180; %Rotate around 



    Ra=[cos(Ang1) -sin(Ang1) 0; sin(Ang1) cos(Ang1) 0; 0 0 1];
    Rb=[cos(Ang2) 0 sin(Ang2); 0 1 0; -sin(Ang2) 0 cos(Ang2)];
    Rc=[cos(Ang3) -sin(Ang3) 0; sin(Ang3) cos(Ang3) 0; 0 0 1];
    R=Ra*Rb*Rc;
    figure
    hold all
    plot3(modelpoints(:,1),modelpoints(:,2),modelpoints(:,3),'.b')
    modelpoints=modelpoints*R';
    plot3(modelpoints(:,1),modelpoints(:,2),modelpoints(:,3),'.r')

    view(90,90);
    xlabel('x');
    ylabel('x');
    zlabel('x');

    axis equal

    qw=sqrt(1+R(1,1)+R(2,2)+R(3,3))/2
    qx=(R(3,2)-R(2,3))/(4*qw)
    qy=(R(1,3)-R(3,1))/(4*qw)
    qz=(R(2,1)-R(1,2))/(4*qw)

    q=[qw qx qy qz]
end