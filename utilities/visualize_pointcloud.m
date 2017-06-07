function output = visualize_pointcloud(file_name,fig_h)
%     addpath('~/catkin_ws/src/affordance_learning/simple_robot_vision/models/ObjectModeling/Utilities/');
%     file_path = '~/catkin_ws/src/affordance_learning/simple_robot_vision/models/';
    data_file_path = '~/catkin_ws/data/';
    
    [modelpoints modelcolors]=loadCloudColor([data_file_path,file_name]);
        
    output = modelpoints;
    
    %close(figure(fig_h));    
    %fig_h = figure;
    figure(fig_h);
    scatter3(modelpoints(:,1),modelpoints(:,2),modelpoints(:,3),1,modelcolors);
    view(-90,0);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal;
    axis([-1 1 -0.3 0.3 -0.1 0.5]);

    
%     for i=1:size(modelpoints,1)
%         plot3(modelpoints(i,1),modelpoints(i,2),modelpoints(i,3),'Color',[modelcolors(i,1) modelcolors(i,2) modelcolors(i,3)],'Marker','.','Linestyle','none');   
%     end
    
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