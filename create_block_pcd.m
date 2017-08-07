function [modelpoints normpoints] = create_block_pcd(x_size,y_size,z_size,granual,obj_cent)
%     granual = 0.01;
%     x_size = 0.05;
%     y_size = 0.05;
%     z_size = 0.2;
%     obj_cent = [0 0 0];
    %num_point = (x_size/granual)*(y_size/granual)*(z_size/granual);
    %modelpoints = zeros(num_point,3);
    modelpoints = [];
    normpoints = [];

    % build 6 surfaces
    for i=1:x_size/granual+1
        for j=1:y_size/granual+1
            x_val = obj_cent(1,1) - x_size/2; 
            y_val = obj_cent(1,2) - y_size/2;
            z_val = z_size/2 + obj_cent(1,3);
            surf_points(round((i-1)*(y_size/granual+1)+j),:)=[x_val+(i-1)*granual y_val+(j-1)*granual z_val];
            norm_points(round((i-1)*(y_size/granual+1)+j),:)=[0 0 1];
        end
    end
    modelpoints = [modelpoints;surf_points];
    normpoints = [normpoints;norm_points];
    for i=1:x_size/granual+1
        for j=1:y_size/granual+1
            x_val = obj_cent(1,1) - x_size/2; 
            y_val = obj_cent(1,2) - y_size/2;
            z_val = -z_size/2 + obj_cent(1,3);
            surf_points(round((i-1)*(y_size/granual+1)+j),:)=[x_val+(i-1)*granual y_val+(j-1)*granual z_val];
            norm_points(round((i-1)*(y_size/granual+1)+j),:)=[0 0 -1];
        end
    end
    modelpoints = [modelpoints;surf_points];
    normpoints = [normpoints;norm_points];

    for i=1:x_size/granual+1
        for j=1:z_size/granual+1
            x_val = obj_cent(1,1) - x_size/2; 
            y_val = obj_cent(1,2) - y_size/2;
            z_val = obj_cent(1,3) - z_size/2;
            surf_points(round((i-1)*(z_size/granual+1)+j),:)=[x_val+(i-1)*granual y_val z_val+(j-1)*granual];
            norm_points(round((i-1)*(z_size/granual+1)+j),:)=[0 -1 0];
        end
    end
    modelpoints = [modelpoints;surf_points];
    normpoints = [normpoints;norm_points];
    for i=1:x_size/granual+1
        for j=1:z_size/granual+1
            x_val = obj_cent(1,1) - x_size/2; 
            y_val = obj_cent(1,2) + y_size/2;
            z_val = obj_cent(1,3) - z_size/2;
            surf_points(round((i-1)*(z_size/granual+1)+j),:)=[x_val+(i-1)*granual y_val z_val+(j-1)*granual];
            norm_points(round((i-1)*(z_size/granual+1)+j),:)=[0 1 0];
        end
    end
    modelpoints = [modelpoints;surf_points];
    normpoints = [normpoints;norm_points];

    for i=1:y_size/granual+1
        for j=1:z_size/granual+1
            x_val = obj_cent(1,1) - x_size/2; 
            y_val = obj_cent(1,2) - y_size/2;
            z_val = obj_cent(1,3) - z_size/2;
            surf_points(round((i-1)*(z_size/granual+1)+j),:)=[x_val y_val+(i-1)*granual z_val+(j-1)*granual];
            norm_points(round((i-1)*(z_size/granual+1)+j),:)=[-1 0 0];
        end
    end
    modelpoints = [modelpoints;surf_points];
    normpoints = [normpoints;norm_points];
    for i=1:y_size/granual+1
        for j=1:z_size/granual+1
            x_val = obj_cent(1,1) + x_size/2; 
            y_val = obj_cent(1,2) - y_size/2;
            z_val = obj_cent(1,3) - z_size/2;
            surf_points(round((i-1)*(z_size/granual+1)+j),:)=[x_val y_val+(i-1)*granual z_val+(j-1)*granual];
            norm_points(round((i-1)*(z_size/granual+1)+j),:)=[1 0 0];
        end
    end
    modelpoints = [modelpoints;surf_points];
    normpoints = [normpoints;norm_points];

%     figure;
%     plot3(modelpoints(:,1),modelpoints(:,2),modelpoints(:,3),'Color',[1 0 0],'Marker','.','Linestyle','none');
%     hold on;
%     %plotCoord(crnt_obj_pose(1:3,4)',crnt_obj_pose(1:3,1:3),0.025);
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     axis equal

end

