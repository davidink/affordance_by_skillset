function [features] = compute_surface_feature(pcd_ee, norm_ee, pcd_obj, norm_obj, cont_frame)
    
    %convert inputs into ref frame
    [pcd_ee_ref norm_ee_ref] = convert_pcd_frame(pcd_ee,norm_ee,cont_frame^-1);
    [pcd_obj_ref norm_obj_ref] = convert_pcd_frame(pcd_obj,norm_obj,cont_frame^-1);
    
%     figure;
%     hold on;
%     plot3(pcd_ee_ref(:,1),pcd_ee_ref(:,2),pcd_ee_ref(:,3),'Color',[1 0 0],'Marker','.','Linestyle','none');
%     plot3(pcd_obj_ref(:,1),pcd_obj_ref(:,2),pcd_obj_ref(:,3),'Color',[0 1 0],'Marker','.','Linestyle','none');
% %     quiver3(obj_points_CF(:,1),obj_points_CF(:,2),obj_points_CF(:,3),obj_norms_CF(:,1)/100,obj_norms_CF(:,2)/100,obj_norms_CF(:,3)/100);
%     plotCoord([0 0 0]',eye(3),0.025);
%     view(90,0);
%     axis equal;
%     
    % Assume all 6 axes are independent
    ee_cloud = [pcd_ee_ref norm_ee_ref];
    obj_cloud = [pcd_obj_ref norm_obj_ref];
    
%     ee_surface = [mean(pcd_ee_ref) var(pcd_ee_ref) mean(norm_ee_ref) var(norm_ee_ref)];
%     obj_surface = [mean(pcd_obj_ref) var(pcd_obj_ref) mean(norm_obj_ref) var(norm_obj_ref)];
    
    ee_cov_triu =[];
    ee_cov = cov(ee_cloud);
    for i=1:size(ee_cloud,2)
        ee_cov_triu = [ee_cov_triu ee_cov(i,1:size(ee_cloud,2)-i+1)];
    end    
    ee_surface = [mean(ee_cloud) ee_cov_triu];
    
    obj_cov_triu =[];
    obj_cov = cov(obj_cloud);
    for i=1:size(obj_cloud,2)
        obj_cov_triu = [obj_cov_triu obj_cov(i,1:size(obj_cloud,2)-i+1)];
    end    
    obj_surface = [mean(obj_cloud) obj_cov_triu];
    
    
    features = [ee_surface obj_surface];

%     len_roi = roi_size;
%     grid_size = granual;
%     num_side = (len_roi/grid_size);
%     % Build 3D Voxel grid
%     cont_idx_ee =[];
%     cont_idx_obj = [];
%     for i=1:num_side %z
%         occu_2d = zeros(num_side,num_side);        
%         for j=1:num_side %y
%             for k=1:num_side %x
%                 idx_ee = [];
%                 idx_obj = [];
%                 idx_ee = find(pcd_ee_ref(:,1) > (-len_roi/2 + (k-1)*grid_size) & pcd_ee_ref(:,1) < (-len_roi/2 + k*grid_size) ...
%                          & pcd_ee_ref(:,2) > (-len_roi/2 + (j-1)*grid_size) & pcd_ee_ref(:,2) < (-len_roi/2 + j*grid_size) ...
%                          & pcd_ee_ref(:,3) > (-len_roi/2 + (i-1)*grid_size) & pcd_ee_ref(:,3) < (-len_roi/2 + i*grid_size)); 
%                 idx_obj = find(pcd_obj_ref(:,1) > (-len_roi/2 + (k-1)*grid_size) & pcd_obj_ref(:,1) < (-len_roi/2 + k*grid_size) ...
%                          & pcd_obj_ref(:,2) > (-len_roi/2 + (j-1)*grid_size) & pcd_obj_ref(:,2) < (-len_roi/2 + j*grid_size) ...
%                          & pcd_obj_ref(:,3) > (-len_roi/2 + (i-1)*grid_size) & pcd_obj_ref(:,3) < (-len_roi/2 + i*grid_size)); 
%                 if size(idx_ee,1) ~= 0 && size(idx_obj,1) ~= 0
%                     if size(idx_ee,1) >= size(idx_obj)
%                         occu_2d(k,j) = 1; % End-effector
%                     else
%                         occud_2d(k,j) = 2; % obj
%                     end
%                 end                
%                 if size(idx_ee,1) ~= 0 % if some points are within the grid
%                     occu_2d(k,j) = 1; %                    
%                 end
%                 if size(idx_obj,1) ~= 0 % if some points are within the grid
%                     occu_2d(k,j) = 2;%
%                 end                
%             end
%         end
%         voxel_occu(:,:,i) = occu_2d;
%     end
%     features = voxel_occu(:)';
end

