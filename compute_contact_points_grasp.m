function [bool_cont cont_pcd_EE cont_norm_EE cont_pcd_obj cont_norm_obj] = compute_contact_points_grasp(cur_gripper_points, cur_gripper_norms, cur_obj_modelpoints, cur_obj_normalpoints, cont_frame, dist_th, dot_th, num_th)
    
%     prox_list = [];
%     for i=1:size(cur_gripper_points,1)
%         for j=1:size(cur_obj_modelpoints,1)
%             if sum((cur_gripper_points(i,:)-cur_obj_modelpoints(j,:)).^2) < dist_th
%                 prox_list = [prox_list;i j];
%             end
%         end
%     end

    %convert pcd to contact frame
    [ee_points_CF ee_norms_CF] = convert_pcd_frame(cur_gripper_points, cur_gripper_norms, cont_frame^-1);
%     rev_frame = cont_frame^-1;
%     ee_pcd = rev_frame*[cur_gripper_points'; ones(1,size(cur_gripper_points,1))];
%     ee_norms = [rev_frame(1:3,1:3) [0 0 0]';0 0 0 1]*[cur_gripper_norms'; ones(1,size(cur_gripper_norms,1))];
%     ee_points_CF=ee_pcd(1:3,:)';
%     ee_norms_CF = ee_norms(1:3,:)';

    [obj_points_CF obj_norms_CF] = convert_pcd_frame(cur_obj_modelpoints, cur_obj_normalpoints, cont_frame^-1);
    
%     obj_pcd = rev_frame*[cur_obj_modelpoints'; ones(1,size(cur_obj_modelpoints,1))];
%     obj_norms = [rev_frame(1:3,1:3) [0 0 0]';0 0 0 1]*[cur_obj_normalpoints'; ones(1,size(cur_obj_normalpoints,1))];
%     obj_points_CF=obj_pcd(1:3,:)';
%     obj_norms_CF = obj_norms(1:3,:)';
%     
%     figure;
%     hold on;
%     plot3(ee_points_CF(:,1),ee_points_CF(:,2),ee_points_CF(:,3),'Color',[1 0 0],'Marker','.','Linestyle','none');
%     plot3(obj_points_CF(:,1),obj_points_CF(:,2),obj_points_CF(:,3),'Color',[0 1 0],'Marker','.','Linestyle','none');
% %     quiver3(obj_points_CF(:,1),obj_points_CF(:,2),obj_points_CF(:,3),obj_norms_CF(:,1)/100,obj_norms_CF(:,2)/100,obj_norms_CF(:,3)/100);
%     plotCoord([0 0 0]',eye(3),0.025);
%     axis([-0.3 0.3 -0.3 0.3 -0.3 0.3])
%     axis equal;
%     view(90,0);
    
    ee_cent = mean(ee_points_CF);
    obj_cent = mean(obj_points_CF);
%     if ee_cent(2) > obj_cent(2) 
%         bool_cont= false;
%         cont_pcd_EE = [];
%         cont_norm_EE = [];
%         cont_pcd_obj = [];
%         cont_norm_obj = [];
%         return;        
%     end
    
    %quiver3(ee_points_CF(:,1),ee_points_CF(:,2),ee_points_CF(:,3),ee_norms_CF(:,1)/100,ee_norms_CF(:,2)/100,ee_norms_CF(:,3)/100);

    len_roi = 0.5;
    grid_size = dist_th;
    num_side = (len_roi/grid_size);
    % Build 3D Voxel grid
    cont_idx_ee =[];
    cont_idx_obj = [];
    for i=1:num_side %z
        occu_ee = zeros(num_side,num_side);
        occu_obj = zeros(num_side,num_side);
        for j=1:num_side %y
            for k=1:num_side %x
                idx_ee = [];
                idx_obj = [];
                idx_ee = find(ee_points_CF(:,1) > (-len_roi/2 + (k-1)*grid_size) & ee_points_CF(:,1) < (-len_roi/2 + k*grid_size) ...
                         & ee_points_CF(:,2) > (-len_roi/2 + (j-1)*grid_size) & ee_points_CF(:,2) < (-len_roi/2 + j*grid_size) ...
                         & ee_points_CF(:,3) > (-len_roi/2 + (i-1)*grid_size) & ee_points_CF(:,3) < (-len_roi/2 + i*grid_size)); 
                idx_obj = find(obj_points_CF(:,1) > (-len_roi/2 + (k-1)*grid_size) & obj_points_CF(:,1) < (-len_roi/2 + k*grid_size) ...
                         & obj_points_CF(:,2) > (-len_roi/2 + (j-1)*grid_size) & obj_points_CF(:,2) < (-len_roi/2 + j*grid_size) ...
                         & obj_points_CF(:,3) > (-len_roi/2 + (i-1)*grid_size) & obj_points_CF(:,3) < (-len_roi/2 + i*grid_size)); 
                
%                 if size(idx_ee,1) ~= 0 % if some points are within the grid
%                     %plot3(ee_points_CF(idx_ee,1),ee_points_CF(idx_ee,2),ee_points_CF(idx_ee,3),'Color',[0 0 1],'Marker','o','Linestyle','none');
%                     occu_ee(k,j) = 1;%
%                     occu_idx_ee{k+(j-1)*num_side+(i-1)*num_side*num_side} = idx_ee;
%                 end
%                 if size(idx_obj,1) ~= 0 % if some points are within the grid
%                     %plot3(ee_points_CF(idx,1),ee_points_CF(idx,2),ee_points_CF(idx,3),'Color',[0 0 1],'Marker','o','Linestyle','none');
%                     occu_obj(k,j) = 1;%
%                     occu_idx_obj{k+(j-1)*num_side+(i-1)*num_side*num_side} = idx_obj;
%                 end
                if size(idx_ee,1) ~= 0 && size(idx_obj,1) ~= 0
%                     plot3(ee_points_CF(idx_ee,1),ee_points_CF(idx_ee,2),ee_points_CF(idx_ee,3),'Color',[0 0 1],'Marker','o','Linestyle','none');
%                     plot3(obj_points_CF(idx_obj,1),obj_points_CF(idx_obj,2),obj_points_CF(idx_obj,3),'Color',[0 0 1],'Marker','o','Linestyle','none');  

                    idx_norm = [];
                    [idx_norm(:,1) idx_norm(:,2)] = find(ee_norms_CF(idx_ee,2)*obj_norms_CF(idx_obj,2)' < dot_th);
                    cont_idx_ee = [cont_idx_ee;idx_ee(unique(idx_norm(:,1)))];
                    cont_idx_obj = [cont_idx_obj;idx_obj(unique(idx_norm(:,2)))]; 
                    
%                     cont_idx_ee = [cont_idx_ee;idx_ee];
%                     cont_idx_obj = [cont_idx_obj;idx_obj];

%                     plot3(ee_points_CF(cont_idx_ee,1),ee_points_CF(cont_idx_ee,2),ee_points_CF(cont_idx_ee,3),'Color',[0 0 1],'Marker','o','Linestyle','none');
%                     plot3(obj_points_CF(cont_idx_obj,1),obj_points_CF(cont_idx_obj,2),obj_points_CF(cont_idx_obj,3),'Color',[0 1 1],'Marker','o','Linestyle','none');                    
                end
            end
        end
%         voxel_occu_ee(:,:,i) = occu_ee;
%         voxel_occu_obj(:,:,i) = occu_obj;
    end
    
    % search for overlapping voxel
    %if size(cont_idx_ee,1) ~=0 && size(cont_idx_obj,1) > 5
    if size(cont_idx_ee,1) >5 && size(cont_idx_obj,1) > 5
        bool_cont = true;
%         plot3(ee_points_CF(cont_idx_ee,1),ee_points_CF(cont_idx_ee,2),ee_points_CF(cont_idx_ee,3),'Color',[0 0 1],'Marker','o','Linestyle','none');
%         plot3(obj_points_CF(cont_idx_obj,1),obj_points_CF(cont_idx_obj,2),obj_points_CF(cont_idx_obj,3),'Color',[0 1 1],'Marker','o','Linestyle','none');
        % revert contact points into global frame;
        [cont_pcd_EE cont_norm_EE] = convert_pcd_frame(ee_points_CF(cont_idx_ee,:),ee_norms_CF(cont_idx_ee,:),cont_frame);
        [cont_pcd_obj cont_norm_obj] = convert_pcd_frame(obj_points_CF(cont_idx_obj,:),obj_norms_CF(cont_idx_obj,:),cont_frame);
    else
        bool_cont = false;
        cont_pcd_EE = [];
        cont_norm_EE = [];
        cont_pcd_obj = [];
        cont_norm_obj = [];
    end
%     idx_ovlp = find(voxel_occu_ee+voxel_occu_obj==2);
%     if idx_ovlp ~=0
%         bool_cont = true;
%         cont_idx_ee =[];
%         cont_idx_obj = [];
%         for i=1:size(idx_ovlp,1)
%             cont_idx_ee = [cont_idx_ee; occu_idx_ee{idx_ovlp(i)}];
%             cont_idx_obj = [cont_idx_obj; occu_idx_obj{idx_ovlp(i)}];
%         end
% %         plot3(ee_points_CF(cont_idx_ee,1),ee_points_CF(cont_idx_ee,2),ee_points_CF(cont_idx_ee,3),'Color',[0 0 1],'Marker','o','Linestyle','none');
% %         plot3(obj_points_CF(cont_idx_obj,1),obj_points_CF(cont_idx_obj,2),obj_points_CF(cont_idx_obj,3),'Color',[0 1 1],'Marker','o','Linestyle','none');
% %         quiver3(ee_points_CF(cont_idx_ee,1),ee_points_CF(cont_idx_ee,2),ee_points_CF(cont_idx_ee,3),ee_norms_CF(cont_idx_ee,1)/100,ee_norms_CF(cont_idx_ee,2)/100,ee_norms_CF(cont_idx_ee,3)/100);
% %         quiver3(obj_points_CF(cont_idx_obj,1),obj_points_CF(cont_idx_obj,2),obj_points_CF(cont_idx_obj,3),obj_norms_CF(cont_idx_obj,1)/100,obj_norms_CF(cont_idx_obj,2)/100,obj_norms_CF(cont_idx_obj,3)/100);
%         [cont_idx_norm(:,1) cont_idx_norm(:,2) cont_idx_norm(:,3)] = find(ee_norms_CF(cont_idx_ee,:)*obj_norms_CF(cont_idx_obj,:)' <dot_th);
%         cont_pcd_EE = [];
%         cont_pcd_obj = [];
%     else
%         bool_cont = false;
%         cont_pcd_EE = [];
%         cont_pcd_obj = [];
%     end
%     
%     figure;
%     plot3(cont_pcd_EE(:,1),cont_pcd_EE(:,2),cont_pcd_EE(:,3),'.r')
%     hold on
%     plot3(cont_pcd_obj(:,1),cont_pcd_obj(:,2),cont_pcd_obj(:,3),'.g')

end