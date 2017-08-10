function [bool_cont cont_time_step] = test_contact_through_traj(obj1_trajectory, obj1_pointcloud, obj1_pointcloud_norm, obj2_pointcloud, obj2_pointcloud_norm, obj2_pose)
    step_size = 10;
    bool_cont = false;
    cont_time_step = size(obj1_trajectory,1); % does not contact
    for time_step = step_size:step_size:size(obj1_trajectory,1)% visualize according to the pose
        cur_obj1_pose = [eGetR(obj1_trajectory(time_step,4:6)) obj1_trajectory(time_step,1:3)'; 0 0 0 1];
        
        modelpoints=cur_obj1_pose*[obj1_pointcloud'; ones(1,size(obj1_pointcloud,1))];
        normalpoints = [cur_obj1_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[obj1_pointcloud_norm'; ones(1,size(obj1_pointcloud_norm,1))];
        cur_obj1_points=modelpoints(1:3,:)';
        cur_obj1_norms=normalpoints(1:3,:)';
        
        cur_obj2_pose = [eGetR(obj2_pose(1,4:6)) obj2_pose(1,1:3)'; 0 0 0 1];
        obj2_modelPoints=cur_obj2_pose*[obj2_pointcloud'; ones(1,size(obj2_pointcloud,1))];
        obj2_normalPoints = [cur_obj2_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[obj2_pointcloud_norm'; ones(1,size(obj2_pointcloud_norm,1))];
        cur_obj2_modelpoints=obj2_modelPoints(1:3,:)';
        cur_obj2_normalpoints=obj2_normalPoints(1:3,:)';
        
        dist_th = 0.02;
        dot_th = -0.75;
        num_th = 10;
        
        if time_step == size(obj1_trajectory,1)
            nxt_obj1_pose = cur_obj1_pose; % end of execution
        else
            nxt_obj1_pose = [eGetR(obj1_trajectory(time_step+step_size,4:6)) obj1_trajectory(time_step+step_size,1:3)'; 0 0 0 1];
        end
        cur_obj1_action_pose = nxt_obj1_pose * cur_obj1_pose^-1;
        
        cont_frame = zeros(4,4);
        cont_frame(4,4) = 1;
        obj1_cent = mean(cur_obj1_points);
        cont_frame(1:3,4) = obj1_cent';
        cont_frame(1:3,1:3) = cur_obj1_action_pose(1:3,1:3);
        %[bool_cont cont_pcd_EE cont_norm_EE cont_pcd_obj cont_norm_obj] = compute_contact_points(exp_gripper_points, exp_gripper_norms, cur_obj_modelpoints{i}, cur_obj_normalpoints{i}, cont_frame, dist_th, dot_th, num_th);
        [bool_cont cont_pcd_EE cont_norm_EE cont_pcd_obj cont_norm_obj] = compute_contact_points(cur_obj1_points, cur_obj1_norms, cur_obj2_modelpoints, cur_obj2_normalpoints, cont_frame, dist_th, dot_th, num_th);
        
        if bool_cont == true
            cont_time_step = time_step;
            return;
        end
    end

end