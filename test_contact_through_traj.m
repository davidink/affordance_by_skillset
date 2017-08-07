function bool_cont = test_contact_through_traj(obj1_pointcloud, obj1_pointcloud_norm, obj1_trajectory, obj2_pointcloud, obj2_pointcloud_norm, obj2_pose)
    step_size = 10;
    bool_cont = false;
    for time_step = step_size:step_size:size(obj_trajectory_1,1)% visualize according to the pose
        cur_obj1_pose = [eGetR(obj1_trajectory(time_step,4:6)) obj1_trajectory(time_step,1:3)'; 0 0 0 1];
        
        modelpoints=cur_obj1_pose*[obj1_pointcloud'; ones(1,size(obj1_pointcloud,1))];
        normalpoints = [cur_obj1_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[obj1_pointcloud_norm'; ones(1,size(obj1_pointcloud_norm,1))];
        cur_obj1_points=modelpoints(1:3,:)';
        cur_obj1_norms=normalpoints(1:3,:)';
        
        obj2_pose = [];
                        obj_pose_init{i} = [eGetR(obj_trajectory{i}(time_step,4:6)) obj_trajectory{i}(time_step,1:3)'; 0 0 0 1];
                obj_modelPoints=obj_pose_init{i}*[obj_modelpoints{i}'; ones(1,size(obj_modelpoints{i},1))];
                obj_normalPoints = [obj_pose_init{i}(1:3,1:3) [0 0 0]';0 0 0 1]*[obj_normpoints{i}'; ones(1,size(obj_normpoints{i},1))];
                cur_obj_modelpoints{i}=obj_modelPoints(1:3,:)';
                cur_obj_normalpoints{i}=obj_normalPoints(1:3,:)';
    end

end