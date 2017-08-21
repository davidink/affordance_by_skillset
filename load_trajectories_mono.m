function result = load_trajectories_mono(dataFolder, scene_num, opt_save_kernel_feat, opt_save_feat_n_result)

    addpath('utilities');
    % dataFolder = 'data/trajectories/';
    %scene_num = 4;
    trajectory = load([dataFolder 'trajectory' num2str(scene_num) '.csv']);
    fig_gt = figure;
    step_size = 10;
    disp_trajectory(fig_gt, [dataFolder 'trajectory' num2str(scene_num)], step_size);
    
    % Data loading
    Desired_EE_pose_ref_base = trajectory(:,1:6); % relative to robot_base joint
    %convert into global frame 
    Desired_EE_pose = Desired_EE_pose_ref_base + repmat([0 0 0.25 0 1.57 -1.57],size(Desired_EE_pose_ref_base,1),1);
    Executed_EE_pose = trajectory(:,7:12);
    num_obj = trajectory(:,13);

    for i=1:num_obj
        obj_trajectory{i} = trajectory(:,14+(i-1)*6:14+(i-1)*6+5);
        obj_param{i} = trajectory(:,14+6*num_obj+(i-1)*7:14+6*num_obj+(i-1)*7+6);
    end

    fig_traj = figure;
    hold on;
    plot3(Desired_EE_pose(:,1),Desired_EE_pose(:,2),Desired_EE_pose(:,3),'.r');
    plot3(Executed_EE_pose(:,1),Executed_EE_pose(:,2),Executed_EE_pose(:,3),'.b');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal;

    % show trajectory of end-effector 

    % load the gripper
    load('PR2_gripper.mat');

    gripper_points= modelpoints;
    gripper_norms = normalpoints;
    for i=1:num_obj
        granual = 0.005;
        obj_dim(i,:) = ceil(obj_param{i}(1,4:6)*100)/100;
        %obj_dim(i,:) = floor(obj_param{i}(1,4:6)*100)/100;
    %     x_size = ceil(obj_param{i}(1,4)*100)/100;
    %     y_size = ceil(obj_param{i}(1,5)*100)/100;
    %     z_size = ceil(obj_param{i}(1,6)*100)/100;
        obj_cent = [0.0 0.0 0];
        [obj_modelpoints{i} obj_normpoints{i}] = create_block_pcd(obj_dim(i,1),obj_dim(i,2),obj_dim(i,3),granual,obj_cent);
    end

    % Draw each time step of trajectories
    step_size = 10;
    features = [];
    results = [];
    
    % create gripper pcd
    gripper_pose = [eGetR(Executed_EE_pose(1,4:6)) Executed_EE_pose(1,1:3)'; 0 0 0 1];
    end_gripper_pose = [eGetR(Executed_EE_pose(end,4:6)) Executed_EE_pose(end,1:3)'; 0 0 0 1];
    action_pose = end_gripper_pose * gripper_pose^-1;
    modelpoints=gripper_pose*[gripper_points'; ones(1,size(gripper_points,1))];
    normalpoints = [gripper_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[gripper_norms'; ones(1,size(gripper_norms,1))];
    cur_gripper_points=modelpoints(1:3,:)';
    cur_gripper_norms=normalpoints(1:3,:)';
    % create pcd of objects
    for i=1:num_obj
        obj_pose = [eGetR(obj_trajectory{i}(1,4:6)) obj_trajectory{i}(1,1:3)'; 0 0 0 1];
        obj_modelPoints=obj_pose*[obj_modelpoints{i}'; ones(1,size(obj_modelpoints{i},1))];
        obj_normalPoints = [obj_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[obj_normpoints{i}'; ones(1,size(obj_normpoints{i},1))];
        cur_obj_modelpoints{i}=obj_modelPoints(1:3,:)';
        cur_obj_normalpoints{i}=obj_normalPoints(1:3,:)';
    end
    all_obj_modelpoints= [];
    for i=1:num_obj
        all_obj_modelpoints = [all_obj_modelpoints;cur_obj_modelpoints{i}];
    end
    
    %compute features
    cur_obj1_pose = [eGetR(obj_trajectory{1}(1,4:6)) obj_trajectory{1}(1,1:3)'; 0 0 0 1];
    cur_obj1_pose_af = cur_obj1_pose * gripper_pose^-1;
    cur_obj2_pose = [eGetR(obj_trajectory{2}(1,4:6)) obj_trajectory{2}(1,1:3)'; 0 0 0 1];
    cur_obj2_pose_af = cur_obj2_pose * gripper_pose^-1;
    global_shape_features = compute_grid_feature(cur_gripper_points, all_obj_modelpoints, cur_obj1_pose_af, 0.25, 0.05);
    
    relative_pose1 = cur_obj1_pose*gripper_pose^-1;
    relative_pose1_features = [relative_pose1(1:3,4)' RGete(relative_pose1(1:3,1:3))'];
    %cur_obj2_pose = [eGetR(obj_trajectory{2}(1,4:6)) obj_trajectory{2}(1,1:3)'; 0 0 0 1];
    relative_pose2 = cur_obj2_pose*gripper_pose^-1;
    relative_pose2_features = [relative_pose2(1:3,4)' RGete(relative_pose2(1:3,1:3))'];
    
    action_pose_features = [action_pose(1:3,4)' RGete(action_pose(1:3,1:3))'];
    
    features = [global_shape_features relative_pose1_features relative_pose2_features action_pose_features];
    
    obj1_end_pose = [eGetR(obj_trajectory{1}(end,4:6)) obj_trajectory{1}(end,1:3)'; 0 0 0 1];
    obj1_end_pose_af = obj1_end_pose * gripper_pose^-1;    
    obj1_diff = obj1_end_pose_af * cur_obj1_pose_af^-1;
    result1 = [obj1_diff(1:3,4)' RGete(obj1_diff(1:3,1:3))'];
    
    obj2_end_pose = [eGetR(obj_trajectory{2}(end,4:6)) obj_trajectory{2}(end,1:3)'; 0 0 0 1];
    obj2_end_pose_af = obj2_end_pose * gripper_pose^-1;
    obj2_diff = obj2_end_pose_af * cur_obj2_pose_af^-1;
    result2 = [obj2_diff(1:3,4)' RGete(obj2_diff(1:3,1:3))'];
    
    if opt_save_feat_n_result
       save([dataFolder 'feat_n_result1' num2str(scene_num) '.mat'],'features','result1');
       save([dataFolder 'feat_n_result2' num2str(scene_num) '.mat'],'features','result2');
       disp('done!');
    end
    
end

