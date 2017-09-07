function test_trajectory_prediction_mono(dataFolder,scene_num)
    
    addpath('utilities');
    % display ground truth trajectory
    step_size = 10;
    fig_gt = figure;
    disp_trajectory(fig_gt,[dataFolder 'test_trajectory' num2str(scene_num)],step_size);
    trajectory = load([dataFolder 'test_trajectory' num2str(scene_num) '.csv']);

    % Data loading
    Desired_EE_pose_ref_base = trajectory(:,1:6); % relative to robot_base joint
    %convert into global frame 
    Desired_EE_pose = Desired_EE_pose_ref_base + repmat([0 0 0.25 0 1.57 -1.57],size(Desired_EE_pose_ref_base,1),1);
    Executed_EE_pose = trajectory(:,7:12);
    num_obj = trajectory(1,13);

    for i=1:num_obj
        obj_trajectory{i} = trajectory(:,14+(i-1)*6:14+(i-1)*6+5);
        obj_param{i} = trajectory(:,14+6*num_obj+(i-1)*7:14+6*num_obj+(i-1)*7+6);
    end

    fig_exp = figure;
    hold on;
    plot3(Desired_EE_pose(:,1),Desired_EE_pose(:,2),Desired_EE_pose(:,3),'.b');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal;
    axis([-0.25 0.25 -0.05 0.7 -0.2 0.2]);

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
%     step_size = 20;
    features = [];
    results = [];
    
    %compute features
    gripper_pose = [eGetR(Desired_EE_pose(1,4:6)) Desired_EE_pose(1,1:3)'; 0 0 0 1];
    end_gripper_pose = [eGetR(Desired_EE_pose(end,4:6)) Desired_EE_pose(end,1:3)'; 0 0 0 1];
    modelpoints=gripper_pose*[gripper_points'; ones(1,size(gripper_points,1))];
    normalpoints = [gripper_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[gripper_norms'; ones(1,size(gripper_norms,1))];
    cur_gripper_points=modelpoints(1:3,:)';
    cur_gripper_norms=normalpoints(1:3,:)';
    
    action_pose = end_gripper_pose * gripper_pose^-1;
    
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
    
    cur_obj1_pose = [eGetR(obj_trajectory{1}(1,4:6)) obj_trajectory{1}(1,1:3)'; 0 0 0 1];
    cur_obj1_pose_af = cur_obj1_pose * gripper_pose^-1;
    global_shape_features = compute_grid_feature(cur_gripper_points, all_obj_modelpoints, cur_obj1_pose_af, 0.25, 0.05);
    
    relative_pose1 = cur_obj1_pose*gripper_pose^-1;
    relative_pose1_features = [relative_pose1(1:3,4)' RGete(relative_pose1(1:3,1:3))'];
    cur_obj2_pose = [eGetR(obj_trajectory{2}(1,4:6)) obj_trajectory{2}(1,1:3)'; 0 0 0 1];
    cur_obj2_pose_af = cur_obj2_pose * gripper_pose^-1;
    relative_pose2 = cur_obj2_pose*gripper_pose^-1;
    relative_pose2_features = [relative_pose2(1:3,4)' RGete(relative_pose2(1:3,1:3))'];
    
    action_pose_features = [action_pose(1:3,4)' RGete(action_pose(1:3,1:3))'];
    
    features = [global_shape_features relative_pose1_features relative_pose2_features action_pose_features];
    
    obj1_end_pose = [eGetR(obj_trajectory{1}(end,4:6)) obj_trajectory{1}(end,1:3)'; 0 0 0 1];
    obj1_end_pose_af = gripper_pose*obj1_end_pose;
    obj1_diff = obj1_end_pose_af * cur_obj1_pose^-1;
    result1 = [obj1_diff(1:3,4)' RGete(obj1_diff(1:3,1:3))'];
    
    obj2_end_pose = [eGetR(obj_trajectory{2}(end,4:6)) obj_trajectory{2}(end,1:3)'; 0 0 0 1];
    obj2_end_pose_af = gripper_pose*obj2_end_pose;
    obj2_diff = obj2_end_pose_af * cur_obj2_pose^-1;
    result2 = [obj2_diff(1:3,4)' RGete(obj2_diff(1:3,1:3))'];
    
    %load([dataFolder 'GP_models_feature_ard.mat']);
    modelFolder = 'data/models/';
    load([modelFolder 'GP_models_obj1.mat']);
    gprMdl_obj1 = gprMdl;
    load([modelFolder 'GP_models_obj2.mat']);
    gprMdl_obj2 = gprMdl;
    
    for j=1:6
        [pred_obj1(1,j) pred_obj1_sd(1,j)] = predict(gprMdl_obj1{j},features); 
    end
    for j=1:6
        [pred_obj2(1,j) pred_obj2_sd(1,j)] = predict(gprMdl_obj2{j},features); 
    end
    
    for k=1:5
        for j=1:6
            prediction_obj1(1,j) = mvnrnd(pred_obj1(j),pred_obj1_sd(j)^2);
            prediction_obj2(1,j) = mvnrnd(pred_obj2(j),pred_obj2_sd(j)^2);
        end
        pred_obj1_pose_diff = [eGetR(prediction_obj1(1,4:6)) prediction_obj1(1,1:3)'; 0 0 0 1];            
        pred_obj2_pose_diff = [eGetR(prediction_obj2(1,4:6)) prediction_obj2(1,1:3)'; 0 0 0 1];
        pred_obj_pose_af{1} = pred_obj1_pose_diff * cur_obj1_pose_af;  
        pred_obj_pose{1} = pred_obj_pose_af{1} * gripper_pose;  
        pred_obj_pose_af{2} = pred_obj2_pose_diff * cur_obj2_pose_af;
        pred_obj_pose{2} = pred_obj_pose_af{2} * gripper_pose;  
%         figure(fig_exp);
%         for i=1:num_obj
%             obj_modelPoints=pred_obj_pose{i}*[obj_modelpoints{i}'; ones(1,size(obj_modelpoints{i},1))];
%             pred_obj_modelpoints{i}=obj_modelPoints(1:3,:)';
%             plot3(pred_obj_modelpoints{i}(:,1),pred_obj_modelpoints{i}(:,2),pred_obj_modelpoints{i}(:,3),'Color',[0 1 1],'Marker','.','Linestyle','none');
%         end

        prediction(1,:) = [pred_obj_pose{1}(1:3,4)' RGete(pred_obj_pose{1}(1:3,1:3))'];
        prediction(2,:) = [pred_obj_pose{2}(1:3,4)' RGete(pred_obj_pose{2}(1:3,1:3))'];
        ground_truth(1,:) = [obj1_end_pose(1:3,4)' RGete(obj1_end_pose(1:3,1:3))']
        ground_truth(2,:) = [obj2_end_pose(1:3,4)' RGete(obj2_end_pose(1:3,1:3))']

        save([dataFolder 'mono_prediction_result_' num2str(scene_num) '_' num2str(k) '.mat'],'prediction','ground_truth');
    end
    
   %pred_obj_pose_diff = [vrrotvec2mat(pred(1,4:7)) pred(1,1:3)'; 0 0 0 1];            


end

