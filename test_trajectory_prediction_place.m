function test_trajectory_prediction_place(dataFolder,scene_num)
    
    addpath('utilities');
    % display ground truth trajectory
    step_size = 10;
    fig_gt = figure;
    disp_trajectory_place(fig_gt, [dataFolder 'test_trajectory' num2str(scene_num)], step_size);
    trajectory = load([dataFolder 'test_trajectory' num2str(scene_num) '.csv']);

    % Data loading
    Desired_EE_pose_ref_base = trajectory(:,1:6); % relative to robot_base joint
    %convert into global frame 
    Desired_EE_pose = Desired_EE_pose_ref_base + repmat([0 0 0.3 0 1.57 0],size(Desired_EE_pose_ref_base,1),1);
    Executed_EE_pose = trajectory(:,7:12);
    num_obj = trajectory(1,13);

    for i=1:num_obj
        obj_trajectory{i} = trajectory(:,14+(i-1)*6:14+(i-1)*6+5);
        obj_param{i} = trajectory(:,14+6*num_obj+(i-1)*7:14+6*num_obj+(i-1)*7+6);
    end

    for i=1:num_obj
        granual = 0.005;
        obj_dim(i,:) = ceil(obj_param{i}(1,4:6)*100)/100;
        obj_cent = [0.0 0.0 0];
        [obj_modelpoints{i} obj_normpoints{i}] = create_block_pcd(obj_dim(i,1),obj_dim(i,2),obj_dim(i,3),granual,obj_cent);
    end

    % Draw each time step of trajectories
%     step_size = 20;
    features = [];
    results = [];
    load('data/models/GP_models_place.mat');
   
    figure;
    hold on;
    for i=1:num_obj
        obj_pose = [eGetR(obj_trajectory{i}(1,4:6)) obj_trajectory{i}(1,1:3)'; 0 0 0 1];
        obj_modelPoints=obj_pose*[obj_modelpoints{i}'; ones(1,size(obj_modelpoints{i},1))];
        obj_normalPoints = [obj_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[obj_normpoints{i}'; ones(1,size(obj_normpoints{i},1))];
        cur_obj_modelpoints{i}=obj_modelPoints(1:3,:)';
        cur_obj_normalpoints{i}=obj_normalPoints(1:3,:)';
        if i ==1 
            default_color = [0 0.5 0];
        else if i==2
                default_color = [0 0.5 0.5];
            end
        end
        plot3(cur_obj_modelpoints{i}(:,1),cur_obj_modelpoints{i}(:,2),cur_obj_modelpoints{i}(:,3),'Color',default_color,'Marker','.','Linestyle','none');
        %quiver3(cur_obj_modelpoints(:,1),cur_obj_modelpoints(:,2),cur_obj_modelpoints(:,3),cur_obj_normalpoints(:,1)/100,cur_obj_normalpoints(:,2)/100,cur_obj_normalpoints(:,3)/100,'Color',[0 0 1]);
    end
    
    dist_th = 0.02;
    dot_th = -0.75;
    num_th = 10;

    cont_frame = zeros(4,4);
    cont_frame(4,4) = 1;
    obj1_cent = mean(cur_obj_modelpoints{1});
    obj1_pose = [eGetR(obj_trajectory{1}(1,4:6)) obj_trajectory{1}(1,1:3)'; 0 0 0 1];
    obj2_pose = [eGetR(obj_trajectory{2}(1,4:6)) obj_trajectory{2}(1,1:3)'; 0 0 0 1];
    cont_frame(1:3,4) = obj1_cent';
    cont_frame(1:3,1:3) = obj1_pose(1:3,1:3);

    global_shape_features = compute_grid_feature(cur_obj_modelpoints{1}, cur_obj_modelpoints{2},cont_frame, 0.25, 0.05);
    relative_pose_obj_to_obj = obj2_pose*cont_frame^-1;
    relative_pose_features = [relative_pose_obj_to_obj(1:3,4)' RGete(relative_pose_obj_to_obj(1:3,1:3))'];
    %action_pose_features = [cur_action_pose(1:3,4)' RGete(cur_action_pose(1:3,1:3))'];

%                 local_contact_shape_features = compute_grid_feature(cont_pcd_EE, cont_pcd_obj,cont_local_frame, 0.25, 0.05);
    idx_cont_obj1 = find(cur_obj_normalpoints{1}(:,3)<-0.8);
    cont_pcd_obj1 = cur_obj_modelpoints{1}(idx_cont_obj1,:);
    cont_norm_obj1 = cur_obj_normalpoints{1}(idx_cont_obj1,:);
    idx_cont_obj2 = find(cur_obj_normalpoints{2}(:,3)>0.8);
    cont_pcd_obj2 = cur_obj_modelpoints{2}(idx_cont_obj2,:);
    cont_norm_obj2 = cur_obj_normalpoints{2}(idx_cont_obj2,:);
    
    %local_contact_shape_features = compute_surface_feature(cur_obj_modelpoints{1}, cur_obj_normalpoints{1}, cur_obj_modelpoints{2}, cur_obj_normalpoints{2}, cont_frame);
    local_contact_shape_features = compute_surface_feature(cont_pcd_obj1, cont_norm_obj1, cont_pcd_obj2, cont_norm_obj2, cont_frame);
    all_features = [global_shape_features local_contact_shape_features relative_pose_features];

    for j=1:6
        [pred(1,j) pred_sd(1,j)] = predict(gprMdl{j},all_features); 
    end
    
    obj1_diff = [eGetR(pred(4:6)) pred(1:3)';0 0 0 1];
    
    obj1_pose_af = cont_frame*obj1_pose;
    nxt_obj1_pose_af = obj1_diff * obj1_pose_af;
    nxt_obj1_pose = cont_frame^-1 * nxt_obj1_pose_af;
    
    gt_obj1_pose = [eGetR(obj_trajectory{1}(2,4:6)) obj_trajectory{1}(2,1:3)'; 0 0 0 1];
    pose_err = (nxt_obj1_pose(1:3,4)-gt_obj1_pose(1:3,4));
    
end

