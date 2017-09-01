function result = load_trajectories_place(dataFolder, scene_num, opt_save_kernel_feat, opt_save_feat_n_result)

    addpath('utilities');
    % dataFolder = 'data/trajectories/';
    %scene_num = 4;
    trajectory = load([dataFolder 'trajectory' num2str(scene_num) '.csv']);
    fig_gt = figure;
    step_size = 1;
    disp_trajectory_place(fig_gt, [dataFolder 'trajectory' num2str(scene_num)], step_size);
    
    % Data loading
    Desired_EE_pose_ref_base = trajectory(:,1:6); % relative to robot_base joint
    %convert into global frame 
    Desired_EE_pose = Desired_EE_pose_ref_base + repmat([0 0 0.3 0 1.57 0],size(Desired_EE_pose_ref_base,1),1);
    Executed_EE_pose = trajectory(:,7:12);
    num_obj = trajectory(:,13);

    for i=1:num_obj
        obj_trajectory{i} = trajectory(:,14+(i-1)*6:14+(i-1)*6+5);
        obj_param{i} = trajectory(:,14+6*num_obj+(i-1)*7:14+6*num_obj+(i-1)*7+6);
    end

    % show trajectory of end-effector 
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
    features = [];
    results = [];

    % create pcd of objects
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

    %compute contact points
    % compute contact frame
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

    if opt_save_kernel_feat
    %save features for building up kernel
        save([dataFolder 'kernel' num2str(scene_num) '_' num2str(time_step) '.mat'],'global_shape_features','local_contact_shape_features','relative_pose_features','action_pose_features');
    end

    if opt_save_feat_n_result
        %kernel_feat = compute_kernel_features_kok([dataFolder 'feat_kernel.mat'],all_features);
        %obj_pose = [eGetR(obj_trajectory{1}(1,4:6)) obj_trajectory{i}(time_step,1:3)'; 0 0 0 1];
        obj_pose_af = cont_frame*obj1_pose;

        nxt_obj1_pose = [eGetR(obj_trajectory{1}(2,4:6)) obj_trajectory{1}(2,1:3)'; 0 0 0 1];
        nxt_obj1_pose_af = cont_frame*nxt_obj1_pose;

        obj_diff = nxt_obj1_pose_af * obj_pose_af^-1;
        %obj_diff_ang = vrrotmat2vec(obj_diff(1:3,1:3));

        cur_feat = [global_shape_features local_contact_shape_features relative_pose_features];
        %cur_feat = kernel_feat;
        cur_result = [obj_diff(1:3,4)' RGete(obj_diff(1:3,1:3))'];

        features = [features;cur_feat];
        results = [results;cur_result];                    
    end
       

    if opt_save_feat_n_result
       save([dataFolder 'feat_n_result' num2str(scene_num) '.mat'],'features','results');
       disp('done!');
    end
    result = 1;
    
end

