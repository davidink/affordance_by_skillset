function result = load_trajectories_kok(dataFolder, scene_num, opt_save_kernel_feat, opt_save_feat_n_result)

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

    for time_step = 10:step_size:size(trajectory,1)-step_size % visualize according to the pose
        gripper_pose = [eGetR(Executed_EE_pose(time_step,4:6)) Executed_EE_pose(time_step,1:3)'; 0 0 0 1];
        if time_step == size(trajectory,1)
            nxt_gripper_pose = gripper_pose; % end of execution
        else
            nxt_gripper_pose = [eGetR(Executed_EE_pose(time_step+step_size,4:6)) Executed_EE_pose(time_step+step_size,1:3)'; 0 0 0 1];
        end
        cur_action_pose = nxt_gripper_pose * gripper_pose^-1;
        modelpoints=gripper_pose*[gripper_points'; ones(1,size(gripper_points,1))];
        normalpoints = [gripper_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[gripper_norms'; ones(1,size(gripper_norms,1))];
        cur_gripper_points=modelpoints(1:3,:)';
        cur_gripper_norms=normalpoints(1:3,:)';

        figure(fig_traj);
        plot3(cur_gripper_points(:,1),cur_gripper_points(:,2),cur_gripper_points(:,3),'Color',[1 0 0],'Marker','.','Linestyle','none');
        hold on;
        plotCoord(gripper_pose(1:3,4)',gripper_pose(1:3,1:3),0.025);
        %quiver3(cur_gripper_points(:,1),cur_gripper_points(:,2),cur_gripper_points(:,3),cur_gripper_norms(:,1)/100,cur_gripper_norms(:,2)/100,cur_gripper_norms(:,3)/100,'Color',[0 0 1]);

        % create pcd of objects
        for i=1:num_obj
            obj_pose = [eGetR(obj_trajectory{i}(time_step,4:6)) obj_trajectory{i}(time_step,1:3)'; 0 0 0 1];
            obj_modelPoints=obj_pose*[obj_modelpoints{i}'; ones(1,size(obj_modelpoints{i},1))];
            obj_normalPoints = [obj_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[obj_normpoints{i}'; ones(1,size(obj_normpoints{i},1))];
            cur_obj_modelpoints{i}=obj_modelPoints(1:3,:)';
            cur_obj_normalpoints{i}=obj_normalPoints(1:3,:)';
            color_grad = 0.5/size(10:step_size:size(trajectory,1),2) * ((time_step-10)/step_size);
            if i ==1 
                default_color = [0 0.5 0];
            else if i==2
                    default_color = [0 0.5 0.5];
                end
            end
            plot3(cur_obj_modelpoints{i}(:,1),cur_obj_modelpoints{i}(:,2),cur_obj_modelpoints{i}(:,3),'Color',default_color +[color_grad color_grad color_grad],'Marker','.','Linestyle','none');
            %quiver3(cur_obj_modelpoints(:,1),cur_obj_modelpoints(:,2),cur_obj_modelpoints(:,3),cur_obj_normalpoints(:,1)/100,cur_obj_normalpoints(:,2)/100,cur_obj_normalpoints(:,3)/100,'Color',[0 0 1]);
        end

        %compute contact points
        % compute contact frame
        dist_th = 0.02;
        dot_th = -0.75;
        num_th = 10;
        for i=1:num_obj
            cont_frame = zeros(4,4);
            cont_frame(4,4) = 1;
            ee_cent = mean(cur_gripper_points);
            obj_cent = mean(cur_obj_modelpoints{i});
            %cont_frame(1:3,4) = mean([ee_cent;obj_cent])';
            cont_frame(1:3,4) = ee_cent';
            cont_frame(1:3,1:3) = cur_action_pose(1:3,1:3);
    %         plotCoord(cont_frame(1:3,4)',cont_frame(1:3,1:3),0.025);
            [bool_cont cont_pcd_EE cont_norm_EE cont_pcd_obj cont_norm_obj] = compute_contact_points(cur_gripper_points, cur_gripper_norms, cur_obj_modelpoints{i}, cur_obj_normalpoints{i}, cont_frame, dist_th, dot_th, num_th);        
            if bool_cont
                cont_local_frame = zeros(4,4);
                cont_local_frame(4,4) = 1;
                cont_local_frame(1:3,4) = (mean([cont_pcd_EE;cont_pcd_obj]))';
                cont_local_frame(1:3,1:3) = cur_action_pose(1:3,1:3);
                
                % extract obj shape grid feature
                % w.r.t obj_pose to action frame
                cur_obj_pose_gf = [eGetR(obj_trajectory{i}(time_step,4:6)) obj_trajectory{i}(time_step,1:3)'; 0 0 0 1];
                cur_obj_pose = [eGetR(obj_trajectory{i}(time_step,4:6)) obj_trajectory{i}(time_step,1:3)'; 0 0 0 1];
                cur_obj_pose_gf(1:3,1:3) = cur_action_pose(1:3,1:3);
                %global_shape_features = compute_grid_feature(cur_gripper_points, cur_obj_modelpoints{i},cont_frame, 0.25, 0.05);
                global_shape_features = compute_grid_feature(cur_gripper_points, cur_obj_modelpoints{i},cur_obj_pose_gf, 0.25, 0.05);
                relative_pose_obj_to_obj = cur_obj_pose*cont_frame^-1;
                relative_pose_features = [relative_pose_obj_to_obj(1:3,4)' RGete(relative_pose_obj_to_obj(1:3,1:3))'];
                action_pose_features = [cur_action_pose(1:3,4)' RGete(cur_action_pose(1:3,1:3))'];
                
%                 local_contact_shape_features = compute_grid_feature(cont_pcd_EE, cont_pcd_obj,cont_local_frame, 0.25, 0.05);
                local_contact_shape_features = compute_surface_feature(cont_pcd_EE, cont_norm_EE, cont_pcd_obj, cont_norm_obj, cont_local_frame);
                grid_features = [global_shape_features local_contact_shape_features];
                all_features = [global_shape_features local_contact_shape_features relative_pose_features action_pose_features];

                if opt_save_kernel_feat
                %save features for building up kernel
                    save([dataFolder 'kernel' num2str(scene_num) '_' num2str(time_step) '.mat'],'global_shape_features','local_contact_shape_features','relative_pose_features','action_pose_features');
                end
                
                if opt_save_feat_n_result
                    %kernel_feat = compute_kernel_features_kok([dataFolder 'feat_kernel.mat'],all_features);
                    obj_pose = [eGetR(obj_trajectory{i}(time_step,4:6)) obj_trajectory{i}(time_step,1:3)'; 0 0 0 1];
                    obj_pose_af = cont_frame*obj_pose;

                    nxt_obj_pose = [eGetR(obj_trajectory{i}(time_step+step_size,4:6)) obj_trajectory{i}(time_step+step_size,1:3)'; 0 0 0 1];
                    nxt_obj_pose_af = cont_frame*nxt_obj_pose;
                    
                    obj_diff = nxt_obj_pose_af * obj_pose_af^-1;
                    %obj_diff_ang = vrrotmat2vec(obj_diff(1:3,1:3));

                    cur_feat = [global_shape_features local_contact_shape_features relative_pose_features action_pose_features];
                    %cur_feat = kernel_feat;
                    cur_result = [obj_diff(1:3,4)' RGete(obj_diff(1:3,1:3))'];

                    features = [features;cur_feat];
                    results = [results;cur_result];                    
                end
            end        
        end
    end

    if opt_save_feat_n_result
       save([dataFolder 'feat_n_result' num2str(scene_num) '.mat'],'features','results');
       disp('done!');
    end
    result = 1;
end

