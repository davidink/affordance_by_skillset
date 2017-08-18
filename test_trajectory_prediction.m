function test_trajectory_prediction(dataFolder,scene_num)
    addpath('utilities');
    %scene_num = 9;
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
    axis([-0.25 0.25 -0.05 0.5 -0.2 0.2]);

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
    load([dataFolder 'GP_models_kernel_ard.mat']);
    bool_madecontact = false;

    for time_step = 10:step_size:size(trajectory,1)% visualize according to the pose

        exp_gripper_pose = [eGetR(Desired_EE_pose(time_step,4:6)) Desired_EE_pose(time_step,1:3)'; 0 0 0 1];
        if time_step == size(trajectory,1)
            nxt_exp_gripper_pose = exp_gripper_pose; % end of execution
        else
            nxt_exp_gripper_pose = [eGetR(Desired_EE_pose(time_step+step_size,4:6)) Desired_EE_pose(time_step+step_size,1:3)'; 0 0 0 1];
        end
        cur_exp_action_pose = nxt_exp_gripper_pose * exp_gripper_pose^-1;
        modelpoints=exp_gripper_pose*[gripper_points'; ones(1,size(gripper_points,1))];
        normalpoints = [exp_gripper_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[gripper_norms'; ones(1,size(gripper_norms,1))];
        exp_gripper_points=modelpoints(1:3,:)';
        exp_gripper_norms=normalpoints(1:3,:)';

        figure(fig_exp);
        plot3(exp_gripper_points(:,1),exp_gripper_points(:,2),exp_gripper_points(:,3),'Color',[1 0 0],'Marker','.','Linestyle','none');

        % create pcd of objects

        if bool_madecontact == false
            % load up initial pose of objects
            for i=1:num_obj
                %obj_pose_init{i} = [eGetR(obj_trajectory{i}(time_step,4:6)) obj_trajectory{i}(time_step,1:3)'; 0 0 0 1];
                obj_pose_init{i} = [eGetR(obj_trajectory{i}(1,4:6)) obj_trajectory{i}(1,1:3)'; 0 0 0 1];
                obj_modelPoints=obj_pose_init{i}*[obj_modelpoints{i}'; ones(1,size(obj_modelpoints{i},1))];
                obj_normalPoints = [obj_pose_init{i}(1:3,1:3) [0 0 0]';0 0 0 1]*[obj_normpoints{i}'; ones(1,size(obj_normpoints{i},1))];
                cur_obj_modelpoints{i}=obj_modelPoints(1:3,:)';
                cur_obj_normalpoints{i}=obj_normalPoints(1:3,:)';

                figure(fig_exp);
                plot3(cur_obj_modelpoints{i}(:,1),cur_obj_modelpoints{i}(:,2),cur_obj_modelpoints{i}(:,3),'Color',[0 1 0],'Marker','.','Linestyle','none');

                cur_obj_pose{i} = obj_pose_init{i};            
            end
        else
            % if already made a contact, current obj_pose should be
            % exp_obj_pose of previous step
            for i=1:num_obj
                cur_obj_pose{i} = pred_obj_pose{i};
                obj_modelPoints=cur_obj_pose{i}*[obj_modelpoints{i}'; ones(1,size(obj_modelpoints{i},1))];
                obj_normalPoints = [cur_obj_pose{i}(1:3,1:3) [0 0 0]';0 0 0 1]*[obj_normpoints{i}'; ones(1,size(obj_normpoints{i},1))];
                cur_obj_modelpoints{i}=obj_modelPoints(1:3,:)';
                cur_obj_normalpoints{i}=obj_normalPoints(1:3,:)';            
            end
        end

        %compute contact points
        % compute contact frame
        dist_th = 0.02;
        dot_th = -0.75;
        num_th = 10;
        for i=1:num_obj
            cont_frame = zeros(4,4);
            cont_frame(4,4) = 1;
            ee_cent = mean(exp_gripper_points);
            %obj_cent = mean(cur_obj_modelpoints{i});
            %cont_frame(1:3,4) = mean([ee_cent;obj_cent])';
            cont_frame(1:3,4) = ee_cent';
            cont_frame(1:3,1:3) = cur_exp_action_pose(1:3,1:3);
    %         plotCoord(cont_frame(1:3,4)',cont_frame(1:3,1:3),0.025);
            [bool_cont cont_pcd_EE cont_norm_EE cont_pcd_obj cont_norm_obj] = compute_contact_points(exp_gripper_points, exp_gripper_norms, cur_obj_modelpoints{i}, cur_obj_normalpoints{i}, cont_frame, dist_th, dot_th, num_th);        
            if bool_cont
                cont_local_frame = zeros(4,4);
                cont_local_frame(4,4) = 1;
                cont_local_frame(1:3,4) = (mean([cont_pcd_EE;cont_pcd_obj]))';
                cont_local_frame(1:3,1:3) = cur_exp_action_pose(1:3,1:3);
                
                cur_obj_pose_af = cur_obj_pose{i};                
                cur_obj_pose_af(1:3,1:3) = cur_exp_action_pose(1:3,1:3);

                global_shape_features = compute_grid_feature(exp_gripper_points, cur_obj_modelpoints{i},cur_obj_pose_af, 0.25, 0.05);
                relative_pose_obj_to_obj = cur_obj_pose{i}*cont_frame^-1;
                relative_pose_features = [relative_pose_obj_to_obj(1:3,4)' RGete(relative_pose_obj_to_obj(1:3,1:3))'];
                local_contact_shape_features = compute_surface_feature(cont_pcd_EE, cont_norm_EE, cont_pcd_obj, cont_norm_obj, cont_local_frame);
                grid_features = [global_shape_features local_contact_shape_features];
                action_pose_features = [cur_exp_action_pose(1:3,4)' RGete(cur_exp_action_pose(1:3,1:3))'];

                %kernel_feat = compute_kernel_features([dataFolder 'feat_kernel.mat'],grid_features);
                cur_feat = [global_shape_features local_contact_shape_features relative_pose_features action_pose_features];
                for j=1:6
                    [pred(1,j) pred_sd(1,j)] = predict(gprMdl{j},cur_feat); 
                end

                pred_obj_pose_diff = [vrrotvec2mat(pred(1,4:7)) pred(1,1:3)'; 0 0 0 1];            
                pred_obj_pose{i} = pred_obj_pose_diff * cur_obj_pose{i};  
                obj_modelPoints=pred_obj_pose{i}*[obj_modelpoints{i}'; ones(1,size(obj_modelpoints{i},1))];
                pred_obj_modelpoints{i}=obj_modelPoints(1:3,:)';
                figure(fig_exp);
                plot3(pred_obj_modelpoints{i}(:,1),pred_obj_modelpoints{i}(:,2),pred_obj_modelpoints{i}(:,3),'Color',[0 1 1],'Marker','.','Linestyle','none');

                if bool_madecontact == false
                    bool_madecontact = true;
                end
            else
                % if an object did not make any contact 
                pred_obj_pose{i} = cur_obj_pose{i};
            end
        end
    end

    % save(['feat_n_result' num2str(scene_num) '.mat'],'features','results');
    % disp('done!');
end

