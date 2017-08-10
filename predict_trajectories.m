function [obj1_trajectory obj2_trajectory] = predict_trajectories(Desired_obj1_pose, obj1_pointcloud, obj1_norms, obj2_pointcloud, obj2_norms, obj2_init_pose, dataFolder)
    
    addpath('utilities');
    step_size = 10;

    % Draw each time step of trajectories
    features = [];
    results = [];
    
    load([dataFolder 'GP_models_react.mat']);
    react_gprMdl = gprMdl;
    load([dataFolder 'GP_models.mat']);

    bool_madecontact = false;
    
    step_cnt = 1;
    for time_step = 10:step_size:size(Desired_obj1_pose,1)-step_size% visualize according to the pose
        if bool_madecontact==false
            exp_obj1_pose = [eGetR(Desired_obj1_pose(time_step,4:6)) Desired_obj1_pose(time_step,1:3)'; 0 0 0 1];
        else
            % if already made a contact, current gripper_pose should be
            nxt_desired_obj1_pose = [eGetR(Desired_obj1_pose(time_step,4:6)) Desired_obj1_pose(time_step,1:3)'; 0 0 0 1];
            exp_obj1_pose = pred_obj1_pose_diff * nxt_desired_obj1_pose;
        end
        
        nxt_exp_obj1_pose = [eGetR(Desired_obj1_pose(time_step+step_size,4:6)) Desired_obj1_pose(time_step+step_size,1:3)'; 0 0 0 1];
        cur_exp_obj1_action_pose = nxt_exp_obj1_pose * exp_obj1_pose^-1;
        
        modelpoints=exp_obj1_pose*[obj1_pointcloud'; ones(1,size(obj1_pointcloud,1))];
        normalpoints = [exp_obj1_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[obj1_norms'; ones(1,size(obj1_norms,1))];
        exp_obj1_points=modelpoints(1:3,:)';
        exp_obj1_norms=normalpoints(1:3,:)';

        if bool_madecontact == false
            % load up initial pose of objects
            obj2_pose_bc = [eGetR(obj2_init_pose(1,4:6)) obj2_init_pose(1,1:3)'; 0 0 0 1];
            obj_modelPoints=obj2_pose_bc*[obj2_pointcloud'; ones(1,size(obj2_pointcloud,1))];
            obj_normalPoints = [obj2_pose_bc(1:3,1:3) [0 0 0]';0 0 0 1]*[obj2_norms'; ones(1,size(obj2_norms,1))];
            cur_obj2_modelpoints=obj_modelPoints(1:3,:)';
            cur_obj2_normalpoints=obj_normalPoints(1:3,:)';

            cur_obj2_pose = obj2_pose_bc;            
            
%             if time_step == size(trajectory,1)
%                 nxt_exp_gripper_pose = exp_gripper_pose; % end of execution
%             else
%                 nxt_exp_gripper_pose = [eGetR(Desired_EE_pose(time_step+step_size,4:6)) Desired_EE_pose(time_step+step_size,1:3)'; 0 0 0 1];
%             end
%             cur_exp_action_pose = nxt_exp_gripper_pose * exp_gripper_pose^-1;
%             modelpoints=exp_gripper_pose*[gripper_points'; ones(1,size(gripper_points,1))];
%             normalpoints = [exp_gripper_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[gripper_norms'; ones(1,size(gripper_norms,1))];
%             exp_gripper_points=modelpoints(1:3,:)';
%             exp_gripper_norms=normalpoints(1:3,:)';
        else
            % if already made a contact, current obj2_pose should be
            % exp_obj2_pose of previous step
            cur_obj2_pose = pred_obj2_pose;
            obj_modelPoints=cur_obj2_pose*[obj2_pointcloud'; ones(1,size(obj2_pointcloud,1))];
            obj_normalPoints = [cur_obj2_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[obj2_norms'; ones(1,size(obj2_norms,1))];
            cur_obj2_modelpoints=obj_modelPoints(1:3,:)';
            cur_obj2_normalpoints=obj_normalPoints(1:3,:)';                    
        end

        %compute contact points
        % compute contact frame
        dist_th = 0.02;
        dot_th = -0.75;
        num_th = 10;

        cont_frame = zeros(4,4);
        cont_frame(4,4) = 1;
        obj1_cent = mean(exp_obj1_points);            
        cont_frame(1:3,4) = obj1_cent';
        cont_frame(1:3,1:3) = cur_exp_obj1_action_pose(1:3,1:3);
    
        [bool_cont cont_pcd_obj1 cont_norm_obj1 cont_pcd_obj2 cont_norm_obj2] = compute_contact_points(exp_obj1_points, exp_obj1_norms, cur_obj2_modelpoints, cur_obj2_normalpoints, cont_frame, dist_th, dot_th, num_th);        
        if bool_cont
            if bool_madecontact == false
                % first time of contact
                bool_madecontact = true;
                first_obj1_to_obj2_dist = sqrt(sum((cur_obj2_pose(1:3,4)-obj1_cent').^2));
            end

            cont_local_frame = zeros(4,4);
            cont_local_frame(4,4) = 1;
            cont_local_frame(1:3,4) = (mean([cont_pcd_obj1;cont_pcd_obj2]))';
            cont_local_frame(1:3,1:3) = cur_exp_obj1_action_pose(1:3,1:3);

            global_shape_features = compute_grid_feature(exp_obj1_points, cur_obj2_modelpoints, cont_frame, 0.25, 0.05);
            local_contact_shape_features = compute_surface_feature(cont_pcd_obj1, cont_norm_obj1, cont_pcd_obj2, cont_norm_obj2, cont_local_frame);
            grid_features = [global_shape_features local_contact_shape_features];            
            kernel_feat = compute_kernel_features([dataFolder 'feat_kernel.mat'],grid_features);

            for j=1:7
                [pred(1,j) pred_sd(1,j)] = predict(gprMdl{j},kernel_feat); 
            end

            %global_shape_features = compute_grid_feature(exp_gripper_points, cur_obj_modelpoints{i},cont_frame, 0.25, 0.05);
            local_contact_shape_features = compute_surface_feature(cont_pcd_obj1, cont_norm_obj1, cont_pcd_obj2, cont_norm_obj2, cont_frame);
            grid_features = [global_shape_features local_contact_shape_features];            
            kernel_feat = compute_kernel_features([dataFolder 'feat_kernel_react.mat'],grid_features);

            for j=1:7
                [pred_react(1,j) pred_sd_react(1,j)] = predict(react_gprMdl{j},kernel_feat); 
            end

            % pose difference is predicted according to the action frame                
            pred_obj2_pose_diff = [vrrotvec2mat(pred(1,4:7)) pred(1,1:3)'; 0 0 0 1];           
            cur_obj2_pose_af = cont_frame*cur_obj2_pose;
            pred_obj2_pose_af = pred_obj2_pose_diff * cur_obj2_pose_af;
            pred_obj2_pose = cont_frame^-1*pred_obj2_pose_af;

            % pose difference of end-effector is w.r.t global frame
            pred_obj1_pose_diff = [vrrotvec2mat(pred_react(1,4:7)) pred_react(1,1:3)';0 0 0 1];
            nxt_step_obj1_pose = [eGetR(Desired_obj1_pose(time_step+step_size,4:6)) Desired_obj1_pose(time_step+step_size,1:3)'; 0 0 0 1];
            pred_obj1_pose = pred_obj1_pose_diff*nxt_step_obj1_pose;

%             for k=1:100
%                 for j=1:7
%                     if k==1
%                         pred_sample(k,j) = pred(1,j);
%                         pred_react_sample(k,j) = pred_react(1,j);
%                     else
%                         pred_sample(k,j) = normrnd(pred(1,j),pred_sd(1,j));                        
%                         pred_react_sample(k,j) = normrnd(pred_react(1,j),pred_sd_react(1,j));
%                     end
%                 end
%             end
% 
%             for k=1:100
%                 pred_obj_pose_diff_sample{k} = [vrrotvec2mat(pred_sample(k,4:7)) pred_sample(k,1:3)';0 0 0 1];
%                 pred_obj_pose_af{k} = pred_obj_pose_diff_sample{k} * cur_obj_pose_af{i};
%                 pred_obj_pose_sample{k} = cont_frame^-1*pred_obj_pose_af{k};
%                 pred_sample_gf(k,:) = [pred_obj_pose_sample{k}(1:3,4)' vrrotmat2vec(pred_obj_pose_sample{k}(1:3,1:3));];
% 
%                 pred_ee_pose_diff_sample{k} = [vrrotvec2mat(pred_react_sample(k,4:7)) pred_react_sample(k,1:3)';0 0 0 1];
%                 pred_ee_pose_sample{k} = pred_ee_pose_diff_sample{k} * nxt_step_gripper_pose;
%                 pred_ee_pose_sample_gf(k,:) = [pred_ee_pose_sample{k}(1:3,4)' vrrotmat2vec(pred_ee_pose_sample{k}(1:3,1:3))];                                    
%             end

%             obj_modelPoints=pred_obj_pose{i}*[obj_modelpoints{i}'; ones(1,size(obj_modelpoints{i},1))];
%             pred_obj_modelpoints{i}=obj_modelPoints(1:3,:)';
%             figure(fig_exp);
%             plot3(pred_obj_modelpoints{i}(:,1),pred_obj_modelpoints{i}(:,2),pred_obj_modelpoints{i}(:,3),'Color',[0 1 1],'Marker','.','Linestyle','none');

%             if time_step ~= size(trajectory,1)
%                 nxt_obj_pose_gt = [eGetR(obj_trajectory{i}(time_step+step_size,4:6)) obj_trajectory{i}(time_step+step_size,1:3)'; 0 0 0 1];
%                 nxt_ee_pose_gt = [eGetR(Executed_EE_pose(time_step+step_size,4:6)) Executed_EE_pose(time_step+step_size,1:3)'; 0 0 0 1];
%                 pred_obj_err(step_cnt,1) = sqrt(sum((pred_obj_pose{i}(1:3,4) - nxt_obj_pose_gt(1:3,4)).^2));
%                 pred_ee_err(step_cnt,1) = sqrt(sum(pred_ee_pose(1:3,4)-nxt_ee_pose_gt(1:3,4)).^2);
%                 step_cnt = step_cnt +1;
%             end

        else
            % if an object did not make any contact 
            pred_obj2_pose = obj2_pose_bc;
        end
        obj1_trajectory(step_cnt,:) = [exp_obj1_pose(1:3,4)' RGete(exp_obj1_pose(1:3,1:3))'];
        obj2_trajectory(step_cnt,:) = [cur_obj2_pose(1:3,4)' RGete(cur_obj2_pose(1:3,1:3))'];
        step_cnt = step_cnt +1;
    end
    obj1_trajectory = interpolate_traj(obj1_trajectory);
    obj2_trajectory = interpolate_traj(obj2_trajectory);
end

