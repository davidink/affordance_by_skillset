function test_trajectory_prediction_with_react(dataFolder, scene_num)
    addpath('utilities');
    %scene_num =8;
    % display ground truth trajectory
    step_size = 10;
    fig_gt = figure;
    %dataFolder = '';
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
        obj_cent = [0.0 0.0 0];
        [obj_modelpoints{i} obj_normpoints{i}] = create_block_pcd(obj_dim(i,1),obj_dim(i,2),obj_dim(i,3),granual,obj_cent);
    end

    % Draw each time step of trajectories
    %step_size = 20;
    features = [];
    results = [];
    load([dataFolder 'GP_models_feature_ard_react.mat']);
    react_gprMdl = gprMdl;
    load([dataFolder 'GP_models_feature_ard.mat']);

    bool_madecontact = false;
    
    step_cnt = 1;
    for time_step = 10:step_size:size(trajectory,1)-step_size% visualize according to the pose

        if bool_madecontact==false
            exp_gripper_pose = [eGetR(Desired_EE_pose(time_step,4:6)) Desired_EE_pose(time_step,1:3)'; 0 0 0 1];
        else
            % if already made a contact, current gripper_pose should be
            exp_gripper_pose = [eGetR(Desired_EE_pose(time_step,4:6)) Desired_EE_pose(time_step,1:3)'; 0 0 0 1];
            exp_gripper_pose = pred_ee_pose_diff * exp_gripper_pose;
        end
        
        nxt_exp_gripper_pose = [eGetR(Desired_EE_pose(time_step+step_size,4:6)) Desired_EE_pose(time_step+step_size,1:3)'; 0 0 0 1];
        
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
                obj_pose_init{i} = [eGetR(obj_trajectory{i}(time_step,4:6)) obj_trajectory{i}(time_step,1:3)'; 0 0 0 1];
                obj_modelPoints=obj_pose_init{i}*[obj_modelpoints{i}'; ones(1,size(obj_modelpoints{i},1))];
                obj_normalPoints = [obj_pose_init{i}(1:3,1:3) [0 0 0]';0 0 0 1]*[obj_normpoints{i}'; ones(1,size(obj_normpoints{i},1))];
                cur_obj_modelpoints{i}=obj_modelPoints(1:3,:)';
                cur_obj_normalpoints{i}=obj_normalPoints(1:3,:)';

                figure(fig_exp);
                plot3(cur_obj_modelpoints{i}(:,1),cur_obj_modelpoints{i}(:,2),cur_obj_modelpoints{i}(:,3),'Color',[0 1 0],'Marker','.','Linestyle','none');

                cur_obj_pose{i} = obj_pose_init{i};            
            end
            
            nxt_exp_gripper_pose = [eGetR(Desired_EE_pose(time_step+step_size,4:6)) Desired_EE_pose(time_step+step_size,1:3)'; 0 0 0 1];

            cur_exp_action_pose = nxt_exp_gripper_pose * exp_gripper_pose^-1;
            modelpoints=exp_gripper_pose*[gripper_points'; ones(1,size(gripper_points,1))];
            normalpoints = [exp_gripper_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[gripper_norms'; ones(1,size(gripper_norms,1))];
            exp_gripper_points=modelpoints(1:3,:)';
            exp_gripper_norms=normalpoints(1:3,:)';
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
                if bool_madecontact == false
                    % first time of contact
                    bool_madecontact = true;
                    first_obj_to_ee_dist = sqrt(sum((cur_obj_pose{i}(1:3,4)-ee_cent').^2));
                end
                
                cont_local_frame = zeros(4,4);
                cont_local_frame(4,4) = 1;
                cont_local_frame(1:3,4) = (mean([cont_pcd_EE;cont_pcd_obj]))';
                cont_local_frame(1:3,1:3) = cur_exp_action_pose(1:3,1:3);

                cur_obj_pose_af{i} = cur_obj_pose{i} * cur_exp_action_pose^-1;
                %global_shape_features = compute_grid_feature(cur_gripper_points, cur_obj_modelpoints{i},cont_frame, 0.25, 0.05);
                global_shape_features = compute_grid_feature(exp_gripper_points, cur_obj_modelpoints{i},cur_obj_pose_af{i}, 0.25, 0.05);
                %global_shape_features = compute_grid_feature(exp_gripper_points, cur_obj_modelpoints{i},cont_frame, 0.25, 0.05);
                local_contact_shape_features = compute_surface_feature(cont_pcd_EE, cont_norm_EE, cont_pcd_obj, cont_norm_obj, cont_local_frame);
                relative_pose_obj_to_obj = cur_obj_pose{i}*cont_frame^-1;
                relative_pose_features = [relative_pose_obj_to_obj(1:3,4)' RGete(relative_pose_obj_to_obj(1:3,1:3))'];
                action_pose_features = [cur_exp_action_pose(1:3,4)' RGete(cur_exp_action_pose(1:3,1:3))'];
                
                %grid_features = [global_shape_features local_contact_shape_features];            
                %kernel_feat = compute_kernel_features([dataFolder 'feat_kernel.mat'],grid_features);
                cur_feat = [global_shape_features local_contact_shape_features relative_pose_features action_pose_features];

                for j=1:6
                    [pred(1,j) pred_sd(1,j)] = predict(gprMdl{j},cur_feat); 
                end
                
                % pose difference is predicted according to the action frame                
                
                %pred_obj_pose_diff = [vrrotvec2mat(pred(1,4:7)) pred(1,1:3)'; 0 0 0 1];           
                pred_obj_pose_diff = [eGetR(pred(1,4:6)) pred(1,1:3)'; 0 0 0 1];           
                cur_obj_pose_af{i} = cont_frame*cur_obj_pose{i};
                pred_obj_pose_af{i} = pred_obj_pose_diff * cur_obj_pose_af{i};
                pred_obj_pose{i} = cont_frame^-1*pred_obj_pose_af{i};

                %global_shape_features = compute_grid_feature(exp_gripper_points, cur_obj_modelpoints{i},cont_frame, 0.25, 0.05);
                local_contact_shape_features = compute_surface_feature(cont_pcd_EE, cont_norm_EE, cont_pcd_obj, cont_norm_obj, cont_frame);
                %grid_features = [global_shape_features local_contact_shape_features]; 
                reactive_pose_features = [pred_obj_pose_af{i}(1:3,4)' RGete(pred_obj_pose_af{i}(1:3,1:3))'];
                %kernel_feat = compute_kernel_features([dataFolder 'feat_kernel_react.mat'],grid_features);
                cur_feat = [global_shape_features local_contact_shape_features relative_pose_features action_pose_features reactive_pose_features];

                for j=1:6
                    [pred_react(1,j) pred_sd_react(1,j)] = predict(react_gprMdl{j},cur_feat); 
                end
                
                % pose difference of end-effector is w.r.t global frame
                %pred_ee_pose_diff = [vrrotvec2mat(pred_react(1,4:7)) pred_react(1,1:3)';0 0 0 1];
                pred_ee_pose_diff = [eGetR(pred_react(1,4:6)) pred_react(1,1:3)';0 0 0 1];

                nxt_step_gripper_pose = [eGetR(Desired_EE_pose(time_step+step_size,4:6)) Desired_EE_pose(time_step+step_size,1:3)'; 0 0 0 1];
                
                pred_ee_pose = pred_ee_pose_diff*nxt_step_gripper_pose;
                                
%                 for k=1:100
%                     for j=1:6
%                         if k==1
%                             pred_sample(k,j) = pred(1,j);
%                             pred_react_sample(k,j) = pred_react(1,j);
%                         else
%                             pred_sample(k,j) = normrnd(pred(1,j),pred_sd(1,j));                        
%                             pred_react_sample(k,j) = normrnd(pred_react(1,j),pred_sd_react(1,j));
%                         end
%                     end
%                 end
%                 
%                 for k=1:100
%                     pred_obj_pose_diff_sample{k} = [vrrotvec2mat(pred_sample(k,4:7)) pred_sample(k,1:3)';0 0 0 1];
%                     pred_obj_pose_af{k} = pred_obj_pose_diff_sample{k} * cur_obj_pose_af{i};
%                     pred_obj_pose_sample{k} = cont_frame^-1*pred_obj_pose_af{k};
%                     pred_sample_gf(k,:) = [pred_obj_pose_sample{k}(1:3,4)' vrrotmat2vec(pred_obj_pose_sample{k}(1:3,1:3));];
%                     
%                     pred_ee_pose_diff_sample{k} = [vrrotvec2mat(pred_react_sample(k,4:7)) pred_react_sample(k,1:3)';0 0 0 1];
%                     pred_ee_pose_sample{k} = pred_ee_pose_diff_sample{k} * nxt_step_gripper_pose;
%                     pred_ee_pose_sample_gf(k,:) = [pred_ee_pose_sample{k}(1:3,4)' vrrotmat2vec(pred_ee_pose_sample{k}(1:3,1:3))];                                    
%                 end
%                 
                
%                 figure;
%                 hold on;
%                 plot3(pred_sample_gf(:,1),pred_sample_gf(:,2),pred_sample_gf(:,3),'.r');
%                 plot3(pred_ee_pose_sample_gf(:,1),pred_ee_pose_sample_gf(:,2),pred_ee_pose_sample_gf(:,3),'.b');
%                 plot3(pred_sample_gf(1,1),pred_sample_gf(1,2),pred_sample_gf(1,3),'ro');
%                 plot3(pred_ee_pose_sample_gf(1,1),pred_ee_pose_sample_gf(1,2),pred_ee_pose_sample_gf(1,3),'rx');
%                 %plot3(pred_react_sample(:,1),pred_react_sample(:,2),pred_react_sample(:,3),'.b');
%                 %plot3(pred(1,1),pred(1,2),pred(1,3),'co');
%                 %plot3(pred_react(1,1),pred_react(1,2),pred_react(1,3),'cx');
%                 plot3(cur_obj_pose{i}(1,4),cur_obj_pose{i}(2,4),cur_obj_pose{i}(3,4),'go');
%                 plot3(ee_cent(1),ee_cent(2),ee_cent(3),'gx');
%                 axis equal;
                
                %cur_obj_to_ee_dist = sqrt(sum((cur_obj_pose{i}(1:3,4)-ee_cent').^2));
                
                cand_num = 0;
                for a=1:100
%                     sam_exp_gripper_pose = pred_ee_pose_sample{a};
%                     modelpoints=sam_exp_gripper_pose*[gripper_points'; ones(1,size(gripper_points,1))];
%                     normalpoints = [sam_exp_gripper_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[gripper_norms'; ones(1,size(gripper_norms,1))];
%                     exp_gripper_points=modelpoints(1:3,:)';
%                     exp_gripper_norms=normalpoints(1:3,:)';
%                     
%                     nxt_obj_pose{i} = pred_obj_pose{i};
%                     obj_modelPoints=nxt_obj_pose{i}*[obj_modelpoints{i}'; ones(1,size(obj_modelpoints{i},1))];
%                     obj_normalPoints = [nxt_obj_pose{i}(1:3,1:3) [0 0 0]';0 0 0 1]*[obj_normpoints{i}'; ones(1,size(obj_normpoints{i},1))];
%                     nxt_obj_modelpoints{i}=obj_modelPoints(1:3,:)';
%                     nxt_obj_normalpoints{i}=obj_normalPoints(1:3,:)';  
%                     nxt_nxt_gripper_pose = [eGetR(Desired_EE_pose(time_step+2*step_size,4:6)) Desired_EE_pose(time_step+2*step_size,1:3)'; 0 0 0 1];                    
%                     nxt_action_pose = nxt_nxt_gripper_pose * sam_exp_gripper_pose^-1;
%                     cont_frame(1:3,4) = mean(exp_gripper_points)';
%                     cont_frame(1:3,1:3) = nxt_action_pose(1:3,1:3);
%                     
%                     [bool_cont] = compute_contact_points(exp_gripper_points, exp_gripper_norms, nxt_obj_modelpoints{i}, nxt_obj_normalpoints{i}, cont_frame, dist_th, dot_th, num_th);
%                     if bool_cont
%                         pred_ee_pose_diff = pred_ee_pose_diff_sample{a};
%                         pred_ee_pose = pred_ee_pose_diff*nxt_step_gripper_pose;
%                         break;
%                     end
                    
                    %dist_obj_to_ee(a,1) = sqrt(sum((pred_sample_gf(1,1:3)-pred_ee_pose_sample_gf(a,1:3)).^2));
%                     for b=a:100
%                         dist_pair(a,b) = sqrt(sum((pred_sample_gf(a,1:3) - pred_ee_pose_sample_gf(b,1:3)).^2));
%                         if abs(cur_obj_to_ee_dist - dist_pair(a,b)) < 0.005
% %                             pred_obj_pose{i} = pred_obj_pose_sample{a};
% %                             pred_ee_pose_diff = pred_ee_pose_diff_sample{b};
%                             break;
%                         end
%                         %dist_piar(b,a) = dist_pair(b,a);
%                     end
                end                
%                 pred_ee_pose_diff = pred_ee_pose_diff_sample{find(dist_obj_to_ee==max(dist_obj_to_ee))};
%                 pred_ee_pose = pred_ee_pose_diff*nxt_step_gripper_pose;
                
                
                obj_modelPoints=pred_obj_pose{i}*[obj_modelpoints{i}'; ones(1,size(obj_modelpoints{i},1))];
                pred_obj_modelpoints{i}=obj_modelPoints(1:3,:)';
                figure(fig_exp);
                plot3(pred_obj_modelpoints{i}(:,1),pred_obj_modelpoints{i}(:,2),pred_obj_modelpoints{i}(:,3),'Color',[0 1 1],'Marker','.','Linestyle','none');
                
                if time_step ~= size(trajectory,1)
                    nxt_obj_pose_gt = [eGetR(obj_trajectory{i}(time_step+step_size,4:6)) obj_trajectory{i}(time_step+step_size,1:3)'; 0 0 0 1];
                    nxt_ee_pose_gt = [eGetR(Executed_EE_pose(time_step+step_size,4:6)) Executed_EE_pose(time_step+step_size,1:3)'; 0 0 0 1];
                    pred_obj_err(step_cnt,1) = sqrt(sum((pred_obj_pose{i}(1:3,4) - nxt_obj_pose_gt(1:3,4)).^2));
                    pred_ee_err(step_cnt,1) = sqrt(sum(pred_ee_pose(1:3,4)-nxt_ee_pose_gt(1:3,4)).^2);
                    step_cnt = step_cnt +1;
                end
                
                
            else
                % if an object did not make any contact 
                pred_obj_pose{i} = obj_pose_init{i};
            end
        end
    end

    % save(['feat_n_result' num2str(scene_num) '.mat'],'features','results');
    % disp('done!');
end

