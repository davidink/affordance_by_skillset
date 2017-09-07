function modular_prediction_push(test_trajectory, fwd_gp_model, react_gp_model)
    addpath('utilities');
    step_size = 10;
    fig_gt = figure;
    disp_trajectory(fig_gt,test_trajectory,step_size);
    trajectory = load([test_trajectory '.csv']);

    % Data loading
    Desired_EE_pose_ref_base = trajectory(:,1:6); % relative to robot_base joint
    %convert into global frame 
    Desired_EE_pose = Desired_EE_pose_ref_base + repmat([0 0 0.25 0 1.57 -1.57],size(Desired_EE_pose_ref_base,1),1);
    Executed_EE_pose = trajectory(:,7:12);
    num_obj = trajectory(1,13);

    fig_exp = figure;
    hold on;
    %plot3(Desired_EE_pose(:,1),Desired_EE_pose(:,2),Desired_EE_pose(:,3),'.b');
    plot3(Executed_EE_pose(:,1),Executed_EE_pose(:,2),Executed_EE_pose(:,3),'.r');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal;
    axis([-0.25 0.25 -0.05 1.0 -0.2 0.2]);
    
    for i=1:num_obj
        obj_trajectory{i} = trajectory(:,14+(i-1)*6:14+(i-1)*6+5);
        obj_trajectory_gt{i} = trajectory(:,14+(i-1)*6:14+(i-1)*6+5);
        obj_param{i} = trajectory(:,14+6*num_obj+(i-1)*7:14+6*num_obj+(i-1)*7+6);
        plot3(obj_trajectory{i}(:,1),obj_trajectory{i}(:,2),obj_trajectory{i}(:,3),'.g');
    end

    % Generate obj models
    load('PR2_gripper.mat');
    gripper_points= modelpoints;
    gripper_norms = normalpoints;
    for i=1:num_obj
        granual = 0.005;
        obj_dim(i,:) = ceil(obj_param{i}(1,4:6)*100)/100;
        obj_cent = [0.0 0.0 0];
        [obj_modelpoints{i} obj_normpoints{i}] = create_block_pcd(obj_dim(i,1),obj_dim(i,2),obj_dim(i,3),granual,obj_cent);
    end
    
    % test which obj it interacts
    for i=1:num_obj
        [bool_cont(i) cont_time_step(i)] = test_contact_through_traj(Desired_EE_pose,gripper_points,gripper_norms,obj_modelpoints{i},obj_normpoints{i},obj_trajectory{i}(1,:));
        if bool_cont(i)
            disp(['End effector first contacts object ' num2str(i) ' at time step ' num2str(cont_time_step(i))]);
        end
    end
%     bool_cont(1) = 1;
%     cont_time_step(1) = 50;
%     bool_cont(2) = 0;

    %predict with first object in contact
    for i=1:num_obj
        if bool_cont(i) && cont_time_step(i)==min(cont_time_step)
            %compute expected trajectory of obj1(end-effector) and %obj2(obj1)
            [ee_traj obj_traj] = predict_trajectories(Desired_EE_pose, gripper_points, gripper_norms, obj_modelpoints{i}, obj_normpoints{i}, obj_trajectory{i}(1,:), fwd_gp_model, react_gp_model);
            plot3(ee_traj(:,1),ee_traj(:,2),ee_traj(:,3),'.b');
            plot3(obj_traj(:,1),obj_traj(:,2),obj_traj(:,3),'xb');
            obj_trajectory{i} = obj_traj;
            %exp_ee_traj = interpolate_traj(ee_traj);
        end
    end

    pause(0.5);

    % test for chain reaction
    %load('sample_traj.mat');    
    for i=1:num_obj
        for j=i:num_obj
            if i~=j
                bool_cont_obj = test_contact_through_traj(obj_trajectory{i},obj_modelpoints{i},obj_normpoints{i},obj_modelpoints{j},obj_normpoints{j},obj_trajectory{j}(1,:));
                if bool_cont_obj==1
                    %compute expected trajectory of obj1(end-effector) and %obj2(obj1)
                    [obj1_traj obj2_traj] = predict_trajectories(obj_trajectory{i}, obj_modelpoints{i},obj_normpoints{i}, obj_modelpoints{j},obj_normpoints{j}, obj_trajectory{j}(1,:), fwd_gp_model, react_gp_model);
                    %plot3(ee_traj(:,1),ee_traj(:,2),ee_traj(:,3),'.b');
                    plot3(obj1_traj(:,1),obj1_traj(:,2),obj1_traj(:,3),'oc');
                    plot3(obj2_traj(:,1),obj2_traj(:,2),obj2_traj(:,3),'xc');
                    %exp_ee_traj = interpolate_traj(ee_traj);
                    obj_trajectory{i} = obj1_traj;
                    obj_trajectory{j} = obj2_traj;
                end
            end
        end
    end
    
    testdataFolder = 'data/test/';
    save([testdataFolder 'modular_prediction_result_' test_trajectory(end) '.mat'],'obj_trajectory', 'obj_trajectory_gt');

end

