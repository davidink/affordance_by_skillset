function disp_trajectory_grasp(fig_hd, scene_name,step_size)
    addpath('utilities');
    figure(fig_hd);
    trajectory = load([scene_name '.csv']);

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

    hold on;
    plot3(Desired_EE_pose(:,1),Desired_EE_pose(:,2),Desired_EE_pose(:,3),'.r');
    plot3(Executed_EE_pose(:,1),Executed_EE_pose(:,2),Executed_EE_pose(:,3),'.b');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal;
    axis([-0.25 0.25 -0.05 0.7 -0.2 0.4]);

    % show trajectory of end-effector 

    % load the gripper
    %load('PR2_gripper.mat');
%     gripper_points= modelpoints;
%     gripper_norms = normalpoints;
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

    for time_step = 10:step_size:size(trajectory,1)-step_size% visualize according to the pose
        gripper_pose = [eGetR(Executed_EE_pose(time_step,4:6)) Executed_EE_pose(time_step,1:3)'; 0 0 0 1];
        gripper_width = trajectory(time_step,end);        
        [gripper_l_points gripper_l_norms] = create_block_pcd(0.1,0.025,0.025,0.005,[0.05 -gripper_width/2 0]);
        [gripper_r_points gripper_r_norms] = create_block_pcd(0.1,0.025,0.025,0.005,[0.05 gripper_width/2 0]);
        gripper_points = [gripper_l_points;gripper_r_points];
        gripper_norms = [gripper_l_norms;gripper_r_norms];
    
        nxt_gripper_pose = [eGetR(Executed_EE_pose(time_step+step_size,4:6)) Executed_EE_pose(time_step+step_size,1:3)'; 0 0 0 1];

        cur_action_pose = nxt_gripper_pose * gripper_pose^-1;
        modelpoints=gripper_pose*[gripper_points'; ones(1,size(gripper_points,1))];
        normalpoints = [gripper_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[gripper_norms'; ones(1,size(gripper_norms,1))];
        cur_gripper_points=modelpoints(1:3,:)';
        cur_gripper_norms=normalpoints(1:3,:)';

        figure(fig_hd);
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
            %plot3(obj_pose(:,1),obj_pose(:,2),obj_pose(:,3),'Color',default_color +[color_grad color_grad color_grad],'Marker','.','Linestyle','none');
            %quiver3(cur_obj_modelpoints(:,1),cur_obj_modelpoints(:,2),cur_obj_modelpoints(:,3),cur_obj_normalpoints(:,1)/100,cur_obj_normalpoints(:,2)/100,cur_obj_normalpoints(:,3)/100,'Color',[0 0 1]);
        end

    end
end
% save(['feat_n_result' num2str(scene_num) '.mat'],'features','results');
% disp('done!');

