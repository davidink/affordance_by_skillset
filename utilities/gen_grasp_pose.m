function goal_pose = gen_grasp_pose(table_dim)
    load 'SVM_model_grasp_global_feat.mat';
    %SVMmodle_grasp = SVMModel;
    
    %generate a grid over the config space

    [X,Y] = meshgrid(table_dim(1):0.05:table_dim(2),table_dim(3):0.05:table_dim(4));
    points = [X(:) Y(:)];
    %mean_z = mean(feat_glb_pose(:,3));
    test_pose = [points repmat(table_dim(6),size(points,1),1)];
        
    % generate random quaternion
    for i=1:size(points,1)
        smp_theta(i) = rand * 3.141592 * 2;
        smp_theta_rot = smp_theta(i)+ 3.141592;
        m = [0 -sin(smp_theta_rot) cos(smp_theta_rot); 0 cos(smp_theta_rot) sin(smp_theta_rot); -1 0 0];
        qw = sqrt(1 + m(1,1)+m(2,2)+m(3,3))/2;
        test_pose(i,4) = (m(3,2) - m(2,3))/(4*qw);
        test_pose(i,5) = (m(1,3) - m(3,1))/(4*qw);
        test_pose(i,6) = (m(2,1) - m(1,2))/(4*qw);
        test_pose(i,7) = qw; 
    end

    pred_test_pose = predict(SVMModel,test_pose);
    
    graspable_poses = [];
    for i=1:size(test_pose,1)
        if pred_test_pose(i)==1
            graspable_poses = [graspable_poses;test_pose(i,:)];
        end
    end
    
    graspable_poses_inside_table= [];
    for i=1:size(graspable_poses,1)
        if graspable_poses(i,1) > table_dim(1) && graspable_poses(i,2) > table_dim(3)
            graspable_poses_inside_table= [graspable_poses_inside_table;graspable_poses(i,:)];
        end
    end
    
    graspable_pose_far_y = [];
    for i=1:size(graspable_poses_inside_table,1)
        if graspable_poses_inside_table(i,2) == min(graspable_poses_inside_table(:,2))
            graspable_pose_far_y = [graspable_pose_far_y;graspable_poses_inside_table(i,:)];
        end
    end
    
    goal_idx = randperm(size(graspable_pose_far_y,1));
%     figure;
%     gscatter(test_pose(:,1),test_pose(:,2),pred_test_pose,'rg');
%     hold on;
%     plot3([table_dim(1) table_dim(2) table_dim(2) table_dim(1) table_dim(1)],[table_dim(3) table_dim(3) table_dim(4) table_dim(4) table_dim(3)],[table_dim(6) table_dim(6) table_dim(6) table_dim(6) table_dim(6)],'LineWidth',3);
%     xlabel('x');
%     ylabel('y');
%     axis equal;
%     axis([table_dim(1)-0.15 table_dim(2)+0.15 table_dim(3)-0.15 table_dim(4)+0.15]);
%     camroll(90);
%     legend('Not graspable','Graspable','Location','best');
    
    goal_pose = graspable_pose_far_y(goal_idx(1),:);
    
end