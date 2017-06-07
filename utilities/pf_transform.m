function [pf_transform pointcloud_pf] = pf_transform(pointcloud, push_command)
    
    %pointcloud_pf = pointcloud - repmat(push_command(1,1:3),size(pointcloud,1),1);
    % 2. rotate axis to pushing direction
    % rotate axis to be: pushing direction(z), push off direction(y),
    % into the table(x)
    
    %push_start_pose = [qGetR(push_command(1,4:7)) push_command(1,1:3)'; 0 0 0 1]; % not correct for rotation since push frame is in qx qy qz qw, not qw qx qy qz
    %push_end_pose = [qGetR(push_command(2,4:7)) push_command(2,1:3)'; 0 0 0 1];
    push_vec = push_command(2,1:3)' - push_command(1,1:3)';
    rot_sin = push_vec(1)/sqrt(push_vec(1)^2+push_vec(2)^2);
    rot_cos = push_vec(2)/sqrt(push_vec(1)^2+push_vec(2)^2);
%     push_trans = [0 0 -1 0; -rot_cos rot_sin 0 0; rot_sin rot_cos 0 0; 0 0 0 1];
%     modelpoints_pf = push_trans * [pointcloud_pf'; ones(1,size(pointcloud_pf,1))];
    %pointcloud_pf = modelpoints_pf(1:3,:)';
    
    pf_transform  = [0 0 -1 push_command(1,3); -rot_cos rot_sin 0 rot_cos*push_command(1,1)-rot_sin*push_command(1,2); rot_sin rot_cos 0 -rot_sin*push_command(1,1)-rot_cos*push_command(1,2); 0 0 0 1];
    pointcloud_pf = pf_transform *  [pointcloud'; ones(1,size(pointcloud,1))];
    pointcloud_pf = pointcloud_pf(1:3,:)';
%     % 3. Object pose w.r.t push frame
%         obj_cent_pf{j} = obj_poses_wID(j,2:4) - push_command(1,1:3);
%         obj_cent_point = push_trans* [obj_cent_pf{j}'; 1];
%         obj_poses_pf{j} = obj_cent_point(1:3)';
%         obj_cent_pf_ap{j} = obj_poses_wID_ap(j,2:4) - push_command(1,1:3);
%         obj_cent_point_ap = push_trans* [obj_cent_pf_ap{j}'; 1];
%         obj_poses_pf_ap{j} = obj_cent_point_ap(1:3)';
%         obj_pos_diff{j} = obj_poses_pf_ap{j} - obj_poses_pf{j};
%         
%         obj_rot = qGetR(obj_poses_wID(j,5:8));
%         obj_rot_ap = qGetR(obj_poses_wID_ap(j,5:8));
%         obj_rot_ang = atan2(obj_rot(2,1),obj_rot(1,1));
%         obj_rot_ang_ap = atan2(obj_rot_ap(2,1),obj_rot_ap(1,1));
%         obj_rot_diff{j} = obj_rot_ang_ap - obj_rot_ang;
end