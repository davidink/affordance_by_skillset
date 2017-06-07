%function [feat rslt] = load_test_dataset()
clear all
close all

dataFolder='~/catkin_ws/data/';
dataSaveFolder=[dataFolder 'data_push_result/'];
load([dataFolder 'models/ARmodels/pair_map.mat']);

current_exp_num = size(dir(dataSaveFolder),1)-2;

fig_table = figure;
hold on;

for i=1:current_exp_num
    disp(['loading ' num2str(i) 'th experiment..']);
    load([dataSaveFolder num2str(i) '/push_command.csv']);
    load([dataSaveFolder num2str(i) '/push_result.csv']);
    num_obj = size(push_result,1)/2;
    obj_poses_wID = push_result(1:num_obj,:);
    obj_poses_wID_ap = push_result(num_obj+1:end,:);
    
    %align ids
    if obj_poses_wID_ap(1,1) ~= obj_poses_wID(1,1)
        %swap lines
        temp = obj_poses_wID_ap(1,:);
        obj_poses_wID_ap(1,:) = obj_poses_wID_ap(2,:);
        obj_poses_wID_ap(2,:) = temp;
    end
    
    % obtain pointcloud from model with pose
    [pointcloud_bp pointcloud_norm_bp] = load_pointcloud(obj_poses_wID); %before push
    [pointcloud_ap pointcloud_norm_ap] = load_pointcloud(obj_poses_wID_ap); %after push
    
%     if i==79
    clf(fig_table);
    visualize_objs(fig_table, obj_poses_wID, [0 1 0]); %before push
    visualize_objs(fig_table, obj_poses_wID_ap, [0 0 1]); %after push
    view(0,90)  % XY
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal;
    pause(0.1);
%     end
    % extract features
    for j=1:size(pointcloud_bp,2)
        % rotate pointcloud to pushing frame
        [pf_transform pointcloud_bp_pf{j}] = get_pf_transform(pointcloud_bp{j}, push_command);
        %[pf_transform pointcloud_ap_pf{j}] = get_pf_transform(pointcloud_ap{j}, push_command);
        
        normalpoints = [pf_transform(1:3,1:3) [0 0 0]';0 0 0 1]*[pointcloud_norm_bp{j}'; ones(1,size(pointcloud_norm_bp{j},1))];
        pointcloud_norm_bp_pf{j}=normalpoints(1:3,:)';
        
        
        % 3. Object pose w.r.t push frame
%         push_trans = [[pf_transform(1:3,1:3) [0 0 0]']; 0 0 0 1];
%         obj_cent_pf{j} = obj_poses_wID(j,2:4) - push_command(1,1:3);
%         obj_cent_point = push_trans* [obj_cent_pf{j}'; 1];
%         obj_poses_pf{j} = obj_cent_point(1:3)';
%         obj_cent_pf_ap{j} = obj_poses_wID_ap(j,2:4) - push_command(1,1:3);
%         obj_cent_point_ap = push_trans* [obj_cent_pf_ap{j}'; 1];
%         obj_poses_pf_ap{j} = obj_cent_point_ap(1:3)';
%         obj_pos_diff{j} = obj_poses_pf_ap{j} - obj_poses_pf{j};
         
%         obj_rot = qGetR(obj_poses_wID(j,5:8));
%         obj_rot_ap = qGetR(obj_poses_wID_ap(j,5:8));
%         obj_rot_ang = atan2(obj_rot(2,1),obj_rot(1,1)); % relative to global 'x' axis
%         obj_rot_ang_ap = atan2(obj_rot_ap(2,1),obj_rot_ap(1,1));
%         obj_rot_diff{j} = obj_rot_ang_ap - obj_rot_ang;
        
        obj_poses_pf_bp{j} = pf_transform * [qGetR(obj_poses_wID(j,5:8)) obj_poses_wID(j,2:4)'; 0 0 0 1];
        obj_poses_pf_ap{j} = pf_transform * [qGetR(obj_poses_wID_ap(j,5:8)) obj_poses_wID_ap(j,2:4)'; 0 0 0 1];
                
%         figure;
%         plot3(pointcloud_bp_pf{j}(:,1),pointcloud_bp_pf{j}(:,2),pointcloud_bp_pf{j}(:,3),'Color',[1 0 0 ],'Marker','.','Linestyle','none');
%         %extract_feature(pointcloud_bp_pf, push_command, obj_poses_wID);
%         xlabel('x');
%         ylabel('y');
%         zlabel('z');
%         axis equal
%         hold on;
%         %plotCoord(obj_poses_pf{j}(1:3,4)',obj_poses_pf{j}(1:3,1:3),0.025);
%         quiver3(pointcloud_bp_pf{j}(:,1),pointcloud_bp_pf{j}(:,2),pointcloud_bp_pf{j}(:,3),pointcloud_norm_bp{j}(:,1),pointcloud_norm_bp{j}(:,2),pointcloud_norm_bp{j}(:,3));
        
    end
    % extract features
    features(i,:) = extract_feature(pointcloud_bp_pf, pointcloud_norm_bp_pf, obj_poses_pf_bp, push_command);

    % extract results
    if obj_poses_pf_bp{1}(3,4) < obj_poses_pf_bp{2}(3,4) % which one is closer to the pushing frame
        tar_id = 1;
        obs_id = 2;
    else
        tar_id = 2;
        obs_id = 1;
    end
    tar_displacement(1,1:3) = (obj_poses_pf_ap{tar_id}(1:3,4) - obj_poses_pf_bp{tar_id}(1:3,4))';
    
    tar_rotation = atan2(obj_poses_pf_ap{tar_id}(2,1),obj_poses_pf_ap{tar_id}(3,1)) - atan2(obj_poses_pf_bp{tar_id}(2,1),obj_poses_pf_bp{tar_id}(3,1)); % in 'z' axis
    if tar_rotation > pi
        tar_rotation = tar_rotation - 2*pi;
    else if tar_rotation < -pi
            tar_rotation = tar_rotation + 2*pi;
        end
    end
  
    tar_rot_before = atan2(obj_poses_pf_bp{tar_id}(2,1),obj_poses_pf_bp{tar_id}(3,1))*180/pi;
    tar_rot_after = atan2(obj_poses_pf_ap{tar_id}(2,1),obj_poses_pf_ap{tar_id}(3,1))*180/pi;
    tar_rot_diff_ang = tar_rotation*180/pi;

    obs_displacement(1,1:3) = (obj_poses_pf_ap{obs_id}(1:3,4) - obj_poses_pf_bp{obs_id}(1:3,4))';
    obs_rotation = atan2(obj_poses_pf_ap{obs_id}(2,1),obj_poses_pf_ap{obs_id}(3,1)) - atan2(obj_poses_pf_bp{obs_id}(2,1),obj_poses_pf_bp{obs_id}(3,1)); % in 'z' axis
    if obs_rotation > pi
        obs_rotation = obs_rotation - 2*pi;
    else if obs_rotation < -pi
            obs_rotation = obs_rotation + 2*pi;
        end
    end
    
    obs_rot_before = atan2(obj_poses_pf_bp{obs_id}(2,1),obj_poses_pf_bp{obs_id}(3,1))*180/pi;
    obs_rot_after = atan2(obj_poses_pf_ap{obs_id}(2,1),obj_poses_pf_ap{obs_id}(3,1))*180/pi;
    obs_rot_diff_ang = obs_rotation*180/pi;
    
    if abs(obs_displacement(2)) > 0.02 | abs(obs_displacement(3)) > 0.02 | abs(obs_rotation) > pi/36 % obstacle moved
        interaction = 1;
    else
        interaction = 0;
    end
    
    if tar_rotation > pi/2 | tar_rotation < -pi/2
        something_wrong = 1;
    end
    
    if obs_rotation > pi/2 | obs_rotation < -pi/2
        something_wrong = 1;
    end
    
       
    %save results
    results(i,:) = [tar_displacement tar_rotation obs_displacement obs_rotation];
    results_interaction(i,:) = interaction;
    
    if results(i,2) > 0.09
        stop =1;
    end
   
    
    
    
end

% save features