function features = extract_feature(pointcloud, pointcloud_norm, obj_poses, push_command)
    % object features
    obj_feat = [];
    % pose and orientation
%     figure;
%     hold on;
    for i=1:size(pointcloud,2)
        % orientation relative to 'z' axis of push frame
        obj_ori(i) = atan2(obj_poses{i}(2,1),obj_poses{i}(3,1));
        %obj_ori_pf_ap(i) = atan2(obj_poses_pf_ap{i}(2,1),obj_poses_pf_ap{i}(3,1));
%         plot3(pointcloud{i}(:,1),pointcloud{i}(:,2),pointcloud{i}(:,3),'Color',[0 1 0],'Marker','.','Linestyle','none');
%         xlabel('x');
%         ylabel('y');
%         zlabel('z');
%         axis equal
%         axis([-0.5 0.5 -0.25 0.25 0 0.5]);        
%         view(-90,0);
%         plotCoord(obj_poses{i}(1:3,4)',obj_poses{i}(1:3,1:3),0.025);
        %quiver3(pointcloud{i}(:,1),pointcloud{i}(:,2),pointcloud{i}(:,3),pointcloud_norm{i}(:,1),pointcloud_norm{i}(:,2),pointcloud_norm{i}(:,3));
    end
    
    % find the target pushing object
    if obj_poses{1}(3,4) < obj_poses{2}(3,4)
        if obj_poses{1}(3,4) <0
            tar_id =2; % this object is behind the push
            sec_id =1;
        else
            tar_id = 1;
            sec_id = 2;
        end
    else if obj_poses{2}(3,4) < 0 % this object is behind the push
            tar_id = 1;
            sec_id = 2;
        else
            tar_id = 2;
            sec_id = 1;        
        end
    end
    
    % first obj
    % object pose
    obj_feat = [obj_feat obj_poses{tar_id}(1:3,4)' obj_ori(tar_id)];
    % object dimension
    %obj_feat = [obj_feat max(pointcloud{tar_id})];
    %obj_feat = [obj_feat min(pointcloud{tar_id})];
      
    % second obj
    % object pose
    obj_feat = [obj_feat obj_poses{sec_id}(1:3,4)' obj_ori(sec_id)];
    % object dimension
    %obj_feat = [obj_feat max(pointcloud{sec_id})];
    %obj_feat = [obj_feat min(pointcloud{sec_id})];
  
    % push distance
    %obj_feat = [obj_feat sqrt(sum((push_command(1,1:3)-push_command(2,1:3)).^2))];
    
    % distance until contact & contact shape
    idx1 = find(pointcloud{tar_id}(:,2) > -0.055 & pointcloud{tar_id}(:,2) < -0.04);
    idx2 = find(pointcloud{tar_id}(:,2) > 0.04 & pointcloud{tar_id}(:,2) < 0.055);
    %idx = sort([idx1;idx2]);
    pointcloud_along_push1 = pointcloud{tar_id}(idx1,:);
    pointcloud_along_push2 = pointcloud{tar_id}(idx2,:);
    pointcloud_norm_along_push1 = pointcloud_norm{tar_id}(idx1,:);
    pointcloud_norm_along_push2 = pointcloud_norm{tar_id}(idx2,:);
    idx_twd_contact1 = find(pointcloud_norm_along_push1(:,3) < 0 & abs(pointcloud_norm_along_push1(:,1)) < 0.5);
    idx_twd_contact2 = find(pointcloud_norm_along_push2(:,3) < 0 & abs(pointcloud_norm_along_push2(:,1)) < 0.5);
    pointcloud_contact1 = pointcloud_along_push1(idx_twd_contact1,:);
    pointcloud_contact2 = pointcloud_along_push2(idx_twd_contact2,:);
    pointcloud_norm_contact1 = pointcloud_norm_along_push1(idx_twd_contact1,:);
    pointcloud_norm_contact2 = pointcloud_norm_along_push2(idx_twd_contact2,:);
    
    %plot3(pointcloud_along_push(:,1),pointcloud_along_push(:,2),pointcloud_along_push(:,3),'oc');
    %plot3(pointcloud_contact(:,1),pointcloud_contact(:,2),pointcloud_contact(:,3),'ob');
    %quiver3(pointcloud_contact(:,1),pointcloud_contact(:,2),pointcloud_contact(:,3),pointcloud_norm_contact(:,1),pointcloud_norm_contact(:,2),pointcloud_norm_contact(:,3));
    %options = statset('Display','final');
    %gm = fitgmdist(pointcloud_contact,2,'Options',options);
    %cluster_contact = cluster(gm, pointcloud_contact);
    
    if isempty(idx1) & isempty(idx2)
        dist_contact = 0.5;
        pointcloud_contact1 = [0 0 0.5; 0 0 0.5];
        pointcloud_contact2 = [0 0 0.5; 0 0 0.5];
        pointcloud_norm_contact1 = [0 0 0; 0 0 0];
        pointcloud_norm_contact2 = [0 0 0; 0 0 0];
    else if isempty(idx1)
        dist_contact = min(pointcloud_contact2(:,3));
        pointcloud_contact1 = [0 0 0.5;0 0 0.5];
        pointcloud_norm_contact1 = [0 0 0;0 0 0];
        else if isempty(idx2) | size(pointcloud_contact1,1) == 1
            dist_contact = min(pointcloud_contact1(:,3));
            pointcloud_contact2 = [0 0 0.5;0 0 0.5];
            pointcloud_norm_contact2 = [0 0 0;0 0 0];     
            else
            dist_contact = min(min(pointcloud_contact1(:,3)),min(pointcloud_contact2(:,3)));
            end
        end
    end
    obj_feat = [obj_feat dist_contact];
    
    if size(pointcloud_contact1,1) == 1
        obj_feat = [obj_feat pointcloud_contact1 mean(pointcloud_contact2)];
        obj_feat = [obj_feat pointcloud_norm_contact1 mean(pointcloud_norm_contact2)];
    else if size(pointcloud_contact2,1) == 1
        obj_feat = [obj_feat mean(pointcloud_contact1) pointcloud_contact2];
        obj_feat = [obj_feat mean(pointcloud_norm_contact1) pointcloud_norm_contact2];
        else
        obj_feat = [obj_feat mean(pointcloud_contact1) mean(pointcloud_contact2)]; % mean of contact pose
        obj_feat = [obj_feat mean(pointcloud_norm_contact1) mean(pointcloud_norm_contact2)];
        end
    end
    
    % would be better if we add contact shape features between objects    
    
    % context features
    pointcloud_prj = [];    
    for i=1:size(pointcloud,2)
        if i == tar_id % target object
            pointcloud_prj = [pointcloud_prj;pointcloud{i}(:,2:3) repmat(1,size(pointcloud{i},1),1)];
        else % second obj
            pointcloud_prj = [pointcloud_prj;pointcloud{i}(:,2:3) repmat(2,size(pointcloud{i},1),1)];
        end
    end
    
    % build grid with labels
    len_roi = 0.5;
    grid_size = 0.05;
    num_side = (len_roi/grid_size);
%     figure;
%     axis([-len_roi/2 len_roi/2 0 len_roi]);
%     axis equal;
    for i=1:num_side
        for j=1:num_side
            idx = find(pointcloud_prj(:,1) > (-len_roi/2 + (j-1)*grid_size) & pointcloud_prj(:,1) < (-len_roi/2 + j*grid_size) & pointcloud_prj(:,2) > (i-1)*grid_size & pointcloud_prj(:,2) < i*grid_size); 
            if size(idx,1) ~= 0 % if some points are within the grid
                point_idx = pointcloud_prj(idx,:);
%                 plot(point_idx(:,1),point_idx(:,2),'.r');
                if sum(point_idx(:,3) == 1) > sum(point_idx(:,3) == 2) % grid filled with target obj
                    grid_feat(j,i) = 1;
                else
                    grid_feat(j,i) = 2; % second obj
                end
            else
                grid_feat(j,i) = 3; % free space 
            end
        end
    end
    
    cont_feat = [];
    for i=1:num_side
        for j=1:num_side
            if grid_feat(j,i) == 1
                cont_feat = [cont_feat 1 0 0];
            else if grid_feat(j,i) ==2
                    cont_feat = [cont_feat 0 1 0];
                else cont_feat = [cont_feat 0 0 1];
                end
            end
        end
    end
    
    % context feature from simpler features
%     figure;
%     plot(pointcloud_prj(:,1),pointcloud_prj(:,2),'.r');
%     hold on;
%     axis equal;

    cont_feat_dist = [];
    num_bin = 20;
    bin_size = len_roi/num_bin;
    for j=1:2
        ctxt_feat_dist{j} = [];
        for i=1:num_bin
            idx = find(pointcloud{j}(:,3) > (i-1)*bin_size & pointcloud{j}(:,3) < i*bin_size);
            pointcloud_ins = pointcloud{j}(idx,:);
            if size(idx,1) == 0
                ctxt_feat_dist{j}= [ctxt_feat_dist{j} 0 0];
            else
                ctxt_feat_dist{j}= [ctxt_feat_dist{j} median(pointcloud_ins(:,2)) max(pointcloud_ins(:,2))-min(pointcloud_ins(:,2))];
            end
        %         plot(pointcloud_ins(:,1),pointcloud_ins(:,2),'.g');
        end
    end
    
    cont_feat_dist = [ctxt_feat_dist{tar_id} ctxt_feat_dist{sec_id}];

    
%     for i=1:num_bin
%             idx = find(pointcloud_prj(:,2) > (i-1)*bin_size & pointcloud_prj(:,2) < i*bin_size);
%             pointcloud_ins = pointcloud_prj(idx,:);
%             if size(idx,1) ==0
%                 cont_feat_dist= [cont_feat_dist 0 0];
%             else
%                 cont_feat_dist= [cont_feat_dist median(pointcloud_ins(:,1)) max(pointcloud_ins(:,1))-min(pointcloud_ins(:,1))];
%             end
%         %         plot(pointcloud_ins(:,1),pointcloud_ins(:,2),'.g');
%     end
    

%     k = boundary(pointcloud_prj(:,1),pointcloud_prj(:,2));
%     plot(pointcloud_prj(k,1),pointcloud_prj(k,2),'.b');
%     figure;
%     title('context features');
%     plot(pointcloud_prj(:,1),pointcloud_prj(:,2),'.r');
%     xlabel('y');
%     ylabel('z');
%     axis equal;

    features = [obj_feat cont_feat_dist cont_feat];
    if(size(features,2) ~= 401)
        stop = 1;
    end
    
end