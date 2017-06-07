function [feat rslt] = load_test_dataset()

dataFolder='~/catkin_ws/data/';

features = csvread([dataFolder,'features_all.txt']);
results = csvread([dataFolder,'results_all.txt']);

%features = features(:,1:497);

push_frame = results(1:3:end,:);
obj_before_push = results(2:3:end,:);
obj_after_push = results(3:3:end,:);

push_quat = push_frame(:,4:7);

n = size(push_frame,1);
for i=1:size(push_frame,1)
    rotm = quat2rotm(push_quat(i,:));
    push_transform = [rotm [push_frame(i,1) push_frame(i,2) push_frame(i,3)]'; 0 0 0 1];
    push_transforms(i).transform = push_transform;
    if i==1
        push_ori = [rotm(1,3) rotm(2,3)];
    else
        push_ori = [push_ori; rotm(1,3) rotm(2,3)];
    end
end

for i=1:size(obj_before_push,1)
    rotm = quat2rotm(obj_before_push(i,4:7));
    obj_transform_before = [rotm [obj_before_push(i,1) obj_before_push(i,2) obj_before_push(i,3)]';0 0 0 1];
    obj_transforms_before(i).transform = obj_transform_before;    
end

for i=1:size(obj_after_push,1)
    rotm = quat2rotm(obj_after_push(i,4:7));
    obj_transform_after = [rotm [obj_after_push(i,1) obj_after_push(i,2) obj_after_push(i,3)]';0 0 0 1];
    obj_transforms_after(i).transform = obj_transform_after;    
end

for i=1:size(push_frame,1)
    obj_pose_before_push_frame(i).transform = push_transforms(i).transform\ obj_transforms_before(i).transform;
    obj_pose_after_push_frame(i).transform = push_transforms(i).transform\ obj_transforms_after(i).transform;
end

for i=1:size(push_frame,1)
    obj_position_before(i,1) = obj_pose_before_push_frame(i).transform(1,4);
    obj_position_before(i,2) = obj_pose_before_push_frame(i).transform(2,4);
    obj_position_before(i,3) = obj_pose_before_push_frame(i).transform(3,4);
    obj_position_after(i,1) = obj_pose_after_push_frame(i).transform(1,4);
    obj_position_after(i,2) = obj_pose_after_push_frame(i).transform(2,4);
    obj_position_after(i,3) = obj_pose_after_push_frame(i).transform(3,4);
    %consider orientation on y-z plane = theta
    obj_orientation_before(i,1) = 0;
    obj_orientation_before(i,2) = obj_pose_before_push_frame(i).transform(2,3);
    obj_orientation_before(i,3) = obj_pose_before_push_frame(i).transform(3,3);
    obj_orientation_after(i,1) = 0;
    obj_orientation_after(i,2) = obj_pose_after_push_frame(i).transform(2,3);
    obj_orientation_after(i,3) = obj_pose_after_push_frame(i).transform(3,3);
    obj_pose_before_push_frame_quat(i,:) = [obj_pose_before_push_frame(i).transform(1:3,4)' rotm2quat(obj_pose_before_push_frame(i).transform(1:3,1:3))];
    obj_pose_after_push_frame_quat(i,:) = [obj_pose_after_push_frame(i).transform(1:3,4)' rotm2quat(obj_pose_after_push_frame(i).transform(1:3,1:3))];
end

%flip direction
% for i=1:n
%     if obj_orientation_before(i,3) < 0
%         obj_orientation_before(i,2) = -obj_orientation_before(i,2);
%         obj_orientation_before(i,3) = -obj_orientation_before(i,3);
%     end
%     if obj_orientation_after(i,3) < 0
%         obj_orientation_after(i,2) = -obj_orientation_after(i,2);
%         obj_orientation_after(i,3) = -obj_orientation_after(i,3);
%     end    
% end

for i=1:n
    obj_theta_before(i,1) = atan2(obj_orientation_before(i,2),obj_orientation_before(i,3));
    obj_theta_after(i,1) = atan2(obj_orientation_after(i,2),obj_orientation_after(i,3));
    if obj_theta_after(i,1) - obj_theta_before(i,1) > pi
         obj_theta_after(i,1) = obj_theta_after(i,1) - 2*pi;
    else if (obj_theta_after(i,1) - obj_theta_before(i,1)) < -pi
         obj_theta_after(i,1) = obj_theta_after(i,1) + 2*pi;
         end
    end
end



%filter out infeasible data

num_outlier=0;
for i=1:n
    if abs(obj_position_after(i,1)) > 0.02 | obj_position_after(i,3) < -0.02 | obj_position_after(i,3) - obj_position_before(i,3) < -0.02 | abs(obj_position_after(i,2) - obj_position_before(i,2)) > 0.08 | abs(obj_theta_after(i,1) - obj_theta_before(i,1)) > pi/2  
        num_outlier= num_outlier+1; 
        outlier(num_outlier) = i;
    end
end

if num_outlier ~=0
    for i=1:num_outlier
        remov = outlier(num_outlier-i+1);
        features(remov,:) = [];
        obj_position_before(remov,:) = [];
        obj_position_after(remov,:) = [];
        obj_orientation_before(remov,:) = [];
        obj_orientation_after(remov,:) = [];
        obj_theta_before(remov,:) = [];
        obj_theta_after(remov,:) = [];
        push_frame(remov,:) = [];
        obj_before_push(remov,:) = [];
        obj_after_push(remov,:) = [];
        obj_pose_before_push_frame_quat(remov,:) = [];
        obj_pose_after_push_frame_quat(remov,:) = [];
    end
end

% fig_pose = figure;
% fig_push = figure;
% for i=1:size(push_frame,1)
%     clf(fig_pose);
%     figure(fig_pose);
%     plot3(push_frame(i,1),push_frame(i,2),push_frame(i,3),'*k','MarkerSize',10);
%     hold on;
%     modelpoints = visualize_model('Book1',obj_before_push(i,:),'r');
%     plot3(modelpoints(:,1),modelpoints(:,2),modelpoints(:,3),'Color','r','Marker','.','Linestyle','none');
%     modelpoints = visualize_model('Book1',obj_after_push(i,:),'r');
%     plot3(modelpoints(:,1),modelpoints(:,2),modelpoints(:,3),'Color','b','Marker','.','Linestyle','none');
%     axis equal;
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     view([0 0 1]);
%     camroll(90);
%     legend('Pushing position','Before','After','Location','bestoutside');
%     figure(fig_push);
%     clf(fig_push);
%     hold on;
%     modelpoints = visualize_model('Book1',obj_pose_before_push_frame_quat(i,:),'r');
%     plot3(modelpoints(:,1),modelpoints(:,2),modelpoints(:,3),'Color','r','Marker','.','Linestyle','none');
%     modelpoints = visualize_model('Book1',obj_pose_after_push_frame_quat(i,:),'b');
%     plot3(modelpoints(:,1),modelpoints(:,2),modelpoints(:,3),'Color','b','Marker','.','Linestyle','none');
%     vec = [obj_position_before(i,:);obj_position_before(i,:)+obj_orientation_before(i,:)];
%     plot3(vec(:,1),vec(:,2),vec(:,3),'-r','LineWidth',3);
%     vec = [obj_position_after(i,:);obj_position_after(i,:)+obj_orientation_after(i,:)];
%     plot3(vec(:,1),vec(:,2),vec(:,3),'-b','LineWidth',3);
%     
%     axis equal;
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     view([-1 0 0]);
%     %camroll(90);
%     axis([-0.1 0.1 -0.2 0.2 0 0.5]);
%     
%     delta_theta = (obj_theta_after(i,:)-obj_theta_before(i,:))/pi*180    
% end

feat = features;
rslt = [obj_position_after(:,2)-obj_position_before(:,2) obj_position_after(:,3)-obj_position_before(:,3) obj_theta_after(:,1)-obj_theta_before(:,1)];

% fig_push = figure;
% fig_cont = figure;
% for i=1:size(rslt,1)
%     draw_contextual_features(feat(i,:),fig_cont);
%     figure(fig_push);
%     hold on;
%     plot(rslt(i,1),rslt(i,2),'.r');
% end

end