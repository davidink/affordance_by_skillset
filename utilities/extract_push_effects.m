close all
clear all

dataFolder='~/catkin_ws/data/data_push';

all_files = dir(dataFolder);
all_dir = all_files([all_files(:).isdir]);
num_dir = numel(all_dir);

for i=1:num_dir-2
    sub_dir = [dataFolder,'/data_push_',num2str(i)];
    cnt = dir([sub_dir,'/obj*']);
    for j=1:length(cnt)
        m = csvread([sub_dir,'/obj_and_tool_frame_',num2str(j),'.txt']);
        if i==1
            results = m;
        else
            results = [results; m];
        end        
    end
end

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
end

%flip direction
for i=1:n
    if obj_orientation_before(i,3) < 0
        obj_orientation_before(i,2) = -obj_orientation_before(i,2);
        obj_orientation_before(i,3) = -obj_orientation_before(i,3);
    end
    if obj_orientation_after(i,3) < 0
        obj_orientation_after(i,2) = -obj_orientation_after(i,2);
        obj_orientation_after(i,3) = -obj_orientation_after(i,3);
    end    
end

for i=1:n
    obj_theta_before(i,1) = atan2(obj_orientation_before(i,2),obj_orientation_before(i,3));
    obj_theta_after(i,1) = atan2(obj_orientation_after(i,2),obj_orientation_after(i,3));
end

%filter out infeasible data

num_outlier=0;
for i=1:n
    if abs(obj_position_after(i,2)) > 0.2 | obj_position_after(i,3) < 0  | obj_position_after(i,3)-obj_position_before(i,3) > 0.12 | obj_position_after(i,3)-obj_position_before(i,3) < 0.05
        num_outlier= num_outlier+1; 
        outlier(num_outlier) = i;
    end
end

if num_outlier ~=0
    for i=1:num_outlier
        remov = outlier(num_outlier-i+1);
        obj_position_before(remov,:) = [];
        obj_position_after(remov,:) = [];
        obj_orientation_before(remov,:) = [];
        obj_orientation_after(remov,:) = [];
        obj_theta_before(remov,:) = [];
        obj_theta_after(remov,:) = [];
    end
end

%plot(obj_theta_before(:,1), obj_theta_after(:,1),'.r');
%hold on

X1 = [obj_theta_before obj_theta_after];
X2 = [obj_theta_before obj_position_after(:,2)-obj_position_before(:,2)];
X3 = [obj_theta_before obj_position_after(:,3)-obj_position_before(:,3)];

subplot(3,1,1)
plot(X1(:,1),X1(:,2),'.r');
ylabel('theta_after');
subplot(3,1,2)
plot(X2(:,1),X2(:,2),'.g');
ylabel('y_after');
subplot(3,1,3)
plot(X3(:,1),X3(:,2),'.b');
ylabel('z_after');
xlabel('theta_before');

cluster_num = 2;
%options = statset('Display','final');
%gmmodel = fitgmdist(X,2,'Options',options);
load GMM_2.mat
cluster_theta = cluster(gmmodel, X1);

figure
subplot(3,1,1);
gscatter(X1(:,1),X1(:,2),cluster_theta);
hold on
subplot(3,1,2);
gscatter(X2(:,1),X2(:,2),cluster_theta);
hold on
subplot(3,1,3);
gscatter(X3(:,1),X3(:,2),cluster_theta);
hold on

% grouping 
for i=1:cluster_num
    group_cnt(i) = 0;
end
for i=1:size(X1,1)
    if cluster_theta(i) == 2 %belong to group second
        group_cnt(2) = group_cnt(2) +1;
        if group_cnt(2) == 1
            group_sec = [X1(i,1) X1(i,2) X2(i,2) X3(i,2)];
        else
            group_sec = [group_sec; X1(i,1) X1(i,2) X2(i,2) X3(i,2)];
        end
    else % belong to group first
        group_cnt(1) = group_cnt(1) +1;
        if group_cnt(1) ==1
            group_fir = [X1(i,1) X1(i,2) X2(i,2) X3(i,2)];
        else
            group_fir = [group_fir; X1(i,1) X1(i,2) X2(i,2) X3(i,2)];
        end 
    end
end

small_sigma_squared = 0.02; %esitmated variance of y (theta) distribution
eta_squared = 0.01; %estimated variance of w (weight) distribution
[mu{1} lambda{1}] = bayesian_regression([ones(size(group_sec,1),1) group_sec(:,1)],group_sec(:,2),small_sigma_squared,eta_squared);
small_sigma_squared = 0.001;
[mu{2} lambda{2}] = bayesian_regression([ones(size(group_sec,1),1) group_sec(:,1)],group_sec(:,3),small_sigma_squared,eta_squared);
small_sigma_squared = 0.001;
eta_squared = 0.005;
[mu{3} lambda{3}] = bayesian_regression([ones(size(group_sec,1),1) group_sec(:,1)],group_sec(:,4),small_sigma_squared,eta_squared);

x1 = linspace(-2,2);
for i=1:size(x1,2)
    for j=1:3
        a(j,:) = mvnrnd(mu{j},lambda{j});
        y(j,i) = a(j,2) * x1(i) + a(j,1);
    end
end

subplot(3,1,1);
plot(x1,y(1,:),'.g');
subplot(3,1,2);
plot(x1,y(2,:),'.g');
subplot(3,1,3);
plot(x1,y(3,:),'.g');

sss = 0.01;
es = 0.01;
[mu_ori lambda_ori] = bayesian_regression(obj_before_push, obj_after_push,sss,es);


% x_input = [obj_position_before(:,2) obj_position_before(:,3) obj_theta_before(:,1)];
% y_output = [obj_position_after(:,2) obj_position_after(:,3) obj_theta_after(:,1)];
% 
% figure;
% plot3(x_input(:,1), x_input(:,2), x_input(:,3),'.g');
% hold on
% plot3(y_output(:,1), y_output(:,2), y_output(:,3),'.b');
% xlabel('y');
% ylabel('z');
% zlabel('theta');
% 
% options = statset('Display','final');
% gmmodel_y = fitgmdist(y_output,2,'Options',options);
% cluster_y = cluster(gmmodel_y, y_output);
% 
% cnt_one = 0;
% cnt_two= 0;
% for i=1:size(cluster_y,1)
%     if cluster_y(i) == 1;
%         cnt_one = cnt_one + 1;
%         group_one(cnt_one,:) = y_output(i,:);
%     else cluster_y(i) == 2  ;      
%         cnt_two = cnt_two + 1;
%         group_two(cnt_two,:) = y_output(i,:);
%     end
% end
% 
% figure;
% plot3(group_one(:,1), group_one(:,2), group_one(:,3),'.g');
% hold on;
% plot3(group_two(:,1), group_two(:,2), group_two(:,3),'.b');



%h2 = ezcontour(@(x,y,z)pdf(gmmodel_y,[x y z]),[-1 1],[-1 1],[-2 2]);

% plot3(obj_position_before(:,2), obj_position_before(:,3), obj_theta_before(:,1),'.');
% xlabel('y');
% ylabel('z');
% zlabel('theta');
% hold all;
% plot3(obj_position_after(:,2), obj_position_after(:,3), obj_theta_after(:,1),'.');
%plot3(obj_position_after(:,2)-obj_position_before(:,2),obj_position_after(:,3)-obj_position_before(:,3),obj_theta_after(:,1)-obj_theta_before(:,1),'.');

% subplot(3,1,1);
% subplot(3,1,2);

% plot3(obj_position_before(:,1),obj_position_before(:,2),obj_position_before(:,3),'.r');
% hold on
% plot3(obj_position_after(:,1),obj_position_after(:,2),obj_position_after(:,3),'.b');
% 
% xlabel('x');
% ylabel('y');
% zlabel('z');

%quiver(obj_pose_before_push_frame(:).transform(2,4),obj_pose_before_push_frame(:).transform(3,4),obj_pose_before_push_frame(:).transform(2,3),obj_pose_before_push_frame(:).transform(3,3));
% quiver3(obj_position_before(:,1),obj_position_before(:,2),obj_position_before(:,3),obj_orientation_before(:,1),obj_orientation_before(:,2),obj_orientation_before(:,3));

% for i=1:size(push_frame,1)
%     obj_position_before_push_frame(i,1) = obj_pose_before_push_frame(i).transform(1,4);
%     obj_position_before_push_frame(i,2) = obj_pose_before_push_frame(i).transform(2,4);
%     obj_position_before_push_frame(i,3) = obj_pose_before_push_frame(i).transform(3,4);
% %     obj_pose_before_push_frame(i).transform(2,4) obj_pose_before_push_frame(i).transform(3,4)]
% end
% 
% push_sample=cell(n,1)
% for i=1:size(push_frame,1)
%     push_sample{i} = [push_frame(i,1) push_frame(i,2)]
%    
% end
% 
% push_sample{i}(1,1)

% plot3(obj_position_before_push_frame(:,1), obj_position_before_push_frame(:,2),obj_position_before_push_frame(:,3),'.r');
% xlabel('x');
% ylabel('y');
% zlabel('z');

%quiver(push_frame(:,1),push_frame(:,2),push_ori(:,1),push_ori(:,2));
% quiver(zeros(size(push_frame,1),1),zeros(size(push_frame,1),1),push_ori(:,1),push_ori(:,2));
% xlabel('x');
% ylabel('y');
% view(-90,90);


%scatter(push_frame(:,1), push_frame(:,2));

% FID = fopen([dataFolder,'trackerparam_temp.txt'],'w');
% for i=1:length(modelNames)
%    [modelpoints modelcolors]=loadCloudColor([modelFolder,modelNames{i},'.pcd']);
% Ang1=modelOri(i,1)*pi/180; %Rotate around table normal
% Ang2=modelOri(i,2)*pi/180; %Rotate around local y-axis
% Ang3=modelOri(i,3)*pi/180; %Rotate around 
% 
% 
% 
% 
% 
% Ra=[cos(Ang1) -sin(Ang1) 0; sin(Ang1) cos(Ang1) 0; 0 0 1];
% Rb=[cos(Ang2) 0 sin(Ang2); 0 1 0; -sin(Ang2) 0 cos(Ang2)];
% Rc=[cos(Ang3) -sin(Ang3) 0; sin(Ang3) cos(Ang3) 0; 0 0 1];
% R=Ra*Rb*Rc;
% modelpoints=bsxfun(@plus, modelpoints*R', modelPos(i,:));
% plot3(modelpoints(:,1),modelpoints(:,2),modelpoints(:,3),'.')
% 
% qw=sqrt(1+R(1,1)+R(2,2)+R(3,3))/2
% qx=(R(3,2)-R(2,3))/(4*qw)
% qy=(R(1,3)-R(3,1))/(4*qw)
% qz=(R(2,1)-R(1,2))/(4*qw)
% 
% q=[qw qx qy qz];
% 
% 
% fwrite(FID,['Tracker', num2str(i),', ',modelNames{i},'_wNorm, ',num2str(modelPos(i,1)),', ',num2str(modelPos(i,2)),', ',num2str(modelPos(i,3)),', ',num2str(qw),', ',num2str(qx),', ',num2str(qy),', ',num2str(qz), char(10) ])
%     
% end