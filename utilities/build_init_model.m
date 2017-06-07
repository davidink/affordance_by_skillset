close all
clear all

dataFolder='~/catkin_hydro_ws/data/';

features = csvread([dataFolder,'features_init.txt']);
results = csvread([dataFolder,'results_init.txt']);

%features = features(:,1:38);
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
    if (obj_theta_after(i,1) - obj_theta_before(i,1)) > 3.141592/2
    obj_theta_after(i,1) = obj_theta_after(i,1) - 3.141592/2;
    else if (obj_theta_after(i,1) - obj_theta_before(i,1)) < -3.141592/2
        obj_theta_after(i,1) = obj_theta_after(i,1) + 3.141592/2;
        end
    end
end

%filter out infeasible data

num_outlier=0;
for i=1:n
    if abs(obj_position_after(i,1)) > 0.02 | obj_position_after(i,3) < -0.02 | obj_position_after(i,3) - obj_position_before(i,3) < -0.02  
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
    end
end

%% clustering on YZ plane
X = [obj_position_after(:,2)-obj_position_before(:,2) obj_position_after(:,3)-obj_position_before(:,3)];
cluster_num = 2;
%options = statset('Display','final');
%gmmodel = fitgmdist(X,cluster_num,'Options',options);
load 'GMM_3.mat'
cluster_yz = cluster(gmmodel, X);

figure
gscatter(X(:,1),X(:,2),cluster_yz);
hold on

% grouping 
for i=1:cluster_num
    group_cnt(i) = 0;
end
for i=1:size(X,1)
    if cluster_yz(i) == 2 % belongs to group second
        group_cnt(2) = group_cnt(2) +1;
        if group_cnt(2) == 1
            group_b = [obj_position_after(i,2)-obj_position_before(i,2) obj_position_after(i,3)-obj_position_before(i,3) obj_theta_after(i,1)-obj_theta_before(i,1)];
            features_b = [features(i,:)];
        else
            group_b = [group_b; obj_position_after(i,2)-obj_position_before(i,2) obj_position_after(i,3)-obj_position_before(i,3) obj_theta_after(i,1)-obj_theta_before(i,1)];
            features_b = [features_b;features(i,:)];
        end
    else % belong to group first
        group_cnt(1) = group_cnt(1) +1;
        if group_cnt(1) ==1
            group_a = [obj_position_after(i,2)-obj_position_before(i,2) obj_position_after(i,3)-obj_position_before(i,3) obj_theta_after(i,1)-obj_theta_before(i,1)];
            features_a = [features(i,:)];
        else
            group_a = [group_a; obj_position_after(i,2)-obj_position_before(i,2) obj_position_after(i,3)-obj_position_before(i,3) obj_theta_after(i,1)-obj_theta_before(i,1)];
            features_a = [features_a;features(i,:)];
        end 
    end
end
%%

%%No grouping
Y1 = [obj_position_after(:,2)-obj_position_before(:,2) obj_position_after(:,3)-obj_position_before(:,3) obj_theta_after(:,1)-obj_theta_before(:,1)];

small_sigma_squared = 0.0005; %esitmated variance of y distribution
eta_squared = 0.0001; %estimated variance of w (weight) distribution
%[mu lambda] = bayesian_regression([ones(size(features,1),1) features],Y1(:,1),small_sigma_squared,eta_squared);
%[mu lambda] = bayesian_regression(X1,Y1(:,1),small_sigma_squared,eta_squared);
%[mu lambda] = bayesian_regression([ones(size(features,1),1) features],Y1(:,1),small_sigma_squared,eta_squared);
[mu{1} lambda{1}] = bayesian_regression(features,Y1(:,1),small_sigma_squared,eta_squared);

small_sigma_squared = 0.0005; %esitmated variance of y distribution
eta_squared = 0.0001; %estimated variance of w (weight) distribution
[mu{2} lambda{2}] = bayesian_regression(features,Y1(:,2),small_sigma_squared,eta_squared);

small_sigma_squared = 0.0005; %esitmated variance of y distribution
eta_squared = 0.0001; %estimated variance of w (weight) distribution
[mu{3} lambda{3}] = bayesian_regression(features,Y1(:,3),small_sigma_squared,eta_squared);

for k=1:3
    for i=1:size(lambda{k},1)
        for j=i:size(lambda{k},1)
        lambda{k}(j,i) = lambda{k}(i,j);        
        end
    end
end

for i=1:3
    w{i} = mvnrnd(mu{i},lambda{i});
    E_Y{i} = w{i} * features';    
end

E_Yt = [E_Y{1}' E_Y{2}' E_Y{3}'];
plot(E_Yt(:,1),E_Yt(:,2),'xb');

for i=1:size(E_Yt,1)
    plot([Y1(i,1);E_Yt(i,1)],[Y1(i,2);E_Yt(i,2)],'-y');
end

RMSE = rms(Y1 - E_Yt)

% save as a file for weights
weight_sample_num = 100;
for i=1:3
    for j=1:weight_sample_num
        weight = mvnrnd(mu{i},lambda{i});
        fileID = fopen([dataFolder,'weight_samples_',num2str(i),'.txt'],'a+');
        for k=1:size(weight,2)
            fprintf(fileID,'%f,',weight(k));
        end
        fprintf(fileID,'\n');
        fclose(fileID);
    end
end

%%

%% grouping

%group_a
small_sigma_squared = 0.0005; %esitmated variance of y distribution
eta_squared = 0.0001; %estimated variance of w (weight) distribution
[mu_a{1} lambda_a{1}] = bayesian_regression(features_a,group_a(:,1),small_sigma_squared,eta_squared);

small_sigma_squared = 0.0005; %esitmated variance of y distribution
eta_squared = 0.0001; %estimated variance of w (weight) distribution
[mu_a{2} lambda_a{2}] = bayesian_regression(features_a,group_a(:,2),small_sigma_squared,eta_squared);

small_sigma_squared = 0.0005; %esitmated variance of y distribution
eta_squared = 0.0001; %estimated variance of w (weight) distribution
[mu_a{3} lambda_a{3}] = bayesian_regression(features_a,group_a(:,3),small_sigma_squared,eta_squared);

for k=1:3
    for i=1:size(lambda_a{k},1)
        for j=i:size(lambda_a{k},1)
        lambda_a{k}(j,i) = lambda_a{k}(i,j);        
        end
    end
end

for i=1:3
    w_a{i} = mvnrnd(mu_a{i},lambda_a{i});
    E_Y_a{i} = w_a{i} * features_a';    
end

figure;
gscatter(X(:,1),X(:,2),cluster_yz);
hold on

E_Yt_a = [E_Y_a{1}' E_Y_a{2}' E_Y_a{3}'];
plot(E_Yt_a(:,1),E_Yt_a(:,2),'xr');

for i=1:size(E_Yt_a,1)
    plot([group_a(i,1);E_Yt_a(i,1)],[group_a(i,2);E_Yt_a(i,2)],'-y');
end

RMSE_a = rms(group_a - E_Yt_a)

%group_b
small_sigma_squared = 0.0015; %esitmated variance of y distribution
eta_squared = 0.0005; %estimated variance of w (weight) distribution
[mu_b{1} lambda_b{1}] = bayesian_regression(features_b,group_b(:,1),small_sigma_squared,eta_squared);

small_sigma_squared = 0.0005; %esitmated variance of y distribution
eta_squared = 0.0005; %estimated variance of w (weight) distribution
[mu_b{2} lambda_b{2}] = bayesian_regression(features_b,group_b(:,2),small_sigma_squared,eta_squared);

small_sigma_squared = 0.0005; %esitmated variance of y distribution
eta_squared = 0.0005; %estimated variance of w (weight) distribution
[mu_b{3} lambda_b{3}] = bayesian_regression(features_b,group_b(:,3),small_sigma_squared,eta_squared);

for k=1:3
    for i=1:size(lambda_b{k},1)
        for j=i:size(lambda_b{k},1)
        lambda_b{k}(j,i) = lambda_b{k}(i,j);        
        end
    end
end

for i=1:3
    w_b{i} = mvnrnd(mu_b{i},lambda_b{i});
    E_Y_b{i} = w_b{i} * features_b';    
end

E_Yt_b = [E_Y_b{1}' E_Y_b{2}' E_Y_b{3}'];
plot(E_Yt_b(:,1),E_Yt_b(:,2),'xc');

for i=1:size(E_Yt_b,1)
    plot([group_b(i,1);E_Yt_b(i,1)],[group_b(i,2);E_Yt_b(i,2)],'-y');
end

RMSE_b = rms(group_b - E_Yt_b)


