close all
clear all

[feat_test rslt_test] = load_test_dataset();

dataFolder='~/catkin_ws/data/';

features = csvread([dataFolder,'features_init.txt']);
results = csvread([dataFolder,'results_init.txt']);

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
    if obj_theta_after(i,1) - obj_theta_before(i,1) > 3.141592/2
        obj_theta_after(i,1) = obj_theta_after(i,1) - 3.141592/2;
    else if obj_theta_after(i,1) - obj_theta_before(i,1) < -3.141592/2
        obj_theta_after(i,1) = obj_theta_after(i,1) + 3.141592/2;
        end
    end
end

%filter out infeasible data

num_outlier=0;
for i=1:n
    if abs(obj_position_after(i,1)) > 0.02 | obj_position_after(i,3) < -0.02 | obj_position_after(i,3) - obj_position_before(i,3) < -0.02 | abs(obj_position_after(i,2) - obj_position_before(i,2)) > 0.15  
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

Y1 = [obj_position_after(:,2)-obj_position_before(:,2) obj_position_after(:,3)-obj_position_before(:,3) obj_theta_after(:,1)-obj_theta_before(:,1)];

%% clustering on YZ plane
X = [rslt_test(:,1) rslt_test(:,2)];
cluster_num = 2;
%options = statset('Display','final');
%gmmodel = fitgmdist(X,cluster_num,'Options',options);
load 'GMM_4.mat'
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
            group_b = [rslt_test(i,1) rslt_test(i,2) rslt_test(i,3)];
            features_b = [feat_test(i,:)];
        else
            group_b = [group_b; rslt_test(i,1) rslt_test(i,2) rslt_test(i,3)];
            features_b = [features_b;feat_test(i,:)];
        end
    else % belong to group first
        group_cnt(1) = group_cnt(1) +1;
        if group_cnt(1) ==1
            group_a = [rslt_test(i,1) rslt_test(i,2) rslt_test(i,3)];
            features_a = [feat_test(i,:)];
        else
            group_a = [group_a; rslt_test(i,1) rslt_test(i,2) rslt_test(i,3)];
            features_a = [features_a;feat_test(i,:)];
        end 
    end
end
%%
%generating data sets

[feat_rnd rslt_rnd feat_test rslt_test] = generate_train_test_dataset(features_a,group_a,15);

%random learning
for smp_cnt=3:size(feat_rnd,1)
    feat_usd = feat_rnd(1:smp_cnt,:);
    rslt_usd = rslt_rnd(1:smp_cnt,:);
    
    small_sigma_squared = 0.001; %esitmated variance of y distribution
    eta_squared = 0.001; %estimated variance of w (weight) distribution
    [mu{1} lambda{1}] = bayesian_regression(feat_usd,rslt_usd(:,1),small_sigma_squared,eta_squared);
    [mu{2} lambda{2}] = bayesian_regression(feat_usd,rslt_usd(:,2),small_sigma_squared,eta_squared);
    [mu{3} lambda{3}] = bayesian_regression(feat_usd,rslt_usd(:,3),small_sigma_squared,eta_squared);

    for k=1:3
        for i=1:size(lambda{k},1)
            for j=i:size(lambda{k},1)
            lambda{k}(j,i) = lambda{k}(i,j);        
            end
        end
    end

    for itr=1:50
        for i=1:3
            w{i} = mvnrnd(mu{i},lambda{i});
            %w{i} = mu{i}';
            E_Y{i} = w{i} * feat_test';
        end
        E_Yt = [E_Y{1}' E_Y{2}' E_Y{3}'];
        diff = rslt_test - E_Yt;
        for j=1:size(diff,1)
            if diff(j,3) < -3.141592/2
                diff(j,3) = diff(j,3) + 3.141592/2;
            else if diff(j,3) > 3.141592/2
                    diff(j,3) = diff(j,3) - 3.141592/2;
                end
            end
        end
        RMSE(itr,:) = rms(diff);
    end
    RMSE_rnd_a(smp_cnt-2,:) = mean(RMSE);
end

% feat_rnd = features_b(1:13,:);
% rslt_rnd = group_b(1:13,:);
% 
% feat_test = features_b(14:end,:);
% rslt_test = group_b(14:end,:);

[feat_rnd rslt_rnd feat_test rslt_test] = generate_train_test_dataset(features_b,group_b,15);

%group b
for smp_cnt=3:size(feat_rnd,1)
    feat_usd = feat_rnd(1:smp_cnt,:);
    rslt_usd = rslt_rnd(1:smp_cnt,:);
    
    small_sigma_squared = 0.001; %esitmated variance of y distribution
    eta_squared = 0.001; %estimated variance of w (weight) distribution
    [mu{1} lambda{1}] = bayesian_regression(feat_usd,rslt_usd(:,1),small_sigma_squared,eta_squared);
    [mu{2} lambda{2}] = bayesian_regression(feat_usd,rslt_usd(:,2),small_sigma_squared,eta_squared);
    [mu{3} lambda{3}] = bayesian_regression(feat_usd,rslt_usd(:,3),small_sigma_squared,eta_squared);

    for k=1:3
        for i=1:size(lambda{k},1)
            for j=i:size(lambda{k},1)
            lambda{k}(j,i) = lambda{k}(i,j);        
            end
        end
    end

    for itr=1:50
        for i=1:3
            w{i} = mvnrnd(mu{i},lambda{i});
            %w{i} = mu{i}';
            E_Y{i} = w{i} * feat_test';
        end
        E_Yt = [E_Y{1}' E_Y{2}' E_Y{3}'];
        diff = rslt_test - E_Yt;
        for j=1:size(diff,1)
            if diff(j,3) < -3.141592/2
                diff(j,3) = diff(j,3) + 3.141592/2;
            else if diff(j,3) > 3.141592/2
                    diff(j,3) = diff(j,3) - 3.141592/2;
                end
            end
        end
        RMSE(itr,:) = rms(diff);
    end
    RMSE_rnd_b(smp_cnt-2,:) = mean(RMSE);
end

figure;
subplot(3,1,1);
plot(RMSE_rnd_a(:,1),'-r');
hold on
plot(RMSE_rnd_b(:,1),'-b');
subplot(3,1,2);
plot(RMSE_rnd_a(:,2),'-r');
hold on
plot(RMSE_rnd_b(:,2),'-b');
subplot(3,1,3);
plot(RMSE_rnd_a(:,3),'-r');
hold on
plot(RMSE_rnd_b(:,3),'-b');



