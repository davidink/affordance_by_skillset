close all
clear all

[feat_all rslt_all] = load_test_dataset();

%% Grouping
X = [rslt_all(:,1) rslt_all(:,2)];
cluster_num = 2;
%options = statset('Display','final');
%gmmodel = fitgmdist(X,cluster_num,'Options',options);
load 'GMM_4.mat'
cluster_yz = cluster(gmmodel, X);

figure
gscatter(X(:,1),X(:,2),cluster_yz,'rb');
hold on

% grouping 
for i=1:cluster_num
    group_cnt(i) = 0;
end
for i=1:size(X,1)
    if cluster_yz(i) == 2 % belongs to group second
        group_cnt(2) = group_cnt(2) +1;
        if group_cnt(2) == 1
            group_b = [rslt_all(i,1) rslt_all(i,2) rslt_all(i,3)];
            features_b = [feat_all(i,:)];
        else
            group_b = [group_b; rslt_all(i,1) rslt_all(i,2) rslt_all(i,3)];
            features_b = [features_b;feat_all(i,:)];
        end
    else % belong to group first
        group_cnt(1) = group_cnt(1) +1;
        if group_cnt(1) ==1
            group_a = [rslt_all(i,1) rslt_all(i,2) rslt_all(i,3)];
            features_a = [feat_all(i,:)];
        else
            group_a = [group_a; rslt_all(i,1) rslt_all(i,2) rslt_all(i,3)];
            features_a = [features_a;feat_all(i,:)];
        end 
    end
end

%% train/test generation & build a base model
[feat_a_train rslt_a_train feat_a_test rslt_a_test] = generate_train_test_dataset(features_a,group_a,20);
[feat_b_train rslt_b_train feat_b_test rslt_b_test] = generate_train_test_dataset(features_b,group_b,15);

feat_a_base_train = feat_a_train(1:5,:);
feat_b_base_train = feat_b_train(1:5,:);
rslt_a_base_train = rslt_a_train(1:5,:);
rslt_b_base_train = rslt_b_train(1:5,:);

mu_a = [mean(rslt_a_base_train(:,1)) mean(rslt_a_base_train(:,2))];
mu_b = [mean(rslt_b_base_train(:,1)) mean(rslt_b_base_train(:,2))];

for i=1:3
    sigma0 = 0.2;
    kparams0 = [5, 0.1];
    gprMdl_a_base{i} = fitrgp(feat_a_base_train, rslt_a_base_train(:,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
    gprMdl_b_base{i} = fitrgp(feat_b_base_train, rslt_b_base_train(:,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
    loss_a_base(i) = sqrt(loss(gprMdl_a_base{i},feat_a_test,rslt_a_test(:,i)));
    loss_b_base(i) = sqrt(loss(gprMdl_b_base{i},feat_b_test,rslt_b_test(:,i)));
end

%% grouping sampling
feat_cand = [feat_a_train(6:end,:);feat_b_train(6:end,:)];
rslt_cand = [rslt_a_train(6:end,:);rslt_b_train(6:end,:)];

rnd_idx = randperm(size(feat_cand,1));

for i=1:size(rnd_idx,2)
    feat_cand_rnd(i,:) = feat_cand(rnd_idx(i),:);
    rslt_cand_rnd(i,:) = rslt_cand(rnd_idx(i),:);
end

feat_cand = feat_cand_rnd;
rslt_cand = rslt_cand_rnd;

feat_a_train = feat_a_base_train;
feat_b_train = feat_b_base_train;
rslt_a_train = rslt_a_base_train;
rslt_b_train = rslt_b_base_train;

for smp_cnt=1:size(feat_cand,1)
    X_cand = feat_cand(smp_cnt,:);
    Y_cand = rslt_cand(smp_cnt,:);
    if smp_cnt == 1
        gprMdl_a = gprMdl_a_base;
        gprMdl_b = gprMdl_b_base;
    end
%     for i=1:3
%         %fit into base model
%         rslt_pred_a(i) = predict(gprMdl_a{i},X_cand);
%         rslt_pred_b(i) = predict(gprMdl_b{i},X_cand);
%     end
    %try to find the group
    dist_a = sqrt( (Y_cand(1)-mu_a(1,1))^2 + (Y_cand(2)-mu_a(1,2))^2);
    dist_b = sqrt( (Y_cand(1)-mu_b(1,1))^2 + (Y_cand(2)-mu_b(1,2))^2);
    if(dist_a < dist_b) % belongs to a
        feat_a_train = [feat_a_train; X_cand];
        rslt_a_train = [rslt_a_train; Y_cand];
    else
        feat_b_train = [feat_b_train; X_cand];
        rslt_b_train = [rslt_b_train; Y_cand];
    end
    %update the GP model
    for i=1:3
        gprMdl_a{i} = fitrgp(feat_a_train, rslt_a_train(:,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
        gprMdl_b{i} = fitrgp(feat_b_train, rslt_b_train(:,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
        loss_a(smp_cnt,i) = sqrt(loss(gprMdl_a{i},feat_a_test,rslt_a_test(:,i)));
        loss_b(smp_cnt,i) = sqrt(loss(gprMdl_b{i},feat_b_test,rslt_b_test(:,i)));
    end
    %update the GM model
    mu_a = [mean(rslt_a_train(:,1)) mean(rslt_a_train(:,2))];
    mu_b = [mean(rslt_b_train(:,1)) mean(rslt_b_train(:,2))];
end

%% No grouping

feat_train = [feat_a_train; feat_b_train];
rslt_train = [rslt_a_train; rslt_b_train];
feat_test = [feat_a_test; feat_b_test];
rslt_test = [rslt_a_test; rslt_b_test];
for smp_cnt=1:size(feat_train,1)
    for i=1:3
        sigma0 = 0.2;
        kparams0 = [5, 0.1];
        gprMdl_init{i} = fitrgp(feat_train(1:smp_cnt,:), rslt_train(1:smp_cnt,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
        loss_init(smp_cnt,i) = sqrt(loss(gprMdl_init{i},feat_test,rslt_test(:,i)));
    end
    %rmse_init(smp_cnt,:) = rms(rslt_test-pred_init);    
end

%% drawing
figure;
subplot(3,1,1);
title('Sample selection');
plot(loss_a(:,1),'r');
hold on;
plot(loss_b(:,1),'b');
hold on;
plot(loss_init(:,1),'g');
ylabel('y');

subplot(3,1,2);
plot(loss_a(:,2),'r');
hold on;
plot(loss_b(:,2),'b');
hold on;
plot(loss_init(:,2),'g');
ylabel('z');

subplot(3,1,3);
plot(loss_a(:,3),'r');
hold on
plot(loss_b(:,3),'b');
hold on;
plot(loss_init(:,3),'g');
ylabel('theta');

xlabel('Number of samples to train model');
legend('Group A','Group B','No group','Location','NorthEast');


