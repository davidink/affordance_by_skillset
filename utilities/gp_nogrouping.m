close all
clear all

[feat_all rslt_all] = load_test_dataset();
feat_all = feat_all(:,1:56);

%% Grouping
X = [rslt_all(:,1) rslt_all(:,2)];
cluster_num = 2;
% options = statset('Display','final');
% gmmodel = fitgmdist(X,cluster_num,'Options',options);
%load 'GMM_5.mat'
load 'GMM_4.mat'
cluster_yz = cluster(gmmodel, X);

% grouping 
for i=1:cluster_num
    group_cnt(i) = 0;
end
for i=1:size(X,1)
    if cluster_yz(i) == 2 % belongs to group second
        group_cnt(2) = group_cnt(2) +1;
        if group_cnt(2) == 1
            rslt_b = [rslt_all(i,1) rslt_all(i,2) rslt_all(i,3)];
            features_b = [feat_all(i,:)];
        else
            rslt_b = [rslt_b; rslt_all(i,1) rslt_all(i,2) rslt_all(i,3)];
            features_b = [features_b;feat_all(i,:)];
        end
    else % belong to group first
        group_cnt(1) = group_cnt(1) +1;
        if group_cnt(1) ==1
            rslt_a = [rslt_all(i,1) rslt_all(i,2) rslt_all(i,3)];
            features_a = [feat_all(i,:)];
        else
            rslt_a = [rslt_a; rslt_all(i,1) rslt_all(i,2) rslt_all(i,3)];
            features_a = [features_a;feat_all(i,:)];
        end 
    end
end

%% group a gp
% [gp_model_cont cont_rmse_mean cont_rmse_std] = gp_cv(features_a,rslt_a,5);
% 
% 
% %% group b gp
% [gp_model_free free_rmse_mean free_rmse_std] = gp_cv(features_b,rslt_b,5);
% 
% stop = 1;

feat_a_part = features_a(1:2:end,:);
rslt_a_part = rslt_a(1:2:end,:);
feat_b_part = features_b(1:2:end,:);
rslt_b_part = rslt_b(1:2:end,:);

% feat_nogroup = [feat_a_part;feat_b_part];
% rslt_nogroup = [rslt_a_part;rslt_b_part];

feat_nogroup = [features_a;features_b];
rslt_nogroup = [rslt_a;rslt_b];


%% No grouping
[gp_model_nogroup nogroup_rmse_mean nogroup_rmse_std] = gp_cv(feat_nogroup,rslt_nogroup,5);
% feat_train = [feat_a_train; feat_b_train];
% rslt_train = [rslt_a_train; rslt_b_train];
% feat_test = [feat_a_test; feat_b_test];
% rslt_test = [rslt_a_test; rslt_b_test];
% for smp_cnt=1:size(feat_train,1)
%     for i=1:3
%         sigma0 = 0.2;
%         kparams0 = [5, 0.1];
%         gprMdl_init{i} = fitrgp(feat_train(1:smp_cnt,:), rslt_train(1:smp_cnt,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
%         loss_init(smp_cnt,i) = sqrt(loss(gprMdl_init{i},feat_test,rslt_test(:,i)));
%     end
%     %rmse_init(smp_cnt,:) = rms(rslt_test-pred_init);    
% end
% 
% %% drawing
% figure;
% subplot(3,1,1);
% title('Sample selection');
% plot(loss_a(:,1),'r');
% hold on;
% plot(loss_b(:,1),'b');
% hold on;
% plot(loss_init(:,1),'g');
% ylabel('y');
% 
% subplot(3,1,2);
% plot(loss_a(:,2),'r');
% hold on;
% plot(loss_b(:,2),'b');
% hold on;
% plot(loss_init(:,2),'g');
% ylabel('z');
% 
% subplot(3,1,3);
% plot(loss_a(:,3),'r');
% hold on
% plot(loss_b(:,3),'b');
% hold on;
% plot(loss_init(:,3),'g');
% ylabel('theta');
% 
% xlabel('Number of samples to train model');
% legend('Group A','Group B','No group','Location','NorthEast');
% 
% 
