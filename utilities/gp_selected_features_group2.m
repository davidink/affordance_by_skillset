close all
clear all

[feat_all rslt_all] = load_test_dataset();

% Grouping
X = [rslt_all(:,1) rslt_all(:,2)];
cluster_num = 2;
% options = statset('Display','final');
% gmmodel = fitgmdist(X,cluster_num,'Options',options);
%load 'GMM_5.mat'
load 'GMM_4.mat'
cluster_yz = cluster(gmmodel, X);

%grouping 
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

feat_a_obj = features_a(:,1:14); %object features
feat_a_objshp = features_a(:,1:56); %object features + shape
feat_a_ctxt = features_a(:,57:end); %context feat only

feat_b_obj = features_b(:,1:14); %object features
feat_b_objshp = features_b(:,1:56); %object features + shape
feat_b_ctxt = features_b(:,57:end); %context feat only

%% group a gp
% cross validation
disp('obj feat');
[gp_model_feat_a_obj rmse_mean_feat_a_obj rmse_std_feat_a_obj] = gp_cv(feat_a_obj,rslt_a,5);
disp('obj feat + shape');
[gp_model_feat_a_objshp rmse_mean_feat_a_objshp rmse_std_feat_a_objshp] = gp_cv(feat_a_objshp,rslt_a,5);
disp('context only');
[gp_model_feat_a_ctxt rmse_mean_feat_a_ctxt rmse_std_feat_a_ctxt] = gp_cv(feat_a_ctxt,rslt_a,5);
disp('All feat');
[gp_model_feat_a_all rmse_mean_feat_a_all rmse_std_feat_a_all] = gp_cv(features_a,rslt_a,5);
 
%[gp_model_cont cont_rmse_mean cont_rmse_std] = gp_cv(features_a,rslt_a,5);
%gprModel_cnst = build_gp_model(features_a,rslt_a);
 
%% group b gp
% cross validate
disp('obj feat');
[gp_model_feat_b_obj rmse_mean_feat_b_obj rmse_std_feat_b_obj] = gp_cv(feat_b_obj,rslt_b,5);
disp('obj feat + shape');
[gp_model_feat_b_objshp rmse_mean_feat_b_objshp rmse_std_feat_b_objshp] = gp_cv(feat_b_objshp,rslt_b,5);
disp('context only');
[gp_model_feat_b_ctxt rmse_mean_feat_b_ctxt rmse_std_feat_b_ctxt] = gp_cv(feat_b_ctxt,rslt_b,5);
disp('All feat');
[gp_model_feat_b_all rmse_mean_feat_b_all rmse_std_feat_b_all] = gp_cv(features_b,rslt_b,5);
% [gp_model_free free_rmse_mean free_rmse_std] = gp_cv(features_b,rslt_b,5);
% 
% gprModel_free = build_gp_model(features_b,rslt_b);
% 
%% No grouping
%gprModel_nogroup = build_gp_model(feat_all,rslt_all);
% disp('obj feat');
% [gp_model_feat_obj rmse_mean_feat_obj rmse_std_feat_obj] = gp_cv(feat_obj,rslt_all,5);
% disp('obj feat + shape');
% [gp_model_feat_objshp rmse_mean_feat_objshp rmse_std_feat_objshp] = gp_cv(feat_objshp,rslt_all,5);
% disp('context only');
% [gp_model_feat_ctxt rmse_mean_feat_ctxt rmse_std_feat_ctxt] = gp_cv(feat_ctxt,rslt_all,5);
% disp('All feat');
% [gp_model_feat_all rmse_mean_feat_all rmse_std_feat_all] = gp_cv(feat_all,rslt_all,5);
stop = 1;

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
