close all
clear all

[feat_all rslt_all] = load_test_dataset();

% Grouping
X = [rslt_all(:,1) rslt_all(:,2)];
cluster_num = 4;
%options = statset('Display','final');
%gmmodel = fitgmdist(X,cluster_num,'Options',options);
%load 'GMM_5.mat'
%load 'GMM_4.mat'
load 'GMM_8.mat'
cluster_yz = cluster(gmmodel, X);

%grouping 
for i=1:cluster_num
    group_cnt(i) = 0;
end

for i=1:size(X,1)
    for j=1:cluster_num
        if cluster_yz(i) == j % belongs to first group
            group_cnt(j) = group_cnt(j) +1;
            if group_cnt(j) == 1
                rslt{j} = [rslt_all(i,1) rslt_all(i,2) rslt_all(i,3)];
                feat{j} = [feat_all(i,:)];
            else
                rslt{j} = [rslt{j}; rslt_all(i,1) rslt_all(i,2) rslt_all(i,3)];
                feat{j} = [feat{j};feat_all(i,:)];
            end
        end
    end
end

%selected featuers for learning
for i=1:cluster_num
    disp(['group_' num2str(i)]);
    disp('obj feat');
    feat_obj = feat{i}(:,1:14);
    [gp_model_feat_obj{i} rmse_mean_feat_obj{i} rmse_std_feat_obj{i}] = gp_cv(feat_obj,rslt{i},5);
    disp('obj+shape feat');
    feat_objshp = feat{i}(:,1:56);
    [gp_model_feat_objshp{i} rmse_mean_feat_objshp{i} rmse_std_feat_objshp{i}] = gp_cv(feat_objshp,rslt{i},5);
    disp('ctxt feat');
    feat_ctxt = feat{i}(:,57:end);
    [gp_model_feat_ctxt{i} rmse_mean_feat_ctxt{i} rmse_std_feat_ctxt{i}] = gp_cv(feat_ctxt,rslt{i},5);
    disp('all feat');
    [gp_model_feat_all{i} rmse_mean_feat_all{i} rmse_std_feat_all{i}] = gp_cv(feat{i},rslt{i},5);
end

stop = 1;
