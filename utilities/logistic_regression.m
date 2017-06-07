clear all
close all

load 'gp_a.mat'
load 'gp_b.mat'
load 'gp_nogroup.mat'
load 'GMM_4.mat'

[feat_all rslt_all] = load_test_dataset();

cluster_num = 2;
X = [rslt_all(:,1) rslt_all(:,2)];
cluster_yz = cluster(gmmodel, X);

feat_unary = feat_all(:,1:56);
feat_context = feat_all(:,57:end);

%[B,dev,stats] = mnrfit(feat_unary,cluster_yz);

SVMModel = fitcsvm(feat_all,cluster_yz);