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

figure
gscatter(X(:,1),X(:,2),cluster_yz,'rb');
legend('Constraint Push','Free-space Push','Location','northeast');
xlabel('y (horizontal) (m) ');
ylabel('z (vertical) (m)');
hold on
axis equal;
set(gca,'Xdir','reverse')

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
disp('computing GP model for constraint push');
gprModel_cnst = build_gp_model(features_a,rslt_a);
%[gp_model_cont cont_rmse_mean cont_rmse_std] = gp_cv(features_a,rslt_a,5);


%% group b gp
disp('computing GP model for free push');
gprModel_free = build_gp_model(features_b,rslt_b);

%% no grouping
disp('computing GP model for no grouping');
gprModel_all = build_gp_model(feat_all,rslt_all);

%% Compute prediction for all input features
%load('gp_grouping_initial_model.mat');
for i=1:3
    [pred_cnst_mean(:,i) pred_cnst_sd(:,i)] = predict(gprModel_cnst{i},features_a);
    [pred_free_mean(:,i) pred_free_sd(:,i)] = predict(gprModel_free{i},features_b);
end

figure;
hold on;
plot(rslt_a(:,1),rslt_a(:,2),'.r');
plot(rslt_b(:,1),rslt_b(:,2),'.b');

plot(pred_cnst_mean(:,1),pred_cnst_mean(:,2),'or');
plot(pred_free_mean(:,1),pred_free_mean(:,2),'ob');

for i=1:size(features_a,1)
    plot([rslt_a(i,1) pred_cnst_mean(i,1)],[rslt_a(i,2) pred_cnst_mean(i,2)],'y');    
end

for i=1:size(features_b,1)
    plot([rslt_b(i,1) pred_free_mean(i,1)],[rslt_b(i,2) pred_free_mean(i,2)],'y');
end

for i=1:3
    [pred_cnst_cross_mean(:,i) pred_cnst_cross_sd(:,i)] = predict(gprModel_cnst{i},features_b);
    [pred_free_cross_mean(:,i) pred_free_cross_sd(:,i)] = predict(gprModel_free{i},features_a);
end

plot(pred_cnst_cross_mean(:,1),pred_cnst_cross_mean(:,2),'xr');
plot(pred_free_cross_mean(:,1),pred_free_cross_mean(:,2),'xb');

for i=1:size(features_b,1)
    plot([rslt_b(i,1) pred_cnst_cross_mean(i,1)],[rslt_b(i,2) pred_cnst_cross_mean(i,2)],'y');    
end

for i=1:size(features_a,1)
    plot([rslt_a(i,1) pred_free_cross_mean(i,1)],[rslt_a(i,2) pred_free_cross_mean(i,2)],'y');
end

err_cnst = sum(((pred_cnst_mean(:,1:2) - rslt_a(:,1:2)).^2)');
err_cnst_cross = sum(((pred_free_cross_mean(:,1:2) - rslt_a(:,1:2)).^2)');

err_free = sum(((pred_free_mean(:,1:2) - rslt_b(:,1:2)).^2)');
err_free_cross = sum(((pred_cnst_cross_mean(:,1:2) - rslt_b(:,1:2)).^2)');

for i=1:size(err_cnst,2)
    if err_cnst(1,i) > err_cnst_cross(1,i)
        disp('say hello');
    end
end

for i=1:size(err_free,2)
    if err_free(1,i) > err_free_cross(1,i)
        disp('say hello');
    end
end

figure;
hold on;
plot(rslt_all(:,1),rslt_all(:,2),'.r');

for i=1:3
    [pred_all_mean(:,i) pred_all_sd(:,i)] = predict(gprModel_all{i},feat_all);    
end

plot(pred_all_mean(:,1),pred_all_mean(:,2),'or');

for i=1:size(feat_all,1)
    plot([rslt_all(i,1) pred_all_mean(i,1)],[rslt_all(i,2) pred_all_mean(i,2)],'y');
end


