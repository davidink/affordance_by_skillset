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

figure;
%%
%group_b with env cons
%generating data sets
%[feat_train rslt_train feat_test rslt_test] = generate_train_test_dataset(features_a,group_a,15);
[feat_train rslt_train feat_test rslt_test] = generate_train_test_dataset(features_b,group_b,15);
alpha_inv = 0.01; % esitmated variance of w distribution
beta_inv = 0.0001; % estimated noise variance of y distribution
[mu{1} sigsq{1}] = bayesian_regression_update(feat_train,rslt_train(:,1),alpha_inv,beta_inv);
[mu{2} sigsq{2}] = bayesian_regression_update(feat_train,rslt_train(:,2),alpha_inv,beta_inv);
alpha_inv = 0.01;
beta_inv = 0.1;
[mu{3} sigsq{3}] = bayesian_regression_update(feat_train,rslt_train(:,3),alpha_inv,beta_inv);

for step=1:size(mu{1},2)
    for itr=1:50
        for i=1:3
            w{i} = mvnrnd(mu{i}{step},sigsq{i}{step});
            E_Y{i} = w{i} * feat_test';
            mu_y{i} = mu{i}{step}' * feat_test';
            for k=1:size(feat_test,1)
                var_y{i}(k) = beta_inv + feat_test(k,:) * sigsq{i}{step} * feat_test(k,:)';
            end
            E_Yd{i} = mvnrnd(mu_y{i},var_y{i}); 
            E_Ym{i} = mu_y{i};
        end
        E_Yt = [E_Y{1}' E_Y{2}' E_Y{3}'];
        E_Ytd = [E_Yd{1}' E_Yd{2}' E_Yd{3}'];
        E_Ytm = [E_Ym{1}' E_Ym{2}' E_Ym{3}'];
%         plot(rslt_test(:,1),rslt_test(:,2),'or');
%         hold on;
%         plot(E_Yt(:,1),E_Yt(:,2),'xb');
        diff = rslt_test - E_Yt;
        for j=1:size(diff,1)
            if diff(j,3) < -3.141592/2
                diff(j,3) = diff(j,3) + 3.141592/2;
            else if diff(j,3) > 3.141592/2
                    diff(j,3) = diff(j,3) - 3.141592/2;
                end
            end
        end
        diff_d = rslt_test - E_Ytd;
        for j=1:size(diff_d,1)
            if diff_d(j,3) < -3.141592/2
                diff_d(j,3) = diff_d(j,3) + 3.141592/2;
            else if diff_d(j,3) > 3.141592/2
                    diff_d(j,3) = diff_d(j,3) - 3.141592/2;
                end
            end
        end
        diff_m = rslt_test - E_Ytm;
        for j=1:size(diff_m,1)
            if diff_m(j,3) < -3.141592/2
                diff_m(j,3) = diff_m(j,3) + 3.141592/2;
            else if diff_m(j,3) > 3.141592/2
                    diff_m(j,3) = diff_m(j,3) - 3.141592/2;
                end
            end
        end
        RMSE(itr,:) = rms(diff);
        RMSE_d(itr,:) = rms(diff_d);
        RMSE_m(itr,:) = rms(diff_m);
    end    
    RMSE_rnd_a(step,:) = mean(RMSE);
    RMSE_rnd_ad(step,:) = mean(RMSE_d);
    RMSE_rnd_am(step,:) = mean(RMSE_m);
end

%group_b_without_environmental constraint
[feat_train rslt_train feat_test rslt_test] = generate_train_test_dataset(features_b,group_b,15);
feat_train = feat_train(:,1:13);
feat_test = feat_test(:,1:13);
alpha_inv = 0.01; % esitmated variance of w distribution
beta_inv = 0.0001; % estimated noise variance of y distribution
[mu{1} sigsq{1}] = bayesian_regression_update(feat_train,rslt_train(:,1),alpha_inv,beta_inv);
[mu{2} sigsq{2}] = bayesian_regression_update(feat_train,rslt_train(:,2),alpha_inv,beta_inv);
alpha_inv = 0.01;
beta_inv = 0.1;
[mu{3} sigsq{3}] = bayesian_regression_update(feat_train,rslt_train(:,3),alpha_inv,beta_inv);

for step=1:size(mu{1},2)
    for itr=1:50
        for i=1:3
            w{i} = mvnrnd(mu{i}{step},sigsq{i}{step});
            E_Y{i} = w{i} * feat_test';
            mu_y{i} = mu{i}{step}' * feat_test';
            var_y{i}=[];
            for k=1:size(feat_test,1)
                var_y{i}(k) = beta_inv + feat_test(k,:) * sigsq{i}{step} * feat_test(k,:)';
            end
            E_Yd{i} = mvnrnd(mu_y{i},var_y{i});            
        end
        E_Yt = [E_Y{1}' E_Y{2}' E_Y{3}'];
        E_Ytd = [E_Yd{1}' E_Yd{2}' E_Yd{3}'];
        plot(rslt_test(:,1),rslt_test(:,2),'or');
        hold on;
        plot(E_Yt(:,1),E_Yt(:,2),'xc');
        diff = rslt_test - E_Yt;
        for j=1:size(diff,1)
            if diff(j,3) < -3.141592/2
                diff(j,3) = diff(j,3) + 3.141592/2;
            else if diff(j,3) > 3.141592/2
                    diff(j,3) = diff(j,3) - 3.141592/2;
                end
            end
        end
        diff_d = rslt_test - E_Ytd;
        for j=1:size(diff_d,1)
            if diff_d(j,3) < -3.141592/2
                diff_d(j,3) = diff_d(j,3) + 3.141592/2;
            else if diff_d(j,3) > 3.141592/2
                    diff_d(j,3) = diff_d(j,3) - 3.141592/2;
                end
            end
        end
        RMSE(itr,:) = rms(diff);
        RMSE_d(itr,:) = rms(diff_d);
    end    
    RMSE_rnd_b(step,:) = mean(RMSE);
    RMSE_rnd_bd(step,:) = mean(RMSE_d);
end

% %No grouping
% [feat_train rslt_train feat_test rslt_test] = generate_train_test_dataset(feat_all,rslt_all,15);
% alpha_inv = 0.001; % esitmated variance of w distribution
% beta_inv = 0.001; % estimated noise variance of y distribution
% [mu{1} sigsq{1}] = bayesian_regression_update(feat_train,rslt_train(:,1),alpha_inv,beta_inv);
% [mu{2} sigsq{2}] = bayesian_regression_update(feat_train,rslt_train(:,2),alpha_inv,beta_inv);
% alpha_inv = 0.01;
% beta_inv = 0.01;
% [mu{3} sigsq{3}] = bayesian_regression_update(feat_train,rslt_train(:,3),alpha_inv,beta_inv);
% 
% for step=1:size(mu{1},2)
%     for itr=1:50
%         for i=1:3
%             w{i} = mvnrnd(mu{i}{step},sigsq{i}{step});
%             E_Y{i} = w{i} * feat_test';
%             mu_y{i} = mu{i}{step}' * feat_test';
%             var_y{i}=[];
%             for k=1:size(feat_test,1)
%                 var_y{i}(k) = beta_inv + feat_test(k,:) * sigsq{i}{step} * feat_test(k,:)';
%             end
%             E_Yd{i} = mvnrnd(mu_y{i},var_y{i});            
%         end
%         E_Yt = [E_Y{1}' E_Y{2}' E_Y{3}'];
%         E_Ytd = [E_Yd{1}' E_Yd{2}' E_Yd{3}'];
%         diff = rslt_test - E_Yt;
%         for j=1:size(diff,1)
%             if diff(j,3) < -3.141592/2
%                 diff(j,3) = diff(j,3) + 3.141592/2;
%             else if diff(j,3) > 3.141592/2
%                     diff(j,3) = diff(j,3) - 3.141592/2;
%                 end
%             end
%         end
%         diff_d = rslt_test - E_Ytd;
%         for j=1:size(diff_d,1)
%             if diff_d(j,3) < -3.141592/2
%                 diff_d(j,3) = diff_d(j,3) + 3.141592/2;
%             else if diff_d(j,3) > 3.141592/2
%                     diff_d(j,3) = diff_d(j,3) - 3.141592/2;
%                 end
%             end
%         end
%         RMSE(itr,:) = rms(diff);
%         RMSE_d(itr,:) = rms(diff_d);
%     end    
%     RMSE_rnd_n(step,:) = mean(RMSE);
%     RMSE_rnd_nd(step,:) = mean(RMSE_d);
% end


% %%
% %ploting
% 
% 
figure;
subplot(3,1,1);
plot(RMSE_rnd_a(:,1),'-r');
hold on
plot(RMSE_rnd_b(:,1),'-b');
ylabel('y');
legend('EC','withoutEC','Location','NorthEast');
% hold on
% plot(RMSE_rnd_am(:,1),'-g');
subplot(3,1,2);
plot(RMSE_rnd_a(:,2),'-r');
hold on
plot(RMSE_rnd_b(:,2),'-b');
ylabel('z');
legend('EC','withoutEC','Location','NorthEast');
% hold on
% plot(RMSE_rnd_am(:,2),'-g');
subplot(3,1,3);
plot(RMSE_rnd_a(:,3),'-r');
hold on
plot(RMSE_rnd_b(:,3),'-b');
ylabel('theta');
legend('EC','withoutEC','Location','NorthEast');
% hold on
% plot(RMSE_rnd_am(:,3),'-g');



