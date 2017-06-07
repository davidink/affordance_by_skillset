close all
clear all

[feat_all rslt_all] = load_test_dataset();

%%
%generating data sets
[feat_train rslt_train feat_test rslt_test] = generate_train_test_dataset(feat_all,rslt_all,30);
alpha_inv = 0.0001; % esitmated variance of w distribution
beta_inv = 0.0001; % estimated noise variance of y distribution
[mu{1} sigsq{1}] = bayesian_regression_update(feat_train,rslt_train(:,1),alpha_inv,beta_inv);
[mu{2} sigsq{2}] = bayesian_regression_update(feat_train,rslt_train(:,2),alpha_inv,beta_inv);
alpha_inv = 0.01;
beta_inv = 0.01;
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

%%
%sort the input features to be more informative

feat_obsvd = []; rslt_obsvd = [];
feat_unobsvd = feat_train; rslt_unobsvd = rslt_train;
iNumRegressors = size(feat_train,2);
for i=1:size(feat_train,1)
    alpha_inv = 0.0001;
    beta_inv = 0.01;
    
    if size(feat_obsvd,1) == 0; % starting
        w_var_obsvd = inv(alpha_inv * eye(iNumRegressors));
    else
        w_var_obsvd = bayesian_regression_varonly(feat_obsvd, alpha_inv, beta_inv);
    end
    
    MI = [];
    ent = [];
    for j=1:size(feat_unobsvd,1)
        feat_cand = feat_unobsvd(j,:);
        feat_ex_cand = feat_unobsvd;
        feat_ex_cand(j,:) = [];
        y_var_obsvd = beta_inv + feat_cand * w_var_obsvd * feat_cand';
        w_var_unobsvd = bayesian_regression_varonly(feat_ex_cand, alpha_inv, beta_inv);
        y_var_unobsvd = beta_inv + feat_cand * w_var_unobsvd * feat_cand';
        MI(j) = 0.5 * log(y_var_obsvd/y_var_unobsvd);
        ent(j) = log(y_var_obsvd);
    end
    %[max_mi idx] = max(MI);
    [max_ent idx] = min(ent);
    feat_obsvd = [feat_obsvd; feat_unobsvd(idx,:)];
    rslt_obsvd = [rslt_obsvd; rslt_unobsvd(idx,:)];
    feat_unobsvd(idx,:) = [];
    rslt_unobsvd(idx,:) = [];
end

feat_train = feat_obsvd;
rslt_train = rslt_obsvd;
%[feat_train rslt_train feat_test rslt_test] = generate_train_test_dataset(feat_all,rslt_all,30);
alpha_inv = 0.0001; % esitmated variance of w distribution
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
    RMSE_rnd_b(step,:) = mean(RMSE);
    RMSE_rnd_bd(step,:) = mean(RMSE_d);
    RMSE_rnd_bm(step,:) = mean(RMSE_m);
end

% %%
% %ploting
% 
% 
figure;
subplot(3,1,1);
plot(RMSE_rnd_a(:,1),'-r');
hold on
plot(RMSE_rnd_b(:,1),'-b');
% hold on
% plot(RMSE_rnd_am(:,1),'-g');
subplot(3,1,2);
plot(RMSE_rnd_a(:,2),'-r');
hold on
plot(RMSE_rnd_b(:,2),'-b');
% hold on
% plot(RMSE_rnd_am(:,2),'-g');
subplot(3,1,3);
plot(RMSE_rnd_a(:,3),'-r');
hold on
plot(RMSE_rnd_b(:,3),'-b');
% hold on
% plot(RMSE_rnd_am(:,3),'-g');

figure;
plot(rslt_test(:,3),E_Yt(:,3),'xr')
xlabel('Ground Truth');
ylabel('Expected output theta');


