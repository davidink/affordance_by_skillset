close all
clear all

[feat_all rslt_all] = load_test_dataset();
%feat_all = feat_all(:,1:56);
[feat_train rslt_train feat_test rslt_test] = generate_train_test_dataset(feat_all,rslt_all,80);


% feat_train = feat_all(1:30,:);
% rslt_train = rslt_all(1:30,:);
% feat_test = feat_all(31:end,:);
% rslt_test = rslt_all(31:end,:);

%no init
% for i=1:3
%     gprMdl_noinit{i} = fitrgp(feat_train, rslt_train(:,i),'KernelFunction','squaredexponential');
%     [pred_noinit(:,i) pred_sd_noinit(:,i)] = predict(gprMdl_noinit{i},feat_test);
% end
% rmse_noinit = rms(rslt_test-pred_noinit);

%hyperparameter optimization
feat_prior = feat_all(1:10:end,:);
rslt_prior = rslt_all(1:10:end,:);
for i=1:3
    gprMdl_prior{i} = fitrgp(feat_prior,rslt_prior(:,i),'KernelFunction','squaredexponential');
    kernel_params(i,:) = gprMdl_prior{i}.KernelInformation.KernelParameters';
end

% init
for smp_cnt=1:size(feat_train,1)
    for i=1:3
        sigma0 = 0.2;
        kparams0 = [3.5, 6.2];
        %gprMdl_init{i} = fitrgp(feat_train(1:smp_cnt,:), rslt_train(1:smp_cnt,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
        %kparams0 = kernel_params(i,:)';
        gprMdl_init{i} = fitrgp(feat_train(1:smp_cnt,:), rslt_train(1:smp_cnt,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
        [pred_init(:,i) pred_sd_init(:,i)] = predict(gprMdl_init{i},feat_test);
        loss_init(smp_cnt,i) = sqrt(loss(gprMdl_init{i},feat_test,rslt_test(:,i)));
        %gprMdl_init{i}.KernelInformation.KernelParameters
    end
    rmse_init(smp_cnt,:) = rms(rslt_test-pred_init);    
end

%No grouping
%[feat_train rslt_train feat_test rslt_test] = generate_train_test_dataset(feat_all,rslt_all,15);
alpha_inv = 0.00001; % esitmated variance of w distribution
beta_inv = 0.00001; % estimated noise variance of y distribution
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
            var_y{i}=[];
            for k=1:size(feat_test,1)
                var_y{i}(k) = beta_inv + feat_test(k,:) * sigsq{i}{step} * feat_test(k,:)';
            end
            E_Yd{i} = mvnrnd(mu_y{i},var_y{i});            
        end
        E_Yt = [E_Y{1}' E_Y{2}' E_Y{3}'];
%         E_Ytd = [E_Yd{1}' E_Yd{2}' E_Yd{3}'];
        diff = rslt_test - E_Yt;
        for j=1:size(diff,1)
            if diff(j,3) < -3.141592/2
                diff(j,3) = diff(j,3) + 3.141592/2;
            else if diff(j,3) > 3.141592/2
                    diff(j,3) = diff(j,3) - 3.141592/2;
                end
            end
        end
%         diff_d = rslt_test - E_Ytd;
%         for j=1:size(diff_d,1)
%             if diff_d(j,3) < -3.141592/2
%                 diff_d(j,3) = diff_d(j,3) + 3.141592/2;
%             else if diff_d(j,3) > 3.141592/2
%                     diff_d(j,3) = diff_d(j,3) - 3.141592/2;
%                 end
%             end
%         end
        RMSE(itr,:) = rms(diff);
%         RMSE_d(itr,:) = rms(diff_d);
    end    
    RMSE_rnd_n(step,:) = mean(RMSE);
%     RMSE_rnd_nd(step,:) = mean(RMSE_d);
end

figure;
subplot(3,1,1);
title('No grouping');
plot(loss_init(:,1),'r');
hold on;
plot(RMSE_rnd_n(:,1),'b');
ylabel('y');

subplot(3,1,2);
plot(loss_init(:,2),'r');
hold on;
plot(RMSE_rnd_n(:,2),'b');
ylabel('z');

subplot(3,1,3);
plot(loss_init(:,3),'r');
hold on
plot(RMSE_rnd_n(:,3),'b');
ylabel('theta');

xlabel('Number of samples to train model');
legend('GP','BLR','Location','NorthEast');

% %ard
% for i=1:3
%     sigma0 = std(rslt_train);
%     sigmaF0 = sigma0;
%     d = size(feat_train,2);
%     sigmaM0 = 10*ones(d,1);
%     gprMdl_ard{i} = fitrgp(feat_train,rslt_train(:,i),'Basis','constant','FitMethod','exact',...
%     'Predictmethod','exact','KernelFunction','ardsquaredexponential',...
%     'KernelParameters',[sigmaM0;sigmaF0(i)],'Sigma',sigma0(i),'Standardize',1);
%     pred_ard(:,i) = predict(gprMdl_ard{i},feat_test);
% end
% rmse_ard = rms(rslt_test-pred_ard)

% %cross-val
% for i=1:3
%     gprMdl_cv{i} = fitrgp(feat_all, rslt_all(:,i),'Holdout',0.25);
%     pred_cv(:,i) = kfoldPredict(gprMdl_cv{i});
% end
% rmse_cv = rms(rslt_all-pred_cv)

%kL = kfoldLoss(cvgprMdl_y)
%L = loss(gprMdl_y,feat_test,rslt_test(:,1))
 
% pred_y = resubPredict(gprMdl_y);
% pred_z = resubPredict(gprMdl_z);
% pred_theta = resubPredict(gprMdl_theta);
% [pred_y pred_var_y] = predict(gprMdl_y,feat_test);
% [pred_z pred_var_z] = predict(gprMdl_z,feat_test);
% [pred_theta pred_var_theta] = predict(gprMdl_theta,feat_test);

% pred_y_all = predict(gprMdl_y_all,feat_test);
% pred_z_all = predict(gprMdl_z_all,feat_test);
% pred_theta_all = predict(gprMdl_theta_all,feat_test);
% pred_y_all = kfoldPredict(cvgprMdl_y);
% pred_z_all = kfoldPredict(cvgprMdl_z);
% pred_theta_all = kfoldPredict(cvgprMdl_theta);
% 
% diff = [rslt_test(:,1) - pred_y rslt_test(:,2) - pred_z rslt_test(:,3)-pred_theta];
% diff_all = [rslt_test(:,1)-pred_y_all rslt_test(:,2)-pred_z_all rslt_test(:,3)-pred_theta_all];
% rmse = rms(diff)
% rmse_all = rms(diff_all)
% 
% figure;
% plot(rslt_test(:,1),rslt_test(:,2),'ob');
% hold on
% plot(pred_y,pred_z,'xr');
% plot(pred_y_all,pred_z_all,'*g');
% 
% figure;
% plot(rslt_test(:,3), pred_theta, '*r');