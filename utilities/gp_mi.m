close all
clear all

[feat_all rslt_all] = load_test_dataset();
%[feat_train rslt_train feat_test rslt_test] = generate_train_test_dataset(feat_all,rslt_all,40);

%hyperparameter optimization
feat_prior = feat_all(1:5:end,:);
rslt_prior = rslt_all(1:5:end,:);

for i=1:3
    gprMdl_prior{i} = fitrgp(feat_prior,rslt_prior(:,i),'KernelFunction','squaredexponential');
    kernel_params(i,:) = gprMdl_prior{i}.KernelInformation.KernelParameters';
    [mean_prior(:,i) sd_prior(:,i)] = resubPredict(gprMdl_prior{i}); 
end

feat_train = feat_all(1:30,:);
rslt_train = rslt_all(1:30,:);
feat_test = feat_all(31:end,:);
rslt_test = rslt_all(31:end,:);


%% no sample selection(=sequential learning)
% init
for smp_cnt=1:size(feat_train,1)
    for i=1:3
        sigma0 = 0.2;
        kparams0 = [5, 0.1];
        %gprMdl_init{i} = fitrgp(feat_train(1:smp_cnt,:), rslt_train(1:smp_cnt,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
        %kparams0 = kernel_params(i,:)';
        gprMdl_init{i} = fitrgp(feat_train(1:smp_cnt,:), rslt_train(1:smp_cnt,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
        %[pred_init(:,i) pred_sd_init(:,i)] = predict(gprMdl_init{i},feat_test);
        loss_init(smp_cnt,i) = sqrt(loss(gprMdl_init{i},feat_test,rslt_test(:,i)));
        %gprMdl_init{i}.KernelInformation.KernelParameters
    end
    %rmse_init(smp_cnt,:) = rms(rslt_test-pred_init);    
end

% rnd_idx = randperm(size(feat_train,1));
% for i=1:size(feat_train,1)
%     feat_train_rnd(i,:) = feat_train(rnd_idx(i),:);
%     rslt_train_rnd(i,:) = rslt_train(rnd_idx(i),:);
% end
% feat_train = feat_train_rnd;
% rslt_train = rslt_train_rnd;
% 
% for smp_cnt=1:size(feat_train,1)
%     for i=1:3
%         sigma0 = 0.2;
%         kparams0 = [3.5, 6.2];
%         %gprMdl_init{i} = fitrgp(feat_train(1:smp_cnt,:), rslt_train(1:smp_cnt,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
%         %kparams0 = kernel_params(i,:)';
%         gprMdl_init{i} = fitrgp(feat_train(1:smp_cnt,:), rslt_train(1:smp_cnt,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
%         %[pred_init(:,i) pred_sd_init(:,i)] = predict(gprMdl_init{i},feat_test);
%         loss_rnd(smp_cnt,i) = sqrt(loss(gprMdl_init{i},feat_test,rslt_test(:,i)));
%         %gprMdl_init{i}.KernelInformation.KernelParameters
%     end
%     %rmse_init(smp_cnt,:) = rms(rslt_test-pred_init);    
% end

%% random sample selection
% rnd_idx = randperm(size(feat_train,1));
% for i=1:size(rnd_idx,2)
%     feat_train_rnd(i,:) = feat_train(rnd_idx(i),:);
%     rslt_train_rnd(i,:) = rslt_train(rnd_idx(i),:);
% end
% 
% for smp_cnt=1:size(feat_train,1)
%     for i=1:3
%         sigma0 = 0.2;
%         %kparams0 = [3.5, 6.2];
%         %gprMdl_init{i} = fitrgp(feat_train(1:smp_cnt,:), rslt_train(1:smp_cnt,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
%         kparams0 = kernel_params(i,:)';
%         gprMdl_init_rnd{i} = fitrgp(feat_train_rnd(1:smp_cnt,:), rslt_train_rnd(1:smp_cnt,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
%         %[pred_init_rnd(:,i) pred_sd_init(:,i)] = predict(gprMdl_init{i},feat_test);
%         loss_init_rnd(smp_cnt,i) = sqrt(loss(gprMdl_init_rnd{i},feat_test,rslt_test(:,i)));
%         %gprMdl_init{i}.KernelInformation.KernelParameters
%     end
%     %rmse_init(smp_cnt,:) = rms(rslt_test-pred_init);    
% end

%% Informative sample selection

kernel_params(1,:) = kparams0;
%kernel_params(1,2) = 0.1;
%sample selection
X_obsvd = [];
Y_obsvd = [];
X_unobsvd = feat_train;
Y_unobsvd = rslt_train;
for smp_cnt=1:size(feat_train,1)
    %select sample
    mi = [];
    for j=1:size(X_unobsvd,1)
        X_cand = X_unobsvd(j,:);
        if smp_cnt==1
            Y_var_obsvd = 1;
        else
            Y_var_obsvd = gp_varonly(X_cand,X_obsvd,kernel_params(1,1),kernel_params(1,2));
        end
        X_unobsvd_cand = X_unobsvd;
        X_unobsvd_cand(j,:) = [];
        if smp_cnt==size(feat_train,1)
            Y_var_unobsvd = 1;
        else
            Y_var_unobsvd = gp_varonly(X_cand, X_unobsvd_cand,kernel_params(1,1),kernel_params(1,2));
        end
        mi(j) = 0.5 * log( Y_var_obsvd / Y_var_unobsvd);
    end
    [max_mi idx] = max(mi);
    if smp_cnt ==1
        idx = 1;
    end
    %compute the model
    X_obsvd = [X_obsvd; X_unobsvd(idx,:)];
    Y_obsvd = [Y_obsvd; Y_unobsvd(idx,:)];
    for i=1:3
        sigma0 = 0.2;
        kparams0 = kernel_params(i,:)';
        gprMdl_AL{i} = fitrgp(X_obsvd,Y_obsvd(:,i),'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
        loss_AL(smp_cnt,i) = sqrt(loss(gprMdl_AL{i},feat_test,rslt_test(:,i)));
    end
    %update the obsvd/unobsvd sets
    X_unobsvd(idx,:) = [];
    Y_unobsvd(idx,:) = [];
end

figure;
subplot(3,1,1);
title('Sample selection');
plot(loss_init(:,1),'r');
hold on;
plot(loss_AL(:,1),'b');
ylabel('y');

subplot(3,1,2);
plot(loss_init(:,2),'r');
hold on;
plot(loss_AL(:,2),'b');
ylabel('z');

subplot(3,1,3);
plot(loss_init(:,3),'r');
hold on
plot(loss_AL(:,3),'b');
ylabel('theta');

xlabel('Number of samples to train model');
legend('Sequential','Affordance Learning','Location','NorthEast');

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