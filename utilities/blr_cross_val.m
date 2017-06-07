close all
clear all

[feat_all rslt_all] = load_test_dataset();

% randomly sort the features
rnd_idx = randperm(size(feat_all,1),size(feat_all,1));

for i=1:size(rnd_idx,2);
    feat_rnd(i,:) = feat_all(rnd_idx(1,i),:);
    rslt_rnd(i,:) = rslt_all(rnd_idx(1,i),:);
end
feat_all = feat_rnd;
rlst_all = rslt_rnd;

% cross validation
k = 5;
fold_size = floor(size(feat_all,1)/k);

%figure;
hold on;
pred_all = [];
for f=1:k
    disp(['Computing for fold ' num2str(f)]);
    %divide train/test
    if f ~= k
        feat_test = feat_all((f-1)*fold_size+1:f*fold_size,:);
        rslt_test = rslt_all((f-1)*fold_size+1:f*fold_size,:);
        feat_train = feat_all;
        rslt_train = rslt_all;
        feat_train((f-1)*fold_size+1:f*fold_size,:) = [];
        rslt_train((f-1)*fold_size+1:f*fold_size,:) = [];
    else
        feat_test = feat_all((f-1)*fold_size+1:end,:);
        rslt_test = rslt_all((f-1)*fold_size+1:end,:);
        feat_train = feat_all;
        rslt_train = rslt_all;
        feat_train((f-1)*fold_size+1:end,:) = [];
        rslt_train((f-1)*fold_size+1:end,:) = [];
    end
    
    % build model
    pred=[];
    pred_sd=[];
    pred_train = [];
    %% GP model
%     % GP with ARD exponential
%     length_scale_init_val = 10;
%     %sigma0 = std(rslt_train);
%     sigma0 = [0.01 0.01 0.05];
%     sigmaF0 = sigma0;
%     dim = size(feat_train,2);
%     sigmaM0 = length_scale_init_val*ones(dim,1);
%     for i=1:3
% %         sigma0 = 0.2;
% %         kparams0 = [3.5, 6.2];
% %         %gprMdl_init{i} = fitrgp(feat_train(1:smp_cnt,:), rslt_train(1:smp_cnt,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
% %         %kparams0 = kernel_params(i,:)';
% %         gprMdl{i} = fitrgp(feat_train, rslt_train(:,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
%         gprMdl{i} = fitrgp(feat_train, rslt_train(:,i), 'Basis', 'constant', 'FitMethod','exact',...
%             'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
%             'KernelParameters',[sigmaM0;sigmaF0(1,i)],'Sigma',sigma0(1,i),'Standardize',1);
%         [pred(:,i) pred_sd(:,i)] = predict(gprMdl{i},feat_test);        
%         pred_train(:,i) = resubPredict(gprMdl{i});
%     end
    %% Bayesian Linear Regression
    %No grouping
    %[feat_train rslt_train feat_test rslt_test] = generate_train_test_dataset(feat_all,rslt_all,15);
    alpha_inv = 0.0001; % esitmated variance of w distribution
    beta_inv = 0.0001; % estimated noise variance of y distribution
    [mu{1} sigsq{1}] = bayesian_regression(feat_train,rslt_train(:,1),alpha_inv,beta_inv);
    [mu{2} sigsq{2}] = bayesian_regression(feat_train,rslt_train(:,2),alpha_inv,beta_inv);
    alpha_inv = 0.01;
    beta_inv = 0.01;
    [mu{3} sigsq{3}] = bayesian_regression(feat_train,rslt_train(:,3),alpha_inv,beta_inv);

    for itr=1:50
        for i=1:3
            w{i} = mvnrnd(mu{i},sigsq{i});
            E_Y{i} = w{i} * feat_test';
        end
        E_Yt = [E_Y{1}' E_Y{2}' E_Y{3}'];
        diff = rslt_test - E_Yt;
        RMSE(itr,:) = rms(diff);
    end
    
    rmse(f,:) = mean(RMSE);
%     plot(rslt_train(:,1),rslt_train(:,2),'or');
%     plot(pred(:,1),pred(:,2),'xb');    
    %pred_all = [pred_all;pred];    
end

figure;
plot(pred_all(:,1),pred_all(:,2),'xb');
hold on;
plot(rslt_all(:,1),rslt_all(:,2),'or');

diff = pred_all - rslt_all;
rmse_all = rms(diff);

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

