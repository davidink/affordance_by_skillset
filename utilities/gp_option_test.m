close all
clear all

[feat_all rslt_all] = load_test_dataset();
%feat_all = feat_all(:,1:56);
%feat_all = feat_all(60:end,:);
%rslt_all = rslt_all(60:end,:);
[feat_train rslt_train feat_test rslt_test] = generate_train_test_dataset(feat_all,rslt_all,80);

% feat_train = feat_all(1:80,:);
% rslt_train = rslt_all(1:80,:);
% feat_test = feat_all(81:end,:);
% rslt_test = rslt_all(81:end,:);

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

% cross validation model building
% for i=1:3
%     cvgprMdl{i} = fitrgp(feat_all,rslt_all(:,i),'Standardize',1,'Holdout',0.25);
%     loss(1,i) = kfoldLoss(cvgprMdl{i});
%     ypred = kfoldPredict(cvgprMdl{i});
% end

% GP with ARD exponential
length_scale_init_val = 10;
%sigma0 = std(rslt_train);
sigma0 = [0.01 0.01 0.2];
sigmaF0 = sigma0;
dim = size(feat_train,2);
sigmaM0 = length_scale_init_val*ones(dim,1);

for i=1:3
    kparams0 = kernel_params(i,:)';
    %kparams0 = [3.5, 6.2];
    
    %choose model option for gp
    %sigma0 = 0.05;
    %gprMdl{i} = fitrgp(feat_train, rslt_train(:,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
    %gprMdl{i} = fitrgp(feat_train, rslt_train(:,i), 'KernelFunction','squaredexponential');
    % GP with ARD exponential
%     gprMdl{i} = fitrgp(feat_train, rslt_train(:,i), 'Basis', 'constant', 'FitMethod','exact',...
%         'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
%         'KernelParameters',[sigmaM0;sigmaF0(1,i)],'Sigma',sigma0(1,i),'Standardize',1);
    gprMdl{i} = fitrgp(feat_all, rslt_all(:,i), 'Basis', 'constant', 'FitMethod','exact',...
        'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
        'KernelParameters',[sigmaM0;sigmaF0(1,i)],'Sigma',sigma0(1,i),'Standardize',1);
   
    %cvgprMdl{i} = fitrgp(feat_all,rslt_all,'Standardize',1,'Holdout',0.25);

    pred_train(:,i) = resubPredict(gprMdl{i});
    [pred_gp(:,i) pred_gp_sd(:,i)] = predict(gprMdl{i},feat_test);
    %loss(i) = sqrt(loss(gprMdl{i},feat_test,rslt_test(:,i)));
    %gprMdl_init{i}.KernelInformation.KernelParameters
end
loss_train = rms(pred_train -rslt_train)
loss_test = rms(pred_gp -rslt_test)
%fig_plot = figure('Position', [1000, 300, 1000, 500]);
figure1 = figure;
%plot(rslt_train(:,1),rslt_train(:,2),'or');
plot(rslt_all(:,1),rslt_all(:,2),'or');
hold on;
plot(pred_train(:,1),pred_train(:,2),'xb');
title('Train');
axis equal;
figure2 = figure('Position',[1243 678 560 420]);
plot(rslt_test(:,1),rslt_test(:,2),'or');
hold on
plot(pred_gp(:,1),pred_gp(:,2),'xb');
title('Test');
axis equal;

rmse_init(smp_cnt,:) = rms(rslt_test-pred_init);
for smp_cnt=1:size(feat_train,1)
    for i=1:3
        sigma0 = 0.2;
        %kparams0 = [3.5, 6.2];
        %gprMdl_init{i} = fitrgp(feat_train(1:smp_cnt,:), rslt_train(1:smp_cnt,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
        kparams0 = kernel_params(i,:)';
        gprMdl_init{i} = fitrgp(feat_train(1:smp_cnt,:), rslt_train(1:smp_cnt,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
        [pred_init(:,i) pred_sd_init(:,i)] = predict(gprMdl_init{i},feat_test);
        loss_init(smp_cnt,i) = sqrt(loss(gprMdl_init{i},feat_test,rslt_test(:,i)));
        %gprMdl_init{i}.KernelInformation.KernelParameters
    end
    rmse_init(smp_cnt,:) = rms(rslt_test-pred_init);    
end