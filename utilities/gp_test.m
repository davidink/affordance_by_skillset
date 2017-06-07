close all
clear all

[feat_all rslt_all] = load_test_dataset();
feat_all = feat_all(70:end,:);
rslt_all = rslt_all(70:end,:);
[feat_train rslt_train feat_test rslt_test] = generate_train_test_dataset(feat_all,rslt_all,30);


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
    pred_prior(:,i) = resubPredict(gprMdl_prior{i});
end

% fig_prior = figure;
% hold on;
% plot(pred_prior(:,1),pred_prior(:,2),'xr');
% plot(rslt_prior(:,1),rslt_prior(:,2),'ob');


pred_all = [];
for i=1:3
    sigma0 = 0.1;
    kparams0 = [3.5, 6.2];
    %kparams0 = kernel_params(i,:)';
    %gprMdl_all{i} = fitrgp(feat_all,rslt_all(:,i),'KernelFunction','squaredexponential');
    gprMdl_all{i} = fitrgp(feat_all,rslt_all(:,i),'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
    pred_all(:,i) = resubPredict(gprMdl_all{i});
    loss_all(:,i) = resubLoss(gprMdl_all{i});
end

fig_pred = figure;
plot(pred_all(:,1),pred_all(:,2),'xr');
hold on;
plot(rslt_all(:,1),rslt_all(:,2),'ob');
num_sample = size(pred_all,1);
title(['Training Sample number: ' num2str(num_sample)]);
% init
for smp_cnt=1:size(feat_train,1)
    for i=1:3
        sigma0 = 0.2;
        %kparams0 = [3.5, 6.2];
        %gprMdl_init{i} = fitrgp(feat_train(1:smp_cnt,:), rslt_train(1:smp_cnt,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
        kparams0 = kernel_params(i,:)';
        %gprMdl_init{i} = fitrgp(feat_train(1:smp_cnt,:), rslt_train(1:smp_cnt,i), 'KernelFunction','squaredexponential','KernelParameters',kparams0,'Sigma',sigma0);
        crnt_feat_train = feat_train(1:smp_cnt,:);
        crnt_rslt_train = rslt_train(1:smp_cnt,:);
        gprMdl_init{i} = fitrgp(crnt_feat_train, crnt_rslt_train(:,i), 'KernelFunction','squaredexponential');
        [pred_init(:,i) pred_sd_init(:,i)] = predict(gprMdl_init{i},feat_test);
        %pred_train(:,i) = resubPredict(gprMdl_init{i});
        %[pred_train(:,i) pred_sd_train(:,i)] = predict(gprMdl_init{i},crnt_feat_train);
        loss_init(smp_cnt,i) = sqrt(loss(gprMdl_init{i},feat_test,rslt_test(:,i)));
        %gprMdl_init{i}.KernelInformation.KernelParameters
    end
    figure(fig_pred);
    clf(fig_pred);
    hold on;
    plot(rslt_test(:,1),rslt_test(:,2),'.r');
    plot(pred_init(:,1),pred_init(:,2),'xb');
    
    plot(crnt_rslt_train(:,1),crnt_rslt_train(:,2),'oc');
    plot(pred_train(:,1),pred_train(:,2),'xg');
    rmse_init(smp_cnt,:) = rms(rslt_test-pred_init);    
end

figure;
subplot(3,1,1);
title('No grouping');
plot(loss_init(:,1),'r');
hold on;
%plot(RMSE_rnd_n(:,1),'b');
ylabel('y');

subplot(3,1,2);
plot(loss_init(:,2),'r');
hold on;
%plot(RMSE_rnd_n(:,2),'b');
ylabel('z');

subplot(3,1,3);
plot(loss_init(:,3),'r');
hold on
%plot(RMSE_rnd_n(:,3),'b');
ylabel('theta');

xlabel('Number of samples to train model');
legend('GP','BLR','Location','NorthEast');