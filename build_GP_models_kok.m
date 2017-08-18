function build_GP_models_kok(dataFolder)
    addpath('utilities');
    %dataFolder = 'data/reactive/';

    bool_end = false;
    i=0;
    while ~bool_end 
        i = i+1;
        if(size(dir([dataFolder 'react_feat_n_result' num2str(i) '*']),1)==0)
            bool_end = true;
            scene_cnt = i -1;
        end
    end

    features_all = [];
    results_all = [];


    for i=1:scene_cnt
        filename = [dataFolder 'react_feat_n_result' num2str(i) '.mat'];
        load(filename);
        features_all = [features_all;features];
        results_all = [results_all;results];
    end

    %filter unavailable value
    [idx1 idx2] = find(isnan(features_all)==1);
    idx1 = unique(idx1);
    for i=size(idx1,1):-1:1
        features_all(idx1(i),:) = [];
        results_all(idx1(i),:) = [];
    end

%     features_all = features_all(1:10:end,:);
%     results_all = results_all(1:10:end,:);
    
    % GP with ARD exponential
    length_scale_init_val = 1;
    %sigma0 = std(rslt_train);
    sigma0 = [0.1 0.1 0.1 0.01 0.01 0.01 0.1];
    sigmaF0 = sigma0;
    dim = size(features_all,2);
    sigmaM0 = length_scale_init_val*ones(dim,1);
%     sigmaM0(end,1) = 0.1;

    for i=1:6
        disp(['Training on ' num2str(i) 'th model']);
        gprMdl{i} = fitrgp(features_all, results_all(:,i), 'Basis', 'constant', 'FitMethod','exact',...
            'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
            'KernelParameters',[sigmaM0;sigmaF0(1,i)],'Sigma',sigma0(1,i),'Standardize',1);
        %[pred(:,i) pred_sd(:,i)] = predict(gprMdl{i},feat_test);        
        pred_train(:,i) = resubPredict(gprMdl{i});
    end

    save([dataFolder 'GP_models_feature_ard_react.mat'],'gprMdl');
end

% mdl_WGP{i,j}.hyp = repmat(1.0,D+2,1);
% %mdl_WGP{i}.hyp = [3.7;3.6;-1.49;-4];
% 
% % compute subset of training if weight is small
% sub_input = [];
% sub_target = [];
% sub_weight = [];
% for k=1:N
%     if pred_label(k,j) > 1e-9 % weight is non-zero, failing inversion of W
%         sub_input = [sub_input;input(k,:)];
%         sub_target = [sub_target;target(k,i)];
%         sub_weight = [sub_weight;pred_label(k,j)];
%     end
% end
% 
% %             weight = pred_label(:,i);
% %             for j=1:N
% %                 if weight(j) ==0
% %                     weight(j) = 1e-6;
% %                 end
% %             end
% 
% % hyperparameter optimize
% %[hyp_new, fX, i] = minimize(hyp, 'gp01lik', 100, x, y);
% [mdl_WGP{i,j}.hyp, mdl_WGP{i,j}.nll, itr] = minimize(mdl_WGP{i,j}.hyp, 'wgp_lik', 10, sub_weight, sub_input, sub_target);
% mdl_WGP{i,j}.weight = sub_weight;
% mdl_WGP{i,j}.input = sub_input;
% mdl_WGP{i,j}.target = sub_target;
% 
% % Gaussian process
% [mdl_WGP{i,j}.pred_mu mdl_WGP{i,j}.pred_ss] = wgp_pred(mdl_WGP{i,j}.hyp, mdl_WGP{i,j}.weight, mdl_WGP{i,j}.input, mdl_WGP{i,j}.target, input);
% mdl_WGP{i,j}.rmse = rms(mdl_WGP{i,j}.pred_mu-target(:,i));