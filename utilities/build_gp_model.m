function gprModel = build_gp_model(features,results)

    % GP with ARD exponential
    length_scale_init_val = 10;
    %sigma0 = std(rslt_train);
    sigma0 = [0.01 0.01 0.05];
    sigmaF0 = sigma0;
    dim = size(features,2);
    sigmaM0 = length_scale_init_val*ones(dim,1);
    
    for i=1:3
        gprMdl{i} = fitrgp(features, results(:,i), 'Basis', 'constant', 'FitMethod','exact',...
            'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
            'KernelParameters',[sigmaM0;sigmaF0(1,i)],'Sigma',sigma0(1,i),'Standardize',1);
        %[pred(:,i) pred_sd(:,i)] = predict(gprMdl{i},feat_test);        
        %pred_train(:,i) = resubPredict(gprMdl{i});
    end
    %rmse(f,:) = rms(pred - rslt_test);
    %pred_all = [pred_all;pred];
    gprModel = gprMdl;
end
