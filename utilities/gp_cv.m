function [gp_model rmse_mean rmse_std] = gp_cv(features,results,k)
    
    feat_all = features;
    rslt_all = results;

    % randomly sort the features
    rnd_idx = randperm(size(feat_all,1),size(feat_all,1));
    for i=1:size(rnd_idx,2);
        feat_rnd(i,:) = feat_all(rnd_idx(1,i),:);
        rslt_rnd(i,:) = rslt_all(rnd_idx(1,i),:);
    end
    
    feat_all = feat_rnd;
    rlst_all = rslt_rnd;

    % cross validation
    fold_size = floor(size(feat_all,1)/k);

%     figure;
%     hold on;
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
        
        % GP with ARD exponential
        length_scale_init_val = 1;
        sigma0 = std(rslt_train);
        %sigma0 = [0.01 0.01 0.05];
        sigmaF0 = sigma0;
        dim = size(feat_train,2);
        sigmaM0 = length_scale_init_val*ones(dim,1);
        for i=1:3
            gprMdl{i} = fitrgp(feat_train, rslt_train(:,i), 'Basis', 'constant', 'FitMethod','exact',...
                'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
                'KernelParameters',[sigmaM0;sigmaF0(1,i)],'Sigma',sigma0(1,i),'Standardize',1);
            [pred(:,i) pred_sd(:,i)] = predict(gprMdl{i},feat_test);        
            pred_train(:,i) = resubPredict(gprMdl{i});
        end
        rmse(f,:) = rms(pred - rslt_test);
        rmse_train(f,:) = rms(pred_train -rslt_train);
        pred_all = [pred_all;pred];
        
%         figure;
%         hold on;
%         plot(rslt_all(:,1),rslt_all(:,2),'.r');
%         plot(pred_train(:,1),pred_train(:,2),'xb');
%         plot(pred(:,1),pred(:,2),'ob');
%         for l=1:size(rslt_train,1)
%             plot([rslt_train(l,1) pred_train(l,1)],[rslt_train(l,2) pred_train(l,2)],'g');
%         end
%         for l=1:size(rslt_test,1)
%             plot([rslt_test(l,1) pred(l,1)],[rslt_test(l,2) pred(l,2)],'c');
%         end
%         text(-0.05,0,['rmse\_test:  ' num2str(rmse(f,:))]);
%         text(-0.05,-0.005,['rmse\_train:  ' num2str(rmse_train(f,:))]);
            
    end
    
    gp_model = gprMdl;
    rmse_mean = mean(rmse);
    rmse_std = std(rmse);

end