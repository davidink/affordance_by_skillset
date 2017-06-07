function [mdl_WLR mdl_WGP] = computeEM_ind(num_submodel, input_WLR, input_WGP, target, init_label)
    
    N = size(input_WLR,1); % number of training examples
    D_LR = size(input_WLR,2); % size of input(feature) vector
    D_GP = size(input_WGP,2); % size of input(feature) vector
    S = size(target,2);
    M = num_submodel;
    
    EM_done = false; EM_itr = 0;    
    
    if nargin == 5 % given initialization
        run_init = false;
    else
        run_init = true;
    end
    
    while ~EM_done
        EM_itr = EM_itr +1;
        disp([num2str(EM_itr) ' round EM iteration..']);
    
        % E step - compute conditional probabilities p(m_i = j | s_i, a_i, s_i')
        if run_init == true && EM_itr ==1
            % Use GMM for initial submodels assignment
            options = statset('Display','final');
            gmm_init = false;
            while ~gmm_init
                try
                    gmmodel = fitgmdist(target,M,'Options',options);
                    gmm_init = true;
                catch ME
                    gmm_init = false;
                end
            end
            submodel_label = cluster(gmmodel, target);
            target_label = zeros(N,M);
            for i=1:N
                target_label(i,submodel_label(i)) = 1;
            end
            pred_label = target_label;
            
            % true initialization
%             target_label = zeros(N,M);
%             for i=1:N
%                 if input(i,1) > input(i,2)-0.2
%                 %if input(i,1) > input(i,2)
%                     target_label(i,1) = 1;
%                 else
%                     target_label(i,2) = 1;
%                 end
%             end
            
            % random initialization
%             target_label = zeros(N,M);
%             for i=1:N
%                 if input(i,1) > 0.0
%                     target_label(i,1) = 1;
%                 else
%                     target_label(i,2) = 1;
%                 end
%             end           
            

        else if run_init == false && EM_itr ==1
                initial_label = zeros(N,M);
                for i=1:N
                    initial_label(i,init_label(i)) = 1;
                end 
                % Do not too much emphasize on initialization
%                 for i=1:N
%                     for j=1:M
%                         if target_label(i,j) == 1
%                             target_label(i,j) = 0.8;
%                         else
%                             target_label(i,j) = 0.2/(M-1);
%                         end
%                     end
%                 end
                pred_label = initial_label;
            else
                %pred_label = zeros(N,M); % p(m=j|s,a,s')
                sample_likhd = ones(N,M);
                for i=1:S % per axis
                    for j=1:M % per submodel
                        for k=1:N % per sample
                            sample_likhd(k,j) = sample_likhd(k,j) * mdl_WGP{i,j}.lik_target(k,1);
                        end
                    end                     
                end
                pred_label = sample_likhd.*submodel_pred./repmat(sum(sample_likhd.*submodel_pred,2),1,M);
                for i=1:N
                    for j=1:M
                        if pred_label(i,j) < 1e-6
                            pred_label(i,j) = 0;
                        end
                    end
                end
            end
        end
        
%         figure;
%         plot3(input(:,1),input(:,2),pred_label(:,1),'r.');
%         hold on;
%         plot3(input(:,1),input(:,2),pred_label(:,2),'g.');        
%         pause(0.1);
        
        % M step
       
        % weighted logistic regression        
        if EM_itr ==1
            mdl_WLR = LogisticRegression(D_LR,M,'linearFeatures');
        else
            %buffer
            if isnan(mdl_WLR.W(1,1))
                suspicious = 1;
            end
            mdl_WLR_prev = mdl_WLR;
            mdl_WLR_prev_weight = mdl_WLR.W;
            mdl_WGP_prev = mdl_WGP;
        end
        WLR_success = mdl_WLR.Train(input_WLR,pred_label,1e-6);
%         if isnan(mdl_WLR.W(1,1))
%             suspicious = 1;
%         end
        if WLR_success == false;
            if EM_itr ==1 % bad GMM initialization
                options = statset('Display','final');
                gmm_init = false;
                while ~gmm_init
                    try
                        gmmodel = fitgmdist(target(:,4:6),M,'Options',options);
                        gmm_init = true;
                    catch ME
                        gmm_init = false;
                    end
                end                
                submodel_label = cluster(gmmodel, target(:,4:6));
                target_label = zeros(N,M);
                for i=1:N
                    target_label(i,submodel_label(i)) = 1;
                end
                pred_label = target_label;
                mdl_WLR.Train(input_WLR,pred_label,1e-6);
            else
                mdl_WLR.W = mdl_WLR_prev_weight;
                mdl_WGP = mdl_WGP_prev;
                %EM_done = true;
            end
        end
        submodel_pred = mdl_WLR.computeLikelihood(input_WLR,[]); %N x M 
        %figure;
        %plot3(input(:,1),input(:,2),submodel_pred(:,1),'.b');
        %scatter(input(:,1),input(:,2),[],[submodel_pred(:,1) 0*submodel_pred(:,1) submodel_pred(:,2)]);
        %title('Prediction of submodel probability using LR');
        
        % weighted Gaussian Process
        
        if EM_done == false;
            for i=1:S % per running output dim = 6 axis
                for j=1:M % per submodel
                    disp(['Running submodel ' num2str(j) ' for ' num2str(i) ' axis']);
                     %hyperparameter initialize
            %         mdl_WGP{i}.hyp = log([repmat(1/(1^2),D,1);            
            %                                 0.1^2; % sigma_f^2
            %                                 0.1^2; % sigma_n^2
            %                                 0.1^2 % lambda
            %                                 ]);
            %        mdl_WGP{i}.hyp = [-7.61414885800113;-8.76972924992119;-8.09572933635572;-16.5780366637640;-21.1402296243598;-3.41497387195901;10.4438039765914;-0.384362931973412;-10.7711798458792];
                    %mdl_WGP{i}.hyp = [-7.61;-8.77;10.44;-10.771];
                    %mdl_WGP{i}.hyp = [1;1;1;1];
  
                    % compute subset of training if weight is small
                    sub_input = [];
                    sub_target = [];
                    sub_weight = [];
                    for k=1:N
                        %if pred_label(k,j) > 1e-9 % weight is non-zero, failing inversion of W
                        if submodel_pred(k,j) > 1e-9 % weight is non-zero, failing inversion of W
                            sub_input = [sub_input;input_WGP(k,:)];
                            sub_target = [sub_target;target(k,i)];
                            sub_weight = [sub_weight;submodel_pred(k,j)];
                        end
                    end
                    
                    if size(sub_input,1) ~= 0 % non samples belongs to this submodel
                        % number of hyperparameter should be same as number of features
                        % plus two hyperparameter for kernel
                        if EM_itr ==1
                            mdl_WGP{i,j}.hyp = repmat(0.1,D_GP+2,1);
                            mdl_WGP{i,j}.hyp(D_GP+2) = -10;
                        else
                            mdl_WGP{i,j}.hyp = repmat(0.1,D_GP+2,1);
                            mdl_WGP{i,j}.hyp(D_GP+2) = -10;
                            %mdl_WGP{i,j}.hyp;
                        end
                        %mdl_WGP{i}.hyp = [3.7;3.6;-1.49;-4];

                        % Run hyperparameter optimization
                        [mdl_WGP{i,j}.hyp, mdl_WGP{i,j}.nll, itr] = minimize(mdl_WGP{i,j}.hyp, 'wgp_lik', 10, sub_weight, sub_input, sub_target);
                        mdl_WGP{i,j}.weight = sub_weight;
                        mdl_WGP{i,j}.input = sub_input;
                        mdl_WGP{i,j}.target = sub_target;
                    end

                    % Gaussian process
                    [mdl_WGP{i,j}.pred_mu mdl_WGP{i,j}.pred_ss] = wgp_pred(mdl_WGP{i,j}.hyp, mdl_WGP{i,j}.weight, mdl_WGP{i,j}.input, mdl_WGP{i,j}.target, input_WGP);
                    mdl_WGP{i,j}.rmse = rms(mdl_WGP{i,j}.pred_mu-target(:,i));
        %             [mdl_WGP{i}.pred_mu_sub mdl_WGP{i}.pred_ss_sub] = wgp_pred(mdl_WGP{i}.hyp, sub_weight, sub_input, sub_target, sub_input);
                    if ~isempty(mdl_WGP{i,j}.pred_ss(mdl_WGP{i,j}.pred_ss <=0))
                        disp('Optimization error at building GP');
                        if EM_itr ==1
                            length_scale_init_val = 10;
                            %sigma0 = std(rslt_train);
                            sigma0 = [0.01 0.01 0.05];
                            sigmaF0 = sigma0;
                            dim = size(sub_input,2);
                            sigmaM0 = length_scale_init_val*ones(dim,1);
                            gprMdl{i,j} = fitrgp(sub_input, sub_target, 'Basis', 'constant', 'FitMethod','exact',...
                                    'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
                                    'KernelParameters',[sigmaM0;sigmaF0(1)],'Sigma',sigma0(1),'Standardize',1);
                            mdl_WGP{i,j}.hyp = [gprMdl{i,j}.ModelParameters.KernelParameters; gprMdl{i,j}.ModelParameters.KernelParameters(end,1)];
                            [mdl_WGP{i,j}.pred_mu mdl_WGP{i,j}.pred_ss] = wgp_pred(mdl_WGP{i,j}.hyp, mdl_WGP{i,j}.weight, mdl_WGP{i,j}.input, mdl_WGP{i,j}.target, input_WGP);
                        
                                %[pred pred_ss] = predict(gprMdl{i,j},sub_input);
                                %[pred(:,i) pred_sd(:,i)] = predict(gprMdl{i},feat_test);        
                                %pred_train(:,i) = resubPredict(gprMdl{i});
                        else 
                            %EM_done = true;
                            mdl_WLR = mdl_WLR_prev;
                            mdl_WGP = mdl_WGP_prev;                        
                        end
                        % predicted variance is negative
                        %[mdl_WGP{i,j}.hyp, mdl_WGP{i,j}.nll, itr] = minimize(mdl_WGP{i,j}.hyp, 'wgp_lik', 10, mdl_WGP{i,j}.weight, mdl_WGP{i,j}.input, mdl_WGP{i,j}.target);
                        %[mdl_WGP{i,j}.pred_mu mdl_WGP{i,j}.pred_ss] = wgp_pred(mdl_WGP{i,j}.hyp, mdl_WGP{i,j}.weight, mdl_WGP{i,j}.input, mdl_WGP{i,j}.target, input_WGP);
                    end
                    %[mu_test ss_test] = wgp_pred(mdl_WGP{i}.hyp, weight, input, target, input(1,:));
                    %mdl_WGP{i}.pd = makedist('Normal',mdl_WGP{i}.pred_mu,sqrt(mdl_WGP{i}.pred_ss));
                    %[mu ss] = wgp_pred(mdl_WGP{i}.hyp, weight, input, target, input);
                    mdl_WGP{i,j}.lik_target = pdf('Normal',target(:,i),mdl_WGP{i,j}.pred_mu,sqrt(mdl_WGP{i,j}.pred_ss));
                    for ck=1:size(mdl_WGP{i,j}.pred_mu,1)
                        if isnan(mdl_WGP{i,j}.lik_target(ck,1)) ==1
                            stop= 1;
                        end
                    end
                end            
                %lik_target(:,i) = mdl_WGP{i,j}.lik_target;
            end
        end
        %verify the quality of current estimator...
%         
        % label designation
%         figure;
%         scatter(input(:,1)',input(:,2)','.',[pred_label(:,1) pred_label(:,2) 0*pred_label(:,1)]')

        %figure;
        %scatter(input(:,1),input(:,2),[],[pred_label(:,1) 0*pred_label(:,1) pred_label(:,2)]);
        %title('Submodel prediction given observations p(m|s,s^,a)');
%         hold on;
%         for i=1:N
%             max_pred = max(pred_label');
%             for j=1:M
%                 if pred_label(i,M) == max_pred(i)
%                     plot(input(i,1),input(i,2),'.r');
%                 else
%                     plot(input(i,1),input(i,2),'.b');
%                 end
%             end
%         end
%         
%         figure
%         plot3(input(:,1),input(:,2),mdl_WGP{1}.pred_mu,'.r')
%         hold on
%         plot3(input(:,1),input(:,2),mdl_WGP{2}.pred_mu,'.b')
%         plot3(input(:,1),input(:,2),target,'.g')
        for i=1:S % per running output dim = 6 axis
                for j=1:M % per submodel
                    mdl_rmse(i,j) = mdl_WGP{i,j}.rmse;
                end
        end
        EM_err(EM_itr,:) = mean(mdl_rmse);
        
        if EM_itr == 5
            EM_done = true;
            sample_likhd = ones(N,M);
            for i=1:S % per axis
                for j=1:M % per submodel
                    for k=1:N % per sample
                        sample_likhd(k,j) = sample_likhd(k,j) * mdl_WGP{i,j}.lik_target(k,1);
                    end
                end                     
            end
            pred_label = sample_likhd.*submodel_pred./repmat(sum(sample_likhd.*submodel_pred,2),1,M);
            for i=1:size(pred_label,1)
                for j=1:M
                    if pred_label(i,j) == max(pred_label(i,:))
                        prdt_label(i,1) = j;
                    end
                end
            end
            for i=1:size(pred_label,1)
                for j=1:S
                    train_err(i,j) = mdl_WGP{j,prdt_label(i)}.pred_mu(i) - target(i,j);
                end                
            end
            for i=1:size(submodel_pred,1)
                for j=1:M
                    if submodel_pred(i,j) == max(submodel_pred(i,:))
                        prdt_sub_label(i,1) = j;
                    end
                end
            end
            for i=1:size(submodel_pred,1)
                for j=1:S
                    train_sub_err(i,j) = mdl_WGP{j,prdt_sub_label(i)}.pred_mu(i) - target(i,j);
                end                
            end
            rmse_err = [rms(train_err); rms(train_sub_err)];
        end
        pause(0.1);
    end
end