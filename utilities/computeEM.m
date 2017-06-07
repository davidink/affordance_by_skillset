function [mdl_WLR mdl_WGP] = computeEM(num_submodel, input, target, init_label)
    
    N = size(input,1); % number of training examples
    D = size(input,2); % size of input(feature) vector
    S = size(target,2);
    M = num_submodel;
    
    EM_done = false; EM_itr = 0;
    
    
    if nargin == 4 % given initialization
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
            gmmodel = fitgmdist(target,M,'Options',options);
            submodel_label = cluster(gmmodel, target);
            target_label = zeros(N,M);
            for i=1:N
                target_label(i,submodel_label(i)) = 1;
            end
            
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
            
            pred_label = target_label;
        else if run_init == false && EM_itr ==1
                target_label = zeros(N,M);
                for i=1:N
                    target_label(i,init_label(i)) = 1;
                end                
                pred_label = target_label;
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
        
%         figure;
%         plot3(input(:,1),input(:,2),pred_label(:,1),'r.');
%         hold on;
%         plot3(input(:,1),input(:,2),pred_label(:,2),'g.');        
%         pause(0.1);
        
        % M step
       
        % weighted logistic regression        
        if EM_itr ==1
            mdl_WLR = LogisticRegression(D,M,'linearFeatures');
        else
            %buffer
            mdl_WLR_prev = mdl_WLR;
            mdl_WGP_prev = mdl_WGP;
        end
        WLR_success = mdl_WLR.Train(input,pred_label,1e-6);
        if WLR_success == false;
            mdl_WLR = mdl_WLR_prev;
            mdl_WGP = mdl_WGP_prev;
            EM_done = true;
        end
        submodel_pred = mdl_WLR.computeLikelihood(input,[]); %N x M 
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
                    % number of hyperparameter should be same as number of features
                    % plus two hyperparameter for kernel
                    mdl_WGP{i,j}.hyp = repmat(1.0,D+2,1);
                    %mdl_WGP{i}.hyp = [3.7;3.6;-1.49;-4];

                    % compute subset of training if weight is small
                    sub_input = [];
                    sub_target = [];
                    sub_weight = [];
                    for k=1:N
                        if pred_label(k,j) > 1e-9 % weight is non-zero, failing inversion of W
                            sub_input = [sub_input;input(k,:)];
                            sub_target = [sub_target;target(k,i)];
                            sub_weight = [sub_weight;pred_label(k,j)];
                        end
                    end

        %             weight = pred_label(:,i);
        %             for j=1:N
        %                 if weight(j) ==0
        %                     weight(j) = 1e-6;
        %                 end
        %             end

                    % hyperparameter optimize
                    %[hyp_new, fX, i] = minimize(hyp, 'gp01lik', 100, x, y);
                    [mdl_WGP{i,j}.hyp, mdl_WGP{i,j}.nll, itr] = minimize(mdl_WGP{i,j}.hyp, 'wgp_lik', 10, sub_weight, sub_input, sub_target);
                    mdl_WGP{i,j}.weight = sub_weight;
                    mdl_WGP{i,j}.input = sub_input;
                    mdl_WGP{i,j}.target = sub_target;

                    % Gaussian process
                    [mdl_WGP{i,j}.pred_mu mdl_WGP{i,j}.pred_ss] = wgp_pred(mdl_WGP{i,j}.hyp, mdl_WGP{i,j}.weight, mdl_WGP{i,j}.input, mdl_WGP{i,j}.target, input);
                    mdl_WGP{i,j}.rmse = rms(mdl_WGP{i,j}.pred_mu-target(:,i));
        %             [mdl_WGP{i}.pred_mu_sub mdl_WGP{i}.pred_ss_sub] = wgp_pred(mdl_WGP{i}.hyp, sub_weight, sub_input, sub_target, sub_input);
                    while ~isempty(mdl_WGP{i,j}.pred_ss(mdl_WGP{i,j}.pred_ss <0))
                        % predicted variance is negative
                        [mdl_WGP{i,j}.hyp, mdl_WGP{i,j}.nll, itr] = minimize(mdl_WGP{i,j}.hyp, 'wgp_lik', 10, mdl_WGP{i,j}.weight, mdl_WGP{i,j}.input, mdl_WGP{i,j}.target);
                        [mdl_WGP{i,j}.pred_mu mdl_WGP{i,j}.pred_ss] = wgp_pred(mdl_WGP{i,j}.hyp, mdl_WGP{i,j}.weight, mdl_WGP{i,j}.input, mdl_WGP{i,j}.target, input);
                    end
                    %[mu_test ss_test] = wgp_pred(mdl_WGP{i}.hyp, weight, input, target, input(1,:));
                    %mdl_WGP{i}.pd = makedist('Normal',mdl_WGP{i}.pred_mu,sqrt(mdl_WGP{i}.pred_ss));
                    %[mu ss] = wgp_pred(mdl_WGP{i}.hyp, weight, input, target, input);
                    mdl_WGP{i,j}.lik_target = pdf('Normal',target(:,i),mdl_WGP{i,j}.pred_mu,sqrt(mdl_WGP{i,j}.pred_ss));
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
        
        if EM_itr == 10
            EM_done = true;
        end
        pause(0.1);
    end
end