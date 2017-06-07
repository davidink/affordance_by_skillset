clear all
close all

%reset random number generator
% s = RandStream('mt19937ar','Seed',0)
% RandStream.setGlobalStream(s);

%dbstop if error
current_result_saving_folder = 'rslts2';
          
load 'features_n_results_interaction4';

colorOrder = [0.8 0 0; 0 0.8 0; 0 0 0.8; 0.8 0.8 0; 0.8 0 0.8; 0 0.8 0.8; 0.8 0.8 0.8];

results_all = [results results_interaction];

% figure;
% plot3(results(:,2),results(:,3),results(:,4),'.r');
% xlabel('y');
% ylabel('z');
% zlabel('theta');
% figure;
% plot3(results(:,6),results(:,7),results(:,8),'.r');
% xlabel('y');
% ylabel('z');
% zlabel('theta');

num_cv = 3;
[feat_grp rslt_grp] = generate_cv_dataset(features,results_all,num_cv);


%[feat_train_s rslt_train_s feat_test_s rslt_test_s] = generate_train_test_dataset(features,results_all,55);

% figure;
% %plot(results(:,2),results(:,3),'.r');
% plot3(results(:,2),results(:,3),results(:,4),'.r');
% hold on;
% result_nointeraction = results(results_interaction==0,:);
% plot3(result_nointeraction(:,2),result_nointeraction(:,3),result_nointeraction(:,4),'ob');
% xlabel('y');
% ylabel('z');
% zlabel('rot');
% %axis equal;
% axis([-0.1 0.1 0 0.3]);
% plot3(results(:,6),results(:,7),results(:,8),'.g');

% num_submodel = 1;

for num_submodel=4:4
    %for set_num=1:8
    for set_num=7:7
        for cv=1:num_cv
            disp(['Running set_num ' num2str(set_num) ' of cv' num2str(cv)]);
            % generate train/test dataset
            feat_test = feat_grp{cv};
            rslt_test = rslt_grp{cv}(:,[2:4 6:8]);

            feat_train=[];
            rslt_train=[];
            for i=1:num_cv
                if cv ~= i
                    feat_train = [feat_train;feat_grp{i}];
                    rslt_train = [rslt_train;rslt_grp{i}(:,[2:4 6:8])];
                end
            end

            feat_all = [feat_train;feat_test];
            rslt_all = [rslt_train;rslt_test];

            rslt_train_init = bsxfun(@rdivide,bsxfun(@minus,rslt_train,mean(rslt_train)),std(rslt_train));

        %     for i=1:size(rslt_train,1)
        %         
        %         for j=1:6
        %         rslt_train(i,3) = rslt_train(i,3)/10;
        %         rslt_train(i,6) = rslt_train(i,6)/10;        
        %         end
        %     end
        %    
            %run EM

            %initialize with GMM
            options = statset('Display','final');
            gmm_init = false;
            while ~gmm_init
                %gmmodel = fitgmdist(rslt_train(:,2:4),num_submodel,'Options',options);
                %submodel_label = cluster(gmmodel, rslt_train(:,2:4));
                try
                    %gmmodel = fitgmdist([rslt_train(:,2:4) rslt_train(:,6:8)],num_submodel,'Options',options);
                    %gmmodel = fitgmdist(rslt_train_init(:,4:5),num_submodel,'Options',options,'Replicates',10);
                    gmmodel = fitgmdist(rslt_train_init,num_submodel,'Options',options,'Replicates',10);
                    gmm_init = true;
                catch ME
                    gmm_init = false;
                end
            end
            %submodel_label = cluster(gmmodel, [rslt_train(:,2:4) rslt_train(:,6:8)]);
            %init_submodel_label = cluster(gmmodel, rslt_train_init(:,4:5));
            init_submodel_label = cluster(gmmodel, rslt_train_init);

            figure;
            title('Initial label from GMM on first object');
            hold on;
            for i=1:num_submodel
                idx = find(init_submodel_label==i);
                plot3(rslt_train(idx,1),rslt_train(idx,2),rslt_train(idx,3),'.','Color',colorOrder(i,:),'MarkerSize',14);
            end
            xlabel('y_{p1}');
            ylabel('z_{p1}');
            zlabel('\theta_{p1}');
            axis equal;
            figure;
            title('Initial label from GMM on second object');
            hold on;
            for i=1:num_submodel
                idx = find(init_submodel_label==i);
                plot3(rslt_train(idx,4),rslt_train(idx,5),rslt_train(idx,6),'.','Color',colorOrder(i,:),'MarkerSize',14);
                %plot3([rslt_train(idx,1);rslt_train(idx,4)],[rslt_train(idx,2);rslt_train(idx,5)],[rslt_train(idx,3);rslt_train(idx,6)],'o');
        %         for j=1:size(idx,1)
        %             plot3([rslt_train(idx(j),1) rslt_train(idx(j),4)],[rslt_train(idx(j),2) rslt_train(idx(j),5)],[rslt_train(idx(j),3) rslt_train(idx(j),6)],'k');
        % %             plot(rslt_train(idx(j),1)
        %         end
                %plot3(rslt_train(idx,1),rslt_train(idx,2),rslt_train(idx,3),'.');
            end
            xlabel('y_{p2}');
            ylabel('z_{p2}');
            zlabel('\theta_{p2}');
            axis equal;

            % Selection of different sets of features

            % 1 set
            if set_num == 1
                feat_WLR_idx = [1:8]; % object features 
                feat_WGP_idx = [1:8]; % object features
            end
            if set_num == 2
                feat_WLR_idx = [1:21]; % obj+contact 
                feat_WGP_idx = [1:21]; 
            end
            if set_num == 3
                feat_WLR_idx = [22:101]; % context only
                feat_WGP_idx = [22:101]; 
            end
            if set_num == 4
                feat_WLR_idx = [1:101]; % obj+contact+context 
                feat_WGP_idx = [1:101]; 
            end
            if set_num == 5
                feat_WLR_idx = [1:21]; % obj+contact 
                feat_WGP_idx = [22:101]; %ctxt only
            end
            if set_num == 6
                feat_WLR_idx = [1:21]; % obj+contact 
                feat_WGP_idx = [1:101]; %all 
            end
            if set_num == 7
                feat_WLR_idx = [22:101]; % obj+contact 
                feat_WGP_idx = [1:21]; %obj only
            end
            if set_num == 8
                feat_WLR_idx = [22:101]; % obj+contact 
                feat_WGP_idx = [1:101]; % all
            end
            % 2 set
        %     feat_WLR_idx = [9:21]; % contact only 
        %     feat_WGP_idx = [9:21]; 
            % 3 set
        %     feat_WLR_idx = [22:101]; % context only 
        %     feat_WGP_idx = [22:101];
            % 4 set
    %         feat_set_idx = 4;
    %         feat_WLR_idx = [1:101]; % obj+cnct+ctxt 
    %         feat_WGP_idx = [1:101];
            % 5 set
        %     feat_WLR_idx = [1:21]; 
        %     feat_WGP_idx = [1:21];
            % 6 set
        %     feat_WLR_idx = [1:21]; 
        %     feat_WGP_idx = [1:101];
            % 7 set
        %     feat_WLR_idx = [22:101]; 
        %     feat_WGP_idx = [1:21];
            % 8 set
        %     feat_WLR_idx = [22:101]; 
        %     feat_WGP_idx = [1:101];

            %feat_WLR_idx = [6:17]; % contact features
            %feat_WLR_idx = [18:40]; %context features(cont)
            %feat_WLR_idx = [41:357];  %context features(discrete)
            %feat_WLR_idx = [1:357]; % all features
            %feat_WLR_idx = [1:40];  % all features

            feat_train_WLR = feat_train(:,feat_WLR_idx);
            feat_test_WLR = feat_test(:,feat_WLR_idx);
            feat_all_WLR = feat_all(:,feat_WLR_idx);

            %feat_WGP_idx = [6:17]; % contact features
            %feat_WGP_idx = [18:40]; %context features(cont)
            %feat_WGP_idx = [41:357];  %context features(discrete)
            %feat_WGP_idx = [1:357]; % all features
            %feat_WGP_idx = [1:40]; % all features
            feat_train_WGP = feat_train(:,feat_WGP_idx);
            feat_test_WGP = feat_test(:,feat_WGP_idx);
            feat_all_WGP = feat_all(:,feat_WGP_idx);

            % build prediction models
            [mdl_WLR mdl_WGP] = computeEM_ind(num_submodel, feat_train_WLR, feat_train_WGP, rslt_train, init_submodel_label);
            pause(0.1);
            WLR_test_pred = mdl_WLR.computeLikelihood(feat_test_WLR,[]);
            WLR_all_pred = mdl_WLR.computeLikelihood(feat_all_WLR,[]);
            WLR_train_pred = mdl_WLR.computeLikelihood(feat_train_WLR,[]);

            for i=1:6
                for j=1:num_submodel
                    [WGP_test_pred{i,j}.pred_mu WGP_test_pred{i,j}.pred_ss] = wgp_pred(mdl_WGP{i,j}.hyp, mdl_WGP{i,j}.weight, mdl_WGP{i,j}.input, mdl_WGP{i,j}.target, feat_test_WGP);
                end
            end

            %find the optimal submodel labels
            pred_test_label = [];
            for i=1:size(feat_test,1)
                for j=1:num_submodel
                    if WLR_test_pred(i,j) == max(WLR_test_pred(i,:))
                        pred_test_label(i,1) = j;
                    end
                end
            end
            pred_train_label = [];
            for i=1:size(feat_train,1)
                for j=1:num_submodel
                    if WLR_train_pred(i,j) == max(WLR_train_pred(i,:))
                        pred_train_label(i,1) = j;
                    end
                end
            end
        %     for i=1:size(feat_test,1)
        %         for j=1:num_submodel
        %             if pred_test_likhd(i,j) == max(pred_test_likhd(i,:))
        %                 pred_test_label(i,1) = j;
        %             end
        %         end
        %     end
            pred_all_label = [];
            for i=1:size(feat_all,1)
                for j=1:num_submodel
                    if WLR_all_pred(i,j) == max(WLR_all_pred(i,:))
                        pred_all_label(i,1) = j;
                    end
                end
            end

            % draw predicted label for train inputs
%             figure;
%             hold on;
%             for i=1:num_submodel
%                 idx = find(pred_train_label==i);
%                 plot3(rslt_train(idx,1),rslt_train(idx,2),rslt_train(idx,3),'.');
%                 %text(rslt_all(idx,1),rslt_all(idx,2),rslt_all(idx,3),[' ' num2str(i)]);
%             end
%             xlabel('y');
%             ylabel('z');
%             zlabel('rot');
%             title('Predicted labels for train inputs for first object');
% 
%             figure;
%             hold on;
%             for i=1:num_submodel
%                 idx = find(pred_train_label==i);
%                 plot3(rslt_train(idx,4),rslt_train(idx,5),rslt_train(idx,6),'.');
%                 %text(rslt_all(idx,4),rslt_all(idx,5),rslt_all(idx,6),[' ' num2str(i)]);
%             end
%             xlabel('y');
%             ylabel('z');
%             zlabel('rot');
%             title('Predicted labels for train inputs for second object');

            %find the prediction value from selected labels
            for i=1:size(feat_test,1)
                for k=1:6
                        pred_test_mu(i,k) = WGP_test_pred{k,pred_test_label(i)}.pred_mu(i);
                        pred_test_ss(i,k) = WGP_test_pred{k,pred_test_label(i)}.pred_ss(i);
%                         pred_test_mu(i,k) = WGP_test_pred{k,1}.pred_mu(i);
%                         pred_test_ss(i,k) = WGP_test_pred{k,1}.pred_ss(i);
                end    
            end

            %visualize final predicted submodel labels
%             figure;
%             hold on;
%             for i=1:num_submodel
%                 idx = find(pred_test_label==i);
%                 plot3(pred_test_mu(idx,1),pred_test_mu(idx,2),pred_test_mu(idx,3),'o');
%                 %text(pred_test_mu(idx,1),pred_test_mu(idx,2),pred_test_mu(idx,3),[' ' num2str(i)]);
%                 plot3(rslt_test(idx,1),rslt_test(idx,2),rslt_test(idx,3),'bx');
%             end
%             for i=1:size(rslt_test,1)
%                 plot3([pred_test_mu(i,1);rslt_test(i,1)], [pred_test_mu(i,2);rslt_test(i,2)],[pred_test_mu(i,3);rslt_test(i,3)],'c');
%             end
%             xlabel('y');
%             ylabel('z');
%             zlabel('rot');
%             title('Predicted values for first object');
% 
%             figure;
%             hold on;
%             for i=1:num_submodel
%                 idx = find(pred_test_label==i);
%                 plot3(pred_test_mu(idx,4),pred_test_mu(idx,5),pred_test_mu(idx,6),'o');
%                 %text(pred_test_mu(idx,4),pred_test_mu(idx,5),pred_test_mu(idx,6),[' ' num2str(i)]);
%                 plot3(rslt_test(idx,4),rslt_test(idx,5),rslt_test(idx,6),'bx');
%             end
%             for i=1:size(rslt_test,1)
%                 plot3([pred_test_mu(i,4);rslt_test(i,4)], [pred_test_mu(i,5);rslt_test(i,5)],[pred_test_mu(i,6);rslt_test(i,6)],'c');
%             end
%             xlabel('y');
%             ylabel('z');
%             zlabel('rot');
%             title('Predicted values for second object');

            %compute error of prediction for each 6 axes
            %rslt_test_axis = [rslt_test(:,2:4) rslt_test(:,6:8)];
            pred_test_err(cv,:) = rms(rslt_test - pred_test_mu);

%             for i=1:num_submodel
%                 for j=1:size(feat_test,1)
%                     for k=1:6
%                         pred_test_by_label{i}(j,k) = WGP_test_pred{k,i}.pred_mu(j);
%                     end
%                 end
%             end
%             for i=1:num_submodel
%                 pred_test_err_by_label{i} = abs(rslt_test - pred_test_by_label{1,i});
%             end
%             for i=1:6
%                 for j=1:size(feat_test,1)
%                     minimum(j,i) = 100000;
%                     for k=1:num_submodel
%                         if pred_test_err_by_label{k}(j,i) < minimum(j,i)
%                             minimum(j,i) = pred_test_err_by_label{k}(j,i);
%                             pred_test_label_by_err(j,i) = k;
%                         end
%                     end
%                 end        
%             end
%             pred_label_all = [pred_test_label_by_err pred_test_label];
            for i=1:size(rslt_test,1)
                for k=1:6
                    pred_test_ll(i,k) = log(pdf('Normal',rslt_test(i,k),pred_test_mu(i,k),sqrt(pred_test_ss(i,k))));
                    %1/sqrt(2*pi*pred_test_ss(1,1))*exp((rslt_test_axis(1,1)-pred_test_mu(1,1))^2*-1/(2*pred_test_ss(1,1)))
                end
            end
            pred_test_ll_mean(cv,:) = mean(pred_test_ll);
            
            for i=1:num_submodel
                idx = find(pred_test_label==i);
                if isempty(idx)
                    pred_test_err_by_label(i,:) = [-1 -1 -1 -1 -1 -1];
                else
                    pred_test_err_by_label(i,:) = rms(rslt_test(idx,:) - pred_test_mu(idx,:));                
                end
            end

        end
        current_save_file_name = [current_result_saving_folder '/sub_' num2str(num_submodel) '_featset' num2str(set_num) '.mat'];
        save(current_save_file_name,'pred_test_err','pred_test_ll_mean');
    end
    %avg_pred_error(num_submodel,:) = mean(pred_test_err);
    %std_pred_error(num_submodel,:) = std(pred_test_err);
    %data_pred_error{num_submodel} = pred_test_err;

end
