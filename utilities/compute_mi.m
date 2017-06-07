close all
clear all

dataFolder='~/catkin_hydro_ws/data/';

features_obsvd = csvread([dataFolder,'features_init.txt']);
%results = csvread([dataFolder,'results_init.txt']);

%features = features(:,1:38);
%features_obsvd = features_obsvd(:,1:497);

num_cand = 10;
num_horizon = 3;
for i=1:num_cand
    m = csvread([dataFolder,'features_candidate_',num2str(i),'.txt']);
    if i==1
        feat_cand = m;
    else
        feat_cand = [feat_cand; m];
    end        
end

small_sigma_squared = 0.0005; %esitmated variance of y distribution
eta_squared = 0.0001; %estimated variance of w (weight) distribution

%w_var_obsvd = bayesian
big_sigma = small_sigma_squared * eye(size(features_obsvd,1)); %a
big_omega = eta_squared * eye(size(features_obsvd,2));
w_var_obsvd = inv(features_obsvd' * inv(big_sigma) * features_obsvd + inv(big_omega));
remain_feat_cand = [];
for i=1:num_cand
    sample_feat_cand = feat_cand(num_horizon*(i-1)+1:num_horizon*(i-1)+num_horizon,:);
    if(i~=1) remain_feat_cand = feat_cand(1:num_horizon*(i-2)+num_horizon,:);
    end
    if(i~=num_cand) remain_feat_cand = [remain_feat_cand; feat_cand(num_horizon*i+1:end,:)];
    end
    for j=1:num_horizon
        current_sample_feat = sample_feat_cand(j,:);
        sample_unobsvd = sample_feat_cand;
        sample_unobsvd(j,:) = [];
        features_unobsvd = [remain_feat_cand;sample_unobsvd(:,:)];
        w_var_obsvd = bayesian_regression_varonly(features_obsvd,small_sigma_squared,eta_squared);
        w_var_unobsvd = bayesian_regression_varonly(features_unobsvd,small_sigma_squared,eta_squared);
        y_var_obsvd = small_sigma_squared + current_sample_feat * w_var_obsvd * current_sample_feat';
        y_var_unobsvd = small_sigma_squared + current_sample_feat * w_var_unobsvd * current_sample_feat';
        MI(i,j) = 0.5 * log(y_var_obsvd/y_var_unobsvd);
    end
end

MI_score = MI(:,1)+MI(:,2)+MI(:,3)
[max_MI idx] = max(MI_score);

fileID = fopen([dataFolder,'push_selected.txt'],'w');
fprintf(fileID,'%d',idx);
fclose(fileID);


% %%No grouping
% Y1 = [obj_position_after(:,2)-obj_position_before(:,2) obj_position_after(:,3)-obj_position_before(:,3) obj_theta_after(:,1)-obj_theta_before(:,1)];
% 
% small_sigma_squared = 0.0005; %esitmated variance of y distribution
% eta_squared = 0.0001; %estimated variance of w (weight) distribution
% %[mu lambda] = bayesian_regression([ones(size(features,1),1) features],Y1(:,1),small_sigma_squared,eta_squared);
% %[mu lambda] = bayesian_regression(X1,Y1(:,1),small_sigma_squared,eta_squared);
% %[mu lambda] = bayesian_regression([ones(size(features,1),1) features],Y1(:,1),small_sigma_squared,eta_squared);
% [mu{1} lambda{1}] = bayesian_regression(features,Y1(:,1),small_sigma_squared,eta_squared);
% 
% small_sigma_squared = 0.0005; %esitmated variance of y distribution
% eta_squared = 0.0001; %estimated variance of w (weight) distribution
% [mu{2} lambda{2}] = bayesian_regression(features,Y1(:,2),small_sigma_squared,eta_squared);
% 
% small_sigma_squared = 0.0005; %esitmated variance of y distribution
% eta_squared = 0.0001; %estimated variance of w (weight) distribution
% [mu{3} lambda{3}] = bayesian_regression(features,Y1(:,3),small_sigma_squared,eta_squared);
% 
% for k=1:3
%     for i=1:size(lambda{k},1)
%         for j=i:size(lambda{k},1)
%         lambda{k}(j,i) = lambda{k}(i,j);        
%         end
%     end
% end
% 
% for i=1:3
%     w{i} = mvnrnd(mu{i},lambda{i});
%     E_Y{i} = w{i} * features';    
% end
% 
% E_Yt = [E_Y{1}' E_Y{2}' E_Y{3}'];
% plot(E_Yt(:,1),E_Yt(:,2),'xb');
% 
% for i=1:size(E_Yt,1)
%     plot([Y1(i,1);E_Yt(i,1)],[Y1(i,2);E_Yt(i,2)],'-y');
% end
% 
% RMSE = rms(Y1 - E_Yt)
% 
% % save as a file for weights
% weight_sample_num = 100;
% for i=1:3
%     for j=1:weight_sample_num
%         weight = mvnrnd(mu{i},lambda{i});
%         fileID = fopen([dataFolder,'weight_samples_',num2str(i),'.txt'],'a+');
%         for k=1:size(weight,2)
%             fprintf(fileID,'%f,',weight(k));
%         end
%         fprintf(fileID,'\n');
%         fclose(fileID);
%     end
% end
% 
% %%
% 
% %% grouping
% 
% %group_a
% small_sigma_squared = 0.0005; %esitmated variance of y distribution
% eta_squared = 0.0001; %estimated variance of w (weight) distribution
% [mu_a{1} lambda_a{1}] = bayesian_regression(features_a,group_a(:,1),small_sigma_squared,eta_squared);
% 
% small_sigma_squared = 0.0005; %esitmated variance of y distribution
% eta_squared = 0.0001; %estimated variance of w (weight) distribution
% [mu_a{2} lambda_a{2}] = bayesian_regression(features_a,group_a(:,2),small_sigma_squared,eta_squared);
% 
% small_sigma_squared = 0.0005; %esitmated variance of y distribution
% eta_squared = 0.0001; %estimated variance of w (weight) distribution
% [mu_a{3} lambda_a{3}] = bayesian_regression(features_a,group_a(:,3),small_sigma_squared,eta_squared);
% 
% for k=1:3
%     for i=1:size(lambda_a{k},1)
%         for j=i:size(lambda_a{k},1)
%         lambda_a{k}(j,i) = lambda_a{k}(i,j);        
%         end
%     end
% end
% 
% for i=1:3
%     w_a{i} = mvnrnd(mu_a{i},lambda_a{i});
%     E_Y_a{i} = w_a{i} * features_a';    
% end
% 
% figure;
% gscatter(X(:,1),X(:,2),cluster_yz);
% hold on
% 
% E_Yt_a = [E_Y_a{1}' E_Y_a{2}' E_Y_a{3}'];
% plot(E_Yt_a(:,1),E_Yt_a(:,2),'xr');
% 
% for i=1:size(E_Yt_a,1)
%     plot([group_a(i,1);E_Yt_a(i,1)],[group_a(i,2);E_Yt_a(i,2)],'-y');
% end
% 
% RMSE_a = rms(group_a - E_Yt_a)
% 
% %group_b
% small_sigma_squared = 0.0015; %esitmated variance of y distribution
% eta_squared = 0.0005; %estimated variance of w (weight) distribution
% [mu_b{1} lambda_b{1}] = bayesian_regression(features_b,group_b(:,1),small_sigma_squared,eta_squared);
% 
% small_sigma_squared = 0.0005; %esitmated variance of y distribution
% eta_squared = 0.0005; %estimated variance of w (weight) distribution
% [mu_b{2} lambda_b{2}] = bayesian_regression(features_b,group_b(:,2),small_sigma_squared,eta_squared);
% 
% small_sigma_squared = 0.0005; %esitmated variance of y distribution
% eta_squared = 0.0005; %estimated variance of w (weight) distribution
% [mu_b{3} lambda_b{3}] = bayesian_regression(features_b,group_b(:,3),small_sigma_squared,eta_squared);
% 
% for k=1:3
%     for i=1:size(lambda_b{k},1)
%         for j=i:size(lambda_b{k},1)
%         lambda_b{k}(j,i) = lambda_b{k}(i,j);        
%         end
%     end
% end
% 
% for i=1:3
%     w_b{i} = mvnrnd(mu_b{i},lambda_b{i});
%     E_Y_b{i} = w_b{i} * features_b';    
% end
% 
% E_Yt_b = [E_Y_b{1}' E_Y_b{2}' E_Y_b{3}'];
% plot(E_Yt_b(:,1),E_Yt_b(:,2),'xc');
% 
% for i=1:size(E_Yt_b,1)
%     plot([group_b(i,1);E_Yt_b(i,1)],[group_b(i,2);E_Yt_b(i,2)],'-y');
% end
% 
% RMSE_b = rms(group_b - E_Yt_b)


