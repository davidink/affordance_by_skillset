clear all
close all

%load data

%load gprdata;
load dataset1;

% sigma0 = std(ytrain);
% sigmaF0 = sigma0;
% d = size(Xtrain,2);
% sigmaM0 = 10*ones(d,1);
% 
% gprMdl = fitrgp(Xtrain,ytrain,'Basis','constant','FitMethod','exact',...
% 'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
% 'KernelParameters',[sigmaM0;sigmaF0],'Sigma',sigma0,'Standardize',1);
% 
% L = loss(gprMdl,Xtest,ytest)
% 
% ypred = predict(gprMdl,Xtest);
% 
% figure;
% plot(ytest,'r');
% hold on;
% plot(ypred,'b');
% legend('True response','GPR predicted values','Location','Best');
% hold off

num_submodel = 2;
%run EM

%[mdl_wlr mdl_wgp] = computeEM(num_submodel, Xtrain, ytrain);

[mdl_WLR mdl_WGP] = computeEM(num_submodel, X, Y);

% load EM_models_1.mat

submodel_pred = mdl_WLR.computeLikelihood(X_all,[]); %N x M p(m=j|s,a)
for i=1:num_submodel
    [submodel{i}.pred_mu submodel{i}.pred_ss] = wgp_pred(mdl_WGP{i}.hyp, mdl_WGP{i}.weight, mdl_WGP{i}.input, mdl_WGP{i}.target, X_all);
end

%verify the quality of current estimator...
%         
% label designation
% figure;
% title('submodel prediction on each datapoint');
% hold on;
% grp_1 =[];
% grp_2 =[];
% for i=1:size(X_all,1)
%     max_pred = max(submodel_pred');
%     for j=1:num_submodel
%         if submodel_pred(i,1) == max_pred(i)
%             grp_1 = [grp_1;X_all(i,1) X_all(i,2)];
%             %plot(X_all(i,1),X_all(i,2),'.r');
%         else
%             grp_2 = [grp_2;X_all(i,1) X_all(i,2)];
%             %plot(X_all(i,1),X_all(i,2),'.b');
%         end
%     end
% end
% plot(grp_1(:,1),grp_1(:,2),'.r');
% plot(grp_2(:,1),grp_2(:,2),'.b');

figure
surf(reshape(X_all(:,1),101,101),reshape(X_all(:,2),101,101),reshape(submodel{1}.pred_mu,101,101));
figure
surf(reshape(X_all(:,1),101,101),reshape(X_all(:,2),101,101),reshape(submodel{2}.pred_mu,101,101));
figure
surf(reshape(X_all(:,1),101,101),reshape(X_all(:,2),101,101),reshape(Y_all,101,101));
% figure
% plot3(X_all(:,1),X_all(:,2),submodel{1}.pred_mu,'.r');
% hold on
% plot3(X_all(:,1),X_all(:,2),submodel{2}.pred_mu,'.b');
% plot3(X_all(:,1),X_all(:,2),Y_all,'.g');
% 
% figure
% surf(reshape(X_all(:,1),101,101),reshape(X_all(:,2),101,101),reshape(submodel{1}.pred_mu,101,101));


% compute for new input