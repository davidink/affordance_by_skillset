close all
clear all

[feat_all rslt_all] = load_test_dataset();

%% Grouping
real_push(:) = 0.2 - feat_all(:,12);
%rslt_all(:,2) = rslt_all(:,2) - feat_all(:,14);
X = [rslt_all(:,1) rslt_all(:,2)];
cluster_num = 3;
%options = statset('Display','final');
%gmmodel = fitgmdist(X,cluster_num,'Options',options);
load 'GMM_8.mat'
%load 'GMM_4.mat'
cluster_yz = cluster(gmmodel, X);

figure
%gscatter(X(:,1),X(:,2),cluster_yz,'rbg');

gscatter(X(:,1),X(:,2),cluster_yz,'rbgc');
legend('Location','northeast');
xlabel('y_p (m) ');
ylabel('z_p (m)');
hold on
axis equal;
set(gca,'Xdir','reverse')

% figure;
% plot(real_push,rslt_all(:,2)','.r');
% hold on;
% x= 0:0.01:0.15;
% plot(x,x,'b');
% axis equal