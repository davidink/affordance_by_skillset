clear all
close all

load gprdata1

gprMdl1 = fitrgp(x,y,'KernelFunction','squaredexponential');

%hyperparameter initialize
%hyp = [1/(18.3539^2);3.1559^2;0.7162^2];
% hyp = [1/(2.99569^2);3.9593^2;0.1907^2];
hyp = [1/(2.99569^2);3.9593^2;0.1^2];
hyp = log(hyp);
weight = rand(size(x,1),1);
%weight(x>10 & x<40) = 1e-6;

% hyperparameter optimize
%[hyp_new, fX, i] = minimize(hyp, 'gp01lik', 100, x, y);
[hyp_new, fX, i] = minimize(hyp, 'wgp_lik', 50, weight, x, y);

% Gaussian process
[mu S2] = wgp_pred(hyp_new, weight, x, y, x);

figure();
hold on;
plot(x,y,'r.');
plot(x,mu,'b');



% iter = 0;
% [fX, dfX] = gp01lik(hyp, x, y)
% 
% sigma0 = 0.2;
% kparams0 = [3.5, 6.2];
% gprMdl2 = fitrgp(x,y,'KernelFunction','squaredexponential',...
%      'KernelParameters',kparams0,'Sigma',sigma0);
%  
%  
% hyp = [1/(2.99569^2);3.9593^2;0.1907^2];
% hyp = log(hyp);
% 
% [fX, dfX] = gp01lik(hyp, x, y)
 
% ypred1 = resubPredict(gprMdl1);
% ypred2 = resubPredict(gprMdl2);
% 
% figure();
% plot(x,y,'r.');
% hold on
% plot(x,ypred1,'b');
% plot(x,ypred2,'g');
% xlabel('x');
% ylabel('y');
% legend({'data','default kernel parameters',...
% 'kparams0 = [3.5,6.2], sigma0 = 0.2'},...
% 'Location','Best');
% title('Impact of initial kernel parameter values');
% hold off
 
 
