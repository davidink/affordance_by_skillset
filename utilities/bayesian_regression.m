function [w_mu w_sigsq] = bayesian_regression(X,Y,alpha_inv,beta_inv)

iNumRegressors   = size(X,2);

% p(w|t) = N(w,w_mu,w_sq), p(w|a) = N(w|0,a_inv*I)
% alpha = initial variance of weight w
alpha = inv(alpha_inv * eye(iNumRegressors)); 

% p(y|x,w,b) = N(t|y(x,w),b_inv)
% beta = precision(inverse variance) of noise of target output t
beta = inv(beta_inv * eye(iNumRegressors)); 

% bayesian linear regression
w_sigsq = inv(alpha + beta * X' * X);
w_mu = beta * w_sigsq * X' * Y;    

for j=1:size(w_sigsq,1)
    for k=j:size(w_sigsq,1)
        w_sigsq(k,j) = w_sigsq(j,k);
    end
end
        
% iNumMeasurements = size(X,1);
% iNumRegressors   = size(X,2);
% 
% big_sigma = small_sigma_squared * eye(iNumMeasurements); %a = 0.01 * eye(5);
% big_omega = eta_squared * eye(iNumRegressors); %b = 0.01 * eye(5);
% 
% %disp('covariance matrix and mean vector of posterior distribution')
% lambda = inv(X' * inv(big_sigma) * X + inv(big_omega));
% mu     = lambda * X' * inv(big_sigma) * Y;

end
