function w_sigsq = bayesian_regression_varonly(X,alpha_inv,beta_inv)

iNumRegressors   = size(X,2);

% p(w|t) = N(w,w_mu,w_sq), p(w|a) = N(w|0,a_inv*I)
% alpha = initial variance of weight w
alpha = inv(alpha_inv * eye(iNumRegressors)); 

% p(y|x,w,b) = N(t|y(x,w),b_inv)
% beta = precision(inverse variance) of noise of target output t
beta = inv(beta_inv * eye(iNumRegressors)); 

w_sigsq = inv(alpha + beta * X' * X);

end
