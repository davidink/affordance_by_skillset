function [w_mu w_sigsq] = bayesian_regression_update(X,Y,alpha_inv,beta_inv)

num_sample = size(Y,1);
iNumRegressors   = size(X,2);

for i=1:num_sample
    if i==1 % initial modeling
        X_in = X(i,:);
        Y_in = Y(i,:);
        
        % p(w|t) = N(w,w_mu,w_sq), p(w|a) = N(w|0,a_inv*I)
        % alpha = initial variance of weight w
        alpha = inv(alpha_inv * eye(iNumRegressors)); 

        % p(y|x,w,b) = N(t|y(x,w),b_inv)
        % beta = precision(inverse variance) of noise of target output t
        beta = inv(beta_inv * eye(iNumRegressors)); 

        % bayesian linear regression
        w_sigsq{i} = inv(alpha + beta * X_in' * X_in);
        w_mu{i}     = beta * w_sigsq{i} * X_in' * Y_in;    
    else
        X_in = X(1:i,:);
        Y_in = Y(1:i,:);
        w_sigsq_prior = w_sigsq; % S_0
        w_mu_prior = w_mu;
        w_sigsq{i} = inv( inv(w_sigsq{i-1}) + beta * X_in' * X_in);
        w_mu{i} = w_sigsq{i}*(inv(w_sigsq{i-1})*w_mu{i-1} + beta * X_in' * Y_in);
    end
    
    for j=1:size(w_sigsq{i},1)
        for k=j:size(w_sigsq{i},1)
            w_sigsq{i}(k,j) = w_sigsq{i}(j,k);
        end
    end
end

end
