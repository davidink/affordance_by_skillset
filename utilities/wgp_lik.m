function [fX, dfX] = wgp_lik(hyp, weight, input, target);

% Compute minus log likelihood and its derivatives with respect to
% hyperparameters for a Gaussian Process for regression.
%
% Implementation of Gaussian Process for Regression, using a
% simlpe Gaussian covariance function and a noise term. The covariance
% function is controlled by hyperparameters.
%
% usage: [fX dfX] = wgp_lik(hyp, weight, input, target)
%
% where:
%
% hyp    is a (column) vector (of size D+3) of hyperparameters
% weight is a column vector of size N of inputs
% input  is a N by D matrix of training inputs
% target is a (column) vector (of size N) of targets
% fX     is the returned value of minus log likelihood
% dfX    is a (column) vector (of size D+3) of partial derivatives
%        of minus the log likelihood wrt each of the hyperparameters
%
% The form of the covariance function is
%
% C(x^p,x^q) = v1 * exp( -0.5 * sum_{d=1..D} w_d * (x^p_d - x^q_d)^2 )
%            + v0 * W^-1
%            
%
% where the first term is Gaussian, the second term with the kronecker
% delta is the noise contribution, and thrid term is weights. In this function, the hyperparameters w_d,
% v2, v1, and v0 are collected in the vector X as follows:
%
% hyp = [ log(w_1) % legnth scale of input vector 1
%       log(w_2) 
%        .
%       log(w_D)
%       log(v1)  % kernel amplitutde sigma_f^2 
%       log(v0)  % noise signal variance lambda
%     ]
%
% Note: the reason why the log of the parameters are used in hyp is that we may
% then used unconstrained optimisation, while the hyperparameters themselves
% must, naturally, always be positive.
%
% This function can conveniently be used with the "minimize" function to train
% a Gaussian process:
%
% [hyp_new, fX, i] = minimize(hyp, 'wgp_lik', length, weight, input, target)
%
% See also: wgp_pred
%      
% (C) Copyright 1999, 2000 & 2001, Carl Edward Rasmussen (2001-05-25).
% Editted by David Kim 2016

[N, D] = size(input);       % number of examples and dimension of input space
exphyp = exp(hyp);          % exponentiate the hyperparameters once and for all
W = diag(weight);       


% first, we write out the covariance matrix C_N

Z = zeros(N,N);
for d = 1:D
  Z = Z + (repmat(input(:,d),1,N)-repmat(input(:,d)',N,1)).^2*exphyp(d);
end
Z = exphyp(D+1)*exp(-0.5*Z);
C_N = Z + exphyp(D+2)*W^-1;        % Gaussian term, weight, and noise

% then, we compute the negative log likelihood ...

invCN = inv(C_N);
invCNt = invCN*target;
logdetCN = 2*sum(log(diag(chol(C_N))));            % don't compute det(CN) directly
fX = 0.5*logdetCN + 0.5*target'*invCNt + 0.5*N*log(2*pi);

% ... and its partial derivatives

dfX = zeros(D+2,1);                     % set the size of the derivative vector

for d = 1:D % for every length scale
  V = (repmat(input(:,d),1,N)-repmat(input(:,d)',N,1)).^2.*Z;
  dfX(d) = exphyp(d)*(invCNt'*V*invCNt - sum(sum(invCN.*V)))/4;
end 
dfX(D+1) = 0.5*sum(sum(invCN.*Z)) - 0.5*invCNt'*Z*invCNt; % derivative of v1, sigma_f^2
dfX(D+2) = 0.5*trace(invCN*W^-1)*exphyp(D+2) - 0.5*invCNt'*W^-1*invCNt*exphyp(D+2); % derivative of v0, sigma_n^2 = lambda
%dfX(D+3) = 0.5*trace(invCN*W^-1)*exphyp(D+3) - 0.5*invCNt'*W^-1*invCNt*exphyp(D+3); % derivative of v0, lambda