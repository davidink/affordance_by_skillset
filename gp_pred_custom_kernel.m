function [mu, S2] = gp_pred_custom_kernel(hyp, input, target, test);

% where: 
%
%   hyp     is a (column) vector (of size D+3) of hyperparameters
%   weight  is a column vector of size N of inputs
%   input   is a n by D matrix of training inputs
%   target  is a (column) vector (of size n) of targets
%   test    is a nn by D matrix of test inputs
%   mu      is a (column) vector (of size nn) of prediced means
%   S2      is a (column) vector (of size nn) of predicted variances

% Note that the reported variances in S2 are for the noise-free signal; to get
% the noisy variance, simply add the noise variance log(X(D+3)).
%
% C(x^p,x^q) = v1 * exp( -0.5 * sum_{d=1..D} w_d * (x^p_d - x^q_d)^2 )
%            + v0 * W^-1
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
%       log(v0)  % noise signal variance sigma_n^2
%       ]
% 

[N, D] = size(input);   % number of training cases and dimension of input space
[nn, D] = size(test);       % number of test cases and dimension of input space
exphyp = exp(hyp);              % exponentiate the hyperparameters once and for all
%W = diag(weight);

% first, we write out the covariance matrix C_N for the training inputs ...

C_N = zeros(N,N);                                         % create and zero space
for i=1:N
    for j=i:N
        C_N(i,j) = compute_kernel_custom(input(i,:),input(j,:),exphyp(1:4));
        C_N(j,i) = C_N(i,j);
    end
end
C_N = C_N + exphyp(5)*eye(N);
% for d = 1:D                                           % non-linear contribution
%   C_N = C_N + exphyp(d)*(repmat(input(:,d),1,N)-repmat(input(:,d)',N,1)).^2;
% end
% C_N = exphyp(D+1)*exp(-0.5*C_N) + exphyp(D+2)*W^-1;


% ... then we compute the covariance between training and test inputs ...

k = zeros(N, nn);                                       % create and zero space
for i=1:N
    for j=1:nn
        k(i,j) = compute_kernel_custom(input(i,:),test(j,:),exphyp(1:4));
    end
end
% for d = 1:D
%   %k = k + exphyp(d)*(repmat(input(:,d),1,nn)-repmat(test(:,d)',N,1)).^2;
% end
% k = exphyp(D+1)*exp(-0.5*k);

% ... and covariance between the test input and themselves 

k_s = compute_kernel_custom(test(1,:),test(1,:),exphyp(1:4));
%k_s = exphyp(D+1);

% Now, write out the desired terms

if nargout == 1
  mu = k'*(C_N\target);     % don't compute invQ explicitly if we only need means
else
  invCN = inv(C_N);
  invCNt = invCN*target;
  mu = k'*invCNt;                                              % predicted means
  S2 = k_s + exphyp(5) - sum(k.*(invCN*k),1)';                            % predicted variance
end
% 
% if nargout > 2
%   deriv = zeros(nn,D);                                           % create space
%   if nargout == 4
%     S2deriv = zeros(nn,D);
%   elseif nargout == 5
%     S2deriv = zeros(nn,D,D); dummy = [];        % assign dummy to avoid warning
%   end
%   for d = 1:D
%     c = a.*(repmat(input(:,d),1,nn)-repmat(test(:,d)',n,1));
%     deriv(1:nn,d) = expX(d)*c'*invQt;                      % derivative of mean
%     if nargout == 4
%       S2deriv(1:nn,d) = expX(d)*(expX(D+1)-expX(d)*sum(c.*(invQ*c),1)');    
%     elseif nargout == 5
%       ainvQc = a.*(invQ*c);
%       for e = 1:D
%         S2deriv(1:nn,d,e) = -expX(d)*expX(e)* ...
%               sum(ainvQc.*(repmat(input(:,e),1,nn)-repmat(test(:,e)',n,1)),1)';
%       end
%       S2deriv(1:nn,d,d) = S2deriv(1:nn,d,d) + expX(d)*expX(D+1);
%     end
%   end
% end