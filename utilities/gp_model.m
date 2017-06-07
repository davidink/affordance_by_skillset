function [mean sigsq] = gp_model(X,A,Y,l_scale,signal_sd,sigma)

% X = new input
% A = training input
% Y = training output
% l_scale, signal_sd = hyperparameter for kernel
% sigma = noise sd for output

for i=1:size(A,1)
    for j=i:size(A,1)
        K_AA(i,j) = signal_sd^2 * exp ( (A(i,:)-A(j,:))*(A(i,:)-A(j,:))'/(-2 * l_scale^2) );
        K_AA(j,i) = K_AA(i,j);
    end
end

for i=1:size(A,1)
    K_XA(i,1) = signal_sd^2 * exp(( (X-A(i,:))*(X-A(i,:))' )/(-2*l_scale^2));
end

K_XX = signal_sd^2;

mean = 
sigsq = K_XX - K_XA' * (K_AA^-1 + sigma*eye(size(A,1))) * K_XA;
end
