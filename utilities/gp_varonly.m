function sigsq = gp_varonly(X,A,l_scale,signal_sd)

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

sigsq = K_XX - K_XA' * K_AA^-1 * K_XA;
end
