clear

A = [-1.5,-1.0,-0.75,-0.4,-0.25,0.00];
X = 0.2;
l_scale = 1;
signal_sd = 1.27;
sigma = 0.3;

for i=1:size(A,2)
    for j=i:size(A,2)
        K_AA(i,j) = signal_sd^2 * exp( (A(i)-A(j))'*(A(i)-A(j)) / (-2*l_scale^2) );
        if i==j
            K_AA(i,j) = K_AA(i,j) + sigma^2;
        end
        K_AA(j,i) = K_AA(i,j);
    end
end

for i=1:size(A,2)
    K_XA(i) = signal_sd^2 * exp( (X-A(i))'*(X-A(i)) / (-2*l_scale^2) );
end
K_XX = signal_sd^2 + sigma^2;

mean_y = 3;
var_y = K_XX - K_XA * K_AA^-1 * K_XA';

% clear
% 
% cd(matlabroot)
% cd('help/toolbox/stats/examples')
% 
% load gprdata
% 
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
% gprMdl.KernelInformation
% 
% gprMdl.KernelInformation.KernelParameterNames
% 
% sigmaM = gprMdl.KernelInformation.KernelParameters(1:end-1,1)
% sigmaF = gprMdl.KernelInformation.KernelParameters(end)
% sigma  = gprMdl.Sigma
% 
% figure()
% plot((1:d)',log(sigmaM),'ro-');
% xlabel('Length scale number');
% ylabel('Log of length scale');
% 
% X = [Xtrain(:,1:3) Xtrain(:,6)];
% sigma0 = std(ytrain);
% sigmaF0 = sigma0;
% d = size(X,2);
% sigmaM0 = 10*ones(d,1);
% 
% gprMdl = fitrgp(X,ytrain,'Basis','constant','FitMethod','exact',...
% 'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
% 'KernelParameters',[sigmaM0;sigmaF0],'Sigma',sigma0,'Standardize',1);
% 
% xtest = [Xtest(:,1:3) Xtest(:,6)];
% L = loss(gprMdl,xtest,ytest)
% 
% ypred = predict(gprMdl,xtest);
% 
% figure;
% plot(ytest,'r');
% hold on;
% plot(ypred,'b');
% legend('True response','GPR predicted values','Location','Best');
% hold off