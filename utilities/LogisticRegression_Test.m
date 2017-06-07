load('fisheriris.mat');

sp = nominal(species);
sp = double(sp);

[B,dev,stats] = mnrfit(meas,sp);

x = [6.2, 3.7, 5.8, 0.2];
pihat = mnrval(B,x);

y_train = zeros(size(sp,1),max(sp));

for i=1:size(y_train)
    y_train(i,sp(i)) = 1;
end

myClassifier=LogisticRegression(size(meas,2),3,'linearFeatures');
myClassifier.Train(meas,y_train,1e-6);
%y_pred = myClassifier.computeExpected(meas(140,:))
%y_sample = myClassifier.Sample(x)
y_pred = myClassifier.computeLikelihood(meas,[]);

%pihat = mnrval(B,meas(2,:))


