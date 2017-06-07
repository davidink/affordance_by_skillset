rng(0,'twister'); % For reproducibility
n = 1000;
x = linspace(-10,10,n)';
y = 1 + x*5e-2 + sin(x)./x + 0.2*randn(n,1);

gprMdl = fitrgp(x,y,'Basis','linear',...
'FitMethod','exact','PredictMethod','exact');
ypred = resubPredict(gprMdl);

x2 = linspace(10,20,n)';
ypred2 = predict(gprMdl,x2);


plot(x,y,'b.');
hold on;
plot(x,ypred,'r','LineWidth',1.5);
plot(x2,ypred2,'b','LineWidth',1.5);
xlabel('x');
ylabel('y');
legend('Data','GPR predictions');
hold off